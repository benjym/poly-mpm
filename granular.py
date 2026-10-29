"""Experimental inertial, rigid-grain MPM with coupled nonsmooth friction.

Inspired by Daviet & Bertails-Descoubes (2016), equations (7)--(14).
This independent small-problem implementation uses the existing nodal bases,
intrinsic 2D symmetric tensors and lumped mass. The default small-problem hybrid
solve has a dense Jacobian and sparse block Gauss-Seidel fallback; a block-only
mode avoids the global factorization. It is NOT a reproduction of Sand6
or a claim of continuum well-posedness. See docs/granular-prototype.md.
The opt-in velocity solver uses compact gradients and the local iteration from
Daviet (2026); see granular_mixed.py and docs/granular-mixed.md.
"""

import numpy as np
from scipy import sparse
from scipy.optimize import brentq, least_squares
from scipy.sparse.linalg import spsolve

from interpolation import stencil
from kinematics import advance
from transfers import gather, nodal_velocity, scatter

# Frobenius-orthonormal intrinsic 2D stress basis, embedded for output.
BASIS = np.array(
    [
        np.diag([1.0, 1.0, 0.0]) / np.sqrt(2),
        np.diag([1.0, -1.0, 0.0]) / np.sqrt(2),
        [[0.0, 1 / np.sqrt(2), 0.0], [1 / np.sqrt(2), 0.0, 0.0], [0.0, 0.0, 0.0]],
    ]
)
COMPONENTS = len(BASIS)
PLANE_BASIS = np.array(
    [
        np.eye(3) / np.sqrt(3),
        np.diag([1.0, -1.0, 0.0]) / np.sqrt(2),
        np.diag([-1.0, -1.0, 2.0]) / np.sqrt(6),
        BASIS[2],
    ]
)


def stress_space(P):
    """Plane strain retains the 3D trace/deviator, with Dzz=Dxz=Dyz=0."""
    return (PLANE_BASIS, 3) if P.granular_stress == "plane_strain" else (BASIS, 2)


def viscosity_matrix(P, G, batch):
    """2 eta sum(Vgrain D(u):D(v)), matching Sand6's trial solve."""
    nodes = canonical_nodes(P, len(G.m))[batch.n]
    gradients = batch.G[:, :, :2]
    weight = P.granular_viscosity * batch.m / batch.material.rho_s
    dot = np.einsum("pia,pja->pij", gradients, gradients)
    rows, cols, data = [], [], []
    for a in range(2):
        for b in range(2):
            value = weight[:, None, None] * (
                (dot if a == b else 0) + gradients[:, :, b, None] * gradients[:, None, :, a]
            )
            rows.append(np.broadcast_to(2 * nodes[:, :, None] + a, value.shape).ravel())
            cols.append(np.broadcast_to(2 * nodes[:, None, :] + b, value.shape).ravel())
            data.append(value.ravel())
    result = sparse.coo_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(2 * len(G.m), 2 * len(G.m)),
    ).tocsr()
    result.eliminate_zeros()
    return result


def cone_projection(x, mu):
    """Euclidean projection and a generalized derivative for |r_t| <= mu*r_n."""
    x = np.asarray(x)
    components = x.shape[-1]
    x = x.reshape(-1, components)
    out = np.zeros_like(x)
    jac = np.zeros((len(x), components, components))
    for i, (value, friction) in enumerate(zip(x, np.broadcast_to(mu, len(x)))):
        n, t = value[0], value[1:]
        length = np.linalg.norm(t)
        if friction == 0:
            out[i, 0] = max(n, 0)
            jac[i, 0, 0] = float(n > 0)
        elif length <= friction * n:
            out[i] = value
            jac[i] = np.eye(components)
        elif n + friction * length > 0:
            direction = t / length
            normal = (n + friction * length) / (1 + friction**2)
            out[i, 0] = normal
            out[i, 1:] = friction * normal * direction
            jac[i, 0, 0] = 1 / (1 + friction**2)
            jac[i, 0, 1:] = jac[i, 1:, 0] = friction * direction / (1 + friction**2)
            outer = np.outer(direction, direction)
            jac[i, 1:, 1:] = friction**2 / (1 + friction**2) * outer + friction * normal / length * (
                np.eye(components - 1) - outer
            )
    return out, jac


def coulomb_mfb(reaction, velocity, mu):
    """SOC Fischer-Burmeister residual with the nonassociated Coulomb shift.

    The reciprocal cone scaling and de Saxce velocity correction reduce
    Signorini-Coulomb to self-dual second-order-cone complementarity.
    mu must be strictly positive. This is used as an independent merit map.
    """
    mu = np.broadcast_to(mu, reaction.shape[:-1])
    a, b = reaction.copy(), velocity.copy()
    a[..., 0] *= mu
    b[..., 0] = b[..., 0] / mu + np.linalg.norm(b[..., 1:], axis=-1)
    normal = np.sum(a * a + b * b, axis=-1)
    tangent = 2 * (a[..., :1] * a[..., 1:] + b[..., :1] * b[..., 1:])
    length = np.linalg.norm(tangent, axis=-1)
    plus, minus = np.sqrt(np.maximum(normal + length, 0)), np.sqrt(np.maximum(normal - length, 0))
    root = np.zeros_like(a)
    root[..., 0] = 0.5 * (plus + minus)
    factor = np.divide(plus - minus, 2 * length, out=np.zeros_like(length), where=length > 0)
    root[..., 1:] = factor[..., None] * tangent
    return a + b - root


def _local_coulomb(W, q, mu, previous):
    """Exact 2D local contact update for nonlinear block Gauss-Seidel."""
    if q[0] >= 0:
        return np.zeros(3)
    stuck = np.linalg.lstsq(W, -q, rcond=1e-12)[0]
    scale = max(np.linalg.norm(q), 1e-30)
    if (
        stuck[0] >= 0
        and np.linalg.norm(stuck[1:]) <= mu * stuck[0]
        and np.linalg.norm(W @ stuck + q) <= 1e-10 * scale
    ):
        return stuck
    if mu == 0:
        if W[0, 0] <= 0:
            raise RuntimeError("packing constraint has no normal support")
        return np.array([-q[0] / W[0, 0], 0.0, 0.0])

    def sliding(angle):
        direction = np.array([np.cos(angle), np.sin(angle)])
        perpendicular = np.array([-direction[1], direction[0]])
        denominator = W[0, 0] - mu * (W[0, 1:] @ direction)
        if denominator <= 0:
            return None
        normal = -q[0] / denominator
        value = np.r_[normal, -mu * normal * direction]
        velocity = W @ value + q
        dn = normal * mu * (W[0, 1:] @ perpendicular) / denominator
        dv = dn * (W[1:, 0] - mu * W[1:, 1:] @ direction) - mu * normal * (W[1:, 1:] @ perpendicular)
        f = perpendicular @ velocity[1:]
        derivative = -direction @ velocity[1:] + perpendicular @ dv
        return value, f, derivative, direction @ velocity[1:]

    direction = (
        -previous[1:] if np.linalg.norm(previous[1:]) else q[1:] - W[1:, 0] * q[0] / max(W[0, 0], 1e-30)
    )
    angle = np.arctan2(direction[1], direction[0])
    for _ in range(15):
        candidate = sliding(angle)
        if candidate is None:
            break
        value, f, df, slip = candidate
        if abs(f) <= 1e-11 * scale and slip >= -1e-11 * scale:
            return value
        if abs(df) <= 1e-30:
            break
        angle -= np.clip(f / df, -0.5, 0.5)
    # Bracket every orientation to select the sliding, not reversed-slip root.
    angles = np.linspace(-np.pi, np.pi, 65)
    for left, right in zip(angles[:-1], angles[1:]):
        a, b = sliding(left), sliding(right)
        if a is None or b is None or a[1] * b[1] > 0:
            continue
        root = brentq(lambda x: sliding(x)[1], left, right, xtol=1e-13)
        value, f, df, slip = sliding(root)
        if slip >= -1e-10 * scale:
            return value
    raise RuntimeError("granular local sliding problem could not be resolved")


def solve_coulomb(W, q, mu, *, tolerance=1e-9, max_iterations=200, initial=None, method="hybrid"):
    """Solve w=W*r+q and non-associated Coulomb complementarity.

    r_n >= 0, w_n >= 0, r_n*w_n=0; |r_t|<=mu*r_n; sliding
    r_t=-mu*r_n*w_t/|w_t|. Normal complementarity and the tangential friction disk are enforced
    separately: an uncorrected joint cone projection would give associated dilation.
    Input friction is the tensor-cone coefficient (equal to tau_eq/p in 2D).
    """
    W = sparse.csr_matrix(W)
    q = np.asarray(q, dtype=float)
    components = 4 if q.ndim == 1 else q.shape[-1]
    q = q.reshape(-1, components)
    mu = np.broadcast_to(np.asarray(mu, dtype=float), len(q))
    if method not in ("hybrid", "block") or (method == "block" and components != 3):
        raise ValueError("method must be hybrid, or block with three contact components")
    if W.shape != (q.size, q.size) or not np.isfinite(W.data).all() or not np.isfinite(q).all():
        raise ValueError("invalid Coulomb problem dimensions or nonfinite data")
    if (
        np.any(~np.isfinite(mu))
        or np.any(mu < 0)
        or not np.isfinite(tolerance)
        or tolerance <= 0
        or not np.isfinite(max_iterations)
        or max_iterations < 1
        or int(max_iterations) != max_iterations
    ):
        raise ValueError("invalid friction or solver controls")
    if not len(q):
        return q.copy(), dict(residual=0.0, evaluations=0, converged=True)
    norm = float(np.max(np.asarray(abs(W).sum(axis=1)), initial=0))
    if norm == 0:
        if np.any(q[:, 0] < 0):
            raise RuntimeError("packing constraint has no free velocity support")
        return np.zeros_like(q), dict(residual=0.0, evaluations=0, converged=True)
    # Dimensionless unknowns and residual, with no constitutive compliance.
    # Large positive packing clearances are already feasible. Including them
    # in the residual scale could hide small compressive violations elsewhere
    # (especially as dt -> 0, when the clearance grows like 1/dt).
    magnitude = max(float(np.max(-q[:, 0], initial=0)), float(np.max(np.abs(q[:, 1:]), initial=0)))
    stress_scale = max(magnitude, 1e-30) / norm
    operator = W / norm
    rhs = q / max(magnitude, 1e-30)
    identity = sparse.eye(q.size, format="csr")

    def evaluate(z, jacobian=False):
        reaction = z.reshape(-1, components)
        w = (operator @ z).reshape(-1, components) + rhs
        argument = reaction - w
        radius = mu * np.maximum(reaction[:, 0], 0)
        length = np.linalg.norm(argument[:, 1:], axis=1)
        factor = np.divide(radius, length, out=np.ones_like(radius), where=length > 0)
        projected = argument.copy()
        projected[:, 0] = np.maximum(argument[:, 0], 0)
        projected[:, 1:] *= np.minimum(1, factor)[:, None]
        residual = reaction - projected
        if not jacobian:
            return residual.ravel()
        derivative = np.zeros((len(q), components, components))
        radius_derivative = np.zeros_like(derivative)
        derivative[:, 0, 0] = argument[:, 0] > 0
        for i in range(len(q)):
            if length[i] <= radius[i]:
                derivative[i, 1:, 1:] = np.eye(components - 1)
            else:
                unit = argument[i, 1:] / length[i]
                derivative[i, 1:, 1:] = (
                    radius[i] / length[i] * (np.eye(components - 1) - np.outer(unit, unit))
                )
                if reaction[i, 0] >= 0:
                    radius_derivative[i, 1:, 0] = mu[i] * unit

        def blocks(values):
            return sparse.bsr_matrix(
                (values, np.arange(len(q)), np.arange(len(q) + 1)), shape=W.shape
            ).tocsr()

        return identity - blocks(derivative) @ (identity - operator) - blocks(radius_derivative)

    start = (
        np.zeros(q.size) if initial is None or magnitude == 0 else np.asarray(initial).ravel() / stress_scale
    )
    candidate = start.copy()
    evaluations = 0
    if method == "hybrid":
        solution = least_squares(
            evaluate,
            start,
            jac=lambda z: evaluate(z, True).toarray(),
            tr_solver="exact",
            ftol=1e-13,
            xtol=1e-13,
            gtol=1e-13,
            max_nfev=max_iterations,
        )
        candidate = solution.x
        evaluations = solution.nfev
    residual = float(np.max(np.abs(evaluate(candidate)), initial=0))
    sweeps = 0
    mfb_evaluations = 0
    if residual > tolerance and components == 4 and np.all(mu > 0):
        # A different equivalent merit map can escape a nonzero minimum of
        # the natural-map least squares in the dependent plane-strain space.
        def mfb(z):
            velocity = (operator @ z).reshape(-1, components) + rhs
            return coulomb_mfb(z.reshape(-1, components), velocity, mu).ravel()

        alternative = least_squares(mfb, start, ftol=1e-13, xtol=1e-13, gtol=1e-13, max_nfev=max_iterations)
        alternative_residual = float(np.max(np.abs(evaluate(alternative.x)), initial=0))
        mfb_evaluations = alternative.nfev
        if alternative_residual < residual:
            candidate, residual = alternative.x, alternative_residual
    if residual > tolerance and components == 3:
        # Newton/least-squares can settle at a nonzero local minimum of a
        # non-associated residual. A block contact solve is an independent
        # fallback, not acceptance of that stationary point.
        candidate = start.copy()
        # Store sparse 3x3 columns. A local change only updates neighbouring
        # contacts, without an O(nodes^2) dense Delassus matrix.
        columns = operator.T.tobsr(blocksize=(3, 3))
        columns.sum_duplicates()
        columns.sort_indices()
        local_blocks = [operator[3 * i : 3 * i + 3, 3 * i : 3 * i + 3].toarray() for i in range(len(q))]
        velocity = (operator @ candidate).reshape(-1, 3) + rhs
        for sweeps in range(1, max_iterations * 10 + 1):
            for i in range(len(q)):
                sl = slice(3 * i, 3 * i + 3)
                local = local_blocks[i]
                other = velocity[i] - local @ candidate[sl]
                updated = _local_coulomb(local, other, mu[i], candidate[sl])
                delta = updated - candidate[sl]
                candidate[sl] = updated
                lo, hi = columns.indptr[i : i + 2]
                velocity[columns.indices[lo:hi]] += columns.data[lo:hi].transpose(0, 2, 1) @ delta
            # Avoid accumulating incremental-update roundoff over many sweeps.
            velocity = (operator @ candidate).reshape(-1, 3) + rhs
            residual = float(np.max(np.abs(evaluate(candidate)), initial=0))
            if residual <= tolerance:
                break

    if not np.isfinite(residual) or residual > tolerance:
        raise RuntimeError(
            f"granular Coulomb solve did not converge: residual={residual:.3g}, "
            f"evaluations={evaluations}, sweeps={sweeps}; reduce dt or inspect particle support"
        )
    reaction = candidate.reshape(-1, components) * stress_scale
    return reaction, dict(
        residual=residual,
        evaluations=evaluations,
        sweeps=sweeps,
        converged=True,
        method=method,
        mfb_evaluations=mfb_evaluations,
        operator_nnz=W.nnz,
        operator_storage_bytes=W.data.nbytes + W.indices.nbytes + W.indptr.nbytes,
    )


def canonical_nodes(P, count):
    nodes = np.arange(count)
    if P.B.cyclic_lr:
        nodes[P.G.nx - 1 :: P.G.nx] -= P.G.nx - 1
    return nodes


def basis_volumes(P, G):
    """Integrate the actual nodal basis, including modified quadratic edges."""
    key = (P.interpolation, P.B.cyclic_lr, P.G.nx, P.G.ny, P.G.thickness, P.G.x_m, P.G.x_M, P.G.y_m, P.G.y_M)
    cached = getattr(G, "_granular_basis_volumes", None)
    if cached is not None and cached[0] == key:
        return cached[1]
    xi, wi = np.polynomial.legendre.leggauss(3)
    axes, weights = [], []
    for lo, hi, n in ((P.G.x_m, P.G.x_M, P.G.nx), (P.G.y_m, P.G.y_M, P.G.ny)):
        edges = np.linspace(lo, hi, 2 * (n - 1) + 1)
        half = np.diff(edges) / 2
        axes.append(((edges[:-1] + half)[:, None] + half[:, None] * xi).ravel())
        weights.append((half[:, None] * wi).ravel())
    xx, yy = np.meshgrid(*axes)
    points = np.column_stack((xx.ravel(), yy.ravel(), np.zeros(xx.size)))
    volume = (weights[1][:, None] * weights[0][None, :]).ravel() * P.G.thickness
    nodes, shape, _, _ = stencil(points, P, G)
    nodes = canonical_nodes(P, len(G.m))[nodes]
    result = np.bincount(nodes.ravel(), weights=(volume[:, None] * shape).ravel(), minlength=len(G.m))
    G._granular_basis_volumes = key, result
    return result


def packing_fraction(P, G, batch):
    """Reconstruct nodal packing from current positions, independently of F/J.

    Call after boundary processing, so positions are inside the grid and any
    periodic coordinates have been wrapped. No particle state is changed.
    """
    nodes, shape, _, _ = stencil(batch.x, P, G)
    nodes = canonical_nodes(P, len(G.m))[nodes]
    occupied = np.bincount(
        nodes.ravel(),
        weights=((batch.m / batch.material.rho_s)[:, None] * shape).ravel(),
        minlength=len(G.m),
    )
    volumes = basis_volumes(P, G)
    return np.divide(occupied, volumes, out=np.zeros_like(occupied), where=volumes > 0)


def coarse_factor(P):
    return int(getattr(P, "granular_stress_coarsening", 2))


def coarse_shape(P):
    """Coarse stress grid: every k-th velocity node in each direction (k = 2 by default)."""
    k = coarse_factor(P)
    if (P.G.nx - 1) % k or (P.G.ny - 1) % k:
        raise ValueError(f"coarse stress basis requires cell counts divisible by {k} in x and y")
    return (P.G.nx - 1) // k + 1, (P.G.ny - 1) // k + 1


def coarse_canonical(P):
    nx, ny = coarse_shape(P)
    nodes = np.arange(nx * ny)
    if P.B.cyclic_lr:
        nodes[nx - 1 :: nx] -= nx - 1
    return nodes


def coarse_stencil(P, G, positions):
    """Bilinear shape functions of the coarse stress grid (spacing k dx, k dy).

    Returns canonical coarse nodes and weights in the corner order of
    interpolation.stencil, including the closed upper edge and periodic wrap.
    """
    nx, ny = coarse_shape(P)
    x = positions[:, 0]
    if P.B.cyclic_lr:
        x = P.G.x_m + (x - P.G.x_m) % (P.G.x_M - P.G.x_m)
    k = coarse_factor(P)
    local_x = (x - P.G.x_m) / (k * G.dx)
    local_y = (positions[:, 1] - P.G.y_m) / (k * G.dy)
    cx = np.minimum(np.floor(local_x).astype(int), nx - 2)
    cy = np.minimum(np.floor(local_y).astype(int), ny - 2)
    fx, fy = local_x - cx, local_y - cy
    base = cx + nx * cy
    nodes = base[:, None] + [0, 1, nx + 1, nx]
    weights = np.stack(((1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy), axis=1)
    return coarse_canonical(P)[nodes], weights


def coarse_volumes(P, G):
    """Exact integrals of the coarse bilinear basis over the domain."""
    nx, ny = coarse_shape(P)
    k = coarse_factor(P)
    wx = np.full(nx, k * G.dx)
    wy = np.full(ny, k * G.dy)
    wy[[0, -1]] = k * G.dy / 2
    if P.B.cyclic_lr:
        wx[-1] = 0  # alias column; its canonical node carries the full support
    else:
        wx[[0, -1]] = k * G.dx / 2
    return (wy[:, None] * wx[None, :]).ravel() * P.G.thickness


def cell_indices(P, G, positions):
    """Q1 cell ownership, including the closed upper edge and periodic wrap."""
    x = positions[:, 0]
    if P.B.cyclic_lr:
        x = P.G.x_m + (x - P.G.x_m) % (P.G.x_M - P.G.x_m)
    ix = np.minimum(np.floor((x - P.G.x_m) / G.dx).astype(int), P.G.nx - 2)
    iy = np.minimum(np.floor((positions[:, 1] - P.G.y_m) / G.dy).astype(int), P.G.ny - 2)
    return ix + (P.G.nx - 1) * iy


def cell_packing_fraction(P, G, batch):
    """Cell-average fixed grain volume; distinct from nodal/F-derived packing."""
    volume = G.dx * G.dy * P.G.thickness
    return (
        np.bincount(
            cell_indices(P, G, batch.x),
            weights=batch.m / batch.material.rho_s,
            minlength=(P.G.nx - 1) * (P.G.ny - 1),
        )
        / volume
    )


def assemble(P, G, batch, *, compact=False):
    """Build b(tau,u)=sum_p V_grain*N_stress*tau:D(u) and packing gap."""
    basis, dimension = stress_space(P)
    components = len(basis)
    canonical = canonical_nodes(P, len(G.m))
    nodes = canonical[batch.n]
    valid = G.valid_nodes[canonical]
    if np.any((~valid[nodes]) & (np.linalg.norm(batch.G, axis=2) > 1e-12)):
        raise ValueError("granular solver requires mass on every derivative-support node")
    solid_volume = batch.m / batch.material.rho_s
    if P.granular_stress_basis == "P0":
        stress_nodes = cell_indices(P, G, batch.x)[:, None]
        stress_weights = np.ones_like(stress_nodes, dtype=float)
        volumes = np.full((P.G.nx - 1) * (P.G.ny - 1), G.dx * G.dy * P.G.thickness)
        unique = np.ones(len(volumes), dtype=bool)
    elif P.granular_stress_basis == "coarse":
        stress_nodes, stress_weights = coarse_stencil(P, G, batch.x)
        volumes = coarse_volumes(P, G)
        unique = coarse_canonical(P) == np.arange(len(volumes))
    else:
        stress_nodes, stress_weights = nodes, batch.N
        volumes = basis_volumes(P, G)
        unique = canonical == np.arange(len(G.m))
    occupied = np.bincount(
        stress_nodes.ravel(), weights=(solid_volume[:, None] * stress_weights).ravel(), minlength=len(volumes)
    )
    active = np.flatnonzero(unique & (occupied > 0))
    if len(active) > P.granular_max_nodes:
        raise ValueError(f"granular prototype supports at most {P.granular_max_nodes} active stress nodes")
    index = np.full(len(volumes), -1, dtype=int)
    index[active] = np.arange(len(active))
    # Each row has all tensor components at one stress node; integration and
    # momentum use exactly the same B and its transpose.
    count, support = nodes.shape
    shape = (count, stress_nodes.shape[1], support)
    rows, cols, data = [], [], []
    for mode in range(0 if compact else components):
        for component in range(2):
            derivative = np.einsum("j,prj->pr", basis[mode, component], batch.G)
            contribution = solid_volume[:, None, None] * stress_weights[:, :, None] * derivative[:, None, :]
            row = np.broadcast_to(components * index[stress_nodes][:, :, None] + mode, shape)
            col = np.broadcast_to(2 * nodes[:, None, :] + component, shape)
            keep = contribution != 0
            rows.append(row[keep])
            cols.append(col[keep])
            data.append(contribution[keep])
    if compact:
        from granular_mixed import GradientOperator

        contribution = (
            solid_volume[:, None, None, None] * stress_weights[:, :, None, None] * batch.G[:, None, :, :2]
        )
        B = GradientOperator(
            np.broadcast_to(index[stress_nodes][:, :, None], shape).ravel(),
            np.broadcast_to(nodes[:, None, :], shape).ravel(),
            contribution.reshape(-1, 2),
            len(active),
            len(G.m),
            basis,
        )
    else:
        B = sparse.coo_matrix(
            (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
            shape=(components * len(active), 2 * len(G.m)),
        ).tocsr()
    free = np.repeat(G.valid_nodes & (canonical == np.arange(len(G.m))), 2)
    free[2 * np.flatnonzero(G.boundary_v)] = False
    free[2 * np.flatnonzero(G.boundary_h) + 1] = False
    if P.B.roughness:
        free[2 * np.flatnonzero(G.boundary_h)] = False
        free[2 * np.flatnonzero(G.boundary_v) + 1] = False
    gap = np.zeros((len(active), components))
    capacity = batch.material.packing_max * volumes[active]
    # Positive-part gap follows the paper's dissipative packing linearisation.
    # It does not geometrically repair pre-existing overpacking.
    clearance = capacity - occupied[active]
    # Do not amplify quadrature/summation roundoff into an artificial opening
    # velocity as dt becomes small. Resolve meaningful physical clearance.
    roundoff = 64 * np.finfo(float).eps * np.maximum(capacity, occupied[active])
    clearance[np.abs(clearance) <= roundoff] = 0
    gap[:, 0] = np.maximum(clearance, 0) / (np.sqrt(dimension) * P.dt)
    return B, free, active, occupied, gap


def step(P, G, L):
    """Advance the experimental rigid-grain branch using end-of-step stresses."""
    from stability import limits
    from diagnostics import inventory

    if not np.isfinite(P.dt) or P.dt <= 0:
        raise ValueError("dt must be finite and positive")
    P.timestep_limits = limits(P, G, L)
    if P.timestep_policy == "enforce" and P.dt > min(P.timestep_limits.values(), default=P.dt) * (1 + 1e-12):
        raise ValueError("granular dt exceeds transport timestep estimate")
    if L.diagnostics is not None:
        L.diagnostics.begin(P, L)
    G.wipe(P)
    P.update_forces()
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    if P.O.check_positions:
        L.recover_position(P, G)
    L.get_nodal_mass_momentum(P, G)
    b = L.materials[0]
    basis, dimension = stress_space(P)
    components = len(basis)
    body = np.array([P.g * np.sin(P.theta), P.g * np.cos(P.theta), 0.0])
    scatter(G.fe, b.n, b.N[:, :, None] * b.m[:, None, None] * body)
    if P.B.cyclic_lr:
        G.make_cyclic(P, G, ["fe"])
    G.BCs(P)
    G.update_momentum(P)
    trial = nodal_velocity(G, P)[:, :2].ravel()
    B, free, active, occupied, gap = assemble(P, G, b, compact=P.granular_solver == "velocity")
    inverse_mass = 1 / np.repeat(G.m, 2)[free]
    viscous_force = np.zeros_like(trial)
    viscous_dissipation = 0.0
    if P.granular_viscosity:
        stiffness = viscosity_matrix(P, G, b)
        mass = 1 / inverse_mass
        system = sparse.diags(mass) + P.dt * stiffness[free][:, free]
        trial[free] = spsolve(system, mass * trial[free])
        if P.B.cyclic_lr:
            periodic_velocity = trial.reshape(-1, 2)
            periodic_velocity[P.G.nx - 1 :: P.G.nx] = periodic_velocity[:: P.G.nx]
        viscous_force = stiffness @ trial
        viscous_dissipation = float(P.dt * trial @ viscous_force)
    old = getattr(G, "granular_contact", np.zeros((len(occupied), components)))
    options = dict(
        tolerance=P.granular_tolerance, max_iterations=P.granular_max_iterations, initial=old[active]
    )
    mu = b.material.mu * np.sqrt(2 / dimension)
    if P.granular_solver == "velocity":
        from granular_mixed import solve_velocity

        full_inverse_mass = np.zeros_like(trial)
        full_inverse_mass[free] = inverse_mass
        contact, stats = solve_velocity(
            B,
            trial,
            full_inverse_mass,
            gap,
            mu,
            P.dt,
            stopping=P.granular_stopping,
            increment=P.granular_increment,
            **options,
        )
    else:
        operator = B[:, free]
        W = (operator.multiply(inverse_mass) @ operator.T) * P.dt
        q = (B @ trial).reshape(-1, components) + gap
        contact, stats = solve_coulomb(W, q, mu, method=P.granular_solver, **options)
    force = np.asarray(B.T @ contact.ravel()).reshape(-1, 2)
    final = trial.copy()
    final[free] += P.dt * inverse_mass * force.ravel()[free]
    velocity = np.column_stack((final.reshape(-1, 2), np.zeros(len(G.m))))
    G.fi[:, :2] = -force + viscous_force.reshape(-1, 2)
    G.q = G.m[:, None] * velocity
    G.q_dot = (G.q - G.q_reference) / P.dt
    if P.B.cyclic_lr:
        for name in ("q", "q_dot", "fi"):
            value = getattr(G, name)
            value[P.G.nx - 1 :: P.G.nx] = value[:: P.G.nx]
    velocity = nodal_velocity(G, P)
    G.granular_contact = np.zeros((len(occupied), components))
    G.granular_contact[active] = contact
    if P.B.cyclic_lr and P.granular_stress_basis == "nodal":
        G.granular_contact[P.G.nx - 1 :: P.G.nx] = G.granular_contact[:: P.G.nx]
    # Recompute physical work/impulse budgets with the final contact response.
    unique = canonical_nodes(P, len(G.m)) == np.arange(len(G.m))
    midpoint = 0.5 * (
        np.divide(G.q_reference, G.m[:, None], out=np.zeros_like(G.q), where=G.m[:, None] > P.M_tol)
        + velocity
    )
    reaction = G.q - G.q_reference - P.dt * (G.fe - G.fi)
    coordinates = np.column_stack((G.X, G.Y, np.zeros_like(G.X)))
    G.constraint_impulse = reaction[unique].sum(axis=0)
    G.constraint_angular_impulse = np.cross(coordinates[unique], reaction[unique]).sum(axis=0)
    G.constraint_work = float(np.sum(reaction[unique] * midpoint[unique]))
    G.external_work = float(P.dt * np.sum(G.fe[unique] * midpoint[unique]))
    # Phase contact stresses are converted to bulk Cauchy stress for the existing
    # output API. This branch does not evolve an elastic pressure history.
    old_stress = b.stress.copy()
    old_pressure = b.pressure.copy()
    b.L = np.einsum("pri,prj->pij", velocity[b.n], b.G)
    b.dstrain = 0.5 * (b.L + b.L_T) * P.dt
    b.strain += b.dstrain
    advance(b, P)
    if P.granular_stress_basis == "P0":
        stress_coefficients = G.granular_contact[cell_indices(P, G, b.x)]
    elif P.granular_stress_basis == "coarse":
        nodes, weights = coarse_stencil(P, G, b.x)
        stress_coefficients = np.einsum("pr,prk->pk", weights, G.granular_contact[nodes])
    else:
        stress_coefficients = gather(G.granular_contact, b.n, b.N, G.valid_nodes)
    phase_stress = -np.einsum("pk,kij->pij", stress_coefficients, basis)
    phase_pressure_2d = -np.trace(phase_stress[:, :2, :2], axis1=1, axis2=2) / 2
    if dimension == 2:
        phase_stress[:, 2, 2] = -phase_pressure_2d
    b.stress = (b.rho / b.material.rho_s)[:, None, None] * phase_stress
    b.dstress = b.stress - old_stress
    b.pressure = np.trace(b.stress, axis1=1, axis2=2) / 3
    b.dp = b.pressure - old_pressure
    b.sigma_kk = 3 * b.pressure
    b.dev_stress = b.stress - b.pressure[:, None, None] * np.eye(3)
    b.dev_dstress = b.dstress - np.trace(b.dstress, axis1=1, axis2=2)[:, None, None] * np.eye(3) / 3
    rate = 0.5 * (b.L + b.L_T)
    dev = rate.copy()
    dev[:, :dimension, :dimension] -= (
        np.trace(rate, axis1=1, axis2=2)[:, None, None] * np.eye(dimension) / dimension
    )
    b.gammadot = np.sqrt(2) * np.linalg.norm(dev, axis=(1, 2))
    phase_pressure = -np.trace(phase_stress, axis1=1, axis2=2) / 3
    size = b.phi @ P.G.s
    b.I = np.divide(
        size * b.gammadot * np.sqrt(b.material.rho_s),
        np.sqrt(np.maximum(phase_pressure, 0)),
        out=np.zeros_like(phase_pressure),
        where=phase_pressure > 0,
    )
    b.mu.fill(b.material.mu)
    b.yieldfunction = np.linalg.norm(b.dev_stress, axis=(1, 2)) / np.sqrt(2) + b.material.mu * b.pressure
    b.material_work += b.V * np.sum(b.stress * b.dstrain, axis=(1, 2))
    b.plastic_work += b.V * np.sum(b.stress * b.dstrain, axis=(1, 2))
    if P.granular_viscosity:
        trial_velocity = np.column_stack((trial.reshape(-1, 2), np.zeros(len(G.m))))
        trial_gradient = np.einsum("pri,prj->pij", trial_velocity[b.n], b.G)
        trial_rate = 0.5 * (trial_gradient + trial_gradient.transpose(0, 2, 1))
        dissipation = (
            2 * P.dt * P.granular_viscosity * b.m / b.material.rho_s * np.sum(trial_rate**2, axis=(1, 2))
        )
        b.viscous_work += dissipation
        b.material_work += dissipation
    b.deposit_diagnostics(P, G)
    if P.B.cyclic_lr:
        G.make_cyclic(P, G, G.diagnostic_fields)
    G.calculate_gammadot(P, G, P.smooth_gamma_dot)
    if P.segregate:
        G.update_pk(P, G)
        G.calculate_phi_increment(P, L)
        L.move_grainsize_on_grid(P, G)
    if P.granular_stress_basis == "P0":
        initial_nodal_packing = float(np.max(packing_fraction(P, G, b), initial=0))
        stats["maximum_stress_packing"] = float(np.max(occupied, initial=0) / (G.dx * G.dy * P.G.thickness))
    elif P.granular_stress_basis == "coarse":
        initial_nodal_packing = float(np.max(packing_fraction(P, G, b), initial=0))
        volumes = coarse_volumes(P, G)
        stats["maximum_stress_packing"] = float(
            np.max(np.divide(occupied, volumes, out=np.zeros_like(occupied), where=volumes > 0), initial=0)
        )
    L.move_material_points(P, G)
    stats.update(
        active_nodes=len(active),
        stress_dimension=dimension,
        stress_basis=P.granular_stress_basis,
        viscosity_trial_dissipation=viscous_dissipation,
        contact_power=float(contact.ravel() @ (B @ velocity[:, :2].ravel())),
        maximum_packing=(
            initial_nodal_packing
            if P.granular_stress_basis in ("P0", "coarse")
            else float(
                np.max(
                    np.divide(
                        occupied,
                        basis_volumes(P, G),
                        out=np.zeros_like(occupied),
                        where=basis_volumes(P, G) > 0,
                    ),
                    initial=0,
                )
            )
        ),
    )
    G.granular_stats = stats
    if L.diagnostics is not None:
        before = inventory(P, L)
    L.process_boundaries(P, G)
    if L.diagnostics is not None:
        packing = packing_fraction(P, G, b)
        stats["maximum_packing_after"] = float(np.max(packing, initial=0))
        solid_volume = float(np.sum(b.m / b.material.rho_s))
        stats["integrated_overpacking_fraction"] = (
            float(np.maximum(packing - b.material.packing_max, 0) @ basis_volumes(P, G)) / solid_volume
            if solid_volume
            else 0.0
        )
        stats["maximum_particle_packing"] = float(np.max(b.rho / b.material.rho_s, initial=0))
        if P.granular_stress_basis == "P0":
            cell_packing = cell_packing_fraction(P, G, b)
            stats["maximum_stress_packing_after"] = float(np.max(cell_packing, initial=0))
            stats["integrated_stress_overpacking_fraction"] = (
                float(
                    np.maximum(cell_packing - b.material.packing_max, 0).sum() * G.dx * G.dy * P.G.thickness
                )
                / solid_volume
                if solid_volume
                else 0.0
            )
        L.diagnostics.finish(P, G, L, before)
