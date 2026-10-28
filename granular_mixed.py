"""Compact, velocity-based Coulomb iteration inspired by Daviet (2026), §5.

This implements the rigid, nonassociated limit. The supplied operator can use
nodal or cellwise constant stresses; finite-strain elastoplasticity is outside
this module's scope.
"""

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres


class GradientOperator:
    """B blocks reconstructed from integral(Vgrain N_stress grad N_velocity).

    Store one 2-vector per stress/velocity-node pair, including periodic sums.
    Padding makes independent stress updates amenable to batched NumPy calls.
    No tensor-expanded B or stress Delassus matrix is retained.
    """

    def __init__(self, rows, nodes, gradients, count, node_count, basis):
        keys, inverse = np.unique(rows * node_count + nodes, return_inverse=True)
        values = np.column_stack(
            [np.bincount(inverse, weights=gradients[:, c], minlength=len(keys)) for c in range(2)]
        )
        keep = np.any(values != 0, axis=1)
        keys, values = keys[keep], values[keep]
        row, node = keys // node_count, keys % node_count
        sizes = np.bincount(row, minlength=count)
        width = int(sizes.max(initial=0))
        self.nodes = np.zeros((count, width), dtype=int)
        self.gradients = np.zeros((count, width, 2))
        offsets = np.cumsum(sizes) - sizes
        slots = np.arange(len(keys)) - offsets[row]
        self.nodes[row, slots] = node
        self.gradients[row, slots] = values
        self.basis = np.asarray(basis)[:, :2, :2]
        self.shape = (count * len(basis), node_count * 2)
        self.count, self.node_count = count, node_count
        self.nnz = len(keys)
        self.storage_bytes = self.nodes.nbytes + self.gradients.nbytes + self.basis.nbytes
        self.T = LinearOperator(self.shape[::-1], matvec=self.transpose, dtype=float)

    def __matmul__(self, velocity):
        u = np.asarray(velocity).reshape(self.node_count, 2)
        return np.einsum("kcd,nsd,nsc->nk", self.basis, self.gradients, u[self.nodes]).ravel()

    def transpose(self, stress):
        stress = np.asarray(stress).reshape(self.count, len(self.basis))
        force = np.einsum("kcd,nsd,nk->nsc", self.basis, self.gradients, stress)
        return np.column_stack(
            [
                np.bincount(self.nodes.ravel(), weights=force[:, :, c].ravel(), minlength=self.node_count)
                for c in range(2)
            ]
        ).ravel()

    def colors(self, inverse_mass):
        """Greedy conflict coloring, including seams and component-wise walls."""
        users = [set() for _ in range(self.node_count)]
        colors = np.zeros(self.count, dtype=int)
        for i, (nodes, gradients) in enumerate(zip(self.nodes, self.gradients)):
            supported = nodes[np.any(gradients != 0, axis=1) & np.any(inverse_mass[nodes] > 0, axis=1)]
            forbidden = set().union(*(users[j] for j in supported))
            color = 0
            while color in forbidden:
                color += 1
            colors[i] = color
            for j in supported:
                users[j].add(color)
        return [np.flatnonzero(colors == c) for c in range(int(colors.max(initial=-1)) + 1)]


def diagonal_coulomb(diagonal, q, mu):
    """Batched exact local solves for positive diagonal metrics.

    r_n=max(-q_n/D_n,0). Sticking is explicit; sliding reduces to the monotone
    scalar equation ||q_t/(D_t+alpha)||=mu*r_n. Safeguarded Newton uses the
    smallest/largest tangential eigenvalues for an analytic bracket (§5.3).
    """
    diagonal, q = np.broadcast_arrays(np.asarray(diagonal, dtype=float), np.asarray(q, dtype=float))
    mu = np.broadcast_to(mu, q.shape[:-1])
    if np.any(diagonal <= 0) or not np.isfinite(diagonal).all():
        raise ValueError("local Coulomb metric must be finite and positive")
    result = np.zeros_like(q)
    result[..., 0] = np.maximum(-q[..., 0] / diagonal[..., 0], 0)
    radius = mu * result[..., 0]
    tangent = -q[..., 1:] / diagonal[..., 1:]
    length = np.linalg.norm(tangent, axis=-1)
    stick = (radius > 0) & (length <= radius)
    result[..., 1:] = np.where(stick[..., None], tangent, 0)
    slip = (radius > 0) & ~stick
    if not np.any(slip):
        return result
    d, b, y = diagonal[..., 1:][slip], q[..., 1:][slip], radius[slip]
    # Divide by the yield radius, so the target is one even at small pressure.
    b = b / y[:, None]
    norm = np.linalg.norm(b, axis=1)
    lower = np.maximum(0, norm - d.max(axis=1))
    upper = np.maximum(lower, norm - d.min(axis=1))
    alpha = lower.copy()
    for _ in range(64):
        denominator = d + alpha[:, None]
        value = b / denominator
        length = np.linalg.norm(value, axis=1)
        error = length - 1
        done = np.abs(error) <= 2e-13
        if np.all(done):
            break
        lower = np.where((error > 0) & ~done, alpha, lower)
        upper = np.where((error < 0) & ~done, alpha, upper)
        derivative = -np.sum(value * value / denominator, axis=1) / length
        newton = alpha - error / derivative
        safeguarded = np.where((newton >= lower) & (newton <= upper), newton, (lower + upper) / 2)
        alpha = np.where(done, alpha, safeguarded)
    else:
        raise RuntimeError("granular scalar sliding solve did not converge")
    result[..., 1:][slip] = -y[:, None] * value
    return result


def solve_velocity(
    operator, trial, inverse_mass, gap, mu, dt, *, tolerance=1e-9, max_iterations=200, initial=None
):
    """Colored Gauss-Seidel with cached tangential metrics and exact coupling.

    The block-diagonal *preconditioner* drops normal/tangent coupling and adds
    proximal regularization. Neither changes B, gap, or the residual equations.
    A conservative norm bound and per-node Euclidean natural residual make the
    acceptance check at least as strict as the reference infinity-norm check.
    Nonconvergence raises; no dense fallback is hidden.
    """
    count, components = gap.shape
    trial = np.asarray(trial).reshape(-1, 2)
    inverse_mass = np.asarray(inverse_mass).reshape(-1, 2)
    mu = np.broadcast_to(mu, count)
    q = (operator @ trial.ravel()).reshape(gap.shape) + gap
    magnitude = max(float(np.max(-q[:, 0], initial=0)), float(np.max(np.abs(q[:, 1:]), initial=0)))
    blocks = np.einsum("kcd,nsd->nskc", operator.basis, operator.gradients)
    local = dt * np.einsum("nskc,nsc,nslc->nkl", blocks, inverse_mass[operator.nodes], blocks)
    # ||B M^-1 B.T||_infinity <= || |B| M^-1 |B.T| 1 ||_infinity.
    absolute_force = np.zeros_like(trial)
    np.add.at(absolute_force, operator.nodes, np.abs(blocks).sum(axis=2))
    bound = dt * np.einsum("nskc,nsc->nk", np.abs(blocks), (inverse_mass * absolute_force)[operator.nodes])
    norm = float(np.max(bound, initial=0))
    if norm == 0 and np.any(q[:, 0] < 0):
        raise RuntimeError("packing constraint has no free velocity support")
    stats = dict(
        residual=0.0,
        evaluations=0,
        sweeps=0,
        converged=True,
        method="velocity",
        operator_nnz=operator.nnz,
        operator_storage_bytes=operator.storage_bytes,
        operator_format="padded gradient blocks",
        assembled_delassus=False,
        residual_norm="conservative block Euclidean",
    )
    if magnitude == 0 or norm == 0:
        return np.zeros_like(q), stats
    scale = magnitude / norm
    local /= norm
    regularizer = 1e-6 * np.maximum(np.trace(local, axis1=1, axis2=2), 1e-12)
    eigenvalues, rotation = np.linalg.eigh(local[:, 1:, 1:])
    diagonal = np.column_stack((local[:, 0, 0], eigenvalues)) + regularizer[:, None]
    Q = np.zeros((count, components, components))
    Q[:, 0, 0], Q[:, 1:, 1:] = 1, rotation
    del blocks, local, bound, absolute_force, eigenvalues, rotation
    colors = operator.colors(inverse_mass)
    stats["colors"] = len(colors)
    stats["preconditioner_storage_bytes"] = diagonal.nbytes + Q.nbytes + sum(group.nbytes for group in colors)
    z = np.zeros_like(q) if initial is None else np.array(initial, dtype=float) / scale
    # Begin with an admissible warmstart; Newton polishing can leave the cone
    # by solve-sized roundoff, which the final natural residual also checks.
    z[:, 0] = np.maximum(z[:, 0], 0)
    length = np.linalg.norm(z[:, 1:], axis=1)
    factor = np.minimum(1, np.divide(mu * z[:, 0], length, out=np.ones(count), where=length > 0))
    z[:, 1:] *= factor[:, None]

    def velocity():
        return trial + dt * scale * inverse_mass * (operator.T @ z.ravel()).reshape(-1, 2)

    def residual(u):
        w = ((operator @ u.ravel()).reshape(q.shape) + gap) / magnitude
        argument = z - w
        projected = argument.copy()
        projected[:, 0] = np.maximum(argument[:, 0], 0)
        radius = mu * np.maximum(z[:, 0], 0)
        length = np.linalg.norm(argument[:, 1:], axis=1)
        factor = np.minimum(1, np.divide(radius, length, out=np.ones(count), where=length > 0))
        projected[:, 1:] *= factor[:, None]
        return float(np.max(np.linalg.norm(z - projected, axis=1), initial=0))

    def feasible(value):
        value = value.copy()
        value[:, 0] = np.maximum(value[:, 0], 0)
        length = np.linalg.norm(value[:, 1:], axis=1)
        factor = np.minimum(1, np.divide(mu * value[:, 0], length, out=np.ones(count), where=length > 0))
        value[:, 1:] *= factor[:, None]
        return value

    u = velocity()
    error = residual(u)
    sweeps = 0
    history = []
    accelerated = 0
    newton_steps = 0
    while error > tolerance and sweeps < max_iterations * 10:
        start = z.copy()
        for group in colors:
            nodes, gradients = operator.nodes[group], operator.gradients[group]
            w = (np.einsum("kcd,nsd,nsc->nk", operator.basis, gradients, u[nodes]) + gap[group]) / magnitude
            old = np.einsum("nki,nk->ni", Q[group], z[group])
            other = np.einsum("nki,nk->ni", Q[group], w) - diagonal[group] * old
            new = diagonal_coulomb(diagonal[group], other, mu[group])
            delta = np.einsum("nki,ni->nk", Q[group], new - old)
            z[group] = np.einsum("nki,ni->nk", Q[group], new)
            force = np.einsum("kcd,nsd,nk->nsc", operator.basis, gradients, delta)
            np.add.at(u, nodes, dt * scale * inverse_mass[nodes] * force)
        sweeps += 1
        u = velocity()  # Periodic refresh avoids accumulating incremental roundoff.
        error = residual(u)
        # Safeguarded limited-memory Anderson acceleration. The current Q1-Q1
        # stress space can converge very slowly under plain colored sweeps.
        # Only accept a feasible extrapolate if its true residual decreases.
        history.append((start.ravel(), (z - start).ravel()))
        history = history[-13:]
        if len(history) > 1 and error > tolerance:
            dx = np.column_stack([b[0] - a[0] for a, b in zip(history[:-1], history[1:])])
            df = np.column_stack([b[1] - a[1] for a, b in zip(history[:-1], history[1:])])
            coefficients = np.linalg.lstsq(df, history[-1][1], rcond=1e-10)[0]
            extrapolate = feasible((z.ravel() - (dx + df) @ coefficients).reshape(z.shape))
            saved = z
            z = extrapolate
            accelerated_u = velocity()
            accelerated_error = residual(accelerated_u)
            if np.isfinite(accelerated_error) and accelerated_error < error:
                u, error = accelerated_u, accelerated_error
                accelerated += 1
            else:
                z = saved
        if sweeps % 25 == 0 and tolerance < error < 1e-3:
            # Polish nearly converged but degenerate stick/slip transitions.
            # J*v uses only B and B.T; GMRES retains a bounded Krylov basis.
            w = ((operator @ u.ravel()).reshape(q.shape) + gap) / magnitude
            argument = z - w
            radius = mu * np.maximum(z[:, 0], 0)
            length = np.linalg.norm(argument[:, 1:], axis=1)
            sliding = length > radius
            unit = np.divide(
                argument[:, 1:],
                length[:, None],
                out=np.zeros_like(argument[:, 1:]),
                where=length[:, None] > 0,
            )
            ratio = np.divide(radius, length, out=np.zeros(count), where=length > 0)
            projected = argument.copy()
            projected[:, 0] = np.maximum(argument[:, 0], 0)
            projected[:, 1:] *= np.minimum(
                1, np.divide(radius, length, out=np.ones(count), where=length > 0)
            )[:, None]
            natural = z - projected

            def jacobian(value):
                value = value.reshape(q.shape)
                du = dt / norm * inverse_mass * (operator.T @ value.ravel()).reshape(-1, 2)
                dw = (operator @ du.ravel()).reshape(q.shape)
                da = value - dw
                out = dw.copy()
                out[:, 0] = np.where(argument[:, 0] > 0, dw[:, 0], value[:, 0])
                tangent = ratio[:, None] * (da[:, 1:] - unit * np.sum(unit * da[:, 1:], axis=1)[:, None])
                tangent += (mu * (z[:, 0] >= 0) * value[:, 0])[:, None] * unit
                out[sliding, 1:] = value[sliding, 1:] - tangent[sliding]
                return out.ravel()

            J = LinearOperator((z.size, z.size), matvec=jacobian, dtype=float)
            correction, _ = gmres(J, -natural.ravel(), rtol=1e-5, atol=0, restart=min(z.size, 64), maxiter=4)
            saved = z
            for fraction in (1.0, 0.5, 0.25, 0.125):
                z = saved + fraction * correction.reshape(z.shape)
                candidate_u = velocity()
                candidate_error = residual(candidate_u)
                if np.isfinite(candidate_error) and candidate_error < error:
                    u, error = candidate_u, candidate_error
                    newton_steps += 1
                    history.clear()
                    break
            else:
                z = saved
    if not np.isfinite(error) or error > tolerance:
        raise RuntimeError(f"granular velocity solve did not converge: residual={error:.3g}, sweeps={sweeps}")
    stats.update(residual=error, sweeps=sweeps, acceleration_steps=accelerated, newton_steps=newton_steps)
    return z * scale, stats
