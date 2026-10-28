"""Mechanics hardware gate: device kernels against independent FP64 CPU references."""

import numpy as np

from . import mechanics
from .fixtures import column
from .topology import split_positions

# FP32 accumulation over at most a few hundred terms per entry.
OPERATOR_TOLERANCE = 2e-6
ADJOINT_TOLERANCE = 1e-6
METRIC_TOLERANCE = 2e-6
LOCAL_TOLERANCE = 1e-5
TRANSFER_TOLERANCE = 2e-6
POSITION_TOLERANCE = 1e-6  # in grid spacings
DEFORMATION_TOLERANCE = 2e-5
PENDING = ()


def _require(condition, message):
    if not condition:
        raise AssertionError(message)


def _upload(rt, state):
    P, G, b = state["P"], state["G"], state["batch"]
    grid = mechanics.device_grid(rt, P.G.nx, P.G.ny, G.dx, G.dy, periodic=P.B.cyclic_lr)
    cells, fractions = split_positions(b.x[:, :2], [P.G.x_m, P.G.y_m], [G.dx, G.dy])
    particles = mechanics.upload_particles(rt, cells, fractions, b.m / b.material.rho_s)
    return grid, particles


def operator_check(rt, spacing, particles_per_cell, periodic, *, edges=False):
    state = column(spacing, particles_per_cell, periodic=periodic, edges=edges)
    B, active = state["operator"], state["active"]
    grid, particles = _upload(rt, state)
    rng = np.random.default_rng(int(1000 * spacing) + periodic)
    u = rng.normal(size=2 * grid.nodes)
    z = rng.normal(size=3 * grid.rows)
    velocity, stress = rt.buffer(u, dtype=np.float32), rt.buffer(z, dtype=np.float32)
    strain = rt.buffer(count=3 * grid.rows, dtype=np.float32)
    force = rt.buffer(count=2 * grid.nodes, dtype=np.float32)
    encoder = rt.device.create_command_encoder()
    operator = mechanics.assemble(rt, encoder, grid, particles)
    mechanics.apply_b(rt, encoder, grid, operator, velocity, strain)
    mechanics.apply_bt(rt, encoder, grid, operator, stress, force)
    rt.submit(encoder)
    _require(int(rt.read(operator.status, np.uint32)[0]) == 0, "unexpected assembly status")
    # Dense (stress node, velocity node) blocks: independent of slot ordering.
    T, N = grid.topology, grid.nodes
    expected = np.zeros((N, N, 2))
    for row, node in enumerate(active):
        np.add.at(expected[node], B.nodes[row], B.gradients[row])
    gradients = rt.read(operator.gradients, np.float32).reshape(grid.rows, 9, 2)
    actual = np.zeros_like(expected)
    for row, node in enumerate(T.nodes):
        valid = T.neighbors[row] >= 0
        actual[node, T.neighbors[row][valid]] += gradients[row][valid]
    scale = np.abs(expected).max()
    block_error = float(np.abs(actual - expected).max() / scale)
    _require(block_error <= OPERATOR_TOLERANCE, f"gradient blocks differ by {block_error:.3g}")
    device_row = np.full(N, -1)
    device_row[T.nodes] = np.arange(grid.rows)
    _require(np.all(device_row[active] >= 0), "active stress node without a device row")
    strain_device = rt.read(strain, np.float32).reshape(-1, 3).astype(np.float64)
    strain_reference = (B @ u).reshape(-1, 3)
    strain_error = float(
        np.abs(strain_device[device_row[active]] - strain_reference).max() / np.abs(strain_reference).max()
    )
    inactive = np.setdiff1d(np.arange(grid.rows), device_row[active])
    _require(np.all(strain_device[inactive] == 0), "inactive stress rows have nonzero strain rates")
    force_reference = B.transpose(z.reshape(-1, 3)[device_row[active]].ravel())
    force_device = rt.read(force, np.float32).astype(np.float64)
    force_error = float(np.abs(force_device - force_reference).max() / np.abs(force_reference).max())
    pairing = abs(strain_device.ravel() @ z - u @ force_device) / (np.abs(strain_device.ravel()) @ np.abs(z))
    _require(strain_error <= OPERATOR_TOLERANCE, f"B u differs by {strain_error:.3g}")
    _require(force_error <= OPERATOR_TOLERANCE, f"B.T z differs by {force_error:.3g}")
    _require(pairing <= ADJOINT_TOLERANCE, f"adjoint pairing error {pairing:.3g}")
    return dict(
        particles=particles.count,
        rows=grid.rows,
        active_rows=len(active),
        periodic=periodic,
        edges=edges,
        block_error=block_error,
        strain_error=strain_error,
        force_error=force_error,
        adjoint_error=float(pairing),
    )


def _inverse_mass(state):
    G, free = state["G"], state["free"]
    inverse = np.zeros(2 * len(G.m))
    inverse[free] = 1 / np.repeat(G.m, 2)[free]
    return inverse.reshape(-1, 2)


def _reference_metric(B, inverse_mass, dt):
    """FP64 copy of the preconditioner set-up in granular_mixed.solve_velocity."""
    blocks = np.einsum("kcd,nsd->nskc", B.basis, B.gradients)
    local = dt * np.einsum("nskc,nsc,nslc->nkl", blocks, inverse_mass[B.nodes], blocks)
    absolute_force = np.zeros_like(inverse_mass)
    np.add.at(absolute_force, B.nodes, np.abs(blocks).sum(axis=2))
    bound = dt * np.einsum("nskc,nsc->nk", np.abs(blocks), (inverse_mass * absolute_force)[B.nodes])
    norm = float(np.max(bound, initial=0))
    local /= norm
    regularizer = 1e-6 * np.maximum(np.trace(local, axis1=1, axis2=2), 1e-12)
    return local, norm, regularizer


def metric_check(rt, spacing, particles_per_cell, periodic):
    state = column(spacing, particles_per_cell, periodic=periodic)
    B, active, dt = state["operator"], state["active"], state["P"].dt
    grid, particles = _upload(rt, state)
    inverse_mass = _inverse_mass(state)
    encoder = rt.device.create_command_encoder()
    operator = mechanics.assemble(rt, encoder, grid, particles)
    pre = mechanics.preconditioner(
        rt, encoder, grid, operator, rt.buffer(inverse_mass, dtype=np.float32), dt, status=operator.status
    )
    rt.submit(encoder)
    _require(int(rt.read(pre.status, np.uint32)[0]) == 0, "unexpected preconditioner status")
    local, norm, regularizer = _reference_metric(B, inverse_mass, dt)
    T = grid.topology
    device_row = np.full(grid.nodes, -1)
    device_row[T.nodes] = np.arange(grid.rows)
    rows = device_row[active]
    norm_device = float(rt.read(pre.norm, np.float32)[0])
    norm_error = abs(norm_device - norm) / norm
    upper = np.triu_indices(3)
    metric = rt.read(pre.metric, np.float32).reshape(-1, 6)[rows].astype(np.float64)
    expected = local[:, upper[0], upper[1]]
    metric_error = float(np.abs(metric - expected).max() / np.abs(expected).max())
    diagonal = rt.read(pre.diagonal, np.float32).reshape(-1, 3)[rows].astype(np.float64)
    angle = rt.read(pre.rotation, np.float32).reshape(-1, 2)[rows].astype(np.float64)
    R = np.stack([np.stack([angle[:, 0], -angle[:, 1]], 1), np.stack([angle[:, 1], angle[:, 0]], 1)], 1)
    orthogonality = float(np.abs(np.einsum("nji,njk->nik", R, R) - np.eye(2)).max())
    tangential = local[:, 1:, 1:]
    reconstructed = np.einsum("nij,nj,nkj->nik", R, diagonal[:, 1:] - regularizer[:, None], R)
    scale = np.abs(local).max()
    eigen_error = float(np.abs(reconstructed - tangential).max() / scale)
    eigenvalues = np.sort(diagonal[:, 1:], axis=1)
    eigenvalue_error = float(
        np.abs(eigenvalues - (np.linalg.eigvalsh(tangential) + regularizer[:, None])).max() / scale
    )
    normal_error = float(np.abs(diagonal[:, 0] - (local[:, 0, 0] + regularizer)).max() / scale)
    for name, value in (
        ("norm", norm_error),
        ("metric", metric_error),
        ("rotation orthogonality", orthogonality),
        ("tangential eigensystem", eigen_error),
        ("eigenvalues", eigenvalue_error),
        ("normal diagonal", normal_error),
    ):
        _require(value <= METRIC_TOLERANCE, f"{name} differs by {value:.3g}")
    _require(np.all(diagonal > 0), "nonpositive regularised diagonal")
    return dict(
        rows=len(rows),
        norm_error=norm_error,
        metric_error=metric_error,
        orthogonality_error=orthogonality,
        eigensystem_error=eigen_error,
        eigenvalue_error=eigenvalue_error,
        normal_diagonal_error=normal_error,
    )


def _local_inputs():
    rng = np.random.default_rng(36)
    count = 4000
    diagonal = 10 ** rng.uniform(-3, 3, (count, 3))
    q = rng.normal(size=(count, 3)) * 10 ** rng.uniform(-8, 8, (count, 1))
    mu = rng.uniform(0, 1, count)
    mu[::7] = 0
    # Explicit branches: opening, zero friction, interior stick, cone boundary,
    # zero tangential trial, pure sliding and an anisotropic sliding metric.
    special_d = np.array([[1, 1, 1]] * 6 + [[1e-3, 1e3, 1e-3]], dtype=float)
    special_q = np.array(
        [[1, 1, 0], [-1, 1, 0], [-1, 0.2, 0], [-1, 0.5, 0], [-1, 0, 0], [-1, 1, -1], [-1, 1, 1]], dtype=float
    )
    special_mu = np.array([0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.5])
    diagonal = np.vstack([special_d, diagonal]).astype(np.float32)
    q = np.vstack([special_q, q]).astype(np.float32)
    mu = np.concatenate([special_mu, mu]).astype(np.float32)
    return diagonal, q, mu


def local_coulomb_check(rt):
    from granular_mixed import diagonal_coulomb

    diagonal, q, mu = _local_inputs()
    count = len(q)
    output = rt.buffer(count=3 * count, dtype=np.float32)
    status = rt.buffer(count=count)
    encoder = rt.device.create_command_encoder()
    mechanics.local_coulomb(
        rt,
        encoder,
        rt.buffer(diagonal, dtype=np.float32),
        rt.buffer(q, dtype=np.float32),
        rt.buffer(mu, dtype=np.float32),
        output,
        status,
        count,
    )
    rt.submit(encoder)
    flags = rt.read(status, np.uint32)
    _require(np.all(flags == 0), f"local solve status set on {int(np.sum(flags != 0))} rows")
    D, Q, M = diagonal.astype(np.float64), q.astype(np.float64), mu.astype(np.float64)
    r = rt.read(output, np.float32).reshape(-1, 3).astype(np.float64)
    expected = diagonal_coulomb(D, Q, M)
    # Errors relative to the unconstrained response magnitude, which bounds r.
    scale = np.linalg.norm(Q / D, axis=1)
    error = float(np.max(np.linalg.norm(r - expected, axis=1) / scale))
    w = D * r + Q
    complementarity = float(np.max(np.abs(r[:, 0] * w[:, 0]) / (scale * np.linalg.norm(Q, axis=1))))
    cone = float(np.max((np.linalg.norm(r[:, 1:], axis=1) - M * r[:, 0]) / scale))
    _require(np.all(r[:, 0] >= 0), "negative normal reaction")
    _require(np.all(r[M == 0, 1:] == 0), "tangential reaction without friction")
    _require(np.all(r[expected[:, 0] == 0] == 0), "opening rows must have zero reaction")
    for name, value in (("reaction", error), ("complementarity", complementarity), ("cone", cone)):
        _require(value <= LOCAL_TOLERANCE, f"local {name} error {value:.3g}")
    branches = dict(
        opening=int(np.sum(expected[:, 0] == 0)),
        frictionless=int(np.sum((M == 0) & (expected[:, 0] > 0))),
        sticking=int(np.sum(np.linalg.norm(-Q[:, 1:] / D[:, 1:], axis=1) <= M * expected[:, 0])),
    )
    branches["sliding"] = count - sum(branches.values())
    return dict(
        rows=count,
        reaction_error=error,
        complementarity_error=complementarity,
        cone_error=cone,
        branches=branches,
    )


def _cpu_grid_step(state):
    """The grid part of granular.step before the contact solve, in FP64."""
    from transfers import nodal_velocity, scatter

    P, G, b = state["P"], state["G"], state["batch"]
    body = np.array([P.g * np.sin(P.theta), P.g * np.cos(P.theta), 0.0])
    scatter(G.fe, b.n, b.N[:, :, None] * b.m[:, None, None] * body)
    if P.B.cyclic_lr:
        G.make_cyclic(P, G, ["fe"])
    G.BCs(P)
    G.update_momentum(P)
    return body, nodal_velocity(G, P)[:, :2]


def transfer_check(rt, spacing, particles_per_cell, periodic, transfer, pic_fraction=0.0):
    from copy import deepcopy
    from transfers import nodal_velocity

    from .topology import join_positions

    state = column(
        spacing,
        particles_per_cell,
        periodic=periodic,
        transfer=transfer,
        pic_fraction=pic_fraction,
        motion=True,
    )
    P, G, b = state["P"], state["G"], state["batch"]
    body, trial_reference = _cpu_grid_step(state)
    grid, particles = _upload(rt, state)
    fixed = rt.buffer(mechanics.fixed_components(G.boundary_v, G.boundary_h, P.B.roughness))
    masses = rt.buffer(b.m, dtype=np.float32)
    velocity = rt.buffer(b.v[:, :2], dtype=np.float32)
    gradient = rt.buffer(b.L[:, :2, :2], dtype=np.float32)
    # Arbitrary end-of-solve nodal velocity on resolved nodes, aliases copied.
    final = np.random.default_rng(3).normal(scale=0.2, size=(grid.nodes, 2))
    final[~G.valid_nodes] = 0
    canonical = grid.topology.canonical
    final = final[canonical]
    encoder = rt.device.create_command_encoder()
    operator = mechanics.assemble(rt, encoder, grid, particles)
    nodal = mechanics.p2g(
        rt,
        encoder,
        grid,
        operator,
        fixed,
        masses,
        velocity,
        gradient,
        dt=P.dt,
        m_tol=P.M_tol,
        body=body,
        affine=transfer == "apic_gradient",
    )
    final_device = rt.buffer(final, dtype=np.float32)
    mechanics.g2p(
        rt,
        encoder,
        grid,
        operator,
        particles,
        nodal,
        final_device,
        velocity,
        gradient,
        dt=P.dt,
        m_tol=P.M_tol,
        transfer=transfer,
        pic_fraction=pic_fraction,
    )
    rt.submit(encoder)
    _require(int(rt.read(operator.status, np.uint32)[0]) == 0, "unexpected assembly status")

    def relative(actual, expected):
        return float(np.abs(actual - expected).max() / np.abs(expected).max())

    inverse_reference = _inverse_mass(state)
    inverse_device = rt.read(nodal.inverse_mass, np.float32).reshape(-1, 2)
    _require(np.array_equal(inverse_device > 0, inverse_reference > 0), "free-component pattern differs")
    errors = dict(
        mass_error=relative(rt.read(nodal.mass, np.float32), G.m),
        reference_momentum_error=relative(
            rt.read(nodal.reference, np.float32).reshape(-1, 2), G.q_reference[:, :2]
        ),
        momentum_error=relative(rt.read(nodal.momentum, np.float32).reshape(-1, 2), G.q[:, :2]),
        trial_velocity_error=relative(rt.read(nodal.trial, np.float32).reshape(-1, 2), trial_reference),
        inverse_mass_error=relative(inverse_device, inverse_reference),
    )
    # CPU grid-to-particle with the same final nodal velocity.
    G.q = np.column_stack((final * G.m[:, None], np.zeros(grid.nodes)))
    L_reference = np.einsum("pri,prj->pij", nodal_velocity(G, P)[b.n], b.G)[:, :2, :2]
    moved = deepcopy(b)
    moved.move_material_points(P, G)
    x_reference = moved.x[:, :2].copy()
    if P.B.cyclic_lr:
        x_reference[:, 0] = P.G.x_m + (x_reference[:, 0] - P.G.x_m) % (P.G.x_M - P.G.x_m)
    cells = np.column_stack(
        (
            rt.read(particles.cell_x, np.uint32).view(np.int32),
            rt.read(particles.cell_y, np.uint32).view(np.int32),
        )
    )
    fractions = np.column_stack(
        (rt.read(particles.fraction_x, np.float32), rt.read(particles.fraction_y, np.float32))
    )
    _require(np.all((fractions >= 0) & (fractions < 1)), "invalid split fraction after advection")
    x_device = join_positions(cells, fractions, [P.G.x_m, P.G.y_m], [G.dx, G.dy])
    if P.B.cyclic_lr:
        x_device[:, 0] = P.G.x_m + (x_device[:, 0] - P.G.x_m) % (P.G.x_M - P.G.x_m)
    errors.update(
        particle_velocity_error=relative(rt.read(velocity, np.float32).reshape(-1, 2), moved.v[:, :2]),
        velocity_gradient_error=relative(rt.read(gradient, np.float32).reshape(-1, 2, 2), L_reference),
        position_error=float(np.abs(x_device - x_reference).max() / min(G.dx, G.dy)),
    )
    for name, value in errors.items():
        limit = POSITION_TOLERANCE if name == "position_error" else TRANSFER_TOLERANCE
        _require(value <= limit, f"{name} {value:.3g}")
    return dict(particles=particles.count, transfer=transfer, pic_fraction=pic_fraction, **errors)


def _deformation_inputs():
    rng = np.random.default_rng(11)
    cases = {
        "tiny": rng.normal(scale=1e-6, size=(64, 2, 2)),
        "moderate": rng.normal(scale=0.3, size=(256, 2, 2)),
        "rotation": np.array([[0.0, -1.0], [1.0, 0.0]]) * rng.uniform(1, 20, (64, 1, 1)),
        "shear": np.array([[0.0, 1.0], [0.0, 0.0]]) * rng.uniform(-10, 10, (64, 1, 1)),
        "extension": np.array([[1.0, 0.0], [0.0, -1.0]]) * rng.uniform(-5, 5, (64, 1, 1)),
        "large": rng.normal(scale=4.0, size=(64, 2, 2)),
    }
    return cases


def deformation_check(rt):
    from scipy.linalg import expm

    results = {}
    worst = 0.0
    for name, A in _deformation_inputs().items():
        count = len(A)
        rng = np.random.default_rng(count)
        F = np.tile(np.eye(3), (count, 1, 1))
        F[:, :2, :2] += rng.normal(scale=0.05, size=(count, 2, 2))
        F[:, 2, 2] = rng.uniform(0.9, 1.1, count)
        A32, F32 = A.astype(np.float32), F.astype(np.float32)
        mass = rng.uniform(0.5, 2.0, count).astype(np.float32)
        V0 = rng.uniform(0.5, 2.0, count).astype(np.float32)
        buffers = dict(
            deformation=rt.buffer(F32, dtype=np.float32),
            jacobian=rt.buffer(np.linalg.det(F32.astype(np.float64)), dtype=np.float32),
            volume=rt.buffer(count=count, dtype=np.float32),
            density=rt.buffer(count=count, dtype=np.float32),
        )
        status = rt.buffer(np.zeros(1, dtype=np.uint32))
        encoder = rt.device.create_command_encoder()
        mechanics.deform(
            rt,
            encoder,
            count,
            rt.buffer(A32, dtype=np.float32),
            rt.buffer(V0, dtype=np.float32),
            rt.buffer(mass, dtype=np.float32),
            buffers["deformation"],
            buffers["jacobian"],
            buffers["volume"],
            buffers["density"],
            status,
            dt=1.0,
        )
        rt.submit(encoder)
        _require(int(rt.read(status, np.uint32)[0]) == 0, f"{name}: unexpected deformation status")
        increment = np.tile(np.eye(3), (count, 1, 1))
        increment[:, :2, :2] = [expm(a) for a in A32.astype(np.float64)]
        expected = increment @ F32.astype(np.float64)
        J = np.linalg.det(expected)
        actual = rt.read(buffers["deformation"], np.float32).reshape(-1, 3, 3)
        scale = np.linalg.norm(expected, axis=(1, 2))
        F_error = float(np.max(np.linalg.norm(actual - expected, axis=(1, 2)) / scale))
        J_error = float(np.max(np.abs(rt.read(buffers["jacobian"], np.float32) - J) / J))
        V_error = float(np.max(np.abs(rt.read(buffers["volume"], np.float32) - V0 * J) / (V0 * J)))
        rho_error = float(
            np.max(np.abs(rt.read(buffers["density"], np.float32) - mass / (V0 * J)) / (mass / (V0 * J)))
        )
        _require(
            np.all(actual[:, 2, :2] == 0) and np.all(actual[:, :2, 2] == 0), f"{name}: out-of-plane coupling"
        )
        _require(np.array_equal(actual[:, 2, 2], F32[:, 2, 2]), f"{name}: F_zz changed")
        case = dict(F_error=F_error, J_error=J_error, V_error=V_error, rho_error=rho_error)
        for label, value in case.items():
            _require(value <= DEFORMATION_TOLERANCE, f"{name}: {label} {value:.3g}")
        worst = max(worst, *case.values())
        results[name] = case
    # Nonfinite gradients and expansions beyond FP32 range must be reported.
    bad = np.array([[[np.nan, 0], [0, 0]], [[200.0, 0], [0, 200.0]]], dtype=np.float32)
    status = rt.buffer(np.zeros(1, dtype=np.uint32))
    deformation = rt.buffer(np.tile(np.eye(3, dtype=np.float32), (2, 1, 1)), dtype=np.float32)
    encoder = rt.device.create_command_encoder()
    ones = rt.buffer(np.ones(2, dtype=np.float32), dtype=np.float32)
    scratch = [rt.buffer(np.ones(2), dtype=np.float32) for _ in range(3)]
    mechanics.deform(
        rt, encoder, 2, rt.buffer(bad, dtype=np.float32), ones, ones, deformation, *scratch, status, dt=1.0
    )
    rt.submit(encoder)
    _require(int(rt.read(status, np.uint32)[0]) & 64, "invalid deformation not reported")
    _require(
        np.array_equal(rt.read(deformation, np.float32).reshape(-1, 3, 3), np.tile(np.eye(3), (2, 1, 1))),
        "invalid F was written",
    )
    return dict(cases=results, worst_relative_error=worst)


def invalid_particle_check(rt):
    state = column(0.25, 2)
    grid, particles = _upload(rt, state)
    cells = np.array([[0, 0], [grid.nx - 1, 0], [-1, 0], [0, grid.ny]], dtype=np.int32)
    fractions = np.array([[0.5, 0.5], [0.25, 0.5], [0.5, 0.5], [0.5, 0.0]], dtype=np.float32)
    particles = mechanics.upload_particles(rt, cells, fractions, np.ones(len(cells)))
    encoder = rt.device.create_command_encoder()
    operator = mechanics.assemble(rt, encoder, grid, particles)
    rt.submit(encoder)
    _require(int(rt.read(operator.status, np.uint32)[0]) & 1, "out-of-grid particle not reported")
    keys = rt.read(operator.keys, np.uint32)
    _require(keys[0] == 0 and np.all(keys[1:] == grid.cells), "invalid particles must have an invalid key")
    return dict(particles=len(cells))


def cases(rt):
    checks = [
        (
            f"operator_h{h}_ppc{n}_{'periodic' if periodic else 'walls'}",
            lambda h=h, n=n, p=periodic: operator_check(rt, h, n, p),
        )
        for h, n, periodic in (
            (0.25, 2, False),
            (0.25, 2, True),
            (0.1, 3, False),
            (0.1, 3, True),
            (0.5, 4, True),
            (1.0, 6, True),
        )
    ]
    checks += [
        ("operator_edges_walls", lambda: operator_check(rt, 0.25, 2, False, edges=True)),
        ("operator_edges_periodic", lambda: operator_check(rt, 0.25, 2, True, edges=True)),
        ("invalid_particles", lambda: invalid_particle_check(rt)),
        ("metric_h0.25_walls", lambda: metric_check(rt, 0.25, 2, False)),
        ("metric_h0.1_periodic", lambda: metric_check(rt, 0.1, 3, True)),
        ("metric_h1.0_periodic", lambda: metric_check(rt, 1.0, 6, True)),
        ("local_coulomb", lambda: local_coulomb_check(rt)),
        ("transfer_flip_walls", lambda: transfer_check(rt, 0.1, 3, False, "flip")),
        ("transfer_flip_blend_periodic", lambda: transfer_check(rt, 0.1, 3, True, "flip", 0.3)),
        ("transfer_pic_walls", lambda: transfer_check(rt, 0.25, 2, False, "pic")),
        ("transfer_apic_gradient_walls", lambda: transfer_check(rt, 0.1, 3, False, "apic_gradient")),
        ("transfer_apic_gradient_periodic", lambda: transfer_check(rt, 0.5, 4, True, "apic_gradient")),
        ("deformation", lambda: deformation_check(rt)),
    ]
    return checks
