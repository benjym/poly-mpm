"""Mechanics hardware gate: device kernels against independent FP64 CPU references."""

import numpy as np

from . import mechanics
from .fixtures import column
from .topology import split_positions

# FP32 accumulation over at most a few hundred terms per entry.
OPERATOR_TOLERANCE = 2e-6
ADJOINT_TOLERANCE = 1e-6
PENDING = (
    "local metric and tangential eigensystem",
    "local Coulomb solves",
    "planar exponential deformation",
)


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
    ]
    return checks
