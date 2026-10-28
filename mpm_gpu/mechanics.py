"""Resident compact-operator kernels for linear Q1 with nodal stresses.

Static incidence is built once on the CPU (topology.q1_topology). Particle
binning and evolving gradient values are computed on the device. Each stress
row stores nine integrated gradient 2-vectors; B and B.T use the same values.
"""

from dataclasses import dataclass

import numpy as np

from .primitives import bins
from .topology import q1_topology

MODULE = "mechanics.slang"
INVALID = np.uint32(0xFFFFFFFF)


def _indices(values):
    values = np.asarray(values, dtype=np.int64)
    return np.where(values < 0, INVALID, values).astype(np.uint32)


@dataclass
class DeviceGrid:
    nx: int
    ny: int
    periodic: bool
    dx: float
    dy: float
    topology: object
    row_nodes: object
    neighbors: object
    transpose_rows: object
    transpose_slots: object

    @property
    def rows(self):
        return len(self.topology.nodes)

    @property
    def nodes(self):
        return self.nx * self.ny

    @property
    def cells(self):
        return (self.nx - 1) * (self.ny - 1)


def device_grid(rt, nx, ny, dx, dy, *, periodic=False):
    T = q1_topology(nx, ny, periodic=periodic)
    return DeviceGrid(
        nx,
        ny,
        bool(periodic),
        float(dx),
        float(dy),
        T,
        rt.buffer(T.nodes.astype(np.uint32), label="row_nodes"),
        rt.buffer(_indices(T.neighbors), label="neighbors"),
        rt.buffer(_indices(T.transpose_rows), label="transpose_rows"),
        rt.buffer(_indices(T.transpose_slots), label="transpose_slots"),
    )


@dataclass
class Particles:
    count: int
    cell_x: object
    cell_y: object
    fraction_x: object
    fraction_y: object
    volume: object


def upload_particles(rt, cells, fractions, volume):
    """Initial upload of split coordinates (topology.split_positions) and grain volume."""
    cells = np.asarray(cells, dtype=np.int32)
    fractions = np.asarray(fractions, dtype=np.float32)
    return Particles(
        len(cells),
        rt.buffer(cells[:, 0].view(np.uint32)),
        rt.buffer(cells[:, 1].view(np.uint32)),
        rt.buffer(fractions[:, 0], dtype=np.float32),
        rt.buffer(fractions[:, 1], dtype=np.float32),
        rt.buffer(np.asarray(volume, dtype=np.float32), dtype=np.float32),
    )


@dataclass
class Operator:
    gradients: object
    keys: object
    local_x: object
    local_y: object
    status: object
    bins: object


def assemble(rt, encoder, grid, particles):
    """Bin particles by cell and gather each stress row's integrated gradients.

    Status bit 1: a particle lies outside the grid. Bit 4: a stencil node is
    missing from the static incidence, which would indicate a topology error.
    """
    n = particles.count
    keys = rt.buffer(count=n)
    local_x = rt.buffer(count=n, dtype=np.float32)
    local_y = rt.buffer(count=n, dtype=np.float32)
    status = rt.buffer(count=1)
    gradients = rt.buffer(count=grid.rows * 18, dtype=np.float32)
    rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=status)
    rt.dispatch(
        encoder,
        "particle_cells",
        n,
        module=MODULE,
        count=n,
        nx=grid.nx,
        ny=grid.ny,
        periodic=int(grid.periodic),
        cell_x=particles.cell_x,
        cell_y=particles.cell_y,
        fraction_x=particles.fraction_x,
        fraction_y=particles.fraction_y,
        keys=keys,
        local_x=local_x,
        local_y=local_y,
        status=status,
    )
    cell_bins = bins(rt, encoder, keys, n, grid.cells)
    rt.dispatch(
        encoder,
        "assemble_gradients",
        grid.rows,
        module=MODULE,
        row_count=grid.rows,
        nx=grid.nx,
        ny=grid.ny,
        periodic=int(grid.periodic),
        dx=grid.dx,
        dy=grid.dy,
        row_nodes=grid.row_nodes,
        neighbors=grid.neighbors,
        cell_offsets=cell_bins.offsets,
        cell_counts=cell_bins.counts,
        permutation=cell_bins.permutation,
        local_x=local_x,
        local_y=local_y,
        volume=particles.volume,
        gradients=gradients,
        status=status,
    )
    return Operator(gradients, keys, local_x, local_y, status, cell_bins)


def apply_b(rt, encoder, grid, operator, velocity, output):
    """Strain-rate coefficients w = B u, three intrinsic components per row."""
    rt.dispatch(
        encoder,
        "apply_b",
        grid.rows,
        module=MODULE,
        row_count=grid.rows,
        neighbors=grid.neighbors,
        gradients=operator.gradients,
        velocity=velocity,
        output=output,
    )


def apply_bt(rt, encoder, grid, operator, stress, output):
    """Nodal forces B.T z, two components per canonical node; aliases receive zero."""
    rt.dispatch(
        encoder,
        "apply_bt",
        grid.nodes,
        module=MODULE,
        node_count=grid.nodes,
        transpose_rows=grid.transpose_rows,
        transpose_slots=grid.transpose_slots,
        gradients=operator.gradients,
        stress=stress,
        output=output,
    )


@dataclass
class Preconditioner:
    metric: object
    bound: object
    norm: object
    diagonal: object
    rotation: object
    status: object
    retained: list


def preconditioner(rt, encoder, grid, operator, inverse_mass, dt, status=None):
    """Normalised local metrics, regularised diagonal and tangential rotation.

    inverse_mass holds two components per canonical node, zero where a
    component is fixed. The global norm bound is reduced on the device and
    read by the normalising dispatch; no host round trip. Status bit 8: zero
    norm, which the solver must treat as unsupported compression or a trivial solve.
    """
    from .primitives import reduce

    retained = []
    if status is None:
        status = rt.buffer(count=1)
        rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=status)
    absolute = rt.buffer(count=2 * grid.nodes, dtype=np.float32)
    metric = rt.buffer(count=6 * grid.rows, dtype=np.float32)
    bound = rt.buffer(count=grid.rows, dtype=np.float32)
    diagonal = rt.buffer(count=3 * grid.rows, dtype=np.float32)
    rotation = rt.buffer(count=2 * grid.rows, dtype=np.float32)
    rt.dispatch(
        encoder,
        "absolute_force",
        grid.nodes,
        module=MODULE,
        node_count=grid.nodes,
        transpose_rows=grid.transpose_rows,
        transpose_slots=grid.transpose_slots,
        gradients=operator.gradients,
        output=absolute,
    )
    rt.dispatch(
        encoder,
        "local_metric",
        grid.rows,
        module=MODULE,
        row_count=grid.rows,
        dt=float(dt),
        neighbors=grid.neighbors,
        gradients=operator.gradients,
        inverse_mass=inverse_mass,
        absolute=absolute,
        metric=metric,
        bound=bound,
    )
    norm, reduce_status, scratch = reduce(rt, encoder, bound, grid.rows, operation="max", status=status)
    retained += [absolute, *scratch]
    rt.dispatch(
        encoder,
        "metric_finish",
        grid.rows,
        module=MODULE,
        row_count=grid.rows,
        norm=norm,
        metric=metric,
        diagonal=diagonal,
        rotation=rotation,
        status=status,
    )
    return Preconditioner(metric, bound, norm, diagonal, rotation, status, retained)


def local_coulomb(rt, encoder, diagonal, q, mu, output, status, count):
    """Batched exact local solves (test entry point for the solver's inner kernel)."""
    rt.dispatch(
        encoder,
        "local_coulomb_test",
        count,
        module=MODULE,
        count=count,
        diagonal=diagonal,
        q=q,
        mu=mu,
        output=output,
        status=status,
    )


TRANSFERS = {"flip": 0, "pic": 1, "apic_gradient": 2}


def fixed_components(boundary_v, boundary_h, roughness):
    """Static per-node wall mask: bit 1 fixes x, bit 2 fixes y (grid.update_momentum)."""
    boundary_v, boundary_h = np.asarray(boundary_v, bool), np.asarray(boundary_h, bool)
    fixed_x = boundary_v | (bool(roughness) & boundary_h)
    fixed_y = boundary_h | (bool(roughness) & boundary_v)
    return (fixed_x * 1 + fixed_y * 2).astype(np.uint32)


@dataclass
class GridState:
    mass: object
    reference: object
    momentum: object
    trial: object
    inverse_mass: object


def p2g(rt, encoder, grid, operator, fixed, particle_mass, velocity, gradient, *, dt, m_tol, body, affine):
    """Nodal mass, raw momentum, masked momentum after dt*m*g, trial velocity, free inverse mass.

    Uses the particle bins and local coordinates from assemble(). velocity has
    two and gradient four (row-major 2x2 L) components per particle.
    """
    n = grid.nodes
    state = GridState(
        rt.buffer(count=n, dtype=np.float32),
        rt.buffer(count=2 * n, dtype=np.float32),
        rt.buffer(count=2 * n, dtype=np.float32),
        rt.buffer(count=2 * n, dtype=np.float32),
        rt.buffer(count=2 * n, dtype=np.float32),
    )
    rt.dispatch(
        encoder,
        "p2g_nodes",
        n,
        module=MODULE,
        node_count=n,
        nx=grid.nx,
        ny=grid.ny,
        periodic=int(grid.periodic),
        dx=grid.dx,
        dy=grid.dy,
        dt=float(dt),
        m_tol=float(m_tol),
        body_x=float(body[0]),
        body_y=float(body[1]),
        affine=int(bool(affine)),
        fixed=fixed,
        cell_offsets=operator.bins.offsets,
        cell_counts=operator.bins.counts,
        permutation=operator.bins.permutation,
        local_x=operator.local_x,
        local_y=operator.local_y,
        particle_mass=particle_mass,
        particle_velocity=velocity,
        particle_gradient=gradient,
        mass=state.mass,
        reference=state.reference,
        momentum=state.momentum,
        trial=state.trial,
        inverse_mass=state.inverse_mass,
    )
    return state


def g2p(
    rt,
    encoder,
    grid,
    operator,
    particles,
    state,
    nodal_velocity,
    velocity,
    gradient,
    *,
    dt,
    m_tol,
    transfer,
    pic_fraction,
):
    """Particle velocity, velocity gradient L and split-coordinate advection (in place)."""
    rt.dispatch(
        encoder,
        "g2p_particles",
        particles.count,
        module=MODULE,
        count=particles.count,
        nx=grid.nx,
        periodic=int(grid.periodic),
        dx=grid.dx,
        dy=grid.dy,
        dt=float(dt),
        m_tol=float(m_tol),
        transfer=TRANSFERS[transfer],
        pic_fraction=float(pic_fraction),
        keys=operator.keys,
        local_x=operator.local_x,
        local_y=operator.local_y,
        mass=state.mass,
        reference=state.reference,
        velocity=nodal_velocity,
        cell_x=particles.cell_x,
        cell_y=particles.cell_y,
        fraction_x=particles.fraction_x,
        fraction_y=particles.fraction_y,
        particle_velocity=velocity,
        particle_gradient=gradient,
    )


def deform(
    rt,
    encoder,
    count,
    gradient,
    reference_volume,
    particle_mass,
    deformation,
    jacobian,
    volume,
    density,
    status,
    *,
    dt,
):
    """F <- blockdiag(exp(L dt), 1) F; J <- J exp(tr L dt); V and rho.

    jacobian holds det(F) on entry. Status bit 64: nonpositive or nonfinite
    J or F, in which case that particle's state is left unchanged.
    """
    rt.dispatch(
        encoder,
        "deform_particles",
        count,
        module=MODULE,
        count=count,
        dt=float(dt),
        particle_gradient=gradient,
        reference_volume=reference_volume,
        particle_mass=particle_mass,
        deformation=deformation,
        jacobian=jacobian,
        volume=volume,
        density=density,
        status=status,
    )
