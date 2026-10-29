"""Resident compact-operator kernels for linear Q1 with nodal stresses.

Static incidence is built once on the CPU (topology.q1_topology). Particle
binning and evolving gradient values are computed on the device. Each stress
row stores nine integrated gradient 2-vectors; B and B.T use the same values.
"""

from dataclasses import dataclass

import numpy as np

from .primitives import Binner, Reduction
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


class Operator:
    """Workspace and result of assemble(): bins, local coordinates, gradients, occupancy."""

    def __init__(self, rt, grid, count):
        self.keys = rt.buffer(count=count)
        self.local_x = rt.buffer(count=count, dtype=np.float32)
        self.local_y = rt.buffer(count=count, dtype=np.float32)
        self.status = rt.buffer(count=1)
        self.gradients = rt.buffer(count=grid.rows * 18, dtype=np.float32)
        self.occupied = rt.buffer(count=grid.rows, dtype=np.float32)
        self.bins = Binner(rt, count, grid.cells)


def assemble(rt, encoder, grid, particles, workspace=None):
    """Bin particles by cell and gather each stress row's integrated gradients.

    Also gathers each row's occupied grain volume sum_p V_p N_row(x_p). Status
    bit 1: a particle lies outside the grid. Bit 4: a stencil node is missing
    from the static incidence, which would indicate a topology error.
    """
    n = particles.count
    ws = Operator(rt, grid, n) if workspace is None else workspace
    rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=ws.status)
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
        keys=ws.keys,
        local_x=ws.local_x,
        local_y=ws.local_y,
        status=ws.status,
    )
    cell_bins = ws.bins.encode(encoder, ws.keys)
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
        local_x=ws.local_x,
        local_y=ws.local_y,
        volume=particles.volume,
        gradients=ws.gradients,
        occupied=ws.occupied,
        status=ws.status,
    )
    return ws


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


class Preconditioner:
    """Workspace and result of preconditioner(); the global norm is norm[0]."""

    def __init__(self, rt, grid):
        self.absolute = rt.buffer(count=2 * grid.nodes, dtype=np.float32)
        self.metric = rt.buffer(count=6 * grid.rows, dtype=np.float32)
        self.bound = rt.buffer(count=grid.rows, dtype=np.float32)
        self.diagonal = rt.buffer(count=3 * grid.rows, dtype=np.float32)
        self.rotation = rt.buffer(count=2 * grid.rows, dtype=np.float32)
        self.reduction = Reduction(rt, grid.rows, operation="max")
        self.norm = self.reduction.result
        self.status = None


def preconditioner(rt, encoder, grid, operator, inverse_mass, dt, status=None, workspace=None):
    """Normalised local metrics, regularised diagonal and tangential rotation.

    inverse_mass holds two components per canonical node, zero where a
    component is fixed. The global norm bound is reduced on the device and
    read by the normalising dispatch; no host round trip. Status bit 8: zero
    norm, which the solver must treat as unsupported compression or a trivial solve.
    """
    ws = Preconditioner(rt, grid) if workspace is None else workspace
    if status is None:
        status = rt.buffer(count=1) if ws.status is None else ws.status
        rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=status)
    ws.status = status
    rt.dispatch(
        encoder,
        "absolute_force",
        grid.nodes,
        module=MODULE,
        node_count=grid.nodes,
        transpose_rows=grid.transpose_rows,
        transpose_slots=grid.transpose_slots,
        gradients=operator.gradients,
        output=ws.absolute,
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
        absolute=ws.absolute,
        metric=ws.metric,
        bound=ws.bound,
    )
    ws.reduction.encode(encoder, ws.bound, status)
    rt.dispatch(
        encoder,
        "metric_finish",
        grid.rows,
        module=MODULE,
        row_count=grid.rows,
        norm=ws.norm,
        metric=ws.metric,
        diagonal=ws.diagonal,
        rotation=ws.rotation,
        status=status,
    )
    return ws


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


def fixed_components(boundary_v, boundary_h, nx, separable=False):
    """Static per-node wall mask for grid.wall_contact.

    Bit 1 marks a wall normal to x (boundary_v), bit 2 a wall normal to y
    (boundary_h). With separable, the domain edges are one-sided: bits 4/8 free
    an x wall while the node moves in +x/-x (left/right edge), bits 16/32 a y
    wall while it moves in +y/-y (bottom/top). Interior walls stay two-sided.
    Roughness is applied on the device from the walls in contact.
    """
    boundary_v, boundary_h = np.asarray(boundary_v, bool), np.asarray(boundary_h, bool)
    mask = boundary_v * 1 + boundary_h * 2
    if separable:
        row, column = np.divmod(np.arange(len(mask)), nx)
        ny = len(mask) // nx
        mask += boundary_v * ((column == 0) * 4 + (column == nx - 1) * 8)
        mask += boundary_h * ((row == 0) * 16 + (row == ny - 1) * 32)
    return mask.astype(np.uint32)


class GridState:
    def __init__(self, rt, grid):
        n = grid.nodes
        self.mass = rt.buffer(count=n, dtype=np.float32)
        self.reference = rt.buffer(count=2 * n, dtype=np.float32)
        self.momentum = rt.buffer(count=2 * n, dtype=np.float32)
        self.trial = rt.buffer(count=2 * n, dtype=np.float32)
        self.inverse_mass = rt.buffer(count=2 * n, dtype=np.float32)
        self.contact = rt.buffer(count=n, dtype=np.uint32)


def p2g(
    rt,
    encoder,
    grid,
    operator,
    fixed,
    particle_mass,
    velocity,
    gradient,
    *,
    dt,
    m_tol,
    body,
    affine,
    rough,
    state=None,
):
    """Nodal mass, raw momentum, masked momentum after dt*m*g, trial velocity, free inverse mass.

    state.contact receives the walls each node is in contact with (bit 1 x
    wall, bit 2 y wall); rough fixes both components at a wall in contact.

    Uses the particle bins and local coordinates from assemble(). velocity has
    two and gradient four (row-major 2x2 L) components per particle.
    """
    n = grid.nodes
    state = GridState(rt, grid) if state is None else state
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
        rough=int(bool(rough)),
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
        contact=state.contact,
    )
    return state


def wall_recheck(rt, encoder, grid, fixed, state, velocity, changed, *, rough):
    """Clamp released wall nodes that the solved velocity drives into their wall.

    Adds them to state.contact, zeroes their trial velocity and inverse mass and
    sets changed[0] when any node was added (granular.step).
    """
    rt.dispatch(
        encoder,
        "wall_recheck",
        grid.nodes,
        module=MODULE,
        node_count=grid.nodes,
        rough=int(bool(rough)),
        fixed=fixed,
        velocity=velocity,
        contact=state.contact,
        trial=state.trial,
        inverse_mass=state.inverse_mass,
        changed=changed,
    )


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
    shift=None,
):
    """Particle velocity, velocity gradient L and split-coordinate advection (in place).

    shift receives each particle's periodic wrap in cells (unwrapped minus
    wrapped cell), for budgets that need the position before wrapping.
    """
    if shift is None:
        shift = rt.buffer(count=particles.count, dtype=np.float32)
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
        shift=shift,
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
