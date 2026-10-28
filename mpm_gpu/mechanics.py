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
