"""Static Q1 incidence and translation-resistant particle coordinates.

Only grid geometry is constructed on the CPU. Particle binning and evolving
operator values belong on the GPU. Canonical periodic aliases are merged here,
including grids narrower than the nine-node support.
"""

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Topology:
    nodes: np.ndarray
    neighbors: np.ndarray
    transpose_rows: np.ndarray
    transpose_slots: np.ndarray
    colors: np.ndarray
    canonical: np.ndarray


def q1_topology(nx, ny, *, periodic=False):
    if int(nx) != nx or int(ny) != ny or nx < 2 or ny < 2:
        raise ValueError("Q1 grid dimensions must be integers >= 2")
    nx, ny = int(nx), int(ny)
    width = nx - 1 if periodic else nx
    canonical = np.arange(nx * ny, dtype=np.uint32).reshape(ny, nx)
    if periodic:
        canonical[:, -1] = canonical[:, 0]
    nodes = canonical[:, :width].ravel().copy()
    row_of = {int(node): row for row, node in enumerate(nodes)}
    neighbors = np.full((len(nodes), 9), -1, dtype=np.int32)
    transpose_rows = np.full((nx * ny, 9), -1, dtype=np.int32)
    transpose_slots = np.full((nx * ny, 9), -1, dtype=np.int32)
    incidence = [[] for _ in range(nx * ny)]
    colors = np.zeros(len(nodes), dtype=np.uint32)
    users = [set() for _ in range(nx * ny)]
    for node, row in row_of.items():
        x, y = node % nx, node // nx
        support = set()
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                xx, yy = x + dx, y + dy
                if periodic:
                    xx %= width
                if 0 <= xx < width and 0 <= yy < ny:
                    support.add(yy * nx + xx)
        support = sorted(support)
        forbidden = set().union(*(users[j] for j in support))
        color = 0
        while color in forbidden:
            color += 1
        colors[row] = color
        for slot, j in enumerate(support):
            neighbors[row, slot] = j
            incidence[j].append((row, slot))
            users[j].add(color)
    for j, values in enumerate(incidence):
        for i, (row, slot) in enumerate(values):
            transpose_rows[j, i], transpose_slots[j, i] = row, slot
    return Topology(nodes, neighbors, transpose_rows, transpose_slots, colors, canonical.ravel())


def split_positions(positions, origin, spacing):
    positions = np.asarray(positions, dtype=np.float64)
    origin, spacing = np.asarray(origin, dtype=np.float64), np.asarray(spacing, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 2 or origin.shape != (2,) or spacing.shape != (2,):
        raise ValueError("positions must be (N,2), with two-component origin and spacing")
    if (
        not np.isfinite(positions).all()
        or not np.isfinite(origin).all()
        or not np.isfinite(spacing).all()
        or np.any(spacing <= 0)
    ):
        raise ValueError("positions/origin/spacing must be finite, with positive spacing")
    local = (positions - origin) / spacing
    cells = np.floor(local)
    fraction = (local - cells).astype(np.float32)
    # Rounding near a cell edge can make a fraction exactly one.
    carry = fraction >= 1
    cells += carry
    fraction[carry] = 0
    if np.any(cells < np.iinfo(np.int32).min) or np.any(cells > np.iinfo(np.int32).max):
        raise ValueError("particle cell coordinates exceed int32")
    return cells.astype(np.int32), fraction


def join_positions(cells, fraction, origin, spacing):
    return np.asarray(origin, dtype=np.float64) + (
        np.asarray(cells, dtype=np.float64) + np.asarray(fraction, dtype=np.float64)
    ) * np.asarray(spacing, dtype=np.float64)
