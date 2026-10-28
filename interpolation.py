"""Complete tensor-product particle/grid stencils on a rectangular grid.

Quadratic cardinal B-splines use linear ghost-value extrapolation at physical
edges. In the outer half-cell this reduces to linear weights, remains C1 at the
join, and preserves nonnegative weights, partition, and linear completeness.
Periodic x support retains unwrapped offsets but deposits to unique seam nodes.
"""

import numpy as np

SIZES = {"linear": 4, "quadratic": 9}


def _quadratic_axis(position, lo, h, count, periodic=False):
    coordinate = (position - lo) / h
    raw = np.floor(coordinate - 0.5).astype(np.intp)[:, None] + np.arange(3)
    distance = coordinate[:, None] - raw
    absolute = np.abs(distance)
    weight = np.where(absolute < 0.5, 0.75 - distance**2, 0.5 * (1.5 - absolute) ** 2)
    derivative = np.where(absolute < 0.5, -2 * distance, -np.sign(distance) * (1.5 - absolute)) / h
    if periodic:
        return raw % (count - 1), weight, derivative, raw * h + lo
    left = coordinate < 0.5
    right = coordinate > count - 1.5
    weight[left] = np.column_stack((np.zeros(left.sum()), 1 - coordinate[left], coordinate[left]))
    derivative[left] = [0, -1 / h, 1 / h]
    raw[left] = [0, 0, 1]
    remaining = count - 1 - coordinate[right]
    weight[right] = np.column_stack((remaining, 1 - remaining, np.zeros(right.sum())))
    derivative[right] = [-1 / h, 1 / h, 0]
    raw[right] = [count - 2, count - 1, count - 1]
    raw = np.clip(raw, 0, count - 1)  # only zero-weight outer entries at exact joins
    return raw, weight, derivative, raw * h + lo


def stencil(position, P, G):
    """Return nodes, weights, spatial gradients, and node-minus-particle offsets."""
    if P.interpolation == "linear":
        cell = np.floor((position[:, :2] - [P.G.x_m, P.G.y_m]) / [G.dx, G.dy]).astype(np.intp)
        cell = np.minimum(cell, [P.G.nx - 2, P.G.ny - 2])
        base = cell[:, 0] + P.G.nx * cell[:, 1]
        nodes = base[:, None] + [0, 1, P.G.nx + 1, P.G.nx]
        x = (position[:, 0] - G.X[base]) / G.dx
        y = (position[:, 1] - G.Y[base]) / G.dy
        weight = np.stack(((1 - x) * (1 - y), x * (1 - y), x * y, (1 - x) * y), axis=1)
        gradient = np.zeros((*nodes.shape, 3))
        gradient[:, :, 0] = np.stack((-(1 - y), 1 - y, y, -y), axis=1) / G.dx
        gradient[:, :, 1] = np.stack((-(1 - x), -x, x, 1 - x), axis=1) / G.dy
        coordinates = np.stack((G.X[nodes], G.Y[nodes], np.zeros_like(weight)), axis=2)
    else:
        ix, wx, gx, xx = _quadratic_axis(position[:, 0], P.G.x_m, G.dx, P.G.nx, P.B.cyclic_lr)
        iy, wy, gy, yy = _quadratic_axis(position[:, 1], P.G.y_m, G.dy, P.G.ny)
        nodes = (ix[:, None, :] + P.G.nx * iy[:, :, None]).reshape(-1, 9)
        weight = (wy[:, :, None] * wx[:, None, :]).reshape(-1, 9)
        gradient = np.zeros((*nodes.shape, 3))
        gradient[:, :, 0] = (wy[:, :, None] * gx[:, None, :]).reshape(-1, 9)
        gradient[:, :, 1] = (gy[:, :, None] * wx[:, None, :]).reshape(-1, 9)
        coordinates = np.zeros((*nodes.shape, 3))
        coordinates[:, :, 0] = np.broadcast_to(xx[:, None, :], (len(position), 3, 3)).reshape(-1, 9)
        coordinates[:, :, 1] = np.broadcast_to(yy[:, :, None], (len(position), 3, 3)).reshape(-1, 9)
    return nodes, weight, gradient, coordinates - position[:, None, :]


def moment(weights, offsets):
    return np.einsum("nr,nri,nrj->nij", weights, offsets, offsets)


def affine_matrix(affine_momentum, weights, offsets):
    """B D^+ in the resolved 2D support; no inverse-mass or arbitrary ridge."""
    D = moment(weights, offsets)
    return affine_momentum @ np.linalg.pinv(D, rcond=1e-12, hermitian=True)
