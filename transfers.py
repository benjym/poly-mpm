"""Repeated-index-safe particle/grid transfers."""
import numpy as np


def scatter(target, nodes, contributions):
    """Add (particles, neighbours, ...) contributions to a nodal array."""
    index = nodes.ravel()
    if target.ndim == 1:
        target += np.bincount(index, weights=contributions.ravel(), minlength=len(target))
    else:
        for component in range(target.shape[1]):
            target[:, component] += np.bincount(
                index, weights=contributions[..., component].ravel(), minlength=len(target))


def divide(numerator, denominator):
    shape = np.broadcast_shapes(np.shape(numerator), np.shape(denominator))
    return np.divide(numerator, denominator, out=np.zeros(shape, dtype=float), where=np.asarray(denominator) != 0)


def nodal_velocity(grid, params, field='q'):
    return np.divide(getattr(grid, field), grid.m[:, None], out=np.zeros_like(grid.q),
                     where=(grid.m > params.M_tol)[:, None])


def gather(values, nodes, weights, valid):
    w = weights * valid[nodes]
    return np.sum(w[(...,) + (None,) * (values.ndim - 1)] * values[nodes], axis=1)
