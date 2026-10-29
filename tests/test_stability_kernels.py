"""Closed-form planar singular values used by the elastic timestep estimate."""

import numpy as np
from numpy.testing import assert_allclose

from stability import _singular_values


def test_planar_singular_values_match_lapack():
    rng = np.random.default_rng(0)
    F = np.tile(np.eye(3), (500, 1, 1))
    F[:, :2, :2] = rng.normal(size=(500, 2, 2)) + 2 * np.eye(2)
    F[0, :2, :2] = 1.3 * np.eye(2)
    F[1, :2, :2] = [[1.0, 1e-9], [0.0, 1.0]]
    expected = np.sort(np.linalg.svd(F, compute_uv=False), axis=1)
    assert_allclose(np.sort(_singular_values(F), axis=1), expected, rtol=1e-9, atol=0)


def test_nonplanar_deformation_uses_lapack():
    F = np.eye(3)[None] + 0.1 * np.arange(9.0).reshape(1, 3, 3)
    assert_allclose(_singular_values(F), np.linalg.svd(F, compute_uv=False), rtol=0, atol=0)
