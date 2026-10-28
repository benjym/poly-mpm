"""Regression checks for gradient masks and state-scoped computation reuse."""

from copy import deepcopy
from types import SimpleNamespace as NS

import numpy as np
import pytest
from numpy.testing import assert_allclose

from conftest import material, state
from diagnostics import inventory
from interpolation import moment
from MPM import run, step
import stability


def masked_gradient(P, G, values, smooth=False):
    values = np.ma.masked_where(G.m <= P.M_tol, values).reshape(P.G.ny, P.G.nx)
    dy, dx = np.gradient(values, G.dy, G.dx)
    if smooth:
        from astropy.convolution import Gaussian2DKernel, convolve

        kernel = Gaussian2DKernel(x_stddev=1, y_stddev=1)
        dy = convolve(dy.filled(np.nan), kernel, boundary="extend")
        dx = convolve(dx.filled(np.nan), kernel, boundary="extend")
    return np.column_stack(
        (
            np.nan_to_num(np.ma.filled(dx, 0)).ravel(),
            np.nan_to_num(np.ma.filled(dy, 0)).ravel(),
            np.zeros_like(G.m),
        )
    )


@pytest.mark.parametrize("shape", [(2, 2), (3, 7), (8, 5)])
@pytest.mark.parametrize("smooth", [False, True])
def test_gradient_matches_masked_reference(shape, smooth):
    ny, nx = shape
    P, G, L = state(G=NS(nx=nx, ny=ny, x_m=-1.0, x_M=2.0, y_m=0.0, y_M=2.0))
    rng = np.random.default_rng(10)
    for fraction in (0.0, 0.3, 1.0):
        G.m[:] = rng.random(len(G.m)) >= fraction
        values = rng.normal(size=len(G.m))
        with np.errstate(invalid="ignore", divide="ignore"):
            assert_allclose(
                G.calculate_gradient(P, G, values, smooth),
                masked_gradient(P, G, values, smooth),
                rtol=0,
                atol=0,
            )


def test_gradient_nonfinite_values_match_masked_reference():
    P, G, L = state()
    G.m[:] = 1
    values = np.arange(len(G.m), dtype=float)
    values[[0, 2, 13]] = [np.nan, np.inf, -np.inf]
    with np.errstate(invalid="ignore"):
        assert_allclose(G.calculate_gradient(P, G, values), masked_gradient(P, G, values), rtol=0, atol=0)


@pytest.mark.parametrize("periodic", [False, True])
def test_apic_reuses_inverse_and_invalidates_changed_geometry(periodic, monkeypatch):
    P, G, L = state(interpolation="quadratic", transfer="apic", B=NS(cyclic_lr=periodic))
    b = L.materials[0]
    b.get_basis_functions(P, G)
    original = np.linalg.pinv
    calls = []

    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(np.linalg, "pinv", counted)
    D, inverse = b.apic_matrices()
    assert_allclose(D, moment(b.N, b.offsets), rtol=0, atol=0)
    for scale in (1.0, 2.0):
        b.affine_momentum[:] = np.diag([0.2, -0.1, 0]) @ D * scale
        actual = inventory(P, L)["affine_kinetic"]
        expected = 0.5 * np.sum(
            b.m * np.einsum("nij,njk,nik->n", b.affine_momentum, inverse, b.affine_momentum)
        )
        assert_allclose(actual, expected, rtol=0, atol=0)
        assert b.apic_matrices()[1] is inverse
    assert len(calls) == 1
    b.x[:, 0] += 0.02
    b.get_basis_functions(P, G)
    assert b.apic_matrices()[1] is not inverse
    assert len(calls) == 2
    b.keep(np.arange(len(b)) % 2 == 0)
    assert b.apic_matrices()[0].shape == (len(b), 3, 3)
    assert len(calls) == 3
    b.append(deepcopy(b))
    assert b.apic_matrices()[0].shape == (len(b), 3, 3)
    assert len(calls) == 4


def test_dynamic_limits_share_state_work_but_reevaluate_material_callbacks(monkeypatch):
    callback_dts = []

    def critical(P):
        callback_dts.append(P.dt)
        return 1e-3 + 0.1 * P.dt

    P, G, L = state(time_stepping="dynamic", materials=material(critical_time=critical))
    original = stability._state_limits
    evaluations = []

    def counted(*args):
        evaluations.append(1)
        return original(*args)

    monkeypatch.setattr(stability, "_state_limits", counted)
    run(P, G, L, max_steps=3, output=False)
    # Initial state plus the state following each of three steps. Material
    # hooks still see each trial timestep, just as in the uncached controller.
    assert len(evaluations) == 4
    assert len(callback_dts) == 12
    cached_dts = callback_dts.copy()
    cached_history = deepcopy(L.diagnostics.history)
    cached_arrays = {name: getattr(L.materials[0], name).copy() for name in L.materials[0].persistent_fields}
    uncached = stability.limits

    def reference_limits(P, G, L, *, cache=None):
        return uncached(P, G, L)

    monkeypatch.setattr(stability, "limits", reference_limits)
    callback_dts.clear()
    reference = state(time_stepping="dynamic", materials=material(critical_time=critical))
    run(*reference, max_steps=3, output=False)
    assert callback_dts == cached_dts
    assert reference[2].diagnostics.history == cached_history
    for name, expected in cached_arrays.items():
        assert_allclose(getattr(reference[2].materials[0], name), expected, rtol=0, atol=0)
    before = len(evaluations)
    step(P, G, L)
    assert len(evaluations) == before + 1


def test_limits_are_fresh_after_output_callbacks_and_between_runs(monkeypatch):
    P, G, L = state(time_stepping="dynamic", interpolation="quadratic", savetime=1e-6)
    seen = []
    original = stability._state_limits

    def counted(P, G, L):
        seen.append(L.materials[0].F.copy())
        return original(P, G, L)

    def change_state(L, P, G):
        b = L.materials[0]
        b.F *= 1.01
        b.J = np.linalg.det(b.F)

    P.O.after_every_savetime = change_state
    monkeypatch.setattr(stability, "_state_limits", counted)
    run(P, G, L, max_steps=1)
    assert_allclose(seen[-1], L.materials[0].F, rtol=0, atol=0)
    L.materials[0].F *= 1.01
    changed = L.materials[0].F.copy()
    run(P, G, L, max_steps=1, output=False)
    assert_allclose(seen[2], changed, rtol=0, atol=0)
    assert not np.array_equal(seen[1], seen[2])
