"""Analytical constitutive checks, independent of the archived scalar formulas."""

from copy import deepcopy
import numpy as np
import pytest
from numpy.testing import assert_allclose
from conftest import material, state
from constit import EYE, RETIRED, update_stress
from MPM import run
from settings import MATERIAL_FIELDS
from initialise import get_parameters

NEW = ("bingham", "linear_mu", "dp", "dp_rate")


def sample(name, n=4, **overrides):
    data = dict(
        mu_s=0.6, beta=0.2, t_star=0.2, yield_stress=2.0, eta_0=100.0, eta=10.0, eta_v=5.0, b=0.4, eta_max=1e6
    )
    if name == "linear_mu":
        data["mu_s"] = 0.3
    data = {key: value for key, value in data.items() if key in MATERIAL_FIELDS[name]}
    data.update(overrides)
    P, G, L = state(materials=material(name, n=n, **data), dt=0.1)
    return P, G, L, L.materials[0]


def apply(b, P, strain):
    b.dstrain[:] = strain
    b.law(b, b.material, P)
    b.stress += b.dstress
    b.check_finite(P)


def shear(gamma):
    result = np.zeros((3, 3))
    result[0, 1] = result[1, 0] = gamma / 2
    return result


@pytest.mark.parametrize("name", NEW)
@pytest.mark.parametrize("n", [0, 1, 4])
def test_empty_zero_and_multiple_particles(name, n):
    P, G, L, b = sample(name, n=n)
    update_stress(b, b.material, P)
    assert_allclose(b.stress, 0)
    run(P, G, L, max_steps=3, output=False)
    b.check_finite(P)


def test_bingham_shear_curve_transition_and_reversal():
    P, _, _, b = sample("bingham", n=7)
    gamma_c = b.material.yield_stress / (b.material.eta_0 - b.material.eta)
    rates = np.array([0, gamma_c * 0.5, gamma_c, gamma_c * (1 + 1e-8), 1.0, -1.0, -gamma_c])
    apply(b, P, np.array([shear(rate * P.dt) for rate in rates]))
    stress = np.sign(rates) * np.minimum(
        b.material.eta_0 * abs(rates), b.material.yield_stress + b.material.eta * abs(rates)
    )
    assert_allclose(b.stress[:, 0, 1], stress)
    assert_allclose(b.gammadot, abs(rates))
    assert np.all(np.sum(b.stress * b.dstrain, axis=(1, 2)) >= 0)
    assert_allclose(np.trace(b.dev_stress, axis1=1, axis2=2), 0)


def test_bingham_newtonian_limit_and_bulk_dissipation():
    P, _, _, b = sample("bingham", yield_stress=0.0, K=0.0)
    rate = np.array([[2.0, 1.0, 0], [1.0, -1.0, 0], [0, 0, 2.0]])
    trace = np.trace(rate)
    apply(b, P, rate * P.dt)
    expected = 2 * b.material.eta * (rate - trace * EYE / 3) + b.material.eta_v * trace * EYE
    assert_allclose(b.stress, np.broadcast_to(expected, b.stress.shape))
    assert np.all(np.sum(b.stress * rate, axis=(1, 2)) >= 0)
    # Stopping removes viscous stress rather than retaining a stale increment.
    apply(b, P, np.zeros((3, 3)))
    assert_allclose(b.stress, 0)


def test_bingham_bulk_modulus_and_viscous_pressure_are_separate():
    P, _, _, b = sample("bingham", eta_v=7.0)
    strain = -0.01 * EYE
    apply(b, P, strain)
    assert_allclose(b.pressure, b.material.K * np.trace(strain))
    assert_allclose(b.p, -b.pressure - 7 * np.trace(strain) / P.dt)
    elastic_pressure = b.pressure.copy()
    apply(b, P, np.zeros((3, 3)))
    assert_allclose(b.p, -elastic_pressure)


def test_linear_mu_friction_curve_composition_and_zero_limits():
    P, _, _, b = sample("linear_mu", n=4)
    # Give each particle its own composition-weighted grain diameter.
    P.G.s = np.array([0.001, 0.003])
    b.phi = np.array([[1, 0], [0, 1], [0.25, 0.75], [0.5, 0.5]])
    b.pressure[:] = [-10, -10, 0, -10]
    rates = np.array([2.0, -2.0, 2.0, 0.0])
    apply(b, P, np.array([shear(rate * P.dt) for rate in rates]))
    expected_I = np.array(
        [2 * 0.001 * np.sqrt(b.material.rho_s / 10), 2 * 0.003 * np.sqrt(b.material.rho_s / 10), 0, 0]
    )
    assert_allclose(b.I, expected_I)
    assert_allclose(b.mu, 0.3 + 0.4 * expected_I)
    assert_allclose(b.stress[:2, 0, 1], [10 * b.mu[0], -10 * b.mu[1]])
    assert_allclose(b.stress[2], 0)
    assert_allclose(b.stress[3], -10 * EYE)
    assert np.all(np.sum(b.dev_stress * b.dstrain, axis=(1, 2)) >= 0)


def test_linear_mu_constant_friction_cap_and_no_tension():
    P, _, _, b = sample("linear_mu", b=0.0, eta_max=1.0)
    b.pressure[:] = [-10, -10, -10, -1]
    strains = np.array([shear(1e-6), shear(0.1), shear(-0.1), 0.1 * EYE])
    apply(b, P, strains)
    assert np.all(b.eta <= 1)
    assert_allclose(b.mu, 0.3)
    assert_allclose(b.stress[:3, 0, 1], [1e-5, 1.0, -1.0])
    assert_allclose(b.stress[3], 0)
    assert b.dp[3] == 1


@pytest.mark.parametrize("name", ["dp", "dp_rate"])
@pytest.mark.parametrize("beta", [0.0, 0.2, 0.6])
def test_dp_yield_flow_and_nonnegative_plastic_work(name, beta):
    P, _, _, b = sample(name, beta=beta)
    initial = -10 * EYE
    b.stress[:] = initial
    strain = shear(1.0)
    apply(b, P, strain)
    delta = b.stress - initial
    tr = np.trace(delta, axis1=1, axis2=2)
    elastic = (delta - tr[:, None, None] * EYE / 3) / (2 * b.material.G) + tr[:, None, None] * EYE / (
        9 * b.material.K
    )
    plastic = strain - elastic
    plastic_trace = np.trace(plastic, axis1=1, axis2=2)
    plastic_dev = plastic - plastic_trace[:, None, None] * EYE / 3
    multiplier = np.sqrt(2 / 3) * np.linalg.norm(plastic_dev, axis=(1, 2))
    assert np.all(multiplier > 0)
    assert_allclose(plastic_trace, beta * multiplier, atol=1e-15)
    assert_allclose(
        plastic_dev, multiplier[:, None, None] * 1.5 * b.dev_stress / b.q[:, None, None], atol=1e-15
    )
    f = b.q - b.material.mu_s * b.p
    viscous = 3 * b.material.G * b.material.t_star if name == "dp_rate" else 0
    assert_allclose(f, viscous * multiplier / P.dt, atol=2e-14)
    assert_allclose(b.yieldfunction, f, atol=2e-14)
    assert_allclose(b.mu, b.material.mu_s)
    assert_allclose(b.dp, -b.p + 10.0)
    work = np.sum(b.stress * plastic, axis=(1, 2))
    assert_allclose(work, (b.q - beta * b.p) * multiplier, atol=1e-13)
    assert np.all(work >= -1e-13)


@pytest.mark.parametrize("name", ["dp", "dp_rate"])
def test_dp_elastic_hydrostatic_response_and_unloading(name):
    P, _, _, b = sample(name)
    b.stress[:] = -10 * EYE
    apply(b, P, -0.01 * EYE)
    assert_allclose(b.p, 13.0)
    assert_allclose(b.q, 0)
    b.stress[:] = -10 * EYE + shear(2 * 6 / np.sqrt(3))
    old = b.stress.copy()
    increment = shear(-0.001)
    apply(b, P, increment)
    assert_allclose(b.stress, old + 2 * b.material.G * increment)
    assert np.all(b.yieldfunction < 0)


@pytest.mark.parametrize("name", ["dp", "dp_rate"])
@pytest.mark.parametrize("beta", [0.0, 0.6])
def test_dp_apex_tension_cutoff_and_zero_confinement(name, beta):
    P, _, _, b = sample(name, beta=beta)
    b.stress[:] = [10 * EYE, np.zeros((3, 3)), -EYE, 10 * EYE + shear(2)]
    apply(b, P, np.zeros((3, 3)))
    assert np.all(b.p >= -1e-14)
    assert_allclose(b.stress[0], 0, atol=1e-14)
    assert_allclose(b.stress[1], 0)
    assert_allclose(b.stress[2], -EYE)
    if name == "dp":
        assert_allclose(b.stress[3], 0, atol=1e-14)
        assert np.all(b.yieldfunction <= 1e-12)
    else:
        # Tension is released instantly; a finite viscosity may carry shear at p=0.
        assert b.q[3] > 0
        assert b.p[3] == 0
        assert b.q[3] < np.sqrt(3)


def test_dp_rate_zero_relaxation_time_is_rate_independent():
    P, _, _, b = sample("dp")
    other_P, _, _, other = sample("dp_rate", t_star=0.0)
    initial = np.array([-10 * EYE, 10 * EYE, -EYE + shear(2), np.zeros((3, 3))])
    b.stress[:], other.stress[:] = initial, initial
    apply(b, P, shear(0.1))
    apply(other, other_P, shear(0.1))
    assert_allclose(b.stress, other.stress)


def test_dp_rate_relaxation_and_timestep_convergence():
    errors = []
    for steps in (20, 40, 80):
        P, _, _, b = sample("dp_rate", n=1, beta=0.0)
        duration, rate = 0.4, 0.1
        P.dt = duration / steps
        b.stress[:] = -100 * EYE + shear(2 * b.material.mu_s * 100 / np.sqrt(3))
        for _ in range(steps):
            apply(b, P, shear(rate * P.dt))
        exact = (
            np.sqrt(3) * b.material.G * rate * b.material.t_star * (1 - np.exp(-duration / b.material.t_star))
        )
        errors.append(abs(b.yieldfunction[0] - exact))
        before = b.q.copy()
        apply(b, P, np.zeros((3, 3)))
        assert_allclose(b.yieldfunction, (before - b.material.mu_s * 100) / (1 + P.dt / b.material.t_star))
    assert errors[1] < errors[0] * 0.55
    assert errors[2] < errors[1] * 0.55


@pytest.mark.parametrize("name", NEW)
def test_tensor_rotation_covariance(name):
    P, _, _, b = sample(name)
    rng = np.random.default_rng(42)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    other = deepcopy(b)
    strain = rng.normal(size=b.dstrain.shape) * 0.001
    strain = (strain + strain.swapaxes(1, 2)) / 2
    initial = np.broadcast_to(-10 * EYE + shear(0.2), b.stress.shape).copy()
    b.stress[:] = initial
    other.stress[:] = rotation @ initial @ rotation.T
    b.pressure[:], other.pressure[:] = -10, -10
    apply(b, P, strain)
    apply(other, P, rotation @ strain @ rotation.T)
    assert_allclose(other.stress, rotation @ b.stress @ rotation.T, atol=1e-13)


@pytest.mark.parametrize("name", ["dp", "dp_rate"])
def test_grid_update_rotates_solid_stress_in_correct_direction(name):
    P, G, L, b = sample(name)
    b.get_reference_nodes(P, G)
    b.get_basis_functions(P, G)
    spin = np.array([[0, -2.0, 0], [2.0, 0, 0], [0, 0, 0]])
    G.m[:] = 1
    G.q[:] = np.column_stack((G.X, G.Y, np.zeros_like(G.X))) @ spin.T
    G.get_valid_nodes(P)
    initial = -10 * EYE + shear(0.2)
    b.stress[:] = initial
    angle = 2 * P.dt
    R = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    b.update_stress_strain(P, G)
    assert_allclose(b.stress, np.broadcast_to(R @ initial @ R.T, b.stress.shape), atol=1e-13)
    assert_allclose(b.stress, b.dstress + initial)
    assert_allclose(b.pressure, -10)
    assert_allclose(b.q, np.sqrt(3) * 0.1)


@pytest.mark.parametrize("name", ["bingham", "linear_mu"])
def test_grid_fluid_stress_does_not_receive_legacy_lie_increment(name):
    P, G, L, b = sample(name)
    b.get_reference_nodes(P, G)
    b.get_basis_functions(P, G)
    gradient = np.array([[0, 2.0, 0], [0, 0, 0], [0, 0, 0]])
    G.m[:] = 1
    G.q[:] = np.column_stack((G.X, G.Y, np.zeros_like(G.X))) @ gradient.T
    G.get_valid_nodes(P)
    b.stress[:] = -10 * EYE
    b.pressure[:] = -10
    b.update_stress_strain(P, G)
    expected = -10 * EYE + shear(4 * b.eta[0])
    assert_allclose(b.stress, np.broadcast_to(expected, b.stress.shape), atol=1e-12)


@pytest.mark.parametrize(
    "name,overrides,match",
    [
        ("dp", {"mu_s": -1}, "mu_s"),
        ("dp", {"beta": -0.1}, "beta"),
        ("dp", {"beta": 1}, "beta"),
        ("dp", {"G": 0}, "G"),
        ("dp", {"K": np.nan}, "K"),
        ("dp_rate", {"t_star": -0.1}, "t_star"),
        ("bingham", {"eta": -1}, "eta"),
        ("bingham", {"eta_0": 10}, "eta_0"),
        ("bingham", {"yield_stress": -1}, "yield_stress"),
        ("bingham", {"eta_v": -1}, "eta_v"),
        ("bingham", {"gamma_c": 10}, "gamma_c"),
        ("linear_mu", {"b": -1}, "b"),
        ("linear_mu", {"eta_max": 0}, "eta_max"),
        ("linear_mu", {"rho_s": np.inf}, "rho_s"),
    ],
)
def test_invalid_physical_parameters(name, overrides, match):
    with pytest.raises(ValueError, match=match):
        sample(name, **overrides)


@pytest.mark.parametrize("name", RETIRED)
def test_retired_laws_explain_reason(name):
    with pytest.raises(ValueError, match="Retired:"):
        state(name)


@pytest.mark.parametrize("example", ["dp_unit_test", "dp_rate_unit_test"])
def test_dp_examples_write_array_histories(example, tmp_path):
    P, G, L = get_parameters([example])
    P.save_dir = str(tmp_path)
    run(P, G, L, max_steps=3)
    with np.load(tmp_path / "data/pq_history.npz") as history:
        assert history["p"].shape == (4, 9)
        assert_allclose(history["time"], [0, P.dt, 2 * P.dt, 3 * P.dt])
        assert_allclose(history["p"][-1], L.materials[0].p)
        assert_allclose(history["q"][-1], L.materials[0].q)
    assert (tmp_path / "pq_history.png").stat().st_size > 1000


@pytest.mark.parametrize("name", ["viscous", "pouliquen", "bingham", "linear_mu"])
@pytest.mark.parametrize("mode", ["compression", "lithostatic"])
def test_fluid_initialisation_needs_no_elastic_shear_modulus(name, mode):
    # These laws read no shear modulus G, so material() leaves it out.
    data = dict(yield_stress=2.0, eta_0=100.0, b=0.4, eta_max=1e6)
    data = {key: value for key, value in data.items() if key in MATERIAL_FIELDS[name]}
    # Keep all nonzero stress-gradient nodes mass-supported; this test targets
    # fluid initialisation, not a particle exactly on an empty cell edge.
    S = material(name, **data)
    assert not hasattr(S, "G")
    S.Y = S.Y + 0.01
    for K in (0.0, 100.0):
        S.K = K
        P, G, L = state(materials=deepcopy(S), initial_stress=mode, initial_pressure=10.0, gravity=-1.0)
        L.materials[0].check_finite(P)
        if mode == "compression":
            assert_allclose(L.materials[0].p, 10.0)
        run(P, G, L, max_steps=1, output=False)


def test_linear_mu_zero_friction_newtonian_limit():
    P, _, _, b = sample("linear_mu", mu_s=0.0, b=0.4)
    b.pressure.fill(-10.0)
    apply(b, P, np.zeros((3, 3)))
    expected = 0.4 * (b.phi @ P.G.s) * np.sqrt(b.material.rho_s * 10)
    assert_allclose(b.eta, expected)
    apply(b, P, shear(0.1))
    assert_allclose(b.stress[:, 0, 1], expected * 0.1 / P.dt)


@pytest.mark.parametrize("name", NEW)
def test_batched_history_is_independent_between_particles(name):
    P, G, L, batch = sample(name)
    rng = np.random.default_rng(13)
    batch.stress[:] = [-10 * EYE + shear(2), 10 * EYE, -EYE, np.zeros((3, 3))]
    batch.pressure[:] = [-10.0, 10.0, -1.0, 0.0]
    strain = rng.normal(size=batch.dstrain.shape) * 0.2
    strain = (strain + strain.swapaxes(1, 2)) / 2
    rows = []
    for i in range(len(batch)):
        row = deepcopy(batch)
        row.keep(np.arange(len(batch)) == i)
        apply(row, P, strain[i])
        rows.append(row)
    apply(batch, P, strain)
    for name in ("stress", "pressure", "gammadot", "eta", "I", "mu", "yieldfunction"):
        assert_allclose(getattr(batch, name), np.concatenate([getattr(row, name) for row in rows]))


def test_linear_mu_near_zero_confinement_is_finite_and_has_vanishing_strength():
    P, _, _, b = sample("linear_mu")
    b.pressure[:] = [-1e-310, -1e-200, -1e-100, 0]
    apply(b, P, shear(0.1))
    assert np.isfinite(b.I).all()
    assert np.all(abs(b.stress[:, 0, 1]) < 1e-49)
