"""Scientific acceptance checks for the opt-in rigid-grain prototype."""

from copy import deepcopy
from types import SimpleNamespace as NS

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import sparse

from granular import assemble, basis_volumes, cone_projection, solve_coulomb
from initialise import get_parameters, prepare
from MPM import run, step
from stability import limits


def column(**options):
    from inputs.granular_column import Params

    P = Params([])
    for key, value in options.items():
        setattr(P, key, value)
    return prepare(P)


@pytest.mark.parametrize(
    "q, expected",
    [
        ([-1, 0.2, 0, 0], [1, -0.2, 0, 0]),
        ([-1, 1, 0, 0], [1, -0.5, 0, 0]),
        ([1, 1, 0, 0], [0, 0, 0, 0]),
        ([-1, -1, 0, 0], [1, 0.5, 0, 0]),
    ],
)
def test_coulomb_static_sliding_and_opening_without_associated_dilation(q, expected):
    reaction, stats = solve_coulomb(sparse.eye(4), q, 0.5)
    assert_allclose(reaction.ravel(), expected, atol=1e-11)
    velocity = reaction.ravel() + q
    assert velocity[0] >= -1e-11
    assert abs(reaction[0, 0] * velocity[0]) < 1e-11
    assert stats["converged"]


def test_coulomb_coupled_spd_problem():
    rng = np.random.default_rng(4)
    A = rng.normal(size=(12, 12))
    W = A @ A.T + np.eye(12)
    q = rng.normal(size=(3, 4))
    q[:, 0] -= 1
    r, stats = solve_coulomb(W, q, 0.4)
    w = (W @ r.ravel()).reshape(-1, 4) + q
    assert np.min(r[:, 0]) > -1e-8
    assert np.min(w[:, 0]) > -1e-8
    assert_allclose(r[:, 0] * w[:, 0], 0, atol=1e-8)
    assert np.all(np.linalg.norm(r[:, 1:], axis=1) <= 0.4 * r[:, 0] + 1e-8)
    for i in range(3):
        if np.linalg.norm(w[i, 1:]) > 1e-7:
            assert_allclose(r[i, 1:], -0.4 * r[i, 0] * w[i, 1:] / np.linalg.norm(w[i, 1:]), atol=1e-8)


def test_large_opening_gap_does_not_hide_compression_at_another_node():
    q = [[-1e-10, 0.2e-10, 0, 0], [1e10, 0, 0, 0]]
    reaction, _ = solve_coulomb(sparse.eye(8), q, 0.5)
    assert_allclose(reaction[0], [1e-10, -0.2e-10, 0, 0], rtol=1e-10, atol=1e-22)
    assert_allclose(reaction[1], 0, atol=1e-22)


def test_cone_projection_derivative_matches_finite_difference():
    x = np.array([[2.0, 0.1, 0.2, 0.3], [-2.0, 0.1, 0.2, 0.3], [0.2, 1.0, 0.3, -0.2]])
    _, jac = cone_projection(x, 0.4)
    for k in range(4):
        delta = np.zeros_like(x)
        delta[:, k] = 1e-6
        derivative = (cone_projection(x + delta, 0.4)[0] - cone_projection(x - delta, 0.4)[0]) / 2e-6
        assert_allclose(derivative, jac[:, :, k], atol=2e-10)


@pytest.mark.parametrize("basis", ["linear"])
@pytest.mark.parametrize("periodic", [False, True])
def test_integrated_basis_volume_and_adjoint_operator(basis, periodic):
    P, G, L = column(interpolation=basis)
    P.B.cyclic_lr = periodic
    if periodic:
        G.boundary_v[:] = False
    volume = basis_volumes(P, G)
    assert_allclose(volume.sum(), 1.5, atol=1e-14)
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    B, free, _, _, _ = assemble(P, G, L.materials[0])
    rng = np.random.default_rng(5)
    u = rng.normal(size=free.sum())
    stress = rng.normal(size=B.shape[0])
    assert_allclose(stress @ (B[:, free] @ u), u @ (B[:, free].T @ stress), atol=1e-15)
    for component in (0, 1):
        translation = np.zeros((len(G.m), 2))
        translation[:, component] = 1
        assert_allclose(B @ translation.ravel(), 0, atol=1e-15)


def test_static_column_at_large_timestep_without_viscosity_creep():
    P, G, L = column(dt=0.02, t_f=0.2)
    initial = L.materials[0].x.copy()
    run(P, G, L, output=False)
    b = L.materials[0]
    assert np.max(np.abs(b.v)) < 1e-7
    assert_allclose(b.x, initial, atol=1e-8)
    assert b.pressure.min() < -100
    assert b.pressure.max() <= 1e-6
    assert b.yieldfunction.max() < 1e-6
    assert G.granular_stats["residual"] < P.granular_tolerance
    assert G.granular_stats["contact_power"] <= 1e-6
    assert abs(L.diagnostics.last["mass_error"]) < 1e-10
    assert np.linalg.norm(L.diagnostics.last["momentum_residual"]) < 1e-9
    assert "elastic_wave" not in P.timestep_limits and "viscosity" not in P.timestep_limits
    assert L.diagnostics.last["granular_solver"]["converged"]


def test_low_packing_free_fall_and_recompaction_are_stress_free():
    P, G, L = column(dt=0.005, t_f=0.01, max_g=-2.0)
    P.B.has_bottom = False
    G.boundary_h[:] = False
    P.B.has_left = P.B.has_right = False
    G.boundary_v[:] = False
    b = L.materials[0]
    b.m *= 0.2
    b.rho *= 0.2
    b.initial_rho *= 0.2
    # Affine contraction is permitted while packing is below capacity.
    b.v[:, 0] = -0.01 * (b.x[:, 0] - 0.5)
    before = b.v[:, 1].copy()
    step(P, G, L)
    assert_allclose(b.pressure, 0, atol=1e-8)
    assert_allclose(b.v[:, 1] - before, -2 * P.dt, atol=1e-12)
    assert np.any(b.L[:, 0, 0] < 0)
    assert G.granular_stats["maximum_packing"] < 0.2


def test_rigid_grain_response_independent_of_elastic_and_viscosity_parameters(monkeypatch):
    import granular

    problems = []
    original = granular.solve_coulomb

    def record(W, q, mu, **options):
        problems.append((W.toarray(), q.copy(), mu))
        return original(W, q, mu, **options)

    monkeypatch.setattr(granular, "solve_coulomb", record)
    P, G, L = column()
    reference = deepcopy((P, G, L))
    P.S[0].K = 1e18
    P.S[0].G = 1e17
    P.S[0].eta_max = 1e20
    assert limits(P, G, L) == limits(*reference)
    step(P, G, L)
    step(*reference)
    for actual, expected in zip(problems[0], problems[1]):
        assert_allclose(actual, expected, rtol=0, atol=0)
    # Static reaction stresses need not be unique; compare velocity to solve
    # accuracy rather than requiring identical null-space stresses.
    assert_allclose(L.materials[0].v, reference[2].materials[0].v, rtol=0, atol=1e-8)


def test_periodic_column_preserves_mass_and_species():
    from inputs.granular_column import Params

    P = Params([])
    P.interpolation = "linear"
    P.transfer = "flip"
    P.B = NS(cyclic_lr=True, has_bottom=True, roughness=True)
    P.G.ns = 2
    P.G.s = [0.001, 0.003]
    P.S.PHI = np.tile([0.25, 0.75], (P.S.n, 1))
    P.segregate = True
    P.c = 1e-4
    P, G, L = prepare(P)
    initial = (L.materials[0].m[:, None] * L.materials[0].phi).sum(axis=0)
    run(P, G, L, max_steps=3, output=False)
    b = L.materials[0]
    assert_allclose((b.m[:, None] * b.phi).sum(axis=0), initial, atol=1e-10)
    assert_allclose(b.phi.sum(axis=1), 1, atol=1e-12)
    assert_allclose(G.q[:: P.G.nx], G.q[P.G.nx - 1 :: P.G.nx], atol=1e-12)


def test_column_example_saves_figures_and_solver_diagnostics(tmp_path):
    P, G, L = get_parameters(["granular_column"])
    P.save_dir = str(tmp_path)
    run(P, G, L, max_steps=1)
    assert (tmp_path / "Continuum/00000.png").is_file()
    assert (tmp_path / "MP/00000.png").is_file()
    assert (tmp_path / "data/diagnostics.json").is_file()


def test_nonconvergence_is_not_silently_accepted():
    with pytest.raises(RuntimeError, match="did not converge"):
        solve_coulomb(np.eye(4), [-1, 1, 0, 0], 0.5, max_iterations=1)


@pytest.mark.parametrize(
    "setting",
    [dict(quadrature="domain"), dict(interpolation="quadratic"), dict(damping=True), dict(solver="invalid")],
)
def test_unsupported_formulations_fail_early(setting):
    with pytest.raises(ValueError):
        column(**setting)


def test_block_contact_fallback_solves_when_global_newton_budget_is_exhausted():
    W = np.eye(6)
    W[:3, 3:] = W[3:, :3] = 0.1 * np.eye(3)
    q = np.array([[-1.0, 1.0, 0.0], [-1.0, -0.8, 0.0]])
    r, stats = solve_coulomb(W, q, 0.5, max_iterations=1)
    w = (W @ r.ravel()).reshape(-1, 3) + q
    assert stats["sweeps"] > 0
    assert_allclose(w[:, 0], 0, atol=1e-9)
    assert_allclose(
        r[:, 1:], -0.5 * r[:, 0, None] * w[:, 1:] / np.linalg.norm(w[:, 1:], axis=1)[:, None], atol=1e-9
    )


def test_tiny_timestep_does_not_hide_static_contact():
    P, G, L = column(dt=1e-6, t_f=2e-6)
    run(P, G, L, output=False)
    assert np.max(np.abs(L.materials[0].v)) < 1e-8
    assert L.materials[0].pressure.min() < -100


def test_collapse_yields_and_contact_work_is_dissipative():
    P, G, L = get_parameters(["granular_collapse"])
    P.dt = 0.01
    initial = L.materials[0].x.copy()
    run(P, G, L, output=False)
    b = L.materials[0]
    assert b.x[:, 0].max() > initial[:, 0].max() + 0.05
    assert np.linalg.norm(b.v, axis=1).max() > 0.1
    work = sum(P.dt * r["granular_solver"]["contact_power"] for r in L.diagnostics.history)
    assert work < 0
    assert_allclose(b.plastic_work.sum(), -work, rtol=1e-9, atol=1e-8)
    assert max(r["granular_solver"]["contact_power"] for r in L.diagnostics.history) <= 1e-7
    assert np.linalg.norm(L.diagnostics.last["momentum_residual"]) < 1e-8
