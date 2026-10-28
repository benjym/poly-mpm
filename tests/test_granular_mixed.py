"""Independent mechanics and reference-solver checks for the compact iteration."""

from copy import deepcopy

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import sparse

from benchmarks.granular_prototype import configuration
from granular import BASIS, PLANE_BASIS, assemble, solve_coulomb
from granular_mixed import GradientOperator, diagonal_coulomb, solve_velocity
from initialise import prepare
from MPM import run, step


@pytest.mark.parametrize("components", [3, 4, 6])
def test_diagonal_local_solver_satisfies_nonassociated_coulomb(components):
    rng = np.random.default_rng(36)
    diagonal = 10 ** rng.uniform(-3, 3, (1000, components))
    q = rng.normal(size=diagonal.shape)
    mu = rng.uniform(0, 1, len(q))
    mu[::7] = 0
    r = diagonal_coulomb(diagonal, q, mu)
    w = diagonal * r + q
    assert np.min(r[:, 0]) >= 0
    assert np.min(w[:, 0]) > -1e-12
    assert_allclose(r[:, 0] * w[:, 0], 0, atol=1e-12)
    assert np.all(np.linalg.norm(r[:, 1:], axis=1) <= mu * r[:, 0] * (1 + 3e-13) + 1e-14)
    moving = np.linalg.norm(w[:, 1:], axis=1) > 1e-8
    expected = (
        -mu[moving, None] * r[moving, :1] * w[moving, 1:] / np.linalg.norm(w[moving, 1:], axis=1)[:, None]
    )
    assert_allclose(r[moving, 1:], expected, rtol=3e-10, atol=2e-11)


@pytest.mark.parametrize(
    "q,expected",
    [
        ([-1, 0.2, 0], [1, -0.2, 0]),
        ([-1, 1, 0], [1, -0.5, 0]),
        ([1, 1, 0], [0, 0, 0]),
        ([-1, 0, 0], [1, 0, 0]),
    ],
)
def test_isotropic_local_branches(q, expected):
    assert_allclose(diagonal_coulomb(np.ones((1, 3)), [q], 0.5)[0], expected, atol=1e-14)


@pytest.mark.parametrize("stress", ["intrinsic", "plane_strain"])
@pytest.mark.parametrize("periodic", [False, True])
def test_compact_assembly_and_transpose_match_original(stress, periodic):
    P = configuration("static", 0.005)
    P.granular_stress = stress
    if periodic:
        P.B.has_left = P.B.has_right = False
        P.B.cyclic_lr = True
    P, G, L = prepare(P)
    b = L.materials[0]
    # Break tensor-product symmetry and exercise particles crossing cells.
    b.x[:, :2] += np.random.default_rng(7).uniform(-0.01, 0.01, (len(b), 2))
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    B, *reference = assemble(P, G, b)
    compact, *actual = assemble(P, G, b, compact=True)
    for a, e in zip(actual, reference):
        assert_allclose(a, e, atol=0, rtol=0)
    rng = np.random.default_rng(6)
    u, r = rng.normal(size=B.shape[1]), rng.normal(size=B.shape[0])
    assert_allclose(compact @ u, B @ u, atol=1e-16)
    assert_allclose(compact.T @ r, B.T @ r, atol=1e-16)
    assert_allclose(r @ (compact @ u), u @ (compact.T @ r), atol=1e-16)
    inverse_mass = np.repeat(np.where(G.valid_nodes, 1.0, 0.0), 2).reshape(-1, 2)
    inverse_mass.ravel()[~actual[0]] = 0
    for group in compact.colors(inverse_mass):
        # Distinct nodes of the same color never write the same free velocity.
        indices = compact.nodes[group][np.any(compact.gradients[group] != 0, axis=-1)]
        indices = indices[np.any(inverse_mass[indices] > 0, axis=1)]
        assert len(np.unique(indices)) == len(indices)


def reference_residual(W, q, r, mu):
    """Original hybrid solver's dimensionless infinity-norm natural map."""
    norm = float(np.asarray(abs(W).sum(axis=1)).max())
    magnitude = max(np.max(-q[:, 0], initial=0), np.max(np.abs(q[:, 1:]), initial=0), 1e-30)
    z = r * norm / magnitude
    w = (W @ r.ravel()).reshape(q.shape) / magnitude + q / magnitude
    projected = z - w
    projected[:, 0] = np.maximum(projected[:, 0], 0)
    radius = mu * np.maximum(z[:, 0], 0)
    length = np.linalg.norm(projected[:, 1:], axis=1)
    projected[:, 1:] *= np.minimum(1, np.divide(radius, length, out=np.ones(len(q)), where=length > 0))[
        :, None
    ]
    return np.max(np.abs(z - projected))


@pytest.mark.parametrize("basis", [BASIS, PLANE_BASIS])
def test_matrix_free_solve_retains_full_coupling_and_reference_tolerance(basis):
    rng = np.random.default_rng(12)
    count, nodes = 5, 8
    rows, columns = np.indices((count, nodes))
    operator = GradientOperator(
        rows.ravel(), columns.ravel(), rng.normal(size=(count * nodes, 2)), count, nodes, basis
    )
    # This dense matrix exists only in the independent small-problem oracle.
    B = np.column_stack([operator @ e for e in np.eye(2 * nodes)])
    mass_inverse = rng.uniform(0.3, 1, 2 * nodes)
    mass_inverse[:2] = 0  # Wall-constrained degrees of freedom.
    trial = rng.normal(size=2 * nodes)
    gap = np.zeros((count, len(basis)))
    gap[:, 0] = rng.uniform(0, 0.1, count)
    dt, mu = 0.02, 0.3
    W = sparse.csr_matrix(dt * (B * mass_inverse) @ B.T)
    q = (B @ trial).reshape(gap.shape) + gap
    reference, _ = solve_coulomb(W, q, mu)
    result, stats = solve_velocity(operator, trial, mass_inverse, gap, mu, dt)
    assert reference_residual(W, q, result, mu) <= 1e-9
    # Stress multipliers can be nonunique; velocity is the physical comparison.
    assert_allclose(
        dt * mass_inverse * (B.T @ result.ravel()), dt * mass_inverse * (B.T @ reference.ravel()), atol=2e-8
    )
    assert stats["assembled_delassus"] is False


@pytest.mark.parametrize("stress", ["intrinsic", "plane_strain"])
def test_evolving_collapse_agrees_with_hybrid(stress):
    P = configuration("collapse", 0.005)
    P.granular_stress = stress
    P.transfer = "apic_gradient"
    P.granular_viscosity = 0.001
    reference = prepare(P)
    P, G, L = deepcopy(reference)
    P.granular_solver = "velocity"
    run(P, G, L, output=False)
    run(*reference, output=False)
    assert_allclose(L.materials[0].x, reference[2].materials[0].x, atol=2e-8)
    assert_allclose(L.materials[0].v, reference[2].materials[0].v, atol=2e-7)
    assert max(r["granular_solver"]["residual"] for r in L.diagnostics.history) <= P.granular_tolerance
    assert max(np.linalg.norm(r["momentum_residual"]) for r in L.diagnostics.history) < 1e-8


@pytest.mark.parametrize("slope", [0.25, 0.75])
def test_inclined_layer_remains_static_or_slides_at_coulomb_rate(slope):
    P = configuration("static", 0.005, solver="velocity")
    P.theta = np.arctan(slope)
    P.B.has_left = P.B.has_right = False
    P.B.cyclic_lr = P.B.roughness = True
    P.t_f = 0.05
    P, G, L = prepare(P)
    run(P, G, L, output=False)
    speed = -L.materials[0].v[:, 0]
    if slope < 0.5:
        assert np.max(np.abs(L.materials[0].v)) < 1e-8
    else:
        expected = 9.81 * (np.sin(P.theta) - 0.5 * np.cos(P.theta)) * P.t
        assert_allclose(speed.max(), expected, rtol=1e-3)
        assert speed.mean() > 0.05


def test_compact_path_does_not_call_global_reference_solver(monkeypatch):
    import granular

    def forbidden(*args, **kwargs):
        pytest.fail("compact path called the assembled global solver")

    monkeypatch.setattr(granular, "solve_coulomb", forbidden)
    P, G, L = prepare(configuration("static", 0.02, solver="velocity"))
    step(P, G, L)
    assert G.granular_stats["assembled_delassus"] is False
    assert G.granular_stats["evaluations"] == 0


def test_compact_solver_rejects_unsatisfied_budget():
    P, G, L = prepare(configuration("static", 0.02, solver="velocity"))
    P.granular_max_iterations = 1
    with pytest.raises(RuntimeError, match="did not converge"):
        step(P, G, L)


def test_large_opening_gap_does_not_mask_compression():
    nodes = np.array([0, 1])
    operator = GradientOperator(nodes, nodes, np.array([[1.0, 0], [1.0, 0]]), 2, 2, BASIS)
    trial = np.array([-1e-10, 0.0, 0.0, 0.0])
    gap = np.zeros((2, 3))
    gap[1, 0] = 1e10
    result, _ = solve_velocity(operator, trial, np.ones(4), gap, 0.5, 1.0)
    final = trial + operator.T @ result.ravel()
    assert abs(final[0]) < 1e-18
    assert_allclose(result[1], 0, atol=1e-22)


def test_frozen_sliding_incline_velocity_agrees_under_tolerance_refinement():
    from benchmarks.granular_mixed import frozen_comparison

    comparison = frozen_comparison("incline_sliding")
    assert comparison["residual_passed"]
    assert comparison["velocity_equivalence_verified"]
    if not comparison["velocity_match"]:
        study = comparison["tolerance_study"]
        assert study[-1]["tolerance"] < comparison["requested_tolerance"]
        assert study[-1]["max_grid_velocity_difference"] < comparison["velocity_acceptance_bound"]
