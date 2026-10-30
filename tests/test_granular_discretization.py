"""Q1-P0 operator checks; spatial accuracy is studied separately in benchmarks."""

from copy import deepcopy

import numpy as np
import pytest
from numpy.testing import assert_allclose

from benchmarks.granular_mixed import parameters
from granular import BASIS, PLANE_BASIS, assemble, cell_indices, cell_packing_fraction
from initialise import prepare
from MPM import run


@pytest.mark.parametrize("stress,basis", [("intrinsic", BASIS), ("plane_strain", PLANE_BASIS)])
@pytest.mark.parametrize("periodic", [False, True])
def test_p0_assembly_matches_independent_particle_sum_and_compact_operator(stress, basis, periodic):
    P = parameters("static", "velocity")
    P.granular_stress_basis, P.granular_stress = "P0", stress
    if periodic:
        P.B.has_left = P.B.has_right = False
        P.B.cyclic_lr = True
    P, G, L = prepare(P)
    b = L.materials[0]
    b.x[:, :2] += np.random.default_rng(9).uniform(-0.01, 0.01, (len(b), 2))
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    B, free, active, occupied, gap = assemble(P, G, b)
    compact, *other = assemble(P, G, b, compact=True)
    for actual, expected in zip(other, (free, active, occupied, gap)):
        assert_allclose(actual, expected, atol=0, rtol=0)
    expected = np.zeros(B.shape)
    for p in range(len(b)):
        ix = min(int((b.x[p, 0] - P.G.x_m) / G.dx), P.G.nx - 2)
        iy = min(int((b.x[p, 1] - P.G.y_m) / G.dy), P.G.ny - 2)
        row = int(np.flatnonzero(active == ix + (P.G.nx - 1) * iy)[0])
        for j, node in enumerate(b.n[p]):
            if periodic and node % P.G.nx == P.G.nx - 1:
                node -= P.G.nx - 1
            for k, tensor in enumerate(basis):
                for c in (0, 1):
                    expected[len(basis) * row + k, 2 * node + c] += (
                        b.m[p] / b.material.rho_s * (tensor[c] @ b.G[p, j])
                    )
    assert_allclose(B.toarray(), expected, atol=2e-17)
    rng = np.random.default_rng(14)
    u, r = rng.normal(size=B.shape[1]), rng.normal(size=B.shape[0])
    assert_allclose(compact @ u, B @ u, atol=1e-16)
    assert_allclose(compact.T @ r, B.T @ r, atol=1e-16)
    # Constant velocity has zero strain even when the pressure basis is P0.
    assert_allclose(compact @ np.tile([0.2, -0.3], len(G.m)), 0, atol=1e-17)


def test_p0_integrates_affine_velocity_gradient_and_preserves_grain_volume():
    P = parameters("static", "velocity")
    P.granular_stress_basis = "P0"
    P, G, L = prepare(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    b = L.materials[0]
    B, _, active, occupied, gap = assemble(P, G, b, compact=True)
    C = np.array([[0.2, -0.1], [0.3, -0.2]])
    u = np.column_stack((G.X, G.Y)) @ C.T
    expected = occupied[active, None] * np.einsum("kij,ij->k", BASIS[:, :2, :2], C)
    assert_allclose((B @ u.ravel()).reshape(-1, 3), expected, atol=1e-16)
    assert_allclose(gap, 0, atol=1e-14)
    phi = cell_packing_fraction(P, G, b)
    assert_allclose(phi.sum() * G.dx * G.dy * P.G.thickness, np.sum(b.m / b.material.rho_s), atol=1e-15)
    assert_allclose(phi[active], b.material.phi_max, atol=1e-15)
    assert len(B.colors(np.ones((len(G.m), 2)))) == 4


@pytest.mark.parametrize("periodic", [False, True])
def test_cell_ownership_at_upper_edge_and_periodic_seam(periodic):
    P = parameters("static", "velocity")
    P.B.cyclic_lr = periodic
    P, G, _ = prepare(P)
    points = np.array([[P.G.x_m, P.G.y_m, 0], [P.G.x_M, P.G.y_M, 0]])
    expected_last = (P.G.nx - 1) * (P.G.ny - 2) + (0 if periodic else P.G.nx - 2)
    assert_allclose(cell_indices(P, G, points), [0, expected_last], atol=0)


@pytest.mark.parametrize("stress", ["intrinsic", "plane_strain"])
def test_p0_static_support_and_pressure_are_finite(stress):
    P = parameters("static", "velocity")
    P.granular_stress_basis, P.granular_stress = "P0", stress
    P, G, L = prepare(P)
    before = L.materials[0].x.copy()
    run(P, G, L, output=False)
    b = L.materials[0]
    assert_allclose(b.x, before, atol=1e-10)
    assert np.max(np.abs(b.v)) < 1e-8
    assert b.pressure.min() < -100
    assert b.pressure.max() <= 1e-6
    assert G.granular_contact.shape == ((P.G.nx - 1) * (P.G.ny - 1), 3 if stress == "intrinsic" else 4)


@pytest.mark.parametrize("stress", ["intrinsic", "plane_strain"])
def test_p0_evolving_collapse_matches_assembled_solver(stress):
    P = parameters("collapse", "hybrid")
    P.granular_stress_basis, P.granular_stress = "P0", stress
    P.t_f = 0.05
    reference = prepare(P)
    P, G, L = deepcopy(reference)
    P.granular_solver = "velocity"
    run(P, G, L, output=False)
    run(*reference, output=False)
    assert_allclose(L.materials[0].x, reference[2].materials[0].x, atol=1e-8)
    assert_allclose(L.materials[0].v, reference[2].materials[0].v, atol=1e-7)
    assert max(np.linalg.norm(r["momentum_residual"]) for r in L.diagnostics.history) < 1e-8
    assert max(r["granular_solver"]["contact_power"] for r in L.diagnostics.history) <= 1e-6


def test_p0_periodic_segregation_conserves_species():
    P = parameters("incline_static", "velocity")
    P.granular_stress_basis = "P0"
    P.G.ns, P.G.s = 2, [0.001, 0.003]
    P.S.phi = np.tile([0.25, 0.75], (P.S.n, 1))
    P.segregate, P.c = True, 1e-4
    P, G, L = prepare(P)
    b = L.materials[0]
    before = (b.m[:, None] * b.phi).sum(axis=0)
    run(P, G, L, max_steps=3, output=False)
    assert_allclose((b.m[:, None] * b.phi).sum(axis=0), before, atol=1e-10)
    assert_allclose(b.phi.sum(axis=1), 1, atol=1e-12)
    assert_allclose(G.q[:: P.G.nx], G.q[P.G.nx - 1 :: P.G.nx], atol=1e-12)


def test_p0_preserves_plot_formats_and_save_cadence(tmp_path):
    P = parameters("static", "velocity")
    P.granular_stress_basis = "P0"
    P.dt, P.t_f, P.savetime = 0.005, 0.02, 0.01
    P.save_dir = str(tmp_path)
    P, G, L = prepare(P)
    run(P, G, L)
    # Saves at .01 and .02, followed by the existing at_end output callback.
    assert len(list((tmp_path / "Continuum").glob("*.png"))) == 3
    assert len(list((tmp_path / "MP").glob("*.png"))) == 3
    assert (tmp_path / "data/diagnostics.json").is_file()
