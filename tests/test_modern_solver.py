"""Independent verification of updated mechanics, stencils, and matched transfers."""

from copy import deepcopy
from types import SimpleNamespace as NS
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.linalg import expm
from conftest import state, material
from constit import update_stress, EYE
from kinematics import advance
from interpolation import stencil, moment, affine_matrix
from diagnostics import inventory
from stability import limits
from MPM import run, step


def mapped(**options):
    P, G, L = state(**options)
    b = L.materials[0]
    b.get_reference_nodes(P, G)
    b.get_basis_functions(P, G)
    return P, G, L, b


@pytest.mark.parametrize("basis", ["linear", "quadratic"])
def test_stencil_complete_at_edges_corners_and_interior(basis):
    rng = np.random.default_rng(7)
    xy = np.vstack((rng.uniform(0, 1, (30, 2)), [[0, 0], [1, 1], [0, 0.125], [1, 0.875], [0.5, 0]]))
    P, G, L, b = mapped(interpolation=basis, materials=material(n=len(xy), X=xy[:, 0], Y=xy[:, 1]))
    assert b.N.shape == (len(xy), 4 if basis == "linear" else 9)
    assert np.all(b.N >= -1e-14)
    assert_allclose(b.N.sum(axis=1), 1, atol=1e-14)
    assert_allclose(b.G.sum(axis=1), 0, atol=1e-13)
    assert_allclose(np.einsum("nr,nri->ni", b.N, b.offsets), 0, atol=1e-14)
    derivative = np.einsum("nri,nrj->nij", b.offsets, b.G)
    assert_allclose(derivative, np.broadcast_to(np.diag([1, 1, 0]), derivative.shape), atol=1e-13)
    L.recover_position(P, G)


@pytest.mark.parametrize("basis", ["linear", "quadratic"])
def test_stencil_gradient_matches_weight_derivative(basis):
    P, G, L, b = mapped(interpolation=basis)
    b.x[:, :2] += 0.0123
    b.get_reference_nodes(P, G)
    b.get_basis_functions(P, G)
    delta = 1e-7
    for axis in (0, 1):
        x1, x2 = b.x.copy(), b.x.copy()
        x1[:, axis] += delta
        x2[:, axis] -= delta
        n1, w1, _, _ = stencil(x1, P, G)
        n2, w2, _, _ = stencil(x2, P, G)
        # Accumulate duplicates at modified physical boundaries before comparing.
        numerical = np.zeros((len(b), len(G.m)))
        expected = np.zeros_like(numerical)
        for i in range(len(b)):
            np.add.at(numerical[i], n1[i], w1[i] / (2 * delta))
            np.add.at(numerical[i], n2[i], -w2[i] / (2 * delta))
            np.add.at(expected[i], b.n[i], b.G[i, :, axis])
        assert_allclose(numerical, expected, atol=2e-9)


def test_quadratic_gradient_is_continuous_across_grid_crossing():
    differences = {}
    for basis in ("linear", "quadratic"):
        P, G, L, b = mapped(interpolation=basis, materials=material(n=1, X=[0.5], Y=[0.37]))
        gradients = []
        for offset in (-1e-8, 1e-8):
            b.x[0, 0] = 0.5 + offset
            b.get_reference_nodes(P, G)
            b.get_basis_functions(P, G)
            g = np.zeros((len(G.m), 3))
            np.add.at(g, b.n[0], b.G[0])
            gradients.append(g)
        differences[basis] = np.linalg.norm(gradients[0] - gradients[1])
    assert differences["quadratic"] < 1e-5
    assert differences["linear"] > 1


@pytest.mark.parametrize(
    "gradient",
    [
        np.diag([-0.2, 0.1, 0]),
        np.array([[0, -2.0, 0], [2.0, 0, 0], [0, 0, 0]]),
        np.array([[0, 0.3, 0], [0, 0, 0], [0, 0, 0]]),
    ],
)
def test_affine_deformation_volume_and_density(gradient):
    P, G, L, b = mapped(dt=0.1)
    b.L[:] = gradient
    initial_mass = b.m.copy()
    for _ in range(10):
        advance(b, P)
    exact = expm(gradient)
    assert_allclose(b.F, np.broadcast_to(exact, b.F.shape), atol=2e-14)
    assert_allclose(b.J, np.exp(np.trace(gradient)))
    assert_allclose(b.V, b.V0 * b.J)
    assert_allclose(b.rho * b.V, initial_mass)
    assert_allclose(b.m, initial_mass)


def test_invalid_deformation_rejected_with_location():
    P, G, L, b = mapped()
    b.F[2, 0, 0] = -1
    with pytest.raises(FloatingPointError, match="particles \\[2\\]"):
        advance(b, P)


def test_elastic_rotation_is_objective_and_zero_energy():
    P, G, L, b = mapped(dt=0.1)
    W = np.array([[0, -3.0, 0], [3.0, 0, 0], [0, 0, 0]])
    b.L[:] = W
    b.dstrain.fill(0)
    for _ in range(20):
        advance(b, P)
        update_stress(b, b.material, P)
    assert_allclose(b.J, 1, atol=1e-14)
    assert_allclose(b.stress, 0, atol=1e-12)
    assert_allclose(b.elastic_energy, 0, atol=1e-13)
    assert_allclose(b.rho, b.initial_rho, atol=1e-10)


def test_elastic_energy_derivative_matches_first_piola_stress():
    P, G, L, b = mapped(materials=material(n=1))
    F = np.array([[1.1, 0.2, 0], [0, 0.9, 0], [0, 0, 1.0]])
    direction = np.array([[0.3, -0.2, 0], [0.1, 0.2, 0], [0, 0, 0]])
    b.F[:] = F
    b.J[:] = np.linalg.det(F)
    b.V = b.V0 * b.J
    update_stress(b, b.material, P)
    first_piola = b.J[0] * b.stress[0] @ np.linalg.inv(F).T
    values = []
    for sign in (-1, 1):
        other = deepcopy(b)
        other.F[:] = F + sign * 1e-6 * direction
        other.J[:] = np.linalg.det(other.F)
        other.V = other.V0 * other.J
        update_stress(other, other.material, P)
        values.append(other.elastic_energy[0] / other.V0[0])
    assert_allclose((values[1] - values[0]) / 2e-6, np.sum(first_piola * direction), rtol=1e-8)


def test_viscous_and_zero_yield_bingham_agree_and_dissipate():
    outcomes = []
    for name in ("viscous", "bingham"):
        P, G, L, b = mapped(materials=material(name, tau_0=0.0, mu_0=100.0))
        b.L[:] = [[-0.2, 1.0, 0], [0.2, -0.1, 0], [0, 0, 0]]
        b.dstrain = 0.5 * (b.L + b.L_T) * P.dt
        advance(b, P)
        update_stress(b, b.material, P)
        outcomes.append(b)
        assert np.all(b.viscous_work > 0)
    assert_allclose(outcomes[0].stress, outcomes[1].stress)
    assert_allclose(outcomes[0].viscous_work, outcomes[1].viscous_work)


def test_j2_return_unloading_and_plastic_work():
    P, G, L, b = mapped(materials=material("von_mises", yield_stress=2.0))
    b.dstrain[:] = [[0, 0.1, 0], [0.1, 0, 0], [0, 0, 0]]
    update_stress(b, b.material, P)
    assert_allclose(b.q, 2.0)
    assert np.all(b.plastic_work > 0)
    before = b.plastic_work.copy()
    b.dstrain *= -0.001
    update_stress(b, b.material, P)
    assert np.all(b.q < 2.0)
    assert_allclose(b.plastic_work, before)


def affine_case():
    rng = np.random.default_rng(9)
    xy = rng.uniform(0.2, 0.8, (50, 2))
    A = np.array([[0.1, -0.3, 0], [0.2, 0.05, 0], [0, 0, 0.0]])
    offset = np.array([0.4, -0.1, 0.0])
    xyz = np.column_stack((xy, np.zeros(len(xy))))
    return (
        mapped(
            interpolation="quadratic",
            transfer="apic",
            materials=material(
                n=len(xy),
                X=xy[:, 0],
                Y=xy[:, 1],
                initial_v=xyz @ A.T + offset,
                initial_C=A,
                rho=rng.uniform(500, 1500, len(xy)),
                A=rng.uniform(0.005, 0.015, len(xy)),
            ),
        ),
        A,
        offset,
    )


def test_apic_affine_roundtrip_and_grid_field():
    (P, G, L, b), A, offset = affine_case()
    L.get_nodal_mass_momentum(P, G)
    expected = np.column_stack((G.X, G.Y, np.zeros_like(G.X))) @ A.T + offset
    assert_allclose(G.q[G.valid_nodes] / G.m[G.valid_nodes, None], expected[G.valid_nodes], atol=1e-14)
    old_v = b.v.copy()
    old_B = b.affine_momentum.copy()
    b.move_material_points(P, G)
    assert_allclose(b.v, old_v, atol=1e-14)
    assert_allclose(b.affine_momentum, old_B, atol=1e-14)


def test_apic_matched_transfers_conserve_full_angular_momentum_and_filter_energy():
    (P, G, L, b), A, offset = affine_case()
    rng = np.random.default_rng(32)
    b.v[:] = rng.normal(size=b.v.shape)
    b.v[:, 2] = 0
    b.affine_momentum[:] = rng.normal(size=b.affine_momentum.shape) * 0.01
    b.affine_momentum[:, :, 2] = 0
    b.affine_momentum[:, 2, :] = 0
    before = inventory(P, L)
    L.get_nodal_mass_momentum(P, G)
    coordinates = np.column_stack((G.X, G.Y, np.zeros_like(G.X)))
    assert_allclose(G.q.sum(axis=0), before["momentum"], atol=1e-12)
    assert_allclose(np.cross(coordinates, G.q).sum(axis=0), before["angular_momentum"], atol=1e-12)
    grid_energy = 0.5 * np.sum(np.sum(G.q**2, axis=1)[G.valid_nodes] / G.m[G.valid_nodes])
    b.move_material_points(P, G)
    after = inventory(P, L)
    assert_allclose(after["momentum"], before["momentum"], atol=1e-12)
    assert_allclose(after["angular_momentum"], before["angular_momentum"], atol=1e-12)
    assert after["kinetic"] + after["affine_kinetic"] <= grid_energy + 1e-12
    assert grid_energy <= before["kinetic"] + before["affine_kinetic"] + 1e-12


@pytest.mark.parametrize("xy", [[0, 0], [0, 0.37], [1, 1], [0.125, 0.875]])
def test_apic_rank_deficient_boundary_support_is_finite_for_resolved_state(xy):
    P, G, L, b = mapped(
        interpolation="quadratic",
        transfer="apic",
        materials=material(n=1, X=[xy[0]], Y=[xy[1]], initial_v=[0.1, 0.2, 0]),
    )
    L.get_nodal_mass_momentum(P, G)
    old = b.v.copy()
    b.move_material_points(P, G)
    assert_allclose(b.v, old)
    assert np.isfinite(b.affine_momentum).all()


def test_apic_refuses_to_discard_unresolved_affine_momentum():
    P, G, L, b = mapped(interpolation="quadratic", transfer="apic", materials=material(n=1, X=[0.0], Y=[0.0]))
    b.affine_momentum[0, 0, 0] = 1.0
    with pytest.raises(FloatingPointError, match="support lost"):
        L.get_nodal_mass_momentum(P, G)


@pytest.mark.parametrize("quadrature", ["point", "domain"])
@pytest.mark.parametrize("transfer", ["flip", "pic", "apic"])
@pytest.mark.parametrize("periodic", [False, True])
def test_new_stencil_transport_preserves_species_and_mass(transfer, periodic, quadrature):
    rng = np.random.default_rng(10)
    xy = rng.uniform(0.15, 0.85, (20, 2))
    phi = rng.dirichlet([1, 1, 1], len(xy))
    P, G, L, b = mapped(
        interpolation="quadratic",
        transfer=transfer,
        quadrature=quadrature,
        ns=3,
        segregate=True,
        c=1e-3,
        materials=material(
            "pouliquen",
            n=20,
            X=xy[:, 0],
            Y=xy[:, 1],
            PHI=phi,
            domain_size=[0.1, 0.1],
            initial_v=np.column_stack((0.1 * xy[:, 1], np.zeros((20, 2)))),
            initial_C=[[0, 0.1, 0], [0, 0, 0], [0, 0, 0]],
        ),
        pressure="compression",
        initial_pressure=10.0,
        B=NS(cyclic_lr=periodic),
    )
    before = (b.m[:, None] * b.phi).sum(axis=0)
    initial_phi = b.phi.copy()
    run(P, G, L, max_steps=20, output=False)
    assert np.max(np.abs(b.phi - initial_phi)) > 1e-12
    assert_allclose((b.m[:, None] * b.phi).sum(axis=0), before, atol=1e-12)
    assert_allclose(b.phi.sum(axis=1), 1, atol=1e-14)
    assert np.all(b.phi >= 0) and np.all(b.phi <= 1)
    assert L.diagnostics.history[-1]["mass_error"] == 0


def test_uniform_acceleration_diagnostics_account_for_work_and_impulse():
    P, G, L, b = mapped(interpolation="quadratic", transfer="apic", max_g=-2.0)
    initial = inventory(P, L)
    run(P, G, L, max_steps=20, output=False)
    row = L.diagnostics.history[-1]
    assert_allclose(b.v[:, 1], -2 * P.t, atol=1e-14)
    assert_allclose(np.asarray(row["momentum"]) - initial["momentum"], row["external_impulse"], atol=1e-13)
    assert_allclose(row["kinetic"], row["external_work"], atol=1e-15)
    assert abs(row["energy_residual"]) < 1e-12
    assert row["energy_budget_closed"]


def test_stability_reports_and_enforcement_are_separate_from_static_control():
    P, G, L, b = mapped(materials=material("viscous", mu_s=1000.0, mu_v=20.0))
    restriction = limits(P, G, L)
    assert restriction["viscosity"] < restriction["elastic_wave"]
    P.dt = restriction["viscosity"] * 2
    P.timestep_policy = "enforce"
    with pytest.raises(ValueError, match="estimated stable"):
        step(P, G, L)


def test_dynamic_first_step_obeys_viscous_limit():
    P, G, L, b = mapped(materials=material("viscous", mu_s=1000.0), time_stepping="dynamic", dt=1.0)
    safe = min(limits(P, G, L).values())
    run(P, G, L, max_steps=1, output=False)
    assert P.t <= safe * (1 + 1e-14)


def test_no_tiny_roundoff_final_step():
    P, G, L = state(dt=1e-4, t_f=0.01)
    run(P, G, L, output=False)
    assert P.tstep == 100
    assert P.t == P.t_f


def test_standing_wave_refinement_and_transfer_dissipation():
    from benchmarks.verification import wave

    rows = [wave(n) for n in (8, 16, 32)]
    # Temporal/spatial errors can cancel at individual levels; do not demand a
    # fabricated monotone order. The fine solution must improve substantially.
    assert rows[-1]["velocity_relative_l2"] < rows[0]["velocity_relative_l2"] * 0.2
    assert rows[-1]["velocity_relative_l2"] < 0.02
    assert abs(rows[-1]["relative_energy_residual"]) < abs(rows[0]["relative_energy_residual"]) * 0.1
    assert max(row["maximum_momentum_error"] for row in rows) < 1e-15
    pic = wave(16, transfer="pic")
    assert abs(rows[1]["relative_energy_residual"]) < abs(pic["relative_energy_residual"]) * 0.1


def test_energy_and_impulse_ledger_tracks_fixed_wall_constraint():
    P, G, L, b = mapped(
        materials=material(n=1, X=[0.4], Y=[0.0], initial_v=[0.0, -1.0, 0.0]), B=NS(has_bottom=True)
    )
    run(P, G, L, max_steps=1, output=False)
    row = L.diagnostics.history[-1]
    assert_allclose(row["momentum_residual"], 0, atol=1e-13)
    assert row["constraint_work"] < 0
    assert_allclose(row["energy_residual"], 0, atol=1e-13)


def test_diagnostics_json_and_sampling(tmp_path):
    P, G, L = state(save_dir=str(tmp_path), diagnostics_stride=2)
    run(P, G, L, max_steps=3)
    import json

    report = json.loads((tmp_path / "data/diagnostics.json").read_text())
    assert report["history"][0]["timestep"] == 2
    assert report["history"][0]["mass_error"] == 0


def test_updated_density_is_independent_of_grid_occupancy():
    P, G, L, b = mapped(materials=material(initial_v=[0.1, 0, 0]))
    G.m.fill(999.0)
    G.q[:] = G.m[:, None] * [0.1, 0, 0]
    G.get_valid_nodes(P)
    before = b.rho.copy()
    b.move_material_points(P, G)
    assert_allclose(b.rho, before)


def test_boundary_injection_initialises_deformation_and_affine_history():
    P, G, L, b = mapped(interpolation="quadratic", transfer="apic")
    b.F[:] = 2 * EYE
    b.affine_momentum[:] = 0.1
    template = L.templates[0]
    template.inlet_rate, template.gap, template.yp = 0.1, [0.01, 0.01], [0.2]
    P.B.inlet_right, P.t = True, 0.1
    L.process_boundaries(P, G)
    assert_allclose(b.F[-1], EYE)
    assert b.J[-1] == 1
    assert_allclose(b.affine_momentum[-1], 0)
    assert b.V[-1] == b.V0[-1]


def test_torque_free_apic_body_preserves_angular_momentum_through_full_steps():
    (P, G, L, b), A, offset = affine_case()
    A = np.array([[0, -0.2, 0], [0.2, 0, 0], [0, 0, 0.0]])
    b.v[:] = (b.x - [0.5, 0.5, 0]) @ A.T
    b.affine_momentum = A @ moment(b.N, b.offsets)
    # Reset diagnostics after prescribing this initial condition.
    from diagnostics import Diagnostics

    L.diagnostics = Diagnostics(P, L)
    initial = inventory(P, L)
    run(P, G, L, max_steps=40, output=False)
    final = inventory(P, L)
    assert_allclose(final["angular_momentum"], initial["angular_momentum"], atol=2e-13)
    assert_allclose(L.diagnostics.history[-1]["angular_momentum_residual"], 0, atol=2e-13)


def test_bulk_fluid_energy_integrates_volume_work():
    P, G, L, b = mapped(materials=material("viscous", K=100.0, mu_s=0.0, mu_v=0.0), dt=0.1)
    b.L[:] = -0.1 * EYE
    b.dstrain[:] = -0.01 * EYE
    advance(b, P)
    update_stress(b, b.material, P)
    J = np.exp(-0.03)
    expected = b.V0 * 100 * (J * np.log(J) - J + 1)
    assert_allclose(b.elastic_energy, expected, atol=1e-13)
    assert_allclose(b.pressure, -3.0)
    assert_allclose(b.viscous_work, 0)


@pytest.mark.parametrize("basis,transfer", [("linear", "flip"), ("quadratic", "apic")])
def test_translation_at_sparse_corner_does_not_generate_strain(basis, transfer):
    P, G, L, b = mapped(
        interpolation=basis,
        transfer=transfer,
        materials=material(n=1, X=[0.0], Y=[0.0], initial_v=[0.1, 0.2, 0.0]),
    )
    step(P, G, L)
    assert_allclose(b.F, EYE[None], atol=1e-14)
    assert_allclose(b.stress, 0, atol=1e-12)
    assert_allclose(b.v, [[0.1, 0.2, 0.0]])
    assert L.diagnostics.last["particles_with_incomplete_gradient_support"] == 1


@pytest.mark.parametrize("nodes", [2, 3, 7])
def test_quadratic_completeness_on_small_grids(nodes):
    P, G, L, b = mapped(
        interpolation="quadratic",
        G=NS(nx=nodes, ny=nodes, x_m=0.0, x_M=1.0, y_m=0.0, y_M=1.0, ns=1, s=[0.001]),
    )
    assert_allclose(b.N.sum(axis=1), 1)
    assert_allclose(np.einsum("nr,nri->ni", b.N, b.offsets), 0, atol=1e-14)
    assert_allclose(
        np.einsum("nri,nrj->nij", b.offsets, b.G),
        np.broadcast_to(np.diag([1, 1, 0]), (len(b), 3, 3)),
        atol=1e-14,
    )


def test_stressed_particle_cannot_discard_massless_support_forces():
    P, G, L, b = mapped(
        materials=material(n=1, X=[0.0], Y=[0.0]), pressure="compression", initial_pressure=1.0
    )
    with pytest.raises(FloatingPointError, match="massless grid support"):
        step(P, G, L)


def test_elastic_timestep_bound_responds_to_large_deformation():
    P, G, L, b = mapped()
    original = limits(P, G, L)["elastic_wave"]
    b.F[:] = np.diag([0.2, 1.0, 1.0])
    b.J[:] = 0.2
    b.V = b.V0 * b.J
    b.rho = b.m / b.V
    assert limits(P, G, L)["elastic_wave"] < original
