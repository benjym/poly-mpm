"""Continuum references and operator-isolated verification acceptance checks."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from benchmarks.broader_verification import (
    BULK,
    SHEAR,
    RHO,
    airy_stress,
    manufactured_motion,
    elastic_reference,
    hydrostatic,
    large_deformation,
    affine_patch,
    force_balance,
    galilean_pair,
    remap,
    temporal_study,
)


def test_airy_reference_is_self_equilibrated_and_has_free_surface_tractions():
    xy = np.array([[0.13, 0.27, 0], [0.4, 0.63, 0], [0.86, 0.91, 0]])
    divergence = np.zeros((len(xy), 3))
    for axis in (0, 1):
        shift = np.zeros(3)
        shift[axis] = 1e-6
        divergence += (airy_stress(xy + shift)[:, :, axis] - airy_stress(xy - shift)[:, :, axis]) / 2e-6
    assert_allclose(divergence, 0, atol=2e-8)
    xy[:, 1] = 0
    assert_allclose(airy_stress(xy)[:, :, 1], 0, atol=1e-13)
    xy[:, 1] = 1
    assert_allclose(airy_stress(xy)[:, :, 1], 0, atol=1e-13)


def test_manufactured_finite_strain_load_satisfies_reference_momentum_balance():
    X = np.array([0.08, 0.21, 0.39, 0.62, 0.9])
    t, epsilon = 0.17, 1e-6
    exact = manufactured_motion(X, t)
    piola = []
    for shift in (-epsilon, epsilon):
        motion = manufactured_motion(X + shift, t)
        stress, _ = elastic_reference(motion["F"])
        piola.append(np.linalg.det(motion["F"]) * stress[:, 0, 0] / motion["F"][:, 0, 0])
    divergence = (piola[1] - piola[0]) / (2 * epsilon)
    # Independently differentiate prescribed positions in time for acceleration.
    acc = (manufactured_motion(X, t + epsilon)["v"] - manufactured_motion(X, t - epsilon)["v"]) / (
        2 * epsilon
    )
    assert_allclose(divergence + RHO * exact["body"], RHO * acc, rtol=1e-8, atol=2e-6)
    assert_allclose(
        exact["F"][:, 0, 0],
        (manufactured_motion(X + epsilon, t)["x"] - manufactured_motion(X - epsilon, t)["x"]) / (2 * epsilon),
        rtol=1e-9,
    )


def test_linear_hydrostatic_column_remains_at_exact_rest():
    row = hydrostatic(8, 32, basis="linear", transfer="flip")
    for key in (
        "initial_force_rms_over_gravity",
        "initial_bottom_reaction_relative_error",
        "velocity_rms",
        "displacement_rms",
        "stress_relative_l2",
    ):
        assert row[key] < 1e-12
    assert row["relative_mass_error"] == 0
    assert row["incomplete_support"] == 0


def test_quadratic_hydrostatic_quadrature_error_is_measured_and_reduced():
    coarse = hydrostatic(8, 64, ppc=2)
    fine = hydrostatic(8, 64, ppc=4)
    assert coarse["initial_force_rms_over_gravity"] > 1e-3
    assert fine["initial_force_rms_over_gravity"] < 0.26 * coarse["initial_force_rms_over_gravity"]
    assert fine["velocity_rms"] < 0.3 * coarse["velocity_rms"]
    assert (
        fine["initial_bottom_reaction_relative_error"]
        < 0.26 * coarse["initial_bottom_reaction_relative_error"]
    )


@pytest.mark.parametrize("case", ["hydrostatic", "translated_stress", "large_deformation"])
def test_analytical_force_balance_with_resolved_gaussian_quadrature(case):
    row = force_balance(case, 8, quadrature="gauss")
    assert row["acceleration_rms"] < 1e-7


def test_distorted_midpoint_quadrature_defect_reduces_with_particle_sampling():
    coarse = force_balance("large_deformation", 16, ppc=2)
    fine = force_balance("large_deformation", 16, ppc=8)
    assert fine["acceleration_rms"] < coarse["acceleration_rms"] * 0.1
    assert fine["force_resultant_error"] < 1e-10


def test_stressed_body_crosses_cells_and_periodic_seam_under_galilean_boost():
    row = galilean_pair(16, 128)
    assert row["grid_cells_crossed"] > 1
    assert row["velocity_rms"] < 4e-6
    assert row["position_rms"] < 5e-7
    assert row["stress_relative_l2"] < 0.002


def test_large_deformation_full_solver_against_manufactured_reference():
    row = large_deformation(16, 256)
    assert row["displacement_relative_l2"] < 0.02
    assert row["stress_relative_l2"] < 0.03
    assert row["deformation_relative_l2"] < 0.006
    assert 0.5 < row["minimum_J"] < 0.7
    assert 1.3 < row["maximum_J"] < 1.5
    assert row["relative_mass_error"] == 0
    assert row["maximum_momentum_residual"] < 1e-10
    assert row["all_steps_within_stability_estimate"]
    assert row["incomplete_support"] == 0


def test_prescribed_affine_patch_separates_gradient_integration_from_advection():
    rows = [affine_patch(n) for n in (8, 16, 32)]
    for coarse, fine in zip(rows, rows[1:]):
        # Smooth midpoint L yields second-order F; particle positions use Euler.
        assert fine["deformation_error"] < coarse["deformation_error"] * 0.26
        assert fine["stress_relative_error"] < coarse["stress_relative_error"] * 0.26
        assert fine["position_error"] < coarse["position_error"] * 0.55
    for row in rows:
        assert row["determinant_error"] < 2e-14
        assert row["density_error"] < 1e-9


def test_fixed_cloud_remapping_isolates_filtering_without_motion_or_forces():
    rows = {t: remap(16, 8, transfer=t) for t in ("flip", "pic", "apic")}
    assert rows["flip"]["velocity_relative_l2"] < 1e-14
    assert rows["apic"]["velocity_relative_l2"] < 0.01
    assert rows["pic"]["velocity_relative_l2"] > 0.2
    assert abs(rows["apic"]["relative_energy_change"]) < abs(rows["pic"]["relative_energy_change"]) * 0.03
    for row in rows.values():
        assert row["momentum_error"] < 1e-12
        assert row["maximum_position_change"] == 0
        assert row["relative_energy_change"] <= 1e-14


def test_temporal_comparison_keeps_grid_and_particle_cloud_fixed():
    rows = temporal_study(hydrostatic, 8, steps=(32, 64, 128))
    assert {r["cells"] for r in rows} == {8}
    assert {r["final_time"] for r in rows} == {0.05}
    assert rows[1]["difference_to_finest"]["v"] < rows[0]["difference_to_finest"]["v"]
    assert rows[-1]["difference_to_finest"] == dict(x=0.0, v=0.0, F=0.0, stress=0.0)
