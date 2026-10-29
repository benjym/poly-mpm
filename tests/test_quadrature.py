"""Finite-domain integration, conservation and cell-crossing checks."""

from types import SimpleNamespace as NS
import numpy as np
import pytest
from numpy.testing import assert_allclose
from conftest import material, state
from quadrature import domain_stencil, _polygon_stencil
from interpolation import moment, affine_matrix
from diagnostics import inventory
from MPM import run
from benchmarks.broader_verification import force_balance, large_deformation


def setup(basis="quadratic", periodic=False, F=None, xy=None):
    xy = np.array([[0.25, 0.25], [0.51, 0.43], [0.8, 0.8]]) if xy is None else np.asarray(xy)
    sizes = np.broadcast_to([0.08, 0.06], (len(xy), 2)).copy()
    P, G, L = state(
        materials=material(n=len(xy), X=xy[:, 0], Y=xy[:, 1], A=np.prod(sizes, axis=1), domain_size=sizes),
        quadrature="domain",
        interpolation=basis,
        B=NS(cyclic_lr=periodic),
        transfer="apic" if basis == "quadratic" else "flip",
    )
    b = L.materials[0]
    if F is not None:
        b.F[:] = F
        b.J[:] = np.linalg.det(b.F)
        b.V = b.V0 * b.J
        b.rho = b.m / b.V
    b.get_reference_nodes(P, G)
    b.get_basis_functions(P, G)
    return P, G, L, b


@pytest.mark.parametrize("basis", ["linear", "quadratic"])
@pytest.mark.parametrize("periodic", [False, True])
@pytest.mark.parametrize(
    "F", [np.eye(3), np.diag([0.6, 1.4, 1]), np.array([[0.8, -0.5, 0], [0.4, 1.1, 0], [0, 0, 1]])]
)
def test_domain_basis_complete_under_compression_shear_rotation_and_periodic_seams(basis, periodic, F):
    xy = [[0.01, 0.4], [0.99, 0.6], [0.5, 0.5]] if periodic else [[0.06, 0.08], [0.5, 0.5], [0.94, 0.92]]
    P, G, L, b = setup(basis, periodic, F, xy)
    assert np.all(b.N >= 0)
    assert_allclose(b.N.sum(axis=1), 1, atol=2e-13)
    assert_allclose(b.G.sum(axis=1), 0, atol=2e-12)
    assert_allclose(np.einsum("nr,nri->ni", b.N, b.offsets), 0, atol=1e-13)
    assert_allclose(
        np.einsum("nri,nrj->nij", b.offsets, b.G),
        np.broadcast_to(np.diag([1, 1, 0]), (len(b), 3, 3)),
        atol=3e-12,
    )
    L.recover_position(P, G)


@pytest.mark.parametrize("basis", ["linear", "quadratic"])
def test_integrated_gradient_is_position_derivative_with_domain_held_fixed(basis):
    F = np.array([[0.9, -0.2, 0], [0.3, 1.1, 0], [0, 0, 1]])
    P, G, L, b = setup(basis, F=F)
    delta = 2e-7
    for axis in (0, 1):
        numerical = np.zeros((len(b), len(G.m)))
        expected = numerical.copy()
        for sign in (-1, 1):
            x = b.x.copy()
            x[:, axis] += sign * delta
            n, w, _, _ = domain_stencil(x, b.F, b.domain_halfsize, P, G)
            for p in range(len(b)):
                np.add.at(numerical[p], n[p], sign * w[p] / (2 * delta))
        for p in range(len(b)):
            np.add.at(expected[p], b.n[p], b.G[p, :, axis])
        assert_allclose(numerical, expected, atol=8e-9)


@pytest.mark.parametrize("basis", ["linear", "quadratic"])
def test_axis_integration_agrees_with_independent_clipped_polygon_rule(basis):
    P, G, L, b = setup(basis, F=np.diag([0.6, 1.4, 1]))
    for p in range(len(b)):
        n, w, g, _ = _polygon_stencil(b.x[p], b.F[p], b.domain_halfsize[p], P, G)
        a = np.zeros((len(G.m), 4))
        c = a.copy()
        np.add.at(a[:, 0], b.n[p], b.N[p])
        np.add.at(a[:, 1:], b.n[p], b.G[p])
        np.add.at(c[:, 0], n, w)
        np.add.at(c[:, 1:], n, g)
        assert_allclose(a, c, atol=2e-12)


def test_tiny_domain_limit_does_not_lose_partition_or_affine_completeness():
    P, G, L, b = setup()
    b.F[:] = np.diag([1e-10, 1e-10, 1])
    b.get_basis_functions(P, G)
    assert_allclose(b.N.sum(axis=1), 1, atol=1e-13)
    assert_allclose(
        np.einsum("nri,nrj->nij", b.offsets, b.G),
        np.broadcast_to(np.diag([1, 1, 0]), (len(b), 3, 3)),
        atol=1e-12,
    )


def test_variable_domain_sizes_and_empty_batch_are_supported():
    sizes = np.array([[0.01, 0.02], [0.2, 0.4], [0.1, 0.1]])
    P, G, L = state(materials=material(n=3, A=np.prod(sizes, axis=1), domain_size=sizes), quadrature="domain")
    b = L.materials[0]
    assert_allclose(b.N.sum(axis=1), 1)
    b.keep(np.zeros(len(b), dtype=bool))
    b.get_basis_functions(P, G)
    assert b.N.shape == (0, 4)
    run(P, G, L, max_steps=1, output=False)


def test_internal_force_and_strain_are_power_conjugate_and_torque_free():
    P, G, L, b = setup(F=np.array([[1, 0.3, 0], [-0.2, 0.9, 0], [0, 0, 1]]))
    rng = np.random.default_rng(3)
    s = rng.normal(size=b.stress.shape)
    b.stress[:] = s + s.swapaxes(1, 2)
    b.stress[:, :2, 2] = 0
    b.stress[:, 2, :2] = 0
    L.get_nodal_mass_momentum(P, G)
    velocity = rng.normal(size=G.q.shape)
    velocity[:, 2] = 0
    G.q[:] = G.m[:, None] * velocity
    b.get_nodal_forces(P, G)
    gradient = np.einsum("nri,nrj->nij", velocity[b.n], b.G)
    assert_allclose(
        np.sum(velocity * G.fi), np.sum(b.V * np.einsum("nij,nij->n", b.stress, gradient)), atol=1e-13
    )
    assert_allclose(G.fi.sum(axis=0), 0, atol=1e-13)
    xyz = np.column_stack((G.X, G.Y, np.zeros_like(G.X)))
    assert_allclose(np.cross(xyz, G.fi).sum(axis=0), 0, atol=1e-13)


def test_matched_domain_apic_transfers_reproduce_affine_motion_and_conserve_momentum():
    P, G, L, b = setup(F=np.array([[1, 0.2, 0], [-0.1, 0.9, 0], [0, 0, 1]]))
    A = np.array([[0.1, -0.3, 0], [0.2, 0.05, 0], [0, 0, 0]])
    b.v[:] = b.x @ A.T + [0.4, -0.1, 0]
    b.m *= np.array([0.5, 1.0, 2.0])
    b.affine_momentum = A @ moment(b.N, b.offsets)
    before = inventory(P, L)
    oldv = b.v.copy()
    oldB = b.affine_momentum.copy()
    L.get_nodal_mass_momentum(P, G)
    b.move_material_points(P, G)
    after = inventory(P, L)
    assert_allclose(b.v, oldv, atol=2e-14)
    assert_allclose(b.affine_momentum, oldB, atol=2e-14)
    for key in ("momentum", "angular_momentum"):
        assert_allclose(after[key], before[key], atol=1e-12)
    assert_allclose(
        after["kinetic"] + after["affine_kinetic"], before["kinetic"] + before["affine_kinetic"], atol=1e-12
    )


def test_rotated_domains_preserve_uniform_translation_through_full_steps():
    angle = 0.6
    c, s = np.cos(angle), np.sin(angle)
    P, G, L, b = setup(F=np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]))
    b.v[:] = [0.1, 0.05, 0]
    run(P, G, L, max_steps=5, output=False)
    assert_allclose(b.v, np.broadcast_to([0.1, 0.05, 0], b.v.shape), atol=1e-13)
    assert_allclose(b.stress, 0, atol=1e-11)
    assert_allclose(b.J, 1, atol=1e-13)


@pytest.mark.parametrize(
    "overrides,message",
    [
        ({"quadrature": "unknown"}, "quadrature must"),
        ({"mechanics": "legacy", "quadrature": "domain"}, "updated mechanics"),
        ({"quadrature": "domain", "domain_max_support": 3}, "domain_max_support"),
    ],
)
def test_invalid_quadrature_options_rejected(overrides, message):
    with pytest.raises(ValueError, match=message):
        state(**overrides)


@pytest.mark.parametrize(
    "domain_size,A,message",
    [(None, 0.01, "domain_size"), ([0.1, -0.1], 0.01, "positive"), ([0.1, 0.2], 0.01, "match")],
)
def test_invalid_material_domain_rejected(domain_size, A, message):
    m = material(A=A)
    if domain_size is not None:
        m.domain_size = domain_size
    with pytest.raises(ValueError, match=message):
        state(materials=m, quadrature="domain")


def test_domain_geometry_and_wall_extrapolation_fail_explicitly_when_invalid():
    P, G, L, b = setup()
    b.F[:, 2, 0] = 0.1
    with pytest.raises(ValueError, match="planar"):
        b.get_basis_functions(P, G)
    b.F[:] = np.eye(3)
    b.x[:, 0] = 0
    b.domain_halfsize[:, 0] = 0.2
    with pytest.raises(ValueError, match="negative integrated weights"):
        b.get_basis_functions(P, G)


def test_domain_geometry_travels_with_particle_rows_and_fresh_inlets():
    P, G, L, b = setup()
    b.keep(np.array([True, False, True]))
    assert b.domain_halfsize.shape == (2, 2)
    b.F[:] = np.diag([0.9, 1.1, 1])
    template = L.templates[0]
    template.domain_size = [0.08, 0.06]
    template.A = 0.08 * 0.06
    template.inlet_rate = 0.1
    template.gap = [0.1, 0.1]
    template.yp = [0.3]
    template.rho = 1000.0
    P.B.inlet_right = True
    P.t = 0.1
    L.process_boundaries(P, G)
    assert_allclose(b.domain_halfsize[-1], [0.04, 0.03])
    assert_allclose(b.F[-1], np.eye(3))
    b.get_basis_functions(P, G)


def test_domain_integration_reduces_distorted_cloud_force_error_across_grid_crossings():
    for offset in (0.0, 0.25 / 16, 0.5 / 16):
        point = force_balance("large_deformation", 16, offset=offset)
        domain = force_balance("large_deformation", 16, integration="domain", offset=offset)
        assert domain["acceleration_rms"] < point["acceleration_rms"] * 0.15
        assert domain["force_resultant_error"] < 1e-10


def test_domain_finite_strain_refinement_does_not_have_point_quadrature_error_growth():
    rows = [large_deformation(n, 256, quadrature="domain") for n in (16, 32)]
    assert rows[1]["stress_relative_l2"] < rows[0]["stress_relative_l2"] * 0.5
    assert rows[1]["stress_relative_l2"] < 0.01
    for r in rows:
        assert r["relative_mass_error"] == 0
        assert r["incomplete_support"] == 0
        assert r["all_steps_within_stability_estimate"]


def test_nonaffine_domain_apic_roundtrip_filters_energy_without_losing_angular_momentum():
    P, G, L, b = setup(F=np.array([[1, 0.2, 0], [-0.1, 0.9, 0], [0, 0, 1]]))
    rng = np.random.default_rng(38)
    b.v[:] = rng.normal(size=b.v.shape)
    b.v[:, 2] = 0
    b.affine_momentum[:] = rng.normal(size=b.affine_momentum.shape) * 0.01
    b.affine_momentum[:, :2, 2] = 0
    b.affine_momentum[:, 2, :] = 0
    before = inventory(P, L)
    L.get_nodal_mass_momentum(P, G)
    grid_energy = 0.5 * np.sum(np.sum(G.q[G.valid_nodes] ** 2, axis=1) / G.m[G.valid_nodes])
    b.move_material_points(P, G)
    after = inventory(P, L)
    assert after["kinetic"] + after["affine_kinetic"] <= grid_energy + 1e-12
    assert grid_energy <= before["kinetic"] + before["affine_kinetic"] + 1e-12
    for name in ("momentum", "angular_momentum"):
        assert_allclose(after[name], before[name], atol=1e-12)


def test_domain_quadrature_rejects_out_of_plane_grid_motion_and_excessive_support():
    P, G, L, b = setup()
    L.get_nodal_mass_momentum(P, G)
    G.q[:, 2] = G.m
    with pytest.raises(ValueError, match="planar grid"):
        b.update_stress_strain(P, G)
    with pytest.raises(ValueError, match="planar grid"):
        b.move_material_points(P, G)
    P.domain_max_support = 4
    with pytest.raises(ValueError, match="domain_max_support"):
        b.get_basis_functions(P, G)


def test_full_domain_column_example_records_integration_mode_and_support():
    from initialise import get_parameters

    P, G, L = get_parameters(["domain_column"])
    run(P, G, L, output=False)
    assert P.t == P.t_f
    assert L.diagnostics.last["quadrature"] == "domain"
    assert L.diagnostics.last["maximum_stencil_entries"] > 9
    assert L.diagnostics.last["mass_error"] == 0
    assert L.diagnostics.last["particles_with_incomplete_support"] == 0


def test_mixed_aligned_and_rotated_domains_share_a_complete_padded_stencil():
    P, G, L, b = setup()
    b.F[1, :2, :2] = [[0.9, 0.3], [-0.2, 1.1]]
    b.get_basis_functions(P, G)
    assert_allclose(b.N.sum(axis=1), 1, atol=1e-13)
    assert_allclose(np.einsum("nr,nri->ni", b.N, b.offsets), 0, atol=1e-13)
    assert_allclose(
        np.einsum("nri,nrj->nij", b.offsets, b.G),
        np.broadcast_to(np.diag([1, 1, 0]), (len(b), 3, 3)),
        atol=1e-12,
    )
    b.F[1, 0, 0] = 10000
    with pytest.raises(ValueError, match="128 grid"):
        b.get_basis_functions(P, G)


def test_periodic_domain_translation_crosses_the_seam_without_strain_or_mass_loss():
    P, G, L, b = setup(periodic=True, xy=[[0.999, 0.4], [0.001, 0.6], [0.3, 0.5]])
    b.v[:] = [1, 0, 0]
    run(P, G, L, max_steps=200, output=False)
    assert b.x[0, 0] < 0.01
    assert_allclose(b.v, np.broadcast_to([1, 0, 0], b.v.shape), atol=1e-12)
    assert_allclose(b.F, np.broadcast_to(np.eye(3), b.F.shape), atol=1e-12)
    assert L.diagnostics.last["mass_error"] == 0


def test_moving_stress_benchmark_uses_the_selected_domain_integrator():
    from benchmarks.broader_verification import translated_stress

    row = translated_stress(4, 16, quadrature="domain")
    assert row["quadrature"] == "domain"
    assert row["stress_relative_l2"] < 0.02
    assert row["relative_mass_error"] == 0


@pytest.mark.parametrize("basis", ["linear", "quadratic"])
@pytest.mark.parametrize("periodic", [False, True])
def test_batched_polygon_integration_matches_the_reference_path(basis, periodic):
    rng = np.random.default_rng(7)
    count = 40
    xy = rng.uniform(0.3, 0.7, size=(count, 2))
    if periodic:
        xy[:5, 0] = rng.uniform(-0.02, 0.02, 5) % 1.0
    angle, shear, stretch = (
        rng.uniform(-0.5, 0.5, count),
        rng.uniform(-0.3, 0.3, count),
        rng.uniform(0.7, 1.4, count),
    )
    F = np.tile(np.eye(3), (count, 1, 1))
    F[:, 0, 0], F[:, 0, 1] = stretch * np.cos(angle), -np.sin(angle) + shear
    F[:, 1, 0], F[:, 1, 1] = stretch * np.sin(angle), np.cos(angle)
    P, G, L, b = setup(basis, periodic, xy=xy)
    P.domain_polygon = "reference"
    reference = domain_stencil(b.x, F, b.domain_halfsize, P, G)
    P.domain_polygon = "batched"
    batched = domain_stencil(b.x, F, b.domain_halfsize, P, G)
    for i in range(count):
        a = {int(n): (w, g) for n, w, g in zip(*(r[i] for r in reference[:3])) if w or np.any(g)}
        c = {int(n): (w, g) for n, w, g in zip(*(r[i] for r in batched[:3])) if w or np.any(g)}
        assert a.keys() == c.keys()
        for node in a:
            assert_allclose(c[node][0], a[node][0], rtol=0, atol=1e-13)
            assert_allclose(c[node][1], a[node][1], rtol=0, atol=1e-13 / G.dx)
