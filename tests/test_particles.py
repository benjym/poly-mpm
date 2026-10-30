from copy import deepcopy
from types import SimpleNamespace as NS
import numpy as np
import pytest
from numpy.testing import assert_allclose
from conftest import state, material
from constit import get_law
from particles import Particles
from transfers import gather, divide
from MPM import run, step


def mapping(P, G, L):
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    return L.materials[0]


def test_shape_functions_and_reconstruction(basic):
    P, G, L = basic
    b = mapping(*basic)
    assert_allclose(b.N.sum(axis=1), 1)
    assert_allclose(b.G.sum(axis=1), 0, atol=1e-12)
    L.recover_position(P, G)


def test_repeated_indices_unequal_mass_and_composition():
    P, G, L = state(
        ns=3,
        materials=material(
            n=3,
            rho=[1000, 1200, 1500],
            A=[0.01, 0.02, 0.03],
            phi=np.eye(3),
            initial_v=[[1, 2, 3], [2, 1, 0], [-1, 3, 2]],
        ),
    )
    b = mapping(P, G, L)
    L.get_nodal_mass_momentum(P, G)
    assert_allclose(G.m.sum(), b.m.sum())
    assert_allclose(G.q.sum(axis=0), (b.m[:, None] * b.v).sum(axis=0))
    assert_allclose(G.phim.sum(axis=0), b.m)
    assert_allclose(b.m, [10, 24, 45])
    b.pressure[:] = [1, 2, 3]
    b.deposit_diagnostics(P, G)
    assert_allclose(G.pressure.sum(), b.m @ b.pressure)


@pytest.mark.parametrize(
    "gradient",
    [
        np.zeros((3, 3)),
        np.array([[0, -2, 0], [2, 0, 0], [0, 0, 0]]),
        np.array([[2, 3, 0], [-1, 4, 0], [0, 0, 0]]),
    ],
)
def test_affine_velocity_gradient(basic, gradient):
    P, G, L = basic
    b = mapping(*basic)
    G.m[:] = 1
    G.q[:] = np.column_stack((G.X, G.Y, np.zeros_like(G.X))) @ gradient.T + [1, 2, 0]
    G.get_valid_nodes(P)
    b.update_stress_strain(P, G)
    assert_allclose(b.L, np.broadcast_to(gradient, b.L.shape), atol=1e-12)
    assert_allclose(
        b.dstrain, np.broadcast_to(0.5 * (gradient + gradient.T) * P.dt, b.dstrain.shape), atol=1e-12
    )
    assert_allclose(b.L_T, b.L.swapaxes(1, 2))
    if not np.allclose(gradient, gradient.T):
        assert not np.allclose(b.L, b.L_T)


def test_force_direction_and_total_body_force(basic):
    P, G, L = basic
    b = mapping(*basic)
    b.stress[:] = [[2, 3, 0], [3, -4, 0], [0, 0, 0]]
    P.g = -9.81
    b.get_nodal_forces(P, G)
    expected = np.zeros_like(G.fi)
    for i in range(len(b)):
        for r, n in enumerate(b.n[i]):
            expected[n] += b.V[i] * b.stress[i] @ b.G[i, r]
    assert_allclose(G.fi, expected, atol=1e-12)
    assert_allclose(G.fi.sum(axis=0), 0, atol=1e-12)
    assert_allclose(G.fe.sum(axis=0), [0, P.g * b.m.sum(), 0])


def test_pouliquen_zero_rate_and_viscosity_cap():
    P, G, L = state("pouliquen")
    b = L.materials[0]
    b.pressure[:] = [-1, 0, 1, -2]
    b.law(b, b.material, P)
    assert np.isfinite(b.dstress).all()
    assert np.all(b.eta_limited <= b.material.eta_max)
    assert np.all(b.pressure <= 0)


def test_domain_edges_and_periodic_multiple_crossings(basic):
    P, G, L = basic
    b = L.materials[0]
    b.x[:, :2] = [[0, 0], [1, 1], [0, 1], [1, 0]]
    mapping(P, G, L)
    assert np.all(b.n >= 0) and np.all(b.n < len(G.m))
    b.x[0, 0] = -0.001
    with pytest.raises(ValueError, match="outside domain"):
        b.get_reference_nodes(P, G)
    P, G, L = state(B=NS(periodic=True))
    b = L.materials[0]
    b.x[:, 0] = [-3.1, 4.2, 1, 0]
    b.get_reference_nodes(P, G)
    assert_allclose(b.x[:, 0], [0.9, 0.2, 0, 0])


def test_periodic_mass_validity_after_merge():
    P, G, L = state(materials=material(n=1, X=[0.99], Y=[0.3]), B=NS(periodic=True))
    mapping(P, G, L)
    L.get_nodal_mass_momentum(P, G)
    assert_allclose(G.m[:: P.G.nx], G.m[P.G.nx - 1 :: P.G.nx])
    assert np.array_equal(G.valid_nodes, G.m > P.mass_tolerance)
    # Duplicated seam endpoints count once in physical totals.
    mass = G.m.reshape(P.G.ny, P.G.nx)
    assert_allclose(mass[:, :-1].sum(), L.materials[0].m.sum())


def test_empty_and_mixed_materials():
    P, G, L = state(
        materials=[
            material(n=0),
            material(Y=np.linspace(0.26, 0.46, 4)),
            material("rigid", initial_v=[0.1, 0, 0]),
        ]
    )
    rigid_x = L.materials[2].x.copy()
    run(P, G, L, max_steps=5, output=False)
    assert L.materials[0].x.shape == (0, 3)
    assert_allclose(L.materials[2].x, rigid_x + np.array([0.1, 0, 0]) * P.dt * 5)
    for b in L.materials:
        b.check_finite(P)


def test_initialisation_shapes_independent_storage():
    velocity = np.array([1.0, 2.0, 3.0])
    P, G, L = state(ns=2, materials=material(initial_v=velocity, phi=[0.2, 0.8]))
    b = L.materials[0]
    assert b.m.shape == (4,) and b.phi.shape == (4, 2)
    b.v[0, 0] = 99
    assert velocity[0] == 1 and b.v[1, 0] == 1


@pytest.mark.parametrize(
    "overrides,match",
    [
        ({"initial_v": [1, 2]}, "initial_v"),
        ({"rho": -1}, "positive"),
        ({"A": 0}, "positive"),
        ({"X": [np.nan] * 4}, "finite"),
        ({"phi": [-0.2, 1.2]}, "fractions"),
    ],
)
def test_invalid_initial_state(overrides, match):
    with pytest.raises(ValueError, match=match):
        state(ns=2, materials=material(**overrides))


@pytest.mark.parametrize("name", ["elastic_time", "ken_kamrin", "ken_simple", "viscous_size", "missing"])
def test_unsupported_models_fail_early(name):
    with pytest.raises(ValueError, match="Unsupported material law"):
        state(name)


def test_remove_and_append_all_fields(basic):
    P, G, L = basic
    b = L.materials[0]
    b.keep([True, False, True, False])
    assert b.material.n == 2
    template = deepcopy(L.templates[0])
    template.n, template.X, template.Y = 1, [0.5], [0.5]
    addition = Particles(P, G, template)
    b.append(addition)
    assert b.material.n == 3
    for key in b.persistent_fields:
        assert len(getattr(b, key)) == 3
    mapping(P, G, L)


def test_outlet_and_inlet():
    S = material(inlet_rate=0.1, inlet_gap=[0.01, 0.01], inlet_y=[0.2, 0.4])
    P, G, L = state(materials=S, B=NS(outlets={"bottom"}, inlets={"right"}))
    b = L.materials[0]
    b.x[:2, 1] = -0.1
    L.process_boundaries(P, G)
    assert len(b) == 2
    P.t = 0.1
    L.process_boundaries(P, G)
    assert len(b) == 4
    assert_allclose(b.x[-2:, :2], [[0.99, 0.2], [0.99, 0.4]])


def test_friction_boundaries_use_masks_without_mutation():
    P, G, L = state(B=NS(walls={"bottom", "left"}, wall_contact="coulomb", wall_friction=0.5))
    before = G.boundary_h.copy()
    G.fe[1] = [10, -4, 0]
    G.fe[P.G.nx + 1] = [10, -4, 0]
    G.update_momentum(P)
    assert_allclose(G.q_dot[1], [8, 0, 0])
    assert_allclose(G.q_dot[P.G.nx + 1], [10, -4, 0])
    assert np.array_equal(G.boundary_h, before)


def test_dynamic_timestep_and_callback_alias(basic):
    P, G, L = basic
    old = P.dt
    P.update_timestep(P, G)
    assert P.dt == old
    G.m[:] = 1
    G.q[:, 0] = 2
    P.update_timestep(P, G)
    assert_allclose(P.dt, P.cfl * 0.1 * min(G.dx, G.dy) / 2)
    P.S.append(NS(critical_time=lambda P: 1e-7))
    P.update_timestep(P, G)
    assert_allclose(P.dt, P.cfl * 1e-7)


def test_bad_timestep_and_finite_error_context(basic):
    P, G, L = basic
    P.dt = 0
    with pytest.raises(ValueError, match="dt"):
        step(P, G, L)
    L.materials[0].v[1, 0] = np.nan
    with pytest.raises(FloatingPointError, match="material 0, timestep 0, v, particles \\[1\\]"):
        L.materials[0].check_finite(P)


def test_empty_material_with_empty_property_arrays():
    P, G, L = state(materials=material(n=0, rho=[], A=[]))
    run(P, G, L, max_steps=3, output=False)
    assert len(L.materials[0]) == 0
    assert_allclose(G.m, 0)
