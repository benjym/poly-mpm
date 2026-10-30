from copy import deepcopy
from types import SimpleNamespace as NS
import numpy as np
import pytest
from numpy.testing import assert_allclose
from conftest import state, material
from initialise import prepare
from MPM import run


@pytest.mark.parametrize("mode", [None, "lithostatic", "compression"])
def test_pressure_initialisation(mode):
    P, G, L = state(initial_stress=mode, initial_pressure=7.0, gravity=-10.0)
    b = L.materials[0]
    assert_allclose(b.pressure, np.trace(b.stress, axis1=1, axis2=2) / 3)
    assert_allclose(b.dev_stress, b.stress - b.pressure[:, None, None] * np.eye(3))
    if mode == "compression":
        assert_allclose(b.pressure, -7)
    if mode == "lithostatic":
        assert_allclose(b.stress[:, 1, 1], -10 * b.rho * (1 - b.x[:, 1]))
    if mode in ("compression", "lithostatic"):
        reconstructed = b.material.K * np.trace(b.strain, axis1=1, axis2=2)[:, None, None] * np.eye(3)
        reconstructed += (
            2
            * b.material.G
            * (b.strain - np.trace(b.strain, axis1=1, axis2=2)[:, None, None] * np.eye(3) / 3)
        )
        assert_allclose(reconstructed, b.stress)


@pytest.mark.parametrize("law", ["viscous", "pouliquen"])
def test_steady_flow_initialisation(law):
    P, G, L = state(
        law, initial_stress="lithostatic", initial_flow="steady", gravity=-10.0, gravity_angle=0.5
    )
    b = L.materials[0]
    assert np.isfinite(b.v).all() and np.isfinite(b.pk).all()
    assert np.any(b.v[:, 0] != 0)
    assert_allclose(b.v[:, 1:], 0)


def test_poiseuille_initialisation():
    P, G, L = state(
        "pouliquen",
        initial_flow="steady_poiseuille",
        gravity=-0.01,
        initial_stress="compression",
        initial_pressure=100.0,
    )
    b = L.materials[0]
    assert np.isfinite(b.v).all()
    assert np.any(b.v[:, 0] != 0)


def test_dynamic_limit_across_materials():
    materials = [material(), material("viscous")]
    materials[0].critical_time = lambda P: 0.01
    materials[1].critical_time = lambda P: 0.0001
    P, G, L = state(materials=materials)
    P.update_timestep(P, G)
    assert_allclose(P.dt, P.cfl * 0.0001)


def test_final_time_no_extra_step_or_overshoot():
    P, G, L = state(t_f=2.5e-5)
    run(P, G, L, output=False)
    assert_allclose(P.t, P.t_f)
    assert P.tstep == 3


def test_prescribed_particle_velocities_survive_removal():
    velocity = np.array([[1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]])
    P, G, L = state(materials=material("rigid", initial_v=velocity))
    b = L.materials[0]
    b.keep([True, False, True, False])
    before = b.x.copy()
    b.move_material_points(P, G)
    assert_allclose(b.x, before + P.dt * velocity[[0, 2]])


def test_inlet_initialises_fresh_stress_history():
    S = material(inlet_rate=0.1, inlet_gap=[0.01, 0.01], inlet_y=[0.2])
    P, G, L = state(materials=S, B=NS(inlets={"right"}))
    b = L.materials[0]
    b.stress[:] = 999
    P.t = 0.1
    L.process_boundaries(P, G)
    assert_allclose(b.stress[-1], 0)
    assert_allclose(b.x[-1], [0.99, 0.2, 0])


def test_grid_traction_and_conveyor():
    P, G, L = state(B=NS(traction_vertical=4.0, traction_horizontal=4.0, conveyor_speed=0.25))
    G.m[:] = 2
    G.BCs(P)
    assert_allclose(G.fe[0, :2], [2 * 4 * 2 / 1000 / G.dx, 2 * 4 * 2 / 1000 / G.dy])
    G.update_momentum(P)
    assert_allclose(G.q[G.boundary_conveyor, 0], 0.5)
