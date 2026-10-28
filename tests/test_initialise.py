from copy import deepcopy
from types import SimpleNamespace as NS
import numpy as np
import pytest
from numpy.testing import assert_allclose
from conftest import state, material
from initialise import prepare
from MPM import run


@pytest.mark.parametrize('mode', [None, 'lithostatic', 'compression', 'non-zero', 'temperature'])
def test_pressure_initialisation(mode):
    P, G, L = state(pressure=mode, initial_pressure=7., max_g=-10.)
    b = L.materials[0]
    assert_allclose(b.pressure, np.trace(b.stress, axis1=1, axis2=2)/3)
    assert_allclose(b.dev_stress, b.stress-b.pressure[:, None, None]*np.eye(3))
    if mode == 'compression':
        assert_allclose(b.pressure, -7)
    if mode == 'lithostatic':
        assert_allclose(b.stress[:, 1, 1], -10*b.rho*(1-b.x[:, 1]))
    if mode in ('compression', 'lithostatic'):
        reconstructed = b.material.K*np.trace(b.strain, axis1=1, axis2=2)[:, None, None]*np.eye(3)
        reconstructed += 2*b.material.G*(b.strain-np.trace(b.strain, axis1=1, axis2=2)[:, None, None]*np.eye(3)/3)
        assert_allclose(reconstructed, b.stress)


@pytest.mark.parametrize('law', ['viscous', 'pouliquen'])
def test_steady_flow_initialisation(law):
    P, G, L = state(law, pressure='lithostatic', initial_flow='steady', max_g=-10., theta=.5)
    b = L.materials[0]
    assert np.isfinite(b.v).all() and np.isfinite(b.pk).all()
    assert np.any(b.v[:, 0] != 0)
    assert_allclose(b.v[:, 1:], 0)


def test_poiseuille_initialisation():
    P, G, L = state('pouliquen', initial_flow='steady_poiseulle', max_g=-.01,
                    pressure='compression', initial_pressure=100.)
    b = L.materials[0]
    assert np.isfinite(b.v).all()
    assert np.any(b.v[:, 0] != 0)


def test_aliases_and_conflicting_segregation_flags():
    assert state(segregate_grid=True)[0].segregate
    with pytest.raises(ValueError, match='disagree'):
        state(segregate=True, segregate_grid=False)


def test_dynamic_limit_across_materials():
    materials = [material(), material('viscous')]
    materials[0].critical_time = lambda P: .01
    materials[1].critical_time = lambda P: .0001
    P, G, L = state(materials=materials)
    P.update_timestep(P, G)
    assert_allclose(P.dt, P.CFL*.0001)


def test_final_time_no_extra_step_or_overshoot():
    P, G, L = state(t_f=2.5e-5)
    run(P, G, L, output=False)
    assert_allclose(P.t, P.t_f)
    assert P.tstep == 3


def test_prescribed_particle_velocities_survive_removal():
    velocity = np.array([[1, 0, 0], [2, 0, 0], [3, 0, 0], [4, 0, 0]])
    P, G, L = state(materials=material('rigid', initial_v=velocity))
    b = L.materials[0]
    b.keep([True, False, True, False])
    before = b.x.copy()
    b.move_material_points(P, G)
    assert_allclose(b.x, before+P.dt*velocity[[0, 2]])


def test_inlet_initialises_fresh_stress_history():
    P, G, L = state()
    b = L.materials[0]
    b.stress[:] = 999
    template = L.templates[0]
    template.inlet_rate, template.gap, template.yp = .1, [.01, .01], [.2]
    P.B.inlet_right, P.t = True, .1
    L.process_boundaries(P, G)
    assert_allclose(b.stress[-1], 0)
    assert_allclose(b.x[-1], [.99, .2, 0])


def test_grid_traction_and_conveyor():
    P, G, L = state(B=NS(vertical_force=True, horizontal_force=True, conveyor=True), max_q=4.)
    G.m[:] = 2
    G.BCs(P)
    assert_allclose(G.fe[0, :2], [2*4*2/1000/G.dx, 2*4*2/1000/G.dy])
    P.v_0 = .25
    G.update_momentum(P)
    assert_allclose(G.q[G.boundary_conveyor, 0], .5)


def test_declared_grain_bounds_are_preserved():
    grid = NS(nx=5, ny=5, x_m=0., x_M=1., y_m=0., y_M=1., ns=2,
              s=[.002, .003], s_m=.001, s_M=.004)
    P, G, L = state(ns=2, G=grid)
    assert P.G.s_m == .001 and P.G.s_M == .004
