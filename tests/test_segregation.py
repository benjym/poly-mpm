from types import SimpleNamespace as NS
import numpy as np
import pytest
from numpy.testing import assert_allclose
from conftest import state, material
from integrators import increment_grainsize, KT
from MPM import run


@pytest.mark.parametrize('ns', [1, 2, 3, 8])
@pytest.mark.parametrize('periodic', [False, True])
def test_uniform_zero_forcing(ns, periodic):
    P, G, L = state(ns=ns, segregate=True, B=NS(cyclic_lr=periodic))
    G.m[:] = 1
    G.get_valid_nodes(P)
    G.phi[:] = 1/ns
    G.s_bar[:] = G.phi@P.G.s
    G.shear_rate = np.zeros_like(G.m)
    G.I[:] = 1
    assert_allclose(increment_grainsize(P, G), 0, atol=1e-12)
    before = L.materials[0].phi.copy()
    run(P, G, L, max_steps=20, output=False)
    assert_allclose(L.materials[0].phi, before, atol=1e-12)


@pytest.mark.parametrize('ns', [2, 3, 8])
@pytest.mark.parametrize('ax', [0, 1])
@pytest.mark.parametrize('periodic', [False, True])
def test_conservative_closed_flux(ns, ax, periodic):
    P, G, L = state(ns=ns, segregate=True, c=.001, B=NS(cyclic_lr=periodic))
    G.m[:] = 1
    G.get_valid_nodes(P)
    rng = np.random.default_rng(42)
    G.phi[:] = rng.uniform(.1, 1, size=G.phi.shape)
    G.phi /= G.phi.sum(axis=1, keepdims=True)
    G.I[:] = 1
    G.s_bar = G.phi@P.G.s
    G.shear_rate = np.ones_like(G.m)
    G.grad_pk[:, ax] = .5
    if periodic:
        G.phi[P.G.nx-1::P.G.nx] = G.phi[::P.G.nx]
        G.s_bar = G.phi@P.G.s
    increment = KT(P, G, ax).reshape(P.G.ny, P.G.nx, ns)
    physical = increment[:, :-1] if periodic else increment
    assert_allclose(physical.sum(axis=(0, 1)), 0, atol=1e-12)
    assert np.isfinite(increment).all()
    if periodic:
        assert_allclose(increment[:, 0], increment[:, -1])


def test_particle_increment_gathers_node_indices():
    P, G, L = state(ns=3, segregate=True)
    b = L.materials[0]
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    G.m[:] = 1
    G.get_valid_nodes(P)
    G.dphi[:] = np.arange(len(G.m))[:, None]*np.array([1, -2, 1])*1e-5
    G.dpk[:] = np.arange(len(G.m))*1e-5
    expected_phi = b.phi.copy()
    expected_pk = b.pk.copy()
    for i in range(len(b)):
        for r, n in enumerate(b.n[i]):
            expected_phi[i] += b.N[i, r]*G.dphi[n]
            expected_pk[i] += b.N[i, r]*G.dpk[n]
    b.move_grainsize_on_grid(P, G)
    assert_allclose(b.phi, expected_phi)
    assert_allclose(b.pk, expected_pk)


def test_normalisation_is_opt_in():
    P, G, L = state(ns=2, segregate=True)
    b = L.materials[0]
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    G.m[:] = 1
    G.get_valid_nodes(P)
    G.dphi[:] = [-1, 1]
    b.move_grainsize_on_grid(P, G)
    assert np.any(b.phi < 0)
    G.dphi.fill(0)
    P.normalise_phi = True
    b.move_grainsize_on_grid(P, G)
    assert_allclose(b.phi, np.broadcast_to([0, 1], b.phi.shape))


def test_periodic_two_node_direction():
    P, G, L = state(ns=3, segregate=True, B=NS(cyclic_lr=True))
    P.G.nx = 2
    from grid import Grid
    G = Grid(P)
    G.m[:] = 1
    G.phi[:] = [0.2, .3, .5]
    G.I[:] = 1
    G.s_bar = G.phi@P.G.s
    G.shear_rate = np.zeros_like(G.m)
    assert_allclose(increment_grainsize(P, G), 0)


@pytest.mark.parametrize('ns', [2, 3, 8])
def test_hundred_step_segregating_example(ns):
    from initialise import get_parameters
    P, G, L = get_parameters(['bi_seg_test', '22', str(ns), '11'])
    assert not P.normalise_phi
    run(P, G, L, max_steps=100, output=False)
    b = L.materials[0]
    assert np.isfinite(b.phi).all()
    assert np.all((b.phi >= 0) & (b.phi <= 1))
    assert P.tstep == 100
    # This scheme has documented simplex/remapping drift. Do not hide it by
    # normalisation; exact conservation is tested independently at each transfer.
