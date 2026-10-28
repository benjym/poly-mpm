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
@pytest.mark.parametrize('steps', [100, 1000])
def test_segregating_example_conserves_particle_state(ns, steps):
    from initialise import get_parameters
    P, G, L = get_parameters(['bi_seg_test', '22', str(ns), '11'])
    assert not P.normalise_phi
    initial_mass = (L.materials[0].m[:, None]*L.materials[0].phi).sum(axis=0)
    run(P, G, L, max_steps=steps, output=False)
    b = L.materials[0]
    assert np.isfinite(b.phi).all()
    assert np.all((b.phi >= -1e-14) & (b.phi <= 1+1e-14))
    assert P.tstep == steps
    assert_allclose(b.phi.sum(axis=1), 1, atol=2e-14)
    assert_allclose((b.m[:, None]*b.phi).sum(axis=0), initial_mass, rtol=1e-13, atol=1e-12)


@pytest.mark.parametrize('ns', [2, 3, 8])
@pytest.mark.parametrize('ax', [0, 1])
@pytest.mark.parametrize('periodic', [False, True])
def test_nonuniform_mass_and_local_class_closure(ns, ax, periodic):
    P, G, L = state(ns=ns, segregate=True, c=.001, B=NS(cyclic_lr=periodic))
    rng = np.random.default_rng(12)
    G.m[:] = rng.uniform(.2, 5, len(G.m))
    G.get_valid_nodes(P)
    G.phi[:] = rng.uniform(.1, 1, size=G.phi.shape)
    G.phi /= G.phi.sum(axis=1, keepdims=True)
    G.I[:] = G.m*.1
    G.shear_rate = np.full_like(G.m, 2.)
    G.grad_pk[:, ax] = rng.uniform(-.5, .5, len(G.m))
    if periodic:
        for value in (G.m, G.phi, G.I, G.grad_pk):
            value[P.G.nx-1::P.G.nx] = value[::P.G.nx]
    change = KT(P, G, ax)
    assert_allclose(change.sum(axis=1), 0, atol=1e-14)
    mass_change = (G.m[:, None]*change).reshape(P.G.ny, P.G.nx, ns)
    if periodic:
        mass_change = mass_change[:, :-1]
    assert_allclose(mass_change.sum(axis=(0, 1)), 0, atol=1e-14)
    assert np.max(np.abs(change)) > 0  # conservation is not a disabled update


@pytest.mark.parametrize('ns', [2, 3, 8])
def test_physical_diffusion_flux_has_no_extra_concentration_factor(ns):
    from integrators import KT2_flux, pad_KT
    P, G, L = state(ns=ns, segregate=True, c=0.)
    G.m[:] = 1
    G.get_valid_nodes(P)
    rng = np.random.default_rng(33)
    G.phi[:] = rng.uniform(.1, 1, size=G.phi.shape)
    G.phi /= G.phi.sum(axis=1, keepdims=True)
    G.I[:] = 1
    G.shear_rate = np.ones_like(G.m)
    padded, dx = pad_KT(G.phi, P, 1)
    size = G.phi@P.G.s
    D = (10*size**3).reshape(P.G.ny, P.G.nx)
    D = np.pad(D, ((2, 2), (0, 0)), mode='edge')
    flux = KT2_flux(padded, G, P, 1)
    assert_allclose(flux, -D[..., None]*np.gradient(padded, dx, axis=0), rtol=1e-12, atol=1e-20)
    assert_allclose(flux.sum(axis=-1), 0, atol=1e-20)


@pytest.mark.parametrize('periodic', [False, True])
def test_coupled_particle_species_mass_with_unequal_masses_and_rigid_wall(periodic):
    rng = np.random.default_rng(18)
    count = 200
    phi = rng.uniform(.2, 1, (count, 3))
    phi /= phi.sum(axis=1, keepdims=True)
    mobile = material(n=count, X=rng.uniform(.05, .95, count), Y=rng.uniform(.05, .95, count),
                      rho=rng.uniform(1000, 2000, count), A=rng.uniform(.001, .003, count), PHI=phi)
    rigid = material('rigid', n=10, X=np.linspace(.1, .9, 10), Y=np.full(10, .5), phi=[.9, .05, .05])
    P, G, L = state(ns=3, materials=[mobile, rigid], segregate=True, c=.001, B=NS(cyclic_lr=periodic))
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    b, wall = L.materials
    initial = (b.m[:, None]*b.phi).sum(axis=0)
    wall_before = wall.phi.copy()
    mobile_mass = G.segregation_m.reshape(P.G.ny, P.G.nx)
    if periodic:
        mobile_mass = mobile_mass[:, :-1]
    assert_allclose(mobile_mass.sum(), b.m.sum())
    G.I[:] = G.m*.1
    G.shear_rate = np.ones_like(G.m)
    G.grad_pk[:] = rng.uniform(-.1, .1, size=G.grad_pk.shape)
    if periodic:
        G.grad_pk[P.G.nx-1::P.G.nx] = G.grad_pk[::P.G.nx]
    G.calculate_phi_increment(P, L)
    L.move_grainsize_on_grid(P, G)
    assert_allclose((b.m[:, None]*b.phi).sum(axis=0), initial, rtol=1e-13, atol=1e-12)
    assert_allclose(b.phi.sum(axis=1), 1, atol=1e-14)
    assert_allclose(wall.phi, wall_before)
    assert np.max(np.abs(G.dphi)) > 0


def test_unequal_mass_update_conserves_across_a_free_surface():
    P, G, L = state(ns=2, segregate=True, c=.001)
    G.m[:] = np.linspace(.5, 3., len(G.m))
    empty = np.arange(P.G.nx*3, len(G.m))
    G.m[empty] = 0
    G.get_valid_nodes(P)
    G.phi[:] = [0.3, .7]
    G.phi[empty] = 0
    G.I[:] = G.m*.1
    G.shear_rate = np.ones_like(G.m)
    G.grad_pk[:, 1] = .5
    change = increment_grainsize(P, G)
    assert_allclose(change[empty], 0)
    assert_allclose(change.sum(axis=1), 0, atol=1e-14)
    assert_allclose(np.sum(G.m[:, None]*change, axis=0), 0, atol=1e-14)


@pytest.mark.parametrize('ns', [2, 3, 8])
@pytest.mark.parametrize('periodic', [False, True])
def test_particle_capacity_limits_joint_xy_flux_without_clipping(ns, periodic):
    rng = np.random.default_rng(106)
    count = 200
    phi = rng.uniform(.01, 1, (count, ns))
    phi /= phi.sum(axis=1, keepdims=True)
    mobile = material(n=count, X=rng.uniform(.01, .99, count), Y=rng.uniform(.01, .99, count),
                      rho=rng.uniform(1000, 2000, count), A=rng.uniform(.001, .003, count), PHI=phi)
    P, G, L = state(ns=ns, materials=mobile, segregate=True, c=1., B=NS(cyclic_lr=periodic))
    P.dt = 100.  # Deliberately exceed the unbounded explicit transport step.
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    G.I[:] = G.m*.1
    G.shear_rate = np.ones_like(G.m)
    G.grad_pk[:] = rng.uniform(-10, 10, size=G.grad_pk.shape)
    if periodic:
        G.grad_pk[P.G.nx-1::P.G.nx] = G.grad_pk[::P.G.nx]
    b = L.materials[0]
    initial = (b.m[:, None]*b.phi).sum(axis=0)
    unlimited = increment_grainsize(P, G)
    G.calculate_phi_increment(P, L)
    assert np.max(np.abs(G.dphi)) < np.max(np.abs(unlimited))
    assert np.max(np.abs(G.dphi)) > 0
    L.move_grainsize_on_grid(P, G)
    assert not P.normalise_phi
    assert np.all((b.phi >= -1e-14) & (b.phi <= 1+1e-14))
    assert_allclose(b.phi.sum(axis=1), 1, atol=1e-13)
    assert_allclose((b.m[:, None]*b.phi).sum(axis=0), initial, rtol=1e-13, atol=1e-12)


def test_pure_layers_mix_conservatively_without_negative_fractions():
    count = 40
    y = np.linspace(.01, .99, count)
    phi = np.column_stack((y < .5, y >= .5)).astype(float)
    mobile = material(n=count, X=np.full(count, .5), Y=y, PHI=phi)
    P, G, L = state(ns=2, materials=mobile, segregate=True, c=0.)
    P.dt = 1000.
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    G.I[:] = G.m*.1
    G.shear_rate = np.full_like(G.m, 100.)
    b = L.materials[0]
    initial = (b.m[:, None]*b.phi).sum(axis=0)
    G.calculate_phi_increment(P, L)
    L.move_grainsize_on_grid(P, G)
    assert np.any((b.phi > 1e-8) & (b.phi < 1-1e-8))
    assert np.all((b.phi >= -1e-14) & (b.phi <= 1+1e-14))
    assert_allclose(b.phi.sum(axis=1), 1, atol=1e-14)
    assert_allclose((b.m[:, None]*b.phi).sum(axis=0), initial, atol=1e-12)


@pytest.mark.parametrize('ns', [3, 8])
def test_flux_does_not_depend_on_grain_class_order(ns):
    P, G, L = state(ns=ns, segregate=True, c=.01, B=NS(cyclic_lr=True))
    rng = np.random.default_rng(71)
    G.m[:] = rng.uniform(.2, 2, len(G.m))
    G.phi[:] = rng.uniform(.1, 1, G.phi.shape)
    G.phi /= G.phi.sum(axis=1, keepdims=True)
    G.I[:] = G.m*.1
    G.shear_rate = np.ones_like(G.m)
    G.grad_pk[:] = rng.uniform(-.1, .1, G.grad_pk.shape)
    for value in (G.m, G.phi, G.I, G.grad_pk):
        value[P.G.nx-1::P.G.nx] = value[::P.G.nx]
    before = increment_grainsize(P, G)
    order = rng.permutation(ns)
    G.phi = G.phi[:, order]
    P.G.s = P.G.s[order]
    after = increment_grainsize(P, G)
    assert_allclose(after, before[:, order], rtol=1e-12, atol=1e-14)


def test_diffusion_exchange_matches_discrete_heat_equation(monkeypatch):
    import integrators
    x = np.arange(8)/8
    phi0 = .5+.2*np.sin(2*np.pi*x)
    material0 = material(n=8, X=x, Y=np.full(8, .5), PHI=np.column_stack((phi0, 1-phi0)))
    grid = NS(nx=9, ny=3, x_m=0., x_M=1., y_m=0., y_M=1., ns=2, s=[.001, .003])
    P, G, L = state(ns=2, materials=material0, G=grid, segregate=True, c=0., B=NS(cyclic_lr=True))
    diffusivity = .01
    monkeypatch.setattr(integrators, '_diffusivity', lambda G, P, phi: np.full(len(G.m), diffusivity))
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    expected = P.dt*diffusivity/G.dx**2*(np.roll(phi0, 1)-2*phi0+np.roll(phi0, -1))
    G.calculate_phi_increment(P, L)
    L.move_grainsize_on_grid(P, G)
    assert_allclose(L.materials[0].phi[:, 0], phi0+expected, atol=1e-14)
    assert_allclose(L.materials[0].phi[:, 1], 1-phi0-expected, atol=1e-14)
