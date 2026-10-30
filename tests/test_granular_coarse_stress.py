"""Coarse Q1 stress basis (stresses on a grid k times the velocity spacing)."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from granular import BASIS, assemble, coarse_stencil, coarse_volumes
from initialise import prepare
from MPM import run


def column(**options):
    from inputs.granular_column import Params

    P = Params([])
    P.granular_solver, P.granular_stress_basis = "velocity", "coarse"
    for key, value in options.items():
        setattr(P, key, value)
    return prepare(P)


def start(P, G, L):
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    return L.materials[0]


@pytest.mark.parametrize("periodic", [False, True])
def test_coarse_basis_partition_volumes_and_adjoint(periodic):
    P, G, L = column()
    if periodic:
        P.B.cyclic_lr = True
        P.B.has_left = P.B.has_right = False
        G.boundary_v[:] = False
    b = start(P, G, L)
    nodes, weights = coarse_stencil(P, G, b.x)
    assert_allclose(weights.sum(axis=1), 1, atol=1e-15)
    area = (P.G.x_max - P.G.x_min) * (P.G.y_max - P.G.y_min) * P.G.thickness
    assert_allclose(coarse_volumes(P, G).sum(), area, rtol=1e-14)
    compact, _, active, occupied, gap = assemble(P, G, b, compact=True)
    assembled, *_ = assemble(P, G, b, compact=False)
    assert_allclose(occupied.sum(), np.sum(b.m / b.material.rho_s), rtol=1e-14)
    rng = np.random.default_rng(3)
    u, z = rng.normal(size=2 * len(G.m)), rng.normal(size=3 * len(active))
    assert_allclose(compact @ u, assembled @ u, atol=1e-12 * np.abs(assembled @ u).max())
    assert abs((compact @ u) @ z - u @ compact.transpose(z)) <= 1e-12 * (np.abs(compact @ u) @ np.abs(z))
    assert compact.count < len(G.m)


def test_coarse_basis_integrates_affine_velocity_gradient():
    P, G, L = column()
    b = start(P, G, L)
    B, _, active, occupied, _ = assemble(P, G, b, compact=True)
    C = np.array([[0.2, -0.1], [0.3, -0.2]])
    u = np.column_stack((G.X, G.Y)) @ C.T
    expected = occupied[active, None] * np.einsum("kij,ij->k", BASIS[:, :2, :2], C)
    assert_allclose((B @ u.ravel()).reshape(-1, 3), expected, atol=1e-14)


def test_coarse_basis_static_support():
    P, G, L = column(dt=0.02, t_f=0.2, granular_tolerance=1e-9)
    before = L.materials[0].x.copy()
    run(P, G, L, output=False)
    b = L.materials[0]
    assert_allclose(b.x, before, atol=1e-8)
    assert np.max(np.abs(b.v)) < 1e-7
    assert b.pressure.min() < -100 and b.pressure.max() <= 1e-6
    assert G.granular_stats["stress_basis"] == "coarse"


def test_coarse_basis_requires_divisible_cell_counts():
    from inputs.granular_column import Params

    P = Params([])
    P.granular_stress_basis = "coarse"
    P.G.nx = 6  # five cells
    with pytest.raises(ValueError, match="divisible by 2"):
        prepare(P)
    P = Params([])
    P.granular_stress_basis, P.granular_stress_coarsening = "coarse", 1
    with pytest.raises(ValueError, match="coarsening"):
        prepare(P)
