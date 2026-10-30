"""Frozen original-Sand6 assembly oracle and source-style option checks."""

import json
from pathlib import Path
import numpy as np
import pytest
from numpy.testing import assert_allclose

from benchmarks.granular_prototype import configuration
from diagnostics import inventory
from granular import assemble, viscosity_matrix, coulomb_mfb, solve_coulomb
from initialise import prepare
from MPM import run, step


def frozen():
    data = json.loads((Path(__file__).parent / "fixtures/sand6_frozen.json").read_text())
    P = configuration("static", 0.005)
    P.G.nx = P.G.ny = 2
    P.G.x_max = P.G.y_max = 1.0
    P.granular_stress = "plane_strain"
    P.granular_viscosity = 0.001
    points = np.array(data["points"])
    P.S.X, P.S.Y = points.T
    P.S.n = 4
    P.S.A = 0.25
    P.S.rho = 750.0
    P, G, L = prepare(P)
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    return P, G, L, data


def test_original_cpp_assembly_matches_plane_strain_reduction():
    P, G, L, data = frozen()
    B, _, _, _, _ = assemble(P, G, L.materials[0])
    # Original chi is orthonormal for (tau:sigma)/2, ours for tau:sigma.
    assert_allclose(np.sqrt(2) * B.toarray(), data["B_source"], atol=1e-16)
    K = viscosity_matrix(P, G, L.materials[0])
    assert_allclose(K.toarray(), 2 * P.granular_viscosity * np.array(data["D_inner_product"]), atol=1e-18)
    assert_allclose(G.m, 187.5, atol=1e-12)


def test_source_gradient_transfer_uses_carried_gradient_and_affine_inventory():
    P = configuration("static", 0.005)
    P.transfer = "apic_gradient"
    C = np.array([[0.2, 0.3, 0.0], [-0.1, 0.4, 0.0], [0.0, 0.0, 0.0]])
    P.S.initial_C = C
    P, G, L = prepare(P)
    b = L.materials[0]
    G.wipe(P)
    b.get_nodal_mass_momentum(P, G)
    reference = np.zeros_like(G.q)
    for p in range(len(b)):
        for j, node in enumerate(b.n[p]):
            reference[node] += b.m[p] * b.N[p, j] * (b.v[p] + C @ b.offsets[p, j])
    assert_allclose(G.q, reference, atol=1e-13)
    assert inventory(P, L)["affine_kinetic"] > 0
    previous = b.L.copy()
    b.x[:, 0] += 0.01
    b.get_reference_nodes(P, G)
    b.get_basis_functions(P, G)
    assert_allclose(b.L, previous, atol=0)


@pytest.mark.parametrize("periodic", [False, True])
@pytest.mark.parametrize("method", ["hybrid", "velocity"])
def test_source_viscosity_stage_is_dissipative_and_work_is_accounted(periodic, method):
    P = configuration("static", 0.005)
    P.granular_solver = method
    P.transfer = "apic_gradient"
    P.granular_viscosity = 10.0
    P.gravity = 0.0
    P.S.rho = 150.0
    if periodic:
        P.B.walls = set(P.B.walls) - {"left", "right"}
        P.B.periodic = True
    P.S.initial_v = np.column_stack((0.1 * np.sin(3 * P.S.Y), np.zeros(P.S.n), np.zeros(P.S.n)))
    P, G, L = prepare(P)
    step(P, G, L)
    work = L.materials[0].viscous_work.sum()
    assert work > 0
    assert_allclose(work, G.granular_stats["viscosity_trial_dissipation"], rtol=1e-12)
    assert np.linalg.norm(L.diagnostics.last["momentum_residual"]) < 1e-9


def test_plane_strain_static_support_and_three_dimensional_yield_surface():
    P = configuration("static", 0.02)
    P.granular_stress = "plane_strain"
    P.transfer = "apic_gradient"
    P.t_f = 0.04
    P, G, L = prepare(P)
    run(P, G, L, output=False)
    b = L.materials[0]
    assert np.max(np.abs(b.v)) < 1e-8
    assert b.yieldfunction.max() < 1e-6
    assert G.granular_contact.shape[1] == 4


@pytest.mark.parametrize(
    "setting",
    [
        dict(granular_stress="invalid"),
        dict(granular_viscosity=-1),
        dict(granular_stress="plane_strain", granular_solver="block"),
    ],
)
def test_invalid_source_options_fail_before_stepping(setting):
    P = configuration("static", 0.005)
    for name, value in setting.items():
        setattr(P, name, value)
    with pytest.raises(ValueError):
        prepare(P)


@pytest.mark.parametrize(
    "reaction,velocity",
    [
        ([1.0, -0.2, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]),
        ([1.0, -0.5, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]),
        ([0.0, 0.0, 0.0, 0.0], [1.0, -3.0, 0.0, 0.0]),
    ],
)
def test_modified_fb_has_the_same_stick_slide_open_roots(reaction, velocity):
    assert_allclose(coulomb_mfb(np.array([reaction]), np.array([velocity]), 0.5), 0, atol=1e-14)
    if reaction[0] and velocity[1]:
        velocity[0] = 0.5  # Associated dilation must not solve nonassociated Coulomb.
        assert np.linalg.norm(coulomb_mfb(np.array([reaction]), np.array([velocity]), 0.5)) > 0.01


def test_plane_strain_merit_fallback_is_checked_against_original_residual(monkeypatch):
    from types import SimpleNamespace
    import granular

    original = granular.least_squares

    def stall_natural_map(fun, start, **kw):
        if fun.__name__ == "evaluate":
            return SimpleNamespace(x=start.copy(), nfev=1)
        return original(fun, start, **kw)

    monkeypatch.setattr(granular, "least_squares", stall_natural_map)
    r, stats = solve_coulomb(np.eye(4), [[-1.0, 1.0, 0.0, 0.0]], 0.5)
    assert stats["mfb_evaluations"] > 0
    assert stats["residual"] < 1e-9
    assert_allclose(r, [[1.0, -0.5, 0.0, 0.0]], atol=1e-9)


def test_source_style_plane_strain_collapse_completes_without_accepting_stalled_root():
    P = configuration("collapse", 0.005)
    P.transfer = "apic_gradient"
    P.granular_stress = "plane_strain"
    P.granular_viscosity = 0.001
    P, G, L = prepare(P)
    run(P, G, L, output=False)
    assert P.t == P.t_f
    assert max(r["granular_solver"]["residual"] for r in L.diagnostics.history) <= P.granular_tolerance
    assert np.linalg.norm(L.diagnostics.last["momentum_residual"]) < 1e-8
