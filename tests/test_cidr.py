"""Acceptance checks for the critical-state CIDR granular step (cidr.py)."""

import numpy as np
import pytest

from cidr import ROOT2, CriticalStateCIDR, Prox
from initialise import prepare
from MPM import run

# Glass beads (Jop et al.), 3 mm grains; phi = 0.6 sits on the cap at twice the base pressure
# of a 1 m column, so the column at rest is denser than the cap for its load.
CIDR = dict(
    phi_max=0.605,
    mu_s=np.tan(np.radians(20.9)),
    mu_d=np.tan(np.radians(32.76)),
    I_0=0.279,
    grain_diameter=0.003,
    phi_min=0.05,
    cidr_deficit=0.005,
    cidr_p_ref=2 * 1500 * 9.81 * 1.0,
)


def cidr_case(module, **options):
    Params = __import__(f"inputs.{module}", fromlist=["Params"]).Params
    P = Params([])
    P.granular_solver = "cidr"
    for key, value in options.items():
        setattr(P, key, value)
    for key, value in CIDR.items():
        setattr(P.S, key, value)
    return prepare(P)


def test_proximal_step_is_stationary_or_jammed():
    model = CriticalStateCIDR(CIDR["mu_s"], CIDR["mu_d"], CIDR["I_0"], 2 * 0.003 * np.sqrt(2500))
    rng = np.random.default_rng(0)
    q = rng.standard_normal((120, 3)) * 10 ** rng.uniform(-4, 0, (120, 1))
    C = np.full(len(q), 4e4) * 10 ** rng.uniform(-0.5, 0.5, len(q))
    r = 1500 / 0.002
    d, jam, _ = Prox(model)(q, r, C)
    assert 0 < jam.sum() < len(q)

    def psi(D, C):
        s = np.maximum(np.linalg.norm(D[:, 1:], axis=1), 1e-300) / ROOT2
        return model.evaluate_C(s, ROOT2 * D[:, 0], C)[0][0]

    free = ~jam
    D, Cf = d[free], C[free]
    h = 1e-7 * np.maximum(np.linalg.norm(D, axis=1), 1e-12)[:, None]
    gradient = np.stack(
        [(psi(D + h * e, Cf) - psi(D - h * e, Cf)) / (2 * h[:, 0]) for e in np.eye(3)], axis=1
    )
    residual = np.linalg.norm(gradient + r * (D - q[free]), axis=1) / np.linalg.norm(r * q[free], axis=1)
    assert residual.max() < 1e-6
    # Jammed nodes: r q lies in the static set, so its support inequality holds in every direction.
    E = rng.standard_normal((4000, 3))
    E /= np.linalg.norm(E, axis=1)[:, None]
    support = model.static_support(
        np.linalg.norm(E[:, 1:], axis=1)[None, :] / ROOT2, ROOT2 * E[None, :, 0], C[jam, None]
    )
    excess = (r * q[jam]) @ E.T - support
    assert np.max(excess / np.linalg.norm(r * q[jam], axis=1)[:, None]) < 1e-6


def test_static_column_jams_without_creep():
    P, G, L = cidr_case("granular_column", dt=0.02, t_f=0.2)
    initial = L.materials[0].x.copy()
    run(P, G, L, output=False)
    b = L.materials[0]
    assert np.max(np.abs(b.v)) < 1e-6
    assert np.max(np.abs(b.x - initial)) < 1e-8
    assert b.pressure.max() < 0
    history = [r["granular_solver"] for r in L.diagnostics.history]
    assert not history[0]["warm"] and all(h["warm"] for h in history[1:])
    assert max(h["newton"] for h in history[1:]) <= 10
    assert np.linalg.norm(L.diagnostics.last["momentum_residual"]) < 1e-9


def test_collapse_flows_and_dissipates():
    P, G, L = cidr_case("granular_collapse", dt=0.01, t_f=0.1)
    initial = L.materials[0].x.copy()
    run(P, G, L, output=False)
    b = L.materials[0]
    assert b.x[:, 0].max() > initial[:, 0].max() + 0.01
    assert np.linalg.norm(b.v, axis=1).max() > 0.1
    history = [r["granular_solver"] for r in L.diagnostics.history]
    assert max(h["contact_power"] for h in history) <= 1e-7
    assert max(h["newton"] for h in history[1:]) <= 40


@pytest.mark.parametrize("missing", ["mu_d", "grain_diameter", "cidr_p_ref"])
def test_cidr_requires_its_material_fields(missing):
    from inputs.granular_column import Params

    P = Params([])
    P.granular_solver = "cidr"
    for key, value in CIDR.items():
        if key != missing:
            setattr(P.S, key, value)
    with pytest.raises(ValueError, match=missing):
        prepare(P)


def test_a_step_that_would_compact_past_phi_max_is_stopped():
    P, G, L = cidr_case("granular_collapse", dt=0.05, t_f=0.3)
    with pytest.raises(FloatingPointError, match="past phi_max"):
        run(P, G, L, output=False)


def test_dynamic_stepping_limits_the_packing_closure_per_step():
    P, G, L = cidr_case("granular_collapse", dt=0.02, t_f=0.3, time_stepping="dynamic", cfl=0.2)
    run(P, G, L, output=False)
    assert "packing" in P.timestep_limits
    history = [r["granular_solver"]["packing_closure"] for r in L.diagnostics.history]
    # The limit uses the previous step's closing rate, so the first step is unconstrained.
    assert max(history[1:]) < 0.5


def test_envelope_hessian_keeps_its_precision_at_small_gamma():
    # At small gamma the envelope's Hessian tends to psi's own; written as
    # r (I - r (H + r I)^-1) it lost about log10(r / |H|) digits and stalled Newton.
    from cidr import Problem

    model = CriticalStateCIDR(CIDR["mu_s"], CIDR["mu_d"], CIDR["I_0"], 2 * 0.003 * np.sqrt(2500))
    C = np.array([50.0, 200.0])
    problem = Problem(np.zeros((6, 1)), np.zeros(6), np.ones(6), np.ones(1), np.zeros(1), model, C)
    x = np.array([[-0.006, 0.006, -3.4], [-0.008, 0.008, -5.6]])
    gamma = 5e-9
    *_, H, jam = problem.envelope(x, gamma)
    assert not jam.any()
    r = 1 / gamma
    d, _, _ = problem.prox(x, r)
    Hpsi = problem.hessian(d, C)
    # r H (H + r I)^-1 = H - H^2 / r + H^3 / r^2 - ..., |H| / r about 1e-6 here.
    series = Hpsi - Hpsi @ Hpsi / r + Hpsi @ Hpsi @ Hpsi / r**2
    assert np.abs(H - series).max() <= 1e-12 * np.abs(Hpsi).max()


def test_closing_rate_leaves_out_nodes_already_at_phi_max():
    from granular import closing_rate

    delta = np.array([-1.0, -1e-6, -1e-6])
    phi = np.array([0.5, 0.605, 0.61])
    assert closing_rate(delta, phi, 0.605) == pytest.approx(1.0 / (0.105 / 0.5))


def test_incline_holds_the_bagnold_profile():
    # The steady state of benchmarks/cidr_incline.py, over a short run on the CPU.
    from benchmarks.cidr_incline import profile
    from initialise import get_parameters

    P, G, L = get_parameters(["cidr_incline"])
    P.t_f = 0.01
    run(P, G, L, output=False)
    result = profile(P, G, L, 18)
    ratio = np.array(result["shear_rate"]) / np.array(result["bagnold_shear_rate"])
    assert np.all(np.abs(ratio[2:15] - 1) < 0.03)


def sliding_layer(degrees, backend="cpu", t_f=0.02):
    """The incline's layer at rest and uniform packing 0.6 on a Coulomb floor (friction 0.2)."""
    import inputs.cidr_incline as incline

    angle = incline.ANGLE
    incline.ANGLE = np.radians(degrees)
    try:
        P = incline.Params(["cidr_incline"])
    finally:
        incline.ANGLE = angle
    P.backend, P.t_f = backend, t_f
    P.B.wall_contact, P.B.wall_friction = "coulomb", 0.2
    P.S.initial_v[:] = 0.0
    P.S.rho = np.full(P.S.n, 0.6 * P.S.rho_s)
    P, G, L = prepare(P)
    run(P, G, L, output=False)
    return P, L.materials[0]


@pytest.mark.parametrize("degrees", [15.0, 10.0])
def test_coulomb_floor_holds_or_slides_a_jammed_layer(degrees):
    # The layer is jammed (tan theta < mu_s); on the floor it slides rigidly at
    # g (sin theta - 0.2 cos theta) when tan theta > 0.2 and otherwise stays at rest.
    P, b = sliding_layer(degrees)
    theta = np.radians(degrees)
    expected = max(9.81 * (np.sin(theta) - 0.2 * np.cos(theta)), 0.0)
    assert np.allclose(-b.v[:, 0] / P.t, expected, atol=1e-4)
    assert np.max(np.abs(b.v[:, 1])) < 1e-6
