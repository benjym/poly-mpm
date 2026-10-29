"""Acceptance checks for the critical-state CIDR granular step (cidr.py)."""

import numpy as np
import pytest

from cidr import ROOT2, CriticalStateCIDR, Prox
from initialise import prepare
from MPM import run

# Glass beads (Jop et al.), 3 mm grains; phi = 0.6 sits on the cap at twice the base pressure
# of a 1 m column, so the column at rest is denser than the cap for its load.
CIDR = dict(
    packing_max=0.605,
    mu=np.tan(np.radians(20.9)),
    mu_d=np.tan(np.radians(32.76)),
    I0=0.279,
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
    model = CriticalStateCIDR(CIDR["mu"], CIDR["mu_d"], CIDR["I0"], 2 * 0.003 * np.sqrt(2500))
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
