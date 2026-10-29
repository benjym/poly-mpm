"""Opt-in increment stopping for the compact velocity solver (Daviet 2026 style)."""

import numpy as np
import pytest

from initialise import prepare
from MPM import run


def collapse(**options):
    from inputs.granular_collapse import Params

    P = Params([])
    P.granular_solver, P.diagnostics = "velocity", True
    for key, value in options.items():
        setattr(P, key, value)
    return prepare(P)


def test_increment_stopping_reports_true_residual_without_raising():
    # A sweep budget far too small for the certified tolerance.
    P, G, L = collapse(
        granular_stopping="increment", granular_tolerance=1e-9, granular_max_iterations=2, t_f=0.02
    )
    run(P, G, L, output=False)
    stats = [r["granular_solver"] for r in L.diagnostics.history]
    assert all(s["stopping"] == "increment" for s in stats)
    assert all(s["sweeps"] <= 20 for s in stats)
    assert any(not s["converged"] and s["residual"] > 1e-9 for s in stats)
    assert all(np.isfinite(s["residual"]) for s in stats)


def test_increment_stopping_ends_on_small_change():
    P, G, L = collapse(granular_stopping="increment", granular_increment=1e-3, t_f=0.02)
    run(P, G, L, output=False)
    for s in (r["granular_solver"] for r in L.diagnostics.history):
        if s["sweeps"] < 10 * P.granular_max_iterations and not s["converged"]:
            assert s["increment_max"] <= 1e-3 and s["increment_rms"] <= 1e-3


def test_residual_stopping_still_raises():
    P, G, L = collapse(granular_tolerance=1e-9, granular_max_iterations=1, t_f=0.02)
    with pytest.raises(RuntimeError, match="did not converge"):
        run(P, G, L, output=False)


@pytest.mark.parametrize(
    "options, message",
    [
        (dict(granular_stopping="sometimes"), "granular_stopping"),
        (dict(granular_stopping="increment", granular_solver="hybrid"), "requires granular_solver"),
        (dict(granular_stopping="increment", granular_increment=0.0), "granular_increment"),
    ],
)
def test_increment_stopping_validation(options, message):
    with pytest.raises(ValueError, match=message):
        collapse(**options)
