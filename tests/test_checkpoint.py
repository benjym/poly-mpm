"""A restarted run continues exactly as an uninterrupted one."""

import contextlib
import io

import numpy as np
import pytest

import checkpoint
from initialise import get_parameters, prepare
from MPM import run


def _state(args, **options):
    with contextlib.redirect_stdout(io.StringIO()):
        if args[0] == "granular_column":
            from inputs.granular_column import Params

            P = Params([])
            for name, value in options.pop("B", {}).items():
                setattr(P.B, name, value)
            for name, value in options.items():
                setattr(P, name, value)
            return prepare(P)
        P, G, L = get_parameters(args)
    return P, G, L


def _same(a, b):
    if isinstance(a, np.ndarray):
        return isinstance(b, np.ndarray) and a.shape == b.shape and np.array_equal(a, b, equal_nan=True)
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return a == b or (a != a and b != b)


CASES = [
    (["ball"], {}),
    (["modern_ball"], {}),
    (["bi_seg_test", "22", "2", "11"], {}),
    (["domain_column"], {}),
    (["granular_column"], dict(dt=0.005, t_f=0.2, granular_solver="velocity", granular_tolerance=1e-9)),
    (
        ["granular_column"],
        dict(dt=0.005, t_f=0.2, max_g=2.0, granular_solver="velocity", B=dict(separable=True)),
    ),
    (
        ["granular_column"],
        dict(dt=0.005, t_f=0.2, theta=0.3, granular_solver="velocity", B=dict(wall_mu=0.2)),
    ),
]


@pytest.mark.parametrize("args, options", CASES, ids=lambda case: str(case)[:40])
def test_restart_matches_an_uninterrupted_run(tmp_path, args, options):
    total, first = 12, 5
    reference = _state(args, **dict(options, B=dict(options.get("B", {}))))
    run(*reference, max_steps=total, output=False)

    interrupted = _state(args, **dict(options, B=dict(options.get("B", {}))))
    run(*interrupted, max_steps=first, output=False)
    path = tmp_path / "checkpoint.pkl"
    checkpoint.save(path, *interrupted)

    restarted = _state(args, **dict(options, B=dict(options.get("B", {}))))
    checkpoint.load(path, *restarted)
    run(*restarted, max_steps=total - first, output=False)

    (P0, G0, L0), (P1, G1, L1) = reference, restarted
    assert P1.tstep == P0.tstep == total
    assert _same(checkpoint._state(P0), checkpoint._state(P1))
    assert _same(checkpoint._state(G0), checkpoint._state(G1))
    for a, b in zip(L0.materials, L1.materials):
        assert _same(checkpoint._state(a), checkpoint._state(b))
    assert _same(checkpoint._state(L0.diagnostics), checkpoint._state(L1.diagnostics))


def test_checkpoint_output_flag_writes_a_loadable_file(tmp_path):
    P, G, L = _state(["ball"])
    P.save_dir = str(tmp_path) + "/"
    P.O.save_checkpoint = True
    P.savetime = 2 * P.dt
    for name in ("plot_continuum", "plot_material_points"):
        setattr(P.O, name, False)
    run(P, G, L, max_steps=4)
    path = tmp_path / "data" / "checkpoint.pkl"
    assert path.exists()
    restarted = _state(["ball"])
    checkpoint.load(path, *restarted)
    assert restarted[0].tstep == 4
    assert np.array_equal(restarted[2].materials[0].x, L.materials[0].x)


def test_checkpoint_rejects_a_different_material_count(tmp_path):
    P, G, L = _state(["ball"])
    path = tmp_path / "checkpoint.pkl"
    L.materials.append(L.materials[0])
    checkpoint.save(path, P, G, L)
    with pytest.raises(ValueError, match="material counts"):
        checkpoint.load(path, *_state(["ball"]))
