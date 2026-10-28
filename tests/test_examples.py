import contextlib
import io
import json
from pathlib import Path
import numpy as np
import pytest
from initialise import get_parameters
from MPM import main, run
from benchmarks.audit_examples import ARGUMENTS

ROWS = json.loads((Path(__file__).resolve().parents[1] / "benchmarks/examples.json").read_text())


@pytest.mark.parametrize(
    "args, message",
    [
        ([], "usage:"),
        (["11"], "usage:"),
        (["11", "2", "extra"], "usage:"),
        (["abc", "2"], "grid_nodes: must be an integer"),
        (["2", "2"], "grid_nodes: must be at least 3"),
        (["11", "abc"], "Froude_number: must be a number"),
        (["11", "0"], "Froude_number: must be positive"),
        (["11", "-1"], "Froude_number: must be positive"),
        (["11", "nan"], "Froude_number: must be finite"),
        (["11", "inf"], "Froude_number: must be finite"),
    ],
)
def test_square_drum_invalid_arguments(args, message, capsys):
    assert main(["bi_square_drum", *args]) == 2
    captured = capsys.readouterr()
    assert message in captured.err
    assert "bi_square_drum" in captured.err
    assert "grid_nodes" in captured.err and "Froude_number" in captured.err
    assert "Traceback" not in captured.err


@pytest.mark.parametrize("row", ROWS, ids=[r["example"] for r in ROWS])
def test_example_compatibility(row):
    with contextlib.redirect_stdout(io.StringIO()):
        if row["status"] == "unsupported material law":
            with pytest.raises(ValueError, match="Unsupported material law"):
                get_parameters(row["args"])
            return
        assert row["status"] == "supported smoke-tested", row
        P, G, L = get_parameters(row["args"])
        P.dt = min(P.dt, 1e-6)
        run(P, G, L, max_steps=5, output=False)
        for b in L.materials:
            b.check_finite(P)
        for key in ("m", "q", "fi", "fe", "phi", "dphi"):
            assert np.isfinite(getattr(G, key)).all()


def test_full_ball():
    P, G, L = get_parameters(["ball"])
    run(P, G, L, output=False)
    assert P.t == P.t_f
    assert P.tstep == 1000
    assert len(L.materials[0]) == 500
