import pytest

from conftest import state
from MPM import run


@pytest.mark.parametrize(
    "options, status, percentage",
    [
        ({"t_f": 2e-5}, "Complete", "100.0%"),
        ({"max_steps": 1}, "Stopped", "0.1%"),
        ({"t_f": 0}, "Complete", "100.0%"),
    ],
)
def test_progress_reports_actual_simulation_time(options, status, percentage, capsys):
    options = options.copy()
    max_steps = options.pop("max_steps", None)
    P, G, L = state(**options)
    result = run(P, G, L, max_steps=max_steps, output=False, progress=True)
    assert result == (P, G, L)
    text = capsys.readouterr().out
    assert status in text
    assert percentage in text
    assert f"t = {P.t:.6g} s" in text


def test_progress_disabled_is_silent(capsys):
    run(*state(), max_steps=1, output=False)
    assert capsys.readouterr().out == ""
