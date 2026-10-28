"""Usage errors must be actionable before simulation setup starts."""

import pytest
from MPM import main
from initialise import get_parameters
from benchmarks.audit_examples import ARGUMENTS


@pytest.mark.parametrize(
    "params",
    [[name, *args[:count]] for name, args in ARGUMENTS.items() for count in range(len(args))],
)
def test_missing_example_arguments(params, capsys):
    assert main(params) == 2
    captured = capsys.readouterr()
    assert f"usage: python MPM.py {params[0]}" in captured.err
    assert "required" in captured.err
    assert not captured.out


@pytest.mark.parametrize(
    "params, message",
    [
        ([], "required"),
        (["not_an_example"], "invalid choice"),
        (["inputs/ball.py"], "invalid choice"),
        (["ball", "unexpected"], "unrecognized arguments"),
        (["ball", "--max-steps", "1"], "unrecognized arguments"),
        (["bi_dam_break", "11", "typo", "0", "2"], "choose from 'top', 'bot', 'mix'"),
        (["bi_dam_break", "11", "mix", "nan", "2"], "slope_degrees: must be finite"),
        (["bi_dam_break", "11", "mix", "0", "1"], "size_classes: must be at least 2"),
        (["bi_drum", "11", "0"], "Froude_number: must be positive"),
        (["bi_seg_test", "22", "1", "11"], "size_classes: must be at least 2"),
        (["bi_seg_test", "22", "2", "1"], "grid_nodes: must be at least 2"),
        (["bi_vertical_chute", "0", "11"], "size_classes: must be at least 2"),
        (["bi_double_insets", "1"], "grid_nodes: must be at least 2"),
        (["francois", "2"], "grid_nodes: must be at least 3"),
        (["poly_seg_test", "inf"], "slope_degrees: must be finite"),
        (["silo", "2"], "grid_nodes: must be at least 3"),
        (["silo", "11", "2"], "unrecognized arguments: 2"),
    ],
)
def test_usage_errors(params, message, capsys):
    assert main(params) == 2
    captured = capsys.readouterr()
    assert "usage: python MPM.py" in captured.err
    assert message in captured.err
    assert "Traceback" not in captured.err
    assert not captured.out


@pytest.mark.parametrize("example", [None, "ball", *ARGUMENTS])
def test_help_without_simulation_setup(example, monkeypatch, capsys):
    def unexpected_setup(*args, **kwargs):
        pytest.fail("help must not construct a simulation")

    monkeypatch.setattr("initialise.prepare", unexpected_setup)
    with pytest.raises(SystemExit) as exc:
        main([example, "--help"] if example else ["--help"])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "usage: python MPM.py" in captured.out
    assert not captured.err
    if example is None:
        assert "bi_dam_break" in captured.out
    elif example in ARGUMENTS:
        assert "positional arguments:" in captured.out


def test_loader_preserves_usage_errors_for_library_callers():
    with pytest.raises(ValueError, match="size_classes: must be at least 2"):
        get_parameters(["bi_vertical_chute", "1", "11"])
