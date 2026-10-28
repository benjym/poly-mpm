"""Example argument definitions shared by the CLI and Python loader."""

import argparse
import math
from pkgutil import iter_modules


class UsageParser(argparse.ArgumentParser):
    def error(self, message):
        # Library callers get an exception; MPM.main renders it on stderr.
        raise ValueError(f"{self.format_usage().strip()}\n{message}")


def integer_at_least(minimum):
    def convert(value):
        try:
            number = int(value)
        except ValueError as exc:
            raise argparse.ArgumentTypeError("must be an integer") from exc
        if number < minimum:
            raise argparse.ArgumentTypeError(f"must be at least {minimum}")
        return number

    return convert


def finite_float(value):
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if not math.isfinite(number):
        raise argparse.ArgumentTypeError("must be finite")
    return number


def positive_float(value):
    number = finite_float(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be positive")
    return number


GRID = ("grid_nodes", integer_at_least(3), "vertical grid nodes (integer, at least 3)")
CHANNEL_GRID = ("grid_nodes", integer_at_least(2), "vertical grid nodes (integer, at least 2)")
SIZES = ("size_classes", integer_at_least(2), "grain size classes (integer, at least 2)")
ANGLE = ("slope_degrees", finite_float, "slope angle in degrees (finite number)")
FROUDE = ("Froude_number", positive_float, "Froude number (finite and positive)")
ARRANGEMENT = ("arrangement", ("top", "bot", "mix"), "small grains on top, on bottom, or mixed")

# Inputs without an entry accept no positional arguments. Register new arguments
# here so help and validation stay consistent with Params(args).
ARGUMENTS = {
    "bi_dam_break": (CHANNEL_GRID, ARRANGEMENT, ANGLE, SIZES),
    "bi_double_insets": (CHANNEL_GRID,),
    "bi_drum": (GRID, FROUDE),
    "bi_seg_test": (ANGLE, SIZES, CHANNEL_GRID),
    "bi_square_drum": (GRID, FROUDE),
    "bi_vertical_chute": (SIZES, CHANNEL_GRID),
    "francois": (GRID,),
    "poly_seg_test": (ANGLE,),
    "silo": (GRID,),
}


def validate_arguments(params):
    """Validate before constructing a simulation; --help exits without setup."""
    examples = sorted(name for _, name, is_package in iter_modules(__path__) if not is_package)
    parser = UsageParser(
        prog="python MPM.py",
        description="Run a Material Point Method example.",
        epilog="Examples: " + ", ".join(examples) + ". Use <example> --help for its arguments.",
    )
    parser.add_argument("example", choices=examples, metavar="<example>")
    parser.add_argument("arguments", nargs=argparse.REMAINDER, help="example-specific positional arguments")
    parsed = parser.parse_args(params)
    parser = UsageParser(
        prog=f"python MPM.py {parsed.example}",
        description=f"Run the {parsed.example} example.",
        epilog="For bounded previews, use MPM.run(..., max_steps=100, output=False) in Python.",
    )
    for name, converter, help_text in ARGUMENTS.get(parsed.example, ()):
        options = {"choices": converter} if isinstance(converter, tuple) else {"type": converter}
        parser.add_argument(name, help=help_text, **options)
    parser.parse_args(parsed.arguments)
