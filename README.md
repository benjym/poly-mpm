# Polydisperse material point method

[![CI](https://github.com/benjym/poly-mpm/actions/workflows/ci.yml/badge.svg)](https://github.com/benjym/poly-mpm/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/benjym/poly-mpm/graph/badge.svg)](https://codecov.io/gh/benjym/poly-mpm)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![PyPI](https://img.shields.io/pypi/v/poly-mpm.svg)](https://pypi.org/project/poly-mpm/)
[![Python 3.14](https://img.shields.io/badge/python-3.14-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A small 2D material point method solver for solid and viscous/granular flow,
including grain-size segregation. Particle state uses NumPy arrays, with one
batch per material and a shared regular Cartesian grid.

## Install from PyPI

With Python 3.14 or later, install in a virtual environment and run an example:

```sh
python -m pip install poly-mpm
poly-mpm ball
```

For development, tests, and benchmark modules, use a source checkout below.

## Source checkout

Use Git and a current version of [uv](https://docs.astral.sh/uv/getting-started/installation/).
uv can install the pinned Python 3.14.10 interpreter. From a terminal:

```sh
git clone https://github.com/benjym/poly-mpm.git
cd poly-mpm
uv sync --extra dev --locked
uv run --locked python MPM.py ball
```

For an existing checkout, start at `uv sync` from the repository root.
The installed command `uv run --locked poly-mpm ball` runs the same example.
For quadratic B-splines, APIC, and dynamic timestepping, use
`uv run --locked python MPM.py modern_ball`.

`ball` simulates one second with 500 elastic particles and 1,000 steps. A Rich
progress bar shows simulation time, elapsed time, estimated time remaining, and
when output is being saved. Figures are saved under **`im/ball/elastic/Continuum/`** and
**`im/ball/elastic/MP/`**. No plot window opens. Run time depends on your machine
and plotting; repeating a run can overwrite its output.

For installation without uv, Windows instructions, and a three-step preview with
its own output folder, see [getting started](https://github.com/benjym/poly-mpm/blob/main/docs/getting_started.rst).
`uv.lock` records the tested dependency set. The pip alternative resolves compatible
versions rather than reproducing that exact environment.

## Run a short preview

The CLI accepts an example name and that example's positional arguments. It has
no generic `--max-steps` option; use Python for bounded runs:

```python
from initialise import get_parameters
from MPM import run

P, G, L = get_parameters(['ball'])
run(P, G, L, max_steps=100, output=False)
print(P.t, P.tstep, L.materials[0].x.shape)
```

Save this as `preview.py` in the checkout and run `uv run --locked python preview.py`.
`output=False` disables all callbacks, so this preview creates no plots or saved
histories. Some examples request millions of steps; preview them first.

For segregation, `bi_seg_test 22 2 11` means 22 degrees, two size classes, and
11 vertical grid nodes. See the [complete example matrix](https://github.com/benjym/poly-mpm/blob/main/docs/array-conversion.md#example-compatibility)
for required arguments and support status. A smoke-tested input is not a guarantee
of stability or physical calibration at every setting.

For the square drum, run `uv run --locked python MPM.py bi_square_drum 11 2`:
`11` is the number of grid nodes per side (at least 3), and `2` is the positive
Froude number. Both arguments are required.

Run `uv run --locked python MPM.py --help` to list examples, or
`uv run --locked python MPM.py bi_dam_break --help` for an example's argument
order, types, and choices. Missing, invalid, and extra arguments produce a usage
error before simulation setup. `silo` takes only the vertical grid node count,
for example `silo 11`; its two grain size classes are fixed in the input.

## User guide

- [Experimental rigid-grain solver](docs/granular-prototype.md): coupled packing/friction prototype for static zones and free surfaces; run `granular_column` or `granular_collapse`.
- [Getting started](https://github.com/benjym/poly-mpm/blob/main/docs/getting_started.rst): setup, first run, previews, and examples.
- [Create a simulation](https://github.com/benjym/poly-mpm/blob/main/docs/tutorial.rst): a complete runnable input, units, geometry,
  boundaries, material parameters, and segregation settings.
- [Numerical methods](https://github.com/benjym/poly-mpm/blob/main/docs/numerical-foundations.md): deformation/volume evolution,
  interpolation, APIC, timestep estimates, diagnostics, and compatibility mode.
- [Finite-domain quadrature](https://github.com/benjym/poly-mpm/blob/main/docs/quadrature.md): particle geometry, cell crossing, and integration controls.
- [Verification](https://github.com/benjym/poly-mpm/blob/main/docs/verification.md): equilibrium, transport, finite strain, and error studies.
- [Implementation checklist](https://github.com/benjym/poly-mpm/blob/main/IMPLEMENTATION_TODO.md): completed acceptance checks and deferred work.
- [Material selection](https://github.com/benjym/poly-mpm/blob/main/docs/materials.rst): supported laws and required parameters.
- [Output and API](https://github.com/benjym/poly-mpm/blob/main/docs/outputs.rst): output flags, file formats, callbacks, and arrays.
- [Troubleshooting](https://github.com/benjym/poly-mpm/blob/main/docs/troubleshooting.rst): installation, missing output, numerical
  errors, and current limitations.
- [Constitutive models](https://github.com/benjym/poly-mpm/blob/main/docs/constitutive-models.md): equations, parameters, and model limits.
- [Particle arrays and examples](https://github.com/benjym/poly-mpm/blob/main/docs/array-conversion.md): state shapes, input interfaces, and support status.
- [Segregation transport](https://github.com/benjym/poly-mpm/blob/main/docs/segregation-conservation.md): equations and verification.
- [Build the HTML manual](https://github.com/benjym/poly-mpm/blob/main/docs/README.md): optional documentation dependencies and checks.

The solver supports `rigid`, `elastic`, `viscous`, `von_mises`, `pouliquen`,
`bingham`, `dp`, `dp_rate`, and `linear_mu`. Retired or absent laws fail explicitly.
Geometry is 2D despite three-component vectors and 3-by-3 stress tensors. Updated
mechanics evolves deformation, volume, and physical density. Linear/quadratic
interpolation and FLIP/PIC/APIC transfers are selectable; `modern_ball` demonstrates
the combined configuration. Mechanics still uses explicit stepping and shared
nodal momentum, with no independent body contact, implicit solve, or general
checkpoint/restart. Finite-domain integration is available in source checkouts through
`quadrature='domain'` with declared particle geometry; the published `1.0.0`
package uses point quadrature. The explicit `mechanics='legacy'` compatibility configuration uses fixed-volume,
incremental mechanics; the user guide describes the default configuration.

## Tests and benchmarks

From the checkout root:

```sh
uv run --locked --extra dev black --check .
uv run --locked --extra dev python -m pytest -q -W error::RuntimeWarning --cov --cov-report=term-missing --cov-report=xml
uv run --locked python -m benchmarks.audit_examples --output /tmp/poly-mpm-examples.json
```

Run `uv run --locked --extra dev black .` to format active Python code.
The `legacy/` source archive is excluded to preserve its provenance. Coverage
includes installed solver modules and inputs, including plotting and diagnostics;
HTML coverage is available with `--cov-report=html`.

GitHub Actions checks formatting, runs the coverage suite, builds the strict HTML
manual, and validates wheel/sdist metadata plus an installed-wheel smoke test.
The test job uploads `coverage.xml` with `CODECOV_TOKEN`; fork pull requests skip
the upload because repository secrets are unavailable. Coverage is informational
until a useful baseline and thresholds have been agreed.
See the [release checklist](https://github.com/benjym/poly-mpm/blob/main/docs/releasing.md) for version-tag and manual publishing instructions.

Choose another writable report path on Windows. The audit runs five bounded
steps per supported example. See [benchmark instructions](https://github.com/benjym/poly-mpm/blob/main/benchmarks/README.md)
for timing and conservation reports. Their default output paths overwrite
committed reports; use separate paths for your own measurements.

The [verification guide](https://github.com/benjym/poly-mpm/blob/main/docs/verification.md) describes analytical references,
refinement controls, measured errors, and scope. Benchmark reports record their
own configuration and environment; performance depends on material, particle
count, diagnostics, and output frequency.

## Project

Originally authored by Benjy Marks, University of Sydney. NumPy, SciPy, Astropy,
and Matplotlib provide the numerical and plotting environment.

Licensed under the [GNU GPL v3](https://github.com/benjym/poly-mpm/blob/main/LICENSE).

### Experimental GPU port

An optional Apple silicon/SlangPy backend runs the granular rigid-grain
timestep on the GPU in FP32 (`P.backend = "metal"`). Its hardware gates passed
on an Apple M4; see [the GPU notes](docs/gpu.md) for supported settings,
results and the [validation command](docs/gpu.md). Matched runs on an M4 were
3 to 26 times faster than the CPU. CPU execution remains the default.
