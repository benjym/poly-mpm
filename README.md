# Polydisperse material point method

A small 2D material point method solver for solid and viscous/granular flow,
including grain-size segregation. Particle state uses NumPy arrays, with one
batch per material and a shared regular Cartesian grid.

## Quick start

Use Git and a current version of [uv](https://docs.astral.sh/uv/getting-started/installation/).
uv can install the pinned Python 3.14.10 interpreter. From a terminal:

```sh
git clone https://github.com/benjym/poly-mpm.git
cd poly-mpm
uv sync --extra dev --locked
uv run --locked python MPM.py ball
```

For an existing checkout, start at `uv sync` from the repository root. Use the
branch containing this README; older branches may contain the scalar solver.
The installed command `uv run --locked poly-mpm ball` runs the same example.

`ball` simulates one second with 500 elastic particles and 1,000 steps. It prints
progress and saves figures under **`im/ball/elastic/Continuum/`** and
**`im/ball/elastic/MP/`**. No plot window opens. Run time depends on your machine
and plotting; repeating a run can overwrite its output.

For installation without uv, Windows instructions, and a three-step preview with
its own output folder, see [getting started](docs/getting_started.rst).
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
histories. Some original examples request millions of steps; preview them first.

For segregation, `bi_seg_test 22 2 11` means 22 degrees, two size classes, and
11 vertical grid nodes. See the [complete example matrix](docs/array-conversion.md#example-compatibility)
for required arguments and support status. A smoke-tested input is not a guarantee
of stability or physical calibration at every original setting.

## User guide

- [Getting started](docs/getting_started.rst): setup, first run, previews, and examples.
- [Create a simulation](docs/tutorial.rst): a complete runnable input, units, geometry,
  boundaries, material parameters, and segregation settings.
- [Material selection](docs/materials.rst): supported laws and required parameters.
- [Output and API](docs/outputs.rst): output flags, file formats, callbacks, and arrays.
- [Troubleshooting](docs/troubleshooting.rst): installation, missing output, numerical
  errors, and current limitations.
- [Constitutive audit](docs/constitutive-models.md): restored equations and retired laws.
- [Array reference and examples](docs/array-conversion.md): shapes, migration, and numerical scope.
- [Segregation conservation](docs/segregation-conservation.md): equations and verification.
- [Build the HTML manual](docs/README.md): optional documentation dependencies and checks.

The solver supports `rigid`, `elastic`, `viscous`, `von_mises`, `pouliquen`,
`bingham`, `dp`, `dp_rate`, and `linear_mu`. Retired or absent laws fail explicitly.
Geometry is 2D despite three-component vectors and 3-by-3 stress tensors. Current
mechanics use bilinear interpolation, explicit stepping, and shared nodal momentum;
there is no independent body contact, deformation-gradient volume update, implicit
solve, or general checkpoint/restart. The original five laws retain legacy
conventions; the four restored laws have separate analytical tests and corrections.

## Tests and benchmarks

From the checkout root:

```sh
uv run --locked --extra dev python -m pytest -q -W error::RuntimeWarning
uv run --locked python -m benchmarks.audit_examples --output /tmp/poly-mpm-examples.json
```

Choose another writable report path on Windows. The audit runs five bounded
steps per supported example. See [benchmark instructions](benchmarks/README.md)
for timing and conservation reports. Their default output paths overwrite
committed reports; use separate paths for your own measurements.

The [performance report](benchmarks/README.md) is a historical Python 3.12.8
comparison, measured before the runtime upgrade and later numerical corrections.
Its timings have not been relabelled as current results. Tests and example checks
run separately on Python 3.14 with the locked dependencies.

## Project

Originally authored by Benjy Marks, University of Sydney. NumPy, SciPy, Astropy,
and Matplotlib provide the numerical and plotting environment.

Licensed under the [GNU GPL v3](LICENSE).
