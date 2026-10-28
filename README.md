# Polydisperse material point method

A basic two-dimensional material point method with grain-size segregation during
granular flow. Particle state is stored in NumPy arrays, with one collection per
material. The supported timestep has no Python loops over individual particles.

## Install and run

Python 3.14 is the supported baseline; `.python-version` pins Python 3.14.10.
`uv.lock` records the complete tested runtime and test dependencies. With a current
version of uv, from a checkout:

```sh
uv sync --extra dev --locked
source .venv/bin/activate
python MPM.py ball
# Or use the installed command:
poly-mpm ball
```

For a pip-based setup, use `python3.14 -m venv .venv`, activate it, and run
`python -m pip install -e '.[dev]'`. This follows the supported minimum versions;
use the lockfile workflow for the exact tested dependency set.

The default output is `im/<example>/<first-material-law>/`, unless an input defines
`supername` or `save_dir`. Examples can run for a long time at their original
settings; the compatibility checks use bounded runs.

Supported material laws are `rigid`, `elastic`, `viscous`, `von_mises`, and
`pouliquen`. Rigid particles follow prescribed velocities; no new contact model is
provided. Experimental and absent laws fail explicitly at initialisation.

```sh
python MPM.py bi_seg_test 22 2 11
python MPM.py bi_seg_test 22 8 11
```

See [the example matrix](docs/array-conversion.md#example-compatibility) for required
arguments and the status of every example. A smoke test establishes compatibility,
not long-run physical validation or a guarantee of stability for its original dt.

## Tests and benchmarks

```sh
python -m pytest -q -W error::RuntimeWarning
python -m benchmarks.audit_examples
python -m benchmarks.segregation_drift
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m benchmarks.compare
```

The benchmark compares identical states against an independent corrected scalar
reference derived from the old solver, rather than timing the broken conversion
or claiming numerical equivalence to every bug in `master`. It reports five-run
medians, stage timings, environment versions, particle-state storage, and traced
peak allocations during an array timestep. The complete default run includes
100,000-particle scalar cases and takes several minutes.

The original performance report was measured on Python 3.12.8 before the runtime
upgrade and segregation correction. It remains a historical comparison; those timings have not been relabelled
as Python 3.14 measurements. The full test suite and example smoke tests also pass
on Python 3.14 with the locked dependencies.

Results: [benchmark report](benchmarks/README.md),
[raw timings](benchmarks/results.json),
[example audit](benchmarks/examples.json), and
[segregation drift](benchmarks/segregation_drift.json).

## Array interfaces and migration

`initialise.get_parameters(args)` returns `(P, G, L)` as before. `L.materials`
contains array-backed `Particles` collections; `P.S` is normalised to a material
list. The old `L.S[p][i]` object interface has been removed. Input geometry may
still be generated with lists; it is converted once when particle state is built.

Use `MPM.step(P, G, L)` for one timestep without advancing the clock or producing
output. Use `MPM.run(P, G, L, max_steps=100, output=False)` for a bounded run.

All output callbacks receive `(L, P, G)`. Legacy plotting flags remain supported,
as do `segregate_grid` and the misspelled `after_ever_timestep` aliases. Custom
callbacks replace default output for their stage. See
[conversion details](docs/array-conversion.md) for shapes, conventions, intentional
numerical corrections, and remaining limitations.

Segregation now transports conservative species-mass exchanges, with coupled
mixture fluxes and a bounded, inventory-aware return transfer. The 100- and
1,000-step checks preserve species mass and fraction sums to roundoff without
clipping or normalisation. See [the source and correction of the segregation
drift](docs/segregation-conservation.md) for the derivation, before/after results,
and numerical limitations.

## Project

Originally authored by Benjy Marks, University of Sydney. NumPy, SciPy, Astropy,
and Matplotlib provide the numerical and plotting environment.

Licensed under the [GNU GPL v3](LICENSE).
