# Array conversion benchmark results

This is a historical performance report, measured before the segregation correction
and constitutive-model restoration. For installation and a first simulation, start
with the [user guide](../docs/getting_started.rst).

Five-run median solver timings after one warmup, on identical initial states.
Initialisation, state copying, plotting, logging, and file output are excluded.
The scalar baseline is the corrected test reference derived from the original
solver, using scalar particle kernels on the same state. It is not a direct
benchmark of unmodified master. Numerical agreement is asserted before recording
each result. Grid and segregation integrators are shared between implementations.

Environment: Python 3.12.8, macOS-26.6.2-x86_64-i386-64bit.
NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2, Astropy 8.0.2.
BLAS and OpenMP threads were limited to one.

| Particles | Model | Grain classes | Segregation | Scalar timestep (ms) | Array timestep (ms) | Speedup |
|---:|---|---:|---|---:|---:|---:|
| 1,000 | elastic | 1 | no | 187.34 | 4.13 | 45.4× |
| 1,000 | viscous | 1 | no | 205.04 | 4.25 | 48.2× |
| 1,000 | von_mises | 1 | no | 202.39 | 4.43 | 45.7× |
| 1,000 | pouliquen | 1 | no | 197.76 | 4.20 | 47.0× |
| 1,000 | pouliquen | 2 | yes | 198.39 | 7.11 | 27.9× |
| 1,000 | pouliquen | 8 | yes | 198.48 | 8.06 | 24.6× |
| 1,000 | rigid | 1 | no | 83.83 | 3.18 | 26.3× |
| 10,000 | elastic | 1 | no | 1865.37 | 22.00 | 84.8× |
| 10,000 | viscous | 1 | no | 1879.44 | 21.78 | 86.3× |
| 10,000 | von_mises | 1 | no | 2008.02 | 24.46 | 82.1× |
| 10,000 | pouliquen | 1 | no | 1947.26 | 22.35 | 87.1× |
| 10,000 | pouliquen | 2 | yes | 1902.02 | 31.30 | 60.8× |
| 10,000 | pouliquen | 8 | yes | 1910.08 | 36.16 | 52.8× |
| 10,000 | rigid | 1 | no | 805.14 | 11.14 | 72.3× |
| 100,000 | elastic | 1 | no | 18273.49 | 232.53 | 78.6× |
| 100,000 | viscous | 1 | no | 18684.74 | 226.81 | 82.4× |
| 100,000 | von_mises | 1 | no | 20113.91 | 226.72 | 88.7× |
| 100,000 | pouliquen | 1 | no | 19216.02 | 252.77 | 76.0× |
| 100,000 | pouliquen | 2 | yes | 19129.48 | 272.95 | 70.1× |
| 100,000 | pouliquen | 8 | yes | 19451.08 | 348.16 | 55.9× |
| 100,000 | rigid | 1 | no | 8161.83 | 110.76 | 73.7× |

The array timestep is faster in every measured case. Elastic speedup is
approximately 85× at 10,000 particles and 79× at 100,000 particles.
These figures describe this environment and these fixtures, not all simulations.

## Stage timings and memory

The [raw report](results.json) separates transfer, constitutive/diagnostic,
force, grid, segregation, and movement costs. At 100,000 particles:

| Model | Grain classes | Persistent particle state (MiB) | Peak timestep allocations (MiB) |
|---|---:|---:|---:|
| elastic | 1 | 67.1 | 49.4 |
| viscous | 1 | 67.1 | 64.6 |
| von_mises | 1 | 67.1 | 59.4 |
| pouliquen | 1 | 67.1 | 66.9 |
| pouliquen | 2 | 67.9 | 66.9 |
| pouliquen | 8 | 72.5 | 96.2 |
| rigid | 1 | 67.1 | 25.8 |

Peak timestep allocations are measured by tracemalloc after state construction;
they exclude the already allocated particle/grid state and are not process RSS.
Persistent state excludes scratch/connectivity arrays. Memory measurements run
separately from timing. No memory reduction relative to the scalar reference is
claimed: both comparison implementations use the same state layout.

## Remaining work after conversion

The largest array costs at 100,000 particles are constitutive updating and
diagnostic deposition, followed by force mapping and movement. Allocation reuse
or compiled kernels can be evaluated against these stages in later work.
No Numba/GPU path or new integration/constitutive scheme was added here.

The original performance table above predates the segregation correction and
Python upgrade. It remains a historical measurement. Current 1,000-particle
Pouliquen timings (three-run medians, Python 3.14.10) are recorded in
[the correction performance report](segregation_correction_performance.json).

The [current conservation report](segregation_drift.json) and
[before-correction report](segregation_drift_before.json) separately measure
composition sums, bounds, and species mass without clipping or normalisation.
See [the source and correction](../docs/segregation-conservation.md) and
[conversion notes](../docs/array-conversion.md) for the derivation and limitations.

Install the dev dependencies with `uv sync --extra dev --locked`, activate the
environment (`source .venv/bin/activate` on macOS/Linux), then reproduce from the
repository root with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m benchmarks.compare --sizes 1000 --repeats 3 --output /tmp/poly-mpm-timings.json
python -m benchmarks.audit_examples --output /tmp/poly-mpm-examples.json
python -m benchmarks.segregation_drift --output /tmp/poly-mpm-segregation.json
```

These commands write fresh reports separately from the committed measurements.
Choose another writable path on Windows. Omitting `--sizes`/`--repeats` runs the
full timing suite, including slow 100,000-particle scalar cases. Omitting `--output`
overwrites the corresponding report in `benchmarks/`. Runtime, tests, and benchmark
modules must be available from an activated source-checkout environment.
