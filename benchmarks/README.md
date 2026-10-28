# Running benchmarks and verification

Run these modules from a source checkout with runtime dependencies installed.
Use `uv sync --locked --extra dev` and prefix commands with `uv run --locked`,
or activate the environment first. Reports identify their configuration and
measurement environment; timings from different environments are not directly
comparable. Use explicit output paths to keep your measurements separate from
committed reports.

## Analytical solver verification

For the opt-in rigid-grain branch, run
`uv run --locked python -m benchmarks.granular_prototype --output /tmp/granular-prototype.json`.
The committed [report](granular_prototype.json) checks static support and
coarse-grid collapse timestep sensitivity. See the [prototype guide](../docs/granular-prototype.md)
for its constitutive assumptions, size limits and remaining validation work.
This changes the physical model and does not measure an equivalent-behavior
speedup against explicit elastic MPM.

```sh
uv run --locked python -m benchmarks.broader_verification --output /tmp/poly-mpm-broader.json --figure /tmp/poly-mpm-broader.png
uv run --locked python -m benchmarks.verification --output /tmp/poly-mpm-wave.json
```

`broader_verification` checks hydrostatic equilibrium, self-equilibrated stress
translation, finite-strain manufactured dynamics, frozen-state force integration,
fixed-cloud transfers and prescribed affine deformation. Its grid, timestep,
particle-quadrature and transfer studies have separate controls. See the
[verification guide](../docs/verification.md) for equations, metrics and limits.

`verification` checks periodic standing waves with linear/FLIP and quadratic
FLIP/PIC/APIC configurations, plus complete `modern_ball` and `dp_rate_unit_test`
runs. Wave errors can be non-monotone. Neither report establishes universal
convergence or experimental validity.

Committed evidence: [broader_verification.json](broader_verification.json) and
[verification.json](verification.json), with source fingerprints and environment
metadata. Their plots are scientific summaries; JSON retains all measured levels.

## Particle quadrature and cell crossing

```sh
uv run --locked python -m benchmarks.quadrature_verification --output /tmp/poly-mpm-quadrature.json --figure /tmp/poly-mpm-quadrature.png
```

This compares point sampling and finite-domain basis integration on identical
grids, particle clouds and timestep sequences, including a frozen-cloud cell-
phase sweep. [quadrature_verification.json](quadrature_verification.json) records
all levels and source hashes; the [quadrature guide](../docs/quadrature.md)
explains results, domain assumptions, computational cost and remaining limits.
The point/domain choice is independent of the cloud's particle-quadrature density.

## Example compatibility

```sh
uv run --locked python -m benchmarks.audit_examples --output /tmp/poly-mpm-examples.json
```

The audit runs five bounded steps per supported input, with an initial timestep
at most `1e-6`, configured dynamic stepping, and output disabled. It reports
unsupported laws explicitly. See the [example matrix](../docs/array-conversion.md#example-compatibility)
for positional arguments. Passing a smoke check does not establish stability at
the example's full duration or default timestep.

## Segregation conservation

```sh
uv run --locked python -m benchmarks.segregation_drift --output /tmp/poly-mpm-segregation.json
```

This measures fraction sums, bounds and mass-weighted species inventories with
normalisation disabled, including multigrain runs. The
[transport guide](../docs/segregation-conservation.md) explains the flux and bounded
return transfer. The committed report is [segregation_drift.json](segregation_drift.json).
Conservation checks do not calibrate the physical segregation coefficients.

## Performance profiling

```sh
uv run --locked python -m benchmarks.compare --sizes 1000 --repeats 3 --output /tmp/poly-mpm-timings.json
```

This driver compares vectorised kernels with scalar reference kernels using
explicit `mechanics='legacy'`, diagnostics disabled, and identical state layouts.
It is a compatibility/performance harness, not a measurement of the default
neo-Hookean/APIC solver. Omitting `--sizes` and `--repeats` includes slow scalar
runs at up to 100,000 particles. Stage timing and allocation measurements are
separate; tracemalloc peak allocations exclude persistent state and are not RSS.

[results.json](results.json) records this reference configuration on Python 3.12.8;
[segregation_correction_performance.json](segregation_correction_performance.json)
records a separate small segregation configuration. Use fresh measurements for
performance decisions about the solver you intend to run.

On Windows, choose writable output paths instead of `/tmp`. Generated reports
should state grid, material, transfer, timestep, diagnostics, and output settings
so numerical work and plotting costs can be interpreted separately.

## Updated-solver hot-path comparison

`solver_performance` measures the default updated mechanics using `ball`,
`modern_ball`, a deterministic 10,000-particle quadratic/APIC elastic cloud, and
1,000 particles with three grain classes and segregation. Setup, plotting and
progress are outside the timed work; diagnostics stay enabled. Each case has one
warmup and three timed runs, with the median reported.

Before changing the solver, retain its numerical snapshots:

```sh
uv run --locked python -m benchmarks.solver_performance --output /tmp/mpm-before.json --snapshot-dir /tmp/mpm-before-state
```

After the change, compare all persistent particle arrays, nodal arrays, final
clock values, and the complete diagnostic history:

```sh
uv run --locked python -m benchmarks.solver_performance --output /tmp/mpm-after.json --snapshot-dir /tmp/mpm-after-state --reference /tmp/mpm-before-state
```

[solver_performance.json](solver_performance.json) records the first optimisation
pass, including sample timings, environment and source fingerprints. Frozen source
snapshots excluded concurrent CLI/quadrature work from this comparison.

| Case | Steps | Before (s) | After (s) | Speedup |
| --- | ---: | ---: | ---: | ---: |
| ball, 500 particles | 200 | 1.415 | 1.073 | 1.32x |
| modern_ball, 500 particles | 200 | 3.015 | 1.796 | 1.68x |
| elastic APIC, 10,000 particles | 40 | 7.993 | 5.942 | 1.35x |
| segregation APIC, 1,000 particles | 40 | 1.293 | 1.166 | 1.11x |

All compared numerical values were bit-for-bit identical. These are bounded CPU
solver measurements, not full-duration timings or expected gains for every input.
The changes replace unsmoothed masked gradients with equivalent explicit masks,
reuse the same APIC moment/pseudoinverse until stencil changes, and share central
stability estimates while the state is unchanged. The pseudoinverse tolerance,
custom material timestep callbacks, smoothing path, constitutive models and
output frequency are retained.
