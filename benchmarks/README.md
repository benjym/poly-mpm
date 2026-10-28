# Running benchmarks and verification

Run these modules from a source checkout with runtime dependencies installed.
Use `uv sync --locked --extra dev` and prefix commands with `uv run --locked`,
or activate the environment first. Reports identify their configuration and
measurement environment; timings from different environments are not directly
comparable. Use explicit output paths to keep your measurements separate from
committed reports.

## Analytical solver verification

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
