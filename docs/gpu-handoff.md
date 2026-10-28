# Handoff: continue the Metal GPU port on Apple silicon

This note is for an agent working on the target Mac. Start by reading this file
and [the GPU implementation plan and acceptance gates](gpu.md). The user has
authorized implementing that plan, keeping GPU execution opt-in, and committing
and pushing the work. Continue implementation after each hardware gate passes;
the task is larger than running the first validation command.

## Checkpoint and scope

- Branch: `codex/metal-backend`.
- Runtime/primitives checkpoint: `4dd4ebf`.
- CPU science baseline: `3d72a69`, before the separate density work.
- Target: Apple silicon only, macOS 26+, native CPython 3.14, Metal via
  SlangPy **0.43.1**. The Intel development Mac is not a job target.
- Keep this backend inside the existing package, under `mpm_gpu`, with the same
  inputs, outputs, plots, Rich progress and save frequency. Segregation will be
  added later through this package; implement mechanics first.
- The user explicitly chose validated FP32 and fully device-side nonlinear
  control, including small QR/SVD and GMRES calculations.
- Do not incorporate the density fix or change the physical model as part of
  this port. Inspect the local checkout before switching branches or editing;
  preserve unrelated changes. The other agent's density changes were left in
  the original checkout and are not included in this branch.

**The MPM GPU timestep does not exist yet.** The current `ValidationPending`
exceptions in `initialise.prepare`, `MPM.step` and `MPM.run` are deliberate.
They prevent selecting Metal and accidentally running CPU mechanics. Replace
these guards only when the corresponding resident implementation exists and
passes its gates. There is no GPU MPM speedup measurement yet.

## First action on the target Mac

Fetch and check out `origin/codex/metal-backend` in a suitable clean checkout.
Then run:

```sh
uname -m
xcrun --find metal
xcrun metal --version
uv sync --locked --extra gpu --extra dev --python 3.14
uv run --locked --extra gpu --extra dev python -m mpm_gpu.validate --stage runtime --benchmark --output metal-validation.json
```

`uname -m` must report `arm64`. A successful `xcrun --find metal` alone is
insufficient: Xcode can install a launcher without the Metal Toolchain
component. If prerequisites are missing, report the exact failure and resolve
setup within the user's authorization; do not substitute a different backend.

Read the JSON even if the command fails. Fix actual Metal compilation, binding,
synchronization, indirect dispatch, atomics or timestamp issues before relying
on those primitives. Re-run failing cases with `--debug` if useful; omit debug
layers for benchmark timings. Preserve failed reports alongside corrected ones.
Save a source-fingerprinted hardware report under `benchmarks/` when the gate
passes, recording the actual device and environment.

After a runtime pass, run the explicit device test and CPU regression suite:

```sh
POLY_MPM_TEST_METAL=1 uv run --locked --extra gpu --extra dev python -m pytest tests/test_gpu_hardware.py -q
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 uv run --locked --extra gpu --extra dev python -m pytest -q
```

`--stage runtime` should exit 0 with status `passed`. At this checkpoint,
`--stage all` runs the runtime gate but then exits 3 and labels mechanics,
solver and runner `not_implemented`. That is an implementation status, not a
hardware failure or permission to mark the whole port validated.

## What is already written

| Files | Purpose |
| --- | --- |
| `mpm_gpu/runtime.py` | Lazy SlangPy import, explicit Metal device, compiler/platform checks, kernel cache, scalar buffers, dispatch/barriers, indirect dispatch, timestamps and transfer accounting |
| `mpm_gpu/shaders/primitives.slang` | Nine entry points: fill, histogram, tiled scan, scan offsets, scatter, reduction, device dispatch arguments, indirect marker and ABI probe |
| `mpm_gpu/primitives.py` | Hierarchical scans/reductions and histogram → scan → scatter; intermediate buffers remain on device |
| `mpm_gpu/topology.py` | Static canonical Q1 neighbors, transpose incidence, conflict colors, split cell/fraction coordinates |
| `mpm_gpu/config.py` | Backend defaults, supported configuration validation and temporary runner guard |
| `mpm_gpu/validate.py` | JSON hardware-gate reports and 10k/100k binning benchmarks |
| `tests/test_gpu_contract.py` | CPU-verifiable packaging, topology, coordinates, configuration and error-report contracts |
| `tests/test_gpu_hardware.py` | Actual Metal runtime gate, opt-in via environment variable; unavailable hardware fails when explicitly requested |

No mechanics shaders, device nonlinear solver, GPU runner, GPU diagnostic
adapter or snapshot integration have been implemented. The validation stages
for those components must be implemented as those components are added.

## Numerical references and next implementation work

Use `granular.py:assemble` and `granular.py:step` as the equations/timestep
reference, and `granular_mixed.py:GradientOperator`, `diagonal_coulomb` and
`solve_velocity` as the compact algorithm reference. Read
[the compact solver notes](granular-mixed.md) and
[the rigid-grain model](granular-prototype.md). The source paper is
[Daviet (2026), sections 5.2–5.3](https://research.nvidia.com/labs/prl/mixed_mpm/files/mixed_mpm.pdf).
Our nodal-stress rigid-grain model is an adaptation, not the paper's complete
constitutive model/discretization. P0 remains outside the first GPU release.

Implement the mechanics gate next, then the solver, then runner integration,
following the order and quantitative acceptance criteria in [gpu.md](gpu.md).
The following details are easy to miss:

- Assemble two integrated gradient components per stress/velocity pair using
  device bins and node-owned gathers. Use the same representation for `B` and
  `B.T`; merge periodic aliases before metric construction. Static CPU topology
  setup is allowed; evolving operator assembly and particle binning stay on GPU.
- Preserve full normal/tangent coupling in the global operator. The cached
  local tangential eigensystem and positive regularizer are a preconditioner,
  not constitutive compliance. Keep exact opening, zero-friction, sticking and
  safeguarded scalar sliding branches with explicit FP32 stagnation handling.
- Exclude positive normal gaps from residual scaling. Use the dimensionless
  maximum block Euclidean natural residual and independent FP64 frozen checks.
  GPU default tolerance is `1e-4`; reject explicit values below `1e-5`.
- Do not simplify the solver to plain colored sweeps: a pre-port CPU feasibility
  probe stalled near `1.06e-4` on static fixtures even after 2,000 sweeps without
  Anderson/Newton. This probe was not a hardware or FP32 validation result.
- Anderson uses 12 history differences (13 states), hierarchical TSQR and a
  device-side small Jacobi SVD, with an FP32 rank cutoff. Project candidates to
  feasibility and accept only a reduction in the full residual.
- Newton uses matrix-free Jv, GMRES restart 64 for up to four cycles, two-pass
  orthogonalization, device Givens/Hessenberg control and backtracking fractions
  `1, .5, .25, .125`. The CPU reference's Newton candidates are unprojected;
  feasibility and complementarity errors remain checked by the full residual.
- Keep solver vectors, histories, convergence/acceptance decisions and small
  linear algebra on device. Host orchestration can submit bounded batches and
  read compact status records. Use multiple dispatches for global dependencies;
  cross-workgroup spin barriers and hidden CPU fallback are prohibited.
- Store particle locations as integer cell plus FP32 local coordinates. For
  deformation, use a stable planar exponential with scaling/squaring and retain
  full 3×3 `F` compatibility: left-multiply by `blockdiag(exp(L_2x2*dt), 1)`.
  Preserve `J`, `V`, `rho`, stress, plastic/material work and diagnostic semantics.
- `run()` should read scalar diagnostic/control records each step and snapshot
  large arrays only at output times and final return. Standalone `step()` must
  synchronize host state. Reject live-NumPy custom callbacks explicitly rather
  than leaving them with stale arrays. Preserve built-in CFL and clock behavior.

Use `kinematics.py`, `diagnostics.py`, `stability.py`, `initialise.py`,
`particles.py`, `transfers.py` and `MPM.py` for the rest of those contracts.
Add backend metadata to existing outputs rather than inventing a parallel I/O
format. Keep CPU as the default and importable without native GPU dependencies.

## Evidence from the Intel development host

- Full suite at the runtime checkpoint: **544 passed, 1 hardware test skipped**;
  two existing Astropy convolution warnings.
- CPU-only packaged validation: `benchmarks/metal_runtime_cpu.json`; explicitly
  labelled `cpu_checks_passed_gpu_unvalidated`.
- Source distribution and wheel built successfully. CPU checks also passed from
  the zipped wheel outside the source checkout; source fingerprints matched.
- All nine shader entry points translated to Metal source using Slang
  **2026.12**, the version selected by SlangPy 0.43.1's tagged build metadata.
  Slang emitted binding-layout warnings E31106/E31107; they were not suppressed.
- Apple Metal compilation could not run: Xcode's launcher reported a missing
  Metal Toolchain component. No Metal device execution was attempted/certified.
- Runtime signatures and entry-point binding were checked against the pinned
  wheel's `.pyi` files and SlangPy's tagged native binding source. Do not assume
  examples for newer SlangPy versions apply unchanged.

Do not carry over temporary `/private/tmp` paths or the Intel managed-worktree
path as prerequisites. All source needed for the first gate is committed.

## Reporting and delivery

Commit and push each tested milestone on this branch, preserving unrelated
work. Update [gpu.md](gpu.md), this handoff and validation-stage implementations
as the work advances. Reports must distinguish CPU checks, offline shader
translation, actual device execution and evolving science acceptance.

Before claiming a GPU MPM speedup, run matched CPU/GPU simulations on the same
Apple silicon Mac, with the same initial state, timestep, diagnostics and
requested tolerance. The existing binning timings do not establish a solver or
simulation speedup. Follow the warmup, five alternating paired samples,
connected/disconnected domain cases and accounting requirements in `gpu.md`.
Retain convergence failures and the baseline's known density/energy limitations.
