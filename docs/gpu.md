# Experimental Apple silicon backend

For continuation on the target Mac, start with the [agent handoff](gpu-handoff.md).

`P.backend = "metal"` runs the granular rigid-grain timestep
(`granular.step` with the compact velocity solver) resident on an Apple
silicon GPU, in FP32. All four hardware gates (runtime, mechanics kernels,
nonlinear solver, runner and science) passed on an Apple M4 (macOS 27.0.1,
SlangPy 0.43.1, Metal toolchain 32023.921); the complete report is
`benchmarks/metal_runner_m4.json`. A Metal run never runs CPU mechanics as a
fallback: unsupported settings are rejected, and a missing device or an
unconverged solve raises.

On the M4, matched runs were 3.3 to 26 times faster than the CPU (see
Performance below). A connected domain of 10,000 particles does not converge on
either backend within the sweep budget.

The baseline is commit `3d72a69`, without the concurrent density experiments.
The scientific model and CPU numerical defaults remain unchanged.

## Example

`granular_square_drum` is the square drum of `bi_square_drum` (rotating
gravity, two grain sizes, grid segregation) with rigid grains in place of the
Pouliquen material. Friction is the constant quasi-static value
mu = tan(20.9 deg), so the rate dependence of mu(I) and the elasticity are
absent. The last argument selects the backend:

```sh
uv sync --extra gpu --python 3.14
uv run python MPM.py granular_square_drum 11 2 metal
```

Time stepping is dynamic. The material time sqrt(h/|g|), the time for a grain
at rest to fall one grid spacing, bounds the step while the bed is at rest;
otherwise the advection and segregation limits apply, each scaled by `CFL`.
Once the bed avalanches (near 32 deg of tilt at ny = 11) the solve needs more
than the default 2,000 sweeps, so the example allows 20,000.

The example uses one-sided walls (`B.separable`, described in
[granular-prototype.md](granular-prototype.md)). With two-sided walls a row of
grains stays inside the first cell along each wall as the drum turns: at
t = 1 s, ny = 11, Fr = 2 (450 particles), 32 particles remain in the bottom
row of cells on both backends, against none with one-sided walls. Re-solving
after a released wall node is pushed back raises the CPU cost from 990 solves
and 14,600 sweeps to 2,560 solves and 40,500 sweeps over that second. On the
M4 the CPU took 1,018 steps in 109 s and Metal 1,013 steps in 69 s, against
38 s and 26 s with two-sided walls. The step counts differ because the dynamic
step sequence depends on FP32 rounding.

## Run the first hardware gate

Use an Apple silicon Mac, macOS 26 or newer, native CPython 3.14, Xcode and its
Metal compiler toolchain. Both `xcrun --find metal` and
`xcrun metal --version` must succeed. The Command Line Tools alone are not
enough: shaders translate, but pipeline creation fails with "pass-through
compiler not found". Xcode may also provide a forwarding executable even when
its Metal Toolchain component is missing (`xcodebuild -downloadComponent
MetalToolchain` installs it). The GPU extra requires SlangPy 0.43.1 or newer;
its Apple silicon wheel requires macOS 26 and has no CPython 3.15 build yet.
The extra is platform marked so CPU installation still works elsewhere.

From this branch:

```sh
uv sync --extra gpu --python 3.14
uv run python -m mpm_gpu.validate --stage runtime --benchmark --output metal-validation.json
```

Or install with `python -m pip install '.[gpu]'` and use the same Python command.
Send back `metal-validation.json`, including failures. Run with `--debug` for
API validation layers if a case fails; omit it for timings.

The runtime gate checks shader-resource packaging and compilation, scalar
uint32/FP32 ABI, integer atomics, histogram/scan/scatter, dispatch barriers
and device-guarded dispatch. Scan cases straddle
workgroup boundaries and include hierarchical recursion. Binning checks dense,
sparse, single-bin and invalid-key distributions. Reductions check sum/min/max,
NaN and infinity handling. Invalid data sets a device status bit instead of
being silently treated as a valid result.

SlangPy 0.43.1 on Metal implements neither indirect compute dispatch nor
timestamp queries. Indirect dispatch logs a warning and executes nothing;
timestamp reads return `SLANG_E_NOT_AVAILABLE` and the device reports a
timestamp frequency of zero. Both are recorded under `capabilities` in the
report and are not gate requirements. Work whose extent is decided on the
device is launched with a host upper bound, and each thread reads the active
extent from a device record (`guarded_marker`). Timings are host wall time
around a waited submission; command buffers are encoded outside the timed
region because Metal command buffers cannot be resubmitted.

The optional 10,000/100,000-element benchmarks warm up once and alternate CPU/GPU
order across five pairs. On the M4 the median resident binning time was 0.78 ms
against 5.1 ms for NumPy at 100,000 keys, and 0.57 ms against 0.34 ms at
10,000 keys, where submission and wait latency dominates. They time **binning primitives only**, with resident
GPU inputs. Both produce the order of a stable sort by key (the GPU sorts
each bin after its atomic scatter). These measurements are not MPM speedups. Setup/upload, compilation,
encoding and readback are reported separately; allocation bytes are cumulative,
not a measured device-memory peak. Hardware/compiler errors retain their tracebacks in the JSON report.

On an Intel development Mac:

```sh
python -m mpm_gpu.validate --cpu-only --output metal-validation-cpu.json
```

That checks imports, resources, backend defaults, periodic topology/coloring,
transpose incidence and coordinate reconstruction. It returns
`cpu_checks_passed_gpu_unvalidated`, never a Metal pass. The regular pytest suite
has the same separation. For the real device test:

```sh
POLY_MPM_TEST_METAL=1 python -m pytest tests/test_gpu_hardware.py -q
```

An explicitly requested but unavailable GPU fails rather than skipping.

`--stage all` runs the runtime, mechanics, solver and runner gates. A single
`--stage` also runs every preceding stage, and a stage runs only after the
previous one has passed. Exit 3 is reserved for requested stages that are not
implemented; none remain.
Exit 0 means the requested implemented gate passed (or explicitly requested CPU
checks passed); exit 1 means failure/unavailable; argument errors use exit 2.
A runtime result cannot certify mechanics or solver convergence.

## Backend contract

CPU remains the default, with contact tolerance `1e-9`. On Metal the default
tolerances are `1e-4` for both the Coulomb and the CIDR solvers.
Explicit settings are preserved; GPU tolerances below `1e-5` are rejected.
No supplied tolerance is silently loosened.

The runner supports one intrinsic 2D granular material, nodal stress, linear
interpolation, point quadrature, FLIP (including PIC
blend), PIC, `apic_gradient`, fixed frictionless/no-slip walls (two-sided or
one-sided), frictional walls (`B.wall_mu`), interior and silo walls, periodic
x, inlets and outlets, built-in static/transport-CFL time stepping and grid
segregation (up to eight grain sizes, `smooth_gamma_dot=False`, optional
`normalise_phi`). P0 and coarse stresses, plane strain, viscosity, density
correction, increment stopping and callbacks that require live NumPy arrays
are rejected. These restrictions do not change CPU support.

Inlets and outlets change the particle count, which the resident buffers
cannot. The device counts the particles that `process_boundaries` would remove
(in the per-step particle reduction), and the host knows from the clock when
an inlet fires. On such a step the runner synchronises the state to the host,
applies `process_boundaries` as the CPU step does, books the boundary flux in
the diagnostics and rebuilds itself for the new count. A silo that drains
every step therefore runs mostly on the host path.

The runtime uses only explicitly selected Metal devices. SlangPy is imported
lazily; normal CPU execution never imports the native GPU module. Scalar
structured buffers have four-byte stride. Package resources are retained for
the runtime lifetime, including zipped installations. Cross-workgroup phases
use distinct dispatches and API barriers. There are no floating-point atomics
or cross-workgroup spin barriers in the primitives.

## Remaining implementation gates

Each gate needs a real Apple silicon result before dependent integration is
called verified. Implementations must retain failures, keep the CPU reference
independent and avoid a CPU solver fallback.

1. **Runtime/primitives (passed on Apple M4).** Scans, integer atomics,
   binning, reductions, barriers and guarded dispatch are validated. Indirect
   dispatch and GPU timestamps are unavailable in the current SlangPy Metal
   backend and must not be used.
2. **Mechanics kernels (passed on Apple M4).** Report:
   `benchmarks/metal_mechanics_m4.json`
   (`python -m mpm_gpu.validate --stage mechanics`). All comparisons are with
   independent FP64 CPU references on wall and periodic columns, including
   one- and two-cell periodic strips and particles on the closed upper edge.
   - Compact operator (`mpm_gpu/mechanics.py`): device cell binning and
     row-owned assembly of nine integrated gradient 2-vectors per stress row,
     with periodic aliases matched by canonical index; `B u` gathers by stress
     row and `B.T z` by velocity node through the static transpose incidence.
     Relative block, `B u` and `B.T z` errors against
     `assemble(..., compact=True)` are at most 2.4e-7; adjoint pairing errors
     at most 5e-8.
   - Local metrics, conservative norm bound (reduced on the device),
     regulariser and tangential eigensystem agree with the `solve_velocity`
     set-up to at most 3.2e-7. The rotation uses the symmetric Jacobi formula,
     since Metal's `atan2(0, 0)` gave a degenerate rotation for isotropic blocks.
   - The exact local Coulomb solve (opening, frictionless, sticking and
     safeguarded sliding, with FP32 bracket-collapse/stagnation exits and an
     unconverged status bit) agrees with `diagonal_coulomb` to 2.1e-7 relative
     to the unconstrained response on 4,007 rows spanning sixteen decades of
     load and six of metric.
   - Transfers: node-owned P2G (mass, raw and masked momentum, `dt m g`, trial
     velocity, free inverse mass; the periodic merge is part of the gather) and
     particle-owned G2P (FLIP with PIC blend, PIC, `apic_gradient`, `L`, and
     advection in split coordinates) agree with the CPU `granular.step` grid
     sequence and `move_material_points` to at most 2.3e-7, and positions to
     6e-8 grid spacings.
   - Deformation: `F <- blockdiag(exp(L dt), 1) F` on the full 3×3 `F`, with a
     closed-form 2×2 exponential after scaling and squaring, agrees with SciPy
     `expm` to at most 4.9e-6 (rotations up to 20 rad and stretches near
     `e^10` in one step). `J` is updated as `J exp(tr L dt)`, which equals
     `det F`: an FP32 3×3 determinant of a strongly stretched `F` was in error by
     6%. `J`, `V` and `rho` agree to 1e-6. Nonfinite or nonpositive results set a
     status bit and leave the particle unchanged.

   Stress, work and diagnostic updates are part of runner integration (gate 4).
3. **Resident nonlinear solver (passed on Apple M4).** Report:
   `benchmarks/metal_solver_m4.json` (`python -m mpm_gpu.validate --stage solver`,
   which runs the runtime and mechanics stages first). The solver
   (`mpm_gpu/solver.py`, `shaders/solver.slang`) follows `solve_velocity`:
   colored Gauss-Seidel sweeps with the exact local solve, an exact velocity
   refresh after every sweep, the dimensionless maximum block Euclidean natural
   residual with positive normal gaps excluded from the scaling, Anderson
   acceleration with 12 differences and matrix-free Newton every 25 sweeps when
   the residual lies between the tolerance and `1e-3`.
   - Control: the host encodes batches of 25 sweeps and reads a 64-byte state
     record after each one. Each kernel reads a device flag word and returns
     when its gate is closed; single-thread kernels make every scalar decision
     (convergence, budget, acceptance, history length, Krylov termination). A
     Newton batch is encoded only when the device has set its request flag.
     Every workspace is allocated when the solver is constructed; a repeated
     solve allocates no device memory.
   - Sweeps: rows of one static Q1 color share no velocity node, so each row
     updates its own nodes' velocities without atomics. Colors come from the
     static topology (9 on walls, up to 12 with periodic x), not from the
     CPU's mass-dependent greedy coloring, so the iteration order differs from
     the CPU.
   - Anderson: the 13 × 13 R factor of `[df | f]` comes from Householder TSQR
     over 256-row tiles, reduced level by level (four levels at 108,300
     unknowns). A one-sided Jacobi SVD of `R_df` gives the minimum-norm fit,
     with singular values below `1e-5` of the largest discarded (the FP64
     reference uses `1e-10`). Candidates are projected to feasibility and
     accepted only on a decrease of the full residual.
   - Newton: GMRES(64) for up to four cycles from `x = 0` with `rtol = 1e-5`,
     two classical Gram-Schmidt passes, Givens rotations and back substitution
     on the device, then backtracking `1, .5, .25, .125` on the unprojected
     candidate. An accepted step clears the Anderson history.
   - Results: on ten wall and periodic columns, static and with random particle
     motion (up to 651 stress rows), every solve reached `1e-4`. The
     independent FP64 frozen natural residual of the device contact stress was
     at most 9.7e-5; velocities agree with the CPU solve at `1e-9` to 1.2e-4 of
     the largest speed. A warm start from the converged stress needs no sweeps.
     The device natural residual and generalized Jacobian-vector product agree
     with the FP64 Newton block to 1.7e-6 and 1.7e-7, with identical frozen
     branches. Anderson coefficients agree with FP64 `lstsq` to at most
     2.9e-5, including rank-deficient and short histories.
   - Newton paths: one accepted step (FP64 residual 7.6e-5 to 4.4e-6) and one
     rejected step whose candidates all raise the FP64 residual, in agreement
     with the device decision. Device corrections have FP64 linear residuals of
     1.4e-5 and 6.7e-5 at the device's frozen branch. GMRES also converges in
     one step for `J = I`, and for `J = 0` records a breakdown, makes no
     correction and rejects the unchanged candidate.
   - Explicit failures: an exhausted sweep budget (status 256), a nonfinite
     trial velocity (status 2) and compression with no free velocity
     (status 128) are reported as unconverged, and `SolveResult.require()`
     raises.

   Binning is deterministic: after the atomic scatter each cell's particle
   indices are sorted (`sort_bins`), so FP32 sums and the whole iteration path
   are identical between runs; the runner gate checks that two runs are
   bit-identical. Before this sort the same moving wall fixture took between
   112 and 155 sweeps over five runs. The packing gap was computed on the CPU
   for these checks; it is computed on the device in gate 4.
4. **Runner and science validation (passed on Apple M4).** Report:
   `benchmarks/metal_runner_m4.json` (`python -m mpm_gpu.validate --stage runner`
   or `--stage all`). `mpm_gpu/runner.py` and `shaders/runner.slang` implement
   the step; `MPM.step` and `MPM.run` select it for `P.backend = "metal"`.
   - Step order follows `granular.step`: start-of-step binning, P2G with
     gravity, wall masks and trial velocity; occupancy and packing gap; the
     preconditioner; the solve, warm-started from the previous contact on rows
     that hold grain volume; final velocity, internal force and the constraint
     and external impulse/work budgets on unique nodes; G2P with advection in
     split coordinates; deformation; stress gathered from the nodal contact
     and scaled by `rho / rho_s`; strain and material/plastic work (Kahan
     summation per particle); the inventory; end-of-step binning, which gives
     the packing diagnostics and is reused as the next step's start-of-step
     binning. Every workspace is allocated once per runner.
   - Per step the host reads the solver's state record between batches and one
     reduced record (30 numbers) at the end. These drive `P.timestep_limits`,
     the dynamic timestep, `Diagnostics.finish` (through new `now`/`support`
     arguments, so the history rows have the CPU format) and
     `G.granular_stats`. Budgets are reduced in FP32 per step and accumulated
     in FP64 on the host. Particle mass is constant on the device, so the mass
     error is zero.
   - `run()` uploads the state once and calls `snapshot()` before each save-time
     callback and on return; standalone `step()` uploads, steps and snapshots.
     A snapshot writes all particle fields and the mechanical grid fields. The
     deposited grid diagnostics (`deposit_diagnostics`, grain fractions,
     `shear_rate`) are formed at snapshot time with the CPU routines from the
     snapshot and the last step's start-of-step stencil, which is how the CPU
     step forms them. They do not feed back into the mechanics.
   - Rejected for Metal: custom `after_every_timestep` callbacks (they would see
     stale arrays), `O.measure_energy`, custom timestep controllers, damping,
     `O.check_positions` and the anisotropy mode. Save-time, start and end
     callbacks are supported and receive snapshots.
   - FP32 thresholds: packing clearances within `64 eps_32` of the capacity are
     set to zero (the CPU uses `64 eps_64`, for the same purpose of not turning
     summation roundoff into an opening gap), and stencil coverage counts as
     incomplete below `1 - 1e-6` (CPU `1 - 1e-10`). The CPU prototype's
     `granular_max_nodes` limit is not applied on Metal.
   - Results against the CPU at the same tolerance `1e-4` and state:
     | Case | Result | Gate |
     | --- | --- | --- |
     | FP64 frozen natural residual of one Metal step (static, collapse, periodic incline) | at most 1.08e-4 | ≤ 2e-4 |
     | Relative mass error, static column | 0 | ≤ 1e-5 |
     | Static column displacement, 0.2 s at `dt = 0.02` | 9.1e-6 grid spacings | ≤ 1e-3 |
     | Periodic incline, slope 0.75: speed against `g (sin θ - μ cos θ) t` | 0.29% (CPU 0.28%) | ≤ 1% |
     | Periodic incline, slope 0.25 | static | — |
     | Collapse positions/runout, FLIP, FLIP+PIC 0.3, PIC, `apic_gradient`, dynamic CFL | at most 6.4e-6 grid spacings | ≤ 1% |
     | Integrated nodal overpacking, Metal minus CPU | at most 3.0e-7 | ≤ 1e-3 |
     | Free fall at low packing | stress free; velocity change within 1e-6 of `g dt` | — |
     | Upward gravity, one-sided or frictional walls: velocity change against `g t` | 1.3e-8 m/s; positions 2.2e-8 grid spacings from the CPU | ≤ 1e-5 relative |
     | Collapse positions/runout, one-sided walls (re-solves when a released wall node is pushed back) | 1.3e-5 grid spacings | ≤ 1% |
     | Collapse positions/runout, frictional walls `wall_mu = 0.3` | 7.5e-5 grid spacings | ≤ 1% |
     | Periodic layer on a frictional floor, `tan θ = 0.4`, `μ = 0.5`, `wall_mu = 0.2`: speed against `g (sin θ - wall_mu cos θ) t` | 0.34% (CPU 0.32%) | ≤ 1% |
     | Same layer, `wall_mu = 0.6` | at rest (2.8e-4 m/s; CPU 1.9e-4) | — |
     | Silo draining through `silo_bottom` into `outlet_bottom`, with and without `inlet_top` (0.3 s) | same particle counts (155, 170 of 240); positions 4.8e-5 grid spacings; mass budget exact | ≤ 1% |
   - Budgets: momentum residuals relative to the external impulse are at most
     1.9e-7 on Metal (CPU 1e-15, FP32 sums). Energy residuals relative to the
     larger of external and constraint work agree with the CPU to 5e-6 in the
     collapses (for example -0.04197 and -0.04198 for FLIP). In the static
     column and static incline the work terms are close to zero, so these
     ratios are dominated by solver-tolerance noise and differ between
     backends (static incline 0.14 CPU, -0.06 Metal).
   - Output: save times, files (figures and `data/diagnostics.json`) and
     snapshot contents at each save match a CPU run; a custom
     `after_every_timestep` is rejected.
   - Reproducibility: two Metal runs of a 20-step collapse are bit-identical
     in positions, velocities, stresses, contact and sweep counts.
   - Scale: a collapse with 3,200 particles and 1,271 stress rows (3 steps)
     agrees with the CPU to 1.1e-6 grid spacings. Neither solver converges to
     `1e-4` within the default 2,000 sweeps on the first step of this state
     (CPU residual 1.4e-4, so this is a property of `solve_velocity`, not of the
     port); with 20,000 both need about 2,000-2,500 sweeps. Timings are in the
     performance section below.

## Performance

Report: `benchmarks/metal_benchmark_m4.json` (`python -m mpm_gpu.benchmark`),
Apple M4. Both backends start from one prepared state with the same `dt`
(0.002 s), tolerance (`1e-4`), sweep budget and diagnostics; one Metal warm-up
run compiles the kernels before timing; paired samples alternate the CPU/Metal
order. Times are host wall time for the complete run, including every waited
submission and the final snapshot. The larger cases have fewer than five pairs
because of the CPU run time; the report records the pairs run. Disconnected
cases repeat one collapsing column (0.5 m wide, 1 m high, grid spacing
0.05 m, 4 particles per cell) side by side.

| Case | Particles | Stress rows | Steps | CPU | Metal | CPU / Metal (median) |
| --- | --- | --- | --- | --- | --- | --- |
| Static column | 6,400 | 651 | 10 | 2.6 s | 0.07-0.11 s | 25 |
| Collapse | 3,200 | 651 | 3 | 14.4 s | 4.4-5.0 s | 3.3 |
| Disconnected collapse | 9,600 | 1,891 | 3 | 28.4 s | 5.9 s | 4.8 |
| Disconnected collapse | 99,200 | 19,251 | 2 | 360 s | 14.0 s | 26 (one pair) |
| Connected collapse | 10,000 | 7,676 | 2 | not converged | not converged | — |

Final positions agree with the CPU to at most 2.2e-6 grid spacings. The
collapse cases are dominated by the first step from rest, which needs about
2,000 sweeps and 60 Newton polishes on either backend. The connected
10,000-particle collapse does not converge on either backend within 20,000
sweeps on its first step (residuals 1.1e-2 CPU, 1.4e-2 Metal): the solver's
iteration count grows with the size of a connected domain, which limits both
backends until a better-scaling solver is available. Device memory for one
runner was 62 MB at 99,200 particles.

Saving output costs the same on both backends apart from the Metal snapshot.
For the disconnected 9,600-particle collapse, one save with the material-point
and continuum figures took 0.34 s on either backend; the Metal snapshot that
precedes it took 8.9 ms (`plotting` in the report,
`python -m mpm_gpu.benchmark --plotting`). Against the warm Metal step of
about 2 s for that case, a save every step would add roughly a sixth; at the
usual save intervals the figures are a small part of the run.

Measured costs on the M4 that shaped the implementation:

- Each dispatch in its own compute pass cost about 27 us of GPU time; within
  one pass, 2-6 us. The runtime now keeps one compute pass per command
  encoder. Dispatches within a pass run in order and see earlier writes; the
  runtime gate and every later gate pass on this path.
- Encoding costs 4-7 us per dispatch in Python, mostly the dispatch call
  itself. An Arnoldi step uses 9 dispatches: the Jacobian product
  (normalising the new basis vector on the way), two classical Gram-Schmidt
  passes and the norm, each inner product with one 64-thread workgroup per
  vector and chunk, which streams the basis at memory bandwidth.
- Up to 1,024 stress rows, a solver iteration runs as two 1,024-thread
  single-workgroup kernels around the TSQR/SVD dispatches (265 us against
  381 us per iteration at 651 rows). Above that the multi-dispatch iteration
  is faster (778 us against 1,492 us at 2,501 rows). Both paths use the same
  control functions and are covered by the solver gate.
- Tried and removed because they were slower: a whole colored sweep in one
  workgroup, a whole GMRES cycle in one workgroup (12.5 ms against 8.8 ms at
  651 rows, 80 ms against 17 ms at 9,801 rows) and Gram-Schmidt kernels fused
  with the following inner products.
- The host sends sweeps 5 at a time and Arnoldi steps 16 at a time, reading
  the 64-byte state record between batches, and stops encoding once the
  device has closed the corresponding gate.

Final science criteria:

| Quantity | Gate |
| --- | --- |
| Independent FP64 frozen natural residual | ≤ 2e-4 |
| Relative mass error | ≤ 1e-5 |
| Static displacement | ≤ 1e-3 grid spacing |
| Incline speed error | ≤ 1% |
| Collapse position/runout difference from CPU | ≤ 1% grid spacing |
| Additional integrated nodal overpacking vs CPU | ≤ 1e-3 |

Compare dimensionless momentum/work budgets to the baseline. The baseline's
known energy-budget residual must not be forced to zero or hidden. The GPU port
must not silently incorporate a density-model change.

The eventual full benchmark will use science fixtures plus connected and
repeated disconnected 10k/100k-particle domains, with identical CPU/GPU states,
time steps, diagnostics and tolerance on the same target Mac. Keep cold startup
and compilation separate from warm timings; record wall time per waited
submission (kernel timestamps are unavailable in the current backend), transfer
bytes/time, memory, iteration counts and plotting cost. The current results
are under Performance above, including the plotting cost; repeated connected
domains at 100k are not measured, because connected refined domains do not
converge on either backend.

## Segregation

`P.segregate = True` runs `integrators.segregation_update` on the device
within each step, in the CPU order: after the stress update and before
advection, with the start-of-step stencil (`mpm_gpu/segregation.py`,
`shaders/segregation.slang`).

- Grid side: mass-weighted grain fractions and the pressure and inertial-number
  deposits are gathered per node. The shear rate and `grad pk` use the masked
  central differences of `Grid.calculate_gradient`, and the diffusivity
  `10 s^3 |shear rate| / sqrt(I)`.
- Faces: the static face lists, wall closures and periodic seams are built
  once on the host with the CPU scheme's own orientation and padding helpers;
  the central-upwind flux, limited reconstruction and harmonic face density
  are computed per face on the device.
- Particles: the inventory-aware limiter and particle exchange of
  `_particle_exchanges` run as node- and particle-owned gathers, with no
  floating-point atomics. Each particle keeps its mass; only its fractions
  change.
- The segregation advection and diffusion timestep limits and the species
  inventory come from the per-step record; snapshots fill `G.phi`, `G.dphi`,
  `G.pk` and `G.grad_pk` as the CPU step leaves them.

Gate cases (in `benchmarks/metal_runner_m4.json`), each against the CPU from
the same state at tolerance `1e-4`:

| Case | Steps | CPU fraction change | Metal minus CPU | Species error |
| --- | ---: | ---: | ---: | ---: |
| Periodic incline, two sizes, `c = 1e-2` | 40 | 0.037 | 1.0e-5 | 4.2e-8 |
| Walled collapse, two sizes, `c = 0.1` | 30 | 0.23 | 7.0e-6 | 3.5e-8 |
| Periodic incline, three sizes, `normalise_phi` | 40 | 0.041 | 1.1e-5 | 4.0e-8 |

Grid fields (`phi`, `dphi`, `pk`, shear rate) agree to at most 4.5e-4
relative and the segregation timestep limits to 5e-4. Species errors are
relative to the species mass and come from FP32 fractions; the per-step
diagnostics inventory is an FP32 sum, so its species error can read zero when
the change is below FP32 resolution.

No separate package or GPU-specific input/output format is needed.

## Critical-state CIDR

`granular_solver = "cidr"` runs on Metal (`mpm_gpu/cidr.py`,
`shaders/cidr.slang`) with the algorithm of the CPU solver (`cidr.py`): the
method of multipliers with semismooth Newton subproblems on the Moreau envelope
of the dissipation potential, whose gradient comes from an exact proximal step
at each stress row. The device version differs as follows:

- Each Newton system is solved matrix-free by conjugate gradients to a relative
  residual of 1e-3, preconditioned by the 2x2 block of the Newton matrix at each
  node. Dot products are summed in shared memory within the kernels that form
  them, so an iteration is six dispatches. On grids of up to 2,048 nodes the
  iterations run in one 1,024-thread workgroup, many per dispatch, with group
  barriers in place of dispatch boundaries (three times faster at 651 nodes,
  twice as slow at 7,676). The host reads a flag after the number of iterations
  the previous solve took, then every 16.
- Line searches use the sign of the directional derivative, which FP32
  resolves, rather than objective values, whose decrease near convergence it
  does not. The objective is convex along the direction, so a step with a
  nonpositive derivative decreases it. A full step is accepted if its derivative
  is below a tenth of the initial slope's magnitude, otherwise regula falsi
  (Illinois) brackets the minimiser. When the accepted step is the last one
  evaluated, its envelope is reused by the next Newton iteration. The local
  Newton in the proximal step accepts a step that reduces the objective or its
  gradient, and stops on the gradient.
- Newton stops when its step is below 1e-6 of the trial velocity in the mass
  norm. A test on the decrement, as on the CPU, depends on the stiffness 1/gamma
  of jammed rows: at the larger smoothing below it let a column at rest drift at
  9e-6 m/s; with the step test it moves by less than 6e-7 m/s.
- `granular_cidr_gamma` defaults to 1.0 (in units of `dt / rho`) rather than the
  CPU's 0.075. With iterative linear solves the smoother envelope needs fewer
  Newton iterations and CG iterations per Newton system: on the connected
  collapse, 8 steps took 1.8 s against 5.9 s at 0.075, with positions within
  4.4e-7 m.
- `alpha(I)` is a cubic Hermite table in `ln I` with slopes from its ODE.

The envelope kernel agrees with the CPU envelope on 400 random rows: the same
jammed rows, gradients to 5e-5 of the stress scale, Hessians to 1e-5 (median);
rows within roundoff of the cap or of zero pressure can take either side of the
Hessian's jump (`tests/test_gpu_cidr.py`).

On both backends a warm step starts from the transferred velocity before the
step's body force: the trial velocity includes one step of gravity, which is
incompatible with jammed material at walls and made the first warm step fail on
the CPU.

Benchmark (`benchmarks/metal_benchmark_cidr_m4.json`; `python -m
mpm_gpu.benchmark --solver cidr --steps 20 --pairs 1`), Apple M4, 20 steps of
the benchmark columns, one pair, host wall time, each backend with its defaults
(tolerance 1e-5 on the CPU, 1e-4 on Metal):

| Case | Particles | CPU s/step | Metal s/step | CPU / Metal | Newton per warm step (Metal, mean) | CG per step (Metal) |
| --- | --- | --- | --- | --- | --- | --- |
| static column | 6,400 | 0.15 | 0.018 | 8.5 | 0.8 | 43 |
| collapse | 3,200 | 0.65 | 0.049 | 13 | 13 | 731 |
| three collapses | 9,600 | 3.5 | 0.079 | 44 | 14 | 783 |
| connected collapse | 10,000 | 20.7 | 0.54 | 38 | 22 | 3,212 |

Before the changes above, Metal took 0.044, 0.186, 0.362 and 2.2 s per step on
these cases. CPU and Metal positions agree to 6e-5 grid spacings after 20 steps,
the difference of their tolerances. The Coulomb solver in
[Performance](#performance) is a different material law and was timed over 3
steps from rest, so the comparison is indicative only: there, the collapse of
3,200 particles took 1.7 s per step on Metal and the connected collapse of
10,000 did not converge on either backend.

The remaining cost is in slow, nearly jammed flow on fine grids, where Newton
needs tens of iterations and each needs hundreds to thousands of CG iterations.
A multigrid preconditioner would address the second; the first is the
globalisation question of
[granular-cidr-design](granular-cidr-design.md#frozen-step-results).

## Sources

- [SlangPy release and wheel metadata](https://pypi.org/project/slangpy/).
- [Low-level compute kernel API](https://slangpy.shader-slang.org/en/stable/src/tutorials/compute_shader.html).
- [SlangPy API reference](https://slangpy.shader-slang.org/en/stable/src/api_reference.html).

Runtime signatures were checked against the `.pyi` API shipped inside the
0.43.1 Apple silicon wheel. A newer SlangPy must re-pass the runtime gate.
