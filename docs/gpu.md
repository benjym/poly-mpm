# Experimental Apple silicon backend

For continuation on the target Mac, start with the [agent handoff](gpu-handoff.md).

The first implementation milestone supplies the optional Metal runtime, resident
integer binning, hierarchical scans/reductions, static Q1 incidence, split
particle coordinates and a hardware validation command. **The GPU MPM timestep
is not implemented yet.** Selecting `P.backend = "metal"` raises an explicit
`ValidationPending` error. It never runs CPU mechanics as a fallback.

The baseline is commit `3d72a69`, without the concurrent density experiments.
The scientific model and CPU numerical defaults remain unchanged.

## Run the first hardware gate

Use an Apple silicon Mac, macOS 26 or newer, native CPython 3.14, Xcode and its
Metal compiler toolchain. Both `xcrun --find metal` and
`xcrun metal --version` must succeed. Xcode may provide a forwarding executable
even when its Metal Toolchain component is missing. The GPU extra pins
SlangPy 0.43.1; its Apple silicon wheel requires macOS 26. The extra is platform
marked so CPU installation still works on Intel and other platforms.

From this branch:

```sh
uv sync --extra gpu --python 3.14
uv run python -m mpm_gpu.validate --stage runtime --benchmark --output metal-validation.json
```

Or install with `python -m pip install '.[gpu]'` and use the same Python command.
Send back `metal-validation.json`, including failures. Run with `--debug` for
API validation layers if a case fails; omit it for timings.

The runtime gate checks shader-resource packaging and compilation, scalar
uint32/FP32 ABI, integer atomics, histogram/scan/scatter, dispatch barriers,
device-controlled indirect dispatch and GPU timestamps. Scan cases straddle
workgroup boundaries and include hierarchical recursion. Binning checks dense,
sparse, single-bin and invalid-key distributions. Reductions check sum/min/max,
NaN and infinity handling. Invalid data sets a device status bit instead of
being silently treated as a valid result.

The optional 10,000/100,000-element benchmarks warm up once and alternate CPU/GPU
order across five pairs. They time **binning primitives only**, with resident
GPU inputs. NumPy uses stable sorting; GPU scatter has unspecified within-bin
order. These measurements are not MPM speedups. Setup/upload, compilation,
readback and waits are reported separately; allocation bytes are cumulative,
not a measured device-memory peak. GPU timer intervals include the recorded
barriers. Hardware/compiler errors retain their tracebacks in the JSON report.

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

`--stage all` is accepted now. It runs the runtime gate and records subsequent
stages as `not_implemented`; its exit code is 3 after a successful runtime gate.
Exit 0 means the requested implemented gate passed (or explicitly requested CPU
checks passed); exit 1 means failure/unavailable; argument errors use exit 2.
A runtime result cannot certify mechanics or solver convergence.

## Backend contract

CPU remains the default, with contact tolerance `1e-9`. The proposed Metal
mechanics defaults are the compact velocity solver and tolerance `1e-4`.
Explicit settings are preserved; GPU tolerances below `1e-5` are rejected.
No supplied tolerance is silently loosened.

The first runner will support one intrinsic 2D granular material, nodal stress,
linear interpolation, point quadrature, updated mechanics, FLIP (including PIC
blend), PIC, `apic_gradient`, fixed frictionless/no-slip walls, periodic x and
built-in static/transport-CFL time stepping. Segregation, P0, plane strain,
viscosity, particle injection/removal, density correction and callbacks that
require live NumPy arrays are outside that release. These restrictions do not
change CPU support.

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

1. **Runtime/primitives (implemented, hardware validation pending).** Run the
   command above before relying on scan, atomics, synchronization or indirect
   dispatch in mechanics.
2. **Mechanics kernels.** Node-owned gathers using device bins; compact paired
   `B`/`B.T`; periodic aliases merged; local metric and analytic tangential
   eigensystem; exact opening/sticking/sliding local solves; planar exponential
   with scaling/squaring and the full compatible 3×3 deformation gradient.
   Compare transfers, adjoint pairing and deformation against independent FP64
   fixtures before integration.
3. **Resident nonlinear solver.** Colored contact sweeps and exact velocity
   refresh; natural residual with positive gaps excluded from scaling; Anderson
   memory 12 with hierarchical TSQR and device-side Jacobi SVD; matrix-free
   Newton/GMRES restart 64, four cycles, two-pass orthogonalization and device
   Givens/Hessenberg control; safeguarded full-residual acceptance. All histories,
   small linear algebra and convergence/acceptance state stay on the device.
   Explicit failure on exhausted budgets, nonfinite values and unsupported
   constraints. Cover rank-deficient histories, GMRES breakdown, candidate
   rejection and nonconvergence, not just successful solves.
4. **Runner and science validation.** Preserve full timestep ordering, stress,
   work and diagnostics; reduce scalar records on device; snapshot NumPy arrays
   at output times and final return, and on standalone `step()`. Preserve plots,
   Rich progress, output formats and save frequency. Cover static/freefall,
   collapse, incline, periodic support, transfers, CFL and output scheduling.
   Reject custom live-array callbacks explicitly.

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
and compilation separate from warm timings; record kernel timestamps, wall
time, transfer bytes/time, memory, iteration counts and plotting cost. There is
no measured GPU MPM speedup yet.

Segregation and packing corrections can extend this package after the mechanics
gates pass; no separate package or GPU-specific input/output format is needed.

## Sources

- [SlangPy 0.43.1 release and wheel metadata](https://pypi.org/project/slangpy/0.43.1/).
- [Low-level compute kernel API](https://slangpy.shader-slang.org/en/stable/src/tutorials/compute_shader.html).
- [SlangPy API reference](https://slangpy.shader-slang.org/en/stable/src/api_reference.html).

Runtime signatures were checked against the `.pyi` API shipped inside the pinned
0.43.1 Apple silicon wheel. Documentation for newer versions is not an execution
test of this pin.
