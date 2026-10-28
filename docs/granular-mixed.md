# Compact velocity-based granular solve

Set `P.granular_solver = "velocity"` before `prepare(P)` to select the opt-in
CPU implementation. `P.solver = "granular"` still selects the rigid-grain
equations; `granular_solver` selects the numerical algorithm. The default
remains `"hybrid"`, which is the small-problem reference solver.

This pass adapts the stress iteration and nonassociated local flow-rule solver
from [Daviet (2026), sections 5.2–5.3](https://research.nvidia.com/labs/prl/mixed_mpm/files/mixed_mpm.pdf).
It retains our current bilinear velocity and **nodal** stress spaces, particle
quadrature, packing gap, transfer options, wall constraints and viscosity split.
It does not yet implement the paper's discontinuous stress spaces or GPU backend.

## Operator and iteration

The equations remain those in the [rigid-grain guide](granular-prototype.md):

```text
u = u_trial + dt M^-1 B^T r
w = B u + gap
r_n >= 0, w_n >= 0, r_n*w_n = 0
|r_t| <= mu*r_n
sliding: r_t = -mu*r_n*w_t/|w_t|
```

The new path stores the two integrated gradient components for each supported
stress/velocity-node pair. Symmetric tensor blocks of `B` are reconstructed
from those vectors and the stress basis. Both `B u` and `B^T r` use this same
representation, including sums across periodic seams. No full tensor-expanded
`B`, full `W = dt B M^-1 B^T`, or dense global Jacobian is retained during the
solve. Tensor blocks are temporarily expanded once to construct local metrics
and a residual normalization bound.

For each stress node, initialization computes its local diagonal block of `W`,
diagonalizes its tangential part, and caches the rotation and diagonal metric.
The metric ignores normal/tangential cross terms and adds a small positive
proximal regularizer (`1e-6` times the normalized local trace, with a floor).
These approximations affect only the iteration: every update and acceptance
check uses the full `B` and mass operator, without adding compliance to the
governing equations.

Stress nodes are greedily colored by shared free velocity support. This handles
periodic seams and constrained velocity components without assuming a fixed
number of colors. Within a color, NumPy batches independent local solves and
velocity updates. Our continuous nodal stress stencil needs more colors than
the paper's Q1–P0 discretization.

In the rotated local metric, pressure is explicit. Sticking is also explicit;
sliding solves the scalar equation

```text
||q_t / (D_t + alpha)|| = mu*r_n, alpha >= 0.
```

The minimum and maximum tangential eigenvalues bound the root. Safeguarded
Newton iteration solves that monotone scalar problem. The implementation
supports intrinsic 2D and the existing plane-strain reduction, including zero
friction and stress-free opening.

## Convergence and safeguards

Plain colored sweeps were insufficient at the existing `1e-9` tolerance on
some Q1–Q1 fixtures. This implementation adds two CPU safeguards beyond the
paper's basic iteration:

- Limited-memory Anderson acceleration, using at most twelve history
  differences. Extrapolated stresses are made admissible and accepted only
  when the full residual decreases.
- Every 25 sweeps, a nearly converged solve can attempt a semismooth Newton
  correction. Its Jacobian-vector product uses only `B` and `B^T`; restarted
  GMRES stores at most 64 Krylov vectors and runs at most four restart cycles
  per attempt. Backtracking accepts a correction only when the full residual
  decreases. This helps degenerate stick/slip transitions in inclined flow.

Local solves and Anderson candidates enforce stress admissibility. Newton
corrections can leave small admissibility errors, which the final natural-map
residual checks alongside normal complementarity and the sliding law.
The original dense reference solver is never invoked as a hidden fallback.
Nonconvergence raises an error.

The residual scale excludes positive normal packing gaps. A large free-surface
clearance therefore cannot hide compression elsewhere. Without forming `W`,
the solver computes the upper bound

```text
L = max(dt |B| M^-1 |B^T| 1) >= ||W||_infinity
```

and tests the maximum **Euclidean** block norm of the dimensionless natural
map. For a fixed reaction, the scaled projection residual norm is monotone in
`L`; this bound and block norm make the check conservative relative to the
reference's infinity norm and exact `||W||_infinity` normalization. Frozen
comparisons additionally evaluate the original reference residual directly.

`granular_max_iterations` allows up to ten times that many colored sweeps,
consistent with the existing block-solver budget. Diagnostics report sweeps,
accepted acceleration and Newton steps, colors and residual. For this path,
`operator_nnz` counts gradient **blocks** and `operator_storage_bytes` counts
their padded storage and tensor basis. These are not scalar nonzeros or storage
of the old `W`. Preconditioner storage is reported separately; neither measure
is total process memory or peak allocation.

## Verification and timing

Run the matched benchmark with:

```sh
uv run --locked python -m benchmarks.granular_mixed --output /tmp/granular-mixed.json
```

It warms each solver, alternates three paired timed runs and reports medians.
Setup, plotting and progress are excluded; diagnostics, grids, particles,
timesteps, transfers and requested tolerances match. The committed
[report](../benchmarks/granular_mixed.json) records the environment and source
hashes. Its separately labeled pre-change baseline used single cold runs;
use the paired warmed measurements for speed comparisons.

On the Intel development Mac, this report measured:

| Fixture | Hybrid median (s) | Velocity median (s) | Hybrid / velocity |
| --- | ---: | ---: | ---: |
| Static column, 10 steps | 0.151 | 0.141 | 1.07x |
| Coarse collapse, 60 steps | 0.789 | 1.005 | 0.78x |
| Refined collapse, 20 steps | 4.991 | 3.299 | 1.51x |
| Static incline, 10 steps | 0.116 | 0.480 | 0.24x |
| Sliding incline, 10 steps | 2.105 | 2.080 | 1.01x |
| Plane-strain collapse, 60 steps | 1.014 | 1.386 | 0.73x |

Small-case overhead and convergence sensitivity matter: tiny floating-point
changes can trigger the reference's block fallback on an inclined layer.
The report retains individual samples and iteration counts. These bounded
measurements support keeping the hybrid default; they establish neither a
general CPU speedup nor Apple-silicon/GPU performance.

Every case also replays the **same frozen last-step problem**, including the
reference warmstart, and checks the original residual and grid-velocity
difference. Evolving runs compare particle positions and velocities. Stresses
may differ in dependent modes even when the velocity solution agrees; this is
not a bitwise-equivalence optimization.

The sliding incline is sensitive to stopping tolerance: at `1e-9`, both solvers
meet the original residual criterion while their frozen grid velocities can
differ by several `1e-7`. The report retains that failed `2e-8` velocity bound
and adds separate `1e-11` and `1e-12` studies of both solvers on the identical
problem. Residual acceptance, matched-tolerance velocity agreement, and
velocity agreement after tolerance refinement are reported separately. Tight
reference agreement does not imply the default `1e-9` residual uniformly bounds
velocity errors by `1e-9`.

Tests cover compact assembly against the original operator and transpose,
periodic coloring, local opening/sticking/sliding, coupled frozen problems,
static support, free fall, both collapse stress spaces, inclined flow,
segregation conservation, viscosity work, saved output and explicit failure
under an insufficient iteration budget.

This pass does not correct the packing drift discussed in the
[source audit](granular-source-audit.md). Q1–P0 is a separate next experiment:
it changes the spatial discretization and requires new accuracy and packing
studies before a GPU port. Particle-carried stress warmstarts also remain a
separate experiment; current grid stress warmstarts are retained for comparison.
