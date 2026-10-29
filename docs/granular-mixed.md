# Compact velocity-based granular solve

Set `P.granular_solver = "velocity"` before `prepare(P)` to select the opt-in
CPU implementation. `P.solver = "granular"` still selects the rigid-grain
equations; `granular_solver` selects the numerical algorithm. The default
remains `"hybrid"`, which is the small-problem reference solver.

This pass adapts the stress iteration and nonassociated local flow-rule solver
from [Daviet (2026), sections 5.2–5.3](https://research.nvidia.com/labs/prl/mixed_mpm/files/mixed_mpm.pdf).
The first pass retains our bilinear velocity and **nodal** stress spaces, particle
quadrature, packing gap, transfer options, wall constraints and viscosity split.
The separate Q1–P0 experiment below changes the stress space. No GPU backend
is implemented yet.

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
hashes for the nodal first pass at commit `c8eefdd`; the subsequent P0 study
has its own source-fingerprinted report. Its separately labeled pre-change baseline used single cold runs;
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
[source audit](granular-source-audit.md). Particle-carried stress warmstarts
remain a separate experiment; current grid stress warmstarts are retained for
comparison.

## Experimental Q1–P0 stress space

Set `P.granular_stress_basis = "P0"` before `prepare(P)` to use one constant
stress per occupied cell. `"nodal"` remains the default. This option is
independent of `granular_stress`, which still chooses intrinsic 2D or plane
strain, and `granular_solver`, which chooses the solution algorithm.

For P0, a particle contributes its full grain volume to its containing cell,
and the stress shape function is one in that cell. The coupling and packing
clearance become

```text
B_cell,j = sum(particles in cell, Vgrain * BASIS : grad N_velocity,j)
gap_n = max(phi_max * cell_volume - sum(particles in cell, Vgrain), 0)
        / (sqrt(dimension) * dt)
```

The volume is the geometric cell volume, including thickness. Empty cells have
no stress unknown. Stresses are transferred to particles as cell constants,
using their positions before advection. Periodic velocities share seam nodes;
cell stresses remain separate unknowns on each side of the domain. Both the
assembled reference and compact path use this same definition. Four colors
suffice for the nonperiodic Q1–P0 interior stencil in 2D; the generic coloring
also handles periodic grids whose dimensions need additional colors.

This is an adaptation to our existing fixed-grain-volume rigid model, not a
reproduction of the complete 2026 constitutive implementation. The pressure
constraint still uses the positive-part gap and does not geometrically repair
overpacking. Cell packing also jumps when point particles cross cell edges;
it is distinct from both nodal reconstruction and deformation-based density.

Existing `maximum_packing` and `maximum_packing_after` diagnostics continue to
refer to nodal reconstruction for comparison across spaces. P0 additionally
reports `maximum_stress_packing`, `maximum_stress_packing_after`, and
`integrated_stress_overpacking_fraction`, using cell grain volumes. Material
and segregation outputs, plot formats and save callbacks remain unchanged.

Run the controlled spatial experiment with:

```sh
uv run --locked python -m benchmarks.granular_discretization --output /tmp/granular-discretization.json
```

The [report](../benchmarks/granular_discretization.json) retains warmup failures,
completed time, residuals, sampling/grid changes, packing diagnostics and
inclined-layer reference errors. Cost ratios compare different discrete
models; they are not equivalent-behavior solver speedups.

Independent scalar assembly checks, affine/constant velocity patches, adjoint
work, mass/packing inventory, static support, both stress reductions and
short evolving comparisons against the assembled solver pass. Periodic
segregation conservation and saved plot cadence are checked for P0 too.
These implementation checks do not establish spatial convergence: the refined
collapse exposes an iterative-solver failure, and coarse inclined flow has a
larger discrepancy with the sliding reference. The report includes an assembled
first-step replay when the compact iteration fails immediately, to distinguish
iteration failure from an infeasible constraint problem.

The recorded pre-density-fix checkpoint gives these bounded costs:

| Case | Nodal median (s) | P0 median (s) | Observation |
| --- | ---: | ---: | --- |
| Static column | 0.140 | 0.067 | Both remain static |
| Collapse, 2 particles per cell axis | 0.987 | 0.736 | Peak nodal packing: 0.621 vs 0.646 |
| Collapse, 4 particles per cell axis | 0.979 | 0.780 | Peak nodal packing: 0.607 vs 0.634 |
| Short coarse collapse | 0.395 | 0.280 | Same time and dt as the refined case |
| Refined collapse | 3.355 | Failed | Assembled P0 first step converges |
| Coarse incline | 3.151 | 0.591 | P0 sliding reference error: 3.72% |
| Refined incline | Failed | Failed | Neither compact solve meets tolerance |

The grid comparisons hold `dt=0.0025` and final time `0.05` fixed. The
particle-count comparison holds grid, `dt=0.005` and final time `0.3` fixed.
The successful assembled replay of the failed P0 collapse first step requires
200 least-squares evaluations plus 1,653 block sweeps, reaching residual
`9.93e-10`. Thus the compact iteration's failure is not evidence that this
frozen constraint problem has no solution. The inclination and packing results
also prevent treating P0's lower coarse-case costs as an accepted physical
accuracy/performance improvement.

P0 remains experimental. Resolve these iteration and accuracy gates, and
consider the paper's dP1 alternative, before selecting a compact stress space
for the GPU port. Neither the paper's compact mixed pairs nor this experiment
establish an inf-sup stability guarantee.

## Coarse Q1 stress space and refinement study

Set `P.granular_stress_basis = "coarse"` to place bilinear stresses on a grid
whose spacing is `P.granular_stress_coarsening` (default 2) times the velocity
grid spacing. The cell counts in x and y must be divisible by that factor. The
velocity space, rigid-grain equations, packing gap, solvers and particle stress
transfer are unchanged; the stress shape functions and their exact basis
volumes are the only differences from `"nodal"`. Stresses are transferred to
particles with the coarse shape functions at their start-of-step positions.
Nodal packing diagnostics keep their nodal reconstruction;
`maximum_stress_packing` uses the coarse basis volumes. This option is CPU
only and experimental.

The [refinement study](../benchmarks/granular_stress_scaling.json)
(`uv run --locked python -m benchmarks.granular_stress_scaling`) freezes the
first contact problem of a column collapse from rest (0.5 m wide, 1 m high,
4 particles per cell, `dt = 0.002`) and runs the compact velocity solver at
tolerance `1e-4` with a 3,000-sweep budget:

| Stress space | Grid spacing | Particles | Stress unknowns | Result |
| --- | ---: | ---: | ---: | --- |
| nodal | 0.05 | 3,200 | 693 | converged, 2,329 sweeps |
| nodal | 0.025 | 12,800 | 2,583 | 1.75e-2 after 3,000 sweeps |
| P0 | 0.05 | 3,200 | 600 | 5.3e-4 after 3,000 sweeps |
| P0 | 0.025 | 12,800 | 2,400 | 1.37e-2 after 3,000 sweeps |
| coarse, k = 2 | 0.05 | 3,200 | 198 | converged, 95 sweeps |
| coarse, k = 2 | 0.025 | 12,800 | 693 | converged at the budget (3,000 sweeps) |
| coarse, k = 2 | 0.0125 | 51,200 | 2,583 | 4.4e-3 after 3,000 sweeps |

The assembled reference solver (`hybrid`) also fails on the nodal problem at
0.025 (1.9e-2 after 200 least-squares evaluations and 2,000 block sweeps).
These are not plateaus: run to 20,000 sweeps, the nodal 0.025 residual keeps
falling slowly (9.6e-3 at 5,000, 5.7e-3 at 10,000, 2.6e-3 at 20,000). The
residuals do not change with `dt`: from rest with zero packing gaps the scaled
first-step problem is independent of the timestep.

Friction sets the difficulty. Without friction the nodal problem at 0.025
converges in 44 sweeps and the coarse problem at 0.0125 in 30; with `mu = 0.1`
they need 1,232 and 1,500 sweeps; with the input's `mu = 0.5` neither
converges within 3,000. The difficulty follows the stress resolution of the
frictional problem, not the ratio of stress to velocity unknowns (nodal at
0.025 and coarse k = 2 at 0.0125 share a 0.025 m stress spacing and behave
alike).

Holding the stress spacing at 0.05 m while the velocity grid refines, every
case converges: 3,000, 1,275 and 1,223 sweeps at 12,800, 51,200 and 204,800
particles. Nested iteration (solve on a coarser stress grid, interpolate the
stress field, which stays admissible because the Coulomb cone is convex, and
use it as the warm start) reduces the 12,800-particle case to 48 + 175 sweeps
(3 s against 46 s) and the 204,800-particle case to 41 + 197 sweeps. It does
not bring the 0.025 m stress spacing within the budget.

A fixed physical stress spacing is a regularisation: it sets a minimum length
for the stress field independent of the velocity grid. Whether it, a rheology
change such as dilatancy or a rate-dependent friction law, or another
regularisation is appropriate is a modelling decision that has not been made.
The coarse space has not been compared with the nodal space for accuracy on
evolving flows.

### Dilatancy

The [dilatancy study](../benchmarks/granular_dilatancy_study.json)
(`uv run --locked python -m benchmarks.granular_dilatancy_study`) repeats the
nodal frozen problems with a dilatant flow rule: while pressed, sliding
separates at `w_n = beta |w_t|`, so normal complementarity applies to
`w_n - beta |w_t|` (`beta = 0` is the current law, `beta = mu` associated
Coulomb flow). The exact local solve becomes the monotone scalar equation
`|r_t(alpha)| (D_n - mu beta alpha) + mu q_n = 0` with
`r_t(alpha) = -q_t / (D_t + alpha)` and `r_n = |r_t| / mu`; the prototype
matches `diagonal_coulomb` at `beta = 0` to 2e-12. The global iteration is the
colored sweeps with Anderson acceleration (no Newton). Residuals after 3,000
sweeps:

| Grid spacing | Particles | beta/mu = 0 | 0.1 | 0.25 | 0.5 | 1 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.05 | 3,200 | converged | converged | 2.6e-4 | 6.0e-4 | 4.3e-4 |
| 0.025 | 12,800 | 1.3e-2 | 1.3e-2 | 4.4e-3 | 4.4e-4 | 4.1e-4 |
| 0.0125 | 51,200 | 4.5e-2 | 3.7e-2 | 1.8e-2 | 6.1e-3 | 1.6e-3 |

At 12,800 particles the associated problem converges after 12,200 sweeps,
against 2.6e-3 after 20,000 for `beta = 0`. Dilatancy lowers the residual
reached on refined grids but slows the coarse case, and even the associated
problem, which is a convex cone complementarity problem, needs a number of
sweeps that grows rapidly with resolution. Dilatancy alone therefore does not
remove the obstacle: the colored-sweep iteration itself does not scale, and a
multilevel solver is needed whichever flow rule is used. No dilatancy law is
implemented in the model; this prototype lives only in the study script.
