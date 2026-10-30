# Design note: a well-posed granular step

The rigid-grain solver cannot refine connected flows: with constant friction, a
hard packing limit and no rate dependence the continuum model is ill-posed, and
the discrete problem stalls every solver tried
([literature review](connected-flows-literature.md),
[regularisation study](granular-mixed.md#regularised-laws-under-refinement)).
This note sets out the model, the time step and the solver that would replace
it, what the code already provides and what must be built or decided first. It
began as a plan; the step now runs in the MPM code on the CPU
([Production solver](#production-solver)). Prototypes of the
step are verified in one dimension ([Shear-cell results](#shear-cell-results))
and on the frozen first step of the column collapse
([Frozen-step results](#frozen-step-results)).

## Model

The candidate is the compressible I-dependent rheology (CIDR) of Barker,
Schaeffer, Shearer and Gray (2017), in one of two published forms:

- **Critical-state CIDR** (Barker et al. 2017; parameters as in Fannon, Moyles
  and Fowler 2019): yield function `Y = alpha(I) p - p^2 / C(phi)` and
  dilatancy function `f = beta(I) - 2 p / C(phi)`. It has a static yield stress,
  which the drum and the column collapse need, and keeps rigid grains.
- **Inertial CIDR** (Schaeffer et al. 2019): `Y = mu(Psi) (I / Psi) p` and
  `f = (mu(Psi) / 4)(I / Psi - Psi / I)`, with `Psi` the inverse of the steady
  `phi = Phi(I)`. It is the better supported model but has no static yield
  stress and covers only packings below the critical fraction.

In both, deviatoric stress and deviatoric strain rate are aligned,
`|tau| = Y(p, phi, I)` while the material deforms, and the flow rule is
`div u = 2 f(p, phi, I) |D|`, with `|T| = sqrt(T_ij T_ij / 2)`, `D` the 2D
deviatoric strain rate and `I = 2 d |D| / sqrt(p / rho_s)`. Barker et al. prove
the equations linearly well-posed in 2D when

```text
(a)  Y_p - (I / 2p) Y_I = f + I f_I
(b)  Y_I > 0
(c)  f_p - (I / 2p) f_I < 0
```

### Questions to settle before implementing

1. **Regime.** Static material needs a yield stress. Critical-state CIDR has
   one; inertial CIDR would need a quasi-static extension that has not been
   published.
2. **Packing.** In critical-state CIDR the packing limit is soft: `C(phi)`
   diverges as `phi` approaches `phi_max`. Its scale `P / b` is free, and the
   scale `rho_s g d` used by Barker et al. is not frame-indifferent; Schaeffer et
   al. (2019) suggest a stiffness scale `k / d`. The current hard packing
   constraint would become a stiff constitutive term.
3. **Variational structure.** The time step below needs the stress to derive
   from a convex dissipation potential. Schaeffer et al. show this for inertial
   CIDR; for critical-state CIDR it is derived below, and holds only once the
   flow rule is completed at the cap of the yield surface (see
   [Derivations](#derivations)).
4. **Discrete well-posedness.** The proofs are for the linearised continuum
   equations. Whether a time step that lags `phi` keeps the property, and
   whether the discrete solutions converge under refinement, are to be tested.

## Time step

Work in velocity. Each step solves

```text
min_u  1/2 |u - u*|^2_(M / dt) + sum_q V_q psi(D_q(u); phi_q^n) + wall terms
```

where `u*` includes gravity, `psi` is the convex dissipation potential
evaluated at quadrature points `q` with `phi` from the start of the step, and
the wall terms are the convex contact constraints already used for one-sided
and frictional walls. The mass term makes the problem strictly convex, so the
velocity is unique even where the stress is not: the static indeterminacy that
defeats the present stress-space iteration does not enter. Stress follows as a
subgradient (or as the dual variables of the solver), unique wherever the
material deforms. The packing fraction is then updated from mass conservation
(`phi = rho / rho_s`, already evolved through `J`).

Rate dependence enters through `D` inside `psi`, so it is implicit; only the
`phi` dependence is lagged.

## Discretisation

Velocity on the background grid as now; `psi` integrated at particles, or per
cell if locking appears near the packing limit. Because stress is not a primal
unknown, velocity–stress compatibility matters less than in the present mixed
form. Volumetric locking near jamming may need volumetric averaging, which the
explicit solver already has (`fbar`).

## Solver

A convex solver whose iteration count does not grow with the grid:

- **Interior point on the CPU.** In the regularisation study, Clarabel solved the
  associated (convex) first step in 25 to 41 iterations from 231 to 13,041 stress
  nodes. Cost per iteration is a sparse factorisation; not a GPU method.
- **Multilevel nonsmooth Newton for scale and the GPU.** Truncated nonsmooth
  Newton multigrid has mesh-independent iteration counts for convex plasticity
  problems (Gräser and Sander 2019; Sander 2017): block Gauss–Seidel smoothing,
  a Newton correction on the smooth part solved by one multigrid V-cycle, and a
  line search.

Both start from the previous step's solution and can use coarse-to-fine nested
iteration.

Jammed regions constrain the choice of smoother. Where a cell or quadrature
point jams, `psi` has a cone at `D = 0`, and `D` is a difference of nodal
velocities, so the kink is not separable in the nodal unknowns. Gauss–Seidel over
nodes then cannot move a jammed group of nodes as one body and can stop short of
the minimum. In the one-dimensional shear cell this is avoided by taking the
cell strain rates as unknowns, with the walls as a linear constraint
([Shear-cell results](#shear-cell-results)). That change of variables has no
direct analogue on a two-dimensional grid, where strain rates must be compatible.
The two-dimensional prototype keeps strain-rate unknowns at the stress nodes with
compatibility as a constraint, and solves either by augmented-Lagrangian
splitting or by the method of multipliers with semismooth Newton subproblems on
the Moreau envelope of `psi`. Both converge; neither has an iteration count
independent of the grid ([Frozen-step results](#frozen-step-results)).
Smoothing over aggregates of nodes and primal–dual methods on the product domain
of the fixed cap are not yet tested.

## Verification

1. One-dimensional gravity-free shear cell against the steady solution and the
   transient DEM results of Schaeffer et al. (2019).
2. The frozen first-step refinement study, extended to the new model.
3. Steady flow on an inclined plane against the Bagnold profile: held for
   0.4 s ([Inclined plane](#inclined-plane)), not yet beyond.
4. Refined column collapse against experiments.
5. The rotating drum, with segregation.

## Inventory

| Piece | Status |
| --- | --- |
| MPM core (updated `F`, `J` and density; linear and quadratic splines; APIC; finite-domain quadrature; budgets; checkpoints) | Available |
| Per-particle packing fraction evolving through `J` | Available; currently used only as a hard cap |
| Walls (one-sided, frictional), periodic edges | One-sided walls by re-solving; Coulomb friction as a Tresca term whose threshold follows the normal reaction between outer steps (`cidr.Friction`, `granular.wall_friction`), on both backends |
| Velocity-space convex step with a cone yield | Prototype (study script, frozen first step) |
| Velocity-space convex step with critical-state CIDR in 2D | Prototype (`benchmarks/cidr_frozen_step.py`, frozen first step): converges under refinement; solvers converge but not in a number of iterations independent of the grid ([Frozen-step results](#frozen-step-results)) |
| Nested coarse-to-fine warm start | Prototype (study) |
| Volumetric averaging (`fbar`) | Available in the explicit solver |
| Segregation coupled to the granular step | Available (CPU and Metal) |
| Dissipation potential for the chosen CIDR model | Derived; checked symbolically and numerically (`benchmarks/cidr_potential.py`) |
| Rigid-grain step with critical-state CIDR | `granular_solver = "cidr"`, CPU only ([Production solver](#production-solver)) |
| Newton solver on the potential with line search | Prototype in 1D (`benchmarks/cidr_shear_cell.py`): Newton for the inertial model, truncated nonsmooth Newton in cell strain rates for the critical-state model |
| Multigrid hierarchy on the MPM grid | Missing |
| GPU convex solver | Missing; the Metal port solves the stress-space problem this replaces, though its operators, binning, transfers and segregation carry over |
| Shear-cell verification against published results | Done for both models: grid convergence and lagged packing ([Shear-cell results](#shear-cell-results)); not yet compared point by point with the published transient curves |
| Inclined-flow verification | Missing |

## Derivations

`benchmarks/cidr_potential.py` checks each result below (sympy and scipy
required). Write `s = |D|`, `delta = div u` and `k = 2 d sqrt(rho_s)`, so that
`I = k s / sqrt(p)`, and consider a function `g(s, p)` of shear rate and
pressure. Because `|T| = sqrt(T_ij T_ij / 2)`, a stress `tau = dg/dD` has
`|tau| = g_s / 2`.

### Well-posedness conditions as a potential

Ask for `g` with `g_s = 2 Y` (the yield condition with alignment) and
`g_p = 2 f s` (the flow rule). Such a `g` exists exactly when the mixed
derivatives agree, which with `I = k s / sqrt(p)` reads

```text
Y_p - (I / 2p) Y_I = f + I f_I,
```

condition (a) of Barker et al. (2017). Its second derivatives are

```text
g_ss = 2 Y_I I / s,      g_pp = 2 s (f_p - (I / 2p) f_I),
```

so conditions (b) and (c) are exactly convexity of `g` in shear rate and
concavity in pressure. With `Y >= 0`, `g` is also nondecreasing in `s`, hence
convex in the tensor `D`. The linear well-posedness conditions of the CIDR
framework therefore coincide with the existence of a saddle potential, convex in
`D` and concave in `p`, from which `tau = dg/dD` and `delta = dg/dp`. The
dissipation potential is its partial conjugate,
`psi(D, delta) = sup_p [g(|D|, p) - p delta]`. All three statements are checked
symbolically for general `Y(p, I)` and `f(p, I)`.

Applied to inertial CIDR, the construction gives
`g = 2 d sqrt(rho_s) M s^2 sqrt(p) / Psi - M Psi p^(3/2) / (6 d sqrt(rho_s))`
with `M = mu(Psi)`, which is eq. (B8) of Schaeffer et al. (2019) (checked
symbolically).

### Critical-state CIDR

With Barker et al.'s functions, `alpha - beta / 2 = mu` and their differential
equation `alpha + (5/2) I alpha' = 2 mu + 2 I mu'`, the potential is

```text
g(s, p) = 2 p s [ (5 alpha(I) - 4 mu(I)) / 3 - p / C(phi) ],
```

checked symbolically against `g_s = 2 Y` and `g_p = 2 f s` for general
`mu(I)`, and by finite differences with the `mu(I)` of Jop et al. for glass beads
(relative errors 5e-8 and 8e-8). At zero shear its slope is
`2 p (2 mu_s - p / C)` (`alpha(0) = 2 mu(0) = 2 mu_s`, checked numerically), so `g` has a
corner at `s = 0`: the static yield condition `|tau| <= 2 mu_s p - p^2 / C`. On
`1e-6 <= I <= 10`, `alpha' > 0` and `beta' > 0`, so `Y_I = alpha' p > 0` and
`f_p - (I / 2p) f_I = -2 / C - (I / 2p) beta' < 0`, as Barker et al. prove.

### The cap

`Y = alpha(I) p - p^2 / C` vanishes at `p = alpha(I) C`, which increases with
shear rate from `alpha(0) C = 2 mu_s C` at rest. Without a bound on `p`, `psi`
is not convex in `D` (it decreases with shear rate where `Y < 0`), and pure
compaction has unbounded dissipation: Barker et al.'s flow rule
`div u = 2 f |D|` allows no compaction without shear, so an unsheared column
could not consolidate under its own weight. The flow rule therefore has to be
completed by a bound on `p`, and the bound has to keep `psi` convex.

The bound used is the fixed cap `p <= alpha(0) C = 2 mu_s C`. Then
`psi(s, delta) = sup_{0 <= p <= 2 mu_s C} [g(s, p) - p delta]` is a supremum
over a fixed set of functions convex in `(s, delta)` (`g_ss = 2 Y_I I / s >= 0`),
so it is convex; and `psi_s = g_s = 2 Y >= 0` because `p <= alpha(0) C <=
alpha(I) C`, so it is nondecreasing in `s` and convex as a function of `D`.

The shear-dependent bound `p <= alpha(I) C`, where `Y = 0`, does not give a
convex potential, although the region it bounds in `(s, p)` is convex. On its
cap branch `psi = g(s, p_cap(s)) - p_cap(s) delta`, whose Hessian in
`(s, delta)` has determinant `-p_cap'(s)^2 < 0`. With the glass-bead `mu(I)`
and `C = 4e4` Pa, `I` stays below about 0.02 and the defect is below the
tolerance of a sampled midpoint test; at grain scale (`C = 1`, `I` of order one)
the midpoint test fails by 0.8% of `psi`, against 8e-6 (the optimiser's
tolerance) for the fixed cap. The shear cell exposed it
([Shear-cell results](#shear-cell-results)).

With the fixed cap, compaction at zero shear proceeds at the cap pressure, as
normality gives at the cap in critical-state soil mechanics, and compaction while
shearing does too, with shear stress `Y(2 mu_s C, I) = 2 mu_s C (alpha(I) -
alpha(0))`. Completing the flow rule this way is an addition to the published
model and has to be stated as such; it differs from the shear-dependent bound
only for states with `p > 2 mu_s C`, which are reached only while shearing and
compacting at once.

### States at fixed packing

With `phi`, hence `C`, fixed (checked numerically with the glass-bead `mu(I)`):

| State | Condition | `C` needed to carry pressure `p` |
| --- | --- | --- |
| Steady isochoric flow | `f = 0`, then `|tau| = mu(I) p` (to 1e-16) | `2 p / beta(I)` |
| Critical state (`I -> 0`) | `f = 0`, `I = 0` | `p / mu_s` |
| Loosest state at rest | on the cap at `s = 0` | `p / (2 mu_s)` |

The critical pressure is half the cap pressure at the same packing, as in
modified Cam-clay. Material on the cap (loose for its load) compacts when
sheared; material denser than critical dilates.

### Scale of `C`

With the form of Fannon et al., `C = (P / b) x`, `x = psi / (1 - psi)` and
`psi = (phi - phi_min) / Delta`, `Delta = phi_max - phi_min`, the packing on
the cap at pressure `p` falls short of `phi_max` by

```text
e(p) = Delta / (1 + x),   x = b p / (2 mu_s P),
```

and the effective bulk stiffness there is `K = phi dp/dphi ~ phi p / e`
(checked to 0.8%). The scale `P / b` is therefore a choice of compressibility:
it fixes how far below `phi_max` the material sits at a given load, and the
stiffness follows. A rule for choosing it is a packing deficit `e_ref` at a
reference pressure `p_ref`:

```text
P / b = p_ref / (2 mu_s (Delta / e_ref - 1)).
```

For `e_ref = 0.005` at the pressure under 0.1 m of grains (1,470 Pa,
`Delta = 0.6`), `P / b = 16 Pa`, and the deficit grows to 0.047 at a tenth of
that pressure: a loose layer forms near the free surface, where `p -> 0`.

This is a physical choice as well as a numerical one. The steady packing of
critical-state CIDR depends on `p / (P / b)`, whereas for rigid grains in the
inertial regime DEM gives a packing that depends on `I` alone (Schaeffer et al.
2019). A pressure scale has to come from somewhere; Schaeffer et al. identify it
with grain stiffness, `k / d`. Here it can be read either way: as the stiffness
of real grains, or as a small, stated compressibility added to rigid grains.

### Stress at zero shear

At `D = 0` the potential is `psi(0, delta) = sup_{0 <= p <= 2 mu_s C} (-p delta)`,
and its subdifferential gives the stress:

- `delta > 0` (separating): `p = 0` and `tau = 0`, stress-free, as in the present
  model's opening;
- `delta < 0` (compacting): `p = 2 mu_s C` (the cap) and `tau = 0`;
- `delta = 0` (static): any `p` in `[0, 2 mu_s C]` and any `tau` with
  `|tau| <= 2 mu_s p - p^2 / C`.

In static zones the stress is set-valued; equilibrium selects it globally and
may not select it uniquely. The velocity is unique regardless. The stress used
for output can be taken from the solver's dual variables, which depend on the
algorithm, or made unique by an added criterion (minimum norm, or a small
elastic compliance). Segregation needs no such choice: its potential uses
`p I`, and `I = 0` where `D = 0`.

### Evaluating `psi`

At a quadrature point with `s > 0`, `psi(s, delta) = g(s, p*) - p* delta`,
where `p*` maximises `g - p delta` on `[0, p_cap]`. Because `g` is concave in
`p`, `p*` is unique and found by a bracketed scalar root of `g_p(s, p) = delta`,
with three branches:

| Branch | Condition | `p*` | `psi_s` (= `2 |tau|`) | `psi_delta` (= `-p`) |
| --- | --- | --- | --- | --- |
| Interior | `g_p(s, 0) > delta > g_p(s, p_cap)` | root of `g_p = delta` | `2 Y(p*, I*)` | `-p*` |
| Cap | `g_p(s, p_cap) >= delta` | `p_cap = 2 mu_s C` | `2 Y(p_cap, I)` | `-p_cap` |
| Zero pressure | `g_p(s, 0) <= delta` | `0` | `0` | `0` |

On the interior branch the Hessian in `(s, delta)` is

```text
psi_ss = g_ss - g_sp^2 / g_pp,   psi_sdelta = g_sp / g_pp,   psi_deltadelta = -1 / g_pp,
```

with determinant `-g_ss / g_pp`. Since `g_ss >= 0` and `g_pp < 0` by conditions
(b) and (c), it is positive definite: `psi` is convex on the interior branch
analytically, not only by sampling. On the cap branch `psi_ss = g_ss >= 0` and
the other second derivatives vanish. Gradients on all three
branches agree with finite differences to 1e-8 or better, and the Hessian to
2e-4 (the remaining difference shrinks with the step, so it is
finite-difference error).

The stress tensor follows from `s = |D|`: `tau_ij = psi_s D_ij / (2 s)`, so that
`|tau| = psi_s / 2`, and `p = -psi_delta`.

### Consequences for the time step

- The step is the convex minimisation of the [Time step](#time-step) section with
  `psi` from the capped potential; `phi`, and hence `C`, is held at its value from
  the start of the step.
- With the fixed cap the saddle form `min_u max_{0 <= p <= 2 mu_s C}` has a
  product domain, so primal–dual methods apply. The primal form in `u` alone,
  with `psi` evaluated at each quadrature point by a one-dimensional concave
  maximisation in `p` (no closed form was found), is the one used in the shear
  cell.
- `psi` is nonsmooth at zero shear (a cone: the static yield), and its Hessian
  jumps where the cap becomes active. The solver must handle the cone directly;
  plain Newton stalls once a region jams
  ([Shear-cell results](#shear-cell-results)). `alpha(I)` has no conic
  representation, so the interior-point route used in the regularisation study
  does not carry over directly.

### Open points

- The value of the packing deficit `e_ref`, and hence `P / b`: a modelling
  choice, as above. The form of `C` near `phi_max` is speculative in the
  published model.
- Whether holding `phi` fixed within a step preserves well-posedness. This is
  not settled analytically. In the shear cell lagged and implicit `phi` differ by
  two to three times the time discretisation error, with no change in behaviour
  ([Shear-cell results](#shear-cell-results)).

## Shear-cell results

`benchmarks/cidr_shear_cell.py` solves the gravity-free shear cell of Schaeffer
et al. (2019) with the time step of this note on a Lagrangian one-dimensional
grid. Case A is their (3.16) (height 30, `phi = 0.78`, a localised
vertical-velocity perturbation of amplitude 0.05, to `t = 20`); case B is their
(5.7) (height 100, sinusoidal perturbations of amplitude 0.16 in `u` and 0.1 in
`w`, to `t = 1000`). The critical-state model uses the fixed cap and a packing
deficit of 0.005. Results are in `benchmarks/cidr_shear_cell.json`.

### Solvers

The inertial model has no static yield, and Newton's method on the nodal
velocities with a backtracking line search takes 1.2 to 1.3 iterations per step.

For the critical-state model plain Newton stalls in case B. By `t = 8` a layer
has compacted to `phi = 0.8161` (46 to 400 cells agree on the packing range to
1e-4), by `t = 10` to 0.8164, within 1e-3 of `phi_max`, and cells in it stop
deforming. There
`psi` sits at the cone of [Stress at zero shear](#stress-at-zero-shear), the
gradient Newton evaluates is one arbitrary subgradient, and on 200 cells the
iteration makes no progress. Taking the cell strain rates `d_c = (u', w')` as
unknowns makes the kink separable, with the walls as the linear constraint
`sum_c dz_c d_c = (1, 0)`. Each iteration then has two parts:

- An exact block Gauss–Seidel pass over slow and jammed cells. Each cell is moved
  together with the fastest cell, which keeps the constraint. A cell is jammed
  if and only if the stress across it lies in the static set: the union over
  `0 <= p <= 2 mu_s C` of the discs of radius `Y_0(p) = 2 mu_s p - p^2 / C`
  centred at `(0, -p)`, tested through the support function
  `h(e) = sup_p [Y_0(p) - p e_delta]`, which agrees with `psi(t e) / t` to order
  `t`. Cells that are not jammed are minimised locally, starting from just
  outside the cone if that is lower: near `D = 0` the Hessian of `psi` grows like
  `1 / |D|`, and a Newton step does not move a cell away from there.
- A Newton step on the free cells with the jammed cells held at zero (a
  truncated Newton step), solved as a saddle-point system with the constraint,
  and a backtracking line search on the whole objective.

Where nothing jams this reproduces Newton on the nodal velocities to 1e-16 (case
A). The shear cell also exposed the nonconvexity of the shear-dependent cap
([The cap](#the-cap)): with it, `psi` had negative curvature along compaction
near `phi_max`.

### Time step and packing limit

With `phi` lagged and `C` growing without bound towards `phi_max`, a step can
compact a cell past `phi_max`. Each step is first limited so that no cell closes
more than half its remaining gap `(phi_max - phi) / phi` at its current
compaction rate, then checked after the solve and halved if any cell closed more
than half its gap, or if an implicit pass produced `phi >= phi_max`. Without the
check, the implicit iteration left the admissible range in case B.

### Grid convergence and lagged packing

Grid convergence compares `max |w|` from 46 and 200 cells at 40 times spaced
logarithmically:

| Model | Case | Largest difference | Signal | Iterations per step (46 / 200 cells) | Final `max |w|` |
| --- | --- | --- | --- | --- | --- |
| Inertial | A | 3.1e-4 | 0.021 | 1.2 / 1.3 | 1.8e-5 |
| Inertial | B | 1.3e-3 | 0.10 | 1.2 / 1.3 | 1.3e-4 |
| Critical-state | A | 5.6e-4 | 0.021 | 1.0 / 1.1 | 7.7e-4 |
| Critical-state | B | 4.3e-3 | 0.10 | 1.5 / 3.1 (at most 6 / 43) | 6.4e-5 |

Inertial case A relaxes to the uniform shear by `t = 20`, consistent with
Schaeffer et al.'s report that the flow is steady by about `t = 10`. Critical-state
case A retains a small residual vertical velocity at `t = 20`. In critical-state
case B the iteration count grows with the grid because more cells jam and
release on the finer grid; the 200-cell run took 154 s.

On case B with 46 cells, against a signal of 0.10 in `max |w|`:

| Model | Lagged against implicit `phi` | Lagged at the full step against half the step |
| --- | --- | --- |
| Inertial | 4.6e-4 | 1.8e-4 |
| Critical-state | 1.2e-3 | 4.5e-4 |

Lagging `phi` adds an error of about two to three times the time discretisation
error, about 1% of the signal for the critical-state model, and changes neither
model's behaviour.

### Consequences for the two-dimensional step

- The step needs the compaction limit and the check after the solve.
- The critical-state model needs a solver that treats the cone at `D = 0`
  directly; see [Solver](#solver) for why the one-dimensional change of
  variables does not carry over.
- The fixed cap is required for convexity.

## Frozen-step results

`benchmarks/cidr_frozen_step.py` solves the first step of the column collapse from
rest (width 0.5 m, height 1 m, `dt = 0.002` s) with critical-state CIDR in two
dimensions, on the operators of the rigid-grain solver: velocities at the grid
nodes, and `psi` at the nodal stress nodes with the strain rate averaged over the
grains in each node's support. The packing is uniform (`phi = 0.6`) and frozen,
and the glass-bead `mu(I)` of Jop et al. is used with a grain diameter of 3 mm.
Results are in `benchmarks/cidr_frozen_step_alg2.json` and
`benchmarks/cidr_frozen_step_newton.json`.

`C` sets how far the column is from the cap. With the base on the cap at rest
(`2 mu_s C = rho g H`), the base compacts at the cap pressure, which only balances
the column's weight, and the column falls through the step as if unsupported (probe
speeds of `g dt`). The results below use a cap pressure of twice the base
pressure, which puts the base at the critical pressure and the whole column
denser than the cap for its load. The column then collapses at the free face,
with probe speeds up to `0.83 g dt`.

### Refinement

Nodal velocities are compared at 1,176 probe points on successive grids, as in the
[regularisation study](granular-mixed.md#regularised-laws-under-refinement):

| Spacing (m) | Stress nodes | Change from coarser grid | Associated cone (earlier study) |
| --- | --- | --- | --- |
| 0.05 | 231 | | |
| 0.025 | 861 | 4.5% | 5.9% |
| 0.0125 | 3,321 | 2.1% | 2.6% |
| 0.00625 | 13,041 | 0.90% | 1.1% |

The change halves with the spacing: critical-state CIDR with the fixed cap
converges at first order in this test, like the associated cone. The fraction of
jammed stress nodes rises from 5% to 23% with refinement (at the ALG2 tolerance).

### Solvers

The cone of `psi` at `D = 0` is not separable in the nodal velocities, and a
Newton polish on a jammed set taken from ALG2 failed: the conditions `D = 0` at
neighbouring jammed nodes are linearly dependent, and Newton stalls at free nodes
that creep close to the cone, where the Hessian of `psi` grows like `1 / |D|`.
Two solvers keep an auxiliary strain rate at each stress node instead, and share
an exact proximal step of `psi` per node (checked against the stationarity
condition to 3e-8 and the static set on 400 random inputs).

- **ALG2**, over-relaxed (1.8) with residual balancing of the penalty, to a
  relative residual of 1e-5 (which changes the probe velocities by 4e-4 against
  1e-6).
- **Method of multipliers with semismooth Newton.** Each outer step minimises
  `f(u) + sum W psi_gamma(D(u) + gamma lambda)` over `u`, where `psi_gamma` is the
  Moreau envelope of `psi`, whose gradient is `(x - prox(x)) / gamma` and whose
  generalised Hessian is `I / gamma` at jammed nodes; then
  `lambda = grad psi_gamma(D(u) + gamma lambda)`. The limit is the step with
  `psi` itself. `gamma` falls by `sqrt(10)` per outer step from `dt / rho` to
  `dt / (1000 rho)`. With `gamma` reduced alone and no multiplier, the error
  falls only in proportion to `gamma`, and Newton fails on round-off before it
  is small.

| Spacing (m) | ALG2 iterations | ALG2 time | Newton iterations (outer steps) | Newton time |
| --- | --- | --- | --- | --- |
| 0.05 | 1,070 | 16 s | 93 (18) | 4 s |
| 0.025 | 2,628 | 60 s | 177 (10) | 18 s |
| 0.0125 | 5,929 | 279 s | 325 (9) | 143 s |
| 0.00625 | 11,725 | 1,591 s | not run | |

The two solvers agree to 2e-5 (0.05 m) and 6e-5 (0.025 m) in the probe velocities
against ALG2 at 1e-6, and give the same changes under refinement. Newton is two to
four times faster, but its iteration count also grows, by about 1.9 per halving
of the spacing against 2.2 for ALG2. In the hardest outer step on the 0.025 m grid
the line search cuts Newton's step to about 3% for some thirty iterations while
the jammed set overshoots (from 245 to 330 nodes) and then settles at 85; after
that Newton converges quadratically. The generalised Hessian jumps from the
curvature of `psi` to `1 / gamma` where a node jams, so each step jams too many
nodes, and the number that change grows with the grid. This is a matter of
globalisation (line search, path in `gamma`), not of the linear solver. Starting
each grid from the coarser grid's velocity shortens only the first outer step,
and on the 0.0125 m grid Newton then stalled at the round-off allowance at the
smallest `gamma`.

### Warm-started steps

In a run each step starts close to its solution. `benchmarks/cidr_warm_steps.py`
repeats the step 20 times on the frozen geometry, a proxy for time stepping: the
trial velocity of each step is the previous velocity plus one step of gravity,
and the multiplier (the stress) carries over. Steps after the first keep `gamma`
fixed at 1e-8, about `dt / (130 rho)`; reducing it further, as the cold path does,
made Newton fail on round-off at the 0.025 and 0.0125 m grids.

| Spacing (m) | First step (Newton iterations, time) | Later steps (Newton iterations) | Mean time of later steps |
| --- | --- | --- | --- |
| 0.05 | 88, 4 s | 2–3 (one 6) | 0.12 s |
| 0.025 | 164, 17 s | 2–4 (one 12) | 0.39 s |
| 0.0125 | 276, 104 s | 3–28 (4–6 in the first seven) | 7.8 s |

Warm steps agree with independently converged solves of the same step to 1e-5
to 2e-4 (0.05 m, against cold Newton solves at 1e-7) and 2e-5 to 5e-5 (0.025 m,
against ALG2 at 1e-7), without accumulating over the steps. ALG2 gains little
from warm starts (about 900 iterations per step at 0.05 m). On the finest grid
the warm count rises in later steps, where the frozen geometry has accelerated
to over 20 `g dt`; whether this happens in a moving run is not yet measured.
Each Newton iteration costs a sparse factorisation and the proximal step at
every node in Python, about 0.6 s at 3,321 nodes.

### Production solver

`granular_solver = "cidr"` (`cidr.py`, `granular.cidr_contact`) uses this solver
in the rigid-grain step on the CPU: packing per stress node from the particles,
the previous step's stress as the multiplier, a cold start only on the first
step. The velocity is the solver's own; rebuilding it from the stress, as the
other solvers do, added Newton's residual (1e-5 m/s per step) and made a static
column drift. A column at rest between walls then moves by less than 2e-10 m
over 20 steps of 0.02 s, with five Newton iterations per step. The collapse of
the 0.5 by 1 m column on the coarse test grid (30 steps of 0.01 s, 4.6 s in all)
takes 31 iterations on its cold first step, 8 to 17 on the next five while the
flow starts, then 3 to 7 (`tests/test_cidr.py`).

The compaction limit of the shear cell carries over: after each solve, the
largest ratio of a node's compaction rate to its relative packing gap
`(phi_max - phi) / phi` sets a `packing` timestep limit for dynamic stepping,
and a step that would close a gap stops with an error. In the collapse at
`dt = 0.02` s one step closes 91% of a gap; at `0.05` s the fourth step would
close 173% and is stopped; dynamic stepping at `cfl = 0.2` keeps every step
after the first below half a gap (the limit uses the previous step's rate).
Nodes at or past `phi_max` are left out of the limit. `C` is clipped there to
about 1e12 times its scale, so their compaction is round-off, and they reach
`phi_max` as particles are advected between nodes rather than by a step's
compaction. On the [incline](#inclined-plane), FP32 round-off at such a node
made the closing rate 1e22 per step before they were left out.

The envelope Hessian is formed as `r H (H + r I)^-1` rather than
`r (I - r (H + r I)^-1)`: at the CPU's smallest `gamma` the second form
subtracted nearly equal matrices, lost about `log10(r / |H|)` digits, and
stalled Newton in the incline's first step.

### GPU solver

`granular_solver = "cidr"` also runs on Metal in FP32 (`mpm_gpu/cidr.py`): the
same method, with Jacobi-preconditioned conjugate gradients for each Newton
system and line searches on the sign of the directional derivative, which FP32
resolves where it cannot resolve the objective's decrease. On the benchmark
collapses it is 8.5 to 44 times faster than the CPU solver and agrees with it to
6e-5 grid spacings; the connected collapse of 10,000 particles, on which the
Coulomb solver does not converge, runs at 0.54 s per step. Details and timings
are in [gpu](gpu.md#critical-state-cidr).

## Inclined plane

`inputs/cidr_incline.py` is a layer 30 grains deep (3 mm glass beads, the
`mu(I)` parameters of the shear cell) on a no-slip base at 25 degrees,
periodic along the slope, 18 cells deep with four particles per cell. In
steady isochoric flow CIDR reduces to `mu(I)`: every depth flows at `I_theta`
with `mu(I_theta) = tan(theta)` (0.133 here), so the shear rate is
`I_theta sqrt(p / rho_s) / d` with `p` the weight above (the Bagnold profile).
The layer starts in this state, at the critical packing of each depth from
`C(phi) = 2 p / beta(I_theta)`, loose only within about a grain of the surface,
and at the Bagnold velocity (1.73 m/s at the surface). A layer at uniform
packing 0.6 instead formed a jammed plug over a basal shear band: near the
surface its static strength `2 mu_s - p / C` exceeds `tan(theta)`.

`benchmarks/cidr_incline.py` runs it on Metal and compares the shear rate in
18 layers, one per cell, with the prediction from the simulated weight above
(`benchmarks/cidr_incline.json`). With these layers the initial state agrees
to 0.7% below the second layer. After 0.4 s (6,888 steps, 43 s) the surface
moves at 1.720 m/s and the shear rate of layers 3 to 15 is on average 0.990 of
the prediction, within 8.9% in every layer; the two layers at the base,
where the no-slip wall holds the velocity, are 17% and 40% above it.

The state does not persist beyond about 0.4 s. The spread of particle packing
in the lower half of the layer, 0.594 to 0.600 at 0.2 s, is 0.538 to 0.651 at
0.4 s, and near 0.45 s a stress node passes `phi_max` and the step is stopped
(the CPU behaves the same: by 0.5 s the densest particle is at 0.78). The
solver controls only the volumetric rate of each stress node, a solid-weighted
average over its particles, and under sustained shear the velocity field
moves particles together and apart below that resolution. Two remedies were
tried and rejected:

- Updating particle volumes from the nodes' volumetric rate interpolated to
  the particles (as B-bar) kept their packing within 0.597 to 0.600, but the
  positions still followed the velocity field: the packing measured from the
  occupied area rose to 1.16 by 1.5 s while the particle volumes stayed at
  0.55, and the flow stopped.
- Packing from the nodes' basis volumes, which sees particle positions,
  failed on the first step: the initial state is set from particle packing,
  and partly filled surface nodes read much looser.

A remedy needs either a stabilisation of these modes or a packing estimate
that follows the particle positions and is consistent at the free surface; it
is not yet chosen.

## References

The papers are listed in the [literature review](connected-flows-literature.md).
