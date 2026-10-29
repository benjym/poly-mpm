# Design note: a well-posed granular step

The rigid-grain solver cannot refine connected flows: with constant friction, a
hard packing limit and no rate dependence the continuum model is ill-posed, and
the discrete problem stalls every solver tried
([literature review](connected-flows-literature.md),
[regularisation study](granular-mixed.md#regularised-laws-under-refinement)).
This note sets out the model, the time step and the solver that would replace
it, what the code already provides and what must be built or decided first. It
is a plan for the MPM code, where none of it is implemented yet; a one-dimensional
prototype of the step is verified in [Shear-cell results](#shear-cell-results).

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
direct analogue on a two-dimensional grid, where strain rates must be compatible;
the options there (strain-rate unknowns with compatibility as a constraint,
smoothing over aggregates of nodes, or a primal–dual method on the product domain
of the fixed cap) are not yet tested.

## Verification

1. One-dimensional gravity-free shear cell against the steady solution and the
   transient DEM results of Schaeffer et al. (2019).
2. The frozen first-step refinement study, extended to the new model.
3. Steady flow on an inclined plane against the Bagnold profile.
4. Refined column collapse against experiments.
5. The rotating drum, with segregation.

## Inventory

| Piece | Status |
| --- | --- |
| MPM core (updated `F`, `J` and density; linear and quadratic splines; APIC; finite-domain quadrature; budgets; checkpoints) | Available |
| Per-particle packing fraction evolving through `J` | Available; currently used only as a hard cap |
| Walls as convex constraints (one-sided, frictional), periodic edges | Available in the stress-space solver; velocity form not yet |
| Velocity-space convex step with a cone yield | Prototype (study script, frozen first step) |
| Nested coarse-to-fine warm start | Prototype (study) |
| Volumetric averaging (`fbar`) | Available in the explicit solver |
| Segregation coupled to the granular step | Available (CPU and Metal) |
| Dissipation potential for the chosen CIDR model | Derived; checked symbolically and numerically (`benchmarks/cidr_potential.py`) |
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
`2 p (2 mu_1 - p / C)` (`alpha(0) = 2 mu(0)`, checked numerically), so `g` has a
corner at `s = 0`: the static yield condition `|tau| <= 2 mu_1 p - p^2 / C`. On
`1e-6 <= I <= 10`, `alpha' > 0` and `beta' > 0`, so `Y_I = alpha' p > 0` and
`f_p - (I / 2p) f_I = -2 / C - (I / 2p) beta' < 0`, as Barker et al. prove.

### The cap

`Y = alpha(I) p - p^2 / C` vanishes at `p = alpha(I) C`, which increases with
shear rate from `alpha(0) C = 2 mu_1 C` at rest. Without a bound on `p`, `psi`
is not convex in `D` (it decreases with shear rate where `Y < 0`), and pure
compaction has unbounded dissipation: Barker et al.'s flow rule
`div u = 2 f |D|` allows no compaction without shear, so an unsheared column
could not consolidate under its own weight. The flow rule therefore has to be
completed by a bound on `p`, and the bound has to keep `psi` convex.

The bound used is the fixed cap `p <= alpha(0) C = 2 mu_1 C`. Then
`psi(s, delta) = sup_{0 <= p <= 2 mu_1 C} [g(s, p) - p delta]` is a supremum
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
shearing does too, with shear stress `Y(2 mu_1 C, I) = 2 mu_1 C (alpha(I) -
alpha(0))`. Completing the flow rule this way is an addition to the published
model and has to be stated as such; it differs from the shear-dependent bound
only for states with `p > 2 mu_1 C`, which are reached only while shearing and
compacting at once.

### States at fixed packing

With `phi`, hence `C`, fixed (checked numerically with the glass-bead `mu(I)`):

| State | Condition | `C` needed to carry pressure `p` |
| --- | --- | --- |
| Steady isochoric flow | `f = 0`, then `|tau| = mu(I) p` (to 1e-16) | `2 p / beta(I)` |
| Critical state (`I -> 0`) | `f = 0`, `I = 0` | `p / mu_1` |
| Loosest state at rest | on the cap at `s = 0` | `p / (2 mu_1)` |

The critical pressure is half the cap pressure at the same packing, as in
modified Cam-clay. Material on the cap (loose for its load) compacts when
sheared; material denser than critical dilates.

### Scale of `C`

With the form of Fannon et al., `C = (P / b) x`, `x = psi / (1 - psi)` and
`psi = (phi - phi_min) / Delta`, `Delta = phi_max - phi_min`, the packing on
the cap at pressure `p` falls short of `phi_max` by

```text
e(p) = Delta / (1 + x),   x = b p / (2 mu_1 P),
```

and the effective bulk stiffness there is `K = phi dp/dphi ~ phi p / e`
(checked to 0.8%). The scale `P / b` is therefore a choice of compressibility:
it fixes how far below `phi_max` the material sits at a given load, and the
stiffness follows. A rule for choosing it is a packing deficit `e_ref` at a
reference pressure `p_ref`:

```text
P / b = p_ref / (2 mu_1 (Delta / e_ref - 1)).
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

At `D = 0` the potential is `psi(0, delta) = sup_{0 <= p <= 2 mu_1 C} (-p delta)`,
and its subdifferential gives the stress:

- `delta > 0` (separating): `p = 0` and `tau = 0`, stress-free, as in the present
  model's opening;
- `delta < 0` (compacting): `p = 2 mu_1 C` (the cap) and `tau = 0`;
- `delta = 0` (static): any `p` in `[0, 2 mu_1 C]` and any `tau` with
  `|tau| <= 2 mu_1 p - p^2 / C`.

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
| Cap | `g_p(s, p_cap) >= delta` | `p_cap = 2 mu_1 C` | `2 Y(p_cap, I)` | `-p_cap` |
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
- With the fixed cap the saddle form `min_u max_{0 <= p <= 2 mu_1 C}` has a
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
  `0 <= p <= 2 mu_1 C` of the discs of radius `Y_0(p) = 2 mu_1 p - p^2 / C`
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

## References

The papers are listed in the [literature review](connected-flows-literature.md).
