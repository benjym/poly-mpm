# Design note: a well-posed granular step

The rigid-grain solver cannot refine connected flows: with constant friction, a
hard packing limit and no rate dependence the continuum model is ill-posed, and
the discrete problem stalls every solver tried
([literature review](connected-flows-literature.md),
[regularisation study](granular-mixed.md#regularised-laws-under-refinement)).
This note sets out the model, the time step and the solver that would replace
it, what the code already provides and what must be built or decided first. It
is a plan; nothing here is implemented.

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
| Newton solver on the potential with line search | Missing |
| Multigrid hierarchy on the MPM grid | Missing |
| GPU convex solver | Missing; the Metal port solves the stress-space problem this replaces, though its operators, binning, transfers and segregation carry over |
| Shear-cell verification against published results | Inertial model done; critical-state model blocked on jammed cells ([Shear-cell results](#shear-cell-results)) |
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

`Y = alpha(I) p - p^2 / C` vanishes at the cap `p = p_cap(s)`, the root of
`p = alpha(k s / sqrt(p)) C`; beyond it no stress state is admissible. With the
glass-bead `mu(I)` and `C = 4e4` Pa, `p_cap` increases with `s` and is concave,
starting from `p_cap(0) = 2 mu_1 C`.

The dissipation potential was evaluated by maximising `g - p delta` over `p`,
with and without the bound `p <= p_cap(s)`, and tested for midpoint convexity
in `(s, delta)` and monotonicity in `s` on 400 random pairs:

| Bound on `p` | Convex in `(s, delta)` | Nondecreasing in `s` | `psi` under pure compaction (`s = 0`, `delta = -1`) |
| --- | --- | --- | --- |
| `0 <= p <= p_cap(s)` | no violation found | no violation found | `2 mu_1 C` (finite) |
| `p >= 0` | no violation found | falls by up to 61% | bounded only by the search interval |

Without the cap, `psi` is not convex in `D` (it decreases with shear rate where
`Y < 0`), and pure compaction has unbounded dissipation: Barker et al.'s flow
rule `div u = 2 f |D|` allows no compaction without shear, so an unsheared
column could not consolidate under its own weight. With the cap, compaction at
zero shear proceeds at the cap pressure, as normality gives at the cap in
critical-state soil mechanics. Completing the flow rule this way is an addition
to the published model and has to be stated as such. Sampled midpoint tests are
evidence, not a proof of convexity.

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
where `p*` maximises `g - p delta` on `[0, p_cap(s)]`. Because `g` is concave in
`p`, `p*` is unique and found by a bracketed scalar root of `g_p(s, p) = delta`,
with three branches:

| Branch | Condition | `p*` | `psi_s` (= `2 |tau|`) | `psi_delta` (= `-p`) |
| --- | --- | --- | --- | --- |
| Interior | `g_p(s, 0) > delta > g_p(s, p_cap)` | root of `g_p = delta` | `2 Y(p*, I*)` | `-p*` |
| Cap | `g_p(s, p_cap) >= delta` | `p_cap(s)` | `(g_p - delta) p_cap'(s)` | `-p_cap(s)` |
| Zero pressure | `g_p(s, 0) <= delta` | `0` | `0` | `0` |

On the interior branch the Hessian in `(s, delta)` is

```text
psi_ss = g_ss - g_sp^2 / g_pp,   psi_sdelta = g_sp / g_pp,   psi_deltadelta = -1 / g_pp,
```

with determinant `-g_ss / g_pp`. Since `g_ss >= 0` and `g_pp < 0` by conditions
(b) and (c), it is positive definite: `psi` is convex on the interior branch
analytically, not only by sampling. On the cap branch the shear stress
`(g_p - delta) p_cap' / 2` is nonzero although `Y = 0` there; this follows from
the cap bound depending on `s` through `I`, and is a consequence of the
completion at the cap rather than of the published model. Gradients on all three
branches agree with finite differences to 1e-8 or better, and the Hessian to
2e-4 (the remaining difference shrinks with the step, so it is
finite-difference error).

The stress tensor follows from `s = |D|`: `tau_ij = psi_s D_ij / (2 s)`, so that
`|tau| = psi_s / 2`, and `p = -psi_delta`.

### Consequences for the time step

- The step is the convex minimisation of the [Time step](#time-step) section with
  `psi` from the capped potential; `phi`, and hence `C`, is held at its value from
  the start of the step.
- The saddle form `min_u max_p` does not have a product domain, because the
  bound on `p` depends on `|D(u)|`. The primal form in `u` alone is the natural
  one, with `psi` evaluated at each quadrature point by a one-dimensional concave
  maximisation in `p`; no closed form was found.
- `psi` is nonsmooth at zero shear (yield) and at the cap. The solver must handle
  both: a nonsmooth Newton method on the primal problem, or a smoothing with a
  continuation parameter. `alpha(I)` has no conic representation, so the
  interior-point route used in the regularisation study does not carry over
  directly.

### Open points

- The value of the packing deficit `e_ref`, and hence `P / b`: a modelling
  choice, as above. The form of `C` near `phi_max` is speculative in the
  published model.
- Whether holding `phi` fixed within a step preserves well-posedness. This is
  not settled analytically. In the shear cell with the inertial model, lagged and
  implicit `phi` differ by less than the time discretisation error
  ([Shear-cell results](#shear-cell-results)); the critical-state comparison
  awaits a solver that handles jammed cells.
- Convexity on the cap branch and across branch boundaries rests on sampling
  (section [The cap](#the-cap)); the interior branch is convex analytically.

## Shear-cell results

`benchmarks/cidr_shear_cell.py` solves the gravity-free shear cell of Schaeffer
et al. (2019) with the time step of this note: a Lagrangian one-dimensional
grid, `phi` lagged and Newton's method with a backtracking line search on the
potential. Case A is their (3.16) (height 30, `phi = 0.78`, a localised
vertical-velocity perturbation of amplitude 0.05, to `t = 20`); case B is their
(5.7) (height 100, sinusoidal perturbations of amplitude 0.16 in `u` and 0.1 in
`w`, to `t = 1000`). The critical-state model uses a packing deficit of 0.005.
Results are in `benchmarks/cidr_shear_cell.json`.

Grid convergence compares `max |w|` from 46 and 200 cells at 40 times spaced
logarithmically:

| Model | Case | Largest difference | Signal | Newton iterations per step | Final `max |w|` |
| --- | --- | --- | --- | --- | --- |
| Inertial | A | 3.1e-4 | 0.021 | 1.2 / 1.3 | 1.8e-5 |
| Inertial | B | 1.3e-3 | 0.10 | 1.2 / 1.3 | 1.3e-4 |
| Critical-state | A | 7.2e-4 | 0.021 | 1.0 / 1.0 | 7.7e-4 |
| Critical-state | B | Newton fails on 200 cells at `t = 10` | 0.10 | | |

Inertial case A relaxes to the uniform shear by `t = 20`, consistent with
Schaeffer et al.'s report that the flow is steady by about `t = 10`. Critical-state
case A retains a small residual vertical velocity at `t = 20`. For inertial case B,
`phi` lagged against `phi` implicit (fixed-point iteration within the step)
differs by 4.6e-4 in `max |w|`, and lagged at the full step against lagged at
half the step by 1.8e-4, against a signal of 0.10. Lagging `phi` costs no more
than the time discretisation error here.

Critical-state case B compacts a layer to `phi = 0.8163`, within 1e-3 of
`phi_max`, by `t = 8`. Up to that time 46, 100, 200 and 400 cells give the same
packing range to 1e-4. Within that layer cells stop deforming (`|D|` about
1e-14), and `psi` there sits at the cone tip of [Stress at zero shear](#stress-at-zero-shear):
the gradient Newton evaluates is one arbitrary subgradient, and on 200 cells
Newton stalls with an error of 0.08 in that gradient. On 46 cells the run
reaches `t = 1000`. The packing-lag comparison for this model fails the same way.
The inertial model has no static yield, so its `psi` is smooth at `D = 0`.

Two findings carry over to the two-dimensional step:

- With `phi` lagged and `C` growing without bound towards `phi_max`, a step can
  compact a cell past `phi_max`. The harness limits the step so that no cell
  closes more than half its remaining gap `(phi_max - phi) / phi` at its current
  compaction rate. Before this limit was added (together with a fix to the Newton
  stopping test), case B left the admissible range of `phi`.
- The critical-state model needs a solver that treats the cone at `D = 0`
  directly. Plain Newton on the primal problem is not sufficient once a region
  jams.

## References

The papers are listed in the [literature review](connected-flows-literature.md).
