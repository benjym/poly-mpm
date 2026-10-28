# Constitutive model audit

The active registry now includes `bingham`, `dp`, `dp_rate`, and `linear_mu`, as
well as the five previously converted laws. These four kernels are array based
and checked against analytical responses and physical invariants. The old scalar
code is historical evidence, not the numerical reference for these models.

All tensors are 3 by 3, including the out-of-plane normal stress in plane strain.
Cauchy stress is tension positive. Below, `p = -tr(stress)/3` is compression
positive, `s = stress + p I`, `q = sqrt(3/2 s:s)`, and the engineering equivalent
shear rate is `gamma = sqrt(2 D':D')`, with `D' = dev(dstrain)/dt`.
`Particles.p` and `.q` compute these invariants from current stress. The existing
`.pressure` diagnostic remains tension positive; for Bingham it stores the
elastic mean part, excluding instantaneous bulk-viscous stress.

## Restored laws

| Law | Parameters (SI units) | Behaviour and corrections |
|---|---|---|
| `bingham` | `K >= 0` (Pa), `mu_s >= 0`, `mu_0 > mu_s`, optional `mu_v >= 0` (Pa s), `tau_0 >= 0` (Pa) | Continuous bi-viscosity regularisation, correct stress-increment interface, engineering shear invariant, dissipative bulk viscosity. |
| `dp` | `K, G > 0` (Pa), `mu > 0`, `0 <= beta <= mu` | Cohesionless, perfectly plastic Drucker–Prager solid with nonassociated dilation and a zero mean-tensile-strength cutoff. Replaces the singular legacy H²/loading-direction formula with a return map. |
| `dp_rate` | Same as `dp`, plus `t_star >= 0` (s) | Linear Perzyna overstress regularisation, integrated by backward Euler. Replaces a rate update missing its timestep factor. |
| `linear_mu` | `K >= 0` (Pa), `mu_0, b >= 0`, `rho_s > 0` (kg/m³), `eta_max > 0` (Pa s), grain sizes and particle compositions | Dense flowing `mu(I) = mu_0 + b I`, finite viscosity cap, updated compression-positive confinement, explicit zero-pressure/rate limits. Removes singular divisions and ad hoc shear-invariant factors. |

Initialization rejects missing, nonfinite, nonscalar, and inadmissible parameters
for these laws. `beta` is a plastic-potential slope, not an angle. `mu` in DP is
the **q/p slope**, not `tan(friction angle)` or the `tau/p` friction coefficient
used by granular flow laws. No Mohr–Coulomb fitting convention is assumed.

### Bingham

The shear stress is `s = 2 eta D'`, where

```
eta = min(mu_0, mu_s + tau_0/gamma)
gamma_c = tau_0/(mu_0-mu_s)
```

At zero rate, `eta = mu_0`. At `tau_0 = 0`, the law reduces to a Newtonian fluid
with viscosity `mu_s`, including at zero rate; `mu_0 >= mu_s` then suffices.
An optional old `gamma_c` input must agree with the derived transition.
For simple shear the stress is continuous and becomes
`abs(stress_xy) = tau_0 + mu_s*abs(du/dy)` above that transition.
This finite-viscosity approximation permits small creep below yield; it does not
represent an exactly rigid plug. The engineering shear invariant and the
underlying Bingham relation are also described in the
[FEBio fluid material documentation](https://febiosoftware.github.io/febio-docs/features/features/fluid_materialprop_bingham/).

The elastic mean stress increment is `K tr(dstrain)`, and the instantaneous bulk
viscous stress is `mu_v tr(D) I`. The bulk-viscous sign ensures nonnegative
viscous dissipation; the old scalar implementation had the opposite sign and
used a mean strain increment as though it were a volumetric one.

### Drucker–Prager and its rate regularisation

The yield function and plastic potential are

```
f = q - mu*p
g = q - beta*p
```

For the elastic trial state, the plastic multiplier on the cone is

```
dlambda = max(f_trial, 0)/(3*G + K*mu*beta + eta_p/dt)
q_new = q_trial - 3*G*dlambda
p_new = p_trial + K*beta*dlambda
```

`eta_p = 0` for `dp`, and `eta_p = 3*G*t_star` for `dp_rate`.
Thus the rate model satisfies `dlambda/dt = max(f_new, 0)/eta_p` on an active
surface. `t_star` is the shear-overstress relaxation time for zero dilation.
`t_star = 0` gives the rate-independent update exactly. Loading below yield and
elastic unloading do not produce plastic flow. Plastic work on the cone is
`(q - beta*p)*dlambda >= 0` within the admitted parameter range.
The use of pressure-dependent yield/flow potentials and a Perzyna flow rule
follows standard viscoplasticity; see the
[Sandia LAMÉ theory](https://www.sandia.gov/files/sierra/SM_LAME_5_22/main/models/viscoplastic_viscoscram.html).
This implementation is the simple linear, perfectly plastic case, without that
manual's additional damage, hardening, or viscoelastic mechanisms.

The mean-tension cutoff is solved together with shear return: if the cone-only
candidate has `p < 0`, impose `p = 0` and use
`dlambda = q_trial/(3*G + eta_p/dt)`. The rate-independent case reaches the cone
apex (`p=q=0`) without division by either invariant. The rate model can carry
finite, dissipative viscous shear overstress at zero confinement; it cannot carry
positive mean tension. This is a mean-stress cutoff, not a maximum-principal-
stress criterion. The cone return, tension surface, and apex require separate
handling, as explained in
[Itasca's Drucker–Prager formulation](https://docs.itascacg.com/itasca930/common/models/drucker/doc/modeldrucker.html)
(which uses a different shear-invariant normalisation).

Both solids use a symmetric spin split: rotate old stress by half the incremental
spin, perform the constitutive update, then rotate by the other half. This is a
corotational hypoelastic integration, exact for pure rigid spin with constant
angular velocity. The old solver's Lie correction is not added to this result.
It is not a finite-strain hyperelastic or multiplicative plasticity model.

### Linear friction

```
I = gamma * (phi @ grain_sizes) * sqrt(rho_s/p)
eta = min(eta_max, mu_0*p/gamma + b*(phi @ grain_sizes)*sqrt(rho_s*p))
s = 2 eta D'
```

The sum for `eta` is evaluated directly so `p=0` and `gamma=0` do not cause
singular arithmetic. At `p=0`, stress and effective viscosity are zero and the
undefined `I` diagnostic is recorded as zero (`mu` then records `mu_0`).
At zero rate and positive confinement the effective viscosity is capped, unless
`mu_0=0`, where the finite Newtonian limit of the second term is retained.
`eta` and `eta_limited` both report the actual effective viscosity.
Pressure evolves by `K tr(dstrain)` and is clamped to compression before computing
friction. There is no arbitrary density threshold. Setting `b=0` gives a capped
constant-friction flowing law.

This is a dense-flow approximation. Its friction grows without bound at high
`I`; it is not a dilute/collisional model or a claim of well-posedness over every
flow regime. Local granular rheology has known continuum limitations; see
[Barker et al.'s stability analysis](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/wellposed-and-illposed-behaviour-of-the-itmui-rheology-for-granular-flow/BABE7C7E3E783B9B71D7D12888F4F1D4).
Bingham and linear-friction fluids prescribe current stress directly, so the
legacy Lie correction is not applied to their algebraic viscous stress.

## Retired implementations

| Old name | Audit finding | Replacement or decision |
|---|---|---|
| `elastic_time` | The active code uses constant initial modulus; the intended softening is commented out. An undefined phase index guards a 5% particle-mass reduction per timestep, unrelated to strain or elapsed time. | Retire. Use `elastic` for constant stiffness. Excavation/mass removal or modulus evolution requires a separately specified model; `tunnel` remains explicitly unsupported. |
| `ken_kamrin` | Unfinished duplicate local granular law: signed compressive pressure can zero its viscosity, while a gravity-based cap can make viscosity negative. Fixed grain diameter and arbitrary rate floor also bypass the intended state. | Retire this implementation; use the existing `pouliquen` law for the supported local friction rheology. This is not a rejection of Kamrin's published models. |
| `ken_simple` | Uses stale shear rate and switches off stress below a hardcoded `rho = 2500 kg/m³`, so changing material density changes the constitutive regime arbitrarily. | Retire. `linear_mu` with `b=0` supplies a tested constant-friction flowing alternative without that cutoff. |

These names fail early with specific retirement reasons rather than silently
mapping to another law. Their scalar source remains untouched in `legacy/` for
provenance and is excluded from the package. `viscous_size` has no implementation;
the separate legacy `P.F` compressible-fluid path is also unimplemented. Neither
was classified as physically invalid simply because it is missing.

## Migration and checks

These are physical corrections and replacements, not numerical parity ports.
Old coefficients may need recalibration: DP now consistently uses 3D `q=sqrt(3 J2)`
and constant `mu`; `dp_rate.t_star` has the explicit relaxation-time meaning above.
The old DP example's unused `mu(I)` parameters were replaced by `mu=0.5`.
Flowing laws use engineering shear rate, standard bulk modulus, and actual
viscosity units. Their former ad hoc invariant factors are not carried across.

The DP loading examples now record growing array histories through the standard
callbacks and save `data/pq_history.npz` plus a labelled `pq_history.png`.
They no longer use `L.S[p][i]` or allocate histories from nominal timestep counts.

`tests/test_constitutive.py` checks analytical shear/volume responses, transition
continuity, flow reversal, dissipation, hydrostatic elasticity, elastic unloading,
plastic flow direction, yield consistency, cone apex, zero confinement, Perzyna
relaxation and timestep convergence, tensor rotation covariance, rigid-spin
integration, empty batches, parameter rejection, and output histories. Existing
particle/grid/segregation regressions remain in the full suite.

The example audit now admits 33 of 37 examples for bounded smoke runs, including
all seven formerly blocked Bingham/DP examples. The full `dp_unit_test` (50,000
steps, `dt=2e-5`, 1 s) and `dp_rate_unit_test` (5,000 steps, `dt=1e-3`, 5 s) also
completed at their original timesteps with finite stress. The rate-independent
case ended on `q=0.5 p`; the rate-dependent case retained the permitted positive
overstress. These checks establish equation and implementation consistency; they do not replace experimental calibration,
spatial convergence studies, or long runs at the old example timesteps.

The five previously supported laws retain their existing equations and stress
integration. Their scalar-parity tests remain compatibility checks. In particular,
legacy `viscous` bulk normalisation/sign and `pouliquen` invariant/cutoff conventions
have not been brought into alignment with the newly restored laws in this audit.
