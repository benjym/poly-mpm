# Constitutive models

The material registry supports prescribed rigid motion, hyperelastic solids,
rate-form plastic solids, Newtonian/Bingham fluids, and local granular flow.
Select parameters and units using [material selection](materials.rst).
All tensors are 3 by 3, including the out-of-plane normal stress in plane strain.
Cauchy stress is tension positive. Below, `p = -tr(stress)/3` is compression
positive, `s = stress + p I`, `q = sqrt(3/2 s:s)`, and the engineering equivalent
shear rate is `gamma = sqrt(2 D':D')`, with `D' = dev(dstrain)/dt`.
`Particles.p` and `.q` compute these invariants from current stress. The existing
`.pressure` diagnostic remains tension positive; for Bingham it stores the
elastic mean part, excluding instantaneous bulk-viscous stress.

## Model parameters

| Law | Parameters (SI units) | Behaviour |
|---|---|---|
| `rigid` | Optional `initial_v` | Prescribed planar motion; no constitutive integration. |
| `elastic` | `K, G > 0` (Pa) | Compressible neo-Hookean elasticity with reference prestress; see [deformation and stress](numerical-foundations.md#stress-and-work-conventions). |
| `von_mises` | `K, G > 0`, `yield_stress > 0` (Pa) | Perfect J2 plasticity with objective spin integration. |
| `viscous` | `K >= 0` (Pa), `eta, eta_v >= 0` (Pa s) | Newtonian deviatoric stress `2 eta D'`, bulk viscous stress `eta_v tr(D) I`, elastic mean-stress increment `K tr(dstrain)`. |
| `pouliquen` | `K >= 0`, `rho_s, I_0, eta_max > 0`, `mu_s >= 0`, `mu_d >= mu_s` | Bounded friction `mu(I)=mu_s+(mu_d-mu_s) I/(I_0+I)` and capped viscosity `mu p/gamma`; mean grain diameter sets `I`. |
| `bingham` | `K >= 0` (Pa), `eta >= 0`, `eta_0 > eta`, optional `eta_v >= 0` (Pa s), `yield_stress >= 0` (Pa) | Continuous bi-viscosity regularisation with engineering shear invariant and dissipative bulk viscosity. |
| `dp` | `K, G > 0` (Pa), `mu_s > 0`, `0 <= beta <= mu_s` | Cohesionless, perfectly plastic Drucker–Prager solid with nonassociated dilation and a zero mean-tensile-strength cutoff. Uses a cone/apex return map. |
| `dp_rate` | Same as `dp`, plus `t_star >= 0` (s) | Linear Perzyna overstress regularisation, integrated by backward Euler. |
| `linear_mu` | `K >= 0` (Pa), `mu_s, b >= 0`, `rho_s > 0` (kg/m³), `eta_max > 0` (Pa s), grain sizes and particle compositions | Dense flowing `mu(I) = mu_s + b I`, finite viscosity cap, updated compression-positive confinement, explicit zero-pressure/rate limits. |

Initialization rejects missing, nonfinite, nonscalar, and inadmissible parameters
for these laws. `beta` is a plastic-potential slope, not an angle. `mu_s` in DP is
the **q/p slope**, not `tan(friction angle)` or the `tau/p` friction coefficient
used by granular flow laws. No Mohr–Coulomb fitting convention is assumed.

### Bingham

The shear stress is `s = 2 eta_eff D'`, where

```
eta_eff = min(eta_0, eta + yield_stress/gamma)
gamma_c = yield_stress/(eta_0-eta)
```

`eta` is the plastic viscosity and `eta_0` the unyielded viscosity; the particle
diagnostic `eta` records `eta_eff`. At zero rate, `eta_eff = eta_0`. At
`yield_stress = 0`, the law reduces to a Newtonian fluid with viscosity `eta`,
including at zero rate; `eta_0 >= eta` then suffices.
An optional `gamma_c` input must agree with the derived transition.
For simple shear the stress is continuous and becomes
`abs(stress_xy) = yield_stress + eta*abs(du/dy)` above that transition.
This finite-viscosity approximation permits small creep below yield; it does not
represent an exactly rigid plug. The engineering shear invariant and the
underlying Bingham relation are also described in the
[FEBio fluid material documentation](https://febiosoftware.github.io/febio-docs/features/features/fluid_materialprop_bingham/).

The elastic mean stress increment is `K tr(dstrain)`, and the instantaneous bulk
viscous stress is `eta_v tr(D) I`. The bulk-viscous contribution dissipates `eta_v tr(D)^2` per unit current volume.

### Drucker–Prager and its rate regularisation

The yield function and plastic potential are

```
f = q - mu_s*p
g = q - beta*p
```

For the elastic trial state, the plastic multiplier on the cone is

```
dlambda = max(f_trial, 0)/(3*G + K*mu_s*beta + eta_p/dt)
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

Both solids use a symmetric spin split: rotate the incoming stress by half the incremental
spin, perform the constitutive update, then rotate by the other half. This is a
corotational hypoelastic integration, exact for pure rigid spin with constant
angular velocity.
It is not a finite-strain hyperelastic or multiplicative plasticity model.

### Linear friction

```
I = gamma * (phi @ grain_sizes) * sqrt(rho_s/p)
eta = min(eta_max, mu_s*p/gamma + b*(phi @ grain_sizes)*sqrt(rho_s*p))
s = 2 eta D'
```

The sum for `eta` is evaluated directly so `p=0` and `gamma=0` do not cause
singular arithmetic. At `p=0`, stress and effective viscosity are zero and the
undefined `I` diagnostic is recorded as zero (`mu` then records `mu_s`).
At zero rate and positive confinement the effective viscosity is capped, unless
`mu_s=0`, where the finite Newtonian limit of the second term is retained.
`eta` and `eta_limited` both report the actual effective viscosity.
Pressure evolves by `K tr(dstrain)` and is clamped to compression before computing
friction. There is no arbitrary density threshold. Setting `b=0` gives a capped
constant-friction flowing law.

This is a dense-flow approximation. Its friction grows without bound at high
`I`; it is not a dilute/collisional model or a claim of well-posedness over every
flow regime. Local granular rheology has known continuum limitations; see
[Barker et al.'s stability analysis](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/wellposed-and-illposed-behaviour-of-the-itmui-rheology-for-granular-flow/BABE7C7E3E783B9B71D7D12888F4F1D4).
Bingham and linear-friction fluids prescribe current stress directly, so the
spin correction is not applied to their algebraic viscous stress.

## Model limits and unsupported names

The plastic laws use corotational hypoelasticity, not multiplicative finite-
strain plasticity. Their quadratic stored-energy diagnostic is an estimate.
Granular flow laws describe dense flow and need calibration for the material and
loading regime. Viscosity caps regularise the calculation and permit creep;
they do not impose an exactly rigid plug or prove continuum well-posedness.

`elastic_time`, `ken_kamrin`, and `ken_simple` are unsupported registry names.
Use `elastic` for constant stiffness, `pouliquen` for the supported bounded local
friction law, or `linear_mu` with `b=0` for capped constant friction. These are
model choices requiring appropriate parameters, not automatic substitutions.
`viscous_size` and the separate `P.F` compressible-fluid path are unimplemented.

## Checks and outputs

`tests/test_constitutive.py` covers analytical shear/volume responses, transition
continuity, flow reversal, dissipation, unloading, plastic flow direction, yield
consistency, cone apex, zero confinement, Perzyna relaxation and timestep
convergence, tensor rotation covariance, rigid-spin integration, empty batches,
and parameter rejection. The [verification guide](verification.md) covers
full-solver equilibrium, translation and finite-strain manufactured solutions.
Implementation verification does not replace experimental calibration.

The DP loading inputs save `data/pq_history.npz` and `pq_history.png` through the
standard output callbacks. Tensor stresses use tension-positive signs; plotted
`p` is compression positive. See [outputs](outputs.rst) for file conventions.
