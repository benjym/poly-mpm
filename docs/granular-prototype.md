# Experimental rigid-grain solver

This opt-in branch explores dense, dry, noncohesive granular motion with static
zones and free surfaces. It replaces elastic pressure evolution and capped
apparent viscosity with a coupled unilateral packing and Coulomb friction solve.
The existing explicit solver remains the default.

The [source-code audit and citing-paper review](granular-source-audit.md) adds
controlled gradient-transfer, plane-strain and viscosity options. Those options
are separate from the baseline described here.

```sh
uv run --locked python MPM.py granular_column
uv run --locked python MPM.py granular_collapse
```

The column is supported by frictionless fixed walls. The collapse uses fixed
no-slip walls. Both examples use linear interpolation, point quadrature and FLIP
transfer, with one material configured as follows:

```python
P.solver = "granular"
P.S.law = "granular"
P.S.rho_s = 2500.0       # intrinsic grain density, kg/m^3
P.S.rho = 1500.0        # initial bulk density, kg/m^3
P.S.packing_max = 0.6
P.S.mu = 0.5            # intrinsic 2D shear/pressure yield ratio
```

These attributes can also be supplied in an input's parameter constructor, as
the examples do. No bulk modulus, shear modulus or viscosity cap is required.
Initial bulk density must not exceed `rho_s * packing_max`.

## Discrete model

The starting point is the nonsmooth formulation in
[Daviet and Bertails-Descoubes (2016), equations 7–14](https://gdaviet.fr/files/mpm/mpmGranular.pdf).
This is an independent, intrinsic **2D** prototype, rather than a reproduction
of their 3D Sand6 implementation. With solid particle volume
`V_grain = mass / rho_s`, the same assembled operator `B` and its transpose
couple nodal velocities and phase contact stresses:

```text
u* = velocity after external force and fixed-wall constraints
u  = u* + dt M^-1 B^T r
w  = B u + k
W  = dt B M^-1 B^T
```

`M` is lumped mass on free velocity degrees of freedom. Stress is expanded in
an orthonormal basis of 2D symmetric tensors, with the first component along
`I / sqrt(2)`. At each occupied stress node, the contact reaction `r` obeys

```text
r_n >= 0, w_n >= 0, r_n*w_n = 0
|r_t| <= mu*r_n
sliding: r_t = -mu*r_n*w_t / |w_t|
```

Normal complementarity is enforced separately from the tangential friction
disk. A single uncorrected cone projection would instead impose associated
dilation. Below yield, zero strain rate and a nonzero admissible stress are
possible without a viscosity cap; open regions carry no contact stress.
The baseline uses constant friction, without a mu(I) or nonlocal closure.

The volumetric gap is
`k_n = max(packing_max * integral(N) - sum(V_grain*N), 0) / (sqrt(2)*dt)`.
Basis volumes are integrated over the computational domain. Differences within
64 machine epsilons of the local capacity/occupation scale are treated as zero
before division by `dt`, to avoid amplifying summation roundoff. This
linearized packing condition does not repair existing geometric overpacking
or enforce an exact end-of-step particle density bound.

The nonlinear solve uses a dimensionless natural-map residual, an analytic
generalized Jacobian, and SciPy least squares with a dense Jacobian factorization.
A nonlinear block Gauss-Seidel contact solve is the fallback if least squares
stalls. Positive feasible packing gaps are excluded from residual scaling so
large free-surface clearances cannot hide small compressive violations.
Failure to meet `granular_tolerance` raises an error; optimizer termination alone
does not count as convergence. Dependent stress modes can make stresses
nonunique even when velocities agree.

## Scope and controls

The current branch requires updated mechanics, linear interpolation, point
quadrature and one granular material. FLIP is the example default. PIC is
available but showed substantial timestep-dependent damping in the collapse
experiment. Quadratic interpolation/APIC and finite-domain quadrature are
excluded pending a suitable mixed discretization and convergence checks.
Fixed frictionless or no-slip walls and horizontal periodicity are supported;
conveyors, traction boundary flags, finite wall friction and damping are rejected.

By default a wall node has its normal velocity (and, for no-slip walls, its
tangential velocity) set to zero in both directions. The mass mapped to wall
nodes then cannot move away from the wall, so grains within one cell of a wall
lag the material that is leaving it: under upward gravity a resting column
stays partly on the floor. `B.separable = True` makes the domain edges
one-sided. A wall node is fixed only while its velocity, after the gathered
momentum and again after gravity, does not move it away from the wall. If the
solve then drives a released node into its wall, the node is fixed again and
the step re-solved, starting from the stresses of the attempt it replaces.
Interior walls stay two-sided. A no-slip wall still fixes the tangential
velocity of a node that is not moving away from it, even without a normal
load.
Any existing segregation operator remains separately explicit and uncalibrated
for this constitutive model.

Controls default to `granular_tolerance=1e-9`,
`granular_solver='hybrid'`, `granular_max_iterations=200`, and `granular_max_nodes=256`.
The iteration control bounds least-squares evaluations and allows up to ten
times that many fallback sweeps. The default hybrid solve retains a dense
Jacobian factorization for reference, with a sparse block fallback.
`granular_solver='block'` skips the global factorization and uses sparse 3x3
block columns for local contact updates; only local diagonal blocks are dense.
This reduces storage, but does not guarantee fast convergence. It is optional
because the inclined-flow acceptance study below fails for the block-only solve.
The node cap remains a guard on research runs; raising it does not establish
scalability.
`granular_solver='velocity'` selects the [compact velocity-based iteration](granular-mixed.md),
with cached local metrics and scalar sliding solves inspired by Daviet (2026).
It supports both stress spaces while retaining the existing discrete equations;
the hybrid reference remains the default.
`granular_stress_basis='P0'` enables the separate cellwise-constant stress
experiment; nodal stresses remain default. Its convergence and flow-accuracy
limits are documented in the same guide.
Elastic-wave and effective-viscosity estimates do not apply to this branch.
Advection, segregation and custom material timestep restrictions still apply.
Timestep accuracy and contact convergence must also be checked.

For compatibility with the existing 3D-shaped output arrays, `sigma_zz = -p`
is supplied as bookkeeping. It is not an independently solved plane-strain
stress. Contact stresses are phase stresses; particle Cauchy stresses are
converted to bulk stresses using `rho/rho_s`. Shear-rate diagnostics use the
intrinsic 2D deviator. No elastic stored energy is introduced. Contact work is
recorded as plastic work; transfer and time discretization can contribute an
energy residual.

## Evidence and next decisions

Run the bounded experiment with plotting disabled:

```sh
uv run --locked python -m benchmarks.granular_prototype --output /tmp/granular-prototype.json
```

The committed [report](../benchmarks/granular_prototype.json) records source
fingerprints, environment, runtime, static drift, timestep sensitivity, mass and
momentum residuals, contact dissipation and solver convergence. Packing is now
measured both at the start and after each step, along with particle density
`rho/rho_s` inferred from deformation. These are distinct measures: the latter
does not affect the fixed solid volumes used for nodal packing. Tests additionally cover analytic
sticking/sliding/opening, coupled complementarity, small timesteps, dilute
freefall, periodic species conservation, solver rejection and saved output.

These are coarse-grid verification experiments. They do not establish spatial
convergence, quantitative column-collapse predictions or a 100x speedup against
an elastic model. The model deliberately changes the physics, so its runtime
is not an equivalent numerical-behavior benchmark against the existing solver.

In the recorded column experiment, maximum speed is below `1e-12 m/s`
at `dt=0.02 s`. For collapse at `t=0.3 s`, the rightmost particle is at
`0.56549`, `0.56248` and `0.56084 m` for timesteps `0.01`, `0.005` and
`0.0025 s`, respectively. Maximum observed nodal packing is about `0.621`,
roughly 3.5% above the specified `0.600` limit, and does not disappear in this
timestep refinement. Tight complementarity residuals therefore do not imply
accurate geometric packing. This is an unresolved discretization issue, and
needs investigation before quantitative use.

## Packing and sparse-solver acceptance study

```sh
uv run --locked python -m benchmarks.granular_packing --output /tmp/granular-packing.json
```

The [recorded study](../benchmarks/granular_packing.json) separates particle
sampling from grid refinement. Each collapse lasts `0.3 s`; the timestep is
halved with grid spacing. Particles per cell below means per coordinate axis.

| Grid spacing | Particles per cell axis | Maximum nodal packing after a step | Maximum deformation-based particle packing | Rightmost particle at 0.3 s |
| ---: | ---: | ---: | ---: | ---: |
| 0.25 | 2 | 0.62122 | 0.80824 | 0.56248 |
| 0.25 | 4 | 0.60661 | 0.87128 | 0.61716 |
| 0.125 | 2 | 0.69347 | 0.89205 | 0.66624 |

Changing particle count also changes the initial rightmost sample location;
centre of mass and displacement are retained in the JSON for comparison.
Neither density reconstruction nor deformation-based density is validated by a
small contact residual. These data fail the spatial-convergence and packing
acceptance gates. More particles improve nodal reconstruction on this fixed
grid while the deformation-based density error grows. A timestep reduction
alone is not a demonstrated fix.

A separate manufactured sampling experiment applies the **exact** map
`diag(exp(0.1), exp(-0.1))` to a uniform particle rectangle. Its determinant is
one, so physical packing is unchanged. Despite exact area preservation and
conserved total solid volume, interior nodal reconstruction has absolute packing
error about `0.0120` at two particles per cell axis. At eight it is about
`0.000845`; this phase-dependent error is not monotonically decreasing at every
sampling level. Thus sampling alias exists independently of time integration.
It is one contributor, rather than a complete explanation of the dynamic drift.

The sparse block solver matches the coarse-grid hybrid collapse velocities
within the regression tolerance. A 256-contact coupled test also checks that
the block path never densifies a global matrix. However, the refined collapse
needs up to 788 sweeps, and the block-only inclined layer above yield fails
the `1e-9` residual tolerance after 2,000 sweeps. In this study it is slower
than hybrid even on the coarse grid. Sparse storage alone is not a performance
solution.

For the hybrid reference, a periodic layer at `tan(theta)=0.25 < mu=0.5`
remains static. At `tan(theta)=0.75`, its maximum downslope speed after
`0.05 s` agrees with the uniform-sliding acceleration reference
`g*(sin(theta)-mu*cos(theta))`, about `0.0981 m/s`. This is a short yield-onset
check, not a calibrated steady-flow profile or static/flowing interface test.

Simply allowing a negative packing gap to force expansion is not adopted. The
report includes an analytic one-mass counterexample: with no external force,
repairing overlap through a velocity impulse creates kinetic energy. A future
geometric correction must explicitly account for its work, rather than hiding
it in plastic dissipation.

## Next implementation gates

1. **Source-paper audit:** first reproduce the paper's choices closely enough to
   separate implementation differences from limitations of the method. The
   present equations are an adaptation, not a verified reproduction. Audit
   transfer, dimensional reduction, viscosity, packing reconstruction and solver
   tolerance before replacing the discretization.
2. **Volume fidelity:** compare a better-integrated, compatible pressure/velocity
   discretization against this baseline. Require both nodal and material-volume
   errors to decrease under independent particle and grid refinement, with
   conserved mass and explicit correction-work accounting. Clamping density
   output does not pass this gate.
3. **Contact solve:** keep the hybrid reference until a sparse method passes
   sticking, sliding, opening, collapse and inclined flow at the same residual
   tolerance. Benchmark work per contact and iteration growth on at least three
   grid levels; matrix storage alone does not pass this gate.
4. **Granular predictions:** once those gates pass, test sustained static/flowing
   interfaces, inclined-layer stopping and experimentally calibrated collapse
   runout. Then select a regularized or nonlocal friction closure for the target
   regime. The present constant-friction baseline cannot establish size effects
   or creep predictions.

### What the source paper does and does not guarantee

The [source paper](https://gdaviet.fr/files/mpm/mpmGranular.pdf), Section 5.2,
explicitly notes that clamping the packing-gap coefficients to preserve
dissipativity can cause volume loss over time. Section 7.5 also identifies
volume loss from enforcing the rheology only at discrete points. Its published
results therefore do not establish an exact geometric density bound.
This observation does not establish that the magnitude of our errors is
acceptable or attributable solely to that limitation.

| Item | Source formulation | Current implementation |
| --- | --- | --- |
| Particle solid volume | Fixed grain volume | Fixed `mass/rho_s`; matches this assumption |
| Spaces and integration | Grid-aligned trilinear velocity/stress bases, particle quadrature | Bilinear 2D counterpart with point quadrature |
| Packing gap | Positive-part coefficients | Same positive-part construction in 2D |
| Particle/grid velocity transfer | APIC, equations 17–19 | FLIP by default; a material difference to audit |
| Stress dimensionality | 3D symmetric tensors | Intrinsic 2D tensors; output `sigma_zz` is bookkeeping |
| Newtonian viscosity | Small nonzero viscosity in the trial momentum solve | Zero viscosity |
| Nonlinear solve | Coulomb solver library, including modified Fischer–Burmeister local solves | Independent natural-map least squares and angular local solves |
| Packing diagnostic | Reconstructed grain volume fraction | Reconstructed nodal packing, plus a separate deformation-based diagnostic |

The deformation-based `rho/rho_s` values in our report are **not** the paper's
packing variable. They derive from the existing MPM compatibility fields
`F/J/V/rho`; neither the paper's packing constraint nor our adapted contact
solve directly constrains them. Their disagreement is a useful consistency
warning, but it is not evidence that the source method predicts those same
values. Source comparisons should primarily use the reconstructed fixed-grain
volume fraction, and independently inspect material-volume consistency.

The next audit should check assembled `B`, mass and gap on a frozen cloud
against independently written equations, then compare a source-style APIC
transfer and a justified dimensional reduction on matched geometry and output
measures. Neither the current regression suite nor the sampling experiment
establishes that the implementation is bug-free. A different continuum closure
should not be justified solely by errors from this unverified adaptation.

Before production use, establish mixed-space stability and mesh convergence,
test inclined layers and sustained static/flowing interfaces, and compare
collapse/runout against experiments. A scalable sparse contact algorithm is
needed for larger grids. Constant local Coulomb friction does not establish
continuum well-posedness, size-dependent stopping or creep predictions; a
regularized or nonlocal constitutive extension remains a separate scientific
decision after the baseline passes these checks.
