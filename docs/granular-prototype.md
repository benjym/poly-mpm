# Experimental rigid-grain solver

This opt-in branch explores dense, dry, noncohesive granular motion with static
zones and free surfaces. It replaces elastic pressure evolution and capped
apparent viscosity with a coupled unilateral packing and Coulomb friction solve.
The existing explicit solver remains the default.

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
Any existing segregation operator remains separately explicit and uncalibrated
for this constitutive model.

Controls default to `granular_tolerance=1e-9`,
`granular_max_iterations=200`, and `granular_max_nodes=256`.
The iteration control bounds least-squares evaluations and allows up to ten
times that many fallback sweeps. The node cap limits the prototype's dense
linear algebra; raising it does not make the method scalable.
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
momentum residuals, contact dissipation and solver convergence. Packing in this
report is measured at the start of each step. Tests additionally cover analytic
sticking/sliding/opening, coupled complementarity, small timesteps, dilute
freefall, periodic species conservation, solver rejection and saved output.

These are coarse-grid verification experiments. They do not establish spatial
convergence, quantitative column-collapse predictions or a 100x speedup against
an elastic model. The model deliberately changes the physics, so its runtime
is not an equivalent numerical-behavior benchmark against the existing solver.

In the recorded column experiment, maximum speed is approximately `2e-15 m/s`
at `dt=0.02 s`. For collapse at `t=0.3 s`, the rightmost particle is at
`0.56549`, `0.56248` and `0.56084 m` for timesteps `0.01`, `0.005` and
`0.0025 s`, respectively. Maximum observed nodal packing is about `0.621`,
roughly 3.5% above the specified `0.600` limit, and does not disappear in this
timestep refinement. Tight complementarity residuals therefore do not imply
accurate geometric packing. This is an unresolved discretization issue, and
needs investigation before quantitative use.

Before production use, establish mixed-space stability and mesh convergence,
test inclined layers and sustained static/flowing interfaces, and compare
collapse/runout against experiments. A scalable sparse contact algorithm is
needed for larger grids. Constant local Coulomb friction does not establish
continuum well-posedness, size-dependent stopping or creep predictions; a
regularized or nonlocal constitutive extension remains a separate scientific
decision after the baseline passes these checks.
