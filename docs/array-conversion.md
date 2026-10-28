# NumPy array conversion

## Design and public interfaces

Each material has a writable float64 structure of arrays, accessed through
`L.materials[p]`. Materials share one grid. The coordinator deposits every
material before combining periodic endpoints, identifying occupied nodes, and
normalising mass-weighted fields.

| State | Shape |
|---|---|
| `x`, `v`, `x0` | `(n, 3)` |
| `m`, `rho`, `initial_rho`, `V`, `pressure`, scalar history | `(n,)` |
| `stress`, `strain`, `dstrain`, `dstress`, tensor history | `(n, 3, 3)` |
| `phi` | `(n, ns)` |
| `n_star` | `(n,)`, integer |
| `n` | `(n, 4)`, integer |
| `N` | `(n, 4)` |
| `G` (shape-function gradients) | `(n, 4, 3)` |
| `L` (velocity gradient), `L_T` (transpose view) | `(n, 3, 3)` |

The four neighbours are bottom-left, bottom-right, top-right, top-left. Particle
and grid vectors retain three components; geometry and transport remain two
dimensional. Constitutive tensor traces and norms reduce tensor axes only.

`Particles(P, G, material, material_id)` constructs one batch.
`ParticleSystem(P, G)` coordinates the material list. Supported laws implement
`law(particles, material_parameters, simulation_parameters)` and update history and
`dstress`. The coordinator applies the existing Lie stress-rate correction once;
rigid materials bypass the stress update and move at their prescribed velocities.
The registry in `constit.py` is the only material dispatch mechanism.

Uniform mass/volume properties, velocities, and grain fractions are copied into
independent writable arrays. Initial velocity must be `(3,)` or `(n, 3)`;
composition must be `(ns,)` or `(n, ns)`. Initial fractions must be nonnegative and
sum to one. `PHI`, when provided, takes precedence over `phi`. Coordinates must
have exactly `n` entries. Density, volume, and grain sizes must be positive and
finite. A zero-particle material is valid.

The grid stores `pressure`, `gammadot`, `mu`, `I`, `eta`, and other constitutive
diagnostics as **mass-weighted sums**. Call `G.physical(name)` to divide once by
nodal mass. `G.phi`, `G.s_bar`, `G.pk`, and `G.shear_rate` are already physical
fields. `G.shear_rate` is the bulk velocity-gradient measure formerly written over
the deposited `gammadot` field. Empty-node quotients contribute zero.

Periodic left/right grid endpoints duplicate one physical node. Sum unique
columns (`array.reshape(ny, nx, ...)[:, :-1]`) when measuring periodic totals.
Merge deposited sums once per stage; do not add periodic endpoints a second time.
Particle positions wrap with modulo arithmetic before cell lookup.

## Timestep and compatibility

The supported timestep keeps the existing strain-before-force order:

1. Reset grid, update loading, determine cells and basis functions.
2. Deposit mass, momentum, composition, and kinetic pressure across all materials.
3. Merge periodic sums, calculate valid nodes, and impose momentum constraints.
4. Gather velocity, update strain/history/stress, and deposit diagnostics.
5. Deposit forces, impose force boundaries, and update grid momentum.
6. Calculate bulk shear rate and, if enabled, segregation increments.
7. Update composition and particle position/velocity; handle particle boundaries.

Mass and volume remain fixed during deformable motion. Reconstructed density is
a diagnostic retained from the scalar solver, including its use in the existing
Pouliquen low-density cutoff. This conversion introduces no density evolution law.

The following input conventions remain supported:

- A single material object or list in `P.S`.
- `segregate_grid` as an alias for `segregate`; contradictory explicit values fail.
- Uniform `initial_v`, or per-particle initial velocities.
- Existing pressure initialisation modes, steady viscous/Pouliquen flow, and the
  existing steady Poiseuille formula within its steady regime.
- Boundary flags, default output flags, output filenames, and particle-data column
  order (`x`, `y`, grain fractions).
- Legacy inlet placement and cadence. Only materials declaring `inlet_rate` enter
  through an enabled inlet. Injection properties must be uniform per material;
  each new batch gets freshly initialised constitutive history.

Outlets remove all persistent fields together. Inlets append full batches. Scratch
connectivity is rebuilt after either operation, and material counts follow the
actual array lengths. Prescribed per-particle rigid velocities travel with their
rows through removal.

Callbacks are `at_start(L, P, G)`, `after_every_timestep(L, P, G)`,
`after_every_savetime(L, P, G)`, and `at_end(L, P, G)`. The old
`after_ever_timestep` spelling remains an alias. A custom callback replaces the
corresponding default stage; it is responsible for any desired plotting, energy
recording, or saving at that stage. Default energy histories grow with actual
steps, so dynamic timesteps do not overflow preallocated buffers.

`run` lands on the requested final time without the legacy extra timestep.
Dynamic timestep selection considers all materials and ignores inactive zero-speed
or zero-segregation-rate limits. Invalid resulting timesteps fail explicitly.
`run(..., output=False)` suppresses every output callback.

## Intentional numerical corrections

- Repeated node indices accumulate with `bincount`; advanced-indexed `+=` no
  longer drops particle contributions.
- Velocity gradient and transpose are no longer the same incremented buffer.
  Analytical rotation and nonsymmetric affine-flow tests cover this inherited bug.
- Internal forces retain the spatial component in the stress-gradient contraction.
- Cell lookup uses floor, validates the domain before indexing, and handles exact
  upper edges using the final cell. Periodic wrapping handles multiple crossings.
- Grid validity is computed after periodic merging and all material deposits.
- Segregation gathers use particle connectivity rather than particle masks as grid
  indices. Padding follows spatial rows and derivatives exclude the grain axis.
- Closed-boundary/periodic KT ghost handling uses consistent spatial orientation;
  periodic endpoint duplicates are not treated as extra cells.
- Shear rate, weighted diagnostics, and physical fields are distinguished so no
  field is overwritten or divided by nodal mass twice.
- Friction uses actual boundary masks without modifying them. Rectangular-box
  indexing and the bottom silo mask no longer affect unintended rows.
- At zero pressure, inertial number is undefined. Its stored diagnostic is zero;
  friction uses the infinite-inertial-number limit for nonzero strain rate.
  At zero strain rate with nonzero pressure, viscosity uses the existing cap.
  This avoids overflowing nodal sums with a machine-maximum diagnostic sentinel.
- At exactly zero inertial number, the existing nonfinite-diffusion/zero-flux limit
  is explicit. Near-zero inertial number can still make the retained diffusion
  formula sensitive; it has not been regularised with a new physical parameter.
- Gravity potential energy uses its Cartesian projection, avoiding division by a
  particle's zero x-coordinate in the old angle-based expression.

The constitutive equations, coefficients, tension convention, viscosity cap,
stress-rate formulation, and KT flux/reconstruction equations are otherwise
retained. Changes are checked against analytical cases and independent corrected
scalar kernels, not blind parity with known defects in `master`.

## Validation and remaining limitations

Tests cover shape functions, affine flow, pure rotation, repeated-node transfers,
unequal particle properties, body/internal forces, constitutive histories, pressure
and flow initialisation, multiple materials, empty batches, domain boundaries,
periodic seams, inlets/outlets, dynamic timestep selection, callbacks, saved array
shapes, and readable nonempty figures. The complete ball run and 100-step
bidisperse/three-class/eight-class segregation runs are included.

Kernel reference comparisons use `rtol=1e-10`, `atol=1e-12`; short trajectories use
`rtol=1e-8`, `atol=1e-12`. Tests run with runtime warnings treated as errors.

Closed-flux KT stages and particle deposition are tested independently for
conservation. The complete segregation/remapping scheme is **not** exactly
composition- or species-mass-conserving. The reproducible report in
`benchmarks/segregation_drift.json` records this without clipping. For the 100-step
slope fixtures, the maximum particle fraction-sum error is approximately
0.00016, 0.00035, and 0.00067 for two, three, and eight classes. Those results are
numerical diagnostics, not validation of the segregation physics. The explicit
legacy `normalise_phi=True` option clips and renormalises; the default is false.

Unsupported experimental laws remain archived outside the active registry. Mixed
materials share the existing momentum grid; this work adds neither contact
resolution nor a separate velocity field per material. Legacy traction scaling
still uses the first material's uniform reference density.

Generated historical HTML under `docs/build` is not regenerated. This document
and the root README describe the current solver.

## Example compatibility

The table below is generated from `python -m benchmarks.audit_examples`: five
bounded timesteps, initial `dt <= 1e-6`, output disabled. Dynamic timestep rules
remain enabled where configured. Original example durations/settings are retained
in source. Unsupported laws are deliberately not substituted with different laws.

| Example command | Status |
|---|---|
| `python MPM.py anisotropy` | supported smoke-tested |
| `python MPM.py avalanche` | supported smoke-tested |
| `python MPM.py ball` | supported smoke-tested |
| `python MPM.py bi_conveyor` | supported smoke-tested |
| `python MPM.py bi_dam_break 11 mix 0 2` | supported smoke-tested |
| `python MPM.py bi_double_insets 11` | supported smoke-tested |
| `python MPM.py bi_drum 11 2` | supported smoke-tested |
| `python MPM.py bi_seg_test 22 2 11` | supported smoke-tested |
| `python MPM.py bi_square_drum 11 2` | supported smoke-tested |
| `python MPM.py bi_vertical_chute 2 11` | supported smoke-tested |
| `python MPM.py bingham_unit_test` | unsupported material law: Unsupported material law 'bingham' |
| `python MPM.py brazillian` | supported smoke-tested |
| `python MPM.py chuteflow` | supported smoke-tested |
| `python MPM.py dam_break` | unsupported material law: Unsupported material law 'dp_rate' |
| `python MPM.py dp_rate_unit_test` | unsupported material law: Unsupported material law 'dp_rate' |
| `python MPM.py dp_unit_test` | unsupported material law: Unsupported material law 'dp' |
| `python MPM.py flow` | supported smoke-tested |
| `python MPM.py fluid` | unsupported material law: Unsupported material law 'compressible' in legacy P.F; use supported P.S materials |
| `python MPM.py francois 11` | supported smoke-tested |
| `python MPM.py gsd` | supported smoke-tested |
| `python MPM.py incline` | supported smoke-tested |
| `python MPM.py macro` | supported smoke-tested |
| `python MPM.py pierre` | unsupported material law: Unsupported material law 'bingham' |
| `python MPM.py poly_avalanche` | unsupported material law: Unsupported material law 'viscous_size' |
| `python MPM.py poly_seg_test 22` | supported smoke-tested |
| `python MPM.py pouliquen_unit_test` | supported smoke-tested |
| `python MPM.py seg_unit_test` | supported smoke-tested |
| `python MPM.py silo 11 2` | supported smoke-tested |
| `python MPM.py slope_stab` | unsupported material law: Unsupported material law 'dp' |
| `python MPM.py stokes` | supported smoke-tested |
| `python MPM.py tri_seg_test` | unsupported material law: Unsupported material law 'viscous_size' |
| `python MPM.py triaxial` | supported smoke-tested |
| `python MPM.py tunnel` | unsupported material law: Unsupported material law 'elastic_time' |
| `python MPM.py vicky` | supported smoke-tested |
| `python MPM.py vicky_latest` | unsupported material law: Unsupported material law 'bingham' |
| `python MPM.py viscous_unit_test` | supported smoke-tested |
| `python MPM.py vm_unit_test` | supported smoke-tested |
