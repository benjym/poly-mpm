# Particle arrays and examples

## Design and public interfaces

Each material has a writable float64 structure of arrays, accessed through
`L.materials[p]`. Materials share one grid. The coordinator deposits every
material before combining periodic endpoints, identifying occupied nodes, and
normalising mass-weighted fields.

| State | Shape |
|---|---|
| `x`, `v`, `x0` | `(n, 3)` |
| `m`, `rho`, `initial_rho`, `V`, `pressure`, scalar history | `(n,)` |
| `stress`, `strain`, `dstrain`, `dstress`, `F`, affine/material history | `(n, 3, 3)` |
| `phi` | `(n, ns)` |
| `domain_halfsize` (reference geometry) | `(n, 2)` |
| `n_star` | `(n,)`, integer |
| `n` | `(n, k)`, integer; k=4 linear / k=9 quadratic for point quadrature; domain support varies |
| `N` | `(n, k)` |
| `G` (shape-function gradients), `offsets` | `(n, k, 3)` |
| `L` (velocity gradient), `L_T` (transpose view) | `(n, 3, 3)` |

The reference `domain_halfsize` is zero in point mode. Domain mode validates
material `domain_size` against represented area; see [quadrature](quadrature.md).

Linear neighbours are bottom-left, bottom-right, top-right, top-left. Quadratic
neighbours are a tensor-product nine-entry stencil, with boundary/periodic mapping
described in [numerical methods](numerical-foundations.md). Particle
and grid vectors retain three components; geometry and transport remain two
dimensional. Constitutive tensor traces and norms reduce tensor axes only.

`Particles(P, G, material, material_id)` constructs one batch.
`ParticleSystem(P, G)` coordinates the material list. Supported laws implement
`law(particles, material_parameters, simulation_parameters)` and update history and
`dstress`. The solver uses neo-Hookean elastic stress, objective rate-form
plasticity, and algebraic fluid stress; it evolves `F`, `J`, `V`, and `rho`.
Rigid materials
bypass deformation/stress evolution and follow prescribed velocity. See the
[numerical methods](numerical-foundations.md) and [constitutive models](constitutive-models.md).
The registry in `constit.py` is the only material dispatch mechanism.

Uniform mass/volume properties, velocities, and grain fractions are copied into
independent writable arrays. Initial velocity must be `(3,)` or `(n, 3)`;
composition must be `(ns,)` or `(n, ns)`. Initial fractions must be nonnegative and
sum to one. Coordinates must
have exactly `n` entries. Density, volume, and grain sizes must be positive and
finite. A zero-particle material is valid.

The grid stores `pressure`, `gammadot`, `mu`, `I`, `eta`, and other constitutive
diagnostics as **mass-weighted sums**. Call `G.physical(name)` to divide once by
nodal mass. `G.phi`, `G.s_bar`, `G.pk`, and `G.shear_rate` are already physical
fields. `G.shear_rate` is the bulk velocity-gradient measure; deposited
`gammadot` is a material constitutive diagnostic. Empty-node quotients contribute zero.
`G.segregation_m`, `G.segregation_phim`, and `G.segregation_phi` hold the deformable
transport inventory separately from rigid material. `G.segregation_particle_dphi`
stages the bounded species-aware increments for the next composition update.

Periodic left/right grid endpoints duplicate one physical node. Sum unique
columns (`array.reshape(ny, nx, ...)[:, :-1]`) when measuring periodic totals.
Merge deposited sums once per stage; do not add periodic endpoints a second time.
Particle positions wrap with modulo arithmetic before cell lookup.

## Timestep and input interfaces

The timestep uses a strain-before-force order:

1. Reset grid, update loading, determine cells and basis functions.
2. Deposit mass, momentum, composition, and kinetic pressure across all materials.
3. Merge periodic sums, calculate valid nodes, and impose momentum constraints.
4. Gather velocity, update deformation/volume/density and strain/history/stress,
   and deposit diagnostics.
5. Deposit forces, impose force boundaries, and update grid momentum.
6. Calculate bulk shear rate and, if enabled, conservative segregation exchanges.
7. Return species-aware composition exchanges, then update particle position and
   velocity and handle particle boundaries.

Particle mass remains fixed during deformable motion; updated volume is `J*V0`
and physical density is `m/V`.

Input conventions:

- Every name set on `P`, `P.G`, `P.B`, `P.O` or a material must be a setting
  listed in the [feature reference](features.md); `settings.normalise` rejects
  any other name before the run and fills in defaults.
- A single material object or list in `P.S`.
- Uniform `initial_v`, or per-particle initial velocities.
- `initial_stress` modes (`lithostatic`, `compression`), steady viscous and
  Pouliquen flow (`initial_flow = "steady"`), and the steady Poiseuille formula
  (`initial_flow = "steady_poiseuille"`) within its steady regime.
- `B.walls`, `B.wall_contact` and the other boundary settings, `O.figures` and
  `O.fields`, output filenames and particle-data column order (`x`, `y`, grain
  fractions).
- Inlet placement and cadence use the material's `inlet_rate`, `inlet_gap` and
  `inlet_y`. Only materials declaring `inlet_rate` enter through an inlet listed
  in `B.inlets`. Injection properties must be uniform per material;
  each new batch gets freshly initialised constitutive history.

Outlets remove all persistent fields together. Inlets append full batches. Scratch
connectivity is rebuilt after either operation, and material counts follow the
actual array lengths. Prescribed per-particle rigid velocities travel with their
rows through removal.

Callbacks are `at_start(L, P, G)`, `after_every_timestep(L, P, G)`,
`after_every_savetime(L, P, G)` and `at_end(L, P, G)`. A custom callback replaces the
corresponding default stage; it is responsible for any desired plotting, energy
recording, or saving at that stage. Default energy histories grow with actual
steps, so dynamic timesteps do not overflow preallocated buffers.

`run` stops at the requested final time, shortening the last step when needed.
Dynamic timestep selection considers all materials and ignores inactive zero-speed
or zero-segregation-rate limits. Invalid resulting timesteps fail explicitly.
`run(..., output=False)` suppresses every output callback.

## Numerical conventions

Repeated node indices accumulate all particle contributions. Cell lookup uses
floor, validates positions, and assigns exact upper edges to the final cell.
Periodic coordinates wrap before lookup; node validity is determined after all
materials deposit and periodic sums are merged. Internal forces contract the
Cauchy stress with spatial shape gradients and current particle volume.

Pressure-dependent flow laws use compression-positive confinement. At zero
pressure, the undefined inertial-number diagnostic is stored as zero; Pouliquen
friction uses its high-inertial-number limit when the shear rate is nonzero.
At zero shear rate, effective viscosity follows the selected model's cap or
Newtonian limit. See [constitutive models](constitutive-models.md).

Segregation conserves mobile species inventories using local equal-mass
exchanges. `G.gammadot` and the constitutive fields remain distinct from the
bulk-grid shear rate. The diffusion coefficient has an explicit zero-inertial-
number convention; very small nonzero inertial numbers can still make it
sensitive. See [segregation transport](segregation-conservation.md).

## Runtime and dependency baseline

The supported baseline is Python 3.14, pinned to 3.14.10 in `.python-version`.
Dependency floors are NumPy 2.5.4, SciPy 1.18.1, Astropy 8.0.2, Matplotlib 3.11.2,
and pytest 9.2. The complete resolved environment is recorded in `uv.lock`.
Use `uv sync --extra dev --locked` with a current uv release to reproduce it.
The regression suite runs on Python 3.14. Benchmark reports record their own
environment and source fingerprints; use measurements matching the configuration
you intend to run.

## Validation and remaining limitations

Tests cover shape functions, affine flow, pure rotation, repeated-node transfers,
unequal particle properties, body/internal forces, constitutive histories, pressure
and flow initialisation, multiple materials, empty batches, domain boundaries,
periodic seams, inlets/outlets, dynamic timestep selection, callbacks, saved array
shapes, and readable nonempty figures. The complete ball run and 100-/1,000-step
bidisperse/three-class/eight-class segregation runs are included.

Kernel reference comparisons use `rtol=1e-10`, `atol=1e-12`; short trajectories use
`rtol=1e-8`, `atol=1e-12`. Tests run with runtime warnings treated as errors.

Segregation checks measure species masses, local fraction sums, bounds, and pure-
layer mixing without clipping or normalisation. Run `python -m
benchmarks.segregation_drift` to reproduce the conservation report. The
`segregation_normalise` option is unnecessary for the conservative operator
and can conceal an invalid configuration.

Materials share a momentum grid: this does not provide separate-body contact or
independent material velocity fields. Traction scaling uses the first material's
uniform reference density. The [numerical methods](numerical-foundations.md)
describe diagnostics and the verification suite; unsupported registry names fail
explicitly rather than selecting a different law.

The user manual builds the current reStructuredText and Markdown sources together.
See [documentation maintenance](README.md) for build instructions and the checked-in
HTML layout.

## Example compatibility

The table below is generated from `python -m benchmarks.audit_examples`: five
bounded timesteps, initial `dt <= 1e-6`, output disabled. Dynamic timestep rules
remain enabled where configured. Durations, loading, and output frequency are configured in each input module. Unsupported laws are deliberately not substituted with different laws.

| Example command | Status |
|---|---|
| `python MPM.py avalanche` | supported smoke-tested |
| `python MPM.py ball` | supported smoke-tested |
| `python MPM.py bi_conveyor` | supported smoke-tested |
| `python MPM.py bi_dam_break 11 mix 0 2` | supported smoke-tested |
| `python MPM.py bi_double_insets 11` | supported smoke-tested |
| `python MPM.py bi_drum 11 2` | supported smoke-tested |
| `python MPM.py bi_seg_test 22 2 11` | supported smoke-tested |
| `python MPM.py bi_square_drum 11 2` | supported smoke-tested |
| `python MPM.py bi_vertical_chute 2 11` | supported smoke-tested |
| `python MPM.py cidr_collapse` | supported smoke-tested |
| `python MPM.py cidr_incline` | supported smoke-tested |
| `python MPM.py dam_break` | supported smoke-tested |
| `python MPM.py domain_column` | supported smoke-tested |
| `python MPM.py dp_rate_unit_test` | supported smoke-tested |
| `python MPM.py dp_unit_test` | supported smoke-tested |
| `python MPM.py francois 11` | supported smoke-tested |
| `python MPM.py granular_collapse` | supported smoke-tested |
| `python MPM.py granular_column` | supported smoke-tested |
| `python MPM.py granular_square_drum 11 2 cpu` | supported smoke-tested |
| `python MPM.py modern_ball` | supported smoke-tested |
| `python MPM.py poly_seg_test 22` | supported smoke-tested |
| `python MPM.py pouliquen_unit_test` | supported smoke-tested |
| `python MPM.py seg_unit_test` | supported smoke-tested |
| `python MPM.py silo 11` | supported smoke-tested |
| `python MPM.py slope_stab` | supported smoke-tested |
| `python MPM.py triaxial` | supported smoke-tested |
