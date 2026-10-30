# Feature reference

Every setting an input can change, with its default, the values it accepts,
which solver and backend support it, and the page that describes it in
detail. Settings live on five objects: `P` (the run), `P.G` (the grid), `P.B`
(boundaries), `P.S` (materials) and `P.O` (output). `initialise.prepare`
first calls `settings.normalise(P)`, which rejects any name that is not a
setting (suggesting the nearest one for a misspelling), fills in the defaults
listed here and checks each value. An unsupported combination raises an error
before the run starts. Names that begin with an underscore are left alone, so
an input can keep its own helpers as `_name`; methods defined on the classes,
such as callbacks and `critical_time`, are not settings.

The run, grid, boundary and output tables below are generated from
`settings.SETTINGS` by `python -m settings`; a test checks that they match.

Support columns use these abbreviations: explicit is the default explicit
solver, granular the rigid-grain solver on the CPU and Metal the rigid-grain
solver on the GPU. "Rejected" means `prepare` raises an error.

## Run settings

<!-- settings:P -->
| Setting | Default | Type | Meaning |
| --- | --- | --- | --- |
| `dt` | required | s | Timestep; the initial one under dynamic stepping. |
| `t_f` | required | s | Final time. |
| `savetime` | required | s | Interval between saved outputs. |
| `solver` | `'explicit'` | choice | `explicit` (constitutive laws) or `granular` (rigid grains). |
| `backend` | `'cpu'` | choice | `cpu`, or `metal` for the granular solver on Apple silicon. |
| `time_stepping` | `'static'` | choice | `dynamic` sets `dt = cfl x min(limits)` from advection, segregation, the elastic and viscous estimates and each material's `critical_time(P)`. |
| `cfl` | `0.1` | number | Courant number for dynamic stepping. |
| `timestep_policy` | `'report'` | choice | `enforce` raises when `dt` exceeds the stability limits. |
| `gravity` | `0.0` | m/s^2 | Gravity, signed (negative is down the y axis at zero angle). |
| `gravity_ramp` | `None` | function | Factor on gravity as a function of time, e.g. to ramp it up. |
| `gravity_angle` | `0.0` | rad or function | Angle of gravity from -y (a slope, or a rotating drum as a function of time). |
| `mass_tolerance` | `1e-10` | kg | Grid nodes with less mass are treated as empty. |
| `force_scale` | `1.0` | number | Factor on the nodal force each step (1: none; numerical damping). |
| `initial_stress` | `None` | choice | `lithostatic` from gravity, depth and slope; `compression` from `initial_pressure`. |
| `initial_pressure` | `None` | Pa | Pressure of `initial_stress = 'compression'`. |
| `initial_flow` | `None` | choice | Steady profiles of a granular or viscous layer. |
| `kinetic_pressure_scale` | `10.0` | number | Scale of the kinetic pressure of `initial_flow`. |
| `interpolation` | `'linear'` | choice | Shape functions. |
| `quadrature` | `'point'` | choice | `domain` integrates over each particle's domain. |
| `domain_max_support` | `256` | integer | Largest stencil of domain quadrature. |
| `domain_polygon` | `'batched'` | choice | Domain clipping implementation. |
| `transfer` | `'flip'` | choice | Particle-grid velocity transfer. |
| `pic_fraction` | `0.0` | number | PIC blended into FLIP. |
| `fbar` | `False` | bool | F-bar volumetric averaging (explicit solver). |
| `segregate` | `False` | bool | Grain-size segregation on the grid. |
| `segregation_rate` | `0.0` | number | Segregation rate coefficient. |
| `segregation_normalise` | `False` | bool | Renormalise grain fractions after each step. |
| `smooth_shear_rate` | `False` | bool | Smooth the shear rate used by segregation. |
| `diagnostics` | `True` | bool | Collect conservation and solver diagnostics. |
| `diagnostics_stride` | `1` | integer | Record diagnostics every this many steps. |
| `output_dir` | `None` | path | Output directory (default `im/<example>/<law>/`). |
| `granular_solver` | `'cidr'` | choice | `cidr`: critical-state CIDR. `velocity`, `hybrid`, `block`: Coulomb friction with a hard packing limit. |
| `granular_tolerance` | `1e-09` | number | Coulomb solvers: natural-residual tolerance (Metal 1e-4). |
| `granular_max_iterations` | `200` | integer | Coulomb solvers: iteration budget. |
| `granular_stopping` | `'residual'` | choice | Coulomb `velocity` solver stopping rule. |
| `granular_increment` | `0.001` | number | Threshold of `increment` stopping. |
| `granular_stress` | `'intrinsic'` | choice | Stress space. |
| `granular_stress_basis` | `'nodal'` | choice | Stress discretisation. |
| `granular_stress_coarsening` | `2` | integer | Coarse stress spacing in velocity cells. |
| `granular_viscosity` | `0.0` | Pa s | Coulomb solvers: optional viscous regularisation. |
| `granular_max_nodes` | `256` | integer | CPU guard on active stress nodes. |
| `granular_cidr_tolerance` | `1e-05` | number | CIDR: velocity change between outer steps (Metal 1e-4). |
| `granular_cidr_gamma` | `None` | number | CIDR: envelope smoothing in units of `dt / rho`; None is 0.075 on the CPU and 1.0 on Metal. |
<!-- /settings:P -->

The body acceleration is `g [sin θ, cos θ, 0]`, with `g` the gravity times the
ramp factor and `θ` the gravity angle, so a slope of 30 degrees is
`gravity_angle = np.radians(30)`. `gravity_ramp`, `gravity_angle` and the boundary tractions may
be functions of time; they are evaluated at the start of each step. A rotating
drum sets `gravity_angle = lambda t: omega * t`, and a gravity that rises
smoothly over `t_c` sets `gravity_ramp = lambda t: 1 - np.exp(-3 * t**2 /
t_c**2)`. An input may not define `update_forces`.

`force_scale` multiplies the nodal force each step; it is a numerical scaling,
not physical damping, and the older `damping = True` is `force_scale = 0.7`.
The granular solver and Metal require `force_scale = 1` and `initial_stress =
None`. `initial_pressure` is also the wall pressure of `initial_flow =
"steady_poiseuille"`; `initial_flow = "steady"` supports `viscous` and
`pouliquen` materials and uses `kinetic_pressure_scale` only with
`pouliquen`, to set the initial segregation potential.

`MPM.run(P, G, L, max_steps=..., output=...)` stops after a number of steps.
`O.save_checkpoint` writes a restart file at each save time; see
[outputs](outputs.rst) for checkpoint and restart.

### Interpolation, transfer and quadrature

| Setting | explicit | granular | Metal |
| --- | --- | --- | --- |
| `interpolation` | both | `linear` | `linear` |
| `transfer` | `flip`, `pic`, `apic` (quadratic only) | `flip`, `pic`, `apic_gradient` | `flip`, `pic`, `apic_gradient` |
| `pic_fraction` (0 to 1) | with `flip` | with `flip` | with `flip` |
| `quadrature` | both (`domain` needs `domain_size` per material) | `point` | `point` |
| `domain_max_support`, `domain_polygon` | `domain` only | — | — |

`apic_gradient` is restricted to the granular solver. `domain_polygon` selects
how rotated or sheared particle domains are integrated: `batched` integrates all
of them together, `reference` one at a time; they agree to roundoff. See
[numerical-foundations](numerical-foundations.md) and [quadrature](quadrature.md).

## Grid

<!-- settings:G -->
| Setting | Default | Type | Meaning |
| --- | --- | --- | --- |
| `x_min` | required | m | Left edge of the grid. |
| `x_max` | required | m | Right edge of the grid. |
| `y_min` | required | m | Bottom edge of the grid. |
| `y_max` | required | m | Top edge of the grid. |
| `nx` | required | integer | Grid nodes in x. |
| `ny` | required | integer | Grid nodes in y. |
| `thickness` | `1.0` | m | Out-of-plane thickness. |
| `grain_sizes` | `(0.001,)` | m | Grain size of each class (segregation, mu(I) laws). |
| `top_gap` | `0.0` | m | Free-surface offset of `initial_stress = 'lithostatic'`. |
<!-- /settings:G -->

`nx` and `ny` count nodes, not cells; the spacing `G.dx`, `G.dy` follows from
them and is set by the grid at run time, so a material's `critical_time(P)` may
read it. The number of size classes is `len(G.grain_sizes)`.

## Materials

`P.S` is one material namespace or a list. Every material needs `law`, `n`,
`X`, `Y`, `rho` (kg/m³) and `A` (area per particle, m²). Any material may also
set `initial_v` (a uniform or per-particle initial velocity), `phi` (the
fractions of each size class), `domain_size` (for domain quadrature),
`initial_C` and the inlet fields `inlet_rate`, `inlet_gap` and `inlet_y` (see
[Boundaries](#boundaries)). A field the law does not read is rejected. See
[tutorial](tutorial.rst) for geometry and [materials](materials.rst) for each
law's parameters.

| `law` | Solver | Fields of the law |
| --- | --- | --- |
| `rigid` | explicit | none |
| `elastic` | explicit | `K`, `G` |
| `von_mises` | explicit | `K`, `G`, `yield_stress` |
| `viscous` | explicit | `K`, `eta`, `eta_v` |
| `bingham` | explicit | `K`, `eta`, `eta_0`, `yield_stress`, `eta_v`, `gamma_c` |
| `pouliquen` | explicit | `K`, `rho_s`, `mu_s`, `mu_d`, `I_0`, `eta_max`, `phi_flow` |
| `linear_mu` | explicit | `K`, `rho_s`, `mu_s`, `b`, `eta_max` |
| `dp` | explicit | `K`, `G`, `mu_s`, `beta` |
| `dp_rate` | explicit | `K`, `G`, `mu_s`, `beta`, `t_star` |
| `granular` | granular, Metal | `rho_s`, `mu_s`, `phi_max`, `mu_d`, `I_0`, `grain_diameter`, `phi_min`, `cidr_deficit`, `cidr_p_ref` |

Only one `granular` material is allowed, and it cannot be combined with other
laws. A `critical_time(P)` method on any material bounds the dynamic timestep.
The Coulomb granular solvers read `rho_s`, `mu_s` and `phi_max` (at most 1).
With `granular_solver = "cidr"` the material follows critical-state CIDR (see
[Rigid-grain solver](#rigid-grain-solver)) and needs the other fields as well:

| Field | Meaning |
| --- | --- |
| `mu_s`, `mu_d`, `I_0` | `mu(I) = mu_s + (mu_d - mu_s) I / (I_0 + I)`, with `mu_s < mu_d` |
| `grain_diameter` | Grain diameter `d` (m); `I = 2 d |D| / sqrt(p / rho_s)` |
| `phi_min`, `phi_max` | Packing range of `C(phi)`, which diverges at `phi_max` |
| `cidr_deficit`, `cidr_p_ref` | Material on the cap at rest at pressure `cidr_p_ref` sits `cidr_deficit` below `phi_max`; this sets the pressure scale of `C` |

## Boundaries

<!-- settings:B -->
| Setting | Default | Type | Meaning |
| --- | --- | --- | --- |
| `walls` | `()` | edges | Domain edges with walls: any of bottom, top, left, right. |
| `periodic` | `False` | bool | Periodic in x. |
| `wall_contact` | `'slip'` | choice | `slip`, `no_slip` (all walls), `no_slip_bottom`, or `coulomb` with `wall_friction`. |
| `wall_friction` | `None` | number | Coulomb wall friction coefficient. |
| `one_sided` | `False` | bool | Domain walls release a node moving away from them. |
| `obstacles` | `()` | names | Interior walls: any of centre_wall, two_walls, box, silo_left, silo_bottom. |
| `outlets` | `()` | edges | Edges that remove material: left, bottom. |
| `inlets` | `()` | edges | Edges that add material (see the material's inlet fields): right, top. |
| `conveyor_speed` | `None` | m/s | Speed of a conveyor along the bottom. |
| `traction_vertical` | `None` | Pa or function | Traction on the bottom and top edges, constant or a function of time. |
| `traction_horizontal` | `None` | Pa or function | Traction on the left and right edges, constant or a function of time. |
<!-- /settings:B -->

Walls are applied at grid nodes; there is no contact between separate bodies or
along curved surfaces. An inlet adds a column of particles at `x_max -
inlet_gap[0]` and heights `inlet_y` (right inlet), or one particle at mid-width
and `y_max - inlet_gap[1]` (top inlet), every `inlet_rate` seconds; the new
particles copy the rest of the material. A left outlet removes particles within
one cell of the left edge and a bottom outlet those below the bottom edge.

| Setting | explicit | granular | Metal |
| --- | --- | --- | --- |
| `walls` | yes | yes | yes |
| `wall_contact = "no_slip"`, `"no_slip_bottom"` | yes | yes | yes |
| `wall_contact = "coulomb"` | explicit friction force | contact rows in the solve; needs `granular_solver = "velocity"` | yes |
| `one_sided` | yes | yes | yes, except with `cidr` |
| `periodic` | yes | yes | yes |
| `obstacles`: `centre_wall`, `two_walls`, `box`, `silo_left` | yes | yes | accepted |
| `obstacles`: `silo_bottom` | yes | yes | checked |
| `outlets`, `inlets` | yes | yes | yes, on the host (see [Metal backend](#metal-backend)) |
| `conveyor_speed` | yes | rejected | rejected |
| `traction_vertical`, `traction_horizontal` | yes, with the legacy density scaling (see [array-conversion](array-conversion.md)) | rejected | rejected |

"Accepted" means the Metal backend does not reject the setting but no Metal
gate case exercises it. `centre_wall` is an interior vertical wall in the lower
half at mid-width, `two_walls` two staggered interior vertical walls, `box` a
fixed interior box (grid nodes 4 to 8), and `silo_left` and `silo_bottom` a wall
with an opening at mid-height or mid-width. The conveyor moves the bottom edge
at `conveyor_speed`, away from the corners.

### Wall conditions

| Wall type | Settings | Normal | Tangential |
| --- | --- | --- | --- |
| Frictionless, two-sided | default (`wall_contact = "slip"`) | fixed both ways | free |
| No-slip, two-sided | `wall_contact = "no_slip"` | fixed both ways | fixed |
| No-slip bottom only | `wall_contact = "no_slip_bottom"` | fixed both ways | fixed on the bottom, free elsewhere |
| Frictionless or no-slip, one-sided | `one_sided = True` (with `slip` or `no_slip`) | fixed only while not moving away | as above, only while not moving away |
| Coulomb friction | `wall_contact = "coulomb"`, `wall_friction = μ_w` | unilateral contact (granular, velocity solver) | sticks below `μ_w` times the normal traction, otherwise slides |

With two-sided walls, mass mapped to wall nodes cannot move away from the wall,
so a row of grains within one cell of a wall can stay there when gravity turns
away from it. One-sided and frictional walls release it. A one-sided step
re-solves when the solve pushes a released node back into its wall (about 2.5
solves per step in the square drum); frictional walls need no re-solve. In the
explicit solver, `coulomb` instead limits the tangential wall force explicitly.
Interior walls (`obstacles`) stay two-sided. Details and checks are in
[granular-prototype](granular-prototype.md).

## Output and diagnostics

<!-- settings:O -->
| Setting | Default | Type | Meaning |
| --- | --- | --- | --- |
| `figures` | `()` | names | Figures at each save: any of continuum, material_points, gsd_material_points, gsd_grid, shear_rate. |
| `fields` | `()` | names | Fields saved at each save: any of s_bar, density, u, phi_material_points. |
| `measure_energy` | `False` | bool | Track elastic energy and plot it at the end. |
| `measure_stiffness` | `False` | bool | Measure the stiffness at the end. |
| `check_positions` | `False` | bool | Check particle positions against the grid each step. |
| `save_diagnostics` | `True` | bool | Write diagnostics.json. |
| `save_checkpoint` | `False` | bool | Write a checkpoint at each save. |
| `continuum_fig_size` | `None` | inches | Size of continuum figures. |
| `mp_fig_size` | `None` | inches | Size of material-point figures. |
| `at_start` | `None` | function | Callback `(L, P, G)` before the first step. |
| `after_every_timestep` | `None` | function | Callback `(L, P, G)` after each step. |
| `after_every_savetime` | `None` | function | Callback `(L, P, G)` at each save. |
| `at_end` | `None` | function | Callback `(L, P, G)` after the last step. |
<!-- /settings:O -->

`figures` and `fields` are written at each `savetime` into `output_dir`.
`measure_energy` is rejected by Metal and so is `check_positions`. `diagnostics`
and `diagnostics_stride` (run settings) control the conservation and work
budgets. The callbacks are described in [outputs](outputs.rst); an input may
also define them as methods of its output class. Metal runs reject a custom
`after_every_timestep`, because particle arrays are synchronised only at output
times.

## Rigid-grain solver

The `granular_*` run settings apply when `solver = "granular"`. The Metal
backend changes two defaults: `granular_tolerance = 1e-4` and
`granular_cidr_tolerance = 1e-4`.

| Setting | Metal | Notes |
| --- | --- | --- |
| `granular_solver` | `cidr`, `velocity` | `cidr` is critical-state CIDR (below), which needs its material fields. The others solve Coulomb friction with a hard packing limit: `velocity` is the compact iteration of [granular-mixed](granular-mixed.md), which Coulomb walls and `increment` stopping need; `hybrid` and `block` are small-problem reference solvers. |
| `granular_tolerance` | ≥ `1e-5` (`velocity`) | Natural-residual tolerance. |
| `granular_max_iterations` | yes | The velocity solver allows ten times this many sweeps. |
| `granular_stopping` | `residual` | `increment` stops on the strain-rate change per sweep, without certifying the residual. |
| `granular_stress` | `intrinsic` | `plane_strain` is not supported by `block`. |
| `granular_stress_basis` | `nodal` | `P0` uses one stress per cell; `coarse` a coarser nodal grid, whose cell counts must divide by `granular_stress_coarsening`. |
| `granular_viscosity` | `0` | Optional viscous regularisation. |
| `granular_max_nodes` | ignored | CPU guard on active stress nodes. |
| `granular_cidr_tolerance` | ≥ `1e-5` | `cidr`: velocity change between outer steps, relative to the larger of the velocity and the trial velocity. |
| `granular_cidr_gamma` | yes | `cidr`: Metal's iterative linear solves favour a smoother envelope, hence its default of 1.0. The converged step does not depend on it. |

The `cidr` solver needs `granular_stress = "intrinsic"` on the `"nodal"` basis and
does not support Coulomb walls; on Metal it does not support one-sided walls
(`B.one_sided`) either. It minimises the kinetic energy of the step plus the CIDR
dissipation potential over the nodal velocities, by the method of multipliers
with semismooth Newton subproblems (`cidr.py`; on Metal `mpm_gpu/cidr.py`, with
conjugate gradients for each Newton system). Jammed nodes are decided exactly,
so jammed material does not creep. The first step starts cold; later steps
start from the previous step's stress and take from none to a few Newton
iterations for material at rest or in fast flow, and several tens in slow,
nearly jammed flow on fine grids (see [gpu](gpu.md)). The packing per stress
node comes from the particles, and the hard packing limit is replaced by the
divergence of `C(phi)`. `C` is evaluated at the packing at the start of each step,
so a long step could compact a node past `phi_max`: the time to close the
largest relative packing gap at the last step's compaction rate is a timestep
limit (`packing`, scaled by `cfl`, used by dynamic stepping), and a step that
would close a gap stops the run with an error. The model, its derivation and the verification are in
[granular-cidr-design](granular-cidr-design.md).

Elastic and viscous timestep estimates do not apply; set `critical_time(P)` on
the material for dynamic stepping (the square drum uses the time for a grain to
fall one grid spacing). Convergence on refined, connected domains is an open
problem; see [gpu-handoff](gpu-handoff.md).

## Grain-size segregation

`segregate = True` transports the size fractions `phi` of each particle on the
grid, at a rate set by `segregation_rate` times the gradient of the segregation
potential. `segregation_normalise` clips each particle's fractions to [0, 1] and
rescales them to sum to one after each transport step. `smooth_shear_rate`
smooths the shear rate before use; Metal requires it to be `False`. The size
classes are `G.grain_sizes` and the initial composition is the material's `phi`.
Metal supports up to eight size classes. Rigid materials do not segregate. See
[segregation-conservation](segregation-conservation.md).

## Metal backend

`P.backend = "metal"` requires the `gpu` extra, Apple silicon and macOS 26 or
newer. It runs the complete granular step on the device in FP32. It requires
`solver = "granular"`, `granular_solver = "velocity"` or `"cidr"`,
`granular_stress = "intrinsic"`, `granular_stress_basis = "nodal"`, linear
interpolation, point quadrature and one `granular` material. Supported:
`flip`, `pic` and `apic_gradient` transfers; static or dynamic time stepping;
two-sided, one-sided and frictional walls; obstacles and silo walls; periodic
x; inlets and outlets; grid segregation. On a step where an inlet fires or a
particle reaches an outlet, the state is synchronised to the host,
`process_boundaries` runs as on the CPU and the runner is rebuilt for the new
particle count, so draining silos run mostly on the host path. Rejected:
tolerances below `1e-5`, `increment` stopping, `granular_viscosity`,
`conveyor_speed`, boundary tractions, `force_scale` other than 1,
`O.check_positions`, `O.measure_energy` and custom timestep controllers.
Results, timings and the validation gates are in [gpu](gpu.md).

## Examples

`python MPM.py --help` lists every example and `python MPM.py <name> --help`
its arguments. The [example matrix](array-conversion.md#example-compatibility)
gives each example's status. Examples for the features above:

| Feature | Example |
| --- | --- |
| Elastic, quadratic B-splines, APIC, dynamic stepping | `ball`, `modern_ball` |
| Segregation on an incline | `bi_seg_test 22 2 11` |
| Rotating drum (Pouliquen material) | `bi_square_drum 11 2` |
| Rigid grains: static column and collapse | `granular_column`, `granular_collapse` |
| Rigid grains in a rotating drum, CPU or Metal | `granular_square_drum 11 2 cpu`, `granular_square_drum 11 2 metal` |
| Finite-domain quadrature | `domain_column` |
