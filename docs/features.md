# Feature reference

Every setting an input can change, with its default, the values it accepts,
which solver and backend support it, and the page that describes it in
detail. Settings live on five objects: `P` (the run), `P.G` (the grid), `P.B`
(boundaries), `P.S` (materials) and `P.O` (output). A missing setting takes
the default listed here, applied by `initialise.prepare`. An unsupported
combination raises an error before the run starts rather than being ignored.

Support columns use these abbreviations: **explicit** is the default explicit
solver, **granular** the rigid-grain solver on the CPU and **Metal** the
rigid-grain solver on the GPU. "Rejected" means `prepare` raises an error.

## Solvers and backends

| Setting | Default | Values | Notes |
| --- | --- | --- | --- |
| `solver` | `"explicit"` | `"explicit"`, `"granular"` | `granular` is the rigid-grain contact solver; it requires one `granular` material. See [granular-prototype](granular-prototype.md). |
| `backend` | `"cpu"` | `"cpu"`, `"metal"` | `metal` runs the granular step on Apple silicon in FP32; see [Metal backend](#metal-backend) and [gpu](gpu.md). |
| `fbar` | `False` | bool | F-bar volumetric averaging against locking of nearly incompressible materials; explicit solver only. See [numerical-foundations](numerical-foundations.md). |

## Time and loading

| Setting | Default | Values | Notes |
| --- | --- | --- | --- |
| `dt` | required | > 0 s | Initial or fixed timestep. |
| `t_f` | required | ≥ `t` | End time (s). |
| `savetime` | required | > 0 s | Interval between saved outputs. |
| `time_stepping` | `"static"` | `"static"`, `"dynamic"` | `dynamic` sets `dt = CFL × min(limits)` from advection, segregation and each material's `critical_time(P)`, plus the elastic and viscous estimates. |
| `CFL` | `0.1` | > 0 | Scales the dynamic limits. |
| `timestep_policy` | `"report"` | `"report"`, `"enforce"` | `enforce` rejects a static `dt` above the estimated limits. |
| `max_g` | `0.0` | m/s² | Gravity magnitude; negative is downward at `theta = 0`. |
| `theta` | `0.0` | radians | Body acceleration is `g [sin θ, cos θ, 0]`. |
| `update_forces` | none | method | Called each step; may change `g` and `theta` from `self.t` (the drum rotates gravity this way). |
| `max_q` | `0.0` | Pa | Traction magnitude for `vertical_force` and `horizontal_force`. |
| `damping` | `False` | bool | Multiplies nodal forces by 0.7; a numerical scaling, not physical damping. Rejected by granular and Metal. |
| `M_tol` | `1e-10` | ≥ 0 kg | Nodes with less mass have no resolved velocity. |

`MPM.run(P, G, L, max_steps=..., output=...)` stops after a number of steps.
`O.save_checkpoint` writes a restart file at each save time; see
[outputs](outputs.rst) for checkpoint and restart.

## Grid

| Setting | Default | Values | Notes |
| --- | --- | --- | --- |
| `G.nx`, `G.ny` | required | ≥ 2 | Node counts, not cell counts. |
| `G.x_m`, `G.x_M`, `G.y_m`, `G.y_M` | required | m | Domain extents; spacing follows from the node counts. |
| `G.thickness` | `1.0` | m | Out-of-plane thickness. |
| `G.ns` | `1` | ≥ 1 | Number of grain-size classes. |
| `G.s` | `[0.001]` | `ns` positive diameters (m) | Representative size of each class. |

## Interpolation, transfer and quadrature

| Setting | Default | Values | explicit | granular | Metal |
| --- | --- | --- | --- | --- | --- |
| `interpolation` | `"linear"` | `"linear"`, `"quadratic"` | both | `linear` | `linear` |
| `transfer` | `"flip"` | `"flip"`, `"pic"`, `"apic"`, `"apic_gradient"` | `flip`, `pic`, `apic` (quadratic only) | `flip`, `pic`, `apic_gradient` | `flip`, `pic`, `apic_gradient` |
| `pic_fraction` | `0.0` | 0 to 1 | with `flip` | with `flip` | with `flip` |
| `quadrature` | `"point"` | `"point"`, `"domain"` | both (`domain` needs `domain_size` per material) | `point` | `point` |
| `domain_max_support` | `256` | integer ≥ 4 | `domain` only | — | — |
| `domain_polygon` | `"batched"` | `"batched"`, `"reference"` | `domain` only | — | — |

`apic_gradient` is restricted to the granular solver. `domain_polygon` selects
how rotated or sheared particle domains are integrated: `batched` integrates all
of them together, `reference` one at a time; they agree to roundoff. See
[numerical-foundations](numerical-foundations.md) and [quadrature](quadrature.md).

## Materials

`P.S` is one material namespace or a list. Every material needs `law`, `n`,
`X`, `Y`, `rho` (kg/m³) and `A` (area per particle, m²); `initial_v` and, for
mixtures, `phi` are optional. See [tutorial](tutorial.rst) for geometry and
[materials](materials.rst) for each law's parameters.

| `law` | Solver | Required parameters |
| --- | --- | --- |
| `rigid` | explicit | optional `initial_v` |
| `elastic` | explicit | `K`, `G` |
| `viscous` | explicit | `K`, `eta`, `eta_v` |
| `von_mises` | explicit | `K`, `G`, `yield_stress` |
| `pouliquen` | explicit | `K`, `rho_s`, `mu_s`, `mu_d`, `I_0`, `eta_max` |
| `bingham` | explicit | `K`, `eta`, `eta_0`, `yield_stress`; optional `eta_v` |
| `dp` | explicit | `K`, `G`, `mu_s`, `beta` |
| `dp_rate` | explicit | as `dp`, plus `t_star` |
| `linear_mu` | explicit | `K`, `rho_s`, `mu_s`, `b`, `eta_max` |
| `granular` | granular, Metal | `rho_s`, `mu_s`, `phi_max` (≤ 1); optional `critical_time(P)` |

Only one `granular` material is allowed, and it cannot be combined with other
laws. A `critical_time(P)` method on any material bounds the dynamic timestep.

With `granular_solver = "cidr"` the `granular` material follows critical-state
CIDR (see [Rigid-grain solver](#rigid-grain-solver)) and needs these fields as
well:

| Field | Meaning |
| --- | --- |
| `mu_s`, `mu_d`, `I_0` | `mu(I) = mu_s + (mu_d - mu_s) I / (I_0 + I)`, with `mu_s < mu_d` |
| `grain_diameter` | Grain diameter `d` (m); `I = 2 d |D| / sqrt(p / rho_s)` |
| `phi_min`, `phi_max` | Packing range of `C(phi)`, which diverges at `phi_max` |
| `cidr_deficit`, `cidr_p_ref` | Material on the cap at rest at pressure `cidr_p_ref` sits `cidr_deficit` below `phi_max`; this sets the pressure scale of `C` |

## Initial state

| Setting | Default | Values | Notes |
| --- | --- | --- | --- |
| `pressure` | `None` | `None`, `"lithostatic"`, `"compression"`, `"non-zero"`, `"temperature"` | `lithostatic` uses gravity, depth and slope; `compression` applies `-initial_pressure` isotropically; the legacy `non-zero` and `temperature` modes start from a small elastic strain. Granular and Metal require `None`. |
| `initial_pressure` | — | Pa | Used by `compression` and `steady_poiseulle`. |
| `initial_flow` | `False` | `False`, `"steady"`, `"steady_poiseulle"` | Analytical velocity profiles; `steady` supports `viscous` and `pouliquen`. |
| `l` | `10.0` | — | Only used by `initial_flow = "steady"` with `pouliquen`, to set the initial segregation potential. |

Material `initial_v` sets a uniform or per-particle initial velocity.

## Boundaries

All flags default to `False`. Walls are applied at grid nodes; there is no
contact between separate bodies or along curved surfaces.

| Flag | Effect | explicit | granular | Metal |
| --- | --- | --- | --- | --- |
| `has_bottom`, `has_top`, `has_left`, `has_right` | Wall on that domain edge: fixes the normal velocity of its nodes | yes | yes | yes |
| `no_slip_bottom` | Also fixes horizontal velocity on the bottom edge | yes | yes | yes |
| `roughness` | No-slip walls: fixes both velocity components on every wall node | yes | yes | yes |
| `separable` | One-sided domain walls: a wall node moving away from its wall is released (see [wall conditions](#wall-conditions)) | a wall node moving away from its wall is released | yes | yes |
| `wall_mu` | Frictional walls with coefficient `wall_mu` (see [wall conditions](#wall-conditions)) | explicit friction force | contact rows in the solve; needs `granular_solver = "velocity"` | yes |
| `cyclic_lr` | Periodic left and right edges | yes | yes | yes |
| `wall` | Interior vertical wall in the lower half at mid-width | yes | yes | accepted |
| `two_walls` | Two staggered interior vertical walls | yes | yes | accepted |
| `box` | Fixed interior box of walls (grid nodes 4 to 8) | yes | yes | accepted |
| `silo_left`, `silo_bottom` | Wall with an opening at mid-height or mid-width | yes | yes | `silo_bottom` checked; `silo_left` accepted |
| `outlet_left`, `outlet_bottom` | Remove particles within one cell of the left edge, or below the bottom edge | yes | yes | yes, on the host (see below) |
| `inlet_right`, `inlet_top` | Add particles every `inlet_rate` seconds (material setting) | yes | yes | yes, on the host (see below) |
| `conveyor` | Horizontal velocity `v_0` (m/s) on the bottom edge, away from the corners | yes | rejected | rejected |
| `vertical_force`, `horizontal_force` | Edge traction `max_q` with the legacy density scaling (see [array-conversion](array-conversion.md)) | yes | rejected | rejected |

"Accepted" means the Metal backend does not reject the flag but no Metal gate
case exercises it.

### Wall conditions

| Wall type | Settings | Normal | Tangential |
| --- | --- | --- | --- |
| Frictionless, two-sided | default | fixed both ways | free |
| No-slip, two-sided | `roughness = True` | fixed both ways | fixed |
| Frictionless or no-slip, one-sided | `separable = True` (with or without `roughness`) | fixed only while not moving away | as above, only while not moving away |
| Coulomb friction | `wall_mu = μ_w` (granular, velocity solver) | unilateral contact | sticks below `μ_w` times the normal traction, otherwise slides |

With two-sided walls, mass mapped to wall nodes cannot move away from the wall,
so a row of grains within one cell of a wall can stay there when gravity turns
away from it. One-sided and frictional walls release it. A one-sided step
re-solves when the solve pushes a released node back into its wall (about 2.5
solves per step in the square drum); frictional walls need no re-solve. In the
explicit solver, `wall_mu` instead limits the tangential wall force explicitly,
in place of the no-slip condition of `roughness`. Interior walls (`wall`, `two_walls`,
`box`, silos) stay two-sided. Details and checks are in
[granular-prototype](granular-prototype.md).

## Rigid-grain solver

These settings apply when `solver = "granular"`. The Metal backend changes two
defaults: `granular_tolerance = 1e-4` and `granular_cidr_tolerance = 1e-4`.

| Setting | Default | Values | Metal | Notes |
| --- | --- | --- | --- | --- |
| `granular_solver` | `"cidr"` | `"cidr"`, `"velocity"`, `"hybrid"`, `"block"` | `cidr`, `velocity` | `cidr` is critical-state CIDR (below), which needs its material fields. The others solve Coulomb friction with a hard packing limit: `velocity` is the compact iteration of [granular-mixed](granular-mixed.md), which `wall_mu` and `increment` stopping need; `hybrid` and `block` are small-problem reference solvers. |
| `granular_tolerance` | `1e-9` | > 0 | ≥ `1e-5` (`velocity`) | Natural-residual tolerance. |
| `granular_max_iterations` | `200` | integer ≥ 1 | yes | The velocity solver allows ten times this many sweeps. |
| `granular_stopping` | `"residual"` | `"residual"`, `"increment"` | `residual` | `increment` stops on the strain-rate change per sweep, without certifying the residual. |
| `granular_increment` | `1e-3` | > 0 | — | Threshold for `increment` stopping. |
| `granular_stress` | `"intrinsic"` | `"intrinsic"`, `"plane_strain"` | `intrinsic` | `plane_strain` is not supported by `block`. |
| `granular_stress_basis` | `"nodal"` | `"nodal"`, `"P0"`, `"coarse"` | `nodal` | `P0` uses one stress per cell; `coarse` a coarser nodal grid. |
| `granular_stress_coarsening` | `2` | integer ≥ 2 | — | Coarse spacing in velocity cells; cell counts must divide by it. |
| `granular_viscosity` | `0.0` | ≥ 0 Pa s | `0` | Optional viscous regularisation. |
| `granular_max_nodes` | `256` | integer ≥ 1 | ignored | CPU guard on active stress nodes. |
| `granular_cidr_tolerance` | `1e-5` | > 0 | ≥ `1e-5` | `cidr`: velocity change between outer steps, relative to the larger of the velocity and the trial velocity. |
| `granular_cidr_gamma` | `None` | > 0 or `None` | yes | `cidr`: smoothing of the Moreau envelope in units of `dt / rho`; `None` is 0.075 on the CPU and 1.0 on Metal, whose iterative linear solves favour a smoother envelope. The converged step does not depend on it. |

The `cidr` solver needs `granular_stress = "intrinsic"` on the `"nodal"` basis and
does not support `wall_mu`; on Metal it does not support one-sided walls
(`B.separable`) either. It minimises the kinetic energy of the step plus the CIDR
dissipation potential over the nodal velocities, by the method of multipliers
with semismooth Newton subproblems (`cidr.py`; on Metal `mpm_gpu/cidr.py`, with
conjugate gradients for each Newton system). Jammed nodes are decided exactly,
so jammed material does not creep. The first step starts cold; later steps
start from the previous step's stress and take from none to a few Newton
iterations for material at rest or in fast flow, and several tens in slow,
nearly jammed flow on fine grids (see [gpu](gpu.md)). The packing per stress node comes from the particles, and
the hard packing limit is replaced by the divergence of `C(phi)`. The model,
its derivation and the verification are in
[granular-cidr-design](granular-cidr-design.md).

Elastic and viscous timestep estimates do not apply; set `critical_time(P)` on
the material for dynamic stepping (the square drum uses the time for a grain to
fall one grid spacing). Convergence on refined, connected domains is an open
problem; see [gpu-handoff](gpu-handoff.md).

## Grain-size segregation

| Setting | Default | Values | Notes |
| --- | --- | --- | --- |
| `segregate` | `False` | bool | Grid segregation; `segregate_grid` is an alias. `segregate_mp` is rejected. |
| `c` | `0.0` | ≥ 0 | Segregation coefficient multiplying the gradient of the segregation potential. |
| `normalise_phi` | `False` | bool | Clip each particle's fractions to [0, 1] and rescale them to sum to one after each transport step. |
| `smooth_gamma_dot` | `False` | bool | Smooth the shear rate before use; Metal requires `False`. |
| `G.ns`, `G.s`, material `phi` | see [Grid](#grid) | | Size classes and initial composition. |

Metal supports up to eight size classes. Rigid materials do not segregate.
See [segregation-conservation](segregation-conservation.md).

## Output and diagnostics

| Setting | Default | Notes |
| --- | --- | --- |
| `save_dir` | `supername`, else `im/<example>/<law>/` | Output folder. |
| `O.plot_continuum`, `O.plot_material_points`, `O.plot_gsd_mp`, `O.plot_gsd_grid` | `False` | Figures at each `savetime`. |
| `O.save_u`, `O.save_density`, `O.save_s_bar`, `O.save_phi_MP` | `False` | NumPy arrays at each `savetime`. |
| `O.measure_energy`, `O.measure_stiffness` | `False` | Energy history and stiffness; Metal rejects `measure_energy`. |
| `O.save_diagnostics` | `True` | Writes `data/diagnostics.json`. |
| `O.save_checkpoint` | `False` | Writes `data/checkpoint.pkl` at each save time and at the end; see [outputs](outputs.rst). |
| `O.check_positions` | `False` | Assert each step that particle positions are consistent with the grid stencil; rejected by Metal. |
| `diagnostics` | `True` | Conservation and work budgets each step. |
| `diagnostics_stride` | `1` | Record every n-th step. |

Callbacks (`at_start`, `after_every_timestep`, `after_every_savetime`,
`at_end`) are described in [outputs](outputs.rst). Metal runs reject a custom
`after_every_timestep`, because particle arrays are synchronised only at output
times.

## Metal backend

`P.backend = "metal"` requires the `gpu` extra, Apple silicon and macOS 26 or
newer. It runs the complete granular step on the device in FP32. It requires
`solver = "granular"`, `granular_solver = "velocity"` or `"cidr"`,
`granular_stress = "intrinsic"`, `granular_stress_basis = "nodal"`, linear
interpolation, point quadrature and one `granular` material. Supported:
`flip`, `pic` and `apic_gradient` transfers; static or dynamic time stepping;
two-sided, one-sided and frictional walls; interior and silo walls; periodic
x; inlets and outlets; grid segregation. On a step where an inlet fires or a
particle reaches an outlet, the state is synchronised to the host,
`process_boundaries` runs as on the CPU and the runner is rebuilt for the new
particle count, so draining silos run mostly on the host path. Rejected:
tolerances below `1e-5`, `increment` stopping, `granular_viscosity`,
conveyors, traction flags, `damping`, `O.check_positions`,
`O.measure_energy` and custom timestep controllers.
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
