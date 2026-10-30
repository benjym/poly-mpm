"""User settings: what an input may set, their defaults and checks, and the solver's internal names.

An input defines a Params object P with namespaces P.G (grid), P.B (boundaries),
P.O (output) and P.S (one material or a list of them). normalise(P) runs first in
initialise.prepare: it rejects any name that is not a setting (with a suggestion for
a likely misspelling), fills in defaults, checks values, and translates the settings
into the attributes the solver reads (for example B.walls into has_bottom/has_left,
gravity into max_g, g and theta each step). Those internal attributes are not
settings. Names beginning with an underscore are left alone, for helpers an input
computes; methods and callbacks defined on the classes are not instance settings.

The tables in docs/features.md are generated from SETTINGS (python -m settings).
"""

import difflib
from dataclasses import dataclass

import numpy as np

REQUIRED = object()
EDGES = ("bottom", "top", "left", "right")
OBSTACLES = ("centre_wall", "two_walls", "box", "silo_left", "silo_bottom")
FIGURES = ("continuum", "material_points", "gsd_material_points", "gsd_grid", "shear_rate")
FIELDS = ("s_bar", "density", "u", "phi_material_points")


@dataclass(frozen=True)
class Setting:
    name: str
    default: object
    kind: str  # shown in the docs
    doc: str
    check: object = None  # value -> error message or None


def _number(lower=None, strict=False, integer=False):
    def check(value):
        if isinstance(value, bool) or not np.isscalar(value) or not np.isfinite(value):
            return "must be a finite number"
        if integer and int(value) != value:
            return "must be an integer"
        if lower is not None and (value <= lower if strict else value < lower):
            return f"must be {'>' if strict else '>='} {lower}"
        return None

    return check


def _boolean(value):
    return None if isinstance(value, (bool, np.bool_)) else "must be True or False"


def _choice(*choices):
    def check(value):
        return None if value in choices else f"must be one of {', '.join(map(repr, choices))}"

    return check


def _subset(choices):
    def check(value):
        if isinstance(value, str) or not all(item in choices for item in value):
            return f"must be a collection drawn from {', '.join(map(repr, choices))}"
        return None

    return check


def _number_or_function(value):
    if callable(value):
        return None
    return _number()(value)


def _optional(check):
    return lambda value: None if value is None else check(value)


def _function(value):
    return None if callable(value) else "must be a function of time"


SETTINGS = {
    "P": [
        Setting(
            "dt", REQUIRED, "s", "Timestep; the initial one under dynamic stepping.", _number(0, strict=True)
        ),
        Setting("t_f", REQUIRED, "s", "Final time.", _number(0)),
        Setting("savetime", REQUIRED, "s", "Interval between saved outputs.", _number(0, strict=True)),
        Setting(
            "solver",
            "explicit",
            "choice",
            "`explicit` (constitutive laws) or `granular` (rigid grains).",
            _choice("explicit", "granular"),
        ),
        Setting(
            "backend",
            "cpu",
            "choice",
            "`cpu`, or `metal` for the granular solver on Apple silicon.",
            _choice("cpu", "metal"),
        ),
        Setting(
            "time_stepping",
            "static",
            "choice",
            "`dynamic` sets `dt = cfl x min(limits)` from advection, segregation, the elastic and viscous "
            "estimates and each material's `critical_time(P)`.",
            _choice("static", "dynamic"),
        ),
        Setting("cfl", 0.1, "number", "Courant number for dynamic stepping.", _number(0, strict=True)),
        Setting(
            "timestep_policy",
            "report",
            "choice",
            "`enforce` raises when `dt` exceeds the stability limits.",
            _choice("report", "enforce"),
        ),
        Setting(
            "gravity", 0.0, "m/s^2", "Gravity, signed (negative is down the y axis at zero angle).", _number()
        ),
        Setting(
            "gravity_ramp",
            None,
            "function",
            "Factor on gravity as a function of time, e.g. to ramp it up.",
            _optional(_function),
        ),
        Setting(
            "gravity_angle",
            0.0,
            "rad or function",
            "Angle of gravity from -y (a slope, or a rotating drum as a function of time).",
            _number_or_function,
        ),
        Setting("mass_tolerance", 1e-10, "kg", "Grid nodes with less mass are treated as empty.", _number(0)),
        Setting(
            "force_scale",
            1.0,
            "number",
            "Factor on the nodal force each step (1: none; numerical damping).",
            _number(0, strict=True),
        ),
        Setting(
            "initial_stress",
            None,
            "choice",
            "`lithostatic` from gravity, depth and slope; `compression` from `initial_pressure`.",
            _choice(None, "lithostatic", "compression"),
        ),
        Setting(
            "initial_pressure",
            None,
            "Pa",
            "Pressure of `initial_stress = 'compression'`.",
            _optional(_number(0)),
        ),
        Setting(
            "initial_flow",
            None,
            "choice",
            "Steady profiles of a granular or viscous layer.",
            _choice(None, "steady", "steady_poiseuille"),
        ),
        Setting(
            "kinetic_pressure_scale",
            10.0,
            "number",
            "Scale of the kinetic pressure of `initial_flow`.",
            _number(0),
        ),
        Setting("interpolation", "linear", "choice", "Shape functions.", _choice("linear", "quadratic")),
        Setting(
            "quadrature",
            "point",
            "choice",
            "`domain` integrates over each particle's domain.",
            _choice("point", "domain"),
        ),
        Setting(
            "domain_max_support",
            256,
            "integer",
            "Largest stencil of domain quadrature.",
            _number(4, integer=True),
        ),
        Setting(
            "domain_polygon",
            "batched",
            "choice",
            "Domain clipping implementation.",
            _choice("batched", "reference"),
        ),
        Setting(
            "transfer",
            "flip",
            "choice",
            "Particle-grid velocity transfer.",
            _choice("flip", "pic", "apic", "apic_gradient"),
        ),
        Setting("pic_fraction", 0.0, "number", "PIC blended into FLIP.", _number(0)),
        Setting("fbar", False, "bool", "F-bar volumetric averaging (explicit solver).", _boolean),
        Setting("segregate", False, "bool", "Grain-size segregation on the grid.", _boolean),
        Setting("segregation_rate", 0.0, "number", "Segregation rate coefficient.", _number(0)),
        Setting(
            "segregation_normalise", False, "bool", "Renormalise grain fractions after each step.", _boolean
        ),
        Setting("smooth_shear_rate", False, "bool", "Smooth the shear rate used by segregation.", _boolean),
        Setting("diagnostics", True, "bool", "Collect conservation and solver diagnostics.", _boolean),
        Setting(
            "diagnostics_stride",
            1,
            "integer",
            "Record diagnostics every this many steps.",
            _number(1, integer=True),
        ),
        Setting(
            "output_dir",
            None,
            "path",
            "Output directory (default `im/<example>/<law>/`).",
            _optional(lambda v: None if isinstance(v, str) else "must be a path"),
        ),
        Setting(
            "granular_solver",
            "cidr",
            "choice",
            "`cidr`: critical-state CIDR. `velocity`, `hybrid`, `block`: Coulomb friction with a hard "
            "packing limit.",
            _choice("cidr", "velocity", "hybrid", "block"),
        ),
        Setting(
            "granular_tolerance",
            1e-9,
            "number",
            "Coulomb solvers: natural-residual tolerance (Metal 1e-4).",
            _number(0, strict=True),
        ),
        Setting(
            "granular_max_iterations",
            200,
            "integer",
            "Coulomb solvers: iteration budget.",
            _number(1, integer=True),
        ),
        Setting(
            "granular_stopping",
            "residual",
            "choice",
            "Coulomb `velocity` solver stopping rule.",
            _choice("residual", "increment"),
        ),
        Setting(
            "granular_increment",
            1e-3,
            "number",
            "Threshold of `increment` stopping.",
            _number(0, strict=True),
        ),
        Setting(
            "granular_stress", "intrinsic", "choice", "Stress space.", _choice("intrinsic", "plane_strain")
        ),
        Setting(
            "granular_stress_basis",
            "nodal",
            "choice",
            "Stress discretisation.",
            _choice("nodal", "P0", "coarse"),
        ),
        Setting(
            "granular_stress_coarsening",
            2,
            "integer",
            "Coarse stress spacing in velocity cells.",
            _number(2, integer=True),
        ),
        Setting(
            "granular_viscosity", 0.0, "Pa s", "Coulomb solvers: optional viscous regularisation.", _number(0)
        ),
        Setting(
            "granular_max_nodes",
            256,
            "integer",
            "CPU guard on active stress nodes.",
            _number(1, integer=True),
        ),
        Setting(
            "granular_cidr_tolerance",
            1e-5,
            "number",
            "CIDR: velocity change between outer steps (Metal 1e-4).",
            _number(0, strict=True),
        ),
        Setting(
            "granular_cidr_gamma",
            None,
            "number",
            "CIDR: envelope smoothing in units of `dt / rho`; None is 0.075 on the CPU and 1.0 on Metal.",
            _optional(_number(0, strict=True)),
        ),
    ],
    "G": [
        Setting("x_min", REQUIRED, "m", "Left edge of the grid.", _number()),
        Setting("x_max", REQUIRED, "m", "Right edge of the grid.", _number()),
        Setting("y_min", REQUIRED, "m", "Bottom edge of the grid.", _number()),
        Setting("y_max", REQUIRED, "m", "Top edge of the grid.", _number()),
        Setting("nx", REQUIRED, "integer", "Grid nodes in x.", _number(2, integer=True)),
        Setting("ny", REQUIRED, "integer", "Grid nodes in y.", _number(2, integer=True)),
        Setting("thickness", 1.0, "m", "Out-of-plane thickness.", _number(0, strict=True)),
        Setting("grain_sizes", (0.001,), "m", "Grain size of each class (segregation, mu(I) laws).", None),
        Setting("top_gap", 0.0, "m", "Free-surface offset of `initial_stress = 'lithostatic'`.", _number()),
    ],
    "B": [
        Setting(
            "walls", (), "edges", "Domain edges with walls: any of bottom, top, left, right.", _subset(EDGES)
        ),
        Setting("periodic", False, "bool", "Periodic in x.", _boolean),
        Setting(
            "wall_contact",
            "slip",
            "choice",
            "`slip`, `no_slip` (all walls), `no_slip_bottom`, or `coulomb` with `wall_friction`.",
            _choice("slip", "no_slip", "no_slip_bottom", "coulomb"),
        ),
        Setting(
            "wall_friction",
            None,
            "number",
            "Coulomb wall friction coefficient.",
            _optional(_number(0, strict=True)),
        ),
        Setting("one_sided", False, "bool", "Domain walls release a node moving away from them.", _boolean),
        Setting(
            "obstacles", (), "names", f"Interior walls: any of {', '.join(OBSTACLES)}.", _subset(OBSTACLES)
        ),
        Setting(
            "outlets", (), "edges", "Edges that remove material: left, bottom.", _subset(("left", "bottom"))
        ),
        Setting(
            "inlets",
            (),
            "edges",
            "Edges that add material (see the material's inlet fields): right, top.",
            _subset(("right", "top")),
        ),
        Setting("conveyor_speed", None, "m/s", "Speed of a conveyor along the bottom.", _optional(_number())),
        Setting(
            "traction_vertical",
            None,
            "Pa or function",
            "Traction on the bottom and top edges, constant or a function of time.",
            _optional(_number_or_function),
        ),
        Setting(
            "traction_horizontal",
            None,
            "Pa or function",
            "Traction on the left and right edges, constant or a function of time.",
            _optional(_number_or_function),
        ),
    ],
    "O": [
        Setting(
            "figures", (), "names", f"Figures at each save: any of {', '.join(FIGURES)}.", _subset(FIGURES)
        ),
        Setting(
            "fields", (), "names", f"Fields saved at each save: any of {', '.join(FIELDS)}.", _subset(FIELDS)
        ),
        Setting("measure_energy", False, "bool", "Track elastic energy and plot it at the end.", _boolean),
        Setting("measure_stiffness", False, "bool", "Measure the stiffness at the end.", _boolean),
        Setting(
            "check_positions", False, "bool", "Check particle positions against the grid each step.", _boolean
        ),
        Setting("save_diagnostics", True, "bool", "Write diagnostics.json.", _boolean),
        Setting("save_checkpoint", False, "bool", "Write a checkpoint at each save.", _boolean),
        Setting("continuum_fig_size", None, "inches", "Size of continuum figures.", None),
        Setting("mp_fig_size", None, "inches", "Size of material-point figures.", None),
        Setting(
            "at_start", None, "function", "Callback `(L, P, G)` before the first step.", _optional(_function)
        ),
        Setting(
            "after_every_timestep",
            None,
            "function",
            "Callback `(L, P, G)` after each step.",
            _optional(_function),
        ),
        Setting(
            "after_every_savetime",
            None,
            "function",
            "Callback `(L, P, G)` at each save.",
            _optional(_function),
        ),
        Setting(
            "at_end", None, "function", "Callback `(L, P, G)` after the last step.", _optional(_function)
        ),
    ],
}

# Material fields, by law (constit.validate_material checks their values).
MATERIAL_COMMON = (
    "law",
    "n",
    "X",
    "Y",
    "rho",
    "A",
    "initial_v",
    "phi",
    "domain_size",
    "initial_C",
    "critical_time",
    "inlet_rate",
    "inlet_gap",
    "inlet_y",
)
MATERIAL_FIELDS = {
    "rigid": (),
    "elastic": ("K", "G"),
    "von_mises": ("K", "G", "yield_stress"),
    "viscous": ("K", "eta", "eta_v"),
    "bingham": ("K", "eta", "eta_0", "yield_stress", "eta_v", "gamma_c"),
    "pouliquen": ("K", "rho_s", "mu_s", "mu_d", "I_0", "eta_max", "phi_flow"),
    "linear_mu": ("K", "rho_s", "mu_s", "b", "eta_max"),
    "dp": ("K", "G", "mu_s", "beta"),
    "dp_rate": ("K", "G", "mu_s", "beta", "t_star"),
    "granular": (
        "rho_s",
        "mu_s",
        "phi_max",
        "mu_d",
        "I_0",
        "grain_diameter",
        "phi_min",
        "cidr_deficit",
        "cidr_p_ref",
    ),
}


def _instance_settings(namespace):
    return {name: value for name, value in vars(namespace).items() if not name.startswith("_")}


def _reject_unknown(where, given, known):
    for name in given:
        if name not in known:
            close = difflib.get_close_matches(name, known, n=1)
            hint = f"; did you mean {close[0]!r}?" if close else ""
            raise ValueError(f"unknown setting {where}.{name}{hint}")


def _apply(where, namespace, settings, skip=()):
    given = {name: v for name, v in _instance_settings(namespace).items() if name not in skip}
    _reject_unknown(where, given, [s.name for s in settings])
    for s in settings:
        # A setting defined on the class (a method such as a callback or gravity_ramp) is kept.
        if s.name not in given and not hasattr(type(namespace), s.name):
            if s.default is REQUIRED:
                raise ValueError(f"{where}.{s.name} is required")
            setattr(namespace, s.name, s.default)
        value = getattr(namespace, s.name)
        if s.check is not None:
            message = s.check(value)
            if message:
                raise ValueError(f"{where}.{s.name} {message}")


class _Namespace:
    pass


def normalise(P):
    """Check P's settings, fill defaults and set the internal attributes the solver reads."""
    for where in ("G", "B", "O"):
        if not hasattr(P, where):
            setattr(P, where, _Namespace())
    if "update_forces" in vars(type(P)) or "update_forces" in vars(P):
        raise ValueError(
            "inputs set gravity, gravity_ramp, gravity_angle and B.traction_* instead of update_forces"
        )
    # args: the command line, recorded by initialise.get_parameters.
    _apply("P", P, SETTINGS["P"], skip=("G", "B", "O", "S", "args"))
    for where in ("G", "B", "O"):
        _apply(where, getattr(P, where), SETTINGS[where])
    materials = P.S if isinstance(P.S, list) else [P.S]
    for index, S in enumerate(materials):
        law = getattr(S, "law", None)
        if law not in MATERIAL_FIELDS:
            continue  # constit.get_law reports unknown laws
        _reject_unknown(f"S[{index}] ({law})", _instance_settings(S), MATERIAL_COMMON + MATERIAL_FIELDS[law])
    _translate(P)


def _translate(P):
    G, B, O = P.G, P.B, P.O
    # Run.
    P.CFL, P.M_tol = P.cfl, P.mass_tolerance
    P.pressure = P.initial_stress
    P.l = P.kinetic_pressure_scale
    P.c, P.normalise_phi, P.smooth_gamma_dot = (
        P.segregation_rate,
        P.segregation_normalise,
        P.smooth_shear_rate,
    )
    if P.output_dir is not None:
        P.save_dir = P.output_dir.rstrip("/") + "/"
    P.max_g = P.gravity
    if P.initial_stress == "compression" and P.initial_pressure is None:
        raise ValueError("initial_stress = 'compression' needs initial_pressure")
    # Grid.
    G.x_m, G.x_M, G.y_m, G.y_M = G.x_min, G.x_max, G.y_min, G.y_max
    G.s = np.asarray(G.grain_sizes, dtype=float).ravel()
    G.ns = len(G.s)
    # Boundaries.
    walls = set(B.walls)
    B.has_bottom, B.has_top = "bottom" in walls, "top" in walls
    B.has_left, B.has_right = "left" in walls, "right" in walls
    B.cyclic_lr = B.periodic
    B.roughness = B.wall_contact == "no_slip"
    B.no_slip_bottom = B.wall_contact == "no_slip_bottom"
    if (B.wall_contact == "coulomb") != (B.wall_friction is not None):
        raise ValueError("B.wall_friction is set exactly when B.wall_contact = 'coulomb'")
    B.wall_mu = B.wall_friction or False
    B.separable = B.one_sided
    obstacles = set(B.obstacles)
    B.wall = "centre_wall" in obstacles
    for name in ("two_walls", "box", "silo_left", "silo_bottom"):
        setattr(B, name, name in obstacles)
    B.outlet_left, B.outlet_bottom = "left" in B.outlets, "bottom" in B.outlets
    B.inlet_right, B.inlet_top = "right" in B.inlets, "top" in B.inlets
    B.conveyor = B.conveyor_speed is not None
    P.v_0 = B.conveyor_speed
    B.vertical_force = B.traction_vertical is not None
    B.horizontal_force = B.traction_horizontal is not None
    # Output.
    figures, fields = set(O.figures), set(O.fields)
    O.plot_continuum = "continuum" in figures
    O.plot_material_points = "material_points" in figures
    O.plot_gsd_mp = "gsd_material_points" in figures
    O.plot_gsd_grid = "gsd_grid" in figures
    O.plot_gamma_dot = "shear_rate" in figures
    O.save_s_bar, O.save_density = "s_bar" in fields, "density" in fields
    O.save_u, O.save_phi_MP = "u" in fields, "phi_material_points" in fields
    for name in ("at_start", "after_every_timestep", "after_every_savetime", "at_end"):
        if vars(O).get(name, False) is None:
            delattr(O, name)
    # Forces evaluated at the current time each step (MPM.step calls P.update_forces).
    P.update_forces = lambda: update_forces(P)


def _at(value, t):
    return value(t) if callable(value) else value


def update_forces(P):
    """Gravity, its angle and the boundary tractions at time P.t."""
    P.g = P.gravity * (P.gravity_ramp(P.t) if P.gravity_ramp is not None else 1.0)
    P.theta = _at(P.gravity_angle, P.t)
    P.q_v = _at(P.B.traction_vertical, P.t) if P.B.traction_vertical is not None else 0.0
    P.q_h = _at(P.B.traction_horizontal, P.t) if P.B.traction_horizontal is not None else 0.0


def markdown(where):
    """The settings table of one namespace, for docs/features.md."""
    rows = ["| Setting | Default | Type | Meaning |", "| --- | --- | --- | --- |"]
    for s in SETTINGS[where]:
        default = "required" if s.default is REQUIRED else f"`{s.default!r}`"
        rows.append(f"| `{s.name}` | {default} | {s.kind} | {s.doc} |")
    return "\n".join(rows)


if __name__ == "__main__":
    for where in SETTINGS:
        print(f"<!-- settings:{where} -->")
        print(markdown(where))
        print(f"<!-- /settings:{where} -->\n")
