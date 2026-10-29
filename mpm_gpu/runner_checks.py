"""Runner and science gate: complete Metal timesteps against CPU runs (docs/gpu.md).

Every case runs the same prepared state on both backends. The CPU uses the
requested Metal tolerance unless a case states otherwise, so differences
measure FP32 and the device solver rather than a tolerance mismatch.
"""

from copy import deepcopy

import numpy as np

from .mechanics_checks import _require

TOLERANCE = 1e-4
FROZEN_RESIDUAL = 2e-4
MASS_ERROR = 1e-5
STATIC_DISPLACEMENT = 1e-3  # grid spacings
INCLINE_SPEED = 1e-2  # relative
COLLAPSE_POSITION = 1e-2  # grid spacings
OVERPACKING = 1e-3
PENDING = ()


def _column(**options):
    from initialise import prepare
    from inputs.granular_column import Params

    P = Params([])
    P.granular_solver, P.granular_tolerance, P.diagnostics = "velocity", TOLERANCE, True
    for key, value in options.items():
        setattr(P, key, value)
    return prepare(P)


def _collapse(**options):
    from initialise import prepare
    from inputs.granular_collapse import Params

    P = Params([])
    P.granular_solver, P.granular_tolerance, P.diagnostics = "velocity", TOLERANCE, True
    for key, value in options.items():
        setattr(P, key, value)
    return prepare(P)


def _pair(state):
    """Independent CPU and Metal copies of one prepared state."""
    cpu = deepcopy(state)
    metal = deepcopy(state)
    metal[0].backend = "metal"
    return cpu, metal


def _run(state, **options):
    from MPM import run

    run(*state, output=False, **options)
    return state


def _spacing(G):
    return min(G.dx, G.dy)


def frozen_residual_check(rt, kind):
    """One Metal step from a CPU state; FP64 frozen natural residual of its contact."""
    from granular import assemble
    from MPM import step
    from transfers import nodal_velocity, scatter

    if kind == "static":
        state = _column(dt=0.02)
    elif kind == "collapse":
        state = _run(_collapse(dt=0.005), max_steps=20)
    else:
        state = _incline(0.75)
        state = _run(state, max_steps=5)
    P, G, L = deepcopy(state)
    # FP64 operator, trial velocity and gap of the start-of-step state.
    b = L.materials[0]
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    body = np.array([P.g * np.sin(P.theta), P.g * np.cos(P.theta), 0.0])
    scatter(G.fe, b.n, b.N[:, :, None] * b.m[:, None, None] * body)
    if P.B.cyclic_lr:
        G.make_cyclic(P, G, ["fe"])
    G.BCs(P)
    G.update_momentum(P)
    trial = nodal_velocity(G, P)[:, :2]
    B, free, active, _, gap = assemble(P, G, b, compact=True)
    inverse = np.zeros(2 * len(G.m))
    inverse[free] = 1 / np.repeat(G.m, 2)[free]
    inverse = inverse.reshape(-1, 2)
    P, G, L = state
    P.backend = "metal"
    step(P, G, L)
    contact = G.granular_contact[active]
    q = (B @ trial.ravel()).reshape(-1, 3) + gap
    magnitude = max(float(np.max(-q[:, 0], initial=0)), float(np.max(np.abs(q[:, 1:]), initial=0)))
    from .mechanics_checks import _reference_metric

    _, norm, _ = _reference_metric(B, inverse, P.dt)
    u = trial + P.dt * inverse * B.transpose(contact.ravel()).reshape(-1, 2)
    if P.B.cyclic_lr:
        u[P.G.nx - 1 :: P.G.nx] = u[:: P.G.nx]
    z = contact / (magnitude / norm)
    w = ((B @ u.ravel()).reshape(z.shape) + gap) / magnitude
    argument = z - w
    radius = b.material.mu * np.maximum(z[:, 0], 0)
    length = np.linalg.norm(argument[:, 1:], axis=1)
    projected = argument.copy()
    projected[:, 0] = np.maximum(argument[:, 0], 0)
    projected[:, 1:] *= np.minimum(1, np.divide(radius, length, out=np.ones(len(z)), where=length > 0))[
        :, None
    ]
    residual = float(np.max(np.linalg.norm(z - projected, axis=1), initial=0))
    _require(residual <= FROZEN_RESIDUAL, f"{kind}: FP64 frozen natural residual {residual:.3g}")
    velocity = nodal_velocity(G, P)[:, :2]
    velocity_error = float(np.abs(velocity - u).max() / max(np.abs(trial).max(), np.abs(u).max()))
    _require(velocity_error <= 1e-5, f"{kind}: final velocity differs by {velocity_error:.3g}")
    return dict(
        kind=kind,
        active_rows=len(active),
        frozen_residual=residual,
        device_residual=G.granular_stats["residual"],
        velocity_error=velocity_error,
    )


def _incline(slope, **options):
    state = _column(theta=np.arctan(slope), t_f=0.05, **options)
    P, G, L = state
    P.B.cyclic_lr = P.B.roughness = True
    G.boundary_v[:] = False
    return state


def static_check(rt):
    cpu, metal = _pair(_column(dt=0.02, t_f=0.2))
    initial = metal[2].materials[0].x.copy()
    _run(cpu)
    _run(metal)
    P, G, L = metal
    b = L.materials[0]
    displacement = float(np.abs(b.x - initial).max() / _spacing(G))
    last = L.diagnostics.last
    mass_error = abs(L.diagnostics.last["mass_error"]) / last["mass"]
    device_mass = abs(metal_runner(L).last["device_mass"] - last["mass"]) / last["mass"]
    _require(displacement <= STATIC_DISPLACEMENT, f"static displacement {displacement:.3g} grid spacings")
    _require(
        mass_error <= MASS_ERROR and device_mass <= MASS_ERROR,
        f"mass error {mass_error:.3g}, {device_mass:.3g}",
    )
    _require(b.pressure.min() < -100 and b.pressure.max() <= 1e-6 * abs(b.pressure.min()), "static pressure")
    return dict(
        steps=P.tstep,
        displacement=displacement,
        max_speed=float(np.abs(b.v).max()),
        cpu_max_speed=float(np.abs(cpu[2].materials[0].v).max()),
        relative_mass_error=mass_error,
        device_mass_error=device_mass,
        **_budgets(cpu, metal),
    )


def metal_runner(L):
    from .runner import _runners

    return _runners[L]


def _budgets(cpu, metal):
    """Dimensionless budget residuals on both backends, scaled by the external impulse/work."""
    c, m = cpu[2].diagnostics.last, metal[2].diagnostics.last
    impulse = max(np.linalg.norm(c["external_impulse"]), 1e-30)
    work = max(abs(c["external_work"]), abs(c["constraint_work"]), 1e-30)
    return dict(
        momentum_residual_cpu=float(np.linalg.norm(c["momentum_residual"]) / impulse),
        momentum_residual_metal=float(np.linalg.norm(m["momentum_residual"]) / impulse),
        energy_residual_cpu=float(c["energy_residual"] / work),
        energy_residual_metal=float(m["energy_residual"] / work),
    )


def freefall_check(rt):
    """Low packing without walls: stress free and exactly ballistic."""
    from MPM import step

    P, G, L = _column(dt=0.005, max_g=-2.0)
    P.backend = "metal"
    P.B.has_bottom = P.B.has_left = P.B.has_right = False
    G.boundary_h[:] = False
    G.boundary_v[:] = False
    b = L.materials[0]
    b.m *= 0.2
    b.rho *= 0.2
    b.initial_rho *= 0.2
    b.v[:, 0] = -0.01 * (b.x[:, 0] - 0.5)
    before = b.v[:, 1].copy()
    step(P, G, L)
    change = float(np.abs(b.v[:, 1] - before + 2 * P.dt).max())
    _require(np.abs(b.pressure).max() <= 1e-8, f"free-fall pressure {np.abs(b.pressure).max():.3g}")
    _require(change <= 1e-6 * 2 * P.dt, f"free-fall velocity change error {change:.3g}")
    _require(np.any(b.L[:, 0, 0] < 0), "affine contraction lost")
    _require(G.granular_stats["maximum_packing"] < 0.2, "packing")
    return dict(velocity_change_error=change, maximum_packing=G.granular_stats["maximum_packing"])


def lift_check(rt, kind="separable"):
    """One-sided or frictional walls: a column under upward gravity leaves the floor ballistically."""
    cpu, metal = _pair(_column(dt=0.005, t_f=0.05, max_g=2.0))
    for P, G, L in (cpu, metal):
        if kind == "separable":
            P.B.separable = True
        else:
            P.B.wall_mu = 0.3
    before = metal[2].materials[0].v[:, 1].copy()
    _run(cpu)
    _run(metal)
    (Pc, Gc, Lc), (P, G, L) = cpu, metal
    b, bc = L.materials[0], Lc.materials[0]
    change = float(np.abs(b.v[:, 1] - before - 2.0 * P.t).max())
    cpu_change = float(np.abs(bc.v[:, 1] - before - 2.0 * Pc.t).max())
    _require(np.abs(b.pressure).max() <= 1e-8, f"lift pressure {np.abs(b.pressure).max():.3g}")
    _require(change <= 1e-5 * 2.0 * P.t, f"lift velocity change error {change:.3g}")
    position = float(np.abs(b.x - bc.x).max() / _spacing(G))
    _require(position <= COLLAPSE_POSITION, f"lift positions differ by {position:.3g} grid spacings")
    return dict(
        kind=kind,
        steps=P.tstep,
        velocity_change_error=change,
        cpu_velocity_change_error=cpu_change,
        position_difference=position,
    )


def wall_slide_check(rt, wall_mu):
    """Frictional floor: a periodic layer on tan(theta) = 0.4 with internal mu = 0.5 slides
    as a block at g (sin - wall_mu cos) t when wall_mu < 0.4 and otherwise rests."""
    state = _column(dt=0.005, t_f=0.05, theta=np.arctan(0.4))
    P, G, L = state
    P.B.wall_mu, P.B.cyclic_lr = wall_mu, True
    P.B.has_left = P.B.has_right = False
    G.boundary_v[:] = False
    cpu, metal = _pair(state)
    _run(cpu)
    _run(metal)
    (Pc, Gc, Lc), (P, G, L) = cpu, metal
    speed, cpu_speed = -L.materials[0].v[:, 0], -Lc.materials[0].v[:, 0]
    result = dict(wall_mu=wall_mu, steps=P.tstep)
    if wall_mu < 0.4:
        expected = -P.g * (np.sin(P.theta) - wall_mu * np.cos(P.theta)) * P.t
        error = float(np.abs(speed - expected).max() / expected)
        cpu_error = float(np.abs(cpu_speed - expected).max() / expected)
        _require(error <= INCLINE_SPEED, f"wall sliding speed error {error:.3g}")
        result.update(expected=float(expected), speed_error=error, cpu_speed_error=cpu_error)
    else:
        moved = float(np.abs(L.materials[0].v).max())
        _require(moved <= 1e-3 * _spacing(G) / P.t_f, f"layer on a rough floor moved {moved:.3g}")
        result.update(max_speed=moved, cpu_max_speed=float(np.abs(Lc.materials[0].v).max()))
    result.update(_budgets(cpu, metal))
    return result


def _silo(inlet):
    """A box draining through a silo_bottom opening into outlet_bottom, optionally fed at the top."""
    from types import SimpleNamespace as NS

    from initialise import prepare

    spacing = 0.05
    x, y = np.meshgrid(np.arange(0.025, 1.0, spacing), np.arange(0.025, 0.6, spacing))
    P = NS(
        solver="granular",
        granular_solver="velocity",
        granular_tolerance=TOLERANCE,
        diagnostics=True,
        dt=0.005,
        t_f=0.3,
        savetime=0.1,
        max_g=-9.81,
        save_dir="im/metal_silo/",
        # Ten nodes keep the top inlet at x = 0.5 off a node line.
        G=NS(nx=10, ny=11, x_m=0.0, x_M=1.0, y_m=0.0, y_M=1.0),
        B=NS(has_left=True, has_right=True, silo_bottom=True, outlet_bottom=True, inlet_top=inlet),
        O=NS(),
    )
    P.S = NS(
        law="granular",
        n=x.size,
        X=x.ravel(),
        Y=y.ravel(),
        A=spacing**2,
        rho=1500.0,
        rho_s=2500.0,
        packing_max=0.6,
        mu=0.5,
    )
    if inlet:
        P.S.inlet_rate, P.S.gap = 0.02, (0.025, 0.025)
    return prepare(P)


def boundary_events_check(rt, inlet):
    """Silo outlet (and inlet): host boundary events keep the particle set equal to the CPU's."""
    state = _silo(inlet)
    start = len(state[2].materials[0])
    cpu, metal = _pair(state)
    _run(cpu)
    _run(metal)
    bc, b = cpu[2].materials[0], metal[2].materials[0]
    _require(len(b) == len(bc), f"particle counts differ: {len(b)} and {len(bc)}")
    position = float(np.abs(b.x - bc.x).max() / _spacing(metal[1]))
    mass_error = abs(metal[2].diagnostics.last["mass_error"]) / metal[2].diagnostics.last["mass"]
    _require(position <= COLLAPSE_POSITION, f"silo positions differ by {position:.3g} grid spacings")
    _require(mass_error <= MASS_ERROR, f"mass error with boundary flux {mass_error:.3g}")
    return dict(
        inlet=inlet,
        steps=metal[0].tstep,
        particles_start=start,
        particles_end=len(b),
        position_difference=position,
        relative_mass_error=mass_error,
    )


def incline_check(rt, slope):
    cpu, metal = _pair(_incline(slope))
    _run(cpu)
    _run(metal)
    P, G, L = metal
    b = L.materials[0]
    speed = -b.v[:, 0]
    result = dict(slope=slope, steps=P.tstep)
    if slope < b.material.mu:
        _require(
            np.abs(b.v).max() <= 1e-3 * _spacing(G) / P.t_f, f"static incline moved {np.abs(b.v).max():.3g}"
        )
        result["max_speed"] = float(np.abs(b.v).max())
    else:
        expected = -P.g * (np.sin(P.theta) - b.material.mu * np.cos(P.theta)) * P.t
        error = abs(speed.max() - expected) / expected
        cpu_error = abs((-cpu[2].materials[0].v[:, 0]).max() - expected) / expected
        _require(error <= INCLINE_SPEED, f"incline speed error {error:.3g}")
        result.update(
            speed=float(speed.max()),
            expected=float(expected),
            speed_error=float(error),
            cpu_speed_error=float(cpu_error),
        )
    result.update(_budgets(cpu, metal))
    return result


def collapse_check(rt, transfer="flip", pic_fraction=0.0, *, dynamic=False, separable=False, wall_mu=0.0):
    """Column collapse: positions, runout and overpacking against the CPU at the same tolerance."""
    options = dict(dt=0.005, transfer=transfer, pic_fraction=pic_fraction)
    if dynamic:
        options.update(time_stepping="dynamic", CFL=0.1)
    state = _collapse(**options)
    state[0].B.separable = separable
    state[0].B.wall_mu = wall_mu
    cpu, metal = _pair(state)
    initial = cpu[2].materials[0].x[:, 0].max()
    _run(cpu)
    _run(metal)
    (Pc, Gc, Lc), (P, G, L) = cpu, metal
    bc, b = Lc.materials[0], L.materials[0]
    h = _spacing(G)
    position = float(np.abs(b.x - bc.x).max() / h)
    runout = float(abs(b.x[:, 0].max() - bc.x[:, 0].max()) / h)
    over_cpu = max(r["granular_solver"]["integrated_overpacking_fraction"] for r in Lc.diagnostics.history)
    over = max(r["granular_solver"]["integrated_overpacking_fraction"] for r in L.diagnostics.history)
    _require(P.tstep == Pc.tstep, f"step counts differ: {P.tstep} and {Pc.tstep}")
    _require(abs(P.t - Pc.t) <= 1e-12 * max(P.t, 1), "end times differ")
    _require(position <= COLLAPSE_POSITION, f"collapse positions differ by {position:.3g} grid spacings")
    _require(runout <= COLLAPSE_POSITION, f"runout differs by {runout:.3g} grid spacings")
    _require(over - over_cpu <= OVERPACKING, f"additional overpacking {over - over_cpu:.3g}")
    return dict(
        transfer=transfer,
        pic_fraction=pic_fraction,
        dynamic=dynamic,
        separable=separable,
        wall_mu=wall_mu,
        steps=P.tstep,
        position_difference=position,
        runout_difference=runout,
        runout=float(b.x[:, 0].max() - initial),
        runout_cpu=float(bc.x[:, 0].max() - initial),
        velocity_difference=float(np.abs(b.v - bc.v).max() / max(np.abs(bc.v).max(), 1e-30)),
        overpacking=float(over),
        overpacking_cpu=float(over_cpu),
        plastic_work=float(b.plastic_work.sum()),
        plastic_work_cpu=float(bc.plastic_work.sum()),
        **_budgets(cpu, metal),
    )


def scale_check(rt, spacing=0.05, particle_spacing=0.0125, steps=3):
    """A finer collapse (3,200 particles) against the CPU, with the wall time of each.

    Both solvers need more than the default 2,000 sweeps on the first step of
    this state, so the sweep budget is raised for both. The times are single
    warm samples, not a benchmark.
    """
    from time import perf_counter

    from initialise import prepare
    from inputs.granular_collapse import Params

    P = Params([])
    P.granular_solver, P.granular_tolerance, P.diagnostics = "velocity", TOLERANCE, True
    P.granular_max_nodes, P.granular_max_iterations, P.dt = 10**7, 2000, 0.002
    P.G.nx = round((P.G.x_M - P.G.x_m) / spacing) + 1
    P.G.ny = round((P.G.y_M - P.G.y_m) / spacing) + 1
    x, y = np.meshgrid(
        np.arange(particle_spacing / 2, 0.5, particle_spacing),
        np.arange(particle_spacing / 2, 1.0, particle_spacing),
    )
    P.S.X, P.S.Y, P.S.n, P.S.A = x.ravel(), y.ravel(), x.size, particle_spacing**2
    cpu, metal = _pair(prepare(P))
    seconds = {}
    for name, state in (("cpu", cpu), ("metal", metal)):
        started = perf_counter()
        _run(state, max_steps=steps)
        seconds[name] = perf_counter() - started
    G, b, bc = metal[1], metal[2].materials[0], cpu[2].materials[0]
    position = float(np.abs(b.x - bc.x).max() / _spacing(G))
    _require(position <= COLLAPSE_POSITION, f"positions differ by {position:.3g} grid spacings")
    history = [r["granular_solver"] for r in metal[2].diagnostics.history]
    return dict(
        particles=len(b),
        stress_rows=int(G.m.size),
        steps=steps,
        position_difference=position,
        velocity_difference=float(np.abs(b.v - bc.v).max() / np.abs(bc.v).max()),
        sweeps=[r["sweeps"] for r in history],
        cpu_sweeps=[r["granular_solver"]["sweeps"] for r in cpu[2].diagnostics.history],
        newton_attempts=[r["newton_attempts"] for r in history],
        cpu_seconds=seconds["cpu"],
        metal_seconds=seconds["metal"],
        timing="single warm samples including the first-step cold solve; not a benchmark",
    )


def reproducibility_check(rt, steps=20):
    """Two Metal runs of one state are bit-identical (deterministic binning and reductions)."""
    runs = []
    for _ in range(2):
        P, G, L = _collapse(dt=0.005)
        P.backend = "metal"
        _run((P, G, L), max_steps=steps)
        b = L.materials[0]
        history = [r["granular_solver"]["sweeps"] for r in L.diagnostics.history]
        runs.append((b.x.copy(), b.v.copy(), b.stress.copy(), G.granular_contact.copy(), history))
    for k, name in enumerate(("positions", "velocities", "stresses", "contact")):
        _require(np.array_equal(runs[0][k], runs[1][k]), f"{name} differ between identical runs")
    _require(runs[0][4] == runs[1][4], "sweep counts differ between identical runs")
    return dict(steps=steps, sweeps=runs[0][4])


SEGREGATION_FRACTION = 1e-2  # Metal-CPU fraction difference, relative to the CPU change
SEGREGATION_SPECIES = 1e-6  # relative species-mass error (FP32 fractions)


def segregation_check(rt, kind, sizes=(0.001, 0.003), normalise=False):
    """Bidisperse or polydisperse segregation against the CPU: fractions, species and grid fields."""
    from initialise import prepare

    ns = len(sizes)
    if kind == "incline":
        from inputs.granular_column import Params

        P = Params([])
        P.theta, P.t_f, c = np.arctan(0.75), 0.2, 1e-2
    else:
        from inputs.granular_collapse import Params

        P = Params([])
        P.dt, P.t_f, c = 0.005, 0.15, 1e-1
    P.granular_solver, P.granular_tolerance, P.diagnostics = "velocity", TOLERANCE, True
    P.G.ns, P.G.s = ns, list(sizes)
    P.S.PHI = np.tile(np.full(ns, 1.0 / ns), (P.S.n, 1))
    P.segregate, P.c, P.normalise_phi = True, c, normalise
    state = prepare(P)
    P, G, L = state
    if kind == "incline":
        P.B.cyclic_lr = P.B.roughness = True
        G.boundary_v[:] = False
    b = L.materials[0]
    cpu, metal = _pair(state)
    initial = (b.m[:, None] * b.phi).sum(axis=0)
    _run(cpu)
    _run(metal)
    bc, bm = cpu[2].materials[0], metal[2].materials[0]
    change = float(np.abs(bc.phi - b.phi).max())
    difference = float(np.abs(bm.phi - bc.phi).max())
    species = float(np.abs((bm.m[:, None] * bm.phi).sum(axis=0) - initial).max() / initial.max())
    unity = float(np.abs(bm.phi.sum(axis=1) - 1).max())
    h = _spacing(metal[1])
    position = float(np.abs(bm.x - bc.x).max() / h)
    _require(change > 1e-3, f"segregation too weak to test ({change:.3g})")
    _require(
        difference <= SEGREGATION_FRACTION * change,
        f"fractions differ by {difference:.3g} (change {change:.3g})",
    )
    _require(species <= SEGREGATION_SPECIES, f"species mass error {species:.3g}")
    _require(unity <= 1e-5, f"fractions sum to one within {unity:.3g}")
    _require(position <= COLLAPSE_POSITION, f"positions differ by {position:.3g} grid spacings")
    fields = {}
    for name in ("phi", "dphi", "pk", "shear_rate"):
        expected, actual = getattr(cpu[1], name), getattr(metal[1], name)
        scale = np.abs(expected).max()
        fields[name] = float(np.abs(actual - expected).max() / scale) if scale else 0.0
    _require(max(fields.values()) <= 1e-2, f"grid fields differ: {fields}")
    limits = {
        name: abs(metal[0].timestep_limits[name] - value) / value
        for name, value in cpu[0].timestep_limits.items()
        if name.startswith("segregation")
    }
    _require(limits and max(limits.values()) <= 1e-2, f"segregation limits differ: {limits}")
    return dict(
        kind=kind,
        species_count=ns,
        normalise_phi=normalise,
        steps=metal[0].tstep,
        fraction_change=change,
        fraction_difference=difference,
        species_error=species,
        fraction_sum_error=unity,
        position_difference=position,
        grid_field_differences=fields,
        limit_differences=limits,
    )


def output_check(rt, directory):
    """Save times, snapshot freshness, plots, diagnostics file and callback rejection."""
    from pathlib import Path

    from MPM import run

    cpu, metal = _pair(_collapse(dt=0.01, t_f=0.1, savetime=0.03))
    records = {}
    for name, state in (("cpu", cpu), ("metal", metal)):
        P, G, L = state
        P.save_dir = str(Path(directory) / name) + "/"
        seen = records[name] = []
        default = P.O.after_every_savetime

        def save(L, P, G, seen=seen, default=default):
            seen.append((P.t, L.materials[0].x.copy(), G.m.copy()))
            default(L, P, G)

        P.O.after_every_savetime = save
        run(P, G, L, output=True)
    times = [t for t, _, _ in records["metal"]]
    _require(times == [t for t, _, _ in records["cpu"]], f"save times differ: {times}")
    h = _spacing(metal[1])
    for (_, x, m), (_, xc, mc) in zip(records["metal"], records["cpu"]):
        _require(np.abs(x - xc).max() <= COLLAPSE_POSITION * h, "stale or wrong snapshot at a save time")
        _require(np.abs(m - mc).max() <= 1e-5 * mc.max(), "stale grid snapshot at a save time")
    files = sorted(
        p.relative_to(Path(directory) / "metal").as_posix()
        for p in (Path(directory) / "metal").rglob("*")
        if p.is_file()
    )
    cpu_files = sorted(
        p.relative_to(Path(directory) / "cpu").as_posix()
        for p in (Path(directory) / "cpu").rglob("*")
        if p.is_file()
    )
    _require(
        files == cpu_files and "data/diagnostics.json" in files,
        f"output files differ: {files} vs {cpu_files}",
    )
    # Per-step callbacks would observe stale arrays and are refused.
    P, G, L = _collapse(dt=0.01, t_f=0.02)
    P.backend = "metal"
    P.save_dir = str(Path(directory) / "rejected") + "/"
    P.O.after_every_timestep = lambda L, P, G: None
    try:
        run(P, G, L, output=True)
    except ValueError as exc:
        rejected = str(exc)
    else:
        raise AssertionError("custom after_every_timestep was accepted")
    return dict(save_times=times, files=len(files), rejected=rejected[:80])


def cases(rt):
    import tempfile

    directory = tempfile.mkdtemp(prefix="metal-runner-")
    return [
        ("frozen_residual_static", lambda: frozen_residual_check(rt, "static")),
        ("frozen_residual_collapse", lambda: frozen_residual_check(rt, "collapse")),
        ("frozen_residual_incline", lambda: frozen_residual_check(rt, "incline")),
        ("static_column", lambda: static_check(rt)),
        ("freefall", lambda: freefall_check(rt)),
        ("lift_separable_walls", lambda: lift_check(rt)),
        ("lift_frictional_walls", lambda: lift_check(rt, "frictional")),
        ("wall_friction_sliding", lambda: wall_slide_check(rt, 0.2)),
        ("wall_friction_resting", lambda: wall_slide_check(rt, 0.6)),
        ("incline_static", lambda: incline_check(rt, 0.25)),
        ("incline_yielding", lambda: incline_check(rt, 0.75)),
        ("collapse_flip", lambda: collapse_check(rt)),
        ("collapse_flip_blend", lambda: collapse_check(rt, "flip", 0.3)),
        ("collapse_pic", lambda: collapse_check(rt, "pic")),
        ("collapse_apic_gradient", lambda: collapse_check(rt, "apic_gradient")),
        ("collapse_dynamic_cfl", lambda: collapse_check(rt, dynamic=True)),
        ("collapse_separable_walls", lambda: collapse_check(rt, separable=True)),
        ("collapse_frictional_walls", lambda: collapse_check(rt, wall_mu=0.3)),
        ("silo_outlet", lambda: boundary_events_check(rt, False)),
        ("silo_outlet_and_inlet", lambda: boundary_events_check(rt, True)),
        ("segregation_incline", lambda: segregation_check(rt, "incline")),
        ("segregation_collapse_walls", lambda: segregation_check(rt, "collapse")),
        (
            "segregation_three_sizes_normalised",
            lambda: segregation_check(rt, "incline", sizes=(0.001, 0.002, 0.004), normalise=True),
        ),
        ("output_scheduling", lambda: output_check(rt, directory)),
        ("reproducibility", lambda: reproducibility_check(rt)),
        ("collapse_3200_particles", lambda: scale_check(rt)),
    ]
