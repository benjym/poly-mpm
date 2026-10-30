"""Matched CPU/Metal timestep benchmark on one machine (docs/gpu.md protocol).

Each case starts both backends from one prepared state with the same
timestep, tolerance, sweep budget and diagnostics, and runs the same number of
steps. A Metal warm-up run compiles every kernel before timing. Paired samples
alternate the CPU/Metal order. Times are host wall time for complete steps,
including every waited submission and the final snapshot; kernel timestamps
are not available in the current backend.

Disconnected cases repeat one collapsing column side by side, so the solver
work per column does not grow with the particle count.

    python -m mpm_gpu.benchmark --output benchmarks/metal_benchmark_m4.json
"""

import argparse
from copy import deepcopy
import json
from pathlib import Path
from time import perf_counter

import numpy as np

TOLERANCE = 1e-4


# Glass beads (Jop et al.), 3 mm grains; tests/test_cidr.py.
CIDR = dict(
    phi_max=0.605,
    mu_s=float(np.tan(np.radians(20.9))),
    mu_d=float(np.tan(np.radians(32.76))),
    I_0=0.279,
    grain_diameter=0.003,
    phi_min=0.05,
    cidr_deficit=0.005,
)


def columns(
    count,
    *,
    spacing=0.05,
    particle_spacing=0.0125,
    height=1.0,
    width=0.5,
    pitch=1.0,
    static=False,
    figures=None,
    solver="velocity",
):
    """count columns of width x height, one per `pitch` of domain, on a grid of `spacing`.

    figures, a directory, turns on the material-point and continuum figures there.
    solver="cidr" uses critical-state CIDR with the glass-bead parameters of
    tests/test_cidr.py, the packing at rest on the cap at twice the base pressure.
    """
    from initialise import prepare
    from inputs.granular_collapse import Params

    P = Params([])
    if figures is not None:
        P.output_dir = str(figures)
        P.O.figures = ("material_points", "continuum")
    P.granular_solver, P.granular_tolerance, P.diagnostics = "velocity", TOLERANCE, True
    P.granular_max_nodes, P.granular_max_iterations, P.dt = 10**8, 2000, 0.002
    P.G.x_max = count * pitch
    P.G.y_max = 1.5 * height
    P.G.nx = round(P.G.x_max / spacing) + 1
    P.G.ny = round(P.G.y_max / spacing) + 1
    s = particle_spacing
    x, y = np.meshgrid(np.arange(s / 2, width, s), np.arange(s / 2, height, s))
    X = np.concatenate([x.ravel() + k * pitch for k in range(count)])
    Y = np.tile(y.ravel(), count)
    P.S.X, P.S.Y, P.S.n, P.S.A = X, Y, X.size, s * s
    if static:
        # Walls on both sides of a single full-width column hold it at rest.
        P.B.walls = set(P.B.walls) | {"left", "right"}
    if solver == "cidr":
        P.granular_solver = "cidr"
        for key, value in CIDR.items():
            setattr(P.S, key, value)
        P.S.cidr_p_ref = 2 * P.S.rho * abs(P.gravity) * height
    return prepare(P)


def _run(state, backend, steps):
    from MPM import run

    P, G, L = state
    P.backend = backend
    started = perf_counter()
    run(P, G, L, max_steps=steps, output=False)
    seconds = perf_counter() - started
    history = [r["granular_solver"] for r in L.diagnostics.history]
    return seconds, history


def case(name, state, *, steps, pairs, rt, backends=("cpu", "metal")):
    """Paired CPU/Metal samples; backends=("metal",) times Metal alone, for cases too large
    for the CPU solver."""
    P0, G0, _ = state
    record = dict(
        name=name,
        particles=len(state[2].materials[0]),
        stress_rows=int(G0.m.size),
        steps=steps,
        dt=P0.dt,
        tolerance=P0.granular_tolerance,
        max_sweeps=10 * P0.granular_max_iterations,
        samples=[],
    )
    # Warm-up: compile every kernel and allocate the runner once.
    warm = deepcopy(state)
    allocated = rt.accounting.allocated_bytes
    started = perf_counter()
    try:
        _run(warm, "metal", 1)
    except RuntimeError as exc:
        # A solve that does not converge is a result to record, not a crash.
        record["metal_warmup_error"] = str(exc)[:200]
    record["metal_warmup_seconds"] = perf_counter() - started
    # Every runner allocates the same workspaces; this is one runner's device memory.
    record["metal_runner_allocated_bytes"] = rt.accounting.allocated_bytes - allocated
    finals = {}
    for pair in range(pairs):
        order = tuple(backends) if pair % 2 == 0 else tuple(reversed(backends))
        sample = dict(pair=pair, order="-".join(order))
        for backend in order:
            run_state = deepcopy(state)
            before = dict(rt.accounting.__dict__)
            try:
                seconds, history = _run(run_state, backend, steps)
                sample[backend] = dict(
                    seconds=seconds,
                    sweeps=[r.get("sweeps") for r in history],
                    newton_attempts=[r.get("newton_attempts", None) for r in history],
                    gmres_iterations=[r.get("gmres_iterations", None) for r in history],
                    newton=[r.get("newton") for r in history],
                    cg_iterations=[r.get("cg_iterations") for r in history],
                    converged=True,
                )
                if backend == "metal":
                    after = rt.accounting.__dict__
                    sample[backend].update(
                        submissions=after["submissions"] - before["submissions"],
                        dispatches=after["dispatches"] - before["dispatches"],
                        readback_bytes=after["readback_bytes"] - before["readback_bytes"],
                        uploaded_bytes=after["uploaded_bytes"] - before["uploaded_bytes"],
                    )
                finals[backend] = run_state[2].materials[0].x.copy()
            except RuntimeError as exc:
                sample[backend] = dict(seconds=None, converged=False, error=str(exc)[:200])
        if sample.get("cpu", {}).get("seconds") and sample["metal"].get("seconds"):
            sample["cpu_over_metal"] = sample["cpu"]["seconds"] / sample["metal"]["seconds"]
        record["samples"].append(sample)
    ratios = [s["cpu_over_metal"] for s in record["samples"] if "cpu_over_metal" in s]
    if ratios:
        record["median_cpu_over_metal"] = float(np.median(ratios))
    if len(finals) == 2:
        record["position_difference_grid_spacings"] = float(
            np.abs(finals["cpu"] - finals["metal"]).max() / min(G0.dx, G0.dy)
        )
    return record


def plotting(rt, count=3, steps=3, repeats=3):
    """Cost of one save (Metal snapshot plus the standard figures) against the steps."""
    import tempfile

    from MPM import run
    from .runner import _runners

    rows = []
    for backend in ("cpu", "metal"):
        with tempfile.TemporaryDirectory() as directory:
            P, G, L = columns(count, figures=directory)
            P.backend = backend
            started = perf_counter()
            run(P, G, L, max_steps=steps, output=False)
            stepping = perf_counter() - started
            snapshot, figures = [], []
            for _ in range(repeats):
                if backend == "metal":
                    started = perf_counter()
                    _runners[L].snapshot(P, G, L)
                    snapshot.append(perf_counter() - started)
                started = perf_counter()
                P.O.after_every_savetime(L, P, G)
                figures.append(perf_counter() - started)
            rows.append(
                dict(
                    backend=backend,
                    particles=len(L.materials[0]),
                    steps=steps,
                    seconds_per_step=stepping / steps,
                    snapshot_seconds=float(np.median(snapshot)) if snapshot else 0.0,
                    figure_seconds=float(np.median(figures)),
                )
            )
    return rows


def main(argv=None):
    from .runner import runtime
    from .runtime import environment

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("metal-benchmark.json"))
    parser.add_argument("--quick", action="store_true", help="Only the small cases")
    parser.add_argument("--cases", nargs="*", help="Run only these case names")
    parser.add_argument("--append", action="store_true", help="Add the cases to an existing report")
    parser.add_argument("--plotting", action="store_true", help="Only measure the cost of one save")
    parser.add_argument("--solver", choices=("velocity", "cidr"), default="velocity")
    parser.add_argument("--steps", type=int, help="Steps per run in every case (default: per case)")
    parser.add_argument("--pairs", type=int, help="Paired samples per case (default: per case)")
    args = parser.parse_args(argv)
    solver = args.solver
    rt = runtime()
    if args.plotting:
        report = json.loads(args.output.read_text()) if args.output.exists() else {}
        report["plotting"] = plotting(rt)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report["plotting"], indent=2), flush=True)
        return 0
    report = dict(
        schema_version=1,
        solver=solver,
        environment=rt.environment,
        host=environment(),
        protocol=(
            "Same prepared state, dt, tolerance and sweep budget on both backends; one Metal warm-up "
            "run before timing; paired samples alternate CPU/Metal order; host wall time per complete "
            "run including waited submissions and the final snapshot. Not every case has five pairs: "
            "the CPU run time of the larger cases limits the number of pairs, recorded per case."
        ),
        cases=[],
    )
    plan = [
        ("static_column_3200", lambda: columns(1, width=1.0, static=True, solver=solver), 10, 3),
        ("collapse_3200", lambda: columns(1, solver=solver), 3, 3),
        ("collapse_disconnected_9600", lambda: columns(3, solver=solver), 3, 2),
    ]
    if not args.quick:
        plan += [
            (
                "collapse_connected_10000",
                lambda: columns(1, spacing=0.02, particle_spacing=0.01, width=1.0, pitch=2.0, solver=solver),
                2,
                1,
            ),
            ("collapse_disconnected_99200", lambda: columns(31, solver=solver), 2, 1),
        ]
        if solver == "cidr":
            # The connected collapse refined 3.2 times each way, on Metal alone.
            plan.append(
                (
                    "collapse_connected_102400",
                    lambda: columns(
                        1, spacing=0.00625, particle_spacing=0.003125, width=1.0, pitch=2.0, solver=solver
                    ),
                    5,
                    1,
                    ("metal",),
                )
            )
    plan = [
        (
            entry[0],
            entry[1],
            args.steps or entry[2],
            args.pairs or entry[3],
            entry[4] if len(entry) > 4 else None,
        )
        for entry in plan
    ]
    if args.cases:
        plan = [entry for entry in plan if entry[0] in args.cases]
    if args.append and args.output.exists():
        previous = json.loads(args.output.read_text())
        report["cases"] = [c for c in previous["cases"] if c["name"] not in {entry[0] for entry in plan}]
    for name, build, steps, pairs, backends in plan:
        print(name, flush=True)
        options = dict(backends=backends) if backends else {}
        report["cases"].append(case(name, build(), steps=steps, pairs=pairs, rt=rt, **options))
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
