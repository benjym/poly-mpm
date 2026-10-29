"""How far do increment-stopped solves (Daviet 2026 style) drift from certified ones?

Evolving column collapse from rest (0.5 m wide, 1 m high, 4 particles per
cell), dt = 0.002 s, 150 steps (0.3 s), nodal stresses, compact velocity solver.
Each variant changes only the stopping rule and sweep budget. Positions are
compared with the first variant listed for each grid: the certified solve at
residual 1e-4 where it is reachable, otherwise the tighter increment rule.
A certified 1e-6 reference was tried and was not reachable (the first step
ended at 1.2e-5 after 50,000 sweeps). Every step's true natural residual is
recorded.

    uv run --locked python -m benchmarks.granular_stopping_study --output benchmarks/granular_stopping_study.json
"""

import argparse
from copy import deepcopy
import json
from pathlib import Path
from time import perf_counter

import numpy as np

STEPS = 150


def run_variant(state, *, stopping, tolerance, increment=1e-3, max_iterations):
    from MPM import run

    P, G, L = deepcopy(state)
    P.granular_stopping, P.granular_tolerance = stopping, tolerance
    P.granular_increment, P.granular_max_iterations = increment, max_iterations
    started = perf_counter()
    try:
        run(P, G, L, max_steps=STEPS, output=False)
        error = None
    except RuntimeError as exc:
        error = str(exc)[:200]
    seconds = perf_counter() - started
    history = [r["granular_solver"] for r in L.diagnostics.history]
    b = L.materials[0]
    return dict(
        stopping=stopping,
        tolerance=tolerance,
        increment=increment if stopping == "increment" else None,
        max_sweeps=10 * max_iterations,
        completed_steps=P.tstep,
        error=error,
        seconds=seconds,
        sweeps=[s["sweeps"] for s in history],
        residuals=[s["residual"] for s in history],
        certified_steps=sum(bool(s["converged"]) for s in history),
        x=b.x[:, :2].copy(),
        runout=float(b.x[:, 0].max()),
        maximum_packing=max((s["maximum_packing"] for s in history), default=None),
    ), (P, G, L)


def study(spacing, variants, reference):
    from mpm_gpu.benchmark import columns

    state = columns(1, spacing=spacing, particle_spacing=spacing / 4)
    state[0].granular_max_nodes = 10**8
    initial_runout = float(state[2].materials[0].x[:, 0].max())
    rows = []
    base = None
    for name, options in ([("reference", reference)] if reference else []) + variants:
        print(spacing, name, flush=True)
        row, _ = run_variant(state, **options)
        row["name"] = name
        if name == "reference":
            base = row
        rows.append(row)
    for row in rows:
        if base is not None and row is not base and row["completed_steps"] == base["completed_steps"]:
            row["position_difference_grid_spacings"] = float(np.abs(row["x"] - base["x"]).max() / spacing)
            row["runout_difference_grid_spacings"] = abs(row["runout"] - base["runout"]) / spacing
    for row in rows:
        row["runout_from_start"] = row["runout"] - initial_runout
        row["max_residual"] = max(row["residuals"], default=None)
        row["total_sweeps"] = int(sum(row["sweeps"]))
        del row["x"]
    return dict(
        spacing=spacing, particles=len(state[2].materials[0]), steps=STEPS, dt=state[0].dt, variants=rows
    )


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("granular-stopping-study.json"))
    args = parser.parse_args(argv)
    paper = dict(stopping="increment", tolerance=1e-4, increment=1e-3, max_iterations=25)
    tight = dict(stopping="increment", tolerance=1e-4, increment=1e-4, max_iterations=250)
    report = dict(
        studies=[
            study(
                0.05,
                [("increment_1e-3_250_sweeps", paper), ("increment_1e-4_2500_sweeps", tight)],
                dict(stopping="residual", tolerance=1e-4, max_iterations=2000),
            ),
            study(0.025, [("increment_1e-3_250_sweeps", paper)], tight),
        ]
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
