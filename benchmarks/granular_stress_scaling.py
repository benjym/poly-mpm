"""Frozen first-step collapse problems under grid refinement (CPU, FP64).

Each case builds the first contact problem of a column collapse from rest
(0.5 m wide, 1 m high, 4 particles per cell, dt = 0.002 s) and runs
granular_mixed.solve_velocity at tolerance 1e-4 with a 3,000-sweep budget,
recording the residual history. Cases vary the stress space (nodal, P0,
coarse Q1 with coarsening k), the friction coefficient, and nested-iteration
warm starts from coarser stress grids.

    uv run --locked python -m benchmarks.granular_stress_scaling --output benchmarks/granular_stress_scaling.json
"""

import argparse
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np

TOLERANCE = 1e-4
MAX_ITERATIONS = 300  # 3,000 sweeps


def problem(spacing, basis="nodal", coarsening=2, mu=None):
    from granular import assemble, coarse_shape
    from mpm_gpu.benchmark import columns
    from transfers import nodal_velocity, scatter

    P, G, L = columns(1, spacing=spacing, particle_spacing=spacing / 4)
    P.granular_stress_basis = basis
    P.granular_stress_coarsening = coarsening
    b = L.materials[0]
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    body = np.array([P.g * np.sin(P.theta), P.g * np.cos(P.theta), 0.0])
    scatter(G.fe, b.n, b.N[:, :, None] * b.m[:, None, None] * body)
    G.BCs(P)
    G.update_momentum(P)
    trial = nodal_velocity(G, P)[:, :2].ravel()
    B, free, active, _, gap = assemble(P, G, b, compact=True)
    inverse = np.zeros_like(trial)
    inverse[free] = 1 / np.repeat(G.m, 2)[free]
    record = dict(
        B=B,
        trial=trial,
        inverse=inverse,
        gap=gap,
        mu=b.material.mu if mu is None else mu,
        dt=P.dt,
        particles=len(b),
        free_velocity=int(free.sum()),
        P=P,
        G=G,
        active=active,
    )
    if basis == "coarse":
        nx, ny = coarse_shape(P)
        index = np.arange(nx * ny)
        record["positions"] = np.column_stack(
            (
                P.G.x_m + (index % nx) * coarsening * G.dx,
                P.G.y_m + (index // nx) * coarsening * G.dy,
                np.zeros(nx * ny),
            )
        )[active]
        record["total"] = nx * ny
    return record


def solve(p, initial=None):
    """solve_velocity, returning the last iterate and its statistics even when unconverged."""
    import granular_mixed

    try:
        contact, stats = granular_mixed.solve_velocity(
            p["B"],
            p["trial"],
            p["inverse"],
            p["gap"],
            p["mu"],
            p["dt"],
            tolerance=TOLERANCE,
            max_iterations=MAX_ITERATIONS,
            initial=initial,
        )
        return contact, dict(converged=True, residual=stats["residual"], sweeps=stats["sweeps"])
    except RuntimeError as exc:
        text = str(exc)
        residual = float(text.split("residual=")[1].split(",")[0])
        return None, dict(converged=False, residual=residual, sweeps=10 * MAX_ITERATIONS)


def case(name, spacing, **options):
    started = perf_counter()
    p = problem(spacing, **options)
    _, result = solve(p)
    return dict(
        name=name,
        spacing=spacing,
        particles=p["particles"],
        stress_unknowns=3 * p["B"].count,
        free_velocity=p["free_velocity"],
        mu=float(p["mu"]),
        seconds=perf_counter() - started,
        **{k: v for k, v in options.items() if k != "mu"},
        **result,
    )


def nested(name, spacing, factors):
    """Solve on successively finer coarse stress grids, each warm-started from the last."""
    from granular import coarse_stencil

    started = perf_counter()
    levels, guess, previous = [], None, None
    for k in factors:
        p = problem(spacing, basis="coarse", coarsening=k)
        initial = None
        if guess is not None:
            full = np.zeros((previous["total"], 3))
            full[previous["active"]] = guess
            nodes, weights = coarse_stencil(previous["P"], previous["G"], p["positions"])
            initial = np.einsum("pr,prk->pk", weights, full[nodes])
        guess, result = solve(p, initial)
        levels.append(dict(coarsening=k, stress_unknowns=3 * p["B"].count, **result))
        previous = p
        if guess is None:
            break
    return dict(name=name, spacing=spacing, levels=levels, seconds=perf_counter() - started)


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("granular-stress-scaling.json"))
    args = parser.parse_args(argv)
    cases = [
        case("nodal_h0.05", 0.05),
        case("nodal_h0.025", 0.025),
        case("p0_h0.05", 0.05, basis="P0"),
        case("p0_h0.025", 0.025, basis="P0"),
        case("coarse2_h0.05", 0.05, basis="coarse", coarsening=2),
        case("coarse2_h0.025", 0.025, basis="coarse", coarsening=2),
        case("coarse2_h0.0125", 0.0125, basis="coarse", coarsening=2),
        case("nodal_h0.025_frictionless", 0.025, mu=0.0),
        case("nodal_h0.025_mu0.1", 0.025, mu=0.1),
        case("coarse2_h0.0125_frictionless", 0.0125, basis="coarse", coarsening=2, mu=0.0),
        case("coarse2_h0.0125_mu0.1", 0.0125, basis="coarse", coarsening=2, mu=0.1),
        # Fixed physical stress spacing 0.05 m while the velocity grid refines.
        case("stress0.05_h0.025", 0.025, basis="coarse", coarsening=2),
        case("stress0.05_h0.0125", 0.0125, basis="coarse", coarsening=4),
        case("stress0.05_h0.00625", 0.00625, basis="coarse", coarsening=8),
    ]
    studies = [
        nested("nested_h0.025", 0.025, (4, 2)),
        nested("nested_h0.0125", 0.0125, (8, 4, 2)),
        nested("nested_stress0.05_h0.00625", 0.00625, (16, 8)),
    ]
    report = dict(
        python=platform.python_version(),
        machine=platform.machine(),
        tolerance=TOLERANCE,
        max_sweeps=10 * MAX_ITERATIONS,
        cases=cases,
        nested=studies,
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
