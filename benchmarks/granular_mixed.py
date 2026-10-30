"""Matched CPU comparison of the hybrid and compact velocity contact solvers.

Run python -m benchmarks.granular_mixed --output /tmp/granular-mixed.json.
Setup, plotting and progress are excluded; diagnostics and timestep sequences
are identical. Frozen last-step problems separately check the original residual.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
from time import perf_counter
from unittest.mock import patch

import numpy as np
import scipy

import granular
from granular_mixed import solve_velocity
from benchmarks.granular_prototype import configuration
from initialise import prepare
from MPM import run


def parameters(name, method):
    P = configuration(
        "static" if name in ("static", "incline_static", "incline_sliding") else "collapse",
        0.02 if name == "static" else 0.0025 if name == "refined" else 0.005,
        grid_spacing=0.125 if name == "refined" else 0.25,
        solver=method,
    )
    if name == "refined":
        P.t_f = 0.05
    elif name.startswith("incline"):
        P.B.walls, P.B.periodic, P.B.wall_contact = {"bottom"}, True, "no_slip"
        P.gravity_angle = np.arctan(0.25 if name == "incline_static" else 0.75)
        P.t_f = 0.05
    elif name == "plane_strain":
        P.granular_stress = "plane_strain"
        P.transfer = "apic_gradient"
        P.granular_viscosity = 0.001
    return P


def natural_residual(W, q, reaction, mu):
    """Reference normalization and infinity norm, not the new norm bound."""
    norm = float(np.max(np.asarray(abs(W).sum(axis=1)), initial=0))
    magnitude = max(np.max(-q[:, 0], initial=0), np.max(np.abs(q[:, 1:]), initial=0), 1e-30)
    z = reaction * norm / magnitude
    w = (W @ reaction.ravel()).reshape(q.shape) / magnitude + q / magnitude
    projected = z - w
    projected[:, 0] = np.maximum(projected[:, 0], 0)
    radius = mu * np.maximum(z[:, 0], 0)
    length = np.linalg.norm(projected[:, 1:], axis=1)
    projected[:, 1:] *= np.minimum(1, np.divide(radius, length, out=np.ones(len(q)), where=length > 0))[
        :, None
    ]
    return float(np.max(np.abs(z - projected), initial=0))


def frozen_comparison(name):
    """Intercept the unchanged hybrid run; compare its exact last-step problem."""
    P, G, L = prepare(parameters(name, "hybrid"))
    assemble, solve = granular.assemble, granular.solve_coulomb
    current, last = {}, {}

    def capture_assembly(P, G, batch, **options):
        result = assemble(P, G, batch, **options)
        current["operator"] = assemble(P, G, batch, compact=True)[0]
        inverse = np.zeros(2 * len(G.m))
        inverse[result[1]] = 1 / np.repeat(G.m, 2)[result[1]]
        current["inverse"] = inverse
        return result

    def capture_solve(W, q, mu, **options):
        result = solve(W, q, mu, **options)
        last.update(
            current, W=W, q=q.copy(), mu=mu, reference=result[0], initial=options["initial"].copy(), dt=P.dt
        )
        return result

    with (
        patch.object(granular, "assemble", capture_assembly),
        patch.object(granular, "solve_coulomb", capture_solve),
    ):
        run(P, G, L, output=False)
    operator, inverse = last["operator"], last["inverse"]
    start = perf_counter()
    result, stats = solve_velocity(
        operator,
        np.zeros_like(inverse),
        inverse,
        last["q"],
        last["mu"],
        last["dt"],
        tolerance=P.granular_tolerance,
        max_iterations=P.granular_max_iterations,
        initial=last["initial"],
    )
    elapsed = perf_counter() - start
    velocity_difference = last["dt"] * inverse * (operator.T @ (result - last["reference"]).ravel())
    residual = natural_residual(last["W"], last["q"], result, last["mu"])
    comparison = dict(
        original_residual=residual,
        requested_tolerance=P.granular_tolerance,
        max_grid_velocity_difference=float(np.max(np.abs(velocity_difference), initial=0)),
        residual_passed=bool(residual <= P.granular_tolerance),
        velocity_match=bool(np.max(np.abs(velocity_difference), initial=0) < 2e-8),
        velocity_acceptance_bound=2e-8,
        solve_seconds=elapsed,
        stats=stats,
    )
    comparison["velocity_equivalence_verified"] = comparison["velocity_match"]
    if not comparison["velocity_match"]:
        # Residual tolerances do not uniformly bound velocity error at degenerate
        # stick/slip transitions. Retain the failed matched-tolerance velocity
        # check, then independently test convergence under tighter solves.
        study = []
        for tolerance in (1e-11, 1e-12):
            reference, reference_stats = solve(
                last["W"],
                last["q"],
                last["mu"],
                tolerance=tolerance,
                initial=last["initial"],
                max_iterations=P.granular_max_iterations,
            )
            reaction, tight_stats = solve_velocity(
                operator,
                np.zeros_like(inverse),
                inverse,
                last["q"],
                last["mu"],
                last["dt"],
                tolerance=tolerance,
                initial=last["initial"],
                max_iterations=P.granular_max_iterations,
            )
            difference = last["dt"] * inverse * (operator.T @ (reaction - reference).ravel())
            study.append(
                dict(
                    tolerance=tolerance,
                    reference_residual=reference_stats["residual"],
                    compact_original_residual=natural_residual(last["W"], last["q"], reaction, last["mu"]),
                    max_grid_velocity_difference=float(np.max(np.abs(difference), initial=0)),
                    compact_sweeps=tight_stats["sweeps"],
                )
            )
        comparison["tolerance_study"] = study
        comparison["velocity_equivalence_verified"] = bool(
            study[-1]["max_grid_velocity_difference"] < 2e-8
            and study[-1]["compact_original_residual"] <= study[-1]["tolerance"]
        )
    return comparison


def timed_run(name, method):
    P, G, L = prepare(parameters(name, method))
    start = perf_counter()
    run(P, G, L, output=False)
    elapsed = perf_counter() - start
    b = L.materials[0]
    rows = L.diagnostics.history
    return (
        elapsed,
        dict(
            steps=P.tstep,
            particles=len(b),
            grid=[P.G.nx, P.G.ny],
            dt=P.dt,
            final_time=P.t,
            tolerance=P.granular_tolerance,
            transfer=P.transfer,
            stress=P.granular_stress,
            viscosity=P.granular_viscosity,
            max_active_nodes=max(r["granular_solver"]["active_nodes"] for r in rows),
            max_residual=max(r["granular_solver"]["residual"] for r in rows),
            max_sweeps=max(r["granular_solver"]["sweeps"] for r in rows),
            max_momentum_residual=max(float(np.linalg.norm(r["momentum_residual"])) for r in rows),
            max_integrated_overpacking=max(
                r["granular_solver"]["integrated_overpacking_fraction"] for r in rows
            ),
            max_packing=max(r["granular_solver"]["maximum_packing_after"] for r in rows),
        ),
        dict(x=b.x.copy(), v=b.v.copy()),
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("/tmp/granular-mixed.json"))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--baseline", type=Path, help="Optional pre-change granular_prototype report")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    root = Path(__file__).resolve().parents[1]
    names = [
        "granular.py",
        "granular_mixed.py",
        "initialise.py",
        "benchmarks/granular_mixed.py",
        "benchmarks/granular_prototype.py",
        "inputs/granular_column.py",
        "inputs/granular_collapse.py",
    ]
    report = dict(
        environment=dict(
            python=sys.version,
            platform=platform.platform(),
            numpy=np.__version__,
            scipy=scipy.__version__,
            thread_settings={key: os.environ.get(key) for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")},
        ),
        source_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in names},
        timing="One warmup per method, then alternating paired runs; medians. Setup/output/progress excluded, diagnostics enabled.",
        scope="Current Q1-Q1 nodal stress space and rigid nonassociated law; no density fix or general GPU speedup claim.",
        cases=[],
    )
    if args.baseline:
        report["prechange_baseline"] = json.loads(args.baseline.read_text())
    for name in ("static", "collapse", "refined", "incline_static", "incline_sliding", "plane_strain"):
        print(f"{name}: warming both solvers", flush=True)
        for method in ("hybrid", "velocity"):
            timed_run(name, method)
        samples = {method: [] for method in ("hybrid", "velocity")}
        summaries, states = {}, {}
        for repeat in range(args.repeats):
            methods = ("hybrid", "velocity") if repeat % 2 == 0 else ("velocity", "hybrid")
            for method in methods:
                elapsed, summaries[method], states[method] = timed_run(name, method)
                samples[method].append(elapsed)
        medians = {method: float(np.median(values)) for method, values in samples.items()}
        differences = {
            key: float(np.max(np.abs(states["velocity"][key] - states["hybrid"][key]), initial=0))
            for key in ("x", "v")
        }
        frozen = frozen_comparison(name)
        row = dict(
            case=name,
            samples_seconds=samples,
            median_seconds=medians,
            speedup=medians["hybrid"] / medians["velocity"],
            configurations=summaries,
            evolving_max_absolute_difference=differences,
            frozen_last_step=frozen,
        )
        report["cases"].append(row)
        print(
            f"{name}: hybrid {medians['hybrid']:.3f}s, velocity {medians['velocity']:.3f}s, "
            f"{row['speedup']:.2f}x; frozen residual {frozen['original_residual']:.3g}",
            flush=True,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    report["acceptance"] = dict(
        all_original_residuals_passed=all(
            row["frozen_last_step"]["residual_passed"] for row in report["cases"]
        ),
        all_same_tolerance_velocities_matched=all(
            row["frozen_last_step"]["velocity_match"] for row in report["cases"]
        ),
        all_velocity_equivalence_verified=all(
            row["frozen_last_step"]["velocity_equivalence_verified"] for row in report["cases"]
        ),
        interpretation="A failed matched-tolerance velocity bound is retained and requires an independent tighter-tolerance study; it is not silently waived.",
    )
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if not (
        report["acceptance"]["all_original_residuals_passed"]
        and report["acceptance"]["all_velocity_equivalence_verified"]
    ):
        raise SystemExit("frozen residual or velocity-equivalence study failed; inspect report")


if __name__ == "__main__":
    main()
