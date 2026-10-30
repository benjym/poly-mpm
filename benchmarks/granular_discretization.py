"""Controlled Q1-Q1/Q1-P0 experiments; failures remain in the report.

The stress discretization changes, so runtime ratios are costs of different
discrete models, not equivalent-behavior solver speedups.
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
from scipy import sparse

import granular_mixed
from benchmarks.granular_prototype import configuration
from granular import cell_packing_fraction, solve_coulomb
from initialise import prepare
from MPM import run


def parameters(name, basis):
    refined = name.endswith("refined")
    h = 0.125 if refined else 0.25
    dt = (
        0.0025
        if refined or name == "collapse_coarse_short" or name.startswith("incline")
        else 0.02 if name == "static" else 0.005
    )
    static = name == "static" or name.startswith("incline")
    P = configuration(
        "static" if static else "collapse",
        dt,
        grid_spacing=h,
        particles_per_cell=4 if name == "collapse_ppc4" else 2,
        solver="velocity",
    )
    P.granular_stress_basis = basis
    if refined or name == "collapse_coarse_short":
        P.t_f = 0.05
    if name.startswith("incline"):
        P.B.walls, P.B.periodic, P.B.wall_contact = {"bottom"}, True, "no_slip"
        P.gravity_angle = np.arctan(0.75)
        P.t_f = 0.05
    return P


def measure(name, basis):
    P, G, L = prepare(parameters(name, basis))
    start = perf_counter()
    error = None
    try:
        run(P, G, L, output=False)
    except (RuntimeError, ValueError) as exc:
        error = str(exc)
    elapsed = perf_counter() - start
    b = L.materials[0]
    rows = L.diagnostics.history
    phi = cell_packing_fraction(P, G, b)
    row = dict(
        status="failed" if error else "completed",
        error=error,
        seconds=elapsed,
        completed_time=P.t,
        requested_time=P.t_f,
        steps=P.tstep,
        particles=len(b),
        grid=[P.G.nx, P.G.ny],
        dt=P.dt,
        tolerance=P.granular_tolerance,
        stress_basis=basis,
        solver=P.granular_solver,
        transfer=P.transfer,
        material=dict(rho_s=b.material.rho_s, phi_max=b.material.phi_max, mu_s=b.material.mu_s),
        extent=[P.G.x_m, P.G.x_M, P.G.y_m, P.G.y_M],
        output=False,
        diagnostics=True,
        maximum_speed=float(np.linalg.norm(b.v, axis=1).max()),
        rightmost_particle=float(b.x[:, 0].max()),
        final_cell_packing_maximum=float(phi.max()),
        final_cell_overpacking=float(
            np.maximum(phi - b.material.phi_max, 0).sum()
            * G.dx
            * G.dy
            * P.G.thickness
            / np.sum(b.m / b.material.rho_s)
        ),
    )
    if rows:
        row.update(
            max_nodal_packing=max(r["granular_solver"]["maximum_packing_after"] for r in rows),
            max_deformation_based_packing=max(r["granular_solver"]["maximum_particle_packing"] for r in rows),
            max_integrated_nodal_overpacking=max(
                r["granular_solver"]["integrated_overpacking_fraction"] for r in rows
            ),
            max_residual=max(r["granular_solver"]["residual"] for r in rows),
            max_active_stress_nodes=max(r["granular_solver"]["active_nodes"] for r in rows),
            max_colors=max(r["granular_solver"].get("colors", 0) for r in rows),
            max_sweeps=max(r["granular_solver"]["sweeps"] for r in rows),
            max_momentum_residual=max(float(np.linalg.norm(r["momentum_residual"])) for r in rows),
        )
    if name == "static" and not error:
        row["static_support_passed"] = bool(row["maximum_speed"] < 1e-8)
    if name.startswith("incline") and not error:
        expected = 9.81 * (np.sin(P.theta) - 0.5 * np.cos(P.theta)) * P.t
        actual = float((-b.v[:, 0]).max())
        row.update(
            uniform_sliding_reference_speed=float(expected),
            maximum_downslope_speed=actual,
            sliding_relative_error=float(abs(actual - expected) / expected),
            sliding_reference_passed=bool(abs(actual - expected) <= 1e-3 * expected),
        )
    return row, b.x.copy(), b.v.copy()


def first_step_reference(name, basis):
    """Replay an initial-step failure through the independent assembled solve."""
    P, G, L = prepare(parameters(name, basis))
    result = {}

    def assembled(operator, trial, inverse, gap, mu, dt, **options):
        B = sparse.csr_matrix(np.column_stack([operator @ column for column in np.eye(operator.shape[1])]))
        W = dt * (B.multiply(inverse) @ B.T)
        q = (B @ trial).reshape(gap.shape) + gap
        start = perf_counter()
        reaction, stats = solve_coulomb(W, q, mu, **options)
        result.update(
            status="completed",
            solve_seconds=perf_counter() - start,
            residual=stats["residual"],
            sweeps=stats["sweeps"],
            evaluations=stats["evaluations"],
        )
        return reaction, stats

    try:
        with patch.object(granular_mixed, "solve_velocity", assembled):
            run(P, G, L, max_steps=1, output=False)
    except RuntimeError as exc:
        result.update(status="failed", error=str(exc))
    return result


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("/tmp/granular-discretization.json"))
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    root = Path(__file__).resolve().parents[1]
    names = [
        "granular.py",
        "granular_mixed.py",
        "initialise.py",
        "benchmarks/granular_discretization.py",
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
        timing=f"Warm each basis, then alternate {args.repeats} paired runs; setup and output excluded, diagnostics enabled. Warmup failures are retained without repeat timing.",
        scope="Different discrete models: compare physical/packing accuracy separately from cost. Neither is established inf-sup stable.",
        cases=[],
    )
    for name in (
        "static",
        "collapse_ppc2",
        "collapse_ppc4",
        "collapse_coarse_short",
        "collapse_refined",
        "incline_coarse",
        "incline_refined",
    ):
        samples = {"nodal": [], "P0": []}
        states = {}
        summaries = {}
        failed = set()
        for basis in ("nodal", "P0"):
            row, x, v = measure(name, basis)
            if row["status"] == "failed":
                summaries[basis] = row
                states[basis] = (x, v)
                failed.add(basis)
                print(
                    f"{name} {basis}: warmup failed at {row['completed_time']:g}: {row['error']}", flush=True
                )
        for repeat in range(args.repeats):
            order = ("nodal", "P0") if repeat % 2 == 0 else ("P0", "nodal")
            for basis in order:
                if basis in failed:
                    continue
                row, x, v = measure(name, basis)
                summaries[basis] = row
                states[basis] = (x, v)
                if row["status"] == "failed":
                    failed.add(basis)
                else:
                    samples[basis].append(row["seconds"])
        medians = {basis: float(np.median(times)) if times else None for basis, times in samples.items()}
        row = dict(case=name, configurations=summaries, samples_seconds=samples, median_seconds=medians)
        if "P0" in failed and summaries["P0"]["completed_time"] == 0:
            row["P0_first_step_assembled_reference"] = first_step_reference(name, "P0")
        if not failed:
            row["runtime_ratio_nodal_to_P0"] = medians["nodal"] / medians["P0"]
            row["different_discretization_max_state_difference"] = {
                key: float(np.max(np.abs(a - b)))
                for key, a, b in zip(("x", "v"), states["nodal"], states["P0"])
            }
        report["cases"].append(row)
        print(f"{name}: medians {medians}; P0 status {summaries['P0']['status']}", flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
