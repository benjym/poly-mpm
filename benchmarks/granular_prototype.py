"""Bounded verification and timestep sensitivity of the rigid-grain prototype.

No speedup against explicit elastic MPM is claimed: the constitutive model is
intentionally different. Run with --output /tmp/granular-prototype.json.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
from initialise import prepare
from MPM import run
from inputs.granular_column import Params as Column
from inputs.granular_collapse import Params as Collapse


def configuration(kind, dt, *, grid_spacing=0.25, particles_per_cell=2, solver="hybrid"):
    P = (Column if kind == "static" else Collapse)([])
    P.G.nx = round((P.G.x_M - P.G.x_m) / grid_spacing) + 1
    P.G.ny = round((P.G.y_M - P.G.y_m) / grid_spacing) + 1
    spacing = grid_spacing / particles_per_cell
    width = 1.0 if kind == "static" else 0.5
    x, y = np.meshgrid(np.arange(spacing / 2, width, spacing), np.arange(spacing / 2, 1.0, spacing))
    P.S.X, P.S.Y = x.ravel(), y.ravel()
    P.S.n, P.S.A = x.size, spacing**2
    P.dt = dt
    P.granular_solver = solver
    P.t_f = 0.2 if kind == "static" else 0.3
    return P


def measure(kind, dt, **options):
    P = configuration(kind, dt, **options)
    P, G, L = prepare(P)
    initial = L.materials[0].x.copy()
    start = perf_counter()
    run(P, G, L, output=False)
    elapsed = perf_counter() - start
    b = L.materials[0]
    rows = L.diagnostics.history
    return dict(
        case=kind,
        dt=dt,
        steps=P.tstep,
        seconds=elapsed,
        particles=len(b),
        grid=dict(nx=P.G.nx, ny=P.G.ny, dx=G.dx, dy=G.dy),
        material=dict(rho_s=b.material.rho_s, packing_max=b.material.packing_max, mu=b.material.mu),
        final_time=P.t,
        output=False,
        diagnostics=True,
        interpolation=P.interpolation,
        transfer=P.transfer,
        contact_solver=P.granular_solver,
        max_displacement=float(np.linalg.norm(b.x - initial, axis=1).max()),
        max_speed=float(np.linalg.norm(b.v, axis=1).max()),
        rightmost_particle=float(b.x[:, 0].max()),
        centre_of_mass=(np.sum(b.m[:, None] * b.x, axis=0) / b.m.sum()).tolist(),
        max_mass_error=max(abs(r["mass_error"]) for r in rows),
        max_momentum_residual=max(float(np.linalg.norm(r["momentum_residual"])) for r in rows),
        max_natural_map_residual=max(r["granular_solver"]["residual"] for r in rows),
        maximum_initial_nodal_packing=max(r["granular_solver"]["maximum_packing"] for r in rows),
        maximum_final_nodal_packing=max(r["granular_solver"]["maximum_packing_after"] for r in rows),
        maximum_particle_packing=max(r["granular_solver"]["maximum_particle_packing"] for r in rows),
        maximum_contact_nodes=max(r["granular_solver"]["active_nodes"] for r in rows),
        maximum_operator_nnz=max(r["granular_solver"]["operator_nnz"] for r in rows),
        maximum_operator_storage_bytes=max(r["granular_solver"]["operator_storage_bytes"] for r in rows),
        maximum_contact_power=max(r["granular_solver"]["contact_power"] for r in rows),
        max_fallback_sweeps=max(r["granular_solver"].get("sweeps", 0) for r in rows),
        accumulated_contact_dissipation=float(b.plastic_work.sum()),
        energy_residual=rows[-1]["energy_residual"],
        timestep_limits=rows[-1]["timestep_limits"],
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("/tmp/granular-prototype.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    files = [
        "granular.py",
        "initialise.py",
        "MPM.py",
        "constit.py",
        "particles.py",
        "stability.py",
        "inputs/granular_column.py",
        "inputs/granular_collapse.py",
    ]
    report = dict(
        environment=dict(python=sys.version, platform=platform.platform()),
        source_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in files},
        scope="Intrinsic 2D rigid-grain prototype; no experimental validation or continuum well-posedness claim.",
        timing="One bounded run per configuration, setup and plotting excluded; diagnostics enabled.",
        cases=[],
    )
    for kind, dt in [("static", 0.02), ("collapse", 0.01), ("collapse", 0.005), ("collapse", 0.0025)]:
        row = measure(kind, dt)
        report["cases"].append(row)
        print(
            f"{kind} dt={dt:g}: {row['steps']} steps, {row['seconds']:.3f}s; "
            f"vmax={row['max_speed']:.6g}, xmax={row['rightmost_particle']:.6g}",
            flush=True,
        )
    collapse = report["cases"][1:]
    differences = [
        abs(collapse[i]["rightmost_particle"] - collapse[i + 1]["rightmost_particle"]) for i in range(2)
    ]
    report["timestep_comparison"] = dict(
        rightmost_particle_differences=differences,
        interpretation="Single coarse-grid temporal sensitivity check; not spatial convergence or experimental validation.",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
