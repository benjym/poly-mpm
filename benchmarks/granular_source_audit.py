"""Matched ablations of source-paper choices, not a 3D paper reproduction."""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
from benchmarks.granular_prototype import configuration
from initialise import prepare
from MPM import run


def measure(transfer, stress, viscosity, *, final_time=0.3, particles_per_cell=2):
    P = configuration("collapse", 0.005, particles_per_cell=particles_per_cell)
    P.transfer, P.granular_stress, P.granular_viscosity = transfer, stress, viscosity
    P.t_f = final_time
    P, G, L = prepare(P)
    b = L.materials[0]
    initial_volume = float(b.V.sum())
    initial_state = hashlib.sha256(b.x.tobytes() + b.m.tobytes() + b.v.tobytes()).hexdigest()
    start = perf_counter()
    error = None
    try:
        run(P, G, L, output=False)
    except RuntimeError as exc:
        error = str(exc)
    rows = L.diagnostics.history
    result = dict(
        transfer=transfer,
        stress=stress,
        viscosity=viscosity,
        requested_time=final_time,
        completed_time=P.t,
        dt=P.dt,
        particles_per_cell_axis=particles_per_cell,
        particles=len(b),
        initial_state_sha256=initial_state,
        seconds=perf_counter() - start,
        status="solver failed" if error else "completed",
        error=error,
    )
    if rows:
        result.update(
            rightmost_particle=float(b.x[:, 0].max()),
            centre_of_mass=(np.sum(b.m[:, None] * b.x, axis=0) / b.m.sum()).tolist(),
            peak_packing_ratio=max(r["granular_solver"]["maximum_packing_after"] for r in rows)
            / b.material.packing_max,
            max_integrated_overpacking_fraction=max(
                r["granular_solver"]["integrated_overpacking_fraction"] for r in rows
            ),
            deformation_volume_change=float(b.V.sum() / initial_volume - 1),
            max_contact_residual=max(r["granular_solver"]["residual"] for r in rows),
            max_mfb_evaluations=max(r["granular_solver"].get("mfb_evaluations", 0) for r in rows),
            max_mass_error=max(abs(r["mass_error"]) for r in rows),
            max_momentum_residual=max(float(np.linalg.norm(r["momentum_residual"])) for r in rows),
            energy_residual=rows[-1]["energy_residual"],
        )
    return result


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("/tmp/granular-source-audit.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    names = [
        "granular.py",
        "initialise.py",
        "particles.py",
        "diagnostics.py",
        "benchmarks/granular_source_audit.py",
        "benchmarks/sand6_oracle.cpp",
        "tests/fixtures/sand6_frozen.json",
    ]
    fixture = json.loads((root / names[-1]).read_text())
    report = dict(
        environment=dict(python=sys.version, platform=platform.platform()),
        source_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in names},
        original_source=dict(archive=fixture["source_url"], sha256=fixture["archive_sha256"]),
        scope="Controlled 2D/plane-strain ablations on our geometry; no Sand6 trajectory or experimental reproduction.",
        volume_metric="Integrated nodal excess / total fixed grain volume. Not the unspecified total-volume estimator used in Rousseau et al. 2023; deformation volume is also separate.",
        cases=[],
    )
    for transfer, stress, eta, ppc in [
        ("flip", "intrinsic", 0.0, 2),
        ("apic_gradient", "intrinsic", 0.0, 2),
        ("apic_gradient", "intrinsic", 0.001, 2),
        ("apic_gradient", "plane_strain", 0.001, 2),
        ("apic_gradient", "intrinsic", 0.001, 6),
    ]:
        row = measure(transfer, stress, eta, particles_per_cell=ppc)
        report["cases"].append(row)
        print(
            f"{transfer} {stress} eta={eta} ppc={ppc}: {row['status']} at {row['completed_time']:.3f}s; "
            f"peak ratio={row.get('peak_packing_ratio', float('nan')):.6f}",
            flush=True,
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
