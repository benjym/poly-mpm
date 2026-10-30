"""Packing and solver acceptance studies; failures are retained as evidence.

Run python -m benchmarks.granular_packing --output /tmp/granular-packing.json.
Runtime comparisons use one bounded run and do not establish a general speedup.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np

from benchmarks.granular_prototype import configuration, measure
from granular import basis_volumes, packing_fraction
from initialise import prepare
from MPM import run


def sampling(particles_per_cell):
    """Exact area-preserving affine motion, without any numerical time solver."""
    P = configuration("static", 0.005, particles_per_cell=particles_per_cell)
    P.G.x_max, P.G.nx = 2.0, 9
    P.S.X += 0.5
    P.S.Y += 0.25
    P, G, L = prepare(P)
    b = L.materials[0]
    stretch = np.exp(np.array([0.1, -0.1]))
    centre = np.array([1.0, 0.75])
    before = packing_fraction(P, G, b)
    b.x[:, :2] = centre + (b.x[:, :2] - centre) * stretch
    after = packing_fraction(P, G, b)
    lower, upper = centre - 0.5 * stretch, centre + 0.5 * stretch
    interior = (
        (G.X - G.dx >= lower[0])
        & (G.X + G.dx <= upper[0])
        & (G.Y - G.dy >= lower[1])
        & (G.Y + G.dy <= upper[1])
    )
    assert interior.any()
    volumes = basis_volumes(P, G)
    return dict(
        particles_per_cell_axis=particles_per_cell,
        particles=len(b),
        affine_determinant=float(np.prod(stretch)),
        exact_packing=b.material.phi_max,
        initial_peak=float(before.max()),
        interior_peak=float(after[interior].max()),
        interior_max_error=float(np.max(np.abs(after[interior] - b.material.phi_max))),
        total_solid_volume_error=float(after @ volumes - b.m.sum() / b.material.rho_s),
    )


def incline(slope, method):
    P = configuration("static", 0.005, solver=method)
    P.B.walls, P.B.periodic, P.B.wall_contact = {"bottom"}, True, "no_slip"
    P.gravity_angle = np.arctan(slope)
    P.t_f = 0.05
    P, G, L = prepare(P)
    start = perf_counter()
    try:
        run(P, G, L, output=False)
    except RuntimeError as exc:
        return dict(
            slope=slope,
            method=method,
            status="solver failed",
            completed_time=P.t,
            seconds=perf_counter() - start,
            error=str(exc),
        )
    b = L.materials[0]
    speed = -b.v[:, 0]
    expected = max(0.0, 9.81 * (np.sin(P.theta) - b.material.mu_s * np.cos(P.theta))) * P.t
    return dict(
        slope=slope,
        method=method,
        status="completed",
        completed_time=P.t,
        seconds=perf_counter() - start,
        max_speed=float(np.linalg.norm(b.v, axis=1).max()),
        maximum_downslope_speed=float(speed.max()),
        uniform_sliding_reference_speed=float(expected),
        max_contact_residual=max(r["granular_solver"]["residual"] for r in L.diagnostics.history),
        max_contact_power=max(r["granular_solver"]["contact_power"] for r in L.diagnostics.history),
    )


def signed_gap_energy():
    """Analytic counterexample: velocity-level geometric repair can add energy.

    One mass, one normal constraint: u*=0, B=M=dt=1, k=-overlap.
    Complementarity gives r=overlap and u=overlap, without external work.
    """
    overlap = 0.01
    return dict(
        overlap=overlap,
        initial_kinetic_energy=0.0,
        final_kinetic_energy=0.5 * overlap**2,
        external_work=0.0,
        contact_power=overlap**2,
        conclusion="Do not adopt signed-gap velocity repair without energy accounting.",
    )


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("/tmp/granular-packing.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    names = [
        "granular.py",
        "initialise.py",
        "benchmarks/granular_prototype.py",
        "benchmarks/granular_packing.py",
        "inputs/granular_column.py",
        "inputs/granular_collapse.py",
    ]
    report = dict(
        environment=dict(python=sys.version, platform=platform.platform()),
        source_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in names},
        scope="Independent particle-density and grid-spacing sensitivity; no spatial convergence claim.",
        sampling=[sampling(n) for n in (2, 4, 8, 16)],
        collapse=[],
        inclined_layers=[],
        signed_gap_energy=signed_gap_energy(),
    )
    for h, ppc, method in [(0.25, 2, "hybrid"), (0.25, 2, "block"), (0.25, 4, "block"), (0.125, 2, "block")]:
        row = measure("collapse", 0.005 * h / 0.25, grid_spacing=h, particles_per_cell=ppc, solver=method)
        row["particles_per_cell_axis"] = ppc
        report["collapse"].append(row)
        print(
            f"collapse h={h} ppc={ppc} {method}: {row['seconds']:.3f}s, "
            f"packing={row['maximum_final_nodal_packing']:.6f}, "
            f"particle packing={row['maximum_particle_packing']:.6f}",
            flush=True,
        )
    for slope, method in [(0.25, "hybrid"), (0.75, "hybrid"), (0.75, "block")]:
        row = incline(slope, method)
        report["inclined_layers"].append(row)
        print(f"incline {slope} {method}: {row['status']}", flush=True)
    report["interpretation"] = (
        "Point-sampling packing errors occur even under exact area-preserving motion. "
        "Collapse changes under particle/grid refinement; no packing or spatial-convergence gate passed. "
        "Sparse block storage avoids a dense global factorization but the block solver fails one "
        "inclined-flow case; retain hybrid as the default reference."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
