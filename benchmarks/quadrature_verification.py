"""Controlled point/domain integration comparisons; all refinement levels retained."""

import argparse
import hashlib
import json
import platform
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import scipy
from .broader_verification import (
    force_balance,
    large_deformation,
    hydrostatic,
    translated_stress,
    manufactured_motion,
)


def timed(function, *args, **kwargs):
    start = perf_counter()
    result = function(*args, **kwargs)
    result["elapsed_seconds"] = perf_counter() - start
    return result


def geometry(cells):
    X = (np.arange(2 * cells) + 0.5) / (2 * cells)
    exact = manufactured_motion(X, 0.25)
    half = exact["F"][:, 0, 0] / (4 * cells)
    gap = np.roll(exact["x"] - half, -1) - (exact["x"] + half)
    gap[-1] += 1
    return dict(
        cells=cells,
        maximum_gap=float(max(0, gap.max())),
        maximum_overlap=float(max(0, -gap.min())),
        gap_rms=float(np.sqrt(np.mean(gap**2))),
    )


def report():
    root = Path(__file__).resolve().parents[1]
    sources = (
        "MPM.py",
        "initialise.py",
        "particles.py",
        "interpolation.py",
        "quadrature.py",
        "constit.py",
        "grid.py",
        "kinematics.py",
        "diagnostics.py",
        "integrators.py",
        "transfers.py",
        "benchmarks/broader_verification.py",
        "benchmarks/quadrature_verification.py",
    )
    result = dict(
        environment=dict(
            python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()
        ),
        source_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in sources},
        controls=dict(
            interpolation="quadratic",
            transfer="apic",
            particles_per_axis=2,
            spatial="same physical final time, 256 timesteps and transfer count; grid/cloud resolution varied",
            temporal="fixed 16-cell grid/cloud; timestep and remap frequency vary together",
            phase="frozen exact finite-strain state, shifted across one grid cell without dynamics",
            geometry="affine F-convected particle rectangles; gaps/overlaps are not artificially removed",
            timing="one wall-clock measurement per trajectory including diagnostics, not a throughput benchmark",
        ),
    )
    print("Frozen force integration and cell phases", flush=True)
    result["frozen"] = {
        case: [force_balance(case, n, integration=q) for q in ("point", "domain") for n in (8, 16, 32)]
        for case in ("hydrostatic", "translated_stress", "large_deformation")
    }
    result["cell_crossing"] = [
        force_balance("large_deformation", 16, integration=q, offset=float(f / 16))
        for q in ("point", "domain")
        for f in np.linspace(0, 1, 9)
    ]
    print("Finite-strain dynamic grid/time/particle studies", flush=True)
    result["spatial"] = [
        timed(large_deformation, n, 256, quadrature=q) for q in ("point", "domain") for n in (8, 16, 32)
    ]
    result["temporal"] = [
        timed(large_deformation, 16, n, quadrature=q) for q in ("point", "domain") for n in (128, 256, 512)
    ]
    result["sampling"] = [
        timed(large_deformation, 16, 256, ppc=n, quadrature=q) for q in ("point", "domain") for n in (2, 4)
    ]
    print("Hydrostatic and moving-stress trajectories", flush=True)
    result["hydrostatic"] = [timed(hydrostatic, 16, 128, quadrature=q) for q in ("point", "domain")]
    result["translated_stress"] = [timed(translated_stress, 8, 64, quadrature=q) for q in ("point", "domain")]
    result["domain_geometry"] = [geometry(n) for n in (8, 16, 32, 64)]
    return result


def figure(result, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for q in ("point", "domain"):
        rows = [r for r in result["frozen"]["large_deformation"] if r["integration"] == q]
        axes[0, 0].loglog([r["cells"] for r in rows], [r["acceleration_rms"] for r in rows], "o-", label=q)
        rows = [r for r in result["cell_crossing"] if r["integration"] == q]
        axes[0, 1].plot(
            [16 * r["offset"] for r in rows], [r["acceleration_rms"] for r in rows], "o-", label=q
        )
        rows = [r for r in result["spatial"] if r["quadrature"] == q]
        axes[1, 0].loglog([r["cells"] for r in rows], [r["stress_relative_l2"] for r in rows], "o-", label=q)
        rows = [r for r in result["temporal"] if r["quadrature"] == q]
        axes[1, 1].loglog([r["steps"] for r in rows], [r["stress_relative_l2"] for r in rows], "o-", label=q)
    axes[0, 0].set(
        title="Frozen finite-strain state", xlabel="Grid cells", ylabel="Acceleration defect RMS (m/s²)"
    )
    axes[0, 1].set(
        title="Frozen state translated across one cell",
        xlabel="Translation / grid spacing",
        ylabel="Acceleration defect RMS (m/s²)",
    )
    axes[1, 0].set(
        title="Finite-strain dynamics: timestep fixed", xlabel="Grid cells", ylabel="Relative stress L2 error"
    )
    axes[1, 1].set(
        title="Finite-strain dynamics: grid fixed",
        xlabel="Number of timesteps",
        ylabel="Relative stress L2 error",
    )
    for ax in axes.ravel():
        ax.grid(True, alpha=0.25)
        ax.legend()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", default="/tmp/poly-mpm-quadrature.json")
    parser.add_argument("--figure")
    args = parser.parse_args()
    result = report()
    Path(args.output).write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if args.figure:
        figure(result, args.figure)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
