"""Reproducible modern-solver verification, separate from historical timings.

Run python -m benchmarks.verification --output /tmp/poly-mpm-verification.json.
Small-amplitude periodic standing wave: v_x=a sin(2 pi x) cos(2 pi c t),
c=sqrt((K+4G/3)/rho). Refinement decreases both dt and h; it is a coupled
convergence study, not a claim of a purely spatial asymptotic convergence rate.
"""

import argparse
import json
import hashlib
import platform
import sys
import scipy
from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
from initialise import prepare
from MPM import run


def wave(cells, interpolation="quadratic", transfer="apic", dt_factor=0.08, fixed_steps=None):
    K, G, rho, amplitude = 1e4, 4e3, 1000.0, 1e-5
    speed = np.sqrt((K + 4 * G / 3) / rho)
    xx = (np.arange(2 * cells) + 0.5) / (2 * cells)
    yy = np.array([0.025, 0.075, 0.125, 0.175])
    X, Y = np.meshgrid(xx, yy)
    X = X.ravel()
    Y = Y.ravel()
    velocity = np.zeros((len(X), 3))
    velocity[:, 0] = amplitude * np.sin(2 * np.pi * X)
    C = np.zeros((len(X), 3, 3))
    C[:, 0, 0] = 2 * np.pi * amplitude * np.cos(2 * np.pi * X)
    S = NS(
        law="elastic", n=len(X), X=X, Y=Y, A=0.2 / len(X), rho=rho, K=K, G=G, initial_v=velocity, initial_C=C
    )
    P = NS(
        CFL=0.5,
        dt=(0.2 / speed / fixed_steps if fixed_steps else dt_factor / (cells * speed)),
        t_f=0.2 / speed,
        savetime=1.0,
        max_g=0.0,
        S=S,
        G=NS(nx=cells + 1, ny=3, x_m=0.0, x_M=1.0, y_m=0.0, y_M=0.2),
        B=NS(cyclic_lr=True, has_bottom=True, has_top=True),
        O=NS(),
        interpolation=interpolation,
        transfer=transfer,
    )
    P, grid, L = prepare(P)
    run(P, grid, L, output=False)
    b = L.materials[0]
    exact = amplitude * np.sin(2 * np.pi * b.x0[:, 0]) * np.cos(2 * np.pi * speed * P.t)
    row = L.diagnostics.history[-1]
    initial = L.diagnostics.initial
    return dict(
        fixed_steps=fixed_steps,
        cells=cells,
        interpolation=interpolation,
        transfer=transfer,
        steps=P.tstep,
        velocity_relative_l2=float(np.linalg.norm(b.v[:, 0] - exact) / np.linalg.norm(exact)),
        relative_energy_residual=float(
            row["energy_residual"] / (initial["kinetic"] + initial["affine_kinetic"])
        ),
        maximum_momentum_error=float(np.max(np.abs(row["momentum"]))),
        minimum_J=row["minimum_J"],
        maximum_J=row["maximum_J"],
    )


def full_examples():
    from initialise import get_parameters

    rows = []
    for example in ("modern_ball", "dp_rate_unit_test"):
        P, G, L = get_parameters([example])
        run(P, G, L, output=False)
        rows.append(
            dict(
                example=example,
                steps=P.tstep,
                time=float(P.t),
                minimum_J=min(float(b.J.min()) for b in L.materials if len(b)),
                maximum_J=max(float(b.J.max()) for b in L.materials if len(b)),
                final_diagnostics=L.diagnostics.last,
            )
        )
    return rows


def report():
    root = Path(__file__).resolve().parents[1]
    sources = [
        "MPM.py",
        "particles.py",
        "grid.py",
        "constit.py",
        "interpolation.py",
        "kinematics.py",
        "diagnostics.py",
        "stability.py",
        "integrators.py",
    ]
    return dict(
        description=__doc__,
        environment=dict(
            python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()
        ),
        source_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in sources},
        full_examples=full_examples(),
        spatial=[wave(n, fixed_steps=256) for n in (8, 16, 32, 64)],
        temporal=[wave(64, fixed_steps=n) for n in (64, 128, 256, 512)],
        waves=[
            wave(n, basis, transfer)
            for basis, transfer in [
                ("linear", "flip"),
                ("quadratic", "flip"),
                ("quadratic", "pic"),
                ("quadratic", "apic"),
            ]
            for n in (8, 16, 32)
        ],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", default="/tmp/poly-mpm-verification.json")
    parser.add_argument("--figure", help="optional PNG path for the refinement comparison")
    args = parser.parse_args()
    results = report()
    Path(args.output).write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    if args.figure:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
        for basis, transfer in [
            ("linear", "flip"),
            ("quadratic", "flip"),
            ("quadratic", "pic"),
            ("quadratic", "apic"),
        ]:
            rows = [r for r in results["waves"] if r["interpolation"] == basis and r["transfer"] == transfer]
            label = f"{basis} / {transfer.upper()}"
            axes[0].loglog(
                [r["cells"] for r in rows], [r["velocity_relative_l2"] for r in rows], "o-", label=label
            )
            axes[1].loglog(
                [r["cells"] for r in rows],
                [abs(r["relative_energy_residual"]) for r in rows],
                "o-",
                label=label,
            )
        axes[0].set(ylabel="Relative velocity L2 error")
        axes[1].set(ylabel="Absolute relative energy residual")
        for ax in axes:
            ax.set(xlabel="Cells along periodic wave direction")
            ax.grid(True, alpha=0.2)
            ax.legend(fontsize=8)
        fig.suptitle("Standing wave: coupled grid/time refinement (not a universal accuracy ranking)")
        fig.savefig(args.figure, dpi=150)
        plt.close(fig)
    for row in results["waves"]:
        print(row)
