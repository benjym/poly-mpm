"""Steady flow down a rough incline with critical-state CIDR, against the Bagnold profile.

In steady, uniform, isochoric flow CIDR reduces to mu(I): every depth flows at the
inertial number I_theta with mu(I_theta) = tan(theta), so the shear rate is

    du/dy = I_theta sqrt(p / rho_s) / d,

with p the weight of the grains above (per unit area, times cos theta). This script
runs inputs/cidr_incline.py, records the surface velocity until it settles, and
compares the shear rate in each layer with this prediction, using the simulated
packing through p. No parameter is fitted.

The layer starts in this state, so the check is that the solver holds it. It does
so for about 0.4 s: the spread of particle packing then grows until a node passes
phi_max (docs/granular-cidr-design.md, incline results), so the default end time is
0.4 s and each record carries the spread.

    python -m benchmarks.cidr_incline --output benchmarks/cidr_incline.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from initialise import get_parameters
from MPM import run


def profile(P, G, L, layers):
    """Layer-mean downslope velocity, shear rate and the Bagnold prediction from the particles."""
    b = L.materials[0]
    S = b.material
    y, u = b.x[:, 1], -b.v[:, 0]  # flow is towards -x
    # The surface is half a particle spacing (a quarter cell) above the top particles.
    top = y.max() + P.G.dy / 4
    edges = np.linspace(0.0, top, layers + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    index = np.clip(np.digitize(y, edges) - 1, 0, layers - 1)
    mean = np.array([u[index == k].mean() for k in range(layers)])
    width = P.G.x_max - P.G.x_min
    # Pressure from the weight above each layer boundary, at the layer centres.
    weight = abs(P.gravity) * np.cos(P.gravity_angle) / (width * P.G.thickness)
    pressure = np.array([weight * b.m[y > c].sum() for c in centres])
    mu_target = np.tan(P.gravity_angle)
    I_theta = S.I_0 * (mu_target - S.mu_s) / (S.mu_d - mu_target)
    rate = np.gradient(mean, centres)
    expected = I_theta * np.sqrt(pressure / S.rho_s) / S.grain_diameter
    packing = np.array([(b.rho[index == k] / S.rho_s).mean() for k in range(layers)])
    return dict(
        y=centres.tolist(),
        velocity=mean.tolist(),
        shear_rate=rate.tolist(),
        bagnold_shear_rate=expected.tolist(),
        pressure=pressure.tolist(),
        packing=packing.tolist(),
        I_theta=float(I_theta),
        surface=float(top),
    )


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("cidr-incline.json"))
    parser.add_argument("--backend", default="metal")
    parser.add_argument("--t-end", type=float, default=0.4)
    parser.add_argument("--chunk", type=float, default=0.1)
    parser.add_argument("--layers", type=int, default=18)  # one per cell, two particle rows each
    args = parser.parse_args(argv)
    P, G, L = get_parameters(["cidr_incline"])
    P.backend = args.backend
    history = []
    started = perf_counter()
    t = 0.0
    while t < args.t_end - 1e-12:
        t = min(t + args.chunk, args.t_end)
        P.t_f = t
        run(P, G, L, output=False)
        b = L.materials[0]
        surface = -b.v[b.x[:, 1] > b.x[:, 1].max() - P.G.dy, 0].mean()
        packing = b.rho / b.material.rho_s
        lower = b.x[:, 1] < 0.5 * b.x[:, 1].max()
        history.append(
            dict(
                t=float(P.t),
                surface_velocity=float(surface),
                steps=int(P.tstep),
                lower_packing=[float(packing[lower].min()), float(packing[lower].max())],
            )
        )
        print(history[-1], f"{perf_counter() - started:.0f}s", flush=True)
    result = profile(P, G, L, args.layers)
    rate, expected = np.array(result["shear_rate"]), np.array(result["bagnold_shear_rate"])
    interior = slice(2, args.layers - 3)  # away from the no-slip base and the free surface
    result.update(
        history=history,
        backend=args.backend,
        seconds=perf_counter() - started,
        interior_relative_error=float(np.max(np.abs(rate[interior] / expected[interior] - 1))),
        interior_mean_ratio=float(np.mean(rate[interior] / expected[interior])),
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"I_theta {result['I_theta']:.4f}; interior shear rate over Bagnold: mean "
        f"{result['interior_mean_ratio']:.3f}, largest deviation {result['interior_relative_error']:.3f}"
    )
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
