"""Volumetric locking: plane-strain cantilever under self-weight, with and without F-bar.

The beam (length L, height h) is fixed at a no-slip left wall and loaded by
gravity from rest. Undamped, its tip oscillates about the static deflection,
estimated as the midpoint of the tip displacement's extremes over at least one
and a half analytical periods. Euler-Bernoulli plane strain gives
delta = 1.5 rho g L^4 / (E' h^2) with E' = E / (1 - nu^2).

    uv run --locked python -m benchmarks.locking_cantilever --output benchmarks/locking_cantilever.json
"""

import argparse
import contextlib
import io
import json
import time
from types import SimpleNamespace as NS

import numpy as np

LENGTH, HEIGHT, RHO, E, GRAVITY = 1.0, 0.125, 1000.0, 1e8, -10.0


class Beam:
    def __init__(self, nu, fbar, cells):
        spacing = HEIGHT / cells
        bottom = 0.5 - HEIGHT / 2
        self.fbar = fbar
        self.dt = 1.0
        self.t_f = 1.0
        self.savetime = 1.0
        self.gravity = GRAVITY
        self.output_dir = "im/locking_cantilever/"
        self.time_stepping = "static"
        nx = int(round(1.25 / spacing)) + 1
        ny = int(round(0.5 / spacing)) + 1
        self.G = NS(
            nx=nx, ny=ny, x_min=0.0, x_max=(nx - 1) * spacing, y_min=0.25, y_max=0.25 + (ny - 1) * spacing
        )
        self.B = NS(walls={"left"}, wall_contact="no_slip")
        self.O = NS()
        quarter = spacing / 4
        xs = np.arange(quarter, LENGTH, spacing / 2)
        ys = np.arange(bottom + quarter, bottom + HEIGHT, spacing / 2)
        x, y = np.meshgrid(xs, ys)
        K = E / (3 * (1 - 2 * nu))
        G = E / (2 * (1 + nu))
        self.S = NS(
            law="elastic", n=x.size, X=x.ravel(), Y=y.ravel(), rho=RHO, A=(spacing / 2) ** 2, K=K, G=G
        )
        wave = np.sqrt((K + 4 * G / 3) / RHO)
        self.dt = 0.3 * spacing / wave


def analytical(nu):
    modulus = E / (1 - nu**2)
    deflection = 1.5 * RHO * abs(GRAVITY) * LENGTH**4 / (modulus * HEIGHT**2)
    stiffness = modulus * HEIGHT**3 / 12
    frequency = 1.875**2 / (2 * np.pi) * np.sqrt(stiffness / (RHO * HEIGHT * LENGTH**4))
    return deflection, 1 / frequency


def cantilever(nu, fbar, cells=4, periods=1.5):
    from initialise import prepare
    from MPM import step

    deflection, period = analytical(nu)
    params = Beam(nu, fbar, cells)
    params.t_f = periods * period
    with contextlib.redirect_stdout(io.StringIO()):
        P, G, L = prepare(params, "locking_cantilever")
    b = L.materials[0]
    tip = b.x[:, 0] > LENGTH - HEIGHT / cells
    start = b.x[tip, 1].mean()
    history = []
    began = time.perf_counter()
    while P.t < P.t_f:
        step(P, G, L)
        P.t += P.dt
        P.tstep += 1
        history.append(start - b.x[tip, 1].mean())
    history = np.asarray(history)
    estimate = 0.5 * (history.max() + history.min())
    return dict(
        nu=nu,
        fbar=fbar,
        cells_per_height=cells,
        steps=P.tstep,
        seconds=time.perf_counter() - began,
        analytical_deflection=deflection,
        static_deflection=float(estimate),
        ratio=float(estimate / deflection),
        peak_deflection=float(history.max()),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default="benchmarks/locking_cantilever.json")
    args = parser.parse_args()
    rows = [cantilever(nu, fbar) for nu in (0.3, 0.49, 0.499) for fbar in (False, True)]
    for row in rows:
        print(f"nu={row['nu']} fbar={row['fbar']}: ratio {row['ratio']:.3f} ({row['steps']} steps)")
    with open(args.output, "w") as handle:
        json.dump(dict(cases=rows), handle, indent=2)


if __name__ == "__main__":
    main()
