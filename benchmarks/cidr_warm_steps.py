"""Warm-started steps of the critical-state CIDR solve on the frozen collapse geometry.

benchmarks/cidr_frozen_step.py solves the first step from rest, the worst case for
an iterative solver: a cold start and the whole jammed set forming at once. In a
run the solution changes little from one step to the next. This study repeats the
step on the same geometry (particles and grid frozen, so it is a proxy for time
stepping, not a collapse): the trial velocity of step n + 1 is the velocity of
step n plus one step of gravity, and each solver starts from its own state at the
end of step n.

- alg2: carries the auxiliary strain rate, the multiplier and the penalty.
- newton: carries the multiplier (the stress); the first step reduces gamma from
  dt / rho to --gamma, and later steps stay at --gamma.

    python -m benchmarks.cidr_warm_steps --output benchmarks/cidr_warm_steps.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from benchmarks.cidr_frozen_step import GRAIN, Step, alg2, frozen, newton
from benchmarks.cidr_potential import Jop
from benchmarks.cidr_shear_cell import CriticalState


def march(p, model, C, steps, solver, *, tol=1e-5, reduction=np.sqrt(10), floor=None, max_outer=200, keep=()):
    """Repeat the frozen step `steps` times with warm starts; returns one record per step."""
    step = Step(p, model, C)
    increment = step.trial.copy()  # one step of gravity from rest (free components)
    dt, rho = p["dt"], p["rho"]
    floor = floor if floor is not None else dt / (1000 * rho)
    norm = lambda v: np.sqrt(v @ (step.mass * v))  # noqa: E731
    u = np.zeros_like(increment)
    d = np.zeros((step.n, 3))
    lam = np.zeros((step.n, 3))
    r = rho / dt
    records = []
    for n in range(steps):
        step.trial = u + increment
        trial = step.trial.copy()
        started = perf_counter()
        if solver == "alg2":
            u, d, lam, jam, r, iterations = alg2(step, u, d, lam, r, tol=tol, max_iterations=200000)
            record = dict(iterations=iterations)
        else:
            gamma = dt / rho if n == 0 else floor
            u = step.trial.copy() if n == 0 else u
            total = 0
            for outer in range(1, max_outer + 1):
                previous = u
                u, jam, lam, iterations = newton(step, u, gamma, lam, max_iterations=200)
                total += iterations
                if gamma <= floor and norm(u - previous) <= tol * norm(u):
                    break
                gamma = max(gamma / reduction, floor)
            else:
                raise RuntimeError(f"step {n}: method of multipliers did not settle")
            record = dict(iterations=total, outer=outer)
        speed = np.linalg.norm(u.reshape(-1, 2), axis=1).max() / (p["g"] * dt)
        record.update(
            step=n,
            seconds=perf_counter() - started,
            jammed=int(jam.sum()),
            max_speed=float(speed),
        )
        if n in keep:
            record["state"] = (trial, u.copy())
        records.append(record)
        print(solver, p["grid"][2], record, flush=True)
    return records


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("cidr-warm-steps.json"))
    parser.add_argument("--spacings", type=float, nargs="*", default=[0.05, 0.025, 0.0125])
    parser.add_argument("--steps", type=int, default=20)
    parser.add_argument("--solvers", nargs="*", default=["newton"])
    parser.add_argument(
        "--gamma", type=float, default=1e-8, help="newton: fixed gamma of steps after the first"
    )
    parser.add_argument("--tol", type=float, default=1e-5)
    args = parser.parse_args(argv)
    jop = Jop()
    report = dict(tol=args.tol, steps=args.steps, consolidation=2.0, cases=[])
    for spacing in args.spacings:
        p = frozen(spacing)
        C = 2.0 * p["rho"] * p["g"] / (2 * jop.mu1)
        model = CriticalState(
            mu_s=jop.mu1, mu_d=jop.mu2, i0=jop.I_0, k=2 * GRAIN * np.sqrt(p["rho"] / p["phi"])
        )
        for solver in args.solvers:
            records = march(p, model, C, args.steps, solver, tol=args.tol, floor=args.gamma)
            report["cases"].append(
                dict(spacing=spacing, stress_nodes=len(p["solid"]), solver=solver, steps=records)
            )
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
