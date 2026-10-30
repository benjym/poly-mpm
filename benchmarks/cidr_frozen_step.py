"""Frozen first collapse step under critical-state CIDR, refined.

The two-dimensional counterpart of the shear cell (benchmarks/cidr_shear_cell.py)
and the CIDR entry of the frozen-step study
(benchmarks/granular_regularisation_study.py): the first time step of the column
collapse from rest, solved at successively finer grids with the convex step of
docs/granular-cidr-design.md,

    min_u  1/2 |u - u*|^2_(M/dt) + sum_a V_a psi(D_a(u); C),

where a runs over the active stress nodes (the nodal basis of the rigid-grain
solver), D_a = (B u)_a / Vs_a is the solid-volume-weighted strain rate there, in
the orthonormal basis of granular.BASIS (component 0 volumetric, so
delta = sqrt(2) D_0; components 1 and 2 deviatoric, so s = |D_dev| / sqrt(2)),
and V_a = Vs_a / phi is the bulk volume. The packing is uniform and frozen. C is
set by the consolidation ratio 2 mu_s C / (rho g H), the cap pressure over the
pressure at the base of the column: at 1 the base sits on the cap and compacts
under the column's weight, which then falls freely through the first step; the
default 2 puts the base at the critical pressure, so the whole column is denser
than the cap for its load.

psi has a cone at D = 0 (the static yield), which is not separable in the nodal
velocities. Both solvers keep an auxiliary strain rate d_a = D_a(u) at each stress
node and use the exact local proximal step of psi, which decides whether a node
is jammed:

- alg2: augmented Lagrangian splitting (over-relaxed, with residual balancing of
  the penalty): a sparse solve for u, the proximal step per node, and a
  multiplier update.
- newton: the method of multipliers, each subproblem minimised over u by
  semismooth Newton on the Moreau envelope psi_gamma (whose gradient is given by
  the proximal step), with gamma reduced along the outer steps.

In both the multiplier converges to the stress.

    python -m benchmarks.cidr_frozen_step --solver alg2 --spacings 0.05 0.025 0.0125 0.00625 \
        --output benchmarks/cidr_frozen_step_alg2.json
    python -m benchmarks.cidr_frozen_step --solver newton --output benchmarks/cidr_frozen_step_newton.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from cidr import ROOT2, Problem, newton
from benchmarks.cidr_potential import Jop
from benchmarks.cidr_shear_cell import CriticalState
from benchmarks.granular_regularisation_study import PROBES, sample

GRAIN = 0.003  # grain diameter (m), as in benchmarks/cidr_potential.py


def frozen(spacing):
    """First step of the collapse from rest: operators, volumes, masses and grid."""
    from granular import assemble, canonical_nodes
    from mpm_gpu.benchmark import columns
    from transfers import nodal_velocity, scatter

    P, G, L = columns(1, spacing=spacing, particle_spacing=spacing / 4)
    b = L.materials[0]
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    body = np.array([P.g * np.sin(P.theta), P.g * np.cos(P.theta), 0.0])
    scatter(G.fe, b.n, b.N[:, :, None] * b.m[:, None, None] * body)
    G.BCs(P)
    G.update_momentum(P)
    trial = nodal_velocity(G, P)[:, :2].ravel()
    _, free, active, _, _ = assemble(P, G, b, compact=True)
    B, *_ = assemble(P, G, b, compact=False)
    canonical = canonical_nodes(P, len(G.m))
    solid = np.bincount(
        canonical[b.n].ravel(), weights=((b.m / b.material.rho_s)[:, None] * b.N).ravel(), minlength=len(G.m)
    )[active]
    return dict(
        B=sparse.csr_matrix(B),
        free=free,
        mass=np.repeat(G.m, 2),
        trial=trial,
        solid=solid,
        phi=b.material.rho / b.material.rho_s,
        rho=b.material.rho,
        dt=P.dt,
        g=abs(P.g),
        grid=(P.G.x_m, P.G.y_m, G.dx, G.dy, P.G.nx, P.G.ny),
        particles=len(b),
    )


def Step(p, model, C):
    """The frozen step as a cidr.Problem: D(u) = A u + a0 per stress node, bulk volumes W."""
    free, dt = p["free"], p["dt"]
    Vs = p["solid"]
    inverse = sparse.diags(np.repeat(1 / Vs, 3))
    return Problem(
        (inverse @ p["B"][:, free]).tocsr(),
        inverse @ (p["B"][:, ~free] @ p["trial"][~free]),
        np.repeat(Vs / p["phi"], 3),
        p["mass"][free] / dt,
        p["trial"][free].copy(),
        model,
        np.full(len(Vs), C),
    )


def alg2(step, u, d, lam, r, *, relax=1.8, tol=1e-6, max_iterations=20000, report=500, history=None):
    """Over-relaxed ADMM with residual balancing, from (d, lam); returns u, d, lam, jam, r, iterations."""
    A, a0, W, mass = step.A, step.a0, step.W, step.mass

    def factorise(r):
        return splu((sparse.diags(mass) + r * (A.T @ sparse.diags(W) @ A)).tocsc())

    factor = factorise(r)
    target = mass * step.trial - r * (A.T @ (W * a0))
    for iteration in range(1, max_iterations + 1):
        u = factor.solve(target - A.T @ (W * (lam.ravel() - r * d.ravel())))
        Du = step.D(u)
        d_old = d
        Dr = relax * Du + (1 - relax) * d_old
        d, jam, _ = step.prox(Dr + lam / r, r)
        lam = lam + r * (Dr - d)
        primal = np.sqrt(np.sum(W * (Du - d).ravel() ** 2) / max(np.sum(W * Du.ravel() ** 2), 1e-300))
        dual = np.sqrt(np.sum(W * (d - d_old).ravel() ** 2) / max(np.sum(W * d.ravel() ** 2), 1e-300))
        if history is not None and iteration % report == 0:
            history.append(
                dict(iteration=iteration, primal=float(primal), dual=float(dual), jammed=int(jam.sum()))
            )
        if primal <= tol and dual <= tol:
            break
        if iteration % 100 == 0 and max(primal, dual) > 10 * min(primal, dual):
            # Residual balancing: a lagging primal residual calls for a stiffer penalty.
            r *= 2.0 if primal > dual else 0.5
            factor = factorise(r)
            target = mass * step.trial - r * (A.T @ (W * a0))
    return u, d, lam, jam, r, iteration


def solve(p, model, C, *, r=None, tol=1e-6, max_iterations=100000, report=500):
    """ALG2 from rest; returns velocity, stress (the multiplier) and diagnostics."""
    step = Step(p, model, C)
    if r is None:
        r = p["rho"] * 1.0**2 / p["dt"]
    n = step.n
    started = perf_counter()
    history = []
    u, d, lam, jam, r, iterations = alg2(
        step,
        step.trial.copy(),
        np.zeros((n, 3)),
        np.zeros((n, 3)),
        r,
        tol=tol,
        max_iterations=max_iterations,
        report=report,
        history=history,
    )
    if iterations == max_iterations:
        raise RuntimeError(f"ALG2 did not converge in {max_iterations} iterations")
    full = p["trial"].copy()
    full[p["free"]] = u
    return (
        full,
        lam,
        dict(
            iterations=iterations,
            jammed_fraction=float(jam.mean()),
            seconds=perf_counter() - started,
            final_r=r,
            history=history,
            max_pressure=float(-lam[:, 0].min() / ROOT2),
        ),
    )


def prolong(coarse, u, fine):
    """Bilinear interpolation of nodal velocity u on grid `coarse` to the nodes of grid `fine`."""
    x_m, y_m, dx, dy, nx, ny = coarse["grid"]
    fx_m, fy_m, fdx, fdy, fnx, fny = fine["grid"]
    X, Y = np.meshgrid(fx_m + fdx * np.arange(fnx), fy_m + fdy * np.arange(fny))
    v = u.reshape(ny, nx, 2)
    fx = np.clip((X.ravel() - x_m) / dx, 0, nx - 1 - 1e-9)
    fy = np.clip((Y.ravel() - y_m) / dy, 0, ny - 1 - 1e-9)
    i, j = np.floor(fx).astype(int), np.floor(fy).astype(int)
    a, c = (fx - i)[:, None], (fy - j)[:, None]
    value = (
        (1 - a) * (1 - c) * v[j, i]
        + a * (1 - c) * v[j, i + 1]
        + a * c * v[j + 1, i + 1]
        + (1 - a) * c * v[j + 1, i]
    )
    return value.ravel()


def solve_newton(
    p, model, C, *, u0=None, gamma0=None, floor=None, reduction=np.sqrt(10), tol=1e-5, max_outer=100
):
    """Method of multipliers with semismooth Newton subproblems.

    Each outer step minimises f(u) + sum W psi_gamma(D(u) + gamma lam) by semismooth Newton
    (psi_gamma the Moreau envelope, whose gradient comes from the exact proximal step),
    then sets lam = grad psi_gamma(D(u) + gamma lam), which converges to the stress. The
    limit is the step with psi itself at any gamma; gamma starts at dt / rho and falls
    tenfold per outer step to dt / (1000 rho), which speeds the outer convergence without
    the ill-conditioning of smaller gamma. Stops at the floor once the velocity changes by
    less than tol (relative, mass norm) between outer steps.
    """
    step = Step(p, model, C)
    gamma = gamma0 if gamma0 is not None else p["dt"] / p["rho"]
    floor = floor if floor is not None else p["dt"] / (1000 * p["rho"])
    u = step.trial.copy() if u0 is None else u0[p["free"]].copy()
    lam = np.zeros((step.n, 3))
    started = perf_counter()
    history = []
    norm = lambda v: np.sqrt(v @ (step.mass * v))  # noqa: E731
    for outer in range(1, max_outer + 1):
        previous = u
        u, jam, lam, iterations = newton(step, u, gamma, lam, max_iterations=200)
        change = norm(u - previous) / norm(u)
        history.append(dict(gamma=gamma, newton=iterations, change=float(change), jammed=int(jam.sum())))
        if gamma <= floor and change <= tol:
            break
        gamma = max(gamma / reduction, floor)
    else:
        raise RuntimeError(f"method of multipliers did not settle in {max_outer} steps")
    full = p["trial"].copy()
    full[p["free"]] = u
    return (
        full,
        dict(
            outer=outer,
            newton_total=int(sum(h["newton"] for h in history)),
            jammed_fraction=float(jam.mean()),
            seconds=perf_counter() - started,
            history=history,
        ),
    )


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("cidr-frozen-step.json"))
    parser.add_argument("--spacings", type=float, nargs="*", default=[0.05, 0.025, 0.0125])
    parser.add_argument("--tol", type=float, default=1e-5)
    parser.add_argument("--solver", choices=("alg2", "newton"), default="newton")
    parser.add_argument(
        "--consolidation",
        type=float,
        default=2.0,
        help="cap pressure 2 mu_s C over the base pressure rho g H",
    )
    args = parser.parse_args(argv)
    jop = Jop()
    report = dict(
        tol=args.tol, solver=args.solver, consolidation=args.consolidation, probes=len(PROBES), cases=[]
    )
    previous = None
    for spacing in args.spacings:
        p = frozen(spacing)
        rho_s = p["rho"] / p["phi"]
        C = args.consolidation * p["rho"] * p["g"] * 1.0 / (2 * jop.mu1)
        model = CriticalState(mu_s=jop.mu1, mu_d=jop.mu2, i0=jop.I_0, k=2 * GRAIN * np.sqrt(rho_s))
        if args.solver == "alg2":
            u, _, info = solve(p, model, C, tol=args.tol)
        else:
            u, info = solve_newton(p, model, C, tol=args.tol)
        probe = sample(p, u)
        row = dict(
            spacing=spacing,
            particles=p["particles"],
            stress_nodes=len(p["solid"]),
            C=C,
            max_probe_speed=float(np.linalg.norm(probe, axis=1).max()),
            **info,
            probe=probe.tolist(),
        )
        if previous is not None:
            row["change_from_coarser"] = float(np.linalg.norm(probe - previous) / np.linalg.norm(probe))
        previous = probe
        print({k: v for k, v in row.items() if k not in ("history", "probe")}, flush=True)
        report["cases"].append(row)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
