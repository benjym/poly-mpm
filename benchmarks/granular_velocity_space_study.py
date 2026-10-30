"""Velocity-space formulation of the frozen contact problem (CPU study).

With De Saxcé's change of variable the non-dilatant Coulomb law is the fixed
point of a sequence of convex problems in the velocity u alone:

    minimise   1/2 (u - u_trial)^T (M / dt) (u - u_trial)
    subject to mu |w_t| <= w_n + s,  w = B u + gap   (every stress node)

with s <- (1 - omega) s + omega mu |w_t| updated between solves (omega = 0.5;
the undamped update oscillates), starting from the associated problem
(s = 0). For fixed s the velocity solution is unique. The cone
multipliers lambda give the contact stress r = (lambda_n, mu lambda_t). Each
convex problem is solved with the Clarabel interior-point solver; the result is
certified with the same FP64 natural residual as solve_velocity.

    uv run --locked --with clarabel python -m benchmarks.granular_velocity_space_study \\
        --output benchmarks/granular_velocity_space_study.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import sparse

TOLERANCE = 1e-4


def frozen(spacing, basis="nodal", mu=None):
    """First contact problem of the collapse from rest, with assembled and compact B."""
    from granular import assemble
    from mpm_gpu.benchmark import columns
    from transfers import nodal_velocity, scatter

    P, G, L = columns(1, spacing=spacing, particle_spacing=spacing / 4)
    P.granular_stress_basis = basis
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
    compact, free, active, _, gap = assemble(P, G, b, compact=True)
    assembled, *_ = assemble(P, G, b, compact=False)
    mass = np.repeat(G.m, 2)
    return dict(
        compact=compact,
        B=sparse.csr_matrix(assembled),
        free=free,
        mass=mass,
        trial=trial,
        gap=gap,
        mu=b.material.mu_s if mu is None else mu,
        dt=P.dt,
        particles=len(b),
    )


def natural_residual(p, contact):
    """solve_velocity's certified measure: max block natural residual, scaled as there."""
    B, gap, dt, mu = p["compact"], p["gap"], p["dt"], p["mu"]
    inverse = np.zeros_like(p["trial"])
    inverse[p["free"]] = 1 / p["mass"][p["free"]]
    inverse, trial = inverse.reshape(-1, 2), p["trial"].reshape(-1, 2)
    q = (B @ trial.ravel()).reshape(gap.shape) + gap
    magnitude = max(float(np.max(-q[:, 0], initial=0)), float(np.max(np.abs(q[:, 1:]), initial=0)))
    blocks = np.einsum("kcd,nsd->nskc", B.basis, B.gradients)
    absolute = np.zeros_like(trial)
    np.add.at(absolute, B.nodes, np.abs(blocks).sum(axis=2))
    norm = float(
        np.max(dt * np.einsum("nskc,nsc->nk", np.abs(blocks), (inverse * absolute)[B.nodes]), initial=0)
    )
    u = trial + dt * inverse * B.transpose(contact.ravel()).reshape(-1, 2)
    z = contact / (magnitude / norm)
    w = ((B @ u.ravel()).reshape(z.shape) + gap) / magnitude
    argument = z - w
    radius = mu * np.maximum(z[:, 0], 0)
    length = np.linalg.norm(argument[:, 1:], axis=1)
    projected = argument.copy()
    projected[:, 0] = np.maximum(argument[:, 0], 0)
    projected[:, 1:] *= np.minimum(1, np.divide(radius, length, out=np.ones(len(z)), where=length > 0))[
        :, None
    ]
    return float(np.max(np.linalg.norm(z - projected, axis=1), initial=0)), u


def conic(p, shift, *, tol=1e-10):
    """One convex problem in the free velocity components; returns contact stress and solver info."""
    import clarabel

    B, free, dt, mu = p["B"], p["free"], p["dt"], p["mu"]
    rows = len(p["gap"])
    Bf = B[:, free].tocsr()
    fixed_part = B[:, ~free] @ p["trial"][~free]
    mass = p["mass"][free] / dt
    P = sparse.diags(mass).tocsc()
    q = -mass * p["trial"][free]
    # Cone variable s = (w_n + shift, mu w_t) in the second-order cone, s = b - A u.
    weight = np.tile([1.0, mu, mu], rows)
    A = -sparse.diags(weight) @ Bf
    b = weight * ((fixed_part + p["gap"].ravel()) + np.repeat(shift[:, None] * [1, 0, 0], 1, axis=0).ravel())
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.tol_gap_abs = settings.tol_gap_rel = tol
    settings.tol_feas = tol
    settings.max_iter = 200
    solver = clarabel.DefaultSolver(P, q, A.tocsc(), b, [clarabel.SecondOrderConeT(3)] * rows, settings)
    solution = solver.solve()
    lam = np.asarray(solution.z).reshape(-1, 3)
    contact = np.column_stack((lam[:, 0], mu * lam[:, 1:]))
    u = p["trial"].copy()
    u[free] = np.asarray(solution.x)
    return contact, u, dict(status=str(solution.status), iterations=int(solution.iterations))


def fixed_point(p, *, max_outer=100, omega=0.5):
    """De Saxcé fixed point on the shift s = mu |w_t|, certified by the natural residual."""
    rows = len(p["gap"])
    shift = np.zeros(rows)
    history = []
    started = perf_counter()
    for outer in range(max_outer):
        contact, u, info = conic(p, shift)
        residual, _ = natural_residual(p, contact)
        w = (p["B"] @ u).reshape(-1, 3) + p["gap"]
        new = p["mu"] * np.linalg.norm(w[:, 1:], axis=1)
        change = float(np.max(np.abs(new - shift)) / max(np.max(np.abs(new)), 1e-300))
        history.append(
            dict(
                outer=outer,
                residual=residual,
                ipm_iterations=info["iterations"],
                status=info["status"],
                shift_change=change,
            )
        )
        if residual <= TOLERANCE:
            break
        shift = (1 - omega) * shift + omega * new
    return dict(
        converged=history[-1]["residual"] <= TOLERANCE,
        outer=len(history),
        omega=omega,
        ipm_iterations=sum(h["ipm_iterations"] for h in history),
        residual=history[-1]["residual"],
        seconds=perf_counter() - started,
        history=history,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("granular-velocity-space-study.json"))
    parser.add_argument("--spacings", type=float, nargs="*", default=[0.05, 0.025, 0.0125])
    args = parser.parse_args(argv)
    report = dict(tolerance=TOLERANCE, cases=[])
    for spacing in args.spacings:
        p = frozen(spacing)
        result = fixed_point(p)
        row = dict(
            spacing=spacing,
            particles=p["particles"],
            stress_unknowns=3 * len(p["gap"]),
            free_velocity=int(p["free"].sum()),
            **result,
        )
        print({k: v for k, v in row.items() if k != "history"}, flush=True)
        report["cases"].append(row)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
