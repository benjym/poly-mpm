"""Model or solver: refinement of the frozen first collapse step under regularised laws.

The rigid-grain law (constant mu, non-dilatant, rate-independent) stalls every
solver on refined connected grids; docs/connected-flows-literature.md argues
that the law itself is ill-posed. This study separates the two by solving the
first contact problem of the column collapse from rest at successively finer
grids under four laws, in velocity form with the Clarabel interior-point
solver (benchmarks/granular_velocity_space_study.py):

    minimise   1/2 |u - u_trial|^2_(M/dt) + sum_a (eta / V_a) |w_dev,a|^2
    subject to mu |w_t| <= w_n + s,   w = B u + gap   (every stress node)

- associated:            s = 0, eta = 0 (normality; one convex solve)
- associated, viscous:   s = 0, eta > 0 (convex stand-in for a well-posed
                         dilatant, rate-dependent critical-state law)
- nonassociated:         s = mu |w_t| by damped fixed point, eta = 0 (current law)
- nonassociated, viscous: as above with eta > 0 (rate dependence alone)

The viscous term is the shear dissipation integral(2 eta D:D) with D_a = w_a / V_a
in the orthonormal stress basis. For each law it records whether the solve
converges and how the solution changes under refinement: nodal velocities are
sampled at fixed probe points in the column and successive grids compared.
A well-posed law should give differences that fall with the grid spacing.

    python -m benchmarks.granular_regularisation_study --output benchmarks/granular_regularisation_study.json
(Clarabel is required: pip install clarabel, or uv run --with clarabel.)
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import sparse

TOLERANCE = 1e-4
PROBES = np.stack(np.meshgrid(np.arange(0.02, 0.49, 0.02), np.arange(0.02, 0.99, 0.02)), axis=-1).reshape(
    -1, 2
)


def frozen(spacing):
    """First contact problem of the collapse from rest: operators, masses and grid."""
    from granular import assemble, basis_volumes
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
    compact, free, active, _, gap = assemble(P, G, b, compact=True)
    assembled, *_ = assemble(P, G, b, compact=False)
    return dict(
        compact=compact,
        B=sparse.csr_matrix(assembled),
        free=free,
        mass=np.repeat(G.m, 2),
        trial=trial,
        gap=gap,
        volume=basis_volumes(P, G)[active],
        mu=b.material.mu,
        dt=P.dt,
        g=abs(P.g),
        grid=(P.G.x_m, P.G.y_m, G.dx, G.dy, P.G.nx, P.G.ny),
        particles=len(b),
    )


def solve(p, shift, eta, *, tol=1e-10):
    """One convex problem in the free velocity components."""
    import clarabel

    B, free, dt, mu = p["B"], p["free"], p["dt"], p["mu"]
    rows = len(p["gap"])
    Bf = B[:, free].tocsr()
    fixed = B[:, ~free] @ p["trial"][~free]
    mass = p["mass"][free] / dt
    H = sparse.diags(mass).tocsc()
    q = -mass * p["trial"][free]
    if eta:
        deviator = np.flatnonzero(np.tile([False, True, True], rows))
        weight = np.repeat(2 * eta / p["volume"], 2)
        Bd = Bf[deviator]
        K = sparse.diags(weight)
        H = (H + 2 * (Bd.T @ K @ Bd)).tocsc()
        q = q + 2 * (Bd.T @ (weight * fixed[deviator]))
    scale = np.tile([1.0, mu, mu], rows)
    A = -sparse.diags(scale) @ Bf
    offset = np.zeros((rows, 3))
    offset[:, 0] = shift
    b = scale * (fixed + p["gap"].ravel() + offset.ravel())
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.tol_gap_abs = settings.tol_gap_rel = tol
    settings.tol_feas = tol
    settings.max_iter = 200
    solution = clarabel.DefaultSolver(
        H, q, A.tocsc(), b, [clarabel.SecondOrderConeT(3)] * rows, settings
    ).solve()
    lam = np.asarray(solution.z).reshape(-1, 3)
    contact = np.column_stack((lam[:, 0], mu * lam[:, 1:]))
    u = p["trial"].copy()
    u[free] = np.asarray(solution.x)
    return contact, u, dict(status=str(solution.status), iterations=int(solution.iterations))


def residual(p, contact, u, associated):
    """Natural residual of the frictional contact at velocity u, scaled as in solve_velocity."""
    B, gap, dt, mu = p["compact"], p["gap"], p["dt"], p["mu"]
    trial = p["trial"].reshape(-1, 2)
    inverse = np.zeros_like(p["trial"])
    inverse[p["free"]] = 1 / p["mass"][p["free"]]
    inverse = inverse.reshape(-1, 2)
    q = (B @ trial.ravel()).reshape(gap.shape) + gap
    magnitude = max(float(np.max(-q[:, 0], initial=0)), float(np.max(np.abs(q[:, 1:]), initial=0)))
    blocks = np.einsum("kcd,nsd->nskc", B.basis, B.gradients)
    absolute = np.zeros_like(trial)
    np.add.at(absolute, B.nodes, np.abs(blocks).sum(axis=2))
    norm = float(
        np.max(dt * np.einsum("nskc,nsc->nk", np.abs(blocks), (inverse * absolute)[B.nodes]), initial=0)
    )
    z = contact / (magnitude / norm)
    w = ((B @ u).reshape(z.shape) + gap) / magnitude
    if associated:
        # Normality: w + mu-scaled shift is the dual cone of the stress cone.
        w = w.copy()
        w[:, 0] -= mu * np.linalg.norm(w[:, 1:], axis=1)
    argument = z - w
    radius = mu * np.maximum(z[:, 0], 0)
    length = np.linalg.norm(argument[:, 1:], axis=1)
    projected = argument.copy()
    projected[:, 0] = np.maximum(argument[:, 0], 0)
    projected[:, 1:] *= np.minimum(1, np.divide(radius, length, out=np.ones(len(z)), where=length > 0))[
        :, None
    ]
    return float(np.max(np.linalg.norm(z - projected, axis=1), initial=0))


def law(p, associated, eta, *, max_outer=60, omega=0.5):
    rows = len(p["gap"])
    shift = np.zeros(rows)
    history = []
    started = perf_counter()
    for outer in range(1 if associated else max_outer):
        contact, u, info = solve(p, shift, eta)
        error = residual(p, contact, u, associated)
        history.append(dict(residual=error, ipm_iterations=info["iterations"], status=info["status"]))
        if associated or error <= TOLERANCE:
            break
        w = (p["B"] @ u).reshape(-1, 3) + p["gap"]
        shift = (1 - omega) * shift + omega * p["mu"] * np.linalg.norm(w[:, 1:], axis=1)
    return u, dict(
        converged=history[-1]["residual"] <= TOLERANCE or (associated and history[-1]["status"] == "Solved"),
        outer=len(history),
        ipm_iterations=[h["ipm_iterations"] for h in history],
        status=history[-1]["status"],
        residual=history[-1]["residual"],
        seconds=perf_counter() - started,
    )


def sample(p, u):
    """Bilinear interpolation of nodal velocity at PROBES, in units of g dt."""
    x_m, y_m, dx, dy, nx, ny = p["grid"]
    v = u.reshape(ny, nx, 2)
    fx, fy = (PROBES[:, 0] - x_m) / dx, (PROBES[:, 1] - y_m) / dy
    i, j = np.floor(fx).astype(int), np.floor(fy).astype(int)
    a, c = (fx - i)[:, None], (fy - j)[:, None]
    value = (
        (1 - a) * (1 - c) * v[j, i]
        + a * (1 - c) * v[j, i + 1]
        + a * c * v[j + 1, i + 1]
        + (1 - a) * c * v[j + 1, i]
    )
    return value / (p["g"] * p["dt"])


LAWS = {
    "associated": (True, 0.0),
    "associated_viscous": (True, None),
    "nonassociated": (False, 0.0),
    "nonassociated_viscous": (False, None),
}


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("granular-regularisation-study.json"))
    parser.add_argument("--spacings", type=float, nargs="*", default=[0.05, 0.025, 0.0125, 0.00625])
    parser.add_argument("--eta", type=float, default=100.0, help="viscosity of the viscous laws (Pa s)")
    parser.add_argument("--laws", nargs="*", default=list(LAWS))
    args = parser.parse_args(argv)
    report = dict(tolerance=TOLERANCE, eta=args.eta, probes=len(PROBES), cases=[])
    samples = {}
    for spacing in args.spacings:
        p = frozen(spacing)
        for name in args.laws:
            associated, eta = LAWS[name]
            eta = args.eta if eta is None else eta
            u, result = law(p, associated, eta)
            probe = sample(p, u)
            row = dict(
                law=name,
                eta=eta,
                spacing=spacing,
                particles=p["particles"],
                stress_rows=len(p["gap"]),
                kinetic_energy=float(0.5 * np.sum(p["mass"] * u**2)),
                max_probe_speed=float(np.linalg.norm(probe, axis=1).max()),
                **result,
            )
            previous = samples.get(name)
            if previous is not None:
                row["change_from_coarser"] = float(np.linalg.norm(probe - previous) / np.linalg.norm(probe))
            samples[name] = probe
            print({k: v for k, v in row.items() if k != "ipm_iterations"}, flush=True)
            report["cases"].append(row)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
