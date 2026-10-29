"""Does dilatancy regularise the refined frictional contact problem? (CPU study)

Sliding with dilatancy beta separates at w_n = beta |w_t| while pressed
(beta = 0 is the current non-dilatant law, beta = mu associated Coulomb flow).
Normal complementarity applies to w_n - beta |w_t|; the friction disk and the
sliding direction are unchanged. The local solve in the cached diagonal metric
reduces to the monotone scalar equation

    f(alpha) = |r_t(alpha)| (D_n - mu beta alpha) + mu q_n = 0,
    r_t(alpha) = -q_t / (D_t + alpha),  r_n = |r_t| / mu,

which this prototype brackets and bisects. The global iteration is
solve_velocity's colored sweeps and safeguarded Anderson acceleration (no
Newton polish), on the frozen first-step problems of
benchmarks.granular_stress_scaling.

    uv run --locked python -m benchmarks.granular_dilatancy_study --output benchmarks/granular_dilatancy_study.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from benchmarks.granular_stress_scaling import problem

TOLERANCE = 1e-4


def local_dilatant(D, q, mu, beta):
    """Exact local solve with dilatancy in a diagonal metric (rows of D, q)."""
    r = np.zeros_like(q)
    qn, qt, Dn, Dt = q[:, 0], q[:, 1:], D[:, 0], D[:, 1:]
    size = np.linalg.norm(qt, axis=1)
    opening = qn - beta * size >= 0
    stick_n = -qn / Dn
    stick_t = -qt / Dt
    stick = ~opening & (stick_n > 0) & (np.linalg.norm(stick_t, axis=1) <= mu * stick_n)
    r[stick, 0], r[stick, 1:] = stick_n[stick], stick_t[stick]
    slide = ~opening & ~stick & (mu > 0)
    if np.any(slide):
        d_n, d_t, b_t, b_n = Dn[slide], Dt[slide], qt[slide], qn[slide]

        def f(alpha):
            rt = np.linalg.norm(b_t / (d_t + alpha[:, None]), axis=1)
            return rt * (d_n - mu * beta * alpha) + mu * b_n

        lower = np.zeros(len(d_n))
        upper = np.maximum(d_t.max(axis=1), d_n / max(mu * beta, 1e-300) if beta > 0 else d_t.max(axis=1))
        for _ in range(200):
            grow = f(upper) > 0
            if not np.any(grow):
                break
            upper[grow] *= 2
        for _ in range(100):
            middle = 0.5 * (lower + upper)
            positive = f(middle) > 0
            lower = np.where(positive, middle, lower)
            upper = np.where(positive, upper, middle)
        alpha = 0.5 * (lower + upper)
        rt = -b_t / (d_t + alpha[:, None])
        r[slide, 1:] = rt
        r[slide, 0] = np.linalg.norm(rt, axis=1) / mu
    return r


def solve(p, beta, *, max_sweeps=3000, anderson=True, initial=None):
    """solve_velocity's sweeps and Anderson acceleration with the dilatant law."""
    B, gap, dt, mu = p["B"], p["gap"], p["dt"], p["mu"]
    trial, inverse = p["trial"].reshape(-1, 2), p["inverse"].reshape(-1, 2)
    q = (B @ trial.ravel()).reshape(gap.shape) + gap
    magnitude = max(float(np.max(-q[:, 0], initial=0)), float(np.max(np.abs(q[:, 1:]), initial=0)))
    blocks = np.einsum("kcd,nsd->nskc", B.basis, B.gradients)
    local = dt * np.einsum("nskc,nsc,nslc->nkl", blocks, inverse[B.nodes], blocks)
    absolute = np.zeros_like(trial)
    np.add.at(absolute, B.nodes, np.abs(blocks).sum(axis=2))
    norm = float(
        np.max(dt * np.einsum("nskc,nsc->nk", np.abs(blocks), (inverse * absolute)[B.nodes]), initial=0)
    )
    scale = magnitude / norm
    local /= norm
    regularizer = 1e-6 * np.maximum(np.trace(local, axis1=1, axis2=2), 1e-12)
    eigenvalues, rotation = np.linalg.eigh(local[:, 1:, 1:])
    diagonal = np.column_stack((local[:, 0, 0], eigenvalues)) + regularizer[:, None]
    Q = np.zeros((len(q), 3, 3))
    Q[:, 0, 0], Q[:, 1:, 1:] = 1, rotation
    colors = B.colors(inverse)
    z = np.zeros_like(q) if initial is None else initial / scale

    def velocity():
        return trial + dt * scale * inverse * B.transpose(z.ravel()).reshape(-1, 2)

    def residual(u, z):
        w = ((B @ u.ravel()).reshape(q.shape) + gap) / magnitude
        bar = w.copy()
        bar[:, 0] -= beta * np.linalg.norm(w[:, 1:], axis=1)
        argument = z - bar
        projected = argument.copy()
        projected[:, 0] = np.maximum(argument[:, 0], 0)
        radius = mu * np.maximum(z[:, 0], 0)
        length = np.linalg.norm(argument[:, 1:], axis=1)
        projected[:, 1:] *= np.minimum(1, np.divide(radius, length, out=np.ones(len(z)), where=length > 0))[
            :, None
        ]
        return float(np.max(np.linalg.norm(z - projected, axis=1), initial=0))

    def feasible(value):
        value = value.copy()
        value[:, 0] = np.maximum(value[:, 0], 0)
        length = np.linalg.norm(value[:, 1:], axis=1)
        value[:, 1:] *= np.minimum(
            1, np.divide(mu * value[:, 0], length, out=np.ones(len(value)), where=length > 0)
        )[:, None]
        return value

    u = velocity()
    error = residual(u, z)
    trace = [(0, error)]
    history = []
    sweeps = 0
    while error > TOLERANCE and sweeps < max_sweeps:
        start = z.copy()
        for group in colors:
            nodes, gradients = B.nodes[group], B.gradients[group]
            w = (np.einsum("kcd,nsd,nsc->nk", B.basis, gradients, u[nodes]) + gap[group]) / magnitude
            old = np.einsum("nki,nk->ni", Q[group], z[group])
            other = np.einsum("nki,nk->ni", Q[group], w) - diagonal[group] * old
            new = local_dilatant(diagonal[group], other, mu, beta)
            delta = np.einsum("nki,ni->nk", Q[group], new - old)
            z[group] = np.einsum("nki,ni->nk", Q[group], new)
            force = np.einsum("kcd,nsd,nk->nsc", B.basis, gradients, delta)
            np.add.at(u, nodes, dt * scale * inverse[nodes] * force)
        sweeps += 1
        u = velocity()
        error = residual(u, z)
        history.append((start.ravel(), (z - start).ravel()))
        history = history[-13:]
        if anderson and len(history) > 1 and error > TOLERANCE:
            dx = np.column_stack([b[0] - a[0] for a, b in zip(history[:-1], history[1:])])
            df = np.column_stack([b[1] - a[1] for a, b in zip(history[:-1], history[1:])])
            c = np.linalg.lstsq(df, history[-1][1], rcond=1e-10)[0]
            candidate = feasible((z.ravel() - (dx + df) @ c).reshape(z.shape))
            saved, z = z, candidate
            candidate_u = velocity()
            candidate_error = residual(candidate_u, z)
            if np.isfinite(candidate_error) and candidate_error < error:
                u, error = candidate_u, candidate_error
            else:
                z = saved
        trace.append((sweeps, error))
    return z * scale, dict(converged=error <= TOLERANCE, residual=error, sweeps=sweeps, trace=trace)


def local_check(count=2000, seed=4):
    """Local solves satisfy the dilatant conditions; beta = 0 matches diagonal_coulomb."""
    from granular_mixed import diagonal_coulomb

    rng = np.random.default_rng(seed)
    D = 10 ** rng.uniform(-2, 2, (count, 3))
    q = rng.normal(size=(count, 3))
    mu = 0.5
    worst = {}
    for beta in (0.0, 0.25, 0.5):
        r = local_dilatant(D, q, mu, beta)
        w = D * r + q
        bar = w[:, 0] - beta * np.linalg.norm(w[:, 1:], axis=1)
        pressed = r[:, 0] > 1e-12
        cone = float(np.max(np.linalg.norm(r[:, 1:], axis=1) - mu * r[:, 0]))
        normal = float(np.max(np.abs(bar[pressed])))
        opening = float(np.min(bar[~pressed], initial=0))
        sliding = pressed & (np.linalg.norm(w[:, 1:], axis=1) > 1e-9)
        cosine = np.sum(r[sliding, 1:] * w[sliding, 1:], axis=1) / (
            np.linalg.norm(r[sliding, 1:], axis=1) * np.linalg.norm(w[sliding, 1:], axis=1)
        )
        worst[beta] = dict(
            cone=cone, normal=normal, opening=opening, direction=float(np.max(1 + cosine, initial=0))
        )
        if beta == 0:
            worst[beta]["against_diagonal_coulomb"] = float(np.abs(r - diagonal_coulomb(D, q, mu)).max())
    return worst


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("granular-dilatancy-study.json"))
    args = parser.parse_args(argv)
    report = dict(tolerance=TOLERANCE, local=local_check(), cases=[])
    for spacing in (0.05, 0.025, 0.0125):
        p = problem(spacing)
        for ratio in (0.0, 0.1, 0.25, 0.5, 1.0):
            started = perf_counter()
            _, result = solve(p, ratio * p["mu"])
            marks = {
                f"{k:g}": next((s for s, e in result["trace"] if e <= k), None) for k in (1e-2, 1e-3, 1e-4)
            }
            row = dict(
                spacing=spacing,
                particles=p["particles"],
                stress_unknowns=3 * p["B"].count,
                beta_over_mu=ratio,
                converged=bool(result["converged"]),
                residual=result["residual"],
                sweeps=result["sweeps"],
                first_sweep_below=marks,
                seconds=perf_counter() - started,
            )
            print(row, flush=True)
            report["cases"].append(row)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
