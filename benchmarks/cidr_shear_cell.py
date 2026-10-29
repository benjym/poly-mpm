"""One-dimensional gravity-free shear cell with the convex CIDR time step.

The test of Schaeffer et al. (2019, sections 3 and 5): granular material between
a fixed bottom wall and a top wall moving at unit speed, no slip, invariant in
the shearing direction. Variables are non-dimensional (grain diameter, grain
density and wall speed are 1). Each cell keeps its mass and the nodes move with
the vertical velocity (Lagrangian), so there is no advection term. Each time step
solves the convex problem of docs/granular-cidr-design.md,

    min_v  sum_i M_i / (2 dt) |v_i - v_i^n|^2 + sum_c dz_c psi(s_c, delta_c; phi_c),

with v = (u, w) at the nodes, s = |D| = sqrt(u'^2 + w'^2) / 2 and delta = w', by
Newton's method with a backtracking line search on a block-tridiagonal system.
phi is held at its start-of-step value ("lagged") or iterated to its end-of-step
value ("implicit").

Models: inertial CIDR (Schaeffer et al. 2019) and critical-state CIDR (Barker et
al. 2017) completed at the cap, with C(phi) from a tunable packing deficit.

    python -m benchmarks.cidr_shear_cell --output benchmarks/cidr_shear_cell.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import integrate
from scipy.interpolate import CubicSpline
from scipy.linalg import solve_banded

# Schaeffer et al. (2019), table 1: DEM fit of mu(I) and Phi(I).
MU_S, MU_D, I0, PHI_C, A_PHI = 0.248, 0.588, 0.145, 0.817, 0.314
K = 2.0  # I = K s / sqrt(p) with d = rho_s = 1


def mu(I):
    return MU_S + (MU_D - MU_S) * I / (I0 + I)


def dmu(I):
    return (MU_D - MU_S) * I0 / (I0 + I) ** 2


class Inertial:
    """iCIDR: g = 2 M s^2 sqrt(p) / Psi - M Psi p^(3/2) / 6, M = mu(Psi), Psi = (phi_c - phi) / a."""

    name = "inertial"

    def evaluate(self, s, delta, phi):
        Psi = np.maximum((PHI_C - phi) / A_PHI, 1e-12)
        M = mu(Psi)
        # Optimal pressure, eq. (4.12), non-dimensional (rho = d = 1).
        p = (2 / (M * Psi)) ** 2 * (np.sqrt(delta**2 + M**2 * s**2) - delta) ** 2
        p = np.maximum(p, 1e-300)
        r = np.sqrt(p)
        g = 2 * M * s**2 * r / Psi - M * Psi * p * r / 6
        g_s = 4 * M * s * r / Psi
        g_ss = 4 * M * r / Psi
        g_sp = 2 * M * s / (Psi * r)
        g_pp = -M * s**2 / (2 * Psi * p * r) - M * Psi / (8 * r)
        psi = g - p * delta
        return _interior(psi, g_s, -p, g_ss, g_sp, g_pp), p


class CriticalState:
    """Critical-state CIDR with Barker et al.'s alpha, beta and C = (P/b) x, x = psi / (1 - psi).

    P/b follows from a packing deficit e_ref on the cap at pressure p_ref:
    P/b = p_ref / (2 mu(0) (Delta / e_ref - 1)), Delta = phi_max - phi_min.
    """

    name = "critical_state"

    def __init__(self, deficit=0.005, p_ref=0.08, phi_min=0.05, phi_max=PHI_C):
        self.phi_min, self.phi_max = phi_min, phi_max
        self.Delta = phi_max - phi_min
        self.scale = p_ref / (2 * mu(0.0) * (self.Delta / deficit - 1))
        grid = np.logspace(-12, 3, 600)
        alpha = [
            0.8 * mu(I) + 0.48 * integrate.quad(lambda t: t ** (-0.6) * mu(t * I), 0, 1)[0] for I in grid
        ]
        self.spline = CubicSpline(np.log(grid), np.array(alpha))

    def C(self, phi):
        q = np.clip((phi - self.phi_min) / self.Delta, 1e-12, 1 - 1e-12)
        return self.scale * q / (1 - q)

    def alpha(self, I):
        return self.spline(np.log(np.clip(I, 1e-12, 1e3)))

    def functions(self, s, p, C):
        I = K * s / np.sqrt(p)
        a, m, dm = self.alpha(I), mu(I), dmu(I)
        da = np.where(I > 1e-10, (2 * m + 2 * I * dm - a) / (2.5 * np.maximum(I, 1e-300)), 0.0)
        b = 2 * (a - m)
        db = 2 * (da - dm)
        return I, a, da, b, db

    def p_cap(self, s, C):
        p = 2 * mu(0.0) * C * np.ones_like(s)
        for _ in range(60):
            I, a, da, _, _ = self.functions(s, p, C)
            h = p - a * C
            dh = 1 + da * C * I / (2 * p)
            p = np.maximum(p - h / dh, 0.5 * p)
        return p

    def evaluate(self, s, delta, phi):
        C = self.C(phi)
        cap = self.p_cap(s, C)

        def gp(p):
            I, a, da, b, db = self.functions(s, p, C)
            return 2 * s * (b - 2 * p / C)

        # Bracketed bisection-Newton for g_p(s, p) = delta on (0, cap); g_p decreases in p.
        lo, hi = np.full_like(s, 1e-300), cap.copy()
        at_cap = gp(cap) - delta >= 0
        p = 0.5 * (lo + hi)
        for _ in range(80):
            value = gp(p) - delta
            lo = np.where(value > 0, p, lo)
            hi = np.where(value <= 0, p, hi)
            p = 0.5 * (lo + hi)
        p = np.where(at_cap, cap, p)
        I, a, da, b, db = self.functions(s, p, C)
        g = 2 * p * s * ((5 * a - 4 * mu(I)) / 3 - p / C)
        Y = a * p - p**2 / C
        g_s, g_p = 2 * Y, 2 * s * (b - 2 * p / C)
        g_ss = 2 * da * p * I / np.maximum(s, 1e-300)
        g_sp = 2 * (a - 2 * p / C) - da * I
        g_pp = 2 * s * (-2 / C - I * db / (2 * p))
        psi = g - p * delta
        interior = _interior(psi, g_s, -p, g_ss, g_sp, g_pp)
        # Cap branch: psi_s = (g_p - delta) p_cap', psi_delta = -p_cap, psi_deltadelta = 0.
        dcap = da * C * I / np.maximum(s, 1e-300) / (1 + da * C * I / (2 * p))
        cap_s = (g_p - delta) * dcap
        cap_ss = (g_sp + g_pp * dcap) * dcap
        out = tuple(
            np.where(at_cap, c, i)
            for i, c in zip(interior, (psi, cap_s, -p, np.maximum(cap_ss, 0), -dcap, np.zeros_like(p)))
        )
        return out, p


def _interior(psi, g_s, psi_d, g_ss, g_sp, g_pp):
    return (psi, g_s, psi_d, g_ss - g_sp**2 / g_pp, g_sp / g_pp, -1 / g_pp)


def cell_terms(model, a, b, phi):
    """psi per unit volume of each cell as a function of (u', w') = (a, b): value, gradient, Hessian."""
    r = np.maximum(np.hypot(a, b), 1e-14)
    s = 0.5 * r
    (psi, ps, pd, pss, psd, pdd), p = model.evaluate(s, b, phi)
    sa, sb = a / (2 * r), b / (2 * r)
    saa, sbb, sab = b**2 / (2 * r**3), a**2 / (2 * r**3), -a * b / (2 * r**3)
    grad = np.stack((ps * sa, ps * sb + pd), axis=-1)
    Haa = pss * sa**2 + ps * saa
    Hab = pss * sa * sb + psd * sa + ps * sab
    Hbb = pss * sb**2 + 2 * psd * sb + pdd + ps * sbb
    return psi, grad, np.stack((np.stack((Haa, Hab), -1), np.stack((Hab, Hbb), -1)), -2), p


class Cell:
    def __init__(self, model, n, height, phi0, u0, w0):
        self.model = model
        self.z = np.linspace(0.0, height, n + 1)
        dz = np.diff(self.z)
        self.mass = phi0(0.5 * (self.z[1:] + self.z[:-1])) * dz
        self.u, self.w = u0(self.z), w0(self.z)
        self.u[0], self.u[-1], self.w[0], self.w[-1] = 0.0, 1.0, 0.0, 0.0

    @property
    def phi(self):
        return self.mass / np.diff(self.z)

    def objective(self, v, v0, nodal_mass, dt, dz, phi):
        u, w = v[0::2], v[1::2]
        a, b = np.diff(u) / dz, np.diff(w) / dz
        psi, grad, hess, p = cell_terms(self.model, a, b, phi)
        f = 0.5 * np.sum(nodal_mass[:, None] * (v - v0).reshape(-1, 2) ** 2) / dt + np.sum(dz * psi)
        g = (nodal_mass[:, None] * (v - v0).reshape(-1, 2) / dt).copy()
        # d(a)/d(u_i) = -1/dz, d(a)/d(u_{i+1}) = 1/dz, likewise for b.
        flux = grad  # per cell, times dz / dz
        g[:-1] -= flux
        g[1:] += flux
        return f, g.ravel(), hess, p

    def step(self, dt, phi):
        """One Newton-solved convex step at fixed phi; returns iterations and pressure."""
        dz = np.diff(self.z)
        n = len(self.z)
        nodal_mass = np.zeros(n)
        nodal_mass[:-1] += 0.5 * self.mass
        nodal_mass[1:] += 0.5 * self.mass
        v0 = np.column_stack((self.u, self.w)).ravel()
        v = v0.copy()
        scale = np.sum(np.repeat(nodal_mass, 2) * v0**2) / dt
        free = np.ones(2 * n, dtype=bool)
        free[[0, 1, -2, -1]] = False
        for iteration in range(100):
            f, g, hess, p = self.objective(v, v0, nodal_mass, dt, dz, phi)
            # Block-tridiagonal Hessian: mass / dt plus per-cell (1/dz) [H, -H; -H, H].
            H = hess / dz[:, None, None]
            size = 2 * n
            full = np.zeros((size, size))
            idx = np.arange(n)
            full[2 * idx, 2 * idx] += nodal_mass / dt
            full[2 * idx + 1, 2 * idx + 1] += nodal_mass / dt
            for c in range(n - 1):
                i, j = 2 * c, 2 * c + 2
                full[i : i + 2, i : i + 2] += H[c]
                full[j : j + 2, j : j + 2] += H[c]
                full[i : i + 2, j : j + 2] -= H[c]
                full[j : j + 2, i : i + 2] -= H[c]
            A = full[np.ix_(free, free)]
            gf = g[free]
            # Banded solve (bandwidth 3 on each side for 2-component nodes).
            band = 3
            ab = np.zeros((2 * band + 1, len(gf)))
            for k in range(-band, band + 1):
                diag = np.diagonal(A, k)
                if k >= 0:
                    ab[band - k, k:] = diag
                else:
                    ab[band - k, : len(diag)] = diag
            direction = np.zeros_like(v)
            direction[free] = -solve_banded((band, band), ab, gf)
            slope = g @ direction
            # Newton decrement against the step's kinetic-energy scale, which does
            # not depend on the size of psi's values.
            if -slope <= 1e-24 + 1e-15 * scale:
                break
            t = 1.0
            while True:
                trial = v + t * direction
                ft = self.objective(trial, v0, nodal_mass, dt, dz, phi)[0]
                if ft <= f + 1e-4 * t * slope or t < 1e-12:
                    break
                t *= 0.5
            v = trial
        else:
            # Plain Newton cannot pass the cone of psi at D = 0 (static yield).
            raise RuntimeError(f"Newton did not converge: decrement {-slope:.3g}, scale {scale:.3g}")
        self.u, self.w = v[0::2].copy(), v[1::2].copy()
        return iteration, p

    def advance(self, dt, implicit=False):
        z0, u0, w0 = self.z.copy(), self.u.copy(), self.w.copy()
        phi = self.phi
        for _ in range(30 if implicit else 1):
            # Each pass restarts from the start-of-step state; only phi changes.
            self.z, self.u, self.w = z0, u0.copy(), w0.copy()
            iterations, p = self.step(dt, phi)
            z_new = z0 + dt * self.w
            phi_new = self.mass / np.diff(z_new)
            change = np.abs(phi_new - phi).max()
            phi = phi_new
            if not implicit or change < 1e-12:
                break
        self.z = z0 + dt * self.w
        return iterations, p


def run(model, n, case, t_end, dt0=1e-3, dt_max=0.25, growth=1.02, implicit=False, outputs=()):
    if case == "A":  # Schaeffer et al. (2019), (3.16), figure 6
        H, phi0, eps = 30.0, 0.78, 0.05
        cell = Cell(
            model,
            n,
            H,
            lambda z: np.full_like(z, phi0),
            lambda z: z / H,
            lambda z: eps * (z - H / 2) * np.exp(-((z - H / 2) ** 2)),
        )
    else:  # (5.7), figure 7
        H, phi0, a, b = 100.0, 0.78, 0.16, 0.1
        cell = Cell(
            model,
            n,
            H,
            lambda z: np.full_like(z, phi0),
            lambda z: z / H + a * np.sin(2 * np.pi * z / H),
            lambda z: -b * np.sin(2 * np.pi * z / H),
        )
    t, dt = 0.0, dt0
    history = [(0.0, float(np.abs(cell.w).max()))]
    snapshots = {}
    pending = sorted(outputs)
    started = perf_counter()
    newton = []
    while t < t_end - 1e-12:
        step = min(dt, t_end - t)
        limit = getattr(model, "phi_max", None)
        if limit is not None:
            # With phi lagged, a step must not compact a cell past phi_max: allow at
            # most half the remaining gap, using the current compaction rate.
            rate = -np.diff(cell.w) / np.diff(cell.z)
            gap = (limit - cell.phi) / cell.phi
            closing = rate > 0
            if np.any(closing):
                step = min(step, 0.5 * float(np.min(gap[closing] / rate[closing])))
        if pending and t + step > pending[0]:
            step = pending[0] - t
        iterations, p = cell.advance(step, implicit)
        newton.append(iterations)
        t += step
        history.append((t, float(np.abs(cell.w).max())))
        if pending and abs(t - pending[0]) < 1e-12:
            snapshots[pending.pop(0)] = dict(
                z=cell.z.tolist(), u=cell.u.tolist(), w=cell.w.tolist(), phi=cell.phi.tolist(), p=p.tolist()
            )
        dt = min(dt * growth, dt_max)
    steady_u = np.abs(cell.u - cell.z / H).max()
    return dict(
        model=model.name,
        cells=n,
        case=case,
        implicit=implicit,
        steps=len(newton),
        newton_mean=float(np.mean(newton)),
        newton_max=int(np.max(newton)),
        seconds=perf_counter() - started,
        history=history,
        snapshots=snapshots,
        final=dict(
            max_w=float(np.abs(cell.w).max()),
            u_minus_linear=float(steady_u),
            phi_range=[float(cell.phi.min()), float(cell.phi.max())],
            mean_phi=float(np.sum(cell.mass) / H),
        ),
    )


def compare(a, b, times):
    """Largest difference of max|w| between two runs at given times (interpolated in log time)."""
    ha, hb = np.array(a["history"]), np.array(b["history"])
    out = []
    for t in times:
        va = np.interp(t, ha[:, 0], ha[:, 1])
        vb = np.interp(t, hb[:, 0], hb[:, 1])
        out.append(abs(va - vb))
    return float(max(out))


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=Path("cidr-shear-cell.json"))
    parser.add_argument("--deficit", type=float, default=0.005)
    parser.add_argument("--t-end-b", type=float, default=1000.0)
    parser.add_argument("--models", nargs="+", default=["inertial", "critical_state"])
    parser.add_argument(
        "--lag-only", action="store_true", help="Only rerun the packing-lag comparison, into --output"
    )
    args = parser.parse_args(argv)
    report = dict(deficit=args.deficit, grid=[], packing_lag=[], cases=[])
    if args.lag_only and args.output.exists():
        report = json.loads(args.output.read_text())
        report["packing_lag"] = []
        report["cases"] = [c for c in report["cases"] if c.get("role") == "grid"]
    models = [m for m in (Inertial(), CriticalState(deficit=args.deficit)) if m.name in args.models]
    plan = [("A", 20.0), ("B", args.t_end_b)]

    def keep(result, **extra):
        thin = result["history"][:: max(1, len(result["history"]) // 400)]
        report["cases"].append(
            {k: v for k, v in result.items() if k != "history"} | dict(history=thin) | extra
        )

    for model in models:
        for case, t_end in [] if args.lag_only else plan:
            probe = np.logspace(-2, np.log10(t_end), 40)
            try:
                coarse = run(model, 46, case, t_end, dt_max=1.0, outputs=(1.0, 10.0))
                fine = run(model, 200, case, t_end, dt_max=1.0, outputs=(1.0, 10.0))
            except RuntimeError as error:
                report["grid"].append(dict(model=model.name, case=case, failed=str(error)))
                print("grid", report["grid"][-1], flush=True)
                args.output.write_text(json.dumps(report) + "\n")
                continue
            keep(coarse, role="grid")
            keep(fine, role="grid")
            row = dict(
                model=model.name,
                case=case,
                max_w_difference=compare(coarse, fine, probe),
                max_w_scale=float(np.array(coarse["history"])[:, 1].max()),
                final_coarse=coarse["final"],
                final_fine=fine["final"],
                newton=[coarse["newton_mean"], fine["newton_mean"]],
            )
            report["grid"].append(row)
            print("grid", row, flush=True)
            args.output.write_text(json.dumps(report) + "\n")
        # Lagged against implicit packing, and lagged at half the step, case B at 46 cells.
        t_end = args.t_end_b
        probe = np.logspace(-2, np.log10(t_end), 40)
        try:
            lagged = run(model, 46, "B", t_end, dt0=1e-3, dt_max=1.0)
            half = run(model, 46, "B", t_end, dt0=5e-4, dt_max=0.5)
            implicit = run(model, 46, "B", t_end, dt0=1e-3, dt_max=1.0, implicit=True)
        except RuntimeError as error:
            report["packing_lag"].append(dict(model=model.name, failed=str(error)))
            print("lag", report["packing_lag"][-1], flush=True)
            args.output.write_text(json.dumps(report) + "\n")
            continue
        for result, role in ((lagged, "lagged"), (half, "lagged_half_step"), (implicit, "implicit")):
            keep(result, role=role)
        row = dict(
            model=model.name,
            lagged_vs_implicit=compare(lagged, implicit, probe),
            lagged_vs_half_step=compare(lagged, half, probe),
            max_w_scale=float(np.array(lagged["history"])[:, 1].max()),
        )
        report["packing_lag"].append(row)
        print("lag", row, flush=True)
        args.output.write_text(json.dumps(report) + "\n")
    print(f"written: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
