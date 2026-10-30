"""One-dimensional gravity-free shear cell with the convex CIDR time step.

The test of Schaeffer et al. (2019, sections 3 and 5): granular material between
a fixed bottom wall and a top wall moving at unit speed, no slip, invariant in
the shearing direction. Variables are non-dimensional (grain diameter, grain
density and wall speed are 1). Each cell keeps its mass and the nodes move with
the vertical velocity (Lagrangian), so there is no advection term. Each time step
solves the convex problem of docs/granular-cidr-design.md,

    min_v  sum_i M_i / (2 dt) |v_i - v_i^n|^2 + sum_c dz_c psi(s_c, delta_c; phi_c),

with v = (u, w) at the nodes, s = |D| = sqrt(u'^2 + w'^2) / 2 and delta = w'.
phi is held at its start-of-step value ("lagged") or iterated to its end-of-step
value ("implicit").

Models: inertial CIDR (Schaeffer et al. 2019), solved by Newton's method on the
nodal velocities, and critical-state CIDR (Barker et al. 2017) completed by the
fixed cap p <= alpha(0) C, with C(phi) from a tunable packing deficit. The latter
has a static yield, so psi has a cone at zero strain rate, and its step is
solved in cell strain rates by a truncated nonsmooth Newton method
(Cell.step_nonsmooth).

    python -m benchmarks.cidr_shear_cell --output benchmarks/cidr_shear_cell.json
"""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import solve_banded

from cidr import CriticalStateCIDR, packing_modulus

# Schaeffer et al. (2019), table 1: DEM fit of mu(I) and Phi(I).
MU_S, MU_D, I_0, PHI_C, A_PHI = 0.248, 0.588, 0.145, 0.817, 0.314
K = 2.0  # I = K s / sqrt(p) with d = rho_s = 1


def mu(I):
    return MU_S + (MU_D - MU_S) * I / (I_0 + I)


def dmu(I):
    return (MU_D - MU_S) * I_0 / (I_0 + I) ** 2


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


class CriticalState(CriticalStateCIDR):
    """Critical-state CIDR (cidr.py) with C = (P/b) x, x = psi / (1 - psi), from the packing.

    P/b follows from a packing deficit e_ref on the cap at pressure p_ref:
    P/b = p_ref / (2 mu(0) (Delta / e_ref - 1)), Delta = phi_max - phi_min.
    """

    name = "critical_state"

    def __init__(
        self, deficit=0.005, p_ref=0.08, phi_min=0.05, phi_max=PHI_C, mu_s=MU_S, mu_d=MU_D, i0=I_0, k=K
    ):
        """mu(I) = mu_s + (mu_d - mu_s) I / (i0 + I) and I = k s / sqrt(p), k = 2 d sqrt(rho_s)."""
        super().__init__(mu_s, mu_d, i0, k)
        self.phi_min, self.phi_max = phi_min, phi_max
        self.deficit, self.p_ref = deficit, p_ref

    def C(self, phi):
        return packing_modulus(phi, self.phi_min, self.phi_max, self.deficit, self.p_ref, self.mu_s)

    def evaluate(self, s, delta, phi):
        return self.evaluate_C(s, delta, self.C(phi))

    def support(self, ea, eb, phi):
        """Support function h(e) of the static set (the subdifferential of psi at D = 0), |e| = 1.

        Near D = 0, psi(t e) / t -> h(e) = sup_{0 <= p <= alpha_0 C} [Y_0(p) - p e_b] with
        Y_0 = alpha_0 p - p^2 / C, so the static set is the union over p of discs of radius
        Y_0(p) centred at (0, -p) in (psi_a, psi_b).
        """
        return self.static_support(0.5, eb, self.C(phi))

    def static_pressure(self, ta, tb, phi):
        """Pressure for a stress (ta, tb) in the static set: the disc it lies deepest inside."""
        C = self.C(phi)
        a0 = self.alpha0
        lo, hi = np.zeros_like(ta), a0 * C
        excess = lambda p: np.hypot(ta, tb + p) - (a0 * p - p**2 / C)  # noqa: E731 (convex in p)
        for _ in range(60):
            m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
            smaller = excess(m1) < excess(m2)
            hi, lo = np.where(smaller, m2, hi), np.where(smaller, lo, m1)
        return 0.5 * (lo + hi)


def static_margin(model, fa, fb, phi, samples=512):
    """max over unit e of f.e - h(e) and its maximiser: f is in the static set iff the margin <= 0."""
    theta = np.linspace(0.0, 2 * np.pi, samples, endpoint=False)
    value = lambda th: fa * np.cos(th) + fb * np.sin(th) - model.support(np.cos(th), np.sin(th), phi)
    k = int(np.argmax(value(theta)))
    lo, hi = theta[k] - 2 * np.pi / samples, theta[k] + 2 * np.pi / samples
    for _ in range(40):
        m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
        if value(m1) < value(m2):
            lo = m1
        else:
            hi = m2
    th = 0.5 * (lo + hi)
    return float(value(th)), np.array([np.cos(th), np.sin(th)])


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


class StepRejected(Exception):
    """An implicit packing pass left the admissible range of phi."""


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

    def step_nonsmooth(self, dt, phi):
        """The convex step for a model with a static yield, in cell strain rates.

        psi has a cone at zero strain rate (a jammed cell). With the strain rates
        d_c = (u', w') of the cells as unknowns that kink is separable, and the walls
        become the linear constraint sum_c dz_c d_c = (1, 0). Each iteration is one
        exact block Gauss-Seidel pass over slow and jammed cells, which jams or
        releases them, then a Newton step on the free cells with the jammed ones held
        at zero (a truncated nonsmooth Newton method), with a line search.
        """
        model = self.model
        dz = np.diff(self.z)
        N = len(dz)
        nodal_mass = np.zeros(N + 1)
        nodal_mass[:-1] += 0.5 * self.mass
        nodal_mass[1:] += 0.5 * self.mass
        md = nodal_mass[1:N] / dt  # interior nodes 1 .. N-1
        V0 = np.column_stack((self.u, self.w))
        d = np.diff(V0, axis=0) / dz[:, None]
        jam = np.hypot(d[:, 0], d[:, 1]) < 1e-13
        d[jam] = 0.0
        scale = np.sum(md * np.sum(V0[1:N] ** 2, axis=1))
        # Kinetic term K(d) = sum_i md_i |V_i(d) - V0_i|^2 / 2 with V_i = sum_{c<i} dz_c d_c:
        # its Hessian is Q_cd = dz_c dz_d sum_{i > max(c, d)} md_i, the same for u and w.
        tail = np.concatenate((np.cumsum(md[::-1])[::-1], [0.0]))  # tail[k] = sum_{i >= k+1} md_i
        Q = np.outer(dz, dz) * tail[np.maximum.outer(np.arange(N), np.arange(N))]

        def kinetic(d):
            V = np.cumsum(dz[:, None] * d, axis=0)[:-1]  # nodes 1 .. N-1
            R = md[:, None] * (V - V0[1:N])
            T = np.concatenate((np.cumsum(R[::-1], axis=0)[::-1], np.zeros((1, 2))))
            return 0.5 * np.sum(R * (V - V0[1:N])), dz[:, None] * T, T

        def terms(cells, x):
            psi, grad, hess, p = cell_terms(model, x[:, 0], x[:, 1], phi[cells])
            return psi * dz[cells], grad * dz[cells, None], hess * dz[cells, None, None], p

        def objective(d, jam):
            free = np.flatnonzero(~jam)
            return kinetic(d)[0] + np.sum(terms(free, d[free])[0])

        def sweep(d, jam):
            """Exact block Gauss-Seidel on slow and jammed cells, each paired with the fastest cell."""
            r = np.hypot(d[:, 0], d[:, 1])
            star = int(np.argmax(r))
            changed = False
            for c in np.flatnonzero(jam | (r < 1e-2 * np.mean(r))):
                if c == star:
                    continue
                kappa = dz[c] / dz[star]
                q = Q[c, c] - 2 * kappa * Q[c, star] + kappa**2 * Q[star, star]
                gK = kinetic(d)[1]
                linear = gK[c] - kappa * gK[star]

                # Smooth part of the local problem in Delta (d_c += Delta, d_star -= kappa Delta).
                def local(delta):
                    x = np.array([d[c] + delta, d[star] - kappa * delta])
                    psi, grad, hess, _ = terms(np.array([c, star]), x)
                    value = linear @ delta + 0.5 * q * delta @ delta + psi.sum()
                    gradient = linear + q * delta + grad[0] - kappa * grad[1]
                    return value, gradient, q * np.eye(2) + hess[0] + kappa**2 * hess[1]

                delta0 = -d[c]
                g_star = terms(np.array([star]), (d[star] + kappa * d[c])[None])[1][0]
                G = linear + q * delta0 - kappa * g_star
                margin, e = static_margin(model, -G[0] / dz[c], -G[1] / dz[c], phi[c])
                if margin <= 1e-12 * max(1.0, np.hypot(*G) / dz[c]):
                    if not jam[c]:
                        d[star] += kappa * d[c]
                        d[c] = 0.0
                        jam[c] = changed = True
                    continue
                # Not jammed: minimise the local problem, starting from the current point or
                # from just outside the cone along the steepest descent direction, whichever
                # is lower (Newton alone cannot leave the neighbourhood of D = 0, where the
                # Hessian grows like 1 / |D|).
                eps = dz[c] * margin / q
                base = (
                    linear @ delta0
                    + 0.5 * q * delta0 @ delta0
                    + terms(np.array([star]), (d[star] + kappa * d[c])[None])[0][0]
                )
                while True:
                    delta = delta0 + eps * e
                    f, g, h = local(delta)
                    if f < base or eps < 1e-300:
                        break
                    eps *= 0.5
                if not jam[c]:
                    current = local(np.zeros(2))
                    if current[0] < f:
                        delta, (f, g, h) = np.zeros(2), current
                for _ in range(20):
                    step = -np.linalg.solve(h, g)
                    slope = g @ step
                    if -slope <= 1e-15 * max(abs(f), 1e-300):
                        break
                    t = 1.0
                    while t > 1e-12:
                        ft, gt, ht = local(delta + t * step)
                        if ft <= f + 1e-4 * t * slope:
                            break
                        t *= 0.5
                    delta, f, g, h = delta + t * step, ft, gt, ht
                d[c] += delta
                d[star] -= kappa * delta
                changed |= bool(jam[c])
                jam[c] = False
            return changed

        lam = np.zeros(2)
        for iteration in range(100):
            changed = sweep(d, jam)
            free = np.flatnonzero(~jam)
            F = len(free)
            _, gK, _ = kinetic(d)
            _, grad, hess, _ = terms(free, d[free])
            g = (gK[free] + grad).T.ravel()  # [a-components, b-components]
            H = np.zeros((2 * F + 2, 2 * F + 2))
            QF = Q[np.ix_(free, free)]
            H[:F, :F] = QF + np.diag(hess[:, 0, 0])
            H[F : 2 * F, F : 2 * F] = QF + np.diag(hess[:, 1, 1])
            H[:F, F : 2 * F] = H[F : 2 * F, :F] = np.diag(hess[:, 0, 1])
            H[2 * F, :F] = H[:F, 2 * F] = -dz[free]
            H[2 * F + 1, F : 2 * F] = H[F : 2 * F, 2 * F + 1] = -dz[free]
            solution = np.linalg.solve(H, np.concatenate((-g, [0.0, 0.0])))
            direction = np.zeros_like(d)
            direction[free, 0], direction[free, 1] = solution[:F], solution[F : 2 * F]
            lam = solution[2 * F :]
            slope = g @ solution[: 2 * F]
            if not changed and -slope <= 1e-24 + 1e-15 * scale:
                break
            f = objective(d, jam)
            t = 1.0
            while t > 1e-12:
                trial = d + t * direction
                if objective(trial, jam) <= f + 1e-4 * t * slope:
                    break
                t *= 0.5
            d = trial
        else:
            raise RuntimeError(
                f"nonsmooth Newton did not converge: decrement {-slope:.3g}, scale {scale:.3g}"
            )
        V = np.vstack((np.zeros((1, 2)), np.cumsum(dz[:, None] * d, axis=0)))
        V[-1] = V0[-1]
        self.u, self.w = V[:, 0].copy(), V[:, 1].copy()
        # Stress per cell: psi' = lam - T_{c+1}; set-valued in jammed cells, where the
        # pressure reported is that of the static disc the stress lies deepest inside.
        p = np.zeros(N)
        free = np.flatnonzero(~jam)
        p[free] = terms(free, d[free])[3]
        stress = lam - kinetic(d)[2][:N]
        cells = np.flatnonzero(jam)
        if len(cells):
            p[cells] = model.static_pressure(stress[cells, 0], stress[cells, 1], phi[cells])
        self.jammed = int(jam.sum())
        return iteration, p

    def advance(self, dt, implicit=False, limit=None):
        z0, u0, w0 = self.z.copy(), self.u.copy(), self.w.copy()
        phi = self.phi
        for _ in range(30 if implicit else 1):
            # Each pass restarts from the start-of-step state; only phi changes.
            self.z, self.u, self.w = z0, u0.copy(), w0.copy()
            if hasattr(self.model, "support"):
                iterations, p = self.step_nonsmooth(dt, phi)
            else:
                iterations, p = self.step(dt, phi)
            z_new = z0 + dt * self.w
            phi_new = self.mass / np.diff(z_new)
            if implicit and limit is not None and np.any(phi_new >= limit):
                raise StepRejected
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
        if limit is None:
            iterations, p = cell.advance(step, implicit)
        else:
            # The estimate above uses the start-of-step rate; check the solved step and
            # halve it if any cell closed more than half its gap, or an implicit pass
            # left the admissible range of phi.
            saved = (cell.z.copy(), cell.u.copy(), cell.w.copy())
            while True:
                phi = cell.phi
                try:
                    iterations, p = cell.advance(step, implicit, limit)
                    if np.all(cell.phi - phi <= 0.5 * (limit - phi)):
                        break
                except StepRejected:
                    pass
                cell.z, cell.u, cell.w = (x.copy() for x in saved)
                step *= 0.5
                dt = step
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
