"""Critical-state CIDR for the rigid-grain granular step.

The compressible I-dependent rheology of Barker et al. (2017), in the critical-state
form completed by the fixed cap p <= alpha(0) C (docs/granular-cidr-design.md):
yield Y = alpha(I) p - p^2 / C and dilatancy f = beta(I) - 2 p / C, with
mu(I) = mu_s + (mu_d - mu_s) I / (I_0 + I), I = k s / sqrt(p) and k = 2 d sqrt(rho_s).
The stress derives from the dissipation potential

    psi(s, delta) = sup_{0 <= p <= alpha(0) C} [g(s, p) - p delta],
    g = 2 p s [(5 alpha - 4 mu) / 3 - p / C],

of the shear rate s = |dev D| and the volumetric rate delta = tr D, which is
convex and has a cone at D = 0: the static yield, where material jams.

The time step minimises 1/2 |u - u*|^2_(M/dt) + sum_a W_a psi(D_a(u); C_a) over the
nodal velocities u, with D_a the strain rate at stress node a in the orthonormal
basis of granular.BASIS (component 0 volumetric, delta = sqrt(2) D_0; components 1
and 2 deviatoric, s = |D_dev| / sqrt(2)). It is solved by the method of multipliers
with semismooth Newton subproblems on the Moreau envelope of psi, whose gradient
comes from an exact proximal step that decides jamming node by node. The limit is
the rigid step itself: jammed nodes do not creep. Warm-started from the previous
step's stress, a step typically takes a few Newton iterations
(benchmarks/cidr_warm_steps.py).
"""

import numpy as np
from scipy import integrate, sparse
from scipy.interpolate import CubicSpline
from scipy.sparse.linalg import splu

ROOT2 = np.sqrt(2.0)


class CriticalStateCIDR:
    """psi and its derivatives for critical-state CIDR with the fixed cap; C is passed per point."""

    def __init__(self, mu_s, mu_d, i0, k):
        """mu(I) = mu_s + (mu_d - mu_s) I / (i0 + I) and I = k s / sqrt(p), k = 2 d sqrt(rho_s)."""
        self.mu_s, self.mu_d, self.i0, self.k = mu_s, mu_d, i0, k
        grid = np.logspace(-12, 3, 600)
        alpha = [
            0.8 * self.mu(I) + 0.48 * integrate.quad(lambda t: t ** (-0.6) * self.mu(t * I), 0, 1)[0]
            for I in grid
        ]
        self.spline = CubicSpline(np.log(grid), np.array(alpha))
        self.alpha0 = float(self.alpha(0.0))

    def mu(self, I):
        return self.mu_s + (self.mu_d - self.mu_s) * I / (self.i0 + I)

    def dmu(self, I):
        return (self.mu_d - self.mu_s) * self.i0 / (self.i0 + I) ** 2

    def alpha(self, I):
        return self.spline(np.log(np.clip(I, 1e-12, 1e3)))

    def functions(self, s, p, C):
        I = self.k * s / np.sqrt(p)
        a, m, dm = self.alpha(I), self.mu(I), self.dmu(I)
        da = np.where(I > 1e-10, (2 * m + 2 * I * dm - a) / (2.5 * np.maximum(I, 1e-300)), 0.0)
        b = 2 * (a - m)
        db = 2 * (da - dm)
        return I, a, da, b, db

    def p_cap(self, s, C):
        """The fixed cap alpha(0) C: the bound p <= alpha(I) C, where Y = 0, depends on s
        and makes psi nonconvex on the cap branch; the fixed bound keeps psi a supremum of
        functions convex in (s, delta), and Y >= 0 below it."""
        return self.alpha0 * C * np.ones_like(s)

    def evaluate_C(self, s, delta, C):
        """psi and its derivatives in (s, delta) at packing modulus C, and the pressure."""
        cap = self.p_cap(s, C)

        def gp(p):
            I, a, da, b, db = self.functions(s, p, C)
            return 2 * s * (b - 2 * p / C)

        # Bracketed bisection for g_p(s, p) = delta on (0, cap); g_p decreases in p.
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            lo, hi = np.full_like(s, 1e-300), cap.copy()
            at_cap = gp(cap) - delta >= 0
            p = 0.5 * (lo + hi)
            for _ in range(80):
                value = gp(p) - delta
                lo = np.where(value > 0, p, lo)
                hi = np.where(value <= 0, p, hi)
                p = 0.5 * (lo + hi)
            p = np.where(at_cap, cap, p)
            return self._interior_and_cap(s, delta, C, p, at_cap)

    def _interior_and_cap(self, s, delta, C, p, at_cap):
        I, a, da, b, db = self.functions(s, p, C)
        g = 2 * p * s * ((5 * a - 4 * self.mu(I)) / 3 - p / C)
        Y = a * p - p**2 / C
        g_s = 2 * Y
        g_ss = 2 * da * p * I / np.maximum(s, 1e-300)
        g_sp = 2 * (a - 2 * p / C) - da * I
        g_pp = 2 * s * (-2 / C - I * db / (2 * p))
        psi = g - p * delta
        interior = (psi, g_s, -p, g_ss - g_sp**2 / g_pp, g_sp / g_pp, -1 / g_pp)
        # Cap branch (p fixed at the cap): psi_s = g_s, psi_ss = g_ss, psi_delta = -p, rest zero.
        zero = np.zeros_like(p)
        out = tuple(np.where(at_cap, c, i) for i, c in zip(interior, (psi, g_s, -p, g_ss, zero, zero)))
        # Zero-pressure branch (dilating faster than g_p(s, 0+) allows): psi and its derivatives vanish.
        tiny = np.full_like(s, 1e-300)
        separating = ~at_cap & (2 * s * (self.functions(s, tiny, C)[3] - 2 * tiny / C) - delta <= 0)
        out = tuple(np.where(separating, 0.0, x) for x in out)
        return out, np.where(separating, 0.0, p)

    def static_support(self, s, delta, C):
        """sup_{0 <= p <= alpha_0 C} [2 s Y_0(p) - p delta]: the limit of psi(t s, t delta) / t."""
        a0 = self.alpha0
        s = np.asarray(s, dtype=float)
        p = np.where(
            s > 0, C * (2 * s * a0 - delta) / (4 * np.maximum(s, 1e-300)), np.where(delta < 0, np.inf, 0.0)
        )
        p = np.clip(p, 0.0, a0 * C)
        return 2 * s * (a0 * p - p**2 / C) - p * delta


def packing_modulus(phi, phi_min, phi_max, deficit, p_ref, mu_s):
    """C(phi) = (P/b) x, x = q / (1 - q), q = (phi - phi_min) / (phi_max - phi_min) (Fannon et al.).

    P/b is set by a packing deficit: material on the cap at rest at pressure p_ref sits
    `deficit` below phi_max, P/b = p_ref / (2 mu_s (Delta / deficit - 1)).
    """
    Delta = phi_max - phi_min
    scale = p_ref / (2 * mu_s * (Delta / deficit - 1))
    q = np.clip((phi - phi_min) / Delta, 1e-12, 1 - 1e-12)
    return scale * q / (1 - q)


class Prox:
    """Exact proximal step of psi / r at each stress node, in the orthonormal components."""

    def __init__(self, model, samples=256):
        self.model = model
        self.theta = np.linspace(-np.pi / 2, np.pi / 2, samples)

    def margin(self, f0, fv, C):
        """max over unit e = (e_dev, e_0), e_dev >= 0, of f.e - h(e), and the maximising angle.

        f is in the static set iff the margin is <= 0; otherwise e is the direction of
        steepest descent out of the cone."""

        def value(th, f0=f0, fv=fv, C=C):
            support = self.model.static_support(np.cos(th) / ROOT2, ROOT2 * np.sin(th), C)
            return fv * np.cos(th) + f0 * np.sin(th) - support

        k = np.argmax(value(self.theta[None, :], f0[:, None], fv[:, None], C[:, None]), axis=1)
        step = self.theta[1] - self.theta[0]
        lo = np.clip(self.theta[k] - step, -np.pi / 2, np.pi / 2)
        hi = np.clip(self.theta[k] + step, -np.pi / 2, np.pi / 2)
        for _ in range(40):
            m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
            lower = value(m1) < value(m2)
            lo, hi = np.where(lower, m1, lo), np.where(lower, hi, m2)
        theta = 0.5 * (lo + hi)
        return value(theta), theta

    def terms(self, sigma, d0, nv, q0, r, C):
        """F = psi(sigma / sqrt 2, sqrt 2 d0) + r/2 [(sigma - nv)^2 + (d0 - q0)^2], gradient, Hessian."""
        (psi, ps, pd, pss, psd, pdd), p = self.model.evaluate_C(sigma / ROOT2, ROOT2 * d0, C)
        F = psi + 0.5 * r * ((sigma - nv) ** 2 + (d0 - q0) ** 2)
        g = np.stack((ps / ROOT2 + r * (sigma - nv), ROOT2 * pd + r * (d0 - q0)), axis=-1)
        H = np.empty(sigma.shape + (2, 2))
        H[..., 0, 0] = 0.5 * pss + r
        H[..., 0, 1] = H[..., 1, 0] = psd
        H[..., 1, 1] = 2 * pdd + r
        return F, g, H, p

    def __call__(self, q, r, C):
        """argmin_d psi(d) + r/2 |d - q|^2 for each row of q; returns d, jammed and pressure."""
        q0, qv = q[:, 0], q[:, 1:]
        nv = np.linalg.norm(qv, axis=1)
        scale = r * np.linalg.norm(q, axis=1)
        margin, theta = self.margin(r * q0, r * nv, C)
        jam = margin <= 1e-12 * np.maximum(scale, 1e-300)
        d0 = np.zeros(len(q))
        sigma = np.zeros(len(q))
        pressure = np.zeros(len(q))
        P = self.model.alpha0 * C
        # Pure volumetric (no deviatoric part): soft threshold, psi(0, delta) = P max(-delta, 0).
        flat = ~jam & (nv <= 1e-14 * np.maximum(np.abs(q0), 1e-300))
        d0[flat] = np.where(q0[flat] >= 0, q0[flat], np.minimum(q0[flat] + ROOT2 * P[flat] / r, 0.0))
        pressure[flat] = np.where(d0[flat] < 0, P[flat], 0.0)
        # Otherwise the minimiser has sigma > 0: damped Newton from (|q_dev|, q_0) or from
        # just outside the cone along the steepest descent direction, whichever is lower
        # (near D = 0 the Hessian grows like 1 / |D| and Newton does not leave there).
        rest = np.flatnonzero(~jam & ~flat)
        sigma[rest], d0[rest] = nv[rest], q0[rest]
        if len(rest):
            exit_ = margin[rest] / r
            s_exit, d_exit = exit_ * np.cos(theta[rest]), exit_ * np.sin(theta[rest])
            usable = s_exit > 0
            F_q = self.terms(nv[rest], q0[rest], nv[rest], q0[rest], r, C[rest])[0]
            F_exit = np.full(len(rest), np.inf)
            F_exit[usable] = self.terms(
                s_exit[usable], d_exit[usable], nv[rest][usable], q0[rest][usable], r, C[rest][usable]
            )[0]
            better = F_exit < F_q
            sigma[rest[better]], d0[rest[better]] = s_exit[better], d_exit[better]
        live = rest.copy()
        for _ in range(60):
            if len(live) == 0:
                break
            F, g, H, p = self.terms(sigma[live], d0[live], nv[live], q0[live], r, C[live])
            pressure[live] = p
            det = H[:, 0, 0] * H[:, 1, 1] - H[:, 0, 1] ** 2
            step = (
                -np.stack(
                    (
                        H[:, 1, 1] * g[:, 0] - H[:, 0, 1] * g[:, 1],
                        H[:, 0, 0] * g[:, 1] - H[:, 0, 1] * g[:, 0],
                    ),
                    axis=-1,
                )
                / det[:, None]
            )
            slope = np.sum(g * step, axis=1)
            done = -slope <= 1e-28 + 1e-15 * (np.abs(F) + r * (nv[live] ** 2 + q0[live] ** 2))
            t = np.ones(len(live))
            todo = ~done
            for _ in range(60):
                if not np.any(todo):
                    break
                idx = np.flatnonzero(todo)
                s_t = sigma[live[idx]] + t[idx] * step[idx, 0]
                d_t = d0[live[idx]] + t[idx] * step[idx, 1]
                ok = s_t > 0
                F_t = np.full(len(idx), np.inf)
                if np.any(ok):
                    nodes = live[idx[ok]]
                    F_t[ok] = self.terms(s_t[ok], d_t[ok], nv[nodes], q0[nodes], r, C[nodes])[0]
                accept = F_t <= F[idx] + 1e-4 * t[idx] * slope[idx]
                todo[idx[accept]] = False
                t[idx[~accept]] *= 0.5
            move = ~done
            sigma[live[move]] += t[move] * step[move, 0]
            d0[live[move]] += t[move] * step[move, 1]
            live = live[move]
        else:
            raise RuntimeError(f"local Newton did not converge at {len(live)} nodes")
        d = np.zeros_like(q)
        d[:, 0] = d0
        direction = np.divide(qv, nv[:, None], out=np.zeros_like(qv), where=nv[:, None] > 0)
        d[:, 1:] = sigma[:, None] * direction
        return d, jam, pressure


class Problem:
    """One step: D(u) = A u + a0 per stress node, weights W, masses / dt, trial velocity, C per node.

    u and trial are the free velocity components; A maps them to strain rates (3 per node).
    """

    def __init__(self, A, a0, W, mass, trial, model, C):
        self.A, self.a0, self.W = A, a0, W
        self.mass, self.trial = mass, trial
        self.model, self.C = model, np.asarray(C, dtype=float)
        self.n = len(self.C)
        self._prox = Prox(model)

    def prox(self, x, r):
        return self._prox(x, r, self.C)

    def D(self, u):
        return (self.A @ u + self.a0).reshape(self.n, 3)

    def hessian(self, d, C):
        """Hessian of psi in the orthonormal components at strain rates d with a deviatoric part."""
        nv = np.linalg.norm(d[:, 1:], axis=1)
        s, delta = nv / ROOT2, ROOT2 * d[:, 0]
        (_, ps, _, pss, psd, pdd), _ = self.model.evaluate_C(s, delta, C)
        n = d[:, 1:] / nv[:, None]
        H = np.zeros((len(d), 3, 3))
        H[:, 0, 0] = 2 * pdd
        H[:, 0, 1:] = H[:, 1:, 0] = psd[:, None] * n
        outer = n[:, :, None] * n[:, None, :]
        H[:, 1:, 1:] = 0.5 * pss[:, None, None] * outer + (ps / (ROOT2 * nv))[:, None, None] * (
            np.eye(2) - outer
        )
        return H

    def envelope(self, x, gamma, derivatives=True):
        """Moreau envelope psi_gamma(x) = min_d psi(d) + |d - x|^2 / (2 gamma) at each node.

        Its gradient is (x - d) / gamma with d = prox_{gamma psi}(x), and its generalised
        Hessian (I - dd/dx) / gamma: r I at jammed nodes (r = 1 / gamma), and
        r (I - r (H + r I)^-1) where d has a deviatoric part, H the Hessian of psi at d.
        """
        r = 1 / gamma
        d, jam, _ = self.prox(x, r)
        nv = np.linalg.norm(d[:, 1:], axis=1)
        smooth = ~jam & (nv > 0)
        flat = ~jam & (nv == 0)
        psi = np.zeros(self.n)
        if np.any(smooth):
            (value, *_), _ = self.model.evaluate_C(nv[smooth] / ROOT2, ROOT2 * d[smooth, 0], self.C[smooth])
            psi[smooth] = value
        psi[flat] = ROOT2 * self.model.alpha0 * self.C[flat] * np.maximum(-d[flat, 0], 0.0)
        value = psi + 0.5 * r * np.sum((d - x) ** 2, axis=1)
        if not derivatives:
            return value
        gradient = r * (x - d)
        H = np.zeros((self.n, 3, 3))
        H[jam] = r * np.eye(3)
        if np.any(smooth):
            Hpsi = self.hessian(d[smooth], self.C[smooth])
            H[smooth] = r * (np.eye(3) - r * np.linalg.inv(Hpsi + r * np.eye(3)))
        # Pure volumetric nodes: the deviatoric part is held at zero (stiff), the volumetric
        # part passes through unless thresholded at the cap.
        H[flat] = r * np.eye(3)
        passing = flat & (d[:, 0] == x[:, 0])
        H[passing, 0, 0] = 0.0
        return value, gradient, H, jam

    def objective(self, u, gamma, shift):
        e = u - self.trial
        x = self.D(u) + gamma * shift
        return 0.5 * e @ (self.mass * e) + np.sum(self.W[::3] * self.envelope(x, gamma, derivatives=False))


def newton(problem, u, gamma, shift, *, max_iterations=200):
    """Semismooth Newton with a backtracking line search on f(u) + sum W psi_gamma(D(u) + gamma shift).

    Returns u, the jammed nodes, grad psi_gamma (the updated multiplier) and the iterations.
    """
    A, W, mass = problem.A, problem.W, problem.mass
    scale = max(float(np.sum(mass * problem.trial**2)), 1e-300)
    for iteration in range(1, max_iterations + 1):
        value, gradient, H, jam = problem.envelope(problem.D(u) + gamma * shift, gamma)
        G = mass * (u - problem.trial) + A.T @ (W * gradient.ravel())
        n = len(H)
        blocks = sparse.bsr_matrix(
            (W[::3, None, None] * H, np.arange(n), np.arange(n + 1)), shape=(3 * n, 3 * n)
        )
        K = (sparse.diags(mass) + A.T @ blocks @ A).tocsc()
        delta = -splu(K).solve(G)
        slope = G @ delta
        # A decrement of 1e-11 of the kinetic scale is a velocity error of about 3e-6; below
        # 1e-9 after many iterations, Newton is held up by round-off.
        if -slope <= 1e-11 * scale or (iteration > 20 and -slope <= 1e-9 * scale):
            return u, jam, gradient, iteration - 1
        f = 0.5 * (u - problem.trial) @ (mass * (u - problem.trial)) + np.sum(W[::3] * value)
        t = 1.0
        while t > 1e-12:
            trial = u + t * delta
            if problem.objective(trial, gamma, shift) <= f + 1e-4 * t * slope:
                break
            t *= 0.5
        else:
            # No decrease representable in floating point: accept if already near the minimum.
            if -slope <= 1e-9 * scale:
                return u, jam, gradient, iteration - 1
        u = trial
    raise RuntimeError(f"Newton did not converge at gamma {gamma:.2e}: decrement {-slope:.2e}")


def solve(problem, u, lam, gamma, *, cold=None, reduction=np.sqrt(10), tol=1e-5, max_outer=200):
    """Method of multipliers from (u, lam): returns u, lam (the stress), jammed nodes and statistics.

    Each outer step minimises f(u) + sum W psi_gamma(D(u) + gamma lam) by semismooth Newton,
    then sets lam = grad psi_gamma(D(u) + gamma lam); the limit is the step with psi itself
    at any gamma. A cold start (cold = the starting gamma, larger than gamma) reduces gamma
    by `reduction` per outer step down to `gamma`; a warm start stays at `gamma`. Stops at
    `gamma` once the velocity changes by less than tol (relative to the larger of the
    velocity and the trial velocity, in the mass norm).
    """
    norm = lambda v: np.sqrt(v @ (problem.mass * v))  # noqa: E731
    current = gamma if cold is None else max(cold, gamma)
    newton_total = 0
    for outer in range(1, max_outer + 1):
        previous = u
        u, jam, lam, iterations = newton(problem, u, current, lam)
        newton_total += iterations
        # Relative to the larger of the velocity and the trial velocity: a jammed body has u = 0.
        if current <= gamma and norm(u - previous) <= tol * max(norm(u), norm(problem.trial), 1e-300):
            break
        current = max(current / reduction, gamma)
    else:
        raise RuntimeError(f"CIDR step did not settle in {max_outer} outer steps")
    return u, lam, jam, dict(outer=outer, newton=newton_total, jammed=int(jam.sum()))
