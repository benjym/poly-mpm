"""Variational structure of the compressible I-dependent rheology (CIDR).

Checks, symbolically where possible, the derivations recorded in
docs/granular-cidr-design.md:

A. For yield and dilatancy functions Y(p, I), f(p, I) with I = k s / sqrt(p)
   (s = |D|, k = 2 d sqrt(rho_s)), a function g(s, p) with g_s = 2 Y and
   g_p = 2 f s exists exactly when Barker et al.'s (2017) condition
   Y_p - (I / 2p) Y_I = f + I f_I holds; then g_ss = 2 Y_I I / s and
   g_pp = 2 s (f_p - (I / 2p) f_I), so their inequalities are convexity in s
   and concavity in p.
B. The construction reproduces the inertial CIDR potential of Schaeffer et al.
   (2019, eq. B8).
C. For critical-state CIDR, g = 2 p s [(5 alpha(I) - 4 mu(I)) / 3 - p / C].
D. Its static limit: g_s -> 2 p (2 mu(0) - p / C) as s -> 0.
E. The bound p <= alpha(I) C, where Y >= 0: its shape in (s, p).
F. Convexity of the dissipation potential psi with the fixed cap
   p <= alpha(0) C (used), the shear-dependent cap p <= alpha(I) C, and none.
G. Steady flow, critical state and the cap at rest.
H. The packing scale of C and the effective compressibility it implies.
I. psi's gradient and Hessian on each branch against finite differences.

    python -m benchmarks.cidr_potential            (sympy and scipy required)
"""

import numpy as np
import sympy as sp
from scipy import integrate, optimize


def general():
    s, p, k = sp.symbols("s p k", positive=True)
    J = sp.Symbol("I", positive=True)
    Y, f = sp.Function("Y"), sp.Function("f")
    I = k * s / sp.sqrt(p)
    gs, gp = 2 * Y(p, I), 2 * f(p, I) * s
    Yp, YI = sp.Derivative(Y(p, J), p), sp.Derivative(Y(p, J), J)
    fp, fI = sp.Derivative(f(p, J), p), sp.Derivative(f(p, J), J)

    def at(expr):
        return sp.simplify(expr.doit().subs(k * s / sp.sqrt(p), J))

    mixed = at(sp.diff(gs, p) - sp.diff(gp, s))
    condition = 2 * ((Yp - J / (2 * p) * YI) - (f(p, J) + J * fI))
    convex = at(sp.diff(gs, s)) - 2 * YI * J / s
    concave = at(sp.diff(gp, p)) - 2 * s * (fp - J / (2 * p) * fI)
    return dict(
        integrability_is_condition_a=sp.simplify(at(mixed.subs(J, J)) - at(condition)) == 0,
        g_ss_is_2_Y_I_I_over_s=sp.simplify(convex.doit().subs(k, J * sp.sqrt(p) / s)) == 0,
        g_pp_is_2s_condition_c=sp.simplify(concave.doit().subs(k, J * sp.sqrt(p) / s)) == 0,
    )


def inertial():
    s, p, d, rho, M, Psi = sp.symbols("s p d rho M Psi", positive=True)
    k = 2 * d * sp.sqrt(rho)
    I = k * s / sp.sqrt(p)
    Y = M * I / Psi * p
    f = M / 4 * (I / Psi - Psi / I)
    g = M * k * s**2 * sp.sqrt(p) / Psi - M * Psi * p ** sp.Rational(3, 2) / (3 * k)
    paper = 2 * d * sp.sqrt(rho) * M / Psi * s**2 * sp.sqrt(p) - M * Psi / (
        6 * d * sp.sqrt(rho)
    ) * p ** sp.Rational(3, 2)
    return dict(
        g_s_is_2Y=sp.simplify(sp.diff(g, s) - 2 * Y) == 0,
        g_p_is_2fs=sp.simplify(sp.diff(g, p) - 2 * f * s) == 0,
        matches_schaeffer_B8=sp.simplify(g - paper) == 0,
    )


def critical_state():
    s, p, C = sp.symbols("s p C", positive=True)
    J = sp.Symbol("I", positive=True)
    mu, alpha = sp.Function("mu"), sp.Function("alpha")
    # Treat I as its own variable; with I = k s / sqrt(p) the chain rule gives
    # d/ds = d_s + (I / s) d_I and d/dp = d_p - (I / 2p) d_I.
    dds = lambda e: sp.diff(e, s) + J / s * sp.diff(e, J)  # noqa: E731
    ddp = lambda e: sp.diff(e, p) - J / (2 * p) * sp.diff(e, J)  # noqa: E731
    # Barker et al. (2017): alpha - beta / 2 = mu and the ODE (27) from condition (a).
    beta = 2 * (alpha(J) - mu(J))
    dalpha = (2 * mu(J) + 2 * J * sp.diff(mu(J), J) - alpha(J)) / (sp.Rational(5, 2) * J)
    Y = alpha(J) * p - p**2 / C
    f = beta - 2 * p / C
    g = 2 * p * s * ((5 * alpha(J) - 4 * mu(J)) / 3 - p / C)

    def reduce(expr):
        return sp.simplify(expr.subs(sp.Derivative(alpha(J), J), dalpha))

    return dict(
        g_s_is_2Y=reduce(dds(g) - 2 * Y) == 0,
        g_p_is_2fs=reduce(ddp(g) - 2 * f * s) == 0,
        static_slope=sp.simplify(dds(g).subs(J, 0)),
    )


class Jop:
    """mu(I) of Jop et al. with the glass-bead values used by Fannon et al. (2019)."""

    mu1, mu2, I0 = np.tan(np.radians(20.9)), np.tan(np.radians(32.76)), 0.279

    def mu(self, I):
        return self.mu1 + (self.mu2 - self.mu1) * I / (self.I0 + I)

    def alpha(self, I):
        # Barker et al. (2017), eq. (25), reparametrised: alpha = 4 mu / 5 + 12/25 int_0^1 s^-3/5 mu(s I) ds.
        integral = integrate.quad(lambda t: t ** (-0.6) * self.mu(t * I), 0, 1)[0]
        return 0.8 * self.mu(I) + 0.48 * integral

    def beta(self, I):
        return 2 * (self.alpha(I) - self.mu(I))


def numeric(seed=0):
    """Finite-difference checks of C and D with the Jop mu(I), and conditions (b), (c)."""
    m = Jop()
    rng = np.random.default_rng(seed)
    k, C = 2 * 0.003 * np.sqrt(2500.0), 4e4  # d = 3 mm, rho_s = 2500 kg/m^3, an arbitrary C
    I_of = lambda s, p: k * s / np.sqrt(p)  # noqa: E731
    Y = lambda s, p: m.alpha(I_of(s, p)) * p - p**2 / C  # noqa: E731
    f = lambda s, p: m.beta(I_of(s, p)) - 2 * p / C  # noqa: E731
    g = lambda s, p: 2 * p * s * ((5 * m.alpha(I_of(s, p)) - 4 * m.mu(I_of(s, p))) / 3 - p / C)  # noqa: E731
    worst_s = worst_p = 0.0
    for _ in range(40):
        s, p = 10 ** rng.uniform(-2, 1.5), 10 ** rng.uniform(1, 4)
        hs, hp = 1e-6 * s, 1e-6 * p
        gs = (g(s + hs, p) - g(s - hs, p)) / (2 * hs)
        gp = (g(s, p + hp) - g(s, p - hp)) / (2 * hp)
        worst_s = max(worst_s, abs(gs - 2 * Y(s, p)) / abs(2 * Y(s, p)))
        worst_p = max(worst_p, abs(gp - 2 * f(s, p) * s) / max(abs(2 * f(s, p) * s), 1e-12))
    alpha0 = m.alpha(1e-12)
    # Conditions (b), (c) on a grid of I (with C, p fixed they reduce to signs of
    # alpha' and of f_p - I/(2p) f_I = -2/C - (I/2p) beta'(I)).
    Is = np.logspace(-6, 1, 60)
    alpha = np.array([m.alpha(x) for x in Is])
    beta = np.array([m.beta(x) for x in Is])
    dalpha = np.gradient(alpha, Is)
    dbeta = np.gradient(beta, Is)
    return dict(
        g_s_relative_error=worst_s,
        g_p_relative_error=worst_p,
        alpha_0_over_mu1=alpha0 / m.mu1,
        min_alpha_prime=float(dalpha.min()),
        min_beta_prime=float(dbeta.min()),
    )


def cap(C=4e4, k=2 * 0.003 * np.sqrt(2500.0)):
    """The bound p_cap(s) solving p = alpha(k s / sqrt(p)) C, where Y = 0."""
    m = Jop()
    ss = np.linspace(0.0, 20.0, 81)
    caps = []
    for s in ss:
        h = lambda p: p - m.alpha(k * s / np.sqrt(p)) * C  # noqa: E731
        caps.append(optimize.brentq(h, 1e-6, 10 * m.mu2 * C))
    caps = np.array(caps)
    second = np.diff(caps, 2)
    return dict(
        p_cap_at_rest_over_C=caps[0] / C,
        p_cap_increasing=bool(np.all(np.diff(caps) >= -1e-9 * C)),
        p_cap_concave=bool(np.all(second <= 1e-9 * C)),
        max_second_difference_over_C=float(second.max() / C),
    )


def dissipation(cap, C=4e4, k=2 * 0.003 * np.sqrt(2500.0), seed=1, trials=400):
    """psi(s, delta) = sup_p [g(s, p) - p delta] over 0 <= p <= bound.

    cap is "fixed" (bound alpha(0) C), "shear" (bound p_cap(s), where Y = 0) or None.
    Tests midpoint convexity in (s, delta) and monotonicity in s on random pairs, half
    of them far apart and half close together (the shear-dependent cap fails only
    locally); psi is a function of the tensor D through s = |D| and is convex in D if
    it is convex and nondecreasing in s. Returns the worst violations, relative to |psi|.
    """
    m = Jop()
    rng = np.random.default_rng(seed)
    grid = np.logspace(-8, 4, 300)
    table = np.array([m.alpha(x) for x in grid])
    alpha = lambda I: np.interp(np.log(np.maximum(I, grid[0])), np.log(grid), table)  # noqa: E731
    mu = m.mu

    def g(s, p):
        I = k * s / np.sqrt(p) if p > 0 else np.inf
        if not np.isfinite(I):
            return 0.0
        return 2 * p * s * ((5 * alpha(I) - 4 * mu(I)) / 3 - p / C)

    def p_cap(s):
        return optimize.brentq(lambda p: p - alpha(k * s / np.sqrt(p)) * C, 1e-9, 10 * m.mu2 * C)

    def psi(s, delta):
        if cap == "fixed" or (cap == "shear" and s == 0):
            upper = alpha(0.0) * C
        elif cap == "shear":
            upper = p_cap(s)
        else:
            upper = 50 * m.mu2 * C
        result = optimize.minimize_scalar(
            lambda p: -(g(s, p) - p * delta),
            bounds=(0.0, upper),
            method="bounded",
            options=dict(xatol=1e-10 * upper),
        )
        return -result.fun, result.x, upper

    convexity = monotone = 0.0
    at_cap = 0
    for _ in range(trials):
        a = np.array([rng.uniform(0, 10), rng.uniform(-5, 5)])
        b = np.array([rng.uniform(0, 10), rng.uniform(-5, 5)])
        if _ % 2:
            b = a + 10 ** rng.uniform(-3, -1) * rng.standard_normal(2) * max(a[0], 0.1)
            b[0] = abs(b[0])
        values = [psi(*a)[0], psi(*b)[0], psi(*(0.5 * (a + b)))[0]]
        scale = max(abs(values[0]), abs(values[1]), 1.0)
        convexity = max(convexity, (values[2] - 0.5 * (values[0] + values[1])) / scale)
        s, delta = a
        value, p, upper = psi(s, delta)
        at_cap += p > 0.999 * upper
        monotone = max(monotone, (value - psi(s + 0.5, delta)[0]) / max(abs(value), 1.0))
    pure_compaction = psi(0.0, -1.0)[0]
    return dict(
        worst_midpoint_excess=float(convexity),
        worst_decrease_in_s=float(monotone),
        fraction_at_bound=at_cap / trials,
        psi_pure_compaction=float(pure_compaction),
    )


class Point:
    """Critical-state CIDR at one quadrature point, with phi (hence C) fixed."""

    def __init__(self, C=4e4, k=2 * 0.003 * np.sqrt(2500.0)):
        self.m, self.C, self.k = Jop(), C, k

    def I(self, s, p):
        return self.k * s / np.sqrt(p)

    def g(self, s, p):
        I = self.I(s, p)
        return 2 * p * s * ((5 * self.m.alpha(I) - 4 * self.m.mu(I)) / 3 - p / self.C)

    def Y(self, s, p):
        return self.m.alpha(self.I(s, p)) * p - p**2 / self.C

    def gp(self, s, p):
        return 2 * (self.m.beta(self.I(s, p)) - 2 * p / self.C) * s

    def p_cap(self, s):
        """The fixed cap alpha(0) C."""
        return self.m.alpha(0.0) * self.C

    def optimum(self, s, delta):
        """p* maximising g - p delta on [0, p_cap], and its branch."""
        cap = self.p_cap(s)
        if self.gp(s, cap) - delta >= 0:
            return cap, "cap"
        if s == 0 or self.gp(s, 1e-12 * cap) - delta <= 0:
            return 0.0, "zero"
        return (
            optimize.brentq(lambda p: self.gp(s, p) - delta, 1e-12 * cap, cap, xtol=1e-14 * cap),
            "interior",
        )

    def psi(self, s, delta):
        p, _ = self.optimum(s, delta)
        return (self.g(s, p) if p > 0 else 0.0) - p * delta


def states():
    """Steady flow, critical state and the cap at rest, at fixed C."""
    m, C = Jop(), 4e4
    Is = np.logspace(-5, 0, 30)
    steady = [
        (m.alpha(I) * (m.beta(I) * C / 2) - (m.beta(I) * C / 2) ** 2 / C) / (m.beta(I) * C / 2) - m.mu(I)
        for I in Is
    ]
    return dict(
        steady_flow_recovers_mu_I=float(np.max(np.abs(steady))),
        critical_pressure_over_C_mu1=float(m.beta(1e-12) * C / 2 / (C * m.mu1)),
        cap_over_critical_pressure=float(m.alpha(1e-12) * C / (m.beta(1e-12) * C / 2)),
    )


def packing_scale(deficit=0.005, p_ref=1500 * 9.81 * 0.1, phi_min=0.05, phi_max=0.65):
    """Fannon et al.'s C = (P/b) x, x = psi/(1-psi): choose P/b from a packing deficit at p_ref (on the cap)."""
    m = Jop()
    Delta = phi_max - phi_min
    c0 = 2 * m.mu1  # cap at rest: p = 2 mu1 C
    x_ref = Delta / deficit - 1
    scale = p_ref / (c0 * x_ref)  # P / b

    def phi_cap(p):
        x = p / (c0 * scale)
        return phi_min + Delta * x / (1 + x)

    h = 1e-6 * p_ref
    dphi = (phi_cap(p_ref + h) - phi_cap(p_ref - h)) / (2 * h)
    K = phi_cap(p_ref) / dphi
    return dict(
        P_over_b=scale,
        deficit_at_p_ref=float(phi_max - phi_cap(p_ref)),
        K_over_phi_p_over_e=float(K / (phi_cap(p_ref) * p_ref / deficit)),
        deficit_at_tenth_p_ref=float(phi_max - phi_cap(0.1 * p_ref)),
    )


def derivatives(seed=2, trials=30):
    """psi's gradient and Hessian against finite differences, on each branch."""
    q = Point()
    rng = np.random.default_rng(seed)
    worst = dict(interior_grad=0.0, interior_hess=0.0, cap_grad=0.0, zero_grad=0.0, interior_det=0.0)
    seen = {"interior": 0, "cap": 0, "zero": 0}

    def fd(fun, x, h):
        return (fun(x + h) - fun(x - h)) / (2 * h)

    for _ in range(trials * 4):
        s, delta = rng.uniform(0.2, 10), rng.uniform(-3, 3)
        p, branch = q.optimum(s, delta)
        if seen[branch] >= trials:
            continue
        seen[branch] += 1
        hs, hd = 1e-4 * s, 1e-4 * max(abs(delta), 1)
        psi_s = fd(lambda x: q.psi(x, delta), s, hs)
        psi_d = fd(lambda x: q.psi(s, x), delta, hd)
        if branch == "interior":
            grad = (2 * q.Y(s, p), -p)
            # g second derivatives by finite differences of g.
            e = 1e-4
            g_ss = (q.g(s * (1 + e), p) - 2 * q.g(s, p) + q.g(s * (1 - e), p)) / (s * e) ** 2
            g_pp = (q.g(s, p * (1 + e)) - 2 * q.g(s, p) + q.g(s, p * (1 - e))) / (p * e) ** 2
            g_sp = (
                q.g(s * (1 + e), p * (1 + e))
                - q.g(s * (1 + e), p * (1 - e))
                - q.g(s * (1 - e), p * (1 + e))
                + q.g(s * (1 - e), p * (1 - e))
            ) / (4 * s * p * e * e)
            hess = np.array([[g_ss - g_sp**2 / g_pp, g_sp / g_pp], [g_sp / g_pp, -1 / g_pp]])
            numeric = np.array(
                [
                    [
                        fd(lambda x: fd(lambda y: q.psi(y, delta), x, hs), s, hs),
                        fd(lambda x: fd(lambda y: q.psi(x, y), delta, hd), s, hs),
                    ],
                    [0, fd(lambda x: fd(lambda y: q.psi(s, y), x, hd), delta, hd)],
                ]
            )
            numeric[1, 0] = numeric[0, 1]
            worst["interior_hess"] = max(
                worst["interior_hess"], float(np.abs(hess - numeric).max() / np.abs(hess).max())
            )
            det = np.linalg.det(hess)
            worst["interior_det"] = max(worst["interior_det"], abs(det - (-g_ss / g_pp)) / abs(det))
            key = "interior_grad"
        elif branch == "cap":
            grad = (2 * q.Y(s, p), -p)
            key = "cap_grad"
        else:
            grad = (0.0, 0.0)
            key = "zero_grad"
        scale = max(abs(grad[0]), abs(grad[1]), 1.0)
        worst[key] = max(worst[key], abs(psi_s - grad[0]) / scale, abs(psi_d - grad[1]) / scale)
    return dict(samples=seen, **{k: float(v) for k, v in worst.items()})


def main():
    for name, check in (
        ("A general", general),
        ("B inertial CIDR", inertial),
        ("C critical-state CIDR", critical_state),
        ("C-D numeric (Jop mu)", numeric),
        ("E cap", cap),
        ("F dissipation, fixed cap", lambda: dissipation("fixed")),
        ("F dissipation, shear-dependent cap", lambda: dissipation("shear")),
        ("F dissipation, no cap", lambda: dissipation(None)),
        # Grain-scale units (d = rho_s = 1, C = 1), where I reaches order one.
        ("F dissipation, fixed cap, I ~ 1", lambda: dissipation("fixed", C=1.0, k=2.0)),
        ("F dissipation, shear-dependent cap, I ~ 1", lambda: dissipation("shear", C=1.0, k=2.0)),
        ("G states", states),
        ("H packing scale", packing_scale),
        ("I derivatives of psi", derivatives),
    ):
        print(name, check(), flush=True)


if __name__ == "__main__":
    main()
