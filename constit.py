"""Array constitutive kernels. Cauchy stress uses tension-positive signs.

Experimental scalar implementations are retained in legacy/constit_scalar.py.
Kernels update history and dstress; update_stress selects the integration scheme.
"""

import numpy as np

EYE = np.eye(3)
SUPPORTED = ("rigid", "elastic", "von_mises", "viscous", "pouliquen", "bingham", "dp", "dp_rate", "linear_mu")
ALGEBRAIC = ("bingham", "linear_mu")
FLUIDS = ("viscous", "bingham", "pouliquen", "linear_mu")
COROTATIONAL = ("dp", "dp_rate")
RETIRED = {
    "elastic_time": "Retired: arbitrary per-step mass deletion is not a constitutive law; use elastic for constant stiffness.",
    "ken_kamrin": "Retired: signed-pressure/viscosity errors in an unfinished duplicate; use pouliquen.",
    "ken_simple": "Retired: stale shear rate and hardcoded density cutoff; use linear_mu with b=0 for constant friction.",
}


def validate_material(material, mechanics="updated"):
    """Physical parameter domain for the newly restored laws (scalar SI inputs)."""
    name = material.law
    if name == "rigid" or (mechanics == "legacy" and name not in ALGEBRAIC + COROTATIONAL):
        return

    def scalar(key, lower=0.0, strict=False, default=None):
        value = getattr(material, key, default)
        if (
            value is None
            or np.ndim(value) != 0
            or not np.isfinite(value)
            or (value <= lower if strict else value < lower)
        ):
            comparison = ">" if strict else ">="
            raise ValueError(f"{name}: {key} must be a finite scalar {comparison} {lower}")
        return value

    scalar("K", strict=name in COROTATIONAL + ("elastic", "von_mises"))
    if name in ("elastic", "von_mises"):
        scalar("G", strict=True)
        if name == "von_mises":
            if hasattr(material, "yield_stress"):
                scalar("yield_stress", strict=True)
            else:
                scalar("k", strict=True)
        return
    if name == "viscous":
        scalar("mu_s")
        scalar("mu_v")
        return
    if name == "pouliquen":
        scalar("rho_s", strict=True)
        scalar("mu_0")
        scalar("delta_mu")
        scalar("I_0", strict=True)
        scalar("eta_max", strict=True)
        return
    if name in COROTATIONAL:
        scalar("G", strict=True)
        mu = scalar("mu", strict=True)
        if scalar("beta") > mu:
            raise ValueError(f"{name}: beta must satisfy 0 <= beta <= mu for nonnegative plastic dissipation")
        if name == "dp_rate":
            scalar("t_star")
    elif name == "bingham":
        mu_s = scalar("mu_s")
        tau = scalar("tau_0")
        mu_0 = scalar("mu_0", lower=mu_s, strict=tau > 0)
        scalar("mu_v", default=0.0)
        if hasattr(material, "gamma_c"):
            expected = tau / (mu_0 - mu_s) if tau else 0.0
            if not np.isclose(scalar("gamma_c"), expected, rtol=1e-10, atol=1e-14):
                raise ValueError(
                    "bingham: gamma_c must equal tau_0/(mu_0-mu_s); omit it to derive the transition"
                )
    else:
        scalar("mu_0")
        scalar("b")
        scalar("rho_s", strict=True)
        scalar("eta_max", strict=True)


def invariants(a):
    trace = np.trace(a, axis1=-2, axis2=-1)
    return trace, a - trace[:, None, None] * EYE / 3


def rigid(mp, material, params):
    mp.dstress.fill(0)


def elastic(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    mp.dstress = material.K * trace[:, None, None] * EYE + 2 * material.G * dev
    mp.gammadot = np.linalg.norm(dev, axis=(1, 2)) / params.dt
    mp.pressure = np.trace(mp.stress, axis1=1, axis2=2) / 3


def viscous(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    mp.de_kk = trace / 3
    mp.de_ij = dev
    mp.dp = material.K * mp.de_kk
    mp.pressure += mp.dp
    mp.dev_stress = 2 * material.mu_s * dev / params.dt
    volumetric = material.mu_v * mp.de_kk[:, None, None] * EYE / params.dt
    mp.dstress = mp.pressure[:, None, None] * EYE + mp.dev_stress - volumetric - mp.stress
    mp.gammadot = np.linalg.norm(dev, axis=(1, 2)) / params.dt


def von_mises(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    work = np.abs(np.sum(mp.dev_stress * dev, axis=(1, 2)))
    norm = np.linalg.norm(mp.dev_stress, axis=(1, 2))
    # For s < 2 the legacy expression is singular at zero stress. Its product
    # with the zero work/stress is defined here by its zero-increment limit.
    factor = np.zeros_like(norm)
    active = norm > 0
    factor[active] = (
        work[active] / (2 * material.k**2) * (norm[active] / (np.sqrt(2) * material.k)) ** (material.s - 2)
    )
    mp.dev_dstress = 2 * material.G * (dev - factor[:, None, None] * mp.dev_stress)
    mp.sigma_kk += 3 * material.K * trace
    mp.dev_stress += mp.dev_dstress
    mp.dstress = mp.dev_stress + mp.sigma_kk[:, None, None] * EYE / 3 - mp.stress
    mp.yieldfunction = norm / (np.sqrt(2) * material.k) - 1
    mp.gammadot = np.linalg.norm(dev, axis=(1, 2)) / params.dt
    mp.pressure = np.trace(mp.stress, axis1=1, axis2=2) / 3


def pouliquen(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    mp.de_kk = trace / 3
    mp.de_ij = dev
    mp.gammadot = np.linalg.norm(2 * dev / params.dt, axis=(1, 2))
    size = mp.phi @ params.G.s
    pressure = np.abs(mp.pressure)
    mp.I = np.zeros_like(pressure)
    positive = pressure > 0
    mp.I[positive] = mp.gammadot[positive] * size[positive] * np.sqrt(material.rho_s / pressure[positive])
    # I is undefined at zero pressure: store zero for that diagnostic and
    # evaluate friction using its infinite-I limit without a huge sentinel.
    ratio = np.divide(mp.I, material.I_0 + mp.I, out=np.zeros_like(mp.I), where=mp.I > 0)
    ratio[(~positive) & (mp.gammadot > 0)] = 1
    mp.mu = material.mu_0 + material.delta_mu * ratio
    numerator = 2 * np.sqrt(2) * mp.mu * pressure
    mp.eta = np.divide(numerator, mp.gammadot, out=np.zeros_like(pressure), where=mp.gammadot > 0)
    mp.eta[(mp.gammadot == 0) & (numerator > 0)] = np.finfo(float).max
    mp.eta_limited = np.minimum(mp.eta, material.eta_max)
    mp.dev_stress = mp.eta_limited[:, None, None] * dev / params.dt
    mp.dp = material.K * mp.de_kk
    mp.pressure += mp.dp
    cutoff = (mp.pressure > 0) | (mp.rho < mp.initial_rho * 0.8)
    mp.dp[cutoff] -= mp.pressure[cutoff]
    mp.pressure[cutoff] = 0
    mp.dstress = mp.pressure[:, None, None] * EYE + mp.dev_stress - mp.stress


def _flow_increment(mp, params):
    trace, dev = invariants(mp.dstrain)
    mp.de_kk = trace / 3
    mp.de_ij = dev
    mp.gammadot = np.sqrt(2.0) * np.linalg.norm(dev, axis=(1, 2)) / params.dt
    return trace, dev / params.dt


def bingham(mp, material, params):
    """Bi-viscosity Bingham fluid; gamma = sqrt(2 D':D'), viscosities in Pa s."""
    trace, rate = _flow_increment(mp, params)
    mp.dp = material.K * trace
    mp.pressure += mp.dp
    mp.eta = np.full_like(mp.gammadot, material.mu_s)
    if material.tau_0:
        mp.eta.fill(material.mu_0)
        flowing = mp.gammadot > material.tau_0 / (material.mu_0 - material.mu_s)
        mp.eta[flowing] = material.mu_s + material.tau_0 / mp.gammadot[flowing]
    mp.eta_limited = mp.eta.copy()
    mp.dev_stress = 2 * mp.eta[:, None, None] * rate
    bulk = getattr(material, "mu_v", 0.0) * trace / params.dt
    mp.dstress = (mp.pressure + bulk)[:, None, None] * EYE + mp.dev_stress - mp.stress


def linear_mu(mp, material, params):
    """Dense flowing linear mu(I) law, with finite effective-viscosity cap."""
    trace, rate = _flow_increment(mp, params)
    old_pressure = mp.pressure.copy()
    mp.pressure = np.minimum(old_pressure + material.K * trace, 0.0)
    mp.dp = mp.pressure - old_pressure
    p = -mp.pressure
    size = mp.phi @ params.G.s
    mp.I = np.zeros_like(p)
    positive = p > 0
    mp.I[positive] = mp.gammadot[positive] * size[positive] * np.sqrt(material.rho_s) / np.sqrt(p[positive])
    mp.mu = material.mu_0 + material.b * mp.I
    # Evaluate mu*p/gamma as a sum to avoid zero-pressure and zero-rate singularities.
    mp.eta = np.divide(
        material.mu_0 * p,
        mp.gammadot,
        out=np.where(material.mu_0 * p > 0, material.eta_max, 0.0),
        where=mp.gammadot > 0,
    )
    mp.eta += material.b * size * np.sqrt(material.rho_s) * np.sqrt(p)
    mp.eta = np.minimum(mp.eta, material.eta_max)
    mp.eta[~positive] = 0.0
    mp.eta_limited = mp.eta.copy()
    mp.dev_stress = 2 * mp.eta[:, None, None] * rate
    mp.dstress = mp.pressure[:, None, None] * EYE + mp.dev_stress - mp.stress


def _drucker_prager(mp, material, params, t_star):
    """Backward Euler for f=q-mu*p and g=q-beta*p; no cohesion or tension.

    p=-tr(sigma)/3, q=sqrt(3/2 s:s). Perzyna viscosity is 3 G t_star.
    A separate associated mean-tension cutoff handles the cone apex without
    dividing by p or q. The rate model permits viscous shear overstress at p=0.
    """
    trace, dev = invariants(mp.dstrain)
    mp.de_kk = trace / 3
    mp.de_ij = dev
    trial = mp.stress + material.K * trace[:, None, None] * EYE + 2 * material.G * dev
    tr, shear = invariants(trial)
    p_trial = -tr / 3
    q_trial = np.sqrt(1.5) * np.linalg.norm(shear, axis=(1, 2))
    viscous = 3 * material.G * t_star / params.dt
    multiplier = np.maximum(q_trial - material.mu * p_trial, 0.0) / (
        3 * material.G + material.K * material.mu * material.beta + viscous
    )
    p = p_trial + material.K * material.beta * multiplier
    tension = p < 0
    multiplier[tension] = q_trial[tension] / (3 * material.G + viscous)
    p[tension] = 0.0
    q = np.maximum(q_trial - 3 * material.G * multiplier, 0.0)
    scale = np.divide(q, q_trial, out=np.zeros_like(q), where=q_trial > 0)
    mp.dev_stress = shear * scale[:, None, None]
    mp.pressure = -p
    mp.dp = mp.pressure - np.trace(mp.stress, axis1=1, axis2=2) / 3
    mp.sigma_kk = -3 * p
    mp.mu.fill(material.mu)
    mp.gammadot = np.sqrt(2.0) * np.linalg.norm(dev, axis=(1, 2)) / params.dt
    mp.yieldfunction = q - material.mu * p
    mp.dstress = -p[:, None, None] * EYE + mp.dev_stress - mp.stress


def dp(mp, material, params):
    """Cohesionless, perfectly plastic Drucker-Prager solid."""
    _drucker_prager(mp, material, params, 0.0)


def dp_rate(mp, material, params):
    """Linear Perzyna Drucker-Prager solid; t_star is a relaxation time (s)."""
    _drucker_prager(mp, material, params, material.t_star)


def _half_rotation(mp, dt):
    """Batched Rodrigues exponential of half the incremental spin."""
    spin = (mp.L - mp.L_T) * dt / 4
    angle = np.linalg.norm(spin, axis=(1, 2)) / np.sqrt(2.0)
    a = np.sinc(angle / np.pi)
    b = 0.5 * np.sinc(angle / (2 * np.pi)) ** 2
    return EYE + a[:, None, None] * spin + b[:, None, None] * (spin @ spin)


def update_stress(mp, material, params):
    """Integrate restored laws consistently with their stress/rate definition.

    Keep the original five models' legacy integration for compatibility. New
    fluids prescribe current stress directly. New solids use a symmetric spin
    split for a corotational hypoelastic/plastic update (exact for pure spin).
    """
    if material.law not in ALGEBRAIC + COROTATIONAL:
        mp.law(mp, material, params)
        mp.stress += mp.dstress - (mp.L @ mp.stress + mp.stress @ mp.L_T) * params.dt
        return
    old = mp.stress.copy()
    old_dev = old - np.trace(old, axis1=1, axis2=2)[:, None, None] * EYE / 3
    if material.law in COROTATIONAL:
        rotation = _half_rotation(mp, params.dt)
        mp.stress = rotation @ old @ rotation.swapaxes(1, 2)
        mp.law(mp, material, params)
        mp.stress = rotation @ (mp.stress + mp.dstress) @ rotation.swapaxes(1, 2)
    elif material.law in ALGEBRAIC:
        mp.law(mp, material, params)
        mp.stress += mp.dstress
    mp.dstress = mp.stress - old
    # Fluid pressure is the elastic mean part: bulk viscosity is not history.
    if material.law in COROTATIONAL:
        mp.pressure = np.trace(mp.stress, axis1=1, axis2=2) / 3
    mean = np.trace(mp.stress, axis1=1, axis2=2) / 3
    mp.sigma_kk = 3 * mean
    mp.dev_stress = mp.stress - mean[:, None, None] * EYE
    mp.dev_dstress = mp.dev_stress - old_dev


LEGACY_REGISTRY = {name: globals()[name] for name in SUPPORTED}
_legacy_update_stress = update_stress


def elastic(mp, material, params):
    """Compressible neo-Hookean Cauchy stress, with reference prestress.

    Reference energy density::

        W0 = G/2 (J^(-2/3) tr(F F^T)-3) + K/2 (ln J)^2
        W0 += S0:(F^T F-I)/2

    The reference prestress S0 is the input stress at F=I. The law recovers
    isotropic linear elasticity for small increments around an unstressed state.
    """
    J = mp.J
    logJ = np.log(J)
    b = mp.F @ mp.F.swapaxes(1, 2)
    trace_b = np.trace(b, axis1=1, axis2=2)
    scale = J ** (-2 / 3)
    tau = (
        material.G * scale[:, None, None] * (b - trace_b[:, None, None] * EYE / 3)
        + material.K * logJ[:, None, None] * EYE
    )
    tau += mp.F @ mp.reference_stress @ mp.F.swapaxes(1, 2)
    stress = tau / J[:, None, None]
    mp.dstress = stress - mp.stress
    green = (mp.F.swapaxes(1, 2) @ mp.F - EYE) / 2
    energy = (
        material.G / 2 * (scale * trace_b - 3)
        + material.K / 2 * logJ**2
        + np.sum(mp.reference_stress * green, axis=(1, 2))
    )
    mp.elastic_energy = mp.V0 * energy
    _, dev = invariants(mp.dstrain)
    mp.gammadot = np.sqrt(2) * np.linalg.norm(dev, axis=(1, 2)) / params.dt


def viscous(mp, material, params):
    """Newtonian fluid with positive bulk dissipation and standard bulk modulus."""
    trace, rate = _flow_increment(mp, params)
    mp.dp = material.K * trace
    mp.pressure += mp.dp
    mp.eta.fill(material.mu_s)
    mp.eta_limited = mp.eta.copy()
    mp.dev_stress = 2 * material.mu_s * rate
    bulk = material.mu_v * trace / params.dt
    mp.dstress = (mp.pressure + bulk)[:, None, None] * EYE + mp.dev_stress - mp.stress


def pouliquen(mp, material, params):
    """Bounded mu(I) flowing law using engineering shear rate and true density.

    No grid-occupancy density cutoff. The finite cap regularises the zero-rate
    limit. At zero confinement, I is undefined (stored as zero), mu uses its
    infinite-I limit for nonzero flow, and frictional stress vanishes.
    """
    trace, rate = _flow_increment(mp, params)
    old_pressure = mp.pressure.copy()
    mp.pressure = np.minimum(old_pressure + material.K * trace, 0.0)
    mp.dp = mp.pressure - old_pressure
    p = -mp.pressure
    size = mp.phi @ params.G.s
    mp.I.fill(0.0)
    positive = p > 0
    mp.I[positive] = mp.gammadot[positive] * size[positive] * np.sqrt(material.rho_s) / np.sqrt(p[positive])
    ratio = np.divide(mp.I, material.I_0 + mp.I, out=np.zeros_like(p), where=mp.I > 0)
    ratio[(~positive) & (mp.gammadot > 0)] = 1.0
    mp.mu = material.mu_0 + material.delta_mu * ratio
    mp.eta = np.divide(
        mp.mu * p, mp.gammadot, out=np.where(mp.mu * p > 0, material.eta_max, 0.0), where=mp.gammadot > 0
    )
    mp.eta = np.minimum(mp.eta, material.eta_max)
    mp.eta_limited = mp.eta.copy()
    mp.dev_stress = 2 * mp.eta[:, None, None] * rate
    mp.dstress = mp.pressure[:, None, None] * EYE + mp.dev_stress - mp.stress


def von_mises(mp, material, params):
    """Perfect J2 plasticity; yield_stress (Pa), or sqrt(3)*legacy k."""
    trace, dev = invariants(mp.dstrain)
    trial = mp.stress + material.K * trace[:, None, None] * EYE + 2 * material.G * dev
    tr, shear = invariants(trial)
    q = np.sqrt(1.5) * np.linalg.norm(shear, axis=(1, 2))
    strength = getattr(material, "yield_stress", np.sqrt(3) * getattr(material, "k", 0.0))
    multiplier = np.maximum(q - strength, 0) / (3 * material.G)
    scale = np.divide(np.minimum(q, strength), q, out=np.zeros_like(q), where=q > 0)
    stress = tr[:, None, None] * EYE / 3 + scale[:, None, None] * shear
    mp.dstress = stress - mp.stress
    mp.yieldfunction = np.minimum(q, strength) - strength
    mp.gammadot = np.sqrt(2) * np.linalg.norm(dev, axis=(1, 2)) / params.dt


def _compliance(stress, material):
    trace, dev = invariants(stress)
    return trace[:, None, None] * EYE / (9 * material.K) + dev / (2 * material.G)


def update_stress(mp, material, params):
    """Objective updated mechanics; historical integration is opt-in legacy mode."""
    if params.mechanics == "legacy":
        _legacy_update_stress(mp, material, params)
        return
    old = mp.stress.copy()
    old_bulk_pressure = mp.pressure.copy()
    _, old_dev = invariants(old)
    plastic = material.law in ("dp", "dp_rate", "von_mises")
    if plastic:
        R = _half_rotation(mp, params.dt)
        mp.stress = R @ old @ R.swapaxes(1, 2)
        mp.plastic_strain = R @ mp.plastic_strain @ R.swapaxes(1, 2)
        before = mp.stress.copy()
    mp.law(mp, material, params)
    result = mp.stress + mp.dstress
    if plastic:
        dep = mp.dstrain - _compliance(result - before, material)
        dissipation = np.sum(result * dep, axis=(1, 2)) * mp.V
        tolerance = 1e-10 * np.maximum(1.0, np.linalg.norm(result, axis=(1, 2)) * mp.V)
        if np.any(dissipation < -tolerance):
            raise FloatingPointError("negative plastic dissipation")
        mp.plastic_work += np.maximum(dissipation, 0)
        mp.plastic_strain += dep
        mp.plastic_strain = R @ mp.plastic_strain @ R.swapaxes(1, 2)
        result = R @ result @ R.swapaxes(1, 2)
        mp.elastic_strain = _compliance(result, material)
        mp.elastic_energy = 0.5 * mp.V * np.sum(result * mp.elastic_strain, axis=(1, 2))
    elif material.law != "elastic":
        trace, dev = invariants(mp.dstrain)
        _, shear = invariants(result)
        bulk_viscosity = getattr(material, "mu_v", 0.0) if material.law in ("viscous", "bingham") else 0.0
        mp.viscous_work += mp.V * (np.sum(shear * dev, axis=(1, 2)) + bulk_viscosity * trace**2 / params.dt)
        # Integrate elastic bulk work p(J) dV exactly along this exponential
        # volume increment. No-tension flow laws have a clipped pressure path.
        old_volume = mp.V * np.exp(-trace)
        if material.K:
            start = np.zeros_like(trace)
            stop = trace.copy()
            if material.law in ("pouliquen", "linear_mu"):
                cutoff = -old_bulk_pressure / material.K
                start = np.minimum(start, cutoff)
                stop = np.minimum(stop, cutoff)

            def potential(t):
                remainder = np.where(
                    np.abs(t) < 1e-4,
                    t**2 / 2 + t**3 / 3 + t**4 / 8 + t**5 / 30 + t**6 / 144,
                    t * np.exp(t) - np.expm1(t),
                )
                return old_bulk_pressure * np.expm1(t) + material.K * remainder

            mp.elastic_energy += old_volume * (potential(stop) - potential(start))
        else:
            mp.elastic_energy += mp.pressure * (mp.V - old_volume)
    mp.stress = result
    mp.dstress = result - old
    mean, dev = invariants(result)
    if material.law in ("elastic", "dp", "dp_rate", "von_mises"):
        mp.pressure = mean / 3
    mp.sigma_kk = mean
    mp.dev_stress = dev
    mp.dev_dstress = dev - old_dev
    mp.material_work += np.sum(0.5 * (old + result) * mp.dstrain, axis=(1, 2)) * mp.V


REGISTRY = {name: globals()[name] for name in SUPPORTED}


def get_law(name, mechanics="updated"):
    try:
        return (LEGACY_REGISTRY if mechanics == "legacy" else REGISTRY)[name]
    except KeyError:
        reason = RETIRED.get(name, "No implemented and tested constitutive law with this name.")
        raise ValueError(
            f"Unsupported material law {name!r}. {reason} Supported: {', '.join(SUPPORTED)}."
        ) from None
