"""Conservative explicit timestep estimates, reported separately by mechanism."""

import numpy as np


def _state_limits(P, G, L):
    """Unscaled restrictions depending on particle/grid state, not the timestep."""
    h = min(G.dx, G.dy)
    inverse_h2 = 1 / G.dx**2 + 1 / G.dy**2
    result = {}
    speed = max((np.linalg.norm(b.v, axis=1).max(initial=0) for b in L.materials), default=0.0)
    if speed:
        result["advection"] = h / speed
    wave = viscosity = 0.0
    for b in L.materials:
        if not len(b) or b.material.law == "rigid":
            continue
        S = b.material
        rho = b.rho.min()
        shear = getattr(S, "G", 0.0) if S.law in ("elastic", "von_mises", "dp", "dp_rate") else 0.0
        if S.law == "elastic" and P.mechanics == "updated":
            # Frobenius/operator-norm bound on the reference neo-Hookean
            # Hessian, pushed to spatial wave speed by ||F||_2^2/rho0.
            # Unlike a small-strain speed this responds to compression/stretch.
            singular = np.linalg.svd(b.F, compute_uv=False)
            f2 = np.sum(singular**2, axis=1)
            inverse2 = np.sum(1 / singular**2, axis=1)
            inverse_operator2 = 1 / np.min(singular, axis=1) ** 2
            hessian = (
                S.G
                * b.J ** (-2 / 3)
                * (1 + 4 / 3 * np.sqrt(f2 * inverse2) + 2 / 9 * f2 * inverse2 + f2 * inverse_operator2 / 3)
            )
            hessian += S.K * (inverse2 + np.abs(np.log(b.J)) * inverse_operator2)
            hessian += np.linalg.norm(b.reference_stress, axis=(1, 2))
            speed2 = hessian * np.max(singular, axis=1) ** 2 / b.initial_rho
            wave = max(wave, np.sqrt(speed2.max()))
        else:
            wave = max(wave, np.sqrt((S.K + 4 * shear / 3) / rho))
        eta = (
            getattr(S, "eta_max", 0.0)
            if S.law in ("pouliquen", "linear_mu")
            else (S.mu_0 if S.law == "bingham" else getattr(S, "mu_s", 0.0) if S.law == "viscous" else 0.0)
        )
        bulk = getattr(S, "mu_v", 0.0) if S.law in ("viscous", "bingham") else 0.0
        viscosity = max(viscosity, (2 * eta + 3 * bulk) / rho)
    if wave:
        result["elastic_wave"] = h / wave
    if viscosity:
        result["viscosity"] = 1 / (2 * viscosity * inverse_h2)
    if P.segregate:
        rate = P.c * (P.G.s_M / P.G.s_m - 1) * np.linalg.norm(G.grad_pk, axis=1).max(initial=0)
        if rate:
            result["segregation_advection"] = h / rate
        from integrators import _diffusivity

        diffusion = _diffusivity(G, P, G.phi).max(initial=0)
        if diffusion:
            result["segregation_diffusion"] = 1 / (2 * diffusion * inverse_h2)
    return result


def limits(P, G, L, *, cache=None):
    """Estimate limits, optionally sharing state work within one runner phase.

    The runner clears its local cache before advancing the state. Standalone
    callers always recompute. Material callbacks are deliberately reevaluated:
    they may depend on the trial timestep, unlike the central estimates.
    """
    if cache is None:
        state = _state_limits(P, G, L)
    else:
        if "state" not in cache:
            cache["state"] = _state_limits(P, G, L)
        state = cache["state"]
    result = {"advection": state["advection"]} if "advection" in state else {}
    for b in L.materials:
        if len(b) and b.material.law != "rigid" and hasattr(b.material, "critical_time"):
            value = float(b.material.critical_time(P))
            if not np.isfinite(value) or value <= 0:
                raise ValueError("material critical_time must be positive and finite")
            result[f"material_{b.material_id}"] = value
    result.update((name, value) for name, value in state.items() if name != "advection")
    if any(not np.isfinite(x) or x <= 0 for x in result.values()):
        raise FloatingPointError("nonfinite or nonpositive explicit timestep estimate")
    return {name: float(value * P.CFL) for name, value in result.items()}
