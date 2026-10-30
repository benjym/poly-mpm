"""Input compatibility, defaults, and construction of the array solver."""

from importlib import import_module
from inputs import validate_arguments
from settings import normalise
import numpy as np
from constit import get_law
from grid import Grid
from particles import ParticleSystem
from transfers import nodal_velocity


def defaults(obj, **values):
    for name, value in values.items():
        if not hasattr(obj, name):
            setattr(obj, name, value)


def update_timestep(P, G):
    distance = min(G.dx, G.dy)
    speed = np.max(np.abs(nodal_velocity(G, P)), initial=0)
    limits = []
    if speed > 0:
        limits.append(0.1 * distance / speed)
    if P.segregate:
        rate = P.c * (P.G.s_M / P.G.s_m - 1) * np.max(np.abs(G.grad_pk), initial=0)
        if rate > 0:
            limits.append(distance / rate)
    for S in P.S:
        if hasattr(S, "critical_time"):
            limits.append(S.critical_time(P))
    candidate = P.CFL * min(limits) if limits else P.dt
    if not np.isfinite(candidate) or candidate <= 0:
        raise ValueError(f"invalid dynamic timestep {candidate}")
    P.dt = candidate


def validate_cidr(P):
    """Critical-state CIDR (cidr.py): material fields and the supported configuration."""
    if P.granular_stress != "intrinsic" or P.granular_stress_basis != "nodal":
        raise ValueError("granular_solver='cidr' requires intrinsic stress on the nodal basis")
    S = P.S[0]
    fields = ("mu_s", "mu_d", "I_0", "grain_diameter", "phi_min", "phi_max", "cidr_deficit", "cidr_p_ref")
    missing = [name for name in fields if not hasattr(S, name)]
    if missing:
        raise ValueError(f"granular_solver='cidr' requires material fields {', '.join(missing)}")
    values = {name: float(getattr(S, name)) for name in fields}
    if not all(np.isfinite(v) and v > 0 for v in values.values()):
        raise ValueError("cidr material fields must be positive and finite")
    if not values["mu_s"] < values["mu_d"]:
        raise ValueError("cidr requires mu_s < mu_d")
    if not values["phi_min"] < values["phi_max"] - values["cidr_deficit"]:
        raise ValueError("cidr requires phi_min < phi_max - cidr_deficit")


def prepare(P, mode="simulation"):
    """Normalise a Params instance; also used by deterministic test fixtures."""
    from mpm_gpu.config import defaults as backend_defaults

    backend_defaults(P)
    normalise(P)
    defaults(P, t=0.0, tstep=0, grid_save=0, mp_save=0)
    if P.transfer == "apic_gradient" and P.solver != "granular":
        raise ValueError("apic_gradient is restricted to the granular source audit")
    if P.pic_fraction > 1:
        raise ValueError("pic_fraction must be between zero and one")
    if P.transfer != "flip" and P.pic_fraction:
        raise ValueError("pic_fraction is only used with flip")
    if P.fbar and P.solver != "explicit":
        raise ValueError("fbar requires the explicit solver")
    if P.transfer == "apic" and P.interpolation != "quadratic":
        raise ValueError(
            "APIC requires quadratic interpolation; linear APIC has degenerate-support instabilities"
        )
    if P.t_f < P.t:
        raise ValueError("t_f must not precede the start time")
    P.mode = mode
    P.S = P.S if isinstance(P.S, list) else [P.S]
    if not P.S:
        raise ValueError("at least one material is required")
    for S in P.S:
        get_law(S.law)
    if P.solver == "granular":
        if P.granular_stress_basis == "coarse":
            from granular import coarse_shape

            coarse_shape(P)
        if P.granular_stress == "plane_strain" and P.granular_solver == "block":
            raise ValueError("plane_strain requires the hybrid or velocity contact solver")
        if P.granular_solver == "cidr":
            validate_cidr(P)
        if P.quadrature != "point":
            raise ValueError("granular prototype requires point quadrature")
        if P.interpolation != "linear":
            raise ValueError("granular prototype currently requires linear interpolation")
        if len(P.S) != 1 or P.S[0].law != "granular":
            raise ValueError("granular prototype requires one granular material")
        if P.initial_stress is not None or P.force_scale != 1:
            raise ValueError("granular prototype solves stress directly and requires force_scale = 1")
        if P.granular_stopping == "increment" and P.granular_solver != "velocity":
            raise ValueError("granular_stopping='increment' requires granular_solver='velocity'")
        if P.B.conveyor:
            raise ValueError("granular prototype does not support conveyors")
        if P.B.wall_mu and P.granular_solver not in ("velocity", "cidr"):
            raise ValueError(
                "Coulomb walls (wall_contact = 'coulomb') require granular_solver='velocity' or 'cidr'"
            )
        if P.B.vertical_force or P.B.horizontal_force:
            raise ValueError("granular prototype does not support boundary tractions")
    elif any(S.law == "granular" for S in P.S):
        raise ValueError("granular material requires solver='granular'")
    from mpm_gpu.config import validate as validate_backend

    validate_backend(P)
    if not np.isfinite(P.G.s).all() or np.any(P.G.s <= 0):
        raise ValueError("G.grain_sizes must be positive and finite")
    P.G.s_m, P.G.s_M, P.G.s_bar_0 = P.G.s.min(), P.G.s.max(), P.G.s.mean()
    defaults(P, save_dir=f"im/{mode}/{P.S[0].law}/")
    P.update_forces()
    G = Grid(P)
    L = ParticleSystem(P, G)
    defaults(P, update_timestep=update_timestep)
    if P.time_stepping == "dynamic":
        from stability import limits

        original_controller = P.update_timestep

        def controlled_step(P, G, *, _limit_cache=None):
            original_controller(P, G)
            restrictions = limits(P, G, L, cache=_limit_cache)
            if restrictions:
                P.dt = min(P.dt, min(restrictions.values()))

        controlled_step._uses_limit_cache = True
        controlled_step._controller = original_controller
        P.update_timestep = controlled_step
    from diagnostics import Diagnostics

    L.diagnostics = Diagnostics(P, L) if P.diagnostics else None
    from plotting import install_callbacks

    install_callbacks(P)
    return P, G, L


def get_parameters(params):
    validate_arguments(params)
    module = import_module("inputs." + params[0])
    P = module.Params(params)
    P.args = list(params)
    return prepare(P, params[0])
