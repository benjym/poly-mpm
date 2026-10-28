"""Input compatibility, defaults, and construction of the array solver."""

from importlib import import_module
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


def prepare(P, mode="simulation"):
    """Normalise a Params instance; also used by deterministic test fixtures."""
    defaults(
        P,
        t=0.0,
        tstep=0,
        grid_save=0,
        mp_save=0,
        M_tol=1e-10,
        max_g=0.0,
        max_q=0.0,
        damping=False,
        theta=0.0,
        pressure=None,
        initial_flow=False,
        smooth_gamma_dot=False,
        smooth_grad2=False,
        diagnostics_stride=1,
        timestep_policy="report",
        mechanics="updated",
        interpolation="linear",
        transfer="flip",
        pic_fraction=0.0,
        diagnostics=True,
        normalise_phi=False,
        time_stepping="static",
        CFL=0.1,
        c=0.0,
        l=10.0,
    )
    if hasattr(P, "segregate") and hasattr(P, "segregate_grid") and P.segregate != P.segregate_grid:
        raise ValueError("segregate and segregate_grid disagree")
    defaults(P, segregate=getattr(P, "segregate_grid", False))
    if getattr(P, "segregate_mp", False):
        raise ValueError("legacy segregate_mp is unsupported; use grid segregation")
    if int(P.diagnostics_stride) != P.diagnostics_stride or P.diagnostics_stride < 1:
        raise ValueError("diagnostics_stride must be a positive integer")
    if P.timestep_policy not in ("report", "enforce"):
        raise ValueError("timestep_policy must be report or enforce")
    if P.mechanics not in ("updated", "legacy"):
        raise ValueError("mechanics must be updated or legacy")
    if P.interpolation not in ("linear", "quadratic"):
        raise ValueError("interpolation must be linear or quadratic")
    if P.transfer not in ("flip", "pic", "apic"):
        raise ValueError("transfer must be flip, pic, or apic")
    if not np.isfinite(P.pic_fraction) or not 0 <= P.pic_fraction <= 1:
        raise ValueError("pic_fraction must be between zero and one")
    if P.transfer != "flip" and P.pic_fraction:
        raise ValueError("pic_fraction is only used with flip")
    if P.transfer == "apic" and P.interpolation != "quadratic":
        raise ValueError(
            "APIC requires quadratic interpolation; linear APIC has degenerate-support instabilities"
        )
    P.mode = mode
    P.S = P.S if isinstance(P.S, list) else [P.S]
    if not P.S:
        raise ValueError("at least one material is required")
    P.phases = len(P.S)
    for S in P.S:
        get_law(S.law)
    if getattr(getattr(P, "F", None), "n", 0):
        raise ValueError(
            f"Unsupported material law {getattr(P.F, 'law', 'legacy fluid')!r} in legacy P.F; use supported P.S materials"
        )
    if (
        not np.isfinite([P.dt, P.t, P.t_f, P.savetime, P.CFL, P.M_tol]).all()
        or P.dt <= 0
        or P.t_f < P.t
        or P.savetime <= 0
        or P.CFL <= 0
        or P.M_tol < 0
    ):
        raise ValueError("dt and savetime must be positive, with t_f >= t")
    defaults(P, nt=int(np.ceil(P.t_f / P.dt)), g=P.max_g, q_v=P.max_q, q_h=P.max_q)
    for name in (
        "wall",
        "two_walls",
        "roughness",
        "wall_mu",
        "no_slip_bottom",
        "cyclic_lr",
        "has_top",
        "has_bottom",
        "has_right",
        "has_left",
        "box",
        "outlet_left",
        "outlet_bottom",
        "inlet_right",
        "inlet_top",
        "force_boundaries",
        "vertical_force",
        "horizontal_force",
        "conveyor",
        "silo_left",
        "silo_bottom",
    ):
        defaults(P.B, **{name: False})
    defaults(P.G, thickness=1.0, ns=1, s=[0.001])
    P.G.s = np.asarray(P.G.s, dtype=float)
    if P.G.s.shape != (P.G.ns,) or not np.isfinite(P.G.s).all() or np.any(P.G.s <= 0):
        raise ValueError("s must contain ns positive finite grain sizes")
    defaults(P.G, s_m=P.G.s.min(), s_M=P.G.s.max())
    defaults(P.G, s_bar_0=P.G.s.mean())
    defaults(P, save_dir=str(getattr(P, "supername", f"im/{mode}/{P.S[0].law}")) + "/")
    for name in (
        "check_positions",
        "measure_stiffness",
        "measure_energy",
        "plot_gsd_mp",
        "plot_gsd_grid",
        "plot_gsd_debug",
        "plot_material_points",
        "plot_continuum",
        "save_s_bar",
        "save_density",
        "save_u",
        "save_phi_MP",
    ):
        defaults(P.O, **{name: False})
    defaults(P.O, save_diagnostics=True)
    defaults(P, update_forces=lambda: None)
    P.update_forces()
    G = Grid(P)
    L = ParticleSystem(P, G)
    defaults(P, update_timestep=update_timestep)
    if P.time_stepping == "dynamic" and P.mechanics == "updated":
        from stability import limits

        original_controller = P.update_timestep

        def controlled_step(P, G):
            original_controller(P, G)
            restrictions = limits(P, G, L)
            if restrictions:
                P.dt = min(P.dt, min(restrictions.values()))

        P.update_timestep = controlled_step
    from diagnostics import Diagnostics

    L.diagnostics = Diagnostics(P, L) if P.diagnostics else None
    from plotting import install_callbacks

    install_callbacks(P)
    return P, G, L


def get_parameters(params):
    if not params:
        raise ValueError("usage: python MPM.py <example> [arguments...]")
    module = import_module("inputs." + params[0])
    return prepare(module.Params(params), params[0])
