"""Backend defaults and supported mechanics, without loading native libraries."""

import math


def defaults(P):
    backend = getattr(P, "backend", "cpu")
    if backend not in ("cpu", "metal"):
        raise ValueError("backend must be cpu or metal")
    P.backend = backend
    if backend == "metal":
        if not hasattr(P, "granular_tolerance"):
            P.granular_tolerance = 1e-4
        if not hasattr(P, "granular_cidr_tolerance"):
            P.granular_cidr_tolerance = 1e-4


def validate(P):
    """Reject unsupported configurations explicitly, before native initialization."""
    if P.backend == "cpu":
        return
    if P.granular_solver not in ("velocity", "cidr"):
        raise ValueError("Metal requires granular_solver='velocity' or 'cidr'")
    if P.granular_solver == "cidr" and P.granular_cidr_tolerance < 1e-5:
        raise ValueError("Metal FP32 requires granular_cidr_tolerance >= 1e-5")
    required = dict(
        solver="granular",
        granular_stress="intrinsic",
        granular_stress_basis="nodal",
        interpolation="linear",
        quadrature="point",
    )
    for name, expected in required.items():
        if getattr(P, name) != expected:
            raise ValueError(f"Metal requires {name}={expected!r}")
    if P.granular_solver == "velocity" and (
        not math.isfinite(P.granular_tolerance) or P.granular_tolerance < 1e-5
    ):
        raise ValueError("Metal FP32 requires granular_tolerance >= 1e-5")
    if len(P.S) != 1 or P.S[0].law != "granular":
        raise ValueError("Metal requires one granular material")
    if P.transfer not in ("flip", "pic", "apic_gradient"):
        raise ValueError("Metal supports flip, pic and apic_gradient transfers")
    if getattr(P, "granular_stopping", "residual") != "residual":
        raise ValueError("Metal supports only granular_stopping='residual'")
    if P.granular_viscosity != 0:
        raise ValueError("Metal does not support granular_viscosity")
    if P.segregate and (P.smooth_shear_rate or P.G.ns > 8):
        raise ValueError("Metal segregation requires smooth_shear_rate=False and at most 8 grain sizes")
    if P.time_stepping not in ("static", "dynamic"):
        raise ValueError("Metal supports static or built-in dynamic time stepping")
    for name in ("conveyor_speed", "traction_vertical", "traction_horizontal"):
        if getattr(P.B, name) is not None:
            raise ValueError(f"Metal does not support B.{name}")
    if P.force_scale != 1:
        raise ValueError("Metal does not support force_scale")
    if P.O.check_positions:
        raise ValueError("Metal does not support O.check_positions")
