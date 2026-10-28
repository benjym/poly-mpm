"""Backend defaults and supported mechanics, without loading native libraries."""

import math

from . import ValidationPending


def defaults(P):
    backend = getattr(P, "backend", "cpu")
    if backend not in ("cpu", "metal"):
        raise ValueError("backend must be cpu or metal")
    P.backend = backend
    if backend == "metal":
        if not hasattr(P, "granular_tolerance"):
            P.granular_tolerance = 1e-4
        if not hasattr(P, "granular_solver"):
            P.granular_solver = "velocity"


def validate(P):
    """Reject unsupported configurations explicitly, before native initialization."""
    if P.backend == "cpu":
        return
    required = dict(
        solver="granular",
        mechanics="updated",
        granular_solver="velocity",
        granular_stress="intrinsic",
        granular_stress_basis="nodal",
        interpolation="linear",
        quadrature="point",
    )
    for name, expected in required.items():
        if getattr(P, name) != expected:
            raise ValueError(f"Metal requires {name}={expected!r}")
    if not math.isfinite(P.granular_tolerance) or P.granular_tolerance < 1e-5:
        raise ValueError("Metal FP32 requires granular_tolerance >= 1e-5")
    if len(P.S) != 1 or P.S[0].law != "granular":
        raise ValueError("Metal requires one granular material")
    if P.transfer not in ("flip", "pic", "apic_gradient"):
        raise ValueError("Metal supports flip, pic and apic_gradient transfers")
    if P.granular_viscosity != 0:
        raise ValueError("Metal does not support granular_viscosity")
    if P.segregate:
        raise ValueError("Metal segregation is not implemented")
    if P.time_stepping not in ("static", "dynamic"):
        raise ValueError("Metal supports static or built-in dynamic time stepping")
    for name in (
        "outlet_left",
        "outlet_bottom",
        "inlet_right",
        "inlet_top",
        "conveyor",
        "force_boundaries",
        "vertical_force",
        "horizontal_force",
        "wall_mu",
        "silo_left",
        "silo_bottom",
    ):
        if getattr(P.B, name, False):
            raise ValueError(f"Metal does not support boundary flag {name}")
    for name in ("packing_correction", "density_correction", "granular_density_correction"):
        if getattr(P, name, False):
            raise ValueError(f"Metal does not support {name}")


def require_runner():
    # Keep backend selection from accidentally running the CPU granular solver.
    raise ValidationPending(
        "The resident Metal MPM runner is not implemented yet. "
        "First validate the runtime on Apple silicon with "
        "python -m mpm_gpu.validate --stage runtime --benchmark --output metal-validation.json; "
        "see docs/gpu.md for the remaining hardware gates. No CPU fallback was used."
    )
