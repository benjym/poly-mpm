"""Analytical equilibrium, transport and finite-strain verification.

All dynamic refinements keep final time and the other discretisation fixed.
Fixed-cloud remaps and prescribed-gradient patches deliberately bypass dynamics;
they diagnose individual operators rather than validate a full simulation.
"""

import argparse
import hashlib
import json
import platform
from pathlib import Path
from types import MethodType, SimpleNamespace as NS
import sys

import numpy as np
import scipy

from diagnostics import inventory
from initialise import prepare
from interpolation import moment
from MPM import run
from transfers import scatter

CONFIGURATIONS = (("linear", "flip"), ("quadratic", "flip"), ("quadratic", "pic"), ("quadratic", "apic"))
RHO, BULK, SHEAR = 1000.0, 1e4, 4e3


def strip(
    cells,
    steps,
    final_time,
    *,
    basis="quadratic",
    transfer="apic",
    ppc=2,
    height=0.2,
    vertical=False,
    square=False,
    quadrature="point",
    gravity=0.0,
    initial_stress=None,
    speed=0.0,
    prepared=True,
):
    """Uniform midpoint quadrature, periodic x and slip/impermeable y walls.

    prepared=False returns the parameters before initialise.prepare."""
    nx, ny = (4, cells) if vertical else (cells, 2)
    height = 1.0 if vertical or square else height
    if square:
        nx = ny = cells
    x = (np.arange(nx * ppc) + 0.5) / (nx * ppc)
    y = height * (np.arange(ny * ppc) + 0.5) / (ny * ppc)
    X, Y = np.meshgrid(x, y)
    S = NS(
        law="elastic",
        n=X.size,
        X=X.ravel(),
        Y=Y.ravel(),
        A=height / X.size,
        rho=RHO,
        K=BULK,
        G=SHEAR,
        initial_v=[speed, 0.0, 0.0],
        domain_size=[1 / (nx * ppc), height / (ny * ppc)],
    )
    P = NS(
        dt=final_time / steps,
        t_f=final_time,
        savetime=1.0,
        S=S,
        gravity=gravity,
        initial_stress=initial_stress,
        interpolation=basis,
        transfer=transfer,
        quadrature=quadrature,
        cfl=0.5,
        timestep_policy="enforce",
        G=NS(nx=nx + 1, ny=ny + 1, x_min=0.0, x_max=1.0, y_min=0.0, y_max=height),
        B=NS(periodic=True, walls={"bottom"} if vertical else {"bottom", "top"}),
        O=NS(),
    )
    return prepare(P) if prepared else P


def _fields(b):
    return {name: getattr(b, name).copy() for name in ("x", "v", "F", "stress")}


def _finish(P, L, metrics, return_fields):
    b = L.materials[0]
    row = L.diagnostics.last
    metrics["quadrature"] = P.quadrature
    metrics.update(
        steps=P.tstep,
        dt=P.dt,
        final_time=P.t,
        minimum_J=float(b.J.min()),
        maximum_J=float(b.J.max()),
        relative_mass_error=float(row["mass_error"] / b.m.sum()),
        maximum_momentum_residual=float(np.max(np.abs(row["momentum_residual"]))),
        incomplete_support=int(row["particles_with_incomplete_support"]),
        minimum_stable_dt_ratio=min(min(r["timestep_limits"].values()) / P.dt for r in L.diagnostics.history),
        all_steps_within_stability_estimate=all(r["timestep_within_estimate"] for r in L.diagnostics.history),
    )
    return (metrics, _fields(b)) if return_fields else metrics


def hydrostatic(
    cells=16,
    steps=128,
    *,
    basis="quadratic",
    transfer="apic",
    ppc=2,
    final_time=0.05,
    return_fields=False,
    quadrature="point",
):
    """Exact rest: sigma=-rho*g*(H-y) I, g=1, traction-free top, supported bottom."""
    P, G, L = strip(
        cells,
        steps,
        final_time,
        basis=basis,
        transfer=transfer,
        ppc=ppc,
        quadrature=quadrature,
        vertical=True,
        gravity=-1.0,
        initial_stress="lithostatic",
    )
    b = L.materials[0]
    expected_stress = b.stress.copy()
    # Measure the initial discrete force defect before any time integration.
    G.wipe(P)
    L.get_nodal_mass_momentum(P, G)
    L.get_nodal_forces(P, G)
    residual = G.fe - G.fi
    free = G.valid_nodes & ~G.boundary_h
    free[P.G.nx - 1 :: P.G.nx] = False
    acceleration = residual[free, 1] / G.m[free]
    force_rms = float(np.sqrt(np.sum(G.m[free] * acceleration**2) / G.m[free].sum()))
    reaction = float(np.sum((G.fi - G.fe)[: P.G.nx - 1, 1]))
    run(P, G, L, output=False)
    return _finish(
        P,
        L,
        dict(
            cells=cells,
            quadrature=quadrature,
            particles_per_cell_axis=ppc,
            basis=basis,
            transfer=transfer,
            initial_force_rms_over_gravity=force_rms,
            initial_bottom_reaction_relative_error=abs(reaction / b.m.sum() - 1),
            velocity_rms=float(np.sqrt(np.mean(np.sum(b.v**2, axis=1)))),
            displacement_rms=float(np.sqrt(np.mean(np.sum((b.x - b.x0) ** 2, axis=1)))),
            stress_relative_l2=float(
                np.linalg.norm(b.stress - expected_stress) / np.linalg.norm(expected_stress)
            ),
        ),
        return_fields,
    )


def airy_stress(xy, amplitude=20.0):
    """Airy potential A sin(2*pi*x) sin(pi*y)^2; divergence zero, free y tractions."""
    x, y = xy[:, 0], xy[:, 1]
    # A=amplitude/(4*pi^2); calculate the three derivatives analytically.
    stress = np.zeros((len(x), 3, 3))
    stress[:, 0, 0] = 0.5 * amplitude * np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y)
    stress[:, 1, 1] = -amplitude * np.sin(2 * np.pi * x) * np.sin(np.pi * y) ** 2
    stress[:, 0, 1] = stress[:, 1, 0] = -0.5 * amplitude * np.cos(2 * np.pi * x) * np.sin(2 * np.pi * y)
    return stress


def translated_stress(
    cells=16,
    steps=128,
    *,
    basis="quadratic",
    transfer="apic",
    ppc=2,
    speed=1.0,
    final_time=0.1,
    return_fields=False,
    quadrature="point",
):
    """Self-equilibrated Airy stress advected horizontally on a periodic square."""
    P, G, L = strip(
        cells,
        steps,
        final_time,
        basis=basis,
        transfer=transfer,
        ppc=ppc,
        square=True,
        speed=speed,
        quadrature=quadrature,
    )
    b = L.materials[0]
    reference = airy_stress(b.x0)
    b.stress[:] = reference
    b.reference_stress[:] = reference
    # Airy tractions vanish on both y boundaries; slip walls are compatible.
    run(P, G, L, output=False)
    displacement = b.x - b.x0 - [speed * P.t, 0, 0]
    displacement[:, 0] = (displacement[:, 0] + 0.5) % 1.0 - 0.5
    return _finish(
        P,
        L,
        dict(
            cells=cells,
            quadrature=quadrature,
            particles_per_cell_axis=ppc,
            basis=basis,
            transfer=transfer,
            speed=speed,
            grid_cells_crossed=abs(speed * P.t) * cells,
            velocity_rms=float(np.sqrt(np.mean(np.sum((b.v - [speed, 0, 0]) ** 2, axis=1)))),
            displacement_rms=float(np.sqrt(np.mean(np.sum(displacement**2, axis=1)))),
            deformation_rms=float(np.linalg.norm(b.F - np.eye(3)) / np.sqrt(len(b))),
            stress_relative_l2=float(np.linalg.norm(b.stress - reference) / np.linalg.norm(reference)),
        ),
        return_fields,
    )


def manufactured_motion(X, time, amplitude=0.4):
    """Periodic 1D finite deformation embedded in 2D plane strain.

    x=X+a(t) sin(2*pi*X)/(2*pi), a=.4 sin(2*pi*t), y=Y.
    Reference balance: rho0*x_tt = dP11/dX + rho0*body.
    P11=(2G/3)(f^(1/3)-f^(-5/3))+K*log(f)/f, f=dx/dX.
    """
    sine, cosine = np.sin(2 * np.pi * X), np.cos(2 * np.pi * X)
    a = amplitude * np.sin(2 * np.pi * time)
    adot = 2 * np.pi * amplitude * np.cos(2 * np.pi * time)
    addot = -((2 * np.pi) ** 2) * a
    f = 1 + a * cosine
    position = X + a * sine / (2 * np.pi)
    velocity = adot * sine / (2 * np.pi)
    tangent = 2 * SHEAR / 9 * (f ** (-2 / 3) + 5 * f ** (-8 / 3)) + BULK * (1 - np.log(f)) / f**2
    force = addot * sine / (2 * np.pi) - tangent * (-2 * np.pi * a * sine) / RHO
    F = np.broadcast_to(np.eye(3), (len(X), 3, 3)).copy()
    F[:, 0, 0] = f
    stress, energy = elastic_reference(F)
    return dict(
        x=position, v=velocity, F=F, body=force, stress=stress, energy_density=energy, C=adot * cosine / f
    )


def elastic_reference(F, bulk=BULK, shear=SHEAR):
    """Closed-form continuum Cauchy stress and reference energy, no solver kernel."""
    J = np.linalg.det(F)
    bbar = J[:, None, None] ** (-2 / 3) * (F @ F.swapaxes(1, 2))
    trace = np.trace(bbar, axis1=1, axis2=2)
    stress = (
        shear * (bbar - trace[:, None, None] * np.eye(3) / 3) + bulk * np.log(J)[:, None, None] * np.eye(3)
    ) / J[:, None, None]
    energy = 0.5 * shear * (trace - 3) + 0.5 * bulk * np.log(J) ** 2
    return stress, energy


def large_deformation(
    cells=16,
    steps=128,
    *,
    basis="quadratic",
    transfer="apic",
    ppc=2,
    final_time=0.25,
    return_fields=False,
    quadrature="point",
):
    """Full timestep with analytically manufactured material-coordinate body load."""
    P, G, L = strip(cells, steps, final_time, basis=basis, transfer=transfer, ppc=ppc, quadrature=quadrature)
    b = L.materials[0]
    initial = manufactured_motion(b.x0[:, 0], 0.0)
    b.v[:, 0] = initial["v"]
    C = np.zeros_like(b.F)
    C[:, 0, 0] = initial["C"]
    b.affine_momentum = C @ moment(b.N, b.offsets)
    original = b.get_nodal_forces

    def manufactured_force(self, P, G):
        original(P, G)
        body = manufactured_motion(self.x0[:, 0], P.t)["body"]
        scatter(G.fe[:, 0], self.n, self.N * self.m[:, None] * body[:, None])

    b.get_nodal_forces = MethodType(manufactured_force, b)
    run(P, G, L, output=False)
    exact = manufactured_motion(b.x0[:, 0], P.t)
    dx = (b.x[:, 0] - exact["x"] + 0.5) % 1.0 - 0.5
    return _finish(
        P,
        L,
        dict(
            cells=cells,
            quadrature=quadrature,
            particles_per_cell_axis=ppc,
            basis=basis,
            transfer=transfer,
            displacement_relative_l2=float(np.linalg.norm(dx) / np.linalg.norm(exact["x"] - b.x0[:, 0])),
            velocity_rms=float(np.sqrt(np.mean((b.v[:, 0] - exact["v"]) ** 2 + b.v[:, 1] ** 2))),
            deformation_relative_l2=float(np.linalg.norm(b.F - exact["F"]) / np.linalg.norm(exact["F"])),
            stress_relative_l2=float(
                np.linalg.norm(b.stress - exact["stress"]) / np.linalg.norm(exact["stress"])
            ),
            density_relative_l2=float(
                np.linalg.norm(b.rho - RHO / exact["F"][:, 0, 0]) / np.linalg.norm(RHO / exact["F"][:, 0, 0])
            ),
            elastic_energy_relative_error=float(
                abs(b.elastic_energy.sum() - np.sum(b.V0 * exact["energy_density"]))
                / np.sum(b.V0 * exact["energy_density"])
            ),
        ),
        return_fields,
    )


def remap(cells=16, remaps=8, *, basis="quadratic", transfer="apic"):
    """Repeated P2G/G2P with fixed particles, no forces, motion or constitutive update."""
    P, G, L = strip(cells, 1, 0.01, basis=basis, transfer=transfer)
    b = L.materials[0]
    b.v[:, 0] = np.sin(2 * np.pi * b.x[:, 0])
    C = np.zeros_like(b.F)
    C[:, 0, 0] = 2 * np.pi * np.cos(2 * np.pi * b.x[:, 0])
    b.affine_momentum = C @ moment(b.N, b.offsets)
    initial = inventory(P, L)
    before = b.v.copy()
    positions = b.x.copy()
    for _ in range(remaps):
        G.wipe(P)
        L.get_reference_nodes(P, G)
        L.get_basis_functions(P, G)
        L.get_nodal_mass_momentum(P, G)
        L.move_material_points(P, G)
        b.x[:] = positions
    final = inventory(P, L)
    E0 = initial["kinetic"] + initial["affine_kinetic"]
    return dict(
        cells=cells,
        remaps=remaps,
        basis=basis,
        transfer=transfer,
        velocity_relative_l2=float(np.linalg.norm(b.v - before) / np.linalg.norm(before)),
        relative_energy_change=float((final["kinetic"] + final["affine_kinetic"] - E0) / E0),
        momentum_error=float(np.linalg.norm(np.asarray(final["momentum"]) - initial["momentum"])),
        maximum_position_change=float(np.max(np.abs(b.x - positions))),
    )


def affine_patch(steps=32, cells=8):
    """Prescribed noncommuting compression/shear/rotation; no force balance claimed."""
    P, G, L = strip(cells, steps, 1.0, vertical=True)
    b = L.materials[0]
    # Keep a small patch away from physical/periodic boundaries throughout motion.
    b.x[:, :2] = 0.5 + 0.1 * (b.x[:, :2] - [0.5, 0.1])
    b.x0[:] = b.x
    a, beta, gamma, omega = np.log(0.45), np.log(1.4), 1.2, 1.2

    def path(t):
        c, s = np.cos(omega * t), np.sin(omega * t)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
        A = np.array(
            [[np.exp(a * t), gamma * t * np.exp(beta * t), 0], [0, np.exp(beta * t), 0], [0, 0, 1.0]]
        )
        base = np.array([[a, gamma * (1 + (beta - a) * t), 0], [0, beta, 0], [0, 0, 0.0]])
        spin = np.array([[0, -omega, 0], [omega, 0, 0], [0, 0, 0.0]])
        return R @ A, spin + R @ base @ R.T

    for n in range(steps):
        L.get_reference_nodes(P, G)
        L.get_basis_functions(P, G)
        _, gradient = path((n + 0.5) * P.dt)
        xyz = np.column_stack((G.X, G.Y, np.zeros_like(G.X))) - [0.5, 0.5, 0]
        G.m[:] = 1.0
        G.valid_nodes[:] = True
        G.q[:] = xyz @ gradient.T
        G.q_reference[:] = G.q
        b.update_stress_strain(P, G)
        b.move_material_points(P, G)
    exact, _ = path(1.0)
    expected_stress, expected_energy = elastic_reference(exact[None])
    expected_x = (b.x0 - [0.5, 0.5, 0]) @ exact.T + [0.5, 0.5, 0]
    # Constitutive work integrated independently by the exact end-point potential.
    return dict(
        cells=cells,
        steps=steps,
        dt=P.dt,
        deformation_error=float(np.linalg.norm(b.F - exact) / np.sqrt(len(b))),
        position_error=float(np.linalg.norm(b.x - expected_x) / np.sqrt(len(b))),
        stress_relative_error=float(
            np.linalg.norm(b.stress - expected_stress) / (np.sqrt(len(b)) * np.linalg.norm(expected_stress))
        ),
        energy_relative_error=float(
            abs(b.elastic_energy.sum() - b.V0.sum() * expected_energy[0]) / (b.V0.sum() * expected_energy[0])
        ),
        determinant_error=float(np.max(np.abs(b.J - 0.63))),
        density_error=float(np.max(np.abs(b.rho - RHO / 0.63))),
    )


def force_balance(
    case, cells=16, *, ppc=2, basis="quadratic", quadrature="midpoint", integration="point", offset=0.0
):
    """Frozen exact states: spatial/particle quadrature residual, without remapping/time integration."""
    if case not in ("hydrostatic", "translated_stress", "large_deformation"):
        raise ValueError(f"Unknown analytical case: {case}")
    if quadrature == "gauss" and integration != "point":
        raise ValueError("Gaussian control uses point integration, not finite particle domains")
    hydro = case == "hydrostatic"
    P = strip(
        cells,
        256,
        0.25,
        basis=basis,
        transfer="flip",
        ppc=ppc,
        quadrature=integration,
        vertical=hydro,
        square=case == "translated_stress",
        gravity=-1.0 if hydro else 0.0,
        initial_stress="lithostatic" if hydro else None,
        prepared=False,
    )
    if quadrature not in ("midpoint", "gauss"):
        raise ValueError(f"Unknown quadrature: {quadrature}")
    if quadrature == "gauss":
        # Independent integration control: align Gaussian subintervals with all
        # half-cell spline knots in CURRENT coordinates, including boundaries.
        def axis(nodes, extent):
            z, w = np.polynomial.legendre.leggauss(8)
            count = 2 * (nodes - 1)
            return ((np.arange(count)[:, None] + 0.5 + 0.5 * z) / count * extent).ravel(), np.tile(
                w * extent / (2 * count), count
            )

        gx, wx = axis(P.G.nx, 1.0)
        gy, wy = axis(P.G.ny, P.G.y_max)
        xx, yy = np.meshgrid(gx, gy)
        volume = (wy[:, None] * wx).ravel()
        reference = xx.ravel().copy()
        if case == "large_deformation":
            for _ in range(12):
                reference -= (reference + 0.4 * np.sin(2 * np.pi * reference) / (2 * np.pi) - xx.ravel()) / (
                    1 + 0.4 * np.cos(2 * np.pi * reference)
                )
            volume /= 1 + 0.4 * np.cos(2 * np.pi * reference)
        P.S.n = xx.size
        P.S.X, P.S.Y, P.S.A = xx.ravel(), yy.ravel(), volume
    P, G, L = prepare(P)
    if quadrature == "gauss" and case == "large_deformation":
        L.materials[0].x0[:, 0] = reference
    b = L.materials[0]
    acceleration = np.zeros_like(b.v)
    body = np.zeros_like(b.v)
    if hydro:
        body[:, 1] = -1.0
    elif case == "translated_stress":
        if quadrature == "gauss":
            b.x0[:, 0] = (b.x[:, 0] - 0.1) % 1.0
        else:
            b.x[:, 0] = (b.x0[:, 0] + 0.1) % 1.0
        b.stress[:] = airy_stress(b.x0)
    else:
        exact = manufactured_motion(b.x0[:, 0], 0.25)
        b.x[:, 0] = exact["x"]
        b.F[:] = exact["F"]
        b.J[:] = np.linalg.det(b.F)
        b.V[:] = b.V0 * b.J
        b.rho[:] = b.m / b.V
        b.stress[:] = exact["stress"]
        body[:, 0] = exact["body"]
        acceleration[:, 0] = -2 * np.pi * 0.4 * np.sin(2 * np.pi * b.x0[:, 0])
    b.x[:, 0] += offset
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    G.wipe(P)
    L.get_nodal_mass_momentum(P, G)
    L.get_nodal_forces(P, G)
    # Assemble analytic body/inertia separately, with no discrete-force cancellation.
    load = np.zeros_like(G.fe)
    scatter(load, b.n, b.N[:, :, None] * b.m[:, None, None] * (body - acceleration)[:, None, :])
    if P.B.cyclic_lr:
        # load is not a named grid field: merge unique periodic endpoint values once.
        merged = load[:: P.G.nx] + load[P.G.nx - 1 :: P.G.nx]
        load[:: P.G.nx] = merged
        load[P.G.nx - 1 :: P.G.nx] = merged
    residual = load - G.fi
    residual[G.boundary_h, 1] = 0.0
    residual[G.boundary_v, 0] = 0.0
    active = G.valid_nodes.copy()
    active[P.G.nx - 1 :: P.G.nx] = False
    nodal_error = residual[active] / G.m[active, None]
    return dict(
        case=case,
        integration=integration,
        offset=offset,
        quadrature=quadrature,
        gaussian_points_per_half_cell=8 if quadrature == "gauss" else None,
        cells=cells,
        basis=basis,
        particles_per_cell_axis=ppc if quadrature == "midpoint" else None,
        acceleration_rms=float(np.sqrt(np.sum(G.m[active, None] * nodal_error**2) / G.m[active].sum())),
        acceleration_max=float(np.linalg.norm(nodal_error, axis=1).max()),
        force_resultant_error=float(np.linalg.norm(residual[active].sum(axis=0))),
    )


def galilean_pair(cells=16, steps=128):
    """Compare the same material labels with and without a unit horizontal boost."""
    moving, boosted = translated_stress(cells, steps, return_fields=True)
    stationary, rest = translated_stress(cells, steps, speed=0.0, return_fields=True)
    position = boosted["x"] - rest["x"] - [moving["final_time"], 0, 0]
    position[:, 0] = (position[:, 0] + 0.5) % 1.0 - 0.5
    velocity = boosted["v"] - rest["v"] - [1.0, 0, 0]
    return dict(
        cells=cells,
        steps=steps,
        grid_cells_crossed=moving["grid_cells_crossed"],
        position_rms=float(np.linalg.norm(position) / np.sqrt(len(position))),
        velocity_rms=float(np.linalg.norm(velocity) / np.sqrt(len(velocity))),
        stress_relative_l2=float(
            np.linalg.norm(boosted["stress"] - rest["stress"]) / np.linalg.norm(airy_stress(rest["x"]))
        ),
        deformation_rms=float(np.linalg.norm(boosted["F"] - rest["F"]) / np.sqrt(len(velocity))),
    )


def temporal_study(case, cells, steps=(32, 64, 128, 256), **kwargs):
    results = [case(cells, n, return_fields=True, **kwargs) for n in steps]
    finest = results[-1][1]
    for row, fields in results:
        differences = {name: fields[name] - finest[name] for name in fields}
        differences["x"][:, 0] = (differences["x"][:, 0] + 0.5) % 1.0 - 0.5
        row["difference_to_finest"] = {
            name: float(np.linalg.norm(value) / np.sqrt(len(value))) for name, value in differences.items()
        }
        row["reference_steps"] = steps[-1]
    return [row for row, _ in results]


def report():
    root = Path(__file__).resolve().parents[1]
    source_names = (
        "MPM.py",
        "particles.py",
        "constit.py",
        "grid.py",
        "interpolation.py",
        "kinematics.py",
        "transfers.py",
        "initialise.py",
        "stability.py",
        "diagnostics.py",
        "benchmarks/broader_verification.py",
    )
    return dict(
        description=__doc__,
        environment=dict(
            python=sys.version, numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform()
        ),
        source_sha256={name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in source_names},
        controls=dict(
            spatial="fixed dt and remap count; grid and midpoint particle quadrature refined together",
            temporal="fixed grid, particle cloud, final time and transfer; dt and remap count vary",
            quadrature="fixed grid, dt and final time; midpoint particles per cell varied",
            transfer="fixed grid/cloud; repeated remaps without motion or force updates",
            affine_patch="prescribed nodal velocity; tests finite-strain kinematics/constitutive integration, not dynamics",
        ),
        hydrostatic=dict(
            spatial=[hydrostatic(n, 256) for n in (8, 16, 32)],
            temporal=temporal_study(hydrostatic, 16),
            quadrature=[hydrostatic(16, 128, ppc=n) for n in (2, 4, 8)],
        ),
        translated_stress=dict(
            spatial=[translated_stress(n, 256) for n in (4, 8, 16, 32)],
            temporal=temporal_study(translated_stress, 8),
            transfers=[translated_stress(8, 128, basis=b, transfer=t) for b, t in CONFIGURATIONS],
            stationary=[translated_stress(n, 256, speed=0.0) for n in (4, 8, 16, 32)],
            galilean=[galilean_pair(n, 128) for n in (8, 16)],
        ),
        large_deformation=dict(
            spatial=[large_deformation(n, 256) for n in (8, 16, 32)],
            temporal=temporal_study(large_deformation, 16, steps=(128, 256, 512, 1024)),
            quadrature=[large_deformation(16, 256, ppc=n) for n in (2, 4, 8)],
            transfers=[large_deformation(16, 256, basis=b, transfer=t) for b, t in CONFIGURATIONS],
        ),
        remapping=[
            remap(n, k, basis=b, transfer=t)
            for b, t in CONFIGURATIONS
            for n in (8, 16, 32)
            for k in (1, 8, 32)
        ],
        affine_patch=[affine_patch(n) for n in (8, 16, 32, 64)],
        frozen_force_balance={
            case: dict(
                spatial=[force_balance(case, n) for n in (8, 16, 32)],
                quadrature=[force_balance(case, 16, ppc=n) for n in (2, 4, 8)],
                resolved_quadrature=[force_balance(case, n, quadrature="gauss") for n in (8, 16, 32)],
            )
            for case in ("hydrostatic", "translated_stress", "large_deformation")
        },
    )


def figure(results, path):
    """Scientific comparison plot; all plotted levels are retained in JSON."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for case, studies in results["frozen_force_balance"].items():
        rows = studies["spatial"]
        axes[0, 0].loglog(
            [r["cells"] for r in rows],
            [r["acceleration_rms"] for r in rows],
            "o-",
            label=case.replace("_", " "),
        )
    axes[0, 0].set(
        title="Frozen exact state: midpoint quadrature",
        xlabel="Cells along refined direction",
        ylabel="Force defect / nodal mass (RMS m/s²)",
    )
    rows = results["affine_patch"]
    for key, label in (("deformation_error", "Deformation gradient"), ("position_error", "Position (m)")):
        axes[0, 1].loglog([r["dt"] for r in rows], [r[key] for r in rows], "o-", label=label)
    axes[0, 1].set(
        title="Prescribed affine patch: time integration", xlabel="Timestep (s)", ylabel="RMS error"
    )
    for transfer in ("pic", "apic", "flip"):
        rows = [
            r
            for r in results["remapping"]
            if r["cells"] == 16 and r["basis"] == "quadratic" and r["transfer"] == transfer
        ]
        axes[1, 0].semilogx(
            [r["remaps"] for r in rows],
            [-r["relative_energy_change"] for r in rows],
            "o-",
            label=transfer.upper(),
        )
    axes[1, 0].set(
        title="Fixed cloud: transfer filtering only",
        xlabel="Number of remaps",
        ylabel="Fraction of represented kinetic energy lost",
    )
    for key, label in (("spatial", "Grid study (dt fixed)"), ("temporal", "Time study (grid fixed)")):
        rows = results["large_deformation"][key]
        axes[1, 1].loglog(
            [r["cells"] if key == "spatial" else r["steps"] for r in rows],
            [r["stress_relative_l2"] for r in rows],
            "o-",
            label=label,
        )
    axes[1, 1].set(
        title="Finite-strain dynamics: non-monotone levels",
        xlabel="Grid cells / number of timesteps",
        ylabel="Relative stress L2 error",
    )
    for ax in axes.ravel():
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", default="/tmp/poly-mpm-broader-verification.json")
    parser.add_argument("--figure", help="optional PNG comparison figure")
    args = parser.parse_args()
    results = report()
    Path(args.output).write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    if args.figure:
        figure(results, args.figure)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
