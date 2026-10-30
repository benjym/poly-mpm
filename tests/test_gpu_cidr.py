"""Metal CIDR step (mpm_gpu/cidr.py) on real hardware. Opt in with POLY_MPM_TEST_METAL=1."""

import os
from copy import deepcopy
from types import SimpleNamespace as NS

import numpy as np
import pytest

from test_cidr import CIDR, cidr_case

pytestmark = pytest.mark.skipif(
    os.environ.get("POLY_MPM_TEST_METAL") != "1",
    reason="requires explicit Apple silicon Metal hardware validation",
)


def test_envelope_kernel_matches_the_cpu_envelope():
    import cidr
    from mpm_gpu import mechanics
    from mpm_gpu.cidr import MODULE, CidrSolver
    from mpm_gpu.runner import runtime

    rt = runtime()
    S = NS(rho_s=2500.0, **CIDR)
    solver = CidrSolver(
        rt, mechanics.device_grid(rt, 5, 5, 0.25, 0.25), S, dt=0.002, gamma=1e-8, tolerance=1e-4
    )
    model = cidr.CriticalStateCIDR(S.mu_s, S.mu_d, S.I_0, 2 * S.grain_diameter * np.sqrt(S.rho_s))
    rng = np.random.default_rng(0)
    n = 400
    for gamma in (1.3e-6, 1e-8):
        q = rng.standard_normal((n, 3)) * 10 ** rng.uniform(-4, 0, (n, 1))
        C = 4e4 * 10 ** rng.uniform(-0.5, 0.5, n)
        _, g_ref, H_ref, jam_ref = cidr.Problem(None, None, None, None, None, model, C).envelope(q, gamma)
        f32 = np.float32
        grad, hessian = rt.buffer(count=3 * n, dtype=f32), rt.buffer(count=9 * n, dtype=f32)
        jammed = rt.buffer(count=n, dtype=np.uint32)
        encoder = rt.device.create_command_encoder()
        rt.dispatch(
            encoder,
            "cidr_envelope_test",
            n,
            module=MODULE,
            count=n,
            r=float(1 / gamma),
            **solver.law,
            table=solver.table,
            q=rt.buffer(q, dtype=f32),
            modulus=rt.buffer(C, dtype=f32),
            grad=grad,
            hessian=hessian,
            jammed=jammed,
        )
        rt.submit(encoder)
        assert np.array_equal(rt.read(jammed, np.uint32).astype(bool), jam_ref)
        # Errors relative to the stress scale r |q|: rows in the zero-pressure branch have a
        # gradient of zero up to roundoff on both sides.
        stress = np.linalg.norm(q, axis=1) / gamma
        g = rt.read(grad, f32).reshape(n, 3)
        assert np.quantile(np.linalg.norm(g - g_ref, axis=1) / stress, 0.99) < 1e-4
        # The Hessian jumps where the cap or zero pressure becomes active, and rows within FP32
        # roundoff of those boundaries may take either side: a median and a looser tail.
        H = rt.read(hessian, f32).reshape(n, 3, 3)
        error = np.linalg.norm(H - H_ref, axis=(1, 2)) / np.linalg.norm(H_ref, axis=(1, 2)).clip(1e-300)
        assert np.median(error) < 1e-5
        assert np.quantile(error, 0.9) < 1e-2


def test_static_column_jams_without_creep_on_metal():
    from MPM import run

    P, G, L = cidr_case("granular_column", dt=0.02, t_f=0.2)
    P.backend = "metal"
    initial = L.materials[0].x.copy()
    run(P, G, L, output=False)
    b = L.materials[0]
    assert np.max(np.abs(b.v)) < 1e-6
    assert np.max(np.abs(b.x - initial)) < 1e-7
    history = [r["granular_solver"] for r in L.diagnostics.history]
    assert all(h["backend"] == "metal" and h["method"] == "cidr" for h in history)


def test_collapse_agrees_with_the_cpu():
    from MPM import run

    state = cidr_case("granular_collapse", dt=0.01, t_f=0.1)
    finals = {}
    for backend in ("cpu", "metal"):
        P, G, L = deepcopy(state)
        P.backend = backend
        run(P, G, L, output=False)
        finals[backend] = L.materials[0].x.copy()
    assert np.max(np.abs(finals["cpu"] - finals["metal"])) < 1e-6
    assert np.max(np.abs(finals["metal"][:, 0])) > np.max(np.abs(state[2].materials[0].x[:, 0])) + 0.01


def test_packing_limit_on_metal():
    from MPM import run

    P, G, L = cidr_case("granular_collapse", dt=0.05, t_f=0.3)
    P.backend = "metal"
    with pytest.raises(FloatingPointError, match="past phi_max"):
        run(P, G, L, output=False)
    P, G, L = cidr_case("granular_collapse", dt=0.02, t_f=0.3, time_stepping="dynamic", cfl=0.2)
    P.backend = "metal"
    run(P, G, L, output=False)
    history = [r["granular_solver"]["packing_closure"] for r in L.diagnostics.history]
    assert max(history[1:]) < 0.5


@pytest.mark.parametrize("module, gravity", [("granular_collapse", None), ("granular_column", 2.0)])
def test_one_sided_walls_agree_with_the_cpu(module, gravity):
    from initialise import prepare
    from MPM import run
    from test_cidr import CIDR

    finals = {}
    for backend in ("cpu", "metal"):
        Params = __import__(f"inputs.{module}", fromlist=["Params"]).Params
        P = Params([])
        P.granular_solver, P.dt, P.t_f, P.backend = "cidr", 0.01, 0.1, backend
        if gravity is not None:
            P.gravity = gravity  # upwards: the column leaves the floor
        P.B.one_sided = True
        for key, value in CIDR.items():
            setattr(P.S, key, value)
        P, G, L = prepare(P)
        start = L.materials[0].x.copy()
        run(P, G, L, output=False)
        finals[backend] = L.materials[0].x.copy()
    # FP32 differences grow through the re-solves (1.25e-1 m spacing).
    assert np.max(np.abs(finals["cpu"] - finals["metal"])) < 2e-5
    if gravity is not None:
        assert np.min(finals["metal"][:, 1] - start[:, 1]) > 0.5 * gravity * 0.1**2 * 0.9


@pytest.mark.parametrize("degrees", [15.0, 10.0])
def test_coulomb_floor_on_metal(degrees):
    from test_cidr import sliding_layer

    P, b = sliding_layer(degrees, backend="metal")
    theta = np.radians(degrees)
    expected = max(9.81 * (np.sin(theta) - 0.2 * np.cos(theta)), 0.0)
    assert np.allclose(-b.v[:, 0] / P.t, expected, atol=1e-4)
    assert np.max(np.abs(b.v[:, 1])) < 1e-6


def test_collapse_with_coulomb_walls_agrees_with_the_cpu():
    from initialise import prepare
    from MPM import run
    from test_cidr import CIDR

    finals = {}
    for backend in ("cpu", "metal"):
        from inputs.granular_collapse import Params

        P = Params([])
        P.granular_solver, P.dt, P.t_f, P.backend = "cidr", 0.01, 0.1, backend
        P.B.wall_contact, P.B.wall_friction = "coulomb", 0.3
        for key, value in CIDR.items():
            setattr(P.S, key, value)
        P, G, L = prepare(P)
        run(P, G, L, output=False)
        finals[backend] = L.materials[0].x.copy()
    assert np.max(np.abs(finals["cpu"] - finals["metal"])) < 2e-5
