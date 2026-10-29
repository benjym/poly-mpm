"""F-bar volumetric averaging of the velocity gradient."""

import contextlib
import io
from types import SimpleNamespace as NS

import numpy as np
import pytest

from initialise import prepare


def block(**options):
    P = NS(
        dt=1e-5,
        t_f=1e-3,
        savetime=1e-3,
        max_g=0.0,
        save_dir="im/fbar_test/",
        G=NS(nx=9, ny=9, x_m=0.0, x_M=1.0, y_m=0.0, y_M=1.0),
        B=NS(),
        O=NS(),
        fbar=True,
    )
    x, y = np.meshgrid(np.arange(0.28125, 0.72, 0.0625), np.arange(0.28125, 0.72, 0.0625))
    P.S = NS(law="elastic", n=x.size, X=x.ravel(), Y=y.ravel(), rho=1000.0, A=0.0625**2, K=1e7, G=1e6)
    for name, value in options.items():
        setattr(P, name, value)
    with contextlib.redirect_stdout(io.StringIO()):
        return prepare(P, "fbar_test")


def gradient(P, G, L, velocity):
    """Velocity gradient from a prescribed nodal velocity field, with and without F-bar."""
    b = L.materials[0]
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    G.q = G.m[:, None] * velocity
    b.update_stress_strain(P, G)
    return b


def nodal(G, function):
    return np.column_stack((*function(G.X, G.Y), np.zeros_like(G.X)))


def test_affine_motion_is_unchanged():
    P, G, L = block()
    rate = np.array([[0.3, 0.1], [-0.2, 0.5]])
    b = gradient(
        P, G, L, nodal(G, lambda x, y: (rate[0, 0] * x + rate[0, 1] * y, rate[1, 0] * x + rate[1, 1] * y))
    )
    np.testing.assert_allclose(b.L[:, :2, :2], np.broadcast_to(rate, (len(b), 2, 2)), atol=1e-12)


def test_total_volume_rate_is_kept_and_only_the_trace_changes():
    P, G, L = block()
    field = nodal(G, lambda x, y: (np.sin(5 * x) * y, np.cos(3 * y) * x**2))
    reference = block(fbar=False)
    plain = gradient(*reference, field.copy())
    averaged = gradient(P, G, L, field.copy())
    trace = lambda b: b.L[:, 0, 0] + b.L[:, 1, 1]
    assert not np.allclose(trace(averaged), trace(plain))
    # Weighted by the volume at the start of the step (V0 here), as the projection is.
    np.testing.assert_allclose(
        np.sum(averaged.V0 * trace(averaged)), np.sum(plain.V0 * trace(plain)), rtol=1e-12
    )
    deviator = lambda b: b.L[:, :2, :2] - 0.5 * trace(b)[:, None, None] * np.eye(2)
    np.testing.assert_allclose(deviator(averaged), deviator(plain), atol=1e-12)


@pytest.mark.parametrize("setting", [dict(mechanics="legacy"), dict(solver="granular")])
def test_fbar_requires_explicit_updated_mechanics(setting):
    with pytest.raises(ValueError, match="fbar"):
        block(**setting)
