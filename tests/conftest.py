import os

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/poly-mpm-mpl")
from types import SimpleNamespace as NS
import numpy as np
import pytest
from initialise import prepare
from settings import MATERIAL_COMMON, MATERIAL_FIELDS

# Default values of every law's fields; material() keeps those its law reads.
MATERIAL_DEFAULTS = dict(
    K=100.0,
    G=40.0,
    yield_stress=np.sqrt(3) * 500.0,
    eta=10.0,
    eta_v=5.0,
    rho_s=1666.0,
    mu_s=0.3,
    mu_d=0.7,
    I_0=0.279,
    eta_max=100.0,
    phi_flow=0.6,
)


def material(law="elastic", n=4, **overrides):
    data = dict(law=law, n=n, X=np.linspace(0.2, 0.4, n), Y=np.linspace(0.25, 0.45, n), rho=1000.0, A=0.01)
    known = MATERIAL_COMMON + MATERIAL_FIELDS.get(law, tuple(MATERIAL_DEFAULTS))
    data.update({name: value for name, value in MATERIAL_DEFAULTS.items() if name in known})
    data.update(overrides)
    return NS(**data)


def grid(nx=5, ny=5, ns=1, **overrides):
    """Grid settings of the unit square, with ns grain sizes."""
    data = dict(
        nx=nx, ny=ny, x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.0, grain_sizes=np.linspace(0.001, 0.003, ns)
    )
    data.update(overrides)
    return NS(**data)


def state(law="elastic", ns=1, materials=None, **overrides):
    P = NS(
        dt=1e-5,
        t_f=0.01,
        savetime=0.01,
        G=grid(ns=ns),
        B=NS(),
        O=NS(),
        S=materials if materials is not None else material(law),
    )
    P.__dict__.update(overrides)
    return prepare(P)


@pytest.fixture
def basic():
    return state()
