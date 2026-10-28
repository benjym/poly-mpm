import os

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/poly-mpm-mpl")
from types import SimpleNamespace as NS
import numpy as np
import pytest
from initialise import prepare


def material(law="elastic", n=4, **overrides):
    data = dict(
        law=law,
        n=n,
        X=np.linspace(0.2, 0.4, n),
        Y=np.linspace(0.25, 0.45, n),
        rho=1000.0,
        A=0.01,
        K=100.0,
        G=40.0,
        k=500.0,
        s=2.5,
        mu_s=10.0,
        mu_v=5.0,
        rho_s=1666.0,
        mu_0=0.3,
        mu_1=0.7,
        delta_mu=0.4,
        I_0=0.279,
        eta_max=100.0,
        packing=0.6,
    )
    data.update(overrides)
    return NS(**data)


def state(law="elastic", ns=1, materials=None, **overrides):
    P = NS(
        dt=1e-5,
        t_f=0.01,
        savetime=0.01,
        max_g=0.0,
        G=NS(nx=5, ny=5, x_m=0.0, x_M=1.0, y_m=0.0, y_M=1.0, ns=ns, s=np.linspace(0.001, 0.003, ns)),
        B=NS(),
        O=NS(),
        S=materials if materials is not None else material(law),
    )
    P.__dict__.update(overrides)
    return prepare(P)


@pytest.fixture
def basic():
    return state()
