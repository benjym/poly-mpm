import os
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/poly-mpm-mpl')
from types import SimpleNamespace as NS
import numpy as np
import pytest
from initialise import prepare


def material(law='elastic', n=4, **overrides):
    data = dict(law=law, n=n, X=np.linspace(.2, .4, n), Y=np.linspace(.25, .45, n),
                rho=1000., A=.01, K=100., G=40., k=500., s=2.5,
                mu_s=10., mu_v=5., rho_s=1666., mu_0=.3, mu_1=.7,
                delta_mu=.4, I_0=.279, eta_max=100., packing=.6)
    data.update(overrides)
    return NS(**data)


def state(law='elastic', ns=1, materials=None, **overrides):
    P = NS(dt=1e-5, t_f=.01, savetime=.01, max_g=0.,
           G=NS(nx=5, ny=5, x_m=0., x_M=1., y_m=0., y_M=1.,
                ns=ns, s=np.linspace(.001, .003, ns)),
           B=NS(), O=NS(), S=materials if materials is not None else material(law))
    P.__dict__.update(overrides)
    return prepare(P)


@pytest.fixture
def basic():
    return state()
