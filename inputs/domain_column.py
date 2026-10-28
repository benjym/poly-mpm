"""Hydrostatic column with finite-domain quadratic/APIC integration."""

from types import SimpleNamespace as NS
import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 2.5e-4
        self.t_f = 0.05
        self.savetime = 0.01
        self.max_g = -1.0
        self.pressure = "lithostatic"
        self.interpolation = "quadratic"
        self.quadrature = "domain"
        self.transfer = "apic"
        self.CFL = 0.5
        self.timestep_policy = "enforce"
        self.G = NS(nx=5, ny=17, x_m=0.0, x_M=1.0, y_m=0.0, y_M=1.0)
        width = np.array([1 / 8, 1 / 32])
        X, Y = np.meshgrid((np.arange(8) + 0.5) * width[0], (np.arange(32) + 0.5) * width[1])
        self.S = NS(
            law="elastic",
            n=X.size,
            X=X.ravel(),
            Y=Y.ravel(),
            A=np.prod(width),
            domain_size=width,
            rho=1000.0,
            K=1e4,
            G=4e3,
        )
        self.B = NS(cyclic_lr=True, has_bottom=True)
        self.O = NS(plot_continuum=True)
