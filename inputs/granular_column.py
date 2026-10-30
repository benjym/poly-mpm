"""Small experimental rigid-grain column; not a validated production model."""

from types import SimpleNamespace as NS
import numpy as np


class Params:
    def __init__(self, args):
        self.solver = "granular"
        # Coulomb friction with a hard packing limit; granular_solver defaults to "cidr".
        self.granular_solver = "velocity"
        self.interpolation = "linear"
        self.transfer = "flip"
        self.dt = 0.005
        self.t_f = 0.05
        self.savetime = 0.01
        self.gravity = -9.81
        self.G = NS(nx=5, ny=7, x_min=0.0, x_max=1.0, y_min=0.0, y_max=1.5)
        self.B = NS(walls={"bottom", "left", "right"})
        self.O = NS(figures={"continuum", "material_points"})
        spacing = 0.125
        x, y = np.meshgrid(np.arange(spacing / 2, 1.0, spacing), np.arange(spacing / 2, 1.0, spacing))
        self.S = NS(
            law="granular",
            n=x.size,
            X=x.ravel(),
            Y=y.ravel(),
            A=spacing**2,
            rho=1500.0,
            rho_s=2500.0,
            phi_max=0.6,
            mu_s=0.5,
        )
