"""Copy to inputs/my_sim.py, then run python MPM.py my_sim from the checkout."""

from types import SimpleNamespace as NS
import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4
        self.t_f = 0.01
        self.savetime = 0.005
        self.max_g = -9.81
        self.save_dir = "im/my_sim/"
        self.G = NS(nx=11, ny=11, x_m=0.0, x_M=1.0, y_m=0.0, y_M=1.0)
        self.B = NS(has_bottom=True)
        self.O = NS(plot_material_points=True, plot_continuum=True, save_phi_MP=True)

        # Two-by-two particles in each occupied grid cell: a 0.2 m square block.
        coordinates = np.arange(0.325, 0.5, 0.05)
        x, y = np.meshgrid(coordinates, coordinates)
        E, nu = 1e5, 0.3
        self.S = NS(
            law="elastic",
            n=x.size,
            X=x.ravel(),
            Y=y.ravel(),
            rho=1000.0,
            A=0.05**2,
            initial_v=[0.0, 0.0, 0.0],
            K=E / (3 * (1 - 2 * nu)),
            G=E / (2 * (1 + nu)),
        )
