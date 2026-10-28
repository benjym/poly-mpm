"""Copy to inputs/my_sim.py, then run python MPM.py my_sim from the checkout."""
from types import SimpleNamespace as NS
import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4
        self.t_f = .01
        self.savetime = .005
        self.max_g = -9.81
        self.save_dir = 'im/my_sim/'
        self.G = NS(nx=11, ny=11, x_m=0., x_M=1., y_m=0., y_M=1.)
        self.B = NS(has_bottom=True)
        self.O = NS(plot_material_points=True, plot_continuum=True, save_phi_MP=True)

        # Two-by-two particles in each occupied grid cell: a 0.2 m square block.
        coordinates = np.arange(.325, .5, .05)
        x, y = np.meshgrid(coordinates, coordinates)
        E, nu = 1e5, .3
        self.S = NS(law='elastic', n=x.size, X=x.ravel(), Y=y.ravel(),
                    rho=1000., A=.05**2, initial_v=[0., 0., 0.],
                    K=E/(3*(1-2*nu)), G=E/(2*(1+nu)))
