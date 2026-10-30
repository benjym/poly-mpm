"""One cell of rate-dependent Drucker-Prager material under boundary tractions.

Both tractions rise together to 10 Pa at t = 2.5 s; after that the horizontal traction keeps
rising while the vertical one falls, shearing the cell.
"""

import numpy as np
from plotting import PQHistory

_MAX_Q = 10.0  # traction at the end of loading (Pa)
_T_C = 2.5  # end of loading (s)


class Params:
    def __init__(self, args):
        self.dt = 1e-3  # (s)
        self.savetime = 0.1  # (s)
        self.t_f = 5.0  # (s)
        self.output_dir = "im/dp_rate_unit_test"
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.S = Solid_Params(self.G)
        self.O = Output_Params()


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 1.0  # (m)
        self.y_min, self.y_max = 0.0, 1.0  # (m)
        self.nx = 2
        self.ny = 2


class Boundary_Params:
    def traction_vertical(self, t):
        return _MAX_Q * t / _T_C if t < _T_C else 2 * _MAX_Q - t / _T_C * _MAX_Q

    def traction_horizontal(self, t):
        return _MAX_Q * t / _T_C


class Solid_Params:
    def __init__(self, G):
        self.law = "dp_rate"
        self.rho = 2650.0  # (kg/m^3)
        E, nu = 1e5, 0.2  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.G = E / (2.0 * (1.0 + nu))  # shear modulus (Pa)
        self.mu_s = 1.0
        self.beta = 1.0 * self.mu_s
        self.t_star = 1e-2  # (s)

        pts_per_cell = 3
        dx = (G.x_max - G.x_min) / (G.nx - 1)
        dy = (G.y_max - G.y_min) / (G.ny - 1)
        nx = (G.nx - 1) * pts_per_cell
        ny = (G.ny - 1) * pts_per_cell
        gap = np.array((dx, dy)) / (2 * pts_per_cell)
        xp = np.linspace(G.x_min + gap[0], G.x_max - gap[0], nx)
        yp = np.linspace(G.y_min + gap[1], G.y_max - gap[1], ny)
        self.X = np.tile(xp, ny)
        self.Y = np.repeat(yp, nx)
        self.n = len(self.X)
        self.A = (G.x_max - G.x_min) * (G.y_max - G.y_min) / self.n  # (m^2)


class Output_Params(PQHistory):
    pass
