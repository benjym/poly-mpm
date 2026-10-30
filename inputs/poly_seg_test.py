"""Segregation of five grain sizes in a viscous layer on a periodic slope.

python MPM.py poly_seg_test <slope_degrees>
"""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-5  # (s)
        self.savetime = 1e2 * self.dt  # (s)
        self.t_f = 10.0  # (s)
        self.gravity = -9.81  # (m/s^2)
        self.gravity_angle = -np.radians(float(args[1]))  # slope angle (rad)
        self.initial_stress = "lithostatic"
        self.segregate = True
        self.segregation_rate = 1.0  # inter-particle drag coefficient
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G)]
        self.output_dir = f"im/seg/{len(self.G.grain_sizes)}/theta_{args[1]}/"


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 0.02  # (m)
        self.y_min, self.y_max = 0.0, 1.0  # (m)
        self.nx = 2
        self.ny = 201
        self.grain_sizes = np.array([0.2, 0.4, 0.6, 0.8, 1.0])  # (m)


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom"}
        self.wall_contact = "no_slip"
        self.periodic = True


class Solid_Params:
    def __init__(self, G):
        self.law = "viscous"
        self.rho = 2700.0  # (kg/m^3)
        self.phi = [0.2, 0.2, 0.2, 0.2, 0.2]
        self.eta = 1e3  # shear viscosity (Pa s)
        self.eta_v = 0.0  # bulk viscosity (Pa s)
        E, nu = 1e6, 0.499  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3 * (1 - 2 * nu))  # bulk modulus (Pa)

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


class Output_Params:
    def __init__(self):
        self.figures = {"continuum", "gsd_grid"}
        self.continuum_fig_size = [10, 4]
