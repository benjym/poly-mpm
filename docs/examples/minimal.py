"""Copy to inputs/my_sim.py, then run python MPM.py my_sim from the checkout."""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s)
        self.t_f = 0.01  # (s)
        self.savetime = 0.005  # (s)
        self.gravity = -9.81  # (m/s^2)
        self.output_dir = "im/my_sim/"
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = Solid_Params()


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 1.0  # (m)
        self.y_min, self.y_max = 0.0, 1.0  # (m)
        self.nx = 11
        self.ny = 11


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom"}


class Output_Params:
    def __init__(self):
        self.figures = {"material_points", "continuum"}
        self.fields = {"phi_material_points"}


class Solid_Params:
    def __init__(self):
        self.law = "elastic"
        self.rho = 1000.0  # (kg/m^3)
        E, nu = 1e5, 0.3  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3 * (1 - 2 * nu))  # bulk modulus (Pa)
        self.G = E / (2 * (1 + nu))  # shear modulus (Pa)

        # Two-by-two particles in each occupied grid cell: a 0.2 m square block.
        coordinates = np.arange(0.325, 0.5, 0.05)
        x, y = np.meshgrid(coordinates, coordinates)
        self.X, self.Y, self.n = x.ravel(), y.ravel(), x.size
        self.A = 0.05**2  # (m^2)
        self.initial_v = [0.0, 0.0, 0.0]  # (m/s)
