"""A viscous layer between two rough walls, periodic in x, with grid segregation of one size."""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-8  # (s)
        self.savetime = 10 * self.dt  # (s)
        self.t_f = 100.0  # (s)
        self.gravity = -9.81  # (m/s^2)
        self.segregate = True
        self.segregation_rate = 1e-3  # inter-particle drag coefficient
        self.time_stepping = "dynamic"
        self.cfl = 0.2
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G)]
        self.output_dir = f"im/seg_unit_test/ny_{self.G.ny}/ns_1/c_{self.segregation_rate}/"


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 0.25  # (m)
        self.y_min, self.y_max = 0.0, 1.0  # (m)
        self.nx = 6
        self.ny = 21


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom", "top"}
        self.wall_contact = "no_slip"
        self.periodic = True


class Solid_Params:
    def __init__(self, G):
        self.law = "viscous"
        self.rho = 1000.0  # (kg/m^3)
        self.eta = 1e5  # shear viscosity (Pa s)
        self.eta_v = 1e5  # bulk viscosity (Pa s)
        self.K = 0.0  # bulk modulus (Pa)

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
        self.figures = {"continuum"}
        self.continuum_fig_size = [10, 6]
