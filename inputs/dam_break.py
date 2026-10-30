"""Collapse of a rate-dependent Drucker-Prager column into a channel with an outlet."""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s)
        self.savetime = 0.0001  # (s)
        self.t_f = 10.0  # (s)
        self.gravity = -9.81  # (m/s^2)
        self.output_dir = "im/dam_break_no_seg/"
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G)]

    def gravity_ramp(self, t):
        """Gravity rises smoothly from zero over about 0.05 s."""
        return 1.0 - np.exp(-3.0 * t**2 / 0.05**2)


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 2.0  # (m)
        self.y_min, self.y_max = 0.0, 1.0  # (m)
        self.ny = 21
        self.nx = (self.ny - 1) * 4 + 1
        self.grain_sizes = np.array([0.5, 1.0])  # (m)


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom", "right"}
        self.outlets = {"left"}


class Solid_Params:
    def __init__(self, G):
        self.law = "dp_rate"
        self.rho = 2650.0  # (kg/m^3)
        E, nu = 1e4, 0.2  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.G = E / (2.0 * (1.0 + nu))  # shear modulus (Pa)
        self.mu_s = 1.0
        self.beta = 1.0 * self.mu_s
        self.t_star = 0.01  # (s)

        # A column in the right quarter of the channel, one cell short of the top.
        pts_per_cell = 3
        dx = (G.x_max - G.x_min) / (G.nx - 1)
        dy = (G.y_max - G.y_min) / (G.ny - 1)
        nx = (G.nx - 1) // 4 * pts_per_cell
        ny = (G.ny - 2) * pts_per_cell
        gap = np.array((dx, dy)) / (2 * pts_per_cell)
        xp = np.linspace(G.x_max - (G.x_min + G.x_max) / 4 + gap[0], G.x_max - gap[0], nx)
        yp = np.linspace(G.y_min + gap[1], G.y_max - dy - gap[1], ny)
        self.X = np.tile(xp, ny)
        self.Y = np.repeat(yp, nx)
        self.n = len(self.X)
        self.A = (G.x_max - G.x_min) * (G.y_max - G.y_min - dy) / 4.0 / self.n  # (m^2)


class Output_Params:
    def __init__(self):
        self.figures = {"continuum", "material_points"}
        self.continuum_fig_size = [18, 6]
        self.mp_fig_size = [18, 4]
