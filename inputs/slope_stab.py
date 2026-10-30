"""Stability of a vertical cut in Drucker-Prager soil as gravity is applied."""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s)
        self.t_f = 5.0  # (s)
        self.savetime = 0.05  # (s)
        self.gravity = -9.81  # (m/s^2)
        self.force_scale = 0.7  # local non-viscous damping
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G)]

    def gravity_ramp(self, t):
        """Gravity rises smoothly from zero over about 0.5 s."""
        return 1.0 - np.exp(-3.0 * t**2 / 0.5**2)


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 12.0  # (m)
        self.y_min, self.y_max = 0.0, 6.0  # (m)
        scale = 5  # grid nodes per metre
        self.nx = int(self.x_max - self.x_min) * scale + 1
        self.ny = int(self.y_max - self.y_min) * scale + 1


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom", "left", "right"}
        self.traction_vertical = 0.0  # (Pa)


class Solid_Params:
    def __init__(self, G):
        self.law = "dp"
        self.rho = 2650.0  # (kg/m^3)
        E, nu = 1e6, 0.2  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.G = E / (2.0 * (1.0 + nu))  # shear modulus (Pa)
        self.mu_s = 1.0
        self.beta = 1.0 * self.mu_s

        # Soil filling the domain except the upper left block above a cut of slope_angle at corner.
        pts_per_cell, slope_angle, corner = 3, 90.0, (3.0, 3.0)
        dx = (G.x_max - G.x_min) / (G.nx - 1)
        dy = (G.y_max - G.y_min) / (G.ny - 1)
        nx = (G.nx - 1) * pts_per_cell
        ny = (G.ny - 1) * pts_per_cell
        gap = np.array((dx, dy)) / (2 * pts_per_cell)
        m = np.tan(slope_angle * np.pi / 180.0)  # slope gradient
        xp = np.linspace(G.x_min + gap[0], G.x_max - gap[0], nx)
        yp = np.linspace(G.y_min + gap[1], G.y_max - gap[1], ny)
        X = np.tile(xp, ny)
        Y = np.repeat(yp, nx)
        right, below = X > corner[0], Y < corner[1]
        keep = (right | below) & ~(right & (Y > m * (X - corner[0]) + corner[1]))
        self.X, self.Y = X[keep], Y[keep]
        self.n = len(self.X)
        total_area = (G.x_max - G.x_min) * (G.y_max - G.y_min)
        missing_corner = (corner[0] - G.x_min) * (G.y_max - corner[1])
        self.A = (total_area - missing_corner) / self.n  # (m^2)


class Output_Params:
    def __init__(self):
        self.figures = {"continuum"}
        self.continuum_fig_size = [10, 6]
