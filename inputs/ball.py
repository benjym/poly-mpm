"""An elastic ball dropped onto a floor."""

import numpy as np
from plotting import Plotting

plot = Plotting()


class Params:
    def __init__(self, args):
        self.dt = 1e-3  # (s)
        self.savetime = 0.1  # (s)
        self.t_f = 1.0  # (s)
        self.gravity = -10.0  # (m/s^2)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = Solid_Params()


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = -1.0, 1.0  # (m)
        self.y_min, self.y_max = 0.0, 2.0  # (m)
        self.nx = 6
        self.ny = 6


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom"}


class Solid_Params:
    def __init__(self):
        self.law = "elastic"
        self.rho = 1000.0  # (kg/m^3)
        E, nu = 1.0e5, 0.3  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.G = E / (2.0 * (1.0 + nu))  # shear modulus (Pa)

        # Rings of particles in a disc of radius r centred at (0, 1).
        rings, around, r, centre = 20, 50, 0.3, (0.0, 1.0)
        X, Y = [], []
        for radius in np.linspace(0, r, rings):
            count = int(np.around(around * radius / r))
            for angle in np.linspace(0, (1.0 - 1.0 / max(count, 1)) * 2 * np.pi, count):
                X.append(centre[0] + radius * np.sin(angle))
                Y.append(centre[1] + radius * np.cos(angle))
        self.X, self.Y, self.n = np.array(X), np.array(Y), len(X)
        self.A = np.pi * r**2 / self.n  # (m^2)

    def critical_time(self, P):
        """Time for an elastic wave to cross one grid spacing."""
        return min(P.G.dx, P.G.dy) / np.sqrt(self.K / self.rho)


class Output_Params:
    def __init__(self):
        self.continuum_fig_size = [10, 6]
        self.mp_fig_size = [10, 10]

    def after_every_savetime(self, L, P, G):
        plot.draw_continuum(G, P)
        plot.draw_material_points(L, P, G)
