"""A von Mises sample loaded by a vertical traction that rises smoothly to 200 Pa."""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-3  # (s)
        self.savetime = 0.1  # (s)
        self.t_f = 5.0  # (s)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = Solid_Params(self.G)


_SCALE = 4  # grid nodes per metre


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 1.0  # (m)
        self.y_min, self.y_max = 0.0, 2.0  # (m)
        self.nx = int(self.x_max - self.x_min) * _SCALE + 1
        self.ny = int(self.y_max - self.y_min) * _SCALE + 1


class Boundary_Params:
    def __init__(self):
        self.wall_contact = "no_slip"

    def traction_vertical(self, t):
        """Rises smoothly from zero to 200 Pa over about 0.5 s."""
        return 200.0 * (1.0 - np.exp(-3.0 * t**2 / 0.5**2))


class Solid_Params:
    def __init__(self, G):
        self.law = "von_mises"
        self.rho = 2650.0  # (kg/m^3)
        E, nu = 1.0e5, 0.3  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.G = E / (2.0 * (1.0 + nu))  # shear modulus (Pa)
        self.yield_stress = np.sqrt(3) * 100.0  # von Mises equivalent stress at yield (Pa)

        # Two particles per cell in each direction, less the first particle of row _SCALE.
        dx = (G.x_max - G.x_min) / (G.nx - 1)
        nx, ny = (G.nx - 1) * 2, (G.ny - 1) * 2
        gap = dx / 4.0
        xp = np.linspace(G.x_min + gap, G.x_max - gap, nx)
        yp = np.linspace(G.y_min + gap, G.y_max - gap, ny)
        keep = np.arange(nx * ny) != _SCALE * nx
        self.X = np.tile(xp, ny)[keep]
        self.Y = np.repeat(yp, nx)[keep]
        self.n = len(self.X)
        self.A = (G.x_max - G.x_min) * (G.y_max - G.y_min) / self.n  # (m^2)


class Output_Params:
    def __init__(self):
        self.figures = {"material_points", "continuum"}
        self.measure_stiffness = True
