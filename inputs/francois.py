"""An elastic disc intruder buried in a periodic granular layer.

python MPM.py francois <grid_nodes>
"""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-5  # (s)
        self.savetime = 0.001  # (s)
        self.t_f = 10.0  # (s)
        self.gravity = -9.81  # (m/s^2)
        self.initial_stress = "lithostatic"
        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G, self), Intruder_Params(self.G)]
        self.output_dir = f"im/francois/ny_{self.G.ny}/"


class Grid_Params:
    def __init__(self, args):
        self.x_min, self.x_max = 0.0, 0.5  # (m)
        self.y_min, self.y_max = 0.0, 0.1  # (m)
        self.ny = int(args[1])
        self.nx = 5 * (self.ny - 1) + 1
        self.top_gap = (self.y_max - self.y_min) / (self.ny - 1)  # one empty cell at the top (m)


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom"}
        self.periodic = True


def _lattice(G, inside, pts_per_cell=3):
    """Particles inside (or outside) the intruder, a disc of radius 0.01 m at mid-depth."""
    dx = (G.x_max - G.x_min) / (G.nx - 1)
    dy = (G.y_max - G.y_min) / (G.ny - 1)
    radius, x_c, y_c = 0.01, (G.x_min + G.x_max) / 2.0, (G.y_min + G.y_max - dy) / 2.0
    nx = (G.nx - 1) * pts_per_cell
    ny = (G.ny - 2) * pts_per_cell
    gap = np.array((dx, dy)) / (2 * pts_per_cell)
    xp = np.linspace(G.x_min + gap[0], G.x_max - gap[0], nx)
    yp = np.linspace(G.y_min + gap[1], G.y_max - gap[1] - G.top_gap, ny)
    X, Y = np.tile(xp, ny), np.repeat(yp, nx)
    keep = (X - x_c) ** 2 + (Y - y_c) ** 2 <= radius**2
    if not inside:
        keep = ~keep
    return X[keep], Y[keep], dx * dy / pts_per_cell**2


class Solid_Params:
    def __init__(self, G, P):
        self.law = "pouliquen"
        self.rho = 2700.0  # (kg/m^3)
        self.phi_flow = 0.6  # packing fraction
        self.rho_s = self.rho / self.phi_flow  # solid density (kg/m^3)
        self.mu_s = np.tan(np.radians(20.9))
        self.mu_d = np.tan(np.radians(32.76))
        self.I_0 = 0.279
        E, nu = 1e4, 0.4  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.eta_max = self.rho * np.sqrt(-P.gravity * (G.y_max - G.y_min) ** 3)
        self.X, self.Y, self.A = _lattice(G, inside=False)
        self.n = len(self.X)


class Intruder_Params:
    def __init__(self, G):
        self.law = "elastic"
        self.rho = 2700.0  # (kg/m^3)
        E, nu = 1e4, 0.4  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.G = E / (2.0 * (1.0 + nu))  # shear modulus (Pa)
        self.X, self.Y, self.A = _lattice(G, inside=True)
        self.n = len(self.X)


class Output_Params:
    def __init__(self):
        self.figures = {"continuum"}
        self.continuum_fig_size = [22, 4]
        self.mp_fig_size = [6, 6]
