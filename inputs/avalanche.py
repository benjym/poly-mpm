"""A viscous layer released on a 45 degree periodic slope with a centre wall."""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s)
        self.savetime = 0.1  # (s)
        self.t_f = 1000.0  # (s)
        self.gravity = -10.0  # (m/s^2)
        self.gravity_angle = np.radians(45.0)  # slope angle (rad)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G)]

    def gravity_ramp(self, t):
        """Gravity rises smoothly from zero over about 0.5 s."""
        return 1.0 - np.exp(-3.0 * t**2 / 0.5**2)


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 4.0  # (m)
        self.y_min, self.y_max = 0.0, 1.0  # (m)
        self.nx = 21
        self.ny = 6


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom"}
        self.wall_contact = "no_slip_bottom"
        self.periodic = True
        self.obstacles = {"centre_wall"}


class Solid_Params:
    def __init__(self, G):
        self.law = "viscous"
        self.rho = 1000.0  # (kg/m^3)
        self.eta = 1e3  # shear viscosity (Pa s)
        self.K = 1e6  # bulk modulus (Pa)
        self.eta_v = 0.0  # bulk viscosity (Pa s)

        # A block filling half the depth of the right quarter of the domain.
        filling_fraction, pts_per_cell = 0.5, 3
        dx = (G.x_max - G.x_min) / (G.nx - 1)
        dy = (G.y_max - G.y_min) / (G.ny - 1)
        nx = (G.nx - 1) * pts_per_cell
        ny = int(round((G.ny - 1) * pts_per_cell * filling_fraction))
        gap = np.array((dx, dy)) / (2 * pts_per_cell)
        xp = np.linspace(3.0 * (G.x_max - G.x_min) / 4.0, G.x_max - gap[0], nx)
        yp = np.linspace(G.y_min + gap[1], gap[1] + filling_fraction * (G.y_max - 2 * gap[1]), ny)
        self.X = np.tile(xp, ny)
        self.Y = np.repeat(yp, nx)
        self.n = len(self.X)
        self.A = (G.x_max - G.x_min) * (G.y_max - G.y_min) / self.n * filling_fraction  # (m^2)

        # Inlet at the right edge; it adds material only with B.inlets = {"right"}.
        self.inlet_rate = 0.1  # (s)
        self.inlet_gap = gap
        self.inlet_y = yp


class Output_Params:
    def __init__(self):
        self.figures = {"continuum", "material_points"}
        self.continuum_fig_size = [18, 6]
        self.mp_fig_size = [6, 6]
