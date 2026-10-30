"""A half-filled rotating circular drum of bidisperse grains, the drum wall a rigid material.

python MPM.py bi_drum <grid_nodes> <Froude_number>
"""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s)
        self.savetime = 1e-1  # (s)
        self.t_f = 10.0  # (s)
        self.gravity = -9.81  # (m/s^2)
        self._Fr = float(args[2])  # Froude number
        self._r = 0.1  # drum radius (m)
        self.segregate = True
        self.segregation_rate = 1e-6  # inter-particle drag coefficient
        self.initial_stress = "lithostatic"
        self.G = Grid_Params(self, args)
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G, self), Drum_Params(self.G, self)]
        self.output_dir = f"im/drum/ny_{self.G.ny}/"

    def gravity_angle(self, t):
        """Gravity turns once per rotation period, which the Froude number sets."""
        period = np.sqrt(4 * np.pi**2 * self._r / (abs(self.gravity) * self._Fr))
        return (t / period) * 2.0 * np.pi


class Grid_Params:
    def __init__(self, P, args):
        self.x_min, self.x_max = -P._r * 1.1, P._r * 1.1  # (m)
        self.y_min, self.y_max = -P._r * 1.1, P._r * 1.1  # (m)
        self.ny = int(args[1])
        self.nx = self.ny
        # 3 mm and 0.3 mm beads (https://journals.aps.org/pre/pdf/10.1103/PhysRevE.62.961).
        self.grain_sizes = np.array([0.003 / 10.0, 0.003])  # (m)


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom", "top", "left", "right"}
        self.wall_contact = "no_slip"


def _lattice(G, rows_below_top, pts_per_cell=3):
    """Regular particle positions over the grid, pts_per_cell per cell in each direction."""
    dx = (G.x_max - G.x_min) / (G.nx - 1)
    dy = (G.y_max - G.y_min) / (G.ny - 1)
    nx = (G.nx - 1) * pts_per_cell
    ny = (G.ny - rows_below_top) * pts_per_cell
    gap = np.array((dx, dy)) / (2 * pts_per_cell)
    xp = np.linspace(G.x_min + gap[0], G.x_max - gap[0], nx)
    yp = np.linspace(G.y_min + gap[1], G.y_max - gap[1], ny)
    return np.tile(xp, ny), np.repeat(yp, nx), dx * dy / pts_per_cell**2


class Solid_Params:
    def __init__(self, G, P):
        self.law = "pouliquen"
        self.rho = 2700.0  # (kg/m^3)
        self.phi_flow = 0.6  # packing fraction
        self.rho_s = self.rho / self.phi_flow  # solid density (kg/m^3)
        self.mu_s = np.tan(np.radians(20.9))
        self.mu_d = np.tan(np.radians(32.76))
        self.I_0 = 0.279
        self.eta_max = 100.0 * self.rho * np.sqrt(-P.gravity * (G.y_max - G.y_min) ** 3)
        E, nu = 1e6, 0.4  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3 * (1 - 2 * nu))  # bulk modulus (Pa)

        # The lower half of the drum, centred at the origin.
        X, Y, _ = _lattice(G, 1)
        fill = 0.5
        inside = (X**2 + Y**2 < P._r**2) & ((Y / P._r + 1.0) / 2.0 < fill)
        self.X, self.Y = X[inside], Y[inside]
        self.n = len(self.X)
        self.A = np.pi * P._r**2 / self.n  # (m^2)


class Drum_Params:
    def __init__(self, G, P):
        self.law = "rigid"
        self.rho = 2700.0  # (kg/m^3)
        X, Y, self.A = _lattice(G, 2)
        outside = X**2 + Y**2 >= P._r**2
        self.X, self.Y = X[outside], Y[outside]
        self.n = len(self.X)


class Output_Params:
    def __init__(self):
        self.figures = {"continuum", "material_points"}
        self.fields = {"s_bar", "u"}
        self.continuum_fig_size = [24, 8]
        self.mp_fig_size = [18, 4]
