"""A half-filled rotating square drum of bidisperse grains with no-slip walls.

python MPM.py bi_square_drum <grid_nodes> <Froude_number>
"""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s)
        self.savetime = 1e-2  # (s)
        self.t_f = 100.0  # (s)
        self.gravity = -9.81  # (m/s^2)
        self._Fr = float(args[2])  # Froude number
        self._r = 1.0  # drum size (m)
        self.segregate = True
        self.segregation_rate = 1e-3  # inter-particle drag coefficient
        self.kinetic_pressure_scale = 1.0  # particle diameters in the segregation diffusion coefficient
        self.initial_stress = "lithostatic"
        self.time_stepping = "dynamic"
        self.segregation_normalise = True
        self.cfl = 0.2
        self.G = Grid_Params(self, args)
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self.G, self)]
        self.output_dir = f"im/drum/ny_{self.G.ny}/Fr_{self._Fr}/"

    def gravity_angle(self, t):
        """Gravity turns once per rotation period, which the Froude number sets."""
        period = np.sqrt(4 * np.pi**2 * self._r / (abs(self.gravity) * self._Fr))
        return -(t / period) * 2.0 * np.pi


class Grid_Params:
    def __init__(self, P, args):
        self.x_min, self.x_max = 0.0, P._r  # (m)
        self.y_min, self.y_max = 0.0, P._r  # (m)
        self.ny = int(args[1])
        self.nx = self.ny
        self.top_gap = 0.5 * (self.y_max - self.y_min)  # (m)
        # 3 mm and 0.3 mm beads (https://journals.aps.org/pre/pdf/10.1103/PhysRevE.62.961).
        self.grain_sizes = np.array([0.003 / 10.0, 0.003])  # (m)


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom", "top", "left", "right"}
        self.wall_contact = "no_slip"


class Solid_Params:
    def __init__(self, G, P):
        self.law = "pouliquen"
        self.rho = 2700.0  # (kg/m^3)
        self.phi_flow = 0.6  # packing fraction
        self.rho_s = self.rho / self.phi_flow  # solid density (kg/m^3)
        self.mu_s = np.tan(np.radians(20.9))
        self.mu_d = np.tan(np.radians(32.76))
        self.I_0 = 0.279
        self.eta_max = self.rho * np.sqrt(-P.gravity * (G.y_max - G.y_min) ** 3)
        E, nu = 1e7, 0.4  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3 * (1 - 2 * nu))  # bulk modulus (Pa)

        # The lower half of the drum.
        pts_per_cell = 3
        dx = (G.x_max - G.x_min) / (G.nx - 1)
        dy = (G.y_max - G.y_min) / (G.ny - 1)
        nx = (G.nx - 1) * pts_per_cell
        ny = (G.ny - 1) // 2 * pts_per_cell
        gap = np.array((dx, dy)) / (2 * pts_per_cell)
        xp = np.linspace(G.x_min + gap[0], G.x_max - gap[0], nx)
        yp = np.linspace(G.y_min + gap[1], G.y_max - gap[1] - G.top_gap, ny)
        self.X = np.tile(xp, ny)
        self.Y = np.repeat(yp, nx)
        self.n = len(self.X)
        self.A = (G.y_max - G.y_min - G.top_gap) * (G.x_max - G.x_min) / self.n  # (m^2)

    def critical_time(self, P):
        distance = min(P.G.dx, P.G.dy)
        t_ela = distance / np.sqrt(self.K / self.rho)  # elasticity
        t_diff = distance**2 / self.eta_max * self.rho  # momentum diffusivity
        return min(t_diff, t_ela)


class Output_Params:
    def __init__(self):
        self.figures = {"continuum", "material_points", "gsd_material_points", "gsd_grid"}
        self.fields = {"s_bar", "u", "density", "phi_material_points"}
        self.continuum_fig_size = [24, 8]
        self.mp_fig_size = [18, 4]
