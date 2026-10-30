"""Steady flow of a polydisperse granular layer down a periodic slope, with segregation.

python MPM.py bi_seg_test <slope_degrees> <size_classes> <grid_nodes>
"""

import numpy as np


class Params:
    def __init__(self, args):
        self.dt = 1e-8  # (s)
        self.savetime = 1e-3  # (s)
        self.t_f = 5.0  # (s)
        self.gravity = -9.81  # (m/s^2)
        self.gravity_angle = -np.radians(float(args[1]))  # slope angle (rad)
        self.initial_stress = "lithostatic"
        self.initial_flow = "steady"
        self.segregate = True
        self.segregation_rate = 1e-3  # inter-particle drag coefficient
        self.kinetic_pressure_scale = 10.0  # particle diameters in the segregation diffusion coefficient
        self.time_stepping = "dynamic"
        self.cfl = 0.2
        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = Solid_Params(self.G, self)
        self.output_dir = (
            f"im/seg/lin/theta_{-np.rad2deg(self.gravity_angle)}/ny_{self.G.ny}/ns_{len(self.G.grain_sizes)}"
            f"/c_{self.segregation_rate}/l_{self.kinetic_pressure_scale}/"
        )


class Grid_Params:
    def __init__(self, args):
        self.x_min, self.x_max = 0.0, 1.0  # (m)
        self.y_min, self.y_max = 0.0, 0.1  # (m)
        self.nx = 2
        self.ny = int(args[3])
        # Size classes evenly spread between 0.3 mm and 3 mm.
        s_max = 3e-3  # (m)
        s_edges = np.linspace(s_max / 10.0, s_max, int(args[2]) + 1)
        self.grain_sizes = (s_edges[1:] + s_edges[:-1]) / 2.0  # (m)


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom"}
        self.wall_contact = "no_slip_bottom"
        self.periodic = True


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
        self.A = dx * dy / pts_per_cell**2  # (m^2)

    def critical_time(self, P):
        distance = min(P.G.dx, P.G.dy)
        t_ela = distance / np.sqrt(self.K / self.rho)  # elasticity
        t_diff = distance**2 / self.eta_max * self.rho  # momentum diffusivity
        return min(t_diff, t_ela)


class Output_Params:
    def __init__(self):
        self.figures = {"continuum"}
        self.fields = {"s_bar", "u", "density", "phi_material_points"}
        self.continuum_fig_size = [12, 8]
        self.mp_fig_size = [6, 6]
