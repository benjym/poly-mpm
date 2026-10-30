"""Steady granular flow down a rough incline with critical-state CIDR.

A layer of glass beads (mu(I) of Jop et al., 3 mm grains), 30 grains deep, on a
no-slip base tilted at 25 degrees, between the static (20.9) and dynamic (32.8)
friction angles, periodic along the slope. In steady isochoric flow CIDR reduces to
mu(I), so every depth flows at the inertial number with mu(I) = tan(theta) and the
velocity follows the Bagnold profile (benchmarks/cidr_incline.py). The layer
starts in that steady state: at its critical packing at each depth, which is
loose only within about a grain of the surface, and at the Bagnold velocity.

    python MPM.py cidr_incline
"""

import numpy as np

GRAIN = 0.003  # (m)
DEPTH = 30 * GRAIN  # (m)
ANGLE = np.radians(25.0)
PHI_0 = 0.6


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s), the first step; dynamic stepping follows
        self.savetime = 0.1  # (s)
        self.t_f = 1.0  # (s)
        self.time_stepping = "dynamic"
        self.cfl = 0.2
        self.gravity = -9.81  # (m/s^2)
        self.gravity_angle = ANGLE
        self.solver = "granular"
        # The Metal smoothing on both backends: at the CPU default (0.075) the shear rates of
        # this flow put Newton's last decrements below round-off in the multiplier.
        self.granular_cidr_gamma = 1.0
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = Solid_Params(self.G, abs(self.gravity))


class Grid_Params:
    def __init__(self):
        cells = 18  # across the layer
        self._spacing = DEPTH / cells
        self.x_min, self.x_max = 0.0, 4 * self._spacing  # (m)
        self.y_min, self.y_max = 0.0, 1.5 * DEPTH  # (m)
        self.nx = 5
        self.ny = int(round(1.5 * cells)) + 1


class Boundary_Params:
    def __init__(self):
        self.walls = {"bottom"}
        self.wall_contact = "no_slip_bottom"
        self.periodic = True


class Solid_Params:
    def __init__(self, G, gravity):
        from cidr import CriticalStateCIDR, packing_modulus

        self.law = "granular"
        self.rho_s = 2500.0  # (kg/m^3)
        self.mu_s = np.tan(np.radians(20.9))
        self.mu_d = np.tan(np.radians(32.76))
        self.I_0 = 0.279
        self.grain_diameter = GRAIN
        self.phi_min, self.phi_max = 0.05, 0.605
        # The packing scale: on the cap at rest at twice the base pressure of a layer at PHI_0.
        self.cidr_deficit = 0.005
        self.cidr_p_ref = 2 * PHI_0 * self.rho_s * gravity * np.cos(ANGLE) * DEPTH
        # Steady flow: every depth at I_theta, mu(I_theta) = tan(theta), and at its critical
        # packing, C(phi) = 2 p / beta(I_theta), with p the weight above (integrated downwards).
        model = CriticalStateCIDR(self.mu_s, self.mu_d, self.I_0, 2 * GRAIN * np.sqrt(self.rho_s))
        I_theta = self.I_0 * (np.tan(ANGLE) - self.mu_s) / (self.mu_d - np.tan(ANGLE))
        beta = 2 * (float(model.alpha(I_theta)) - model.mu(I_theta))
        depth = np.linspace(0.0, DEPTH, 4001)  # below the surface
        pressure, packing = np.zeros_like(depth), np.zeros_like(depth)
        modulus = lambda phi: packing_modulus(  # noqa: E731
            phi, self.phi_min, self.phi_max, self.cidr_deficit, self.cidr_p_ref, self.mu_s
        )
        phis = np.linspace(self.phi_min + 1e-9, self.phi_max - 1e-9, 20001)
        for k in range(len(depth)):
            if k:
                step = depth[k] - depth[k - 1]
                pressure[k] = pressure[k - 1] + self.rho_s * packing[k - 1] * gravity * np.cos(ANGLE) * step
            packing[k] = np.interp(2 * pressure[k] / beta, modulus(phis), phis)
        spacing = G._spacing / 2  # two particles per cell in each direction
        x, y = np.meshgrid(np.arange(spacing / 2, G.x_max, spacing), np.arange(spacing / 2, DEPTH, spacing))
        self.X, self.Y = x.ravel(), y.ravel()
        self.n = self.X.size
        self.A = spacing**2  # (m^2)
        below = DEPTH - self.Y
        self.rho = self.rho_s * np.interp(below, depth, packing)  # (kg/m^3)
        # Bagnold velocity, down the slope (-x): du/dy = I_theta sqrt(p / rho_s) / d.
        rate = I_theta * np.sqrt(pressure / self.rho_s) / GRAIN
        height = DEPTH - depth[::-1]
        speed = np.concatenate(([0.0], np.cumsum(0.5 * (rate[::-1][1:] + rate[::-1][:-1]) * np.diff(height))))
        self.initial_v = np.zeros((self.n, 3))
        self.initial_v[:, 0] = -np.interp(self.Y, height, speed)

    def critical_time(self, P):
        """Time for a grain at rest to fall one grid spacing: rigid grains have no wave speed."""
        return float(np.sqrt(min(P.G.dx, P.G.dy) / abs(P.gravity)))


class Output_Params:
    def __init__(self):
        self.figures = ()
