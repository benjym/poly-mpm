import os
from numpy import pi, array, dot, tan, zeros, eye, cumsum, outer, sum, trace, exp, linspace
from numpy import random
import matplotlib.pyplot as plt


class Params:
    def __init__(self, mode):
        self.dt = 1e-3  # timestep (s)
        self.t = 0.0  # physical time (s)
        self.tstep = 0
        self.savetime = 0.1
        self.t_f = 1.0  # 100*self.dt # final time (s)
        self.nt = int(self.t_f / self.dt)  # number of timesteps
        self.save = 0  # save counter
        self.M_tol = 1e-10  # very small mass (kg)
        self.max_g = 0.0  # gravity (ms^-2)
        self.max_q = 1.0
        self.update_forces()
        self.theta = 0.0 * pi / 180.0  # slope angle (degrees)

        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = Solid_Params(self.G)
        self.damping = True
        self.mode = mode

    def update_forces(self):
        t_c = 0.05
        self.g = self.max_g * (1.0 - exp(-3.0 * self.t**2 / t_c**2))
        self.q_v = self.q_h = self.max_q * (1.0 - exp(-3.0 * self.t**2 / t_c**2))


class Grid_Params:
    def __init__(self):
        self.x_m = 0.0  # (m)
        self.x_M = 1.0  # (m)
        self.y_m = 0.0  # (m)
        self.y_M = 1.0  # (m)
        self.nx = 2
        self.ny = 2


class Boundary_Params:
    def __init__(self):
        self.force_boundaries = True
        self.horizontal_force = True


class Solid_Params:
    def __init__(self, G):
        self.X = []
        self.Y = []
        self.x = 3  # (10*2-1)*2 # particles in x direction
        self.y = 3  # (5*2-1)*2-4 # particles in y direction
        self.n = 0

        #        self.law = 'elastic'
        self.law = "von_mises"
        self.rho = 2650.0  # density (kg/m^3)

        self.E = 1.0e8  # elastic modulus (Pa)
        self.nu = 0.0  # poisson's ratio
        self.K = self.E / (3.0 * (1.0 - 2.0 * self.nu))  # bulk modulus (Pa)
        self.G = self.E / (2.0 * (1.0 + self.nu))  # shear modulus (Pa)

        self.s = 2.5
        self.k = self.E / 100

        self.L = 0.9
        self.W = 0.9
        for x in linspace(0.5 - self.L / 2.0, 0.5 + self.L / 2.0, self.x):
            for y in linspace(0.5 - self.W / 2.0, 0.5 + self.W / 2.0, self.y):
                self.X.append(x)
                self.Y.append(y)
                self.n += 1

        self.A = (G.x_M - G.x_m) * (G.y_M - G.y_m) / self.n  # area (m^2)


class Output_Params:
    def __init__(self, nt):
        self.plot_continuum = False
        self.plot_material_points = True
        self.measure_energy = True
        self.measure_stiffness = True

    def measure_E(self, L, P, G):
        from plotting import Plotting

        Plotting().measure_E(L, P, G)
