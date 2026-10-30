from plotting import PQHistory
import os
from numpy import (
    pi,
    array,
    dot,
    tan,
    zeros,
    eye,
    cumsum,
    outer,
    sum,
    trace,
    exp,
    linspace,
    minimum,
    amax,
    tile,
    repeat,
)
from numpy import random, ones
import matplotlib.pyplot as plt


class Params:
    def __init__(self, mode):
        self.dt = 1e-3  # timestep (s)
        self.savetime = 0.1
        self.t_f = 5.0  # final time (s)
        self.nt = int(self.t_f / self.dt)  # number of timesteps
        self.max_q = 10.0
        self.max_g = 0.0
        self.pressure = "non-zero"
        self.theta = 0.0 * pi / 180.0  # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.S = Solid_Params(self.G)
        self.O = Output_Params(self.nt, self.S.n)
        self.supername = "im/dp_rate_unit_test"

    #         self.damping = True

    def update_forces(self):
        self.g = 0
        t_c = self.t_f / 2.0
        if self.t < t_c:
            self.q_h = self.q_v = self.max_q * self.t / t_c
        else:
            self.q_h = self.t * self.max_q / t_c
            self.q_v = 2 * self.max_q - self.t / t_c * self.max_q


#         self.q_h = self.t*self.max_q/t_c
#         self.q_v = minimum(self.q_h,self.max_q)


class Grid_Params:
    def __init__(self):
        self.x_m = 0.0  # (m)
        self.x_M = 1.0  # (m)
        self.y_m = 0.0  # (m)
        self.y_M = 1.0  # (m)
        self.nx = 2
        self.ny = 2
        self.x = linspace(self.x_m, self.x_M, self.nx)
        self.y = linspace(self.y_m, self.y_M, self.ny)
        self.dx = self.x[1] - self.x[0]  # grid spacing (m)
        self.dy = self.y[1] - self.y[0]  # grid spacing (m)


class Boundary_Params:
    def __init__(self):
        #         self.has_right = True
        #         self.has_left = True
        self.vertical_force = True
        self.horizontal_force = True


class Solid_Params:
    def __init__(self, G):
        self.X = []
        self.Y = []
        self.x = 3  # (10*2-1)*2 # particles in x direction
        self.y = 3  # (5*2-1)*2-4 # particles in y direction
        self.n = 0

        self.law = "dp_rate"
        #         self.law = 'rigid'
        self.rho = 2650.0  # density (kg/m^3)

        self.E = 1e5  # elastic modulus (Pa)
        self.nu = 0.2  # poisson's ratio
        self.K = self.E / (3.0 * (1.0 - 2.0 * self.nu))  # bulk modulus (Pa)
        self.G = self.E / (2.0 * (1.0 + self.nu))  # shear modulus (Pa)

        self.mu_s = 1.0
        self.beta = 1.0 * self.mu_s
        self.t_star = 1e-2

        self.pts_per_cell = 3
        self.x = (G.nx - 1) * self.pts_per_cell  # particles in x direction
        self.y = (G.ny - 1) * self.pts_per_cell  # particles in y direction
        gap = array((G.dx, G.dy)) / (2 * self.pts_per_cell)
        xp = linspace(G.x_m + gap[0], G.x_M - gap[0], self.x)
        yp = linspace(G.y_m + gap[1], G.y_M - gap[1], self.y)
        X = tile(xp, self.y)
        Y = repeat(yp, self.x)
        for i in range(self.x * self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        #         self.A = self.conc*(G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
        self.A = (G.x_M - G.x_m) * (G.y_M - G.y_m) / self.n  # area (m^2)


class Output_Params(PQHistory):
    pass
