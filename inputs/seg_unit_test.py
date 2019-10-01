import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-8 # timestep (s)
        self.savetime = 10*self.dt # 0.1
        self.t_f = 100.0 # 100*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)

        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = [Solid_Params(self.G,self)]

        self.segregate_grid = True
        self.c = 1e-3 #1. # inter-particle drag coefficient
        self.D = 1e-3 #1e-6 # shear strain rate to kinetic stress constant

        self.supername = 'im/seg_unit_test/' + '/ny_' + str(self.G.ny) + '/ns_' + str(self.G.ns) + '/c_' + str(self.c) + '/D_' + str(self.D) + '/'
        self.smooth_gamma_dot = False
        self.time_stepping = 'dynamic' # dynamic or static time steps
        self.CFL = 0.2

    def update_forces(self):
        t_c = 0.1
#         self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.g = self.max_g
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 0.25 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = 6
        self.ny = 21
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.thickness = 1. # (m)

class Boundary_Params():
    def __init__(self):
        self.has_top = True
        self.has_bottom = True
        self.cyclic_lr = True
        self.roughness = True

class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0

        self.law = 'viscous'
        self.rho = 1000. # density (kg/m^3)

        self.mu_s = 1e5 # shear viscosity
        self.mu_v = 1e5
        self.K = 0. # pressure viscosity

        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-1)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
#         P.supername += '/' + str(self.mu_s)

class Output_Params():
    def __init__(self,nt):
        self.plot_gsd = True
        self.plot_continuum = True
        self.continuum_fig_size = [10,6]
