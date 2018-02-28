import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-6 # timestep (s)
        self.savetime = 0.001
        self.t_f = 10.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.pressure = 'lithostatic'

        self.segregate_grid = True
#         self.segregate_grid = False
        self.c = 1e-3 # inter-particle drag coefficient
        self.D = 0 # segregation diffusion coefficient

        self.supername = 'im/silo/'
        print(self.supername)

        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self)
        self.S = [Solid_Params(self.G,self),]

    def update_forces(self):
        t_c = 0.05
        # self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.g = self.max_g
#         self.g = 0.

class Grid_Params():
    def __init__(self):
        self.x_m = -0.1 # (m)
        self.x_M =  0.1 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 0.4 # (m)
        self.nx = 11
        self.ny = 21
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.s = array([0.5,1.0]) # s coordinate
        self.ds = self.s[1]-self.s[0]
        self.ns = len(self.s)

class Boundary_Params():
    def __init__(self):
        self.has_right = True
        self.has_left = True
        self.silo_bottom = True
        self.outlet_bottom = True
#         self.roughness = True

class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700. # density (kg/m^3)
        self.packing = 0.6 # packing fraction
        self.rho_s = self.rho/self.packing # solid density
        self.phi = ones([G.ns])/float(G.ns)
        self.law = 'pouliquen'
        self.mu_0 = tan(deg2rad(20.9)) #0.3
        self.mu_1 = tan(deg2rad(32.76))
        self.delta_mu = self.mu_1 - self.mu_0
        self.I_0 = 0.279

        self.mu_v = 0
        self.E = 1e6
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3.*(1.-2.*self.nu))
        self.G = self.E/(2.*(1.+self.nu))
        self.eta_max = 100.*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)

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
        self.A = G.dx*G.dy/self.pts_per_cell**2 # area (m^2)

        elastic_wave_speed = sqrt(self.K/self.rho)
        distance = minimum(G.dx,G.dy)
        critical_time = distance/elastic_wave_speed
        critical_time /= 2. # safety
        if critical_time < P.dt:
#             P.dt = critical_time
            print('WARNING: timestep should be dt = ' + str(critical_time))

class Output_Params():
    def __init__(self,P):
        self.plot_continuum = True
        self.plot_material_points = True
        # self.plot_gsd_grid = True
        self.continuum_fig_size = [10,10]
        self.mp_fig_size = [10,8]
