import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-5 # timestep (s)
        self.savetime = 1e2*self.dt#0.0001
        self.t_f = 10.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.theta = -float(args[1])*pi/180. # slope angle (degrees)
        self.pressure = 'lithostatic'
        
        self.segregate_grid = True
#         self.segregate_grid = False
        self.c = 1.0 # inter-particle drag coefficient
        self.D = 1e-2 # segregation diffusion coefficient
        
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self)#self.nt)
#         self.S = [Large_Params(self.G,self),Small_Params(self.G,self),]
        self.S = [Solid_Params(self.G,self),]

        self.supername = 'im/seg/' + str(self.G.ns) + '/theta_' + args[1] + '/'
        print(self.supername)

    def update_forces(self):
#         t_c = 0.05
#         self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.g = self.max_g
#         self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))
#         if self.t > 0.01:
#             self.segregate_grid = True
#             self.plot_gsd_grid = True

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 0.02 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = 2
        self.ny = 201
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.s = array([0.2,0.4,0.6,0.8,1.0]) # s coordinate
        self.ds = self.s[1]-self.s[0]
        self.ns = len(self.s)
        
class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.cyclic_lr = True
        self.roughness = True
        
class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700. # density (kg/m^3)
        self.phi = [0.2,0.2,0.2,0.2,0.2]
        self.law = 'viscous'
#         mu_s = 1e3
#         self.mu = [mu_s*P.Mu, mu_s]
#         self.law = 'viscous_size'
#         self.law = 'shear_maxwell'

        self.mu_s = 1e3
#         self.mu_s = self.phi[0]*self.mu[0] + self.phi[1]*self.mu[1]
        self.mu_v = 0
        self.E = 1e6
        self.nu = 0.499 # poissons ratio
        self.K = self.E/(3*(1-2*self.nu))
        self.G = self.E/(2*(1+self.nu))
        
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
#         self.A = self.conc*(G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)

#         elastic_wave_speed = sqrt(self.K/(self.A*self.rho*1.))
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
#         self.plot_material_points = True
        self.plot_gsd_grid = True
        self.continuum_fig_size = [10,4]