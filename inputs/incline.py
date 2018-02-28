import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-5 # timestep (s)
        self.savetime = 500*self.dt
        self.t_f = 5.0 # 3*self.dt # final time (s)
        self.max_g = -10. # gravity (ms^-2)
        self.max_q = 0.
        self.theta = 18.*pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()#self.nt)
        self.S = [Solid_Params(self.G,self),]

    def update_forces(self):
        t_c = 0.05
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
#         self.g = self.max_g
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 0.1 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = 2
        self.ny = 21
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.thickness = 1. # (m)
        
class Boundary_Params():
    def __init__(self):
        self.wall = False
        self.has_top = False
        self.has_bottom = True
        self.has_right = False
        self.has_left = False
        self.outlet_left = False
        self.cyclic_lr = True
        self.force_boundaries = False
        self.vertical_force = False
        self.horizontal_force = False
        self.roughness = True
        self.inlet_right = False
        
class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 1000. # density (kg/m^3)
        
        self.law = 'viscous'
        self.mu_s = 1e4
        self.K = 1e6
        self.mu_v = 0 #1e4

#         self.law = 'bingham'
#         self.rho = 1000. # density (kg/m^3)
#         self.mu_s = 1e2 # shear viscosity
#         self.tau_0 = 100. # yield stress
#         self.mu_0 = 100.*self.mu_s # much steeper viscosity to mimic very rigid part
#         self.gamma_c = self.tau_0/(self.mu_0 - self.mu_s)
#         self.K = 1e6 # bulk modulus
#         self.mu_v = 0 #1e3

#         self.law = 'pouliquen'
#         self.K = 1e4
#         self.mu_s = tan(20.*pi/180.)
#         self.mu_2 = tan(30.*pi/180.)
#         self.mu_v = 1e4
#         self.I_0 = 0.279
        
        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-2)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-G.dy-gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)

        elastic_wave_speed = sqrt(self.K/(self.A*self.rho*G.thickness))
        distance = minimum(G.dx,G.dy)
        critical_time = distance/elastic_wave_speed
        critical_time /= 2. # safety
        if critical_time < P.dt:
            P.dt = critical_time
            print('Reducing timestep to dt = ' + str(P.dt))
        
class Output_Params():
    def __init__(self):
        self.measure_energy = False
        self.plot_continuum = True
        self.plot_material_points = True
        self.measure_stiffness = False
        self.check_positions = False
        self.plot_fluid = False
        self.continuum_fig_size = [10,4]
            
class Fluid_Params():
    def __init__(self):
        self.n = 0

        
class Rigid_Params():
    def __init__(self):
        self.n = 0