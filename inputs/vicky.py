import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-4 # timestep (s)
        self.savetime = 0.1
        self.t_f = 2.#100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.max_g = -10. # gravity (ms^-2)
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params(self,self.G)]
        self.supername = 'im/pile/H_'+str(self.S[0].H)+'/L_'+str(self.S[0].L)+'/ny_' + str(self.G.ny)+'/'
        self.damping = False
        
    def update_forces(self):
        t_c = .5
        self.g=self.max_g

class Grid_Params():
    def __init__(self):
        self.x_m = -1.0 # (m)
        self.x_M = 1.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.ny = 10 # number of grid edges in y direction
        self.nx = 1+2*(self.ny-1) # number of grid edges in x direction
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        
class Solid_Params():
    def __init__(self,P,G):
        self.X = [] # list of material point x coordinates
        self.Y = [] # list of material point y coordinates
        self.n = 0 # number of material points to add
        self.rho = 1000. # density (kg/m^3)
        
        self.law = 'von_mises'
        self.E = 1e6 # elastic modulus (Pa) - use dt = 5e-4
        self.nu = 0.3 # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)
        self.yield_stress = 5e3   # Pa
        self.k = self.yield_stress/sqrt(2)
        self.s = 1.

#         self.law = 'dp'
#         self.E = 1e6 # elastic modulus (Pa) - use dt = 5e-4
#         self.nu = 0.3 # poisson's ratio
#         self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
#         self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)
#         self.mu = 1.
#         self.beta = self.mu #0.5
#         self.s = 1.
        
        self.L = 0.5
        self.H = 1.

        self.pts_per_cell = 3
        self.x = int(self.L/(G.x_M-G.x_m)*(G.nx-1))*self.pts_per_cell # particles in x direction
        self.y = int(self.H/(G.y_M-G.y_m)*(G.ny-1))*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        xp = linspace(-self.L/2.+gap[0],self.L/2.-gap[0],self.x)
        yp = linspace(gap[1],self.H-gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = self.L*self.H/self.n # area (m^2)

#         elastic_wave_speed = sqrt(self.E/self.rho)
#         distance = minimum(G.dx,G.dy)
#         critical_time = distance/elastic_wave_speed
#         critical_time /= 10. # safety
#         if critical_time < P.dt:
#             P.dt = critical_time
#             P.nt = int(P.t_f/P.dt)
#             print('Reducing timestep to dt = ' + str(P.dt))        
        
class Output_Params():
    def __init__(self):
        self.continuum_fig_size = [20,12]
        self.mp_fig_size = [10,10]
        self.plot_continuum = True
        self.plot_material_points = False