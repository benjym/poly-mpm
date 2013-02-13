import os
from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp
from numpy import linspace,tile,repeat
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-3 # timestep (s)
        self.t = 0. # physical time (s)
        self.tstep = 0 # current timestep
        self.t_f = 5.#100*self.dt # final time (s)
        self.savetime = .1
        self.nt = int(self.t_f/self.dt) # number of timesteps - WRONG!
        self.save = 0 # save counter
        self.M_tol = 1e-10 # very small mass (kg)
        self.max_g = -10. # gravity (ms^-2)
        self.max_q = 0.
        self.update_forces()
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.thickness = 1. # (m) into page
        self.roughness = 1. # wall roughness (units?)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = Solid_Params(self.G)
        self.F = Fluid_Params(self.G)
        self.R = Rigid_Params()
        self.has_yielded = False
        self.damping = False # local non-viscous damping
        self.mode = mode
        
    def update_forces(self):
        t_c = .5
        if self.t < 4.*t_c:
            self.theta = 0.*pi/180. # slope angle (degrees)
        else:
            self.theta = 45.*pi/180. # slope angle (degrees)
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.scale = 4 # grid points per m
        self.x_m = 0.0 # (m)
        self.x_M = 1.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = int(self.x_M-self.x_m)*self.scale+1
        self.ny = int(self.y_M-self.y_m)*self.scale+1
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

class Boundary_Params():
    def __init__(self):
        self.wall = False
        self.has_top = False
        self.has_bottom = True
        self.has_right = True
        self.has_left = True
        self.outlet_left = False
        self.force_boundaries = False
        self.vertical_force = False
        self.horizontal_force = False
        self.roughness = False

class Solid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.x = G.scale # particles in x direction
        self.y = G.scale # particles in y direction
        self.n = 0

#        self.law = 'elastic'
        self.law = 'von_mises'
        self.rho = 2650. # density (kg/m^3)
        
        self.E = 1.e7 # elastic modulus (Pa)
        self.nu = 0.3 # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        self.s = 100.
        self.k = self.E/100.

#        for i in linspace(G.x_M-2.*G.dx,G.x_M-G.dx,self.x):
#            for j in linspace(G.y_m+G.dy,G.y_m+2.*G.dy,self.y):
#                self.X.append(i)
#                self.Y.append(j)
#                self.n += 1
        self.A = 1. # area (m^2)

class Output_Params():
    def __init__(self,nt):
        self.measure_energy = True
        self.plot_continuum = False
        self.plot_material_points = False
        self.plot_fluid = True
        self.plot_gamma_dot = False
        self.measure_stiffness = False
        self.check_positions = True
        self.energy = zeros((10*nt+1,4)) # DON'T KNOW NUMBER OF TIMESTEPS!

class Fluid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.x = G.scale # particles in x direction
        self.y = G.scale # particles in y direction
        self.n = 0
        self.rho = 999.
        self.P_0 = 200.
        self.int_energy = 1.
        self.gamma = 1.
        self.mu = 1.
#        for i in linspace(G.x_M-2.*G.dx,G.x_M-G.dx,self.x):
#            for j in linspace(G.y_m,G.y_m+G.dy,self.y):
        for i in linspace(G.x_m+0.01,G.x_M-0.01,G.scale*3):
            for j in linspace(G.y_m+0.01,G.y_M-0.01,G.scale*3):
                self.X.append(i)
                self.Y.append(j)
                self.n += 1
        self.A = 1. # area (m^2)
        self.law = 'compressible'
                
class Rigid_Params():
    def __init__(self):
        self.n = 0

