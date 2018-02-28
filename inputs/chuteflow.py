import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-3 # timestep (s)
        self.t = 0. # physical time (s)
        self.tstep = 0
        self.savetime = 0.1
        self.t_f = 5.0 # 100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.grid_save = 0 # save counter
        self.mp_save = 0 # save counter
        self.M_tol = 1e-10 # very small mass (kg)
        self.max_g = -1. # gravity (ms^-2)
        self.max_q = 0.
        self.update_forces()
        self.theta = 45.*pi/180. # slope angle (degrees)
        self.thickness = 1. # (m) into page
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = []
        self.S.append(Solid_Params(self.G))
        self.has_yielded = False
        self.damping = True
        self.mode = mode
        self.supername = 'chute_flow'

    def update_forces(self):
        t_c = 0.5
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
#         self.g = self.max_g
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 0.1 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = 3
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
        self.inlet_right = False
        self.cyclic_lr = True
        self.force_boundaries = False
        self.vertical_force = False
        self.horizontal_force = False
        self.roughness = True

class Solid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.n = 0

        self.law = 'viscous'
#         self.law = 'HB'
        self.rho = 1000. # density (kg/m^3)
        
        self.mu_s = 1e3 # shear viscosity
        self.mu_v = 0 # pressure viscosity
        self.K = 1e3
#         self.k = 1.
#         self.N = 2.
#         self.tau_0 = 1.
#         self.gamma_dot_0 = 0.1
#         self.mu_0 = self.k*self.gamma_dot_0**(self.N-1.) + self.tau_0/self.gamma_dot_0
        
        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-1)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in xrange(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/2./self.n # area (m^2)
        
        
class Output_Params():
    def __init__(self,nt):
        self.measure_energy = False
        self.plot_continuum = True
        self.plot_material_points = True
        self.measure_stiffness = False
        self.check_positions = False
        self.plot_fluid = False
        self.energy = zeros((nt+1,4)) # energy
        
    def measure_E(self,P,L):
        print 'Measuring macro and micro stress/strain for each material point... '
        for i in xrange(P.S.n):
            original_position = array((P.S.X[i],P.S.Y[i],0))
            macro_strain = (original_position-L.S[i].x)/array((P.S.L,P.S.W,1.)) #original_position
            macro_stress = P.max_conf#/(P.S.L*P.S.W)
            print 'From macroscopic stress/strain:'
            print macro_stress/macro_strain/P.S.E
            print 'From microscopic stress/strain:'
            print L.S[i].dstress/L.S[i].dstrain/P.S.E
