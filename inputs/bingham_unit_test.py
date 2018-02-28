import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.supername = 'bingham_unit_test'
        self.dt = 2e-6 # timestep (s) --- TRY 0.1/mu_s
        self.t = 0. # physical time (s)
        self.tstep = 0
        self.savetime = 1000*self.dt # 0.1
        self.t_f = 10.0 # 100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.grid_save = 0 # save counter
        self.mp_save = 0 # save counter
        self.M_tol = 1e-10 # very small mass (kg)
        self.max_g = -1. # gravity (ms^-2)
        self.max_q = 0.
        self.update_forces()
        self.theta = 90.*pi/180. # slope angle (degrees)
        self.thickness = 1. # (m) into page
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = [Solid_Params(self.G,self)]
        self.F = Fluid_Params()
        self.R = Fluid_Params()
        self.has_yielded = False
        self.damping = True
        self.mode = mode

    def update_forces(self):
        t_c = 100.*self.dt
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
        self.has_top = True
        self.has_bottom = True
        self.has_right = False
        self.has_left = False
        self.outlet_left = False
        self.cyclic_lr = True
        self.force_boundaries = False
        self.vertical_force = False
        self.horizontal_force = False
        self.roughness = True

class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0

        self.law = 'bingham'
        self.rho = 1000. # density (kg/m^3)
        
        self.mu_s = 1e4 # shear viscosity
        self.tau_0 = 100. # yield stress
        self.mu_0 = 100.*self.mu_s # much steeper viscosity to mimic very rigid part
        self.gamma_c = self.tau_0/(self.mu_0 - self.mu_s)
        self.K = 0. # bulk modulus
        
#         self.initial_v = array([-10.,0.,0.])
#         self.pressure = 1e4
        
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
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
        P.supername += ('/' + str(self.mu_s) +
                        '/' + str(self.mu_0))

class Output_Params():
    def __init__(self,nt):
        self.measure_energy = False
        self.plot_continuum = True
        self.plot_material_points = False
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

            
class Fluid_Params():
    def __init__(self):
        self.n = 0

        
class Rigid_Params():
    def __init__(self):
        self.n = 0

