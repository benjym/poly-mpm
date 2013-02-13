from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp,linspace
from numpy import random, cos, sin, ceil, around, tile, repeat
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-3 # timestep (s)
        self.t = 0. # physical time (s)
        self.tstep = 0
        self.savetime = self.dt
        self.t_f = 0.5#100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.save = 0 # save counter
        self.M_tol = 1e-10 # very small mass (kg)
        self.max_g = 0. # gravity (ms^-2)
        self.max_q = 0.
        self.update_forces()
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.thickness = 1. # (m) into page
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = Solid_Params(self.G)
        self.F = Fluid_Params()
        self.R = Fluid_Params()
        self.has_yielded = False
        self.damping = False
        self.mode = mode

    def update_forces(self):
        t_c = .1
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
#        self.g=self.max_g
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.scale = 20 # grid points per m
        self.x_m = 0.0 # (m)
        self.x_M = 2.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 2.0 # (m)
        self.nx = int(self.x_M-self.x_m)*self.scale+1
        self.ny = int(self.y_M-self.y_m)*self.scale/2+1
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        
class Boundary_Params():
    def __init__(self):
        self.wall = False
        self.has_top = False
        self.has_bottom = False
        self.has_right = False
        self.has_left = False
        self.outlet_left = False
        self.force_boundaries = True
        self.vertical_force = False
        self.horizontal_force = False
        self.roughness = False

class Solid_Params():
    def __init__(self,G):
        self.x = (G.nx-1)*2 # particles in x direction
        self.y = (G.nx-1)*2 # particles in y direction
        self.X = []
        self.Y = []
        self.n = 0

        self.law = 'elastic'
        self.rho = 1000. # density (kg/m^3)
        
        self.E = 1.e7 # elastic modulus (Pa)
        self.nu = 0. # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        gap = G.dx/4.
        xp = linspace(G.x_m+gap,G.x_M-gap,self.x)
        yp = linspace(G.y_m+gap,G.y_M-gap,self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in xrange(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
        
        
class Output_Params():
    def __init__(self,nt):
        self.measure_energy = True
        self.plot_continuum = False
        self.plot_material_points = False
        self.measure_stiffness = False
        self.check_positions = False
        self.plot_fluid = False
        self.plot_gamma_dot = True
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

