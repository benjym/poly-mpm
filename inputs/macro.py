import os
from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp,linspace
from numpy import random
import matplotlib.pyplot as plt

class Params():
    def __init__(self):
        self.dt = 1e-3 # timestep (s)
        self.t = 0. # physical time (s)
        self.tstep = 0
        self.savetime = 1.
        self.t_f = 0.2#100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.save = 0 # save counter
        self.M_tol = 1e-10 # very small mass (kg)
        self.max_g = 0. # gravity (ms^-2)
        self.max_q = 1.
        self.update_forces()
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.thickness = 1. # (m) into page
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = Solid_Params(self.G)
        self.F = Fluid_Params()
        self.R = Fluid_Params()
        self.check_positions = False
        self.has_yielded = False
        self.damping = True

    def update_forces(self):
        t_c = .05
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 1.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = 2
        self.ny = 2
        
class Boundary_Params():
    def __init__(self):
        self.wall = False
        self.has_top = False
        self.has_bottom = False
        self.has_right = False
        self.has_left = False
        self.outlet_left = False
        self.outlet_bottom = False
        self.force_boundaries = True
        self.vertical_force = False
        self.horizontal_force = True

class Solid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.x = 2#(10*2-1)*2 # particles in x direction
        self.y = 2#(5*2-1)*2-4 # particles in y direction
        self.n = 0

        self.mode = 'blob' # ball, blob, slope, excavate
        self.law = 'elastic'
#        self.law = 'von_mises'
        self.rho = 2650. # density (kg/m^3)
        
        self.E = 1.e8 # elastic modulus (Pa)
        self.nu = 0. # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        self.s = 2.5
        self.k = self.E/100.

        self.L = 0.9
        self.W = 0.9
        for x in linspace(0.5-self.L/2.,0.5+self.L/2.,self.x):
            for y in linspace(0.5-self.W/2.,0.5+self.W/2.,self.y):
                self.X.append(x)
                self.Y.append(y)
                self.n += 1

        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)

class Output_Params():
    def __init__(self,nt):
        self.measure_energy = True
        self.plot_continuum = False
        self.plot_material_points = False
        self.measure_stiffness = True
        self.check_positions = False
        self.plot_fluid = False
        self.energy = zeros((nt+1,3)) # energy
        
    def measure_E(self,P,L):
        print('Measuring macro and micro stress/strain for each material point... ')
        for i in range(P.S.n):
            original_position = array((P.S.X[i],P.S.Y[i],0))
            macro_strain = (original_position-L.S[i].x)/array((P.S.L,P.S.W,1.)) #original_position
            macro_stress = P.max_q/2.
            print('From macroscopic stress/strain:' +
                  str(macro_stress/macro_strain/P.S.E) +
                  'From microscopic stress/strain:' +
                  str(L.S[i].dstress/L.S[i].dstrain/P.S.E)
                  )
                  
class Fluid_Params():
    def __init__(self):
        self.n = 0
        
class Rigid_Params():
    def __init__(self):
        self.n = 0

