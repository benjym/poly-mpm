import os
from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp,linspace
from numpy import random
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-3 # timestep (s)
        self.t = 0. # physical time (s)
        self.tstep = 0
        self.savetime = 0.1
        self.t_f = 1.#100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.save = 0 # save counter
        self.M_tol = 1e-10 # very small mass (kg)
        self.max_g = 0. # gravity (ms^-2)
        self.max_q = 1.
        self.update_forces()
        self.theta = 0.*pi/180. # slope angle (degrees)

        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = Solid_Params(self.G)
        self.damping = True
        self.mode = mode

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
        self.force_boundaries = True
        self.horizontal_force = True

class Solid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.x = 3#(10*2-1)*2 # particles in x direction
        self.y = 3#(5*2-1)*2-4 # particles in y direction
        self.n = 0

#        self.law = 'elastic'
        self.law = 'von_mises'
        self.rho = 2650. # density (kg/m^3)
        
        self.E = 1.e8 # elastic modulus (Pa)
        self.nu = 0. # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        self.s = 2.5
        self.k = self.E/100

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
        self.plot_continuum = False
        self.plot_material_points = True
        self.measure_energy = True
        self.measure_stiffness = True

    def measure_E(self,P,L):
        print('Measuring macro and micro stress/strain for each material point... ')
        for p in range(P.phases):
            for i in range(P.S[p].n):
                original_position = array((P.S[p].X[i],P.S[p].Y[i],0))
                macro_stress = P.max_q/2.
                macro_strain = (original_position-L.S[p][i].x)/array((P.S[p].L,P.S[p].W,1.)) #original_position
                print('MP ' + str(i))
                print('Macro response:' + str(macro_stress/macro_strain/P.S[p].E) + '*E')
                print('Micro values:' + str(L.S[p][i].dstress/L.S[p][i].dstrain/P.S[p].E) + '*E')