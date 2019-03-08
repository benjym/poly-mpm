import os
from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp,linspace
from numpy import random, tile, repeat
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-3 # timestep (s)
        self.tstep = 0
        self.savetime = 0.1
        self.t_f = 5.#100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.max_g = 0. # gravity (ms^-2)
        self.max_q = 200.
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = Solid_Params(self.G)

    def update_forces(self):
        t_c = .5
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))


class Grid_Params():
    def __init__(self):
        self.scale = 4
        self.x_m = 0.0 # (m)
        self.x_M = 1.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 2.0 # (m)
        self.nx = int(self.x_M-self.x_m)*self.scale+1
        self.ny = int(self.y_M-self.y_m)*self.scale+1
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

class Boundary_Params():
    def __init__(self):
        self.force_boundaries = True
        self.vertical_force = True
        self.horizontal_force = False
        self.roughness = True

class Solid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.x = (G.nx-1)*2 # particles in x direction
        self.y = (G.ny-1)*2 # particles in y direction
        self.n = 0
        self.mode = 'triaxial'
#        self.law = 'elastic'
        self.law = 'von_mises'
#        self.law = 'dp'
        self.rho = 2650. # density (kg/m^3)

        self.E = 1.e5 # elastic modulus (Pa)
        self.nu = 0.3 # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        self.s = 20.
        self.k = 100.

        self.beta = .5 # dilatancy coefficient
        self.mu = 1. # ?

        gap = G.dx/4.
        xp = linspace(G.x_m+gap,G.x_M-gap,self.x)
        yp = linspace(G.y_m+gap,G.y_M-gap,self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            if i is not (G.scale)*self.x:
                self.X.append(X[i])
                self.Y.append(Y[i])
                self.n += 1
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)


class Output_Params():
    def __init__(self,nt):
        # self.measure_energy = True
        self.plot_material_points = True
        self.plot_continuum = True
        self.measure_stiffness = True
        # self.check_positions = True
        self.energy = zeros((10*nt+1,4)) # energy

    def measure_E(self,L,P,G):
        print('Measuring macro and micro stress/strain for each material point... ')
        for i in range(P.S[0].n):
            original_position = array((P.S[0].X[i],P.S[0].Y[i],0))
            macro_strain = (original_position-L.S[0][i].x)/array((P.G.x_M - P.G.x_m,P.G.y_M - P.G.y_m,1.)) #original_position
            macro_stress = P.max_q/2.
            print('From macroscopic stress/strain:')
            print(macro_stress/macro_strain/P.S[0].E)
            print('From microscopic stress/strain:')
            print(L.S[i].dstress/L.S[i].dstrain/P.S[0].E)
