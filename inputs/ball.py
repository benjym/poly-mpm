import os
import numpy as np
import matplotlib.pyplot as plt
from plotting import Plotting
plot = Plotting()

class Params():
    def __init__(self,mode):
        self.dt = 1e-3 # timestep (s)
        self.savetime = 0.1
        self.t_f = 1. #self.dt
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.max_g = -10. # gravity (ms^-2)
        self.theta = 0.*np.pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = Solid_Params()
        # self.time_stepping = 'dynamic'

    def update_forces(self):
        self.g=self.max_g

class Grid_Params():
    def __init__(self):
        self.x_m = -1.0 # (m)
        self.x_M = 1.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 2.0 # (m)
        self.nx = 6 # number of grid edges in x direction
        self.ny = 6 # number of grid edges in y direction
        self.x = np.linspace(self.x_m,self.x_M,self.nx)
        self.y = np.linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True

class Solid_Params():
    def __init__(self):
        self.X = []
        self.Y = []
        self.n = 0

        self.rho = 1000. # density (kg/m^3)

        self.law = 'elastic'
        # self.law = 'von_mises'
#         self.law = 'dp'

        self.E = 1.e5 # elastic modulus (Pa)
        self.nu = 0.3 # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        # for von_mises
        self.k = self.E/100.

        # for dp
        self.s = 5.
        self.beta = 10.
        self.mu = 0.5

        nr = 20 # particles in radial direction
        nphi = 50 # particles around circumference
        r = 0.3 # radius
        c = [0.,1.0] # centre

        for i in np.linspace(0,r,nr):
            dnphi = int(np.around(nphi*i/r)) # number in this ring
            if dnphi > 0:
                for j in np.linspace(0,(1.-1./(dnphi))*2*np.pi,dnphi):
                    self.X.append(c[0]+i*np.sin(j))
                    self.Y.append(c[1]+i*np.cos(j))
                    self.n += 1
        self.A = np.pi*r**2/self.n # area (m^2)

    def critical_time(self,P):
        distance = np.minimum(P.G.dx,P.G.dy)
        t_ela = distance/np.sqrt(self.K/self.rho) # elasticity
        return t_ela

class Output_Params():
    def __init__(self):
        self.continuum_fig_size = [10,6]
        self.mp_fig_size = [10,10]
    def after_every_savetime(self,L,P,G):
        plot.draw_continuum(G,P)
        plot.draw_material_points(L,P,G)
