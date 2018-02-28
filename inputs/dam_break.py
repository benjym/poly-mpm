import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-4 # timestep (s)
        self.savetime = 0.0001
        self.t_f = 10.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.max_q = 0.
        self.theta = 0.*pi/180. # slope angle (degrees)
#         self.pressure = 'lithostatic'
        self.pressure = 'non-zero'
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()#self.nt)
#         self.S = [Large_Params(self.G,self),Small_Params(self.G,self),]
        self.S = [Solid_Params(self.G,self),]
        self.supername = 'im/dam_break_no_seg/'

    def update_forces(self):
        t_c = 0.05
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
#         self.g = self.max_g

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 2.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.ny = 21
        self.nx = (self.ny-1)*4 + 1
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.s = array([0.5,1.0]) # s coordinate
        self.ds = self.s[1]-self.s[0]
        self.ns = len(self.s)
        
class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.has_right = True
        self.outlet_left = True
#         self.roughness = True
        
class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2650. # density (kg/m^3)

        self.law = 'dp_rate'
        self.E = 1e4 # elastic modulus (Pa)
        self.nu = 0.2 # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)
        self.mu = 1.
        self.beta = 1.*self.mu
        self.s = 1.
        self.t_star = 0.01
        
        self.pts_per_cell = 3
        self.x = (G.nx-1)//4*self.pts_per_cell # particles in x direction
        self.y = (G.ny-2)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)

        xp = linspace(G.x_M-(G.x_m+G.x_M)/4+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-G.dy-gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
#         self.A = self.conc*(G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m-G.dy)/4./self.n # area (m^2)

        elastic_wave_speed = sqrt(self.K/(self.A*self.rho*1.))
        distance = minimum(G.dx,G.dy)
        critical_time = distance/elastic_wave_speed
        critical_time /= 2. # safety
        if critical_time < P.dt:
#             P.dt = critical_time
            print('WARNING: timestep should be dt = ' + str(critical_time))

class Output_Params():
    def __init__(self):
        self.plot_continuum = True
        self.plot_material_points = True
        self.continuum_fig_size = [18,6]
        self.mp_fig_size = [18,4]