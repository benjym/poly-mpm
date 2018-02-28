import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.t = 0.
        self.dt = 1e-5 # timestep (s)
        self.savetime = 0.001 #1e1*self.dt#0.01
        self.t_f = 10.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.max_v = array([0.1,0.,0.]) # maximum speed for intruder particle
        self.update_forces()
        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params()#self.nt)
        self.S = [Solid_Params(self.G,self)]
        self.S.append(Intruder_Params(self.G,self))
        self.supername = 'im/francois/ny_' + str(self.G.ny) + '/'
        self.pressure = 'lithostatic'
        print(self.supername)

    def update_forces(self):
        self.g = self.max_g
        t_c = 0.01
        self.intruder_v = self.max_v*(1.-exp(-3.*self.t**2/t_c**2))
        # self.intruder_v = self.max_v

class Grid_Params():
    def __init__(self,args):
        self.x_m = 0.0 # (m)
        self.x_M = 0.5 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 0.1 # (m)
        self.ny = int(args[1])
        self.nx = 5*(self.ny-1) + 1
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.R = 0.01 # radius of intruder
        self.x_c = (self.x_m + self.x_M)/2. # initial x location of intruder
        self.y_c = (self.y_m + self.y_M - self.dy)/2. # initial y location of intruder

        self.top_gap = self.dy # no material at top - HACK - NEED TO UPDATE self.y IF THIS CHANGES!!!

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.cyclic_lr = True

class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700. # density (kg/m^3)
        self.packing = 0.6 # packing fraction
        self.rho_s = self.rho/self.packing # solid density

        self.law = 'pouliquen'
        self.mu_0 = tan(deg2rad(20.9)) #0.3
        self.mu_1 = tan(deg2rad(32.76))
        self.delta_mu = self.mu_1 - self.mu_0
        self.I_0 = 0.279

        self.E = 1e4
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3.*(1.-2.*self.nu))
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)
        self.eta_max = 1.*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)

        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-2)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)

        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-gap[1]-G.top_gap,self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)

        for i in range(self.x*self.y):
            if (X[i]-G.x_c)**2 + (Y[i]-G.y_c)**2 > G.R**2:
                self.X.append(X[i])
                self.Y.append(Y[i])
                self.n += 1

        # self.A = ((G.x_M-G.x_m)*(G.y_M-G.y_m-G.top_gap) - pi*G.R**2)/self.n # area (m^2)
        self.A = G.dy*G.dx/self.pts_per_cell**2

class Intruder_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700.

        self.law = 'elastic'
        self.E = 1e4 # young's modulus (Pa)
        self.nu = 0.4 # poisson's ratio (-)
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-2)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)

        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-gap[1]-G.top_gap,self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)

        self.initial_v = P.intruder_v
        for i in range(self.x*self.y):
            if (X[i]-G.x_c)**2 + (Y[i]-G.y_c)**2 <= G.R**2:
                self.X.append(X[i])
                self.Y.append(Y[i])
                self.n += 1

        # self.A = pi*G.R**2/self.n # area (m^2)
        self.A = G.dy*G.dx/self.pts_per_cell**2

        # Advective terms
        elastic_wave_speed = sqrt(self.K/self.rho)
        max_speed = maximum(elastic_wave_speed,P.max_v[0])
        distance = minimum(G.dx,G.dy)
        critical_adv_time = distance/max_speed
        # Diffusive terms
        critical_diff_time = distance**2/P.S[0].eta_max

        critical_time = minimum(critical_adv_time,critical_diff_time)
        CFL = critical_time/P.dt

        if CFL < 2:
            print('WARNING: STABILITY IS GUARANTEED TO BE POOR')
            print(critical_adv_time)
            print(critical_diff_time)
        print('Current CFL: ' + str(CFL))

class Output_Params():
    def __init__(self):
        self.plot_continuum = True
        # self.plot_material_points = True
        self.continuum_fig_size = [22,4]
        self.mp_fig_size = [6,6]
