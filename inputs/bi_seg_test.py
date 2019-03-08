import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-5#5e-8 # timestep (s)
        self.savetime = 1e-2 # (s)
        self.t_f = 10.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.theta = -deg2rad(float(args[1])) # slope angle (degrees)
        self.pressure = 'lithostatic'
        self.initial_flow = 'steady'

        self.segregate_grid = False
        self.c = 1e-4 # inter-particle drag coefficient
        self.D = 0.#1e-2 # segregation diffusion coefficient

        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params(self)#self.nt)
        self.S = [Solid_Params(self.G,self),]

        self.supername = 'im/seg/lin/theta_' + str(-rad2deg(self.theta)) + '/ns_' + str(self.G.ns) + '/'
        self.smooth_gamma_dot = True
        self.smooth_grad2 = True
        print('Saving to ' + self.supername)

    def update_forces(self):
        self.g = self.max_g

class Grid_Params():
    def __init__(self,args):
        self.ny = 31
        self.y_m = 0.0 # (m)
        self.y_M = 0.1 # (m)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

        self.nx = 2
        self.x_m = 0.0 # (m)
        self.x_M = self.dy*(self.nx-1) + self.x_m # (m)
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)

        self.s_m = 0.0001 # 100 micron
        self.s_M = 0.001 # 1 mm
        self.ns = int(args[2])
        s_edges = linspace(self.s_m,self.s_M,self.ns+1)
        self.s = (s_edges[1:] + s_edges[:-1])/2.
        self.ds = self.s[1]-self.s[0]

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.cyclic_lr = True
        self.no_slip_bottom = True

class Solid_Params():
    def __init__(self,G,P):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700. # density (kg/m^3)
        self.packing = 0.6 # packing fraction
        self.rho_s = self.rho/self.packing # solid density
        self.phi = ones([G.ns])/float(G.ns)
        self.law = 'pouliquen'
        self.mu_0 = tan(deg2rad(20.9)) #0.3
        self.mu_1 = tan(deg2rad(32.76))
        self.delta_mu = self.mu_1 - self.mu_0
        self.I_0 = 0.279

        self.mu_v = 0
        self.E = 1e8
        self.nu = 0.3 # poissons ratio
        self.K = self.E/(3.*(1.-2.*self.nu))
        self.G = self.E/(2.*(1.+self.nu))
        self.eta_max = 100.*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3) # taken a paper somewhere

        G.s_bar_0 = (G.s_m+G.s_M)/2.
        self.v_max = ( sqrt(abs(P.max_g)*G.s_bar_0)*(2./3.)*self.I_0*(tan(abs(P.theta))-self.mu_0)/(self.mu_1-tan(abs(P.theta)))*
                  sqrt(self.packing*cos(abs(P.theta)))*G.y_M**1.5/G.s_bar_0**1.5 ) # Andreotti, Pouliquen and Forterre, page 233, equation (6.17)

        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-1)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
#         self.A = self.conc*(G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
        # self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)
        self.A = G.dx*G.dy/self.pts_per_cell**2

        # Advective terms
        elastic_wave_speed = sqrt(self.K/self.rho)
        distance = minimum(G.dx,G.dy)
        critical_adv_time = distance/elastic_wave_speed
        # Diffusive terms
        l_0 = xp[1] - xp[0] # initial distance between material points
        critical_diff_time = l_0**2*self.rho/self.eta_max

        critical_time = minimum(critical_adv_time,critical_diff_time)
        CFL = critical_time/P.dt

        if CFL < 2:
            print('WARNING: STABILITY IS GUARANTEED TO BE POOR')
        print('CFL from elastic wave speed: ' + str(critical_adv_time/P.dt))
        print('CFL from momentum diffusion: ' + str(critical_diff_time/P.dt))
        print('Current CFL: ' + str(CFL))

class Output_Params():
    def __init__(self,P):
        self.plot_continuum = True
        # self.plot_material_points = True
        # self.plot_gsd_grid = True
        self.save_u = True
        self.save_s_bar = True
        self.continuum_fig_size = [12,8]
        self.mp_fig_size = [6,6]
