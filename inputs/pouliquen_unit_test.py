import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-5 # timestep (s)
        self.savetime = 0.01
        self.t_f = 10.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.theta = -deg2rad(float(args[1])) # slope angle (degrees)
        self.pressure = 'lithostatic'
        self.initial_flow = 'steady'

        self.segregate_grid = False
        self.c = 0. # inter-particle drag coefficient
        self.D = 0. # segregation diffusion coefficient

        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params(self)#self.nt)
        self.S = [Solid_Params(self.G,self),]

        self.supername = 'im/pouliquen_unit_test/theta_' + str(-rad2deg(self.theta)) + '/ny_' + str(self.G.ny) + '/'
        print('Saving to ' + self.supername)

    def update_forces(self):
        self.g = self.max_g

class Grid_Params():
    def __init__(self,args):
        self.s = [0.001]
        self.s_bar_0 = self.s[0]

        self.x_m = 0.0 # (m)
        self.x_M = 0.02 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 50 * self.s_bar_0 #0.1 # (m)
        self.nx = 2
        self.ny = int(args[2])
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)


class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.cyclic_lr = True
        self.roughness = True

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

        self.mu_v = 0
        self.E = 1e4
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3*(1-2*self.nu))
        self.G = self.E/(2*(1+self.nu))
        self.eta_max = 1.*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)

        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-1)*self.pts_per_cell # particles in y direction
        G.gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        xp = linspace(G.x_m+G.gap[0],G.x_M-G.gap[0],self.x)
        yp = linspace(G.y_m+G.gap[1],G.y_M-G.gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = G.dx*G.dy/self.pts_per_cell**2

        self.v_max = ( sqrt(abs(P.max_g)*G.s_bar_0)*(2./3.)*self.I_0*(tan(abs(P.theta))-self.mu_0)/(self.mu_1-tan(abs(P.theta)))*
                  sqrt(self.packing*cos(abs(P.theta)))*G.y_M**1.5/G.s_bar_0**1.5 ) # Andreotti, Pouliquen and Forterre, page 233, equation (6.17)
        # print('Max velocity should be ' + str(self.v_max) + ' m/s')

        # Advective terms
        elastic_wave_speed = sqrt(self.K/self.rho)
        max_speed = maximum(elastic_wave_speed,self.v_max)
        distance = minimum(G.dx,G.dy)
        critical_adv_time = distance/max_speed
        # Diffusive terms
        critical_diff_time = distance**2/self.eta_max

        critical_time = minimum(critical_adv_time,critical_diff_time)
        CFL = critical_time/P.dt

        if CFL < 2:
            print('WARNING: STABILITY IS GUARANTEED TO BE POOR')
            print(critical_adv_time)
            print(critical_diff_time)
        print('Current CFL: ' + str(CFL))


class Output_Params():
    def __init__(self,P):
        self.plot_continuum = True
        # self.plot_material_points = True
        # self.plot_gsd_grid = True
        self.continuum_fig_size = [12,8]
        self.mp_fig_size = [6,6]
