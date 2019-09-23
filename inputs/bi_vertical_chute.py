import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-8#5e-8 # timestep (s)
        self.savetime = 1e-3 # (s)
        self.t_f = 5.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.theta = -deg2rad(90.) # slope angle (degrees)
        self.pressure = 'compression'
        self.initial_pressure = 7.0 # 7 Pa - same as Fan and Hill (ish)
        # self.initial_flow = 'steady'

        self.segregate_grid = True
        self.c = 1e-2 # inter-particle drag coefficient
        self.D = 0.#1e-2 # segregation diffusion coefficient

        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params(self)#self.nt)
        self.S = [Solid_Params(self.G,self),]

        self.supername = 'im/vert_chute/lin/ny_' + str(self.G.ny) + '/ns_' + str(self.G.ns) + '/'
        self.smooth_gamma_dot = False
        self.time_stepping = 'dynamic' # dynamic or static time steps
        self.CFL = 0.2

        print('Saving to ' + self.supername)

    def update_forces(self):
        self.g = self.max_g

class Grid_Params():
    def __init__(self,args):
        self.ny = int(args[2])
        self.y_m = 0.0 # (m)
        self.y_M = 0.05 # (m)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

        self.nx = 2
        self.x_m = 0.0 # (m)
        self.x_M = 1.0#self.dy*(self.nx-1) + self.x_m # (m)
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)

        # self.R = 10.#float(args[3])
        self.s_M = 3e-3 # 3 mm
        self.s_m = 2e-3 #self.s_M/self.R
        self.ns = int(args[1])
        s_edges = linspace(self.s_m,self.s_M,self.ns+1)
        self.s = (s_edges[1:] + s_edges[:-1])/2.
        self.ds = self.s[1]-self.s[0]
        self.s_bar_0 = (self.s_m+self.s_M)/2.

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.has_top = True
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
        self.eta_max = 100.*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)/5e1

        self.E = 1e7
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3*(1-2*self.nu))
        self.G = self.E/(2*(1+self.nu))

        # self.law = constit.pouliquen(E,nu,mu_0,mu_1,I_0,eta_max) # DO I WANT TO MAKE A GENERIC LAW CLASS AND HAVE EACH CONSTITUTIVE MODEL INHERIT THOSE PROPERTIES?!??
        # self.v_max = constit.pouliquen.get_analytical_vmax(P,G)
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

        # # Advective terms
        # elastic_wave_speed = sqrt(self.K/self.rho)
        # distance = minimum(G.dx,G.dy)
        # critical_adv_time = distance/elastic_wave_speed
        # # Diffusive terms
        # l_0 = xp[1] - xp[0] # initial distance between material points
        # critical_diff_time = l_0**2*self.rho/self.eta_max
        #
        # critical_time = minimum(critical_adv_time,critical_diff_time)
        # CFL = critical_time/P.dt
        #
        # if CFL < 2:
        #     print('WARNING: STABILITY IS GUARANTEED TO BE POOR')
        # print('CFL from elastic wave speed: ' + str(critical_adv_time/P.dt))
        # print('CFL from momentum diffusion: ' + str(critical_diff_time/P.dt))
        # print('Current CFL: ' + str(CFL))

    def critical_time(self,P):
        distance = minimum(P.G.dx,P.G.dy)
        t_ela = distance/sqrt(self.K/self.rho) # elasticity
        t_diff = distance**2/self.eta_max*self.rho # momentum diffusivity/viscosity
        return minimum(t_diff,t_ela)


class Output_Params():
    def __init__(self,P):
        self.plot_continuum = True
        # self.plot_material_points = True
        # self.plot_gsd_grid = True
        # self.plot_gsd_mp = True
        # self.plot_gsd_grid = True
        self.save_s_bar = True
        self.save_u = True
        self.save_density = True
        self.save_phi_MP = True
        self.continuum_fig_size = [12,8]
        self.mp_fig_size = [6,6]
