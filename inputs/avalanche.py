import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-4 # timestep (s)
        self.savetime = 0.1
        self.t_f = 1000.0 # 100*self.dt # final time (s)
        self.max_g = -10. # gravity (ms^-2)
        self.max_q = 0.
        self.theta = 45.*pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()#self.nt)
        self.S = [Solid_Params(self.G)]

    def update_forces(self):
        t_c = 0.5
        self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
#         self.g = self.max_g
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 4.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = 21
        self.ny = 6
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.thickness = 1. # (m)

class Boundary_Params():
    def __init__(self):
        self.wall = True
        self.has_bottom = True
        self.no_slip_bottom = True
        self.cyclic_lr = True

class Solid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 1000. # density (kg/m^3)

#         self.law = 'pouliquen'
#         self.mu_s = 0.32
#         self.mu_2 = 0.6
#         self.I_0 = 0.4
#         self.K = 1e6
#         self.mu_v = 0.

        self.law = 'viscous'
        self.mu_s = 1e3
        self.K = 1e6
        self.mu_v = 0

#         self.law = 'bingham'
#         self.mu_s = 1e2 # shear viscosity
#         self.tau_0 = 100. # yield stress
#         self.mu_0 = 1e3*self.mu_s # much steeper viscosity to mimic very rigid part
#         self.gamma_c = self.tau_0/(self.mu_0 - self.mu_s)
#         self.K = 1e7 # bulk modulus
#         self.mu_v = 0. #1e3

        self.inlet_rate = 0.1
        filling_fraction = 0.5
        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = int(round((G.ny-1)*self.pts_per_cell*filling_fraction)) # particles in y direction
        self.gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        self.xp = linspace(3.*(G.x_M-G.x_m)/4.,G.x_M-self.gap[0],self.x)
        self.yp = linspace(G.y_m+self.gap[1],self.gap[1]+filling_fraction*(G.y_M-2*self.gap[1]),self.y)
        X = tile(self.xp,self.y)
        Y = repeat(self.yp,self.x)
        for i in range(int(self.x*self.y)):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n*filling_fraction # area (m^2)
#         self.A = 0.002
#         self.initial_v = array([-0.0001,0,0])

class Output_Params():
    def __init__(self):
        self.continuum_fig_size = [18,6]
        self.mp_fig_size = [6,6]

        def after_every_nth_timestep(self,P,G,L,plot):
        plot.draw_continuum(G,P)
        plot.draw_material_points(L,P,G)
        # plot.draw_gsd_mp(L,P,G)
        # plot.draw_gsd_grid(L,P,G)
        # plot.draw_continuum(G,P)
        # plot.draw_material_points(L,P,G)
        # plot.draw_gamma_dot(L,P,G)
        # P.O.measure_E(L,P,G)
        # plot.save_u(L,P,G)
        # plot.save_s_bar(L,P,G)
        # plot.save_density(L,P,G)
        # plot.save_phi_MP(L,P,G)

    def final_graphs(self,P,G,L,plot):
        plot.draw_continuum(G,P)
        plot.draw_material_points(L,P,G,'final')
        # P.O.measure_E(L,P,G)
        # plot.draw_energy(P)
        # plot.draw_gsd_mp(L,P,G)
        # plot.draw_gsd_grid(L,P,G)
        # plot.save_u(L,P,G)
        # plot.save_s_bar(L,P,G)
        # plot.save_density(L,P,G)
        # plot.save_phi_MP(L,P,G)
        # P.O.draw_p_q(P,G,L,plot,P.tstep)
        # P.O.draw_mu(P,G,L,plot,P.tstep)
