import os
from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp,linspace
from numpy import random, cos, sin, ceil, around
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-3 # timestep (s)
        self.savetime = 0.1
        self.t_f = 5.#100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.max_g = -10. # gravity (ms^-2)
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params()
        self.S = [Solid_Params()]

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

        self.E = 1.e6 # elastic modulus (Pa)
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

        for i in linspace(0,r,nr):
            dnphi = int(around(nphi*i/r)) # number in this ring
            if dnphi > 0:
                for j in linspace(0,(1.-1./(dnphi))*2*pi,dnphi):
                    self.X.append(c[0]+i*sin(j))
                    self.Y.append(c[1]+i*cos(j))
                    self.n += 1
        self.A = pi*r**2/self.n # area (m^2)


class Output_Params():
    def __init__(self):
        self.continuum_fig_size = [10,6]
        self.mp_fig_size = [10,10]

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
