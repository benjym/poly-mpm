import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-4 # timestep (s)
        self.savetime = 1e-1 #1e1*self.dt#0.01
        self.t_f = 100.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.max_q = 0.
        self.theta = 0*pi/180. # slope angle (radians)
        self.Fr = float(args[2]) # Froude number
        self.r = 1.0 # radius of drum
        self.G = Grid_Params(self,args)
        self.B = Boundary_Params()
        self.O = Output_Params(self.G)#self.nt)
        self.S = [Solid_Params(self.G,self,args),]
        self.segregate_grid = True
        self.c = 1e-3 # inter-particle drag coefficient
        self.l = 10. # number of particle diameters for seg diffusion coeff
        # self.D = 1e-3 # segregation diffusion coefficient
        self.supername = 'im/drum/ny_' + str(self.G.ny) + '/Fr_' + str(self.Fr) + '/'
        # self.supername = 'im/drum_rough/wall_mu_' + str(self.B.wall_mu) + '/ny_' + str(self.G.ny) + '/Fr_' + str(self.Fr) + '/c_' +
        # self.supername = 'im/drum_no_slip/ny_' + str(self.G.ny) + '/Fr_' + str(self.Fr) + '/c_' + str(self.c) + '/l_' + str(self.l) + '/'
        self.pressure = 'lithostatic'
        # self.smooth_gamma_dot = True # smooth calculation of gamma_dot
        self.time_stepping = 'dynamic' # dynamic or static time steps
        self.normalise_phi = True
        self.CFL = 0.2 # stability criteria for determining timstep
        print(self.supername)

    def update_forces(self):
        t_c = sqrt(4*pi**2*self.r/(abs(self.max_g)*self.Fr)) # rotation period
        self.theta = -(self.t/t_c)*2.*pi # slope angle (radians)
        self.g = self.max_g

class Grid_Params():
    def __init__(self,P,args):
        self.y_m = 0.0 # (m)
        self.y_M = P.r # (m)
        self.x_m = 0.0 # (m)
        self.x_M = P.r # (m)
        self.ny = int(args[1])
        self.nx = self.ny
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        # self.s = array([0.5,1.0]) # s coordinate
        self.top_gap = 0.5*(self.y_M - self.y_m)

        self.R = 10.#float(args[2])
        self.s_M = 0.003 # 3mm beads, following https://journals.aps.org/pre/pdf/10.1103/PhysRevE.62.961
        # self.s_m = 0.1
        self.s_m = self.s_M/self.R
        self.ns = 2
        self.s = array([self.s_m,self.s_M])
        # s_edges = linspace(self.s_m,self.s_M,self.ns+1)
        # self.s = (s_edges[1:] + s_edges[:-1])/2.
        self.ds = self.s[1]-self.s[0]

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.has_top = True
        self.has_right = True
        self.has_left = True
        self.roughness = False

        # self.roughness = True

        # self.wall_mu = 0.1

class Solid_Params():
    def __init__(self,G,P,args):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700. # density (kg/m^3)
        self.packing = 0.6 # packing fraction
        self.rho_s = self.rho/self.packing # solid density
        self.PHI = []

        self.law = 'pouliquen'
        # self.law = 'ken_simple'
        self.mu_0 = tan(deg2rad(20.9)) #0.3
        self.mu_1 = tan(deg2rad(32.76))
        self.delta_mu = self.mu_1 - self.mu_0
        self.I_0 = 0.279
        self.eta_max = 100.*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)/10.0#/1e2

        self.E = 1e7
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3*(1-2*self.nu))
        self.G = self.E/(2*(1+self.nu))

        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-1)//2*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)

        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-gap[1] - G.top_gap,self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = (G.y_M - G.y_m - G.top_gap)*(G.x_M - G.x_m)/self.n # area (m^2)

    def critical_time(self,P):
        distance = minimum(P.G.dx,P.G.dy)
        t_ela = distance/sqrt(self.K/self.rho) # elasticity
        t_diff = distance**2/self.eta_max*self.rho # momentum diffusivity/viscosity
        return minimum(t_diff,t_ela)

class Output_Params():
    def __init__(self,G):
        self.plot_continuum = True
        # self.plot_material_points = True
        self.plot_gsd_mp = True
        # self.plot_gsd_grid = True
        self.save_s_bar = True
        self.save_u = True
        self.save_density = True
        self.save_phi_MP = True
        self.continuum_fig_size = [24,8]
        self.mp_fig_size = [18,4]
        # self.save_fields = [[G.s_bar,'s_bar'],[G.m/G.V,'density'],[G.q[:,0]/G.m,'u'],[G.q[:,1]/G.m],'v']
