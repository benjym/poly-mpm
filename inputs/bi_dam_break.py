import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-5 # timestep (s)
        self.savetime = 0.01 #1e1*self.dt#0.01
        self.t_f = 1.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.max_q = 0.
        self.theta = float(args[4])*pi/180. # slope angle (degrees)
        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params()#self.nt)
        self.S = [Solid_Params(self.G,self,args),]
        self.segregate_grid = True
        self.c = 1e-3 # inter-particle drag coefficient
        self.D = 0. # segregation diffusion coefficient
        self.supername = 'im/dam_break/ny_' + str(self.G.ny) + '/ns_' + str(self.G.ns) + '/R_' + str(self.G.R) + '/' + args[3] + '/theta_' + args[4] + '/'
        self.pressure = 'lithostatic'
        self.smooth_gamma_dot = False # smooth calculation of gamma_dot
        self.smooth_grad2 = True # smooth the gradient of the shear strain rate
        print(self.supername)

    def update_forces(self):
        t_c = 0.05
        self.g = self.max_g

class Grid_Params():
    def __init__(self,args):
        self.L = 5 # aspect ratio of box
        self.L_1 = 2.0 # aspect ratio of initial pile
        self.y_m = 0.0 # (m)
        self.y_M = 5.0 # (m)
        self.x_m = 0.0 # (m)
        self.x_M = self.y_M*self.L # (m)
        self.ny = int(args[1])
        self.nx = (self.ny-1)*self.L + 1
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        # self.s = array([0.5,1.0]) # s coordinate

        self.R = float(args[2])
        self.s_M = 1.0 # 1m grains! 3mm beads, following https://journals.aps.org/pre/pdf/10.1103/PhysRevE.62.961
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
        self.has_right = True
        self.outlet_left = True
        # self.roughness = True
        self.no_slip_bottom = True

class Solid_Params():
    def __init__(self,G,P,args):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700. # density (kg/m^3)
        self.packing = 0.6 # packing fraction
        self.rho_s = self.rho/self.packing # solid density
        self.heterogeneous = True # non-uniform gsd in space
        self.PHI = []

        # self.law = 'pouliquen'
        # self.mu_0 = tan(deg2rad(20.9)) #0.3
        # self.mu_1 = tan(deg2rad(32.76))
        # self.delta_mu = self.mu_1 - self.mu_0
        # self.I_0 = 0.279
        self.eta_max = 100.*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)

        self.law = 'linear_mu'
        self.mu_0 = 0.5
        self.b = 0.5

        self.E = 1e6
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3*(1-2*self.nu))
        self.G = self.E/(2*(1+self.nu))

        self.pts_per_cell = 3
        self.x = int((G.nx-1)//G.L/G.L_1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-1)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)

        xp = linspace(G.x_m + (G.L - 1. + (G.L_1 - 1.)/G.L_1)/G.L*(G.x_M - G.x_m) + gap[0],G.x_M - gap[0],self.x)
        yp = linspace(G.y_m + gap[1],G.y_M - gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)

        for i in range(self.x*self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            if args[3] == 'top': # small on top
                if Y[i] > (G.y_M - G.y_m)/2: self.PHI.append([1.,0.])
                else:                        self.PHI.append([0.,1.])
            elif args[3] == 'bot': # large on top
                if Y[i] > (G.y_M - G.y_m)/2: self.PHI.append([0.,1.])
                else:                        self.PHI.append([1.,0.])
            elif args[3] == 'mix': # mixed
                self.PHI.append(ones([G.ns])/float(G.ns))

            self.n += 1
        self.A = G.dx*G.dy/self.pts_per_cell**2


class Output_Params():
    def __init__(self):
        self.plot_continuum = True
        # self.plot_material_points = True
        # self.plot_gsd_grid = True
        self.save_s_bar = True
        self.save_u = True
        self.continuum_fig_size = [24,8]
        self.mp_fig_size = [18,4]
