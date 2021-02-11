import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-8 # timestep (s)
        self.savetime = 0.1 # 1e3*self.dt # 0.1 (s)
        self.t_f = 100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.max_q = 0.
        self.theta = 0. # world tilt angle (degrees)
        self.slope_angle = deg2rad(45.0)
        self.G = Grid_Params(args)
        self.B = Boundary_Params()
        self.O = Output_Params()#self.nt)
        self.S = [Solid_Params(self.G,self,args),]
        self.segregate_grid = True
        self.c = 1e-3 # inter-particle drag coefficient
        # self.D = 1e-3 # segregation diffusion coefficient
        self.l = 10. # number of particle diameters for seg diffusion coeff
        self.supername = 'im/PSM_slope_stab/ny_' + args[1] + '/ns_' + args[3] + '/' + args[2] + '/c_' + str(self.c) + '/l_' + str(self.l) + '/'
        self.pressure = 'lithostatic'
        self.smooth_gamma_dot = False # smooth calculation of gamma_dot
        self.normalise_phi = True
        self.time_stepping = 'dynamic' # dynamic or static time steps
        self.CFL = 0.2 # stability criteria for determining timstep
        print(self.supername)

    def update_forces(self):
        t_c = 0.05
        self.g = self.max_g

class Grid_Params():
    def __init__(self,args):
        self.L = 2 # aspect ratio of box
        self.y_m = 0.0 # (m)
        self.y_M = 10.0 # (m)
        self.x_m = 0.0 # (m)
        self.L_1 = 1.0*self.y_M # length of bottom of pile
        self.x_M = 2*self.L_1 # (m)
        self.ny = int(args[1])
        self.nx = (self.ny-1)*self.L + 1
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        # self.s = array([0.5,1.0]) # s coordinate

        self.R = 10.
        self.s_M = 0.1 # 10 cm aggregate
        self.s_m = self.s_M/self.R
        self.ns = int(args[3])
        # self.s = array([self.s_m,self.s_M])
        s_edges = linspace(self.s_m,self.s_M,self.ns+1)
        self.s = (s_edges[1:] + s_edges[:-1])/2.
        # self.s = linspace(self.s_m,self.s_M,self.ns)
        self.ds = self.s[1]-self.s[0]

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.has_left = True
        # self.outlet_left = True
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

        self.law = 'pouliquen'
        self.mu_0 = 2.*tan(deg2rad(20.90)) # double glass beads - about 37.5 deg
        self.mu_1 = 2.*tan(deg2rad(32.76)) # double glass beads
        self.delta_mu = self.mu_1 - self.mu_0
        self.I_0 = 0.279
        self.eta_max = 1e2*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)/1e2

        self.E = 1e7
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3*(1-2*self.nu))
        self.G = self.E/(2*(1+self.nu))

        self.pts_per_cell = 3
        self.ny = (G.ny-1)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)

        yp = linspace(G.y_m + gap[1],G.y_M - gap[1],self.ny)
        for j in range(self.ny):
            current_L = G.L_1 - yp[j]/tan(P.slope_angle) # length of slope at this height
            nx = int(current_L/G.dx*self.pts_per_cell) # particles in x direction
            if nx > 0:
                xp = linspace(G.x_m + gap[0],current_L - gap[0],nx)
                for i in range(nx):
                    self.X.append(xp[i])
                    self.Y.append(yp[j])
                    if args[2] == 'top': # small on top
                        this_phi = zeros([G.ns])
                        H = (G.y_M - G.y_m)
                        this_phi_arg = int(floor((H-yp[j])/H*G.ns))
                        this_phi[this_phi_arg] = 1.
                        self.PHI.append(this_phi)
                    elif args[2] == 'bot': # large on top
                        this_phi = zeros([G.ns])
                        H = (G.y_M - G.y_m)
                        this_phi_arg = int(floor(yp[j]/H*G.ns))
                        this_phi[this_phi_arg] = 1.
                        self.PHI.append(this_phi)
                    elif args[2] == 'mix': # mixed
                        self.PHI.append(ones([G.ns])/float(G.ns))
                    elif args[2] == 'layered':
                        this_phi = zeros([G.ns])
                        dist_from_origin = (yp[j] + xp[i]*tan(P.slope_angle))/sqrt(tan(P.slope_angle)**2 + 1)
                        layerwidth = 1
                        this_phi_arg = int(G.ns/2*cos(dist_from_origin/layerwidth*pi) + G.ns/2)
                        this_phi[this_phi_arg] = 1.
                        self.PHI.append(this_phi)
                    self.n += 1
        self.A = G.dx*G.dy/self.pts_per_cell**2 # FIXME

    def critical_time(self,P):
        distance = minimum(P.G.dx,P.G.dy)
        t_ela = distance/sqrt(self.K/self.rho) # elasticity
        t_diff = distance**2/self.eta_max*self.rho # momentum diffusivity/viscosity
        return minimum(t_diff,t_ela)



class Output_Params():
    def __init__(self):
        self.continuum_fig_size = [24,8]
        self.mp_fig_size = [18,4]

    def after_every_nth_timestep(self,P,G,L,plot):
        # plot.draw_gsd_mp(L,P,G)
        # plot.draw_gsd_grid(L,P,G)
        plot.draw_continuum(G,P)
        # plot.draw_material_points(L,P,G)
        plot.save_u(L,P,G)
        plot.save_s_bar(L,P,G)
        plot.save_density(L,P,G)
        plot.save_phi_MP(L,P,G)

    def final_graphs(self,P,G,L,plot):
        # plot.draw_material_points(L,P,G,'final')
        # plot.draw_gsd_mp(L,P,G)
        # plot.draw_gsd_grid(L,P,G)
        plot.save_u(L,P,G)
        plot.save_s_bar(L,P,G)
        plot.save_density(L,P,G)
        plot.save_phi_MP(L,P,G)
