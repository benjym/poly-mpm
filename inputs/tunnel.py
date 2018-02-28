import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,args):
        self.dt = 1e-2 # timestep (s)
        self.savetime = 10*self.dt
        self.t_f = 10.0 #100.0 # 3*self.dt # final time (s)
        self.max_g = -9.81 # gravity (ms^-2)
        self.theta = 0*pi/180. # slope angle (degrees)
        self.pressure = 'lithostatic'
        self.supername = 'im/tunnel/'
        self.damping = True
        
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self)
        self.S = [Solid_Params(self.G,self,'in'),Solid_Params(self.G,self,'out')] # in and out of tunnel

    def update_forces(self):
#         t_c = 0.05
#         self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.g = self.max_g

class Grid_Params():
    def __init__(self):
        self.x_m = -100 # (m)
        self.x_M = 100 # (m)
        self.y_m = -100 # (m)
        self.y_M = 0 # (m)
        self.nx = 51
        self.ny = 21
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.tunnel_centre = [0,-30] # m
        self.tunnel_radius = 10 # m

class Boundary_Params():
    def __init__(self):
        self.has_bottom = True
        self.cyclic_lr = True
        
class Solid_Params():
    def __init__(self,G,P,portion):
        self.X = []
        self.Y = []
        self.n = 0
        self.rho = 2700. # density (kg/m^3)

        self.law = 'elastic_time'
        self.E_0 = 1e7
        self.E_f = self.E_0/1e4
        self.rho_f = self.rho/1e2
        self.nu = 0.3 # poissons ratio
        self.t_0 = 1. # wait for consolidation before relaxation
        self.t_c = 1. # s, relaxation time for elasticity in tunnel

        self.pts_per_cell = 3
        self.x = (G.nx-1)*self.pts_per_cell # particles in x direction
        self.y = (G.ny-1)*self.pts_per_cell # particles in y direction
        gap = array((G.dx,G.dy))/(2*self.pts_per_cell)
        xp = linspace(G.x_m+gap[0],G.x_M-gap[0],self.x)
        yp = linspace(G.y_m+gap[1],G.y_M-gap[1],self.y)
        X = tile(xp,self.y)
        Y = repeat(yp,self.x)
        n_tot = 0
        for i in range(self.x*self.y):
            if portion == 'out':
                if (X[i]-G.tunnel_centre[0])**2 + (Y[i]-G.tunnel_centre[1])**2 > G.tunnel_radius**2:
                    self.X.append(X[i])
                    self.Y.append(Y[i])
                    self.n += 1
            elif portion == 'in':
                if (X[i]-G.tunnel_centre[0])**2 + (Y[i]-G.tunnel_centre[1])**2 < G.tunnel_radius**2:
                    self.X.append(X[i])
                    self.Y.append(Y[i])
                    self.n += 1
            n_tot += 1
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/n_tot # area (m^2)

#         elastic_wave_speed = sqrt(self.K/(self.A*self.rho*1.))
        elastic_wave_speed = sqrt(self.E_0/(3*(1-2*self.nu))/self.rho)
        distance = minimum(G.dx,G.dy)
        critical_time = distance/elastic_wave_speed
        critical_time /= 2. # safety
        if critical_time < P.dt:
#             P.dt = critical_time
            print('WARNING: timestep should be dt = ' + str(critical_time))

class Output_Params():
    def __init__(self,P):
        self.plot_continuum = True
        self.continuum_fig_size = [10,4]
#         self.plot_material_points = True
#         self.mp_fig_size = [10,8]