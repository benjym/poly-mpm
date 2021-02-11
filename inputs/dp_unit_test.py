import os
from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp,linspace,minimum,amax,tile,repeat
from numpy import random,ones
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-4 # timestep (s)
        self.savetime = 1
        self.t_f = 1 # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.max_q = 1.
        self.max_g = 0.
        self.pressure = 'non-zero'
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.S = [Solid_Params(self.G),]
        self.O = Output_Params(self.nt,self.S[0].n)
        self.supername = 'im/dp_unit_test'
#         self.damping = True

    def update_forces(self):
        self.g = 0
        t_c = self.t_f/10.
        if self.t < t_c:
            self.q_h = self.q_v = self.max_q*self.t/t_c
        else:
            self.q_h = self.t*self.max_q/t_c
            self.q_v = 2*self.max_q - self.t/t_c*self.max_q
#         self.q_h = self.t*self.max_q/t_c
#         self.q_v = minimum(self.q_h,self.max_q)

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 1.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 1.0 # (m)
        self.nx = 2
        self.ny = 2
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

class Boundary_Params():
    def __init__(self):
#         self.has_right = True
#         self.has_left = True
        self.vertical_force = True
        self.horizontal_force = True

class Solid_Params():
    def __init__(self,G):
        self.X = []
        self.Y = []
        self.x = 3#(10*2-1)*2 # particles in x direction
        self.y = 3#(5*2-1)*2-4 # particles in y direction
        self.n = 0

        self.law = 'dp'
#         self.law = 'rigid'
        self.rho = 2650. # density (kg/m^3)
        self.packing = 0.6 # packing fraction
        self.rho_s = self.rho/self.packing # solid density

        self.E = 1e7 # elastic modulus (Pa)
        self.nu = 0.4 # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        # self.mu = 1.
        self.mu_0 = 0.5
        self.mu_1 = 1.
        self.delta_mu = self.mu_1 - self.mu_0
        self.I_0 = 1e-3
        self.beta = 0 #1.*self.mu
        self.s = 2.
        self.eta_max = 1e10

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
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)


class Output_Params():
    def __init__(self,nt,n):
#         self.plot_continuum = True
        self.plot_material_points = True
        self.energy = zeros((nt*10+1,4)) # energy
        self.p = zeros((nt*10+10,n))
        self.q = zeros((nt*10+10,n))

    def store_p_q(self,P,G,L,tstep):
        for i in range(P.S[0].n):
            self.p[tstep,i] = L.S[0][i].p
            self.q[tstep,i] = L.S[0][i].q

    def draw_p_q(self,P,G,L,plot,tstep):
#         x = zeros((P.S[0].n,3))
#         v = zeros((P.S[0].n,3))
        plt.figure()
        for i in range(P.S[0].n):
            plt.clf()
            plt.xlabel(r"$p$")
            plt.ylabel(r"$q$",rotation='horizontal')
            plt.plot([0,amax(self.p)],[0,P.S[0].mu_0*amax(self.p)],'r--')
            plt.plot([0,amax(self.p)],[0,P.S[0].mu_1*amax(self.p)],'g--')
            plt.plot(self.p[:tstep-1,i],self.q[:tstep-1,i],'b-')
            plot.savefig(P,str(i))
        #print (self.q[tstep-1,0]-self.q[1,0])/(self.p[tstep-1,0]-self.p[1,0])
