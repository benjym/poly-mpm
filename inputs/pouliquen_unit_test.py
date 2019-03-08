import os
from numpy import *
import matplotlib.pyplot as plt

class Params():
    def __init__(self,mode):
        self.dt = 1e-4 # timestep (s)
        self.savetime = 0.1
        self.t_f = 0.1 # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.max_q = 1e2
        self.max_g = 0.
        self.pressure = 'non-zero'
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.S = Solid_Params(self,self.G)
        self.O = Output_Params(self.nt,self.S.n)
        self.supername = 'im/pouliquen_unit_test'
#         self.damping = True

    def update_forces(self):
        self.g = 0
        t_c = self.t_f/2.
        if self.t < t_c:
            self.q_h = self.t*self.max_q/t_c
            self.q_v = self.t*self.max_q/t_c
        else:
            self.q_h = self.max_q + (self.t - t_c)*self.max_q/t_c #+ 1e2*self.max_q*((self.t-t_c))#**3#/(self.t_f-t_c))**3
            self.q_v = self.max_q - (self.t - t_c)*self.max_q/t_c
            # self.q_v = 2*self.max_q - self.t/t_c*self.max_q
#         self.q_h = self.t*self.max_q/t_c
#         self.q_v = minimum(self.q_h,self.max_q)

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 0.1 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 0.1 # (m)
        self.nx = 2
        self.ny = 2
        self.x = linspace(self.x_m,self.x_M,self.nx)
        self.y = linspace(self.y_m,self.y_M,self.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)

        self.ns = 1
        self.s = [0.01] # (m)

class Boundary_Params():
    def __init__(self):
#         self.has_right = True
#         self.has_left = True
        self.vertical_force = True
        self.horizontal_force = True

class Solid_Params():
    def __init__(self,P,G):
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
        self.nu = 0.4 # poissons ratio
        self.K = self.E/(3.*(1.-2.*self.nu))
        self.G = self.E/(2.*(1.+self.nu))
        self.eta_max = 1e12#*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)
        # self.eta_max = 1e2*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)

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
        self.A = (G.x_M-G.x_m)*(G.y_M-G.y_m)/self.n # area (m^2)

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
    def __init__(self,nt,n):
        # self.plot_continuum = True
        # self.plot_material_points = True
        self.energy = zeros([nt+1,4]) # energy
        self.mu = zeros([nt+2,n])
        self.I = zeros_like(self.mu)
        self.tau = zeros_like(self.mu)
        self.eta = zeros_like(self.mu)
        self.t = zeros_like(self.mu)
        self.pressure = zeros_like(self.mu)
        self.gammadot = zeros_like(self.mu)
        self.dstrain = zeros([nt+2,n,3,3])
        self.dev_stress = zeros_like(self.dstrain)

    def store_mu(self,P,G,L,tstep):
        if tstep > 1:
            for i in range(P.S[0].n):
                self.mu[tstep,i] = L.S[0][i].mu
                self.dev_stress[tstep,i] = L.S[0][i].dev_stress
                self.I[tstep,i] = L.S[0][i].I
                self.eta[tstep,i] = L.S[0][i].eta
                self.t[tstep,i] = P.t
                self.pressure[tstep,i] = L.S[0][i].pressure
                self.dstrain[tstep,i] = L.S[0][i].dstrain
                self.gammadot[tstep,i] = L.S[0][i].gammadot
                self.tau[tstep,i] = linalg.norm(L.S[0][i].dev_stress)/sqrt(2)

    def draw_mu(self,P,G,L,plot,tstep):
        plt.figure(figsize=[12,6])
        # for i in range(P.S[0].n):
        for i in [0]:
            plt.clf()

            plt.subplot(231)
            plt.plot(self.t[:tstep-1],self.I[:tstep-1,i],'k.')
            plt.xlabel(r'$t$',rotation='horizontal')
            plt.ylabel(r'$I$',rotation='horizontal')

            plt.subplot(232)
            plt.plot(self.t[:tstep-1,i],-self.pressure[:tstep-1,i],'k.')
            plt.plot(self.t[:tstep-1,i],self.tau[:tstep-1,i],'b.')
            plt.xlabel(r'$t$',rotation='horizontal')
            plt.ylabel(r'$\sigma_{ij}$',rotation='horizontal')

            plt.subplot(233)
            plt.semilogx(self.I[:tstep-1,i],self.mu[:tstep-1,i],'k.')
            plt.semilogx(self.I[:tstep-1,i],-self.tau[:tstep-1,i]/self.pressure[:tstep-1,i],'b.')
            plt.xlabel(r'$I$',rotation='horizontal')
            plt.ylabel(r'$\mu$',rotation='horizontal')

            plt.subplot(234)
            plt.plot(self.t[:tstep-1],self.eta[:tstep-1,i],'k.')
            plt.xlabel(r'$t$',rotation='horizontal')
            plt.ylabel(r'$\eta$',rotation='horizontal')

            plt.subplot(235)
            plt.plot(self.t[:tstep-1],self.dstrain[:tstep-1,i,0,0]/P.dt,'r.',label=r'$\dot\varepsilon_{xx}$')
            plt.plot(self.t[:tstep-1],self.dstrain[:tstep-1,i,0,1]/P.dt,'g.',label=r'$\dot\varepsilon_{xy}$')
            plt.plot(self.t[:tstep-1],self.dstrain[:tstep-1,i,1,1]/P.dt,'b.',label=r'$\dot\varepsilon_{yy}$')
            plt.plot(self.t[:tstep-1],self.gammadot[:tstep-1,i],'k.',label=r'$\dot\gamma$')
            plt.xlabel(r'$t$',rotation='horizontal')
            plt.ylabel(r'$\dot\varepsilon_{ij}$',rotation='horizontal')

            plt.subplot(236)
            plt.semilogx(self.t[:tstep-1],self.mu[:tstep-1,i]/(-self.tau[:tstep-1,i]/self.pressure[:tstep-1,i]),'k.')
            plt.xlabel(r'$t$',rotation='horizontal')
            plt.ylabel(r'$\mu/\mu_{applied}$',rotation='horizontal')
            plt.ylim(0.5,1.5)

            plot.savefig(P,str(i))
        #print (self.q[tstep-1,0]-self.q[1,0])/(self.p[tstep-1,0]-self.p[1,0])
