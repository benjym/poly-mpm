import os
from numpy import pi,array,dot,tan,zeros,eye,cumsum,outer,sum,trace,exp,linspace
from numpy import random, cos, sin, ceil, around
import matplotlib.pyplot as plt

class Params():
    def __init__(self):
        self.dt = 1e-3 # timestep (s)
        self.t = 0. # physical time (s)
        self.tstep = 0
        self.savetime = .1
        self.t_f = 5.#100*self.dt # final time (s)
        self.nt = int(self.t_f/self.dt) # number of timesteps
        self.save = 0 # save counter
        self.M_tol = 1e-10 # very small mass (kg)
        self.max_g = -10. # gravity (ms^-2)
        self.max_q = 0.
        self.update_forces()
        self.theta = 0.*pi/180. # slope angle (degrees)
        self.thickness = 1. # (m) into page
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.O = Output_Params(self.nt)
        self.S = Solid_Params()
        self.F = Fluid_Params()
        self.R = Fluid_Params()
        self.has_yielded = False
        self.damping = False

    def update_forces(self):
        t_c = .5
        #self.g = self.max_g*(1.-exp(-3.*self.t**2/t_c**2))
        self.g=self.max_g
        self.q = self.max_q*(1.-exp(-3.*self.t**2/t_c**2))

class Grid_Params():
    def __init__(self):
        self.x_m = 0.0 # (m)
        self.x_M = 1.0 # (m)
        self.y_m = 0.0 # (m)
        self.y_M = 2.0 # (m)
        self.nx = 11
        self.ny = 21
        
class Boundary_Params():
    def __init__(self):
        self.wall = False
        self.has_top = False
        self.has_bottom = True
        self.has_right = False
        self.has_left = False
        self.outlet_left = False
        self.force_boundaries = False
        self.vertical_force = False
        self.horizontal_force = False

class Solid_Params():
    def __init__(self):
        self.X = []
        self.Y = []
        self.n = 0

        self.mode = 'ball' # ball, blob, slope, excavate
        self.law = 'elastic'
#        self.law = 'von_mises'
        self.rho = 1000. # density (kg/m^3)
        
        self.E = 1.e6 # elastic modulus (Pa)
        self.nu = 0. # poisson's ratio
        self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
        self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

        self.s = 2.5
        self.k = self.E/100.

        nr = 10 # particles in radial direction
        nphi = 20 # particles around circumference
        r = 0.3 # radius
        c = [0.5,1.] # centre
        
        for i in linspace(0,r,nr):
            dnphi = around(nphi*i/r) # number in this ring
            for j in linspace(0,(1.-1./(dnphi))*2*pi,dnphi):
                self.X.append(c[0]+i*sin(j))
                self.Y.append(c[1]+i*cos(j))
                self.n += 1
        self.A = pi*0.2**2/self.n # area (m^2)
        
        
class Output_Params():
    def __init__(self,nt):
        self.measure_energy = True
        self.plot_continuum = False
        self.plot_material_points = True
        self.measure_stiffness = False
        self.check_positions = True
        self.plot_fluid = False
        self.energy = zeros((nt+1,4)) # energy
        
    def measure_E(self,P,L):
        print 'Measuring macro and micro stress/strain for each material point... '
        for i in xrange(P.S.n):
            original_position = array((P.S.X[i],P.S.Y[i],0))
            macro_strain = (original_position-L.S[i].x)/array((P.S.L,P.S.W,1.)) #original_position
            macro_stress = P.max_conf#/(P.S.L*P.S.W)
            print 'From macroscopic stress/strain:'
            print macro_stress/macro_strain/P.S.E
            print 'From microscopic stress/strain:'
            print L.S[i].dstress/L.S[i].dstrain/P.S.E

            
class Fluid_Params():
    def __init__(self):
        self.n = 0

        
class Rigid_Params():
    def __init__(self):
        self.n = 0

