from numpy import linspace, sin, cos, pi, zeros, outer, array, dot
from numpy import trunc, arctan, eye, trace, nan_to_num, tensordot
from numpy import sqrt, abs, ones, random#, min, max
from numpy.linalg import norm
from constit import *

class Particle():
    def __init__(self,x,y,P,G):
#        self.rho = P.S.rho # density (kg/m^3)
        self.rho = random.rand()*P.S.rho # density (kg/m^3)
        self.x = array([x, y, 0.]) # position (m)
        self.v = array([0.,0.,0.]) # velocity (m/s)
        self.V = P.S.A*P.G.thickness # volume (m^3)
        self.m = self.rho*self.V # mass (kg)
        self.b = array([P.g*sin(P.theta),P.g*cos(P.theta),0.]) # body forces
        if P.S.law == 'dp':
#            self.strain = -1e-15*eye(3) # isotropic compression
            self.strain = -1e-15*ones((3,3)) # mixed state compression
        else:
            self.strain = zeros((3,3)) # corresponding strains
        self.dstrain = zeros((3,3)) # change in strain
        self.stress = (P.S.K*trace(self.strain)*eye(3) +
                       2.*P.S.G*(self.strain - trace(self.strain)*eye(3)/3.))
        self.sigma_kk = trace(self.stress)
        self.dev_stress = self.stress - self.sigma_kk*eye(3)/3.
        self.gammadot = 0.
        self.pressure = 0.
        self.sigmav = 0.
        self.sigmah = 0.
        self.yieldfunction = 0.
        self.G = zeros((4,3)) # gradient of shape function
        self.N = zeros((4)) # basis function
        self.n_star = 0 # reference node
        try:
            self.phi = P.GSD.phi
            self.dphi = zeros(self.phi.shape)
            self.dm = 0.
            self.neighbours = 0.
            self.seg_v = zeros((self.phi.shape[0],3))
#            self.seg_v[:,1] = linspace(-1,1,self.phi.shape[0])
            self.seg_v[:,0] = random.rand()*linspace(-1,1,self.phi.shape[0])
        except:
            pass

    def update_strain(self,P,G):
        self.dstrain = zeros((3,3))
        for r in xrange(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                v = G.q[n]/G.m[n]
            else:
                v = array([0.,0.,0.])
            self.dstrain += 0.5*(outer(v,self.G[r])
                               + outer(self.G[r],v))*P.dt
        self.strain += self.dstrain

    def update_stress(self,P,G):
        exec(P.S.law + '(self,P,G)')

    def get_reference_node(self,P,G):
        self.n_star = int(trunc((self.x[0] - P.G.x_m)/G.dx)
                          + trunc((self.x[1] - P.G.y_m)/G.dy)*P.G.nx)

    def get_basis_functions(self,P,G):
        self.N, self.G = G.N(self.x-[G.X[self.n_star],G.Y[self.n_star],0.])

    def get_nodal_mass_momentum(self,P,G):                
        for r in xrange(4):
            n = G.nearby_nodes(self.n_star,r,P)
            G.m[n] += self.N[r]*self.m
            G.q[n] += self.N[r]*self.m*self.v

    def get_nodal_gsd(self,P,G):
        for position,s in enumerate(P.GSD.s):
            if s == self.s:
                pos = position
        for r in xrange(4):
            n = G.nearby_nodes(self.n_star,r,P)
            G.gsd[n,pos] += self.N[r]*self.m/self.rho # total volume in size bin

    def move_material_points(self,P,G):
        for r in xrange(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                self.x += P.dt*self.N[r]*G.q[n]/G.m[n]
                self.v += P.dt*self.N[r]*G.q_dot[n]/G.m[n]

    def recover_position(self,P,G):
        x = zeros((3,))
        nodes = zeros((4))
        for r in xrange(4):
            n = G.nearby_nodes(self.n_star,r,P)
            x += self.N[r]*array([G.X[n],G.Y[n],0.])
            nodes[r] = n
        return x, nodes

    def get_nodal_forces(self,P,G):
        for r in xrange(4):
            n = G.nearby_nodes(self.n_star,r,P)
            G.fi[n] += self.V*dot(self.stress,self.G[r])
            G.fe[n] += self.N[r]*self.m*self.b
            # if there are tractions add something here

    def increment_gravity(self,P,G):
        self.b = array([P.g*sin(P.theta),P.g*cos(P.theta),0.]) # body forces

    def energy(self,P):
        k = 0.5*self.m*sum(self.v**2) # kinetic
        h = (sum(self.x**2))**0.5*sin(P.theta+arctan(abs(self.x[1]/self.x[0]))) # height
        g = -self.m*P.g*h # gravitational potential
        s = 0.5*sum(sum(self.stress*self.strain))*self.V # elastic potential energy
        return [k,g,s]
        
    def find_f(self,P,L):
        self.f = P.GSD.s/(sum(P.GSD.s*self.phi))
        self.rho_bar = sum(self.rho*self.phi)

    def find_seg_velocity(self,P,L):
        self.u = zeros((len(self.neighbours),3))
        for k in xrange(len(self.neighbours)):
            j = int(self.neighbours[k])
            n = (L.S[j].x - self.x)
            # find seg velocity u
            r = self.phi*self.rho*P.g
            s = self.phi*dot(L.S[j].stress*L.S[j].f - self.stress*self.f,n)
            m = dot(L.S[j].m*L.S[j].v - self.m*self.v,n) # material derivative of momentum
            self.u[k] = self.gammadot/(self.rho_bar*self.phi*P.S.c)*(r - s - m) # seg velocity
            print self.u[k]
        
    def find_flux(self,P,L):
        self.dQ = zeros((len(self.neighbours)))
        self.dm = zeros((len(self.neighbours)))
        for k in xrange(len(self.neighbours)):
            j = int(self.neighbours[k])
            if j >= 0:
                n = (L.S[j].x - self.x)
                d = norm(n) # distance
                n_hat = n/d
                # use seg velocity as advection term below
                self.dQ[k] = -(P.S.D/d*(self.m - L.S[j].m) - dot(self.m*self.u[k] - L.S[j].m*L.S[j].u[0],n_hat)) # advection of seg velocity only
                self.dQ[k] = min(self.m/P.dt/2.,max(0.,self.dQ[k])) # min and max bounds on flux
                self.dm[k] = float(self.areas[k])*self.dQ[k]*P.dt

    def move_flux(self,P,L):
        for k in xrange(len(self.neighbours)):
            j = int(self.neighbours[k])
            if j >= 0:
                self.m += self.dm[k]
                L.S[j].m -= self.dm[k]
                #print self.rho
