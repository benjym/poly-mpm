from numpy import linspace, sin, cos, pi, zeros, outer, dot, array, ones, sum
from numpy import trunc, arctan, eye, trace, tile, repeat, reshape, tan, abs
from numpy import random
import matplotlib.pyplot as plt
from particle import Particle
from voronoi import voronoi2D
import matplotlib.cm as cm
import sys

class MatPointList():
    """
    A class of material points of a single size.
    Contains convenience functions for operating
    on many material points.
    """
    def __init__(self,P,G):
        self.S = [] # particles
        for i in xrange(P.S.n):
            self.S.append(Particle(P.S.X[i],P.S.Y[i],P,G))
                
    def outlet_left(self,P,G):
        for i in linspace(P.S.n-1,0,P.S.n):
            if self.S[int(i)].x[0] < P.G.x_m+P.G.dx:
                del self.S[int(i)]
                P.S.n -= 1
                
    def get_reference_node(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].get_reference_node(P,G)
            
    def get_basis_functions(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].get_basis_functions(P,G)
            
    def get_nodal_mass_momentum(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].get_nodal_mass_momentum(P,G)
        G.q[:,0] = G.q[:,0]*(1.-G.boundary_v) # 0 at boundary
        G.q[:,1] = G.q[:,1]*(1.-G.boundary_h) # 0 at boundary
        if P.B.roughness:
            G.q[:,0] = G.q[:,0]*(1.-G.boundary_h) # roughness on top/bottom
            G.q[:,1] = G.q[:,1]*(1.-G.boundary_v) # roughness on sidewalls
            
    def update_stress_strain(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].update_strain(P,G)
            self.S[i].update_stress(P,G)
            
    def get_nodal_forces(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].get_nodal_forces(P,G)
            
    def move_material_points(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].move_material_points(P,G)
            
    def update_forces(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].increment_gravity(P,G)
            
    def get_nodal_gsd(self,P,G):
        for i in xrange(P.S.n):
            self.S[i].get_nodal_gsd(P,G)
        G.gsd = G.gsd/repeat(sum(G.gsd,1),P.GSD.ns).reshape(G.gsd.shape)
        G.s_bar = sum(G.gsd*P.GSD.s,1)
                
    def measure_elastic_energy(self,P,G):
        K=0.;G=0.;S=0.
        for i in xrange(P.S.n):
            [k,g,s] = self.S[i].energy(P)
            K += k; G += g; S += s
        return [P.t,K,G,S]
        
    def recover_position(self,P,G):
        for i in xrange(P.S.n):
            x, nodes = self.S[i].recover_position(P,G)
            if max(abs(x-self.S[i].x)) > 1e-10:
                print 'Particle ' + str(i) + ' lost by ',
                print x-self.S[i].x
                print 'Particle was at: ',
                print x
                print 'Surrounded by: '
                for i in xrange(4):
                    print G.X[nodes[i]], G.Y[nodes[i]]
            
    def update_timestep(self,P):
        if not P.has_yielded:
            for i in xrange(P.S.n):
                if self.S[i].yieldfunction > -0.01:
                    P.has_yielded = True
            if P.has_yielded:
                print('FAILED at time ' + str(P.t-P.dt) +
                      ' with g = ' + str(P.g) +
                      ' and q = ' + str(P.q))
                P.dt *= .1
                
    def get_tessellation(self,P,G,plot):
        x = zeros((P.S.n,3))
        x[:,0] = [self.S[i].x[0] for i in xrange(P.S.n)]
        x[:,1] = [self.S[i].x[1] for i in xrange(P.S.n)]
        voronoi2D(x,P,self)

    def move_grainsize(self,P):
        for i in xrange(P.S.n): # find segregation velocity
            self.S[i].find_f(P,self)
        for i in xrange(P.S.n): # find segregation velocity
            self.S[i].find_seg_velocity(P,self)
        for i in xrange(P.S.n): # find change
            self.S[i].find_flux(P,self)
        #mass = 0.
        for i in xrange(P.S.n): # accumulate changes
            self.S[i].move_flux(P,self)
        #    mass += self.S[i].m
        #print 'total mass is ' + str(mass)
