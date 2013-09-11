from numpy import *
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
        self.S = []
        for p in range(P.phases):
            self.S.append([]) # particles
            for i in xrange(P.S[p].n):
                self.S[p].append(Particle(P.S[p].X[i],P.S[p].Y[i],P,G,p))

    def outlet_left(self,P,G):
        for p in range(P.phases):
            for i in linspace(P.S[p].n-1,0,P.S[p].n):
                if self.S[p][int(i)].x[0] < P.G.x_m+P.G.dx:
                    del self.S[p][int(i)]
                    P.S[p].n -= 1

    def inlet_right(self,P,G):
        for p in range(P.phases):
            try:
                if P.t%P.S[p].inlet_rate < P.dt and P.t:
                    for y in P.S[p].yp:
                        self.S[p].append(Particle(P.S[p].xp[-1],y,P,G,p))
                        P.S[p].n += 1
            except: pass
                
    def cyclic_lr(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                if self.S[p][i].x[0] < P.G.x_m:
                    self.S[p][i].x[0] += (P.G.x_M - P.G.x_m)
                if self.S[p][i].x[0] > P.G.x_M:
                    self.S[p][i].x[0] -= (P.G.x_M - P.G.x_m)
                
    def get_reference_node(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                self.S[p][i].get_reference_node(P,G)
            
    def get_basis_functions(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                self.S[p][i].get_basis_functions(P,G)
            
    def get_nodal_mass_momentum(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                self.S[p][i].get_nodal_mass_momentum(P,G)
        G.q[:,0] = G.q[:,0]*(1.-G.boundary_v) # 0 at boundary
        G.q[:,1] = G.q[:,1]*(1.-G.boundary_h) # 0 at boundary
        if P.B.roughness:
            G.q[:,0] = G.q[:,0]*(1.-G.boundary_h) # roughness on top/bottom
            G.q[:,1] = G.q[:,1]*(1.-G.boundary_v) # roughness on sidewalls
        
    def update_stress_strain(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                self.S[p][i].update_strain(P,G)
                self.S[p][i].update_stress(P,G)
        if P.B.cyclic_lr:
            G.apply_cyclic_BCs_stress(P)
            
    def get_nodal_forces(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                self.S[p][i].get_nodal_forces(P,G)
        if P.B.cyclic_lr:
            G.apply_cyclic_nodal_forces(P)
            
    def move_material_points(self,P,G):
#        print self.S[0][25].v[0]
        for p in range(P.phases):
            if P.S[p].law is not 'rigid':
                for i in xrange(P.S[p].n):
                    self.S[p][i].move_material_points(P,G)
            
    def update_forces(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                self.S[p][i].increment_gravity(P,G,p)
            
    def get_nodal_gsd(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                self.S[p][i].get_nodal_gsd(P,G)
        G.gsd = G.gsd/repeat(sum(G.gsd,1),P.GSD.ns).reshape(G.gsd.shape)
        G.s_bar = sum(G.gsd*P.GSD.s,1)
                
    def measure_elastic_energy(self,P,G):
        K=0.;G=0.;S=0.
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                [k,g,s] = self.S[p][i].energy(P)
                K += k; G += g; S += s
        return [P.t,K,G,S]
        
    def recover_position(self,P,G):
        for p in range(P.phases):
            for i in xrange(P.S[p].n):
                x, nodes = self.S[p][i].recover_position(P,G)
                if max(abs(x-self.S[p][i].x)) > 1e-10:
                    print 'Particle ' + str(i) + ' lost by ',
                    print x-self.S[p][i].x
                    print 'Particle was at: ',
                    print x
                    print 'Surrounded by: '
                    for i in xrange(4):
                        print G.X[nodes[i]], G.Y[nodes[i]]
            
    def update_timestep(self,P):
        if not P.has_yielded:
            for p in range(P.phases):
                for i in xrange(P.S[p].n):
                    if self.S[p][i].yieldfunction > -0.01:
                        P.has_yielded = True
        if P.has_yielded:
            print('FAILED at time ' + str(P.t-P.dt) +
                  ' with g = ' + str(P.g) +
                  ' and q = ' + str(P.q))
            P.dt *= .1
                
    def get_tessellation(self,P,G,plot):
        for p in range(P.phases):
            x = zeros((P.S[p].n,3))
            x[:,0] = [self.S[p][i].x[0] for i in xrange(P.S[p].n)]
            x[:,1] = [self.S[p][i].x[1] for i in xrange(P.S[p].n)]
            voronoi2D(x,P,self)

    def move_grainsize(self,P):
        for p in range(P.phases):
            for i in xrange(P.S[p].n): # find segregation velocity
                self.S[p][i].find_f(P,self)
            for i in xrange(P.S[p].n): # find segregation velocity
                self.S[p][i].find_seg_velocity(P,self)
            for i in xrange(P.S[p].n): # find change
                self.S[p][i].find_flux(P,self)
            #mass = 0.
            for i in xrange(P.S[p].n): # accumulate changes
                self.S[p][i].move_flux(P,self)
            #    mass += self.S[p][i].m
            #print 'total mass is ' + str(mass)
