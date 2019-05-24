from numpy import *
import matplotlib.pyplot as plt
from particle import Particle
# from voronoi import voronoi2D
import matplotlib.cm as cm
import sys
# from numba import cuda

class MatPointList():
    """
    A class of material points of a single size.
    Contains convenience functions for operating
    on many material points.
    """
    def __init__(self,P,G):
        """
        This constructor creates and appends all material points defined in parameter file to the material point list.
        """
        self.S = []
        for p in range(P.phases):
            self.S.append([]) # particles
            for i in range(P.S[p].n):
                if P.S[p].heterogeneous:
                    self.S[p].append(Particle(P.S[p].X[i],P.S[p].Y[i],P,G,p,phi=P.S[p].PHI[i]))
                else:
                    self.S[p].append(Particle(P.S[p].X[i],P.S[p].Y[i],P,G,p))

    def outlet_left(self,P,G):
        """Particles are deleted if they move left (:math:`-x` direction) out of the domain.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in linspace(P.S[p].n-1,0,P.S[p].n):
                if self.S[p][int(i)].x[0] < P.G.x_m+P.G.dx:
                    del self.S[p][int(i)]
                    P.S[p].n -= 1

    def outlet_bottom(self,P,G):
        """Particles are deleted if they move downwards (:math:`-y` direction) out of the domain.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            # for i in linspace(P.S[p].n-1,0,P.S[p].n):
            for i in range(P.S[p].n-1,-1,-1):
                if self.S[p][i].x[1] < P.G.y_m:
                    del self.S[p][i]
                    P.S[p].n -= 1

    def inlet_right(self,P,G):
        """Particles are created at right boundary (:math:`+x` direction) of the domain.

        This is not physical at all. Needs reimplementation.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
                if P.t%P.S[p].inlet_rate < P.dt and P.t:
                    for y in P.S[p].yp:
#                         self.S[p].append(Particle(P.S[p].xp[-1],y,P,G,p))
                        self.S[p].append(Particle(P.G.x_M - P.S[p].gap[0],y,P,G,p))
                        P.S[p].n += 1

    def inlet_top(self,P,G):
        """Particles are created at the top (:math:`+y` direction) in the middle of the domain.

        This is not physical at all. Needs reimplementation.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
                if P.t%P.S[p].inlet_rate < P.dt and P.t:
                    # for y in P.S[p].yp:
    #                         self.S[p].append(Particle(P.S[p].xp[-1],y,P,G,p))
                    self.S[p].append(Particle((P.G.x_M - P.G.x_m)/2.,P.G.y_M-P.S[p].gap[1],P,G,p))
                    P.S[p].n += 1

    def cyclic_lr(self,P,G):
        """Particles which leave the domain in the :math:`x` direction are placed back in the other side.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                if self.S[p][i].x[0] < P.G.x_m:
                    self.S[p][i].x[0] += (P.G.x_M - P.G.x_m)
                if self.S[p][i].x[0] > P.G.x_M:
                    self.S[p][i].x[0] -= (P.G.x_M - P.G.x_m)

    def get_reference_node(self,P,G):
        """For every particle, get its reference node.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].get_reference_node(P,G)

    def get_basis_functions(self,P,G):
        """For every particle, get its basis function.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].get_basis_functions(P,G)

    def get_nodal_mass_momentum(self,P,G):
        """For every particle, map mass and momentum to the grid. Then, set momentum at boundaries to zero if applicable.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].get_nodal_mass_momentum(P,G)
            if P.segregate_grid:
                for j in range(P.G.ns):
                    G.phi[:,j] = G.phim[:,j]/G.m

        if P.segregate_grid:
            G.phi = nan_to_num(G.phi)
            G.s_bar = sum(G.S*G.phi,1)
            G.pk = G.pkm/G.m
            # print(G.pk)
            # print(G.phi)
#             print(G.S)
            # print(G.s_bar)

        G.q[:,0] = G.q[:,0]*(1.-G.boundary_v) # 0 at boundary
        G.q[:,1] = G.q[:,1]*(1.-G.boundary_h) # 0 at boundary
        if P.B.roughness:
            G.q[:,0] = G.q[:,0]*(1.-G.boundary_h) # roughness on top/bottom
            G.q[:,1] = G.q[:,1]*(1.-G.boundary_v) # roughness on sidewalls

    def update_stress_strain(self,P,G):
        """For every particle, firstly calculate the incremental strain from nearby nodal velocities. Then use the relevant constitutive law to update the stress at the material point. Finally, if cyclic boundaries are active, update the stress at the boundaries.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].update_strain(P,G)
                self.S[p][i].update_stress(P,G,p)
        if P.B.cyclic_lr: G.make_cyclic(P,G,['pressure', 'gammadot', 'dev_stress', 'yieldfunction', 'mu_s', 'mu', 'I'])

    def get_nodal_forces(self,P,G):
        """For every particle, get nodal forces. If cyclic boundaries are active, update the nodal cyclic forces.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].get_nodal_forces(P,G)
        if P.B.cyclic_lr: G.make_cyclic(P,G,['fi','fe'])

    def move_material_points(self,P,G):
        """For every particle, update the particle location.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            if P.S[p].law is 'rigid':
                for i in range(P.S[p].n):
                    self.S[p][i].move_rigid_points(P,G)
            else:
                for i in range(P.S[p].n):
                    self.S[p][i].move_material_points(P,G)

    def update_forces(self,P,G):
        """For every particle, update the value of the global gravity.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].increment_gravity(P,G,p)

    def get_nodal_gsd(self,P,G):
        """Measure the volume weighted average size at the node.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """

        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].update_nodal_gsd(P,G)
        if P.segregate_grid or P.segregate_mp:
            # if not P.B.cyclic_lr: G.s_bar /= G.m
            G.s_bar /= G.m
            # print(G.s_bar)
#     def get_mp_gsd(self,P,G):
#         """Map mass weighted average size back to the material points.
#
#         :param P: A param.Param instance.
#         :param G: A grid.Grid instance.
#
#         """
#         for p in range(P.phases):
#             for i in range(P.S[p].n):
#                 self.S[p][i].map_gsd_to_mp(P,G)
#
# #         G.gsd = G.gsd/repeat(sum(G.gsd,1),P.GSD.ns).reshape(G.gsd.shape)
# #         G.s_bar = sum(G.gsd*P.GSD.s,1)

    def measure_elastic_energy(self,P,G):
        """I've forgotten what this does.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        K=0.;G=0.;S=0.
        for p in range(P.phases):
            for i in range(P.S[p].n):
                [k,g,s] = self.S[p][i].energy(P)
                K += k; G += g; S += s
        return [P.t,K,G,S]

    def recover_position(self,P,G):
        """Just used for checking if the code is working. Checks that the particle mapping recovers the real particle position.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                x, nodes = self.S[p][i].recover_position(P,G)
                if max(abs(x-self.S[p][i].x)) > 1e-10:
                    print('Particle ' + str(i) + ' lost by ' +
                          x-self.S[p][i].x +
                          'Particle was at: ' +
                          x +
                          'Surrounded by: ',
                          end='\r')
                    for i in range(4):
                        print(G.X[nodes[i]], G.Y[nodes[i]])

    def update_timestep(self,P):
        """Check through every particle, and if any have yielded, reduce the timestep.

        :param P: A param.Param instance.

        """
        if not P.has_yielded:
            for p in range(P.phases):
                for i in range(P.S[p].n):
                    if self.S[p][i].yieldfunction > 0.99:#-0.01:
                        P.has_yielded = True
        if P.has_yielded:
            print('FAILED at time ' + str(P.t-P.dt) +
                  ' with g = ' + str(P.g) +
                  ' and q = ' + str(P.q))
            print('nt was ' + str(P.nt))
            time_reduction = 20.
            nt_rem = P.nt - P.tstep
            P.nt = P.tstep + int(nt_rem*time_reduction)
            P.dt /= time_reduction
            P.savetime /= time_reduction
            print('nt is ' + str(P.nt))


    def move_grainsize_on_grid(self,P,G):
        """Advect grainsize distribution across grid.

        :param P: A param.Param instance.

        """
        for p in range(P.phases):
            for i in range(P.S[p].n):
                self.S[p][i].map_phi_to_mp(P,G)

#         for p in range(P.phases): # CHEATING!!!!!!!
#             for i in range(P.S[p].n):
#                 self.S[p][i].phi[self.S[p][i].phi<0] = 0
#                 self.S[p][i].phi[self.S[p][i].phi>1] = 1
#                 self.S[p][i].phi /= sum(self.S[p][i].phi)
