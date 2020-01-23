# import cupy as cp
import numpy as np
import matplotlib.pyplot as plt
from constit import *

class Particles():
    def __init__(self,P,G):
        """
        This constructor initialises material point arrays from the  input file as necessary.
        """

        self.x = np.zeros([P.S.n,3])
        self.x[:,0] = P.S.X
        self.x[:,1] = P.S.Y
        if hasattr(P.S, 'PHI'): self.phi = P.S.PHI
        else: self.phi = np.ones([P.G.ns])/float(P.G.ns)

        # Shape functions and gradients
        self.N = np.zeros([P.S.n,4])
        self.G = np.zeros([P.S.n,4,3]) # {x,y,z}

        # Reference nodes
        self.n_star = np.empty([P.S.n],dtype='int')
        self.n = np.empty([P.S.n,4],dtype='int')
        self.valid_nodes = np.empty([P.S.n,4],dtype='bool') # just the n's where the mass is over the tolerance

        # Helper functions
        self.eye = np.repeat(np.expand_dims(np.eye(3),0),P.S.n,axis=0) # eye function for each MP

        self.rho = P.S.rho # density (kg/m^3)
        try: self.v = P.S.initial_v
        except: self.v = np.zeros([P.S.n,3]) # velocity (m/s)
        self.V = P.S.A*P.G.thickness # volume (m^3)
        self.m = self.rho*self.V # mass (kg)

        # if P.pressure == 'lithostatic':
        #     K_0 = 1.
        #     if hasattr(P.G, 'top_gap'):
        #         sigma_yy = cos(P.theta)*self.rho*P.max_g*(P.G.y_M - P.G.top_gap - self.x[:,1]) # build in a pressure gradient
        #         sigma_xy = sin(P.theta)*self.rho*P.max_g*(P.G.y_M - P.G.top_gap - self.x[:,1]) # build in a pressure gradient
        #         sigma_xx = K_0*sigma_yy # build in a pressure gradient
        #         sigma_zz = K_0*sigma_yy # build in a pressure gradient
        #     else:
        #         sigma_yy = cos(P.theta)*self.rho*P.max_g*(P.G.y_M - self.x[:,1]) # build in a pressure gradient
        #         sigma_xy = sin(P.theta)*self.rho*P.max_g*(P.G.y_M - self.x[:,1]) # build in a pressure gradient
        #         sigma_xx = K_0*sigma_yy # build in a pressure gradient
        #         sigma_zz = K_0*sigma_yy # build in a pressure gradient
        #     self.stress = np.array([[sigma_xx,sigma_xy,0],[sigma_xy,sigma_yy,0],[0,0,sigma_zz]])
        #     # self.strain = zeros((3,3)) # strain-free
        #     self.strain = MP.eye(3)*np.trace(self.stress)/(9.*P.S[p].K) + (self.stress - MP.eye(3)*np.trace(self.stress)/3.)/(2.*P.S[p].G)
        # elif P.pressure == 'non-zero':
        #     self.strain = -1e-15*np.ones((P.S.n,3,3)) # mixed state compression
        #     self.stress = (P.S[p].K*np.trace(self.strain)*MP.eye(3) +
        #                    2.*P.S[p].G*(self.strain - np.trace(self.strain)*MP.eye(3)/3.))
        # elif P.pressure == 'temperature':
        #     self.strain = -1e-5*MP.eye(3) # pure compression everywhere
        #     self.stress = (P.S[p].K*np.trace(self.strain)*np.eye(3) +
        #                    2.*P.S[p].G*(self.strain - np.trace(self.strain)*np.eye(3)/3.))
        # elif P.pressure == 'compression':
        #     self.stress = np.array([[-P.initial_pressure,0,0],
        #                             [0,-P.initial_pressure,0],
        #                             [0,0,-P.initial_pressure]])
        #     # self.strain = zeros((3,3)) # strain-free
        #     self.strain = MP.eye(3)*np.trace(self.stress)/(9.*P.S[p].K) + (self.stress - MP.eye(3)*np.trace(self.stress)/3.)/(2.*P.S[p].G)
        # else:
        self.strain = np.zeros((P.S.n,3,3)) # corresponding strains
        self.stress = np.zeros((P.S.n,3,3))
        # self.dstrain = np.zeros((P.S.n,3,3)) # change in strain
        # self.dstress = np.zeros((P.S.n,3,3)) # change in strain
        # self.de_ij = np.zeros((P.S.n,3,3))
        # self.sigma_kk = np.trace(self.stress)
        # self.dev_stress = self.stress - self.sigma_kk*MP.eye(3)/3.
        # self.gammadot = 0.
        # self.pressure = self.sigma_kk/3.
        # self.sigmav = self.stress[1,1]
        # self.sigmah = self.stress[0,0]
        # self.yieldfunction = 0.

        # self.pk = 0. # kinetic pressure
        # if phi is not None:
        #     self.phi = array(phi)
        # else:
        #     self.phi = array(P.S[p].phi)
        # if P.initial_flow == 'steady':
        #     if P.S[0].law == 'viscous' or P.S[0].law == 'viscous_size':
        #         self.v = array([self.rho*P.max_g*sin(P.theta)*(P.G.y_M*y - y**2/2.)/P.S[p].mu_s,0.,0.]) # VISCOUS FLOW DOWN A SLOPE
        #     elif P.S[0].law == 'pouliquen' or P.S[0].law == 'pouliquen2D':
        #         s_bar = (P.G.s_m+P.G.s_M)/2.
        #         self.v = array([sqrt(abs(P.max_g)*P.G.s_bar_0)*(2./3.)*P.S[p].I_0*
        #                        (tan(abs(P.theta))-P.S[p].mu_0)/(P.S[p].mu_1-tan(abs(P.theta)))*
        #                        sqrt(P.S[p].packing*cos(abs(P.theta)))*
        #                        (P.G.y_M**1.5 - (P.G.y_M - y)**1.5)/P.G.s_bar_0**1.5,
        #                        0.,0.])
        # elif P.initial_flow == 'steady_poiseulle': # https://link.springer.com/content/pdf/10.1007%2Fs10035-013-0447-3.pdf
        #     H = P.G.y_M - P.G.y_m
        #     y_star = y/H
        #     if y_star > 0.5: y_star = 0.5 - abs(y_star - 0.5) # y_star only defined on LHS
        #     K = abs(P.S[0].rho*P.max_g)
        #     Lambda = H/P.G.s_bar_0
        #     beta = sqrt(abs(K*H/P.initial_pressure))
        #     y_star_c = 0.5 - P.S[p].mu_0/(beta**2)
        #     if y_star < y_star_c:
        #         u_star = (Lambda*P.S[p].I_0)/(beta**3)*\
        #                 ((beta**2)*y_star + (P.S[p].mu_0 - P.S[p].mu_1)*log(-(2.*(beta**2)*y_star + 2*P.S[p].mu_1 - beta**2)/(beta**2 - 2*P.S[p].mu_1)))
        #     else:
        #         u_star = (Lambda*P.S[p].I_0)/(beta**3)*\
        #                  ((beta**2)/2. - P.S[p].mu_0 + (P.S[p].mu_0 - P.S[p].mu_1)*log((P.S[p].mu_1-P.S[p].mu_0)/(P.S[p].mu_1-(beta**2)/2.)))
        #     self.v = array([-u_star*sqrt(K*H/P.S[p].rho_s),0.,0.])
        #
        #     stable = beta<sqrt(2*P.S[p].mu_1) # if True, this should work
        #     if not stable:
        #         print('WARNING: THIS SETUP WILL ACCELERATE FOREVER')
        #         print(sqrt(2*P.S[p].mu_0),beta,sqrt(2*P.S[p].mu_1))

    def cyclic_lr(self,P):
        """Particles which leave the domain in the :math:`x` direction are placed back in the other side.

        :param p_x: An array containing the locations of all of the particles.
        :param G: A grid.Grid instance.

        """
        self.x[self.x[:,0]<P.G.x_m,0] += P.G.x_M - P.G.x_m
        self.x[self.x[:,0]>P.G.x_M,0] -= P.G.x_M - P.G.x_m

    def get_reference_nodes(self,P,G):
        """For every particle, get its reference node (n_star) and nearby nodes (n).

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        self.n_star = (np.trunc((self.x[:,0] - P.G.x_m)/G.dx)
                     + np.trunc((self.x[:,1] - P.G.y_m)/G.dy)*P.G.nx).astype('int')

        self.n[:,0] = self.n_star
        self.n[:,1] = self.n_star + 1
        self.n[:,2] = self.n_star + P.G.nx + 1
        self.n[:,3] = self.n_star + P.G.nx


    def get_basis_functions(self,P,G):
        """For every particle, get its basis function.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        rel_pos = self.x - np.array([G.X[self.n_star],G.Y[self.n_star],[0]*P.S.n]).T

        # Lagrange 4-node system
        self.N[:,0] = (G.dx-rel_pos[:,0])*(G.dy-rel_pos[:,1])/(G.dx*G.dy)
        self.N[:,1] = (rel_pos[:,0])*(G.dy-rel_pos[:,1])/(G.dx*G.dy)
        self.N[:,2] = (rel_pos[:,0])*(rel_pos[:,1])/(G.dx*G.dy)
        self.N[:,3] = (G.dx-rel_pos[:,0])*(rel_pos[:,1])/(G.dx*G.dy)
        self.G[:,0,:] = np.array([-(G.dy-rel_pos[:,1]),-(G.dx-rel_pos[:,0]),[0]*P.S.n]).T/(G.dx*G.dy)
        self.G[:,1,:] = np.array([ (G.dy-rel_pos[:,1]),-rel_pos[:,0],       [0]*P.S.n]).T/(G.dx*G.dy)
        self.G[:,2,:] = np.array([ rel_pos[:,1],        rel_pos[:,0],       [0]*P.S.n]).T/(G.dx*G.dy)
        self.G[:,3,:] = np.array([-rel_pos[:,1],        (G.dx-rel_pos[:,0]),[0]*P.S.n]).T/(G.dx*G.dy)

    def get_nodal_mass_momentum(self,P,G):
        """For every particle, map mass and momentum to the grid. Then, set momentum at boundaries to zero if applicable.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for r in range(4):
            G.m[self.n[:,r]] += self.N[:,r]*self.m
            for i in range(2): G.q[self.n[:,r],i] += self.N[:,r]*self.m*self.v[:,i]
            if P.segregate:
                G.phim[self.n[:,r]] += self.N[:,r]*self.m*self.phi
                G.pkm[self.n[:,r]]  += self.N[:,r]*self.m*self.pk

        if P.segregate:
            for j in range(P.G.ns):
                G.phi[:,j] = G.phim[:,j]/G.m

        if P.segregate:
            G.phi = nan_to_num(G.phi)
            G.s_bar = sum(G.S*G.phi,1)
            G.pk = G.pkm/G.m

        # NOTE: WHY IS THIS HERE??????????? SHOULD ALREADY BE SET BY SETTING q_dot = 0???
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
        # Update strain
        self.dstrain = np.zeros([P.S.n,3,3])
        self.L = self.L_T = np.zeros([P.S.n,3,3])
        for r in range(4):
            v = np.zeros([P.S.n,3])
            for j in range(3):
                v[:,j] = G.q[self.n[:,r],j]/G.m[self.n[:,r]]*G.valid_nodes[self.n[:,r]]
            for i in range(P.S.n): # FIXME - NEED TO REMOVE THIS LOOP
                L = np.outer(v[i],self.G[i,r])
                L_T = np.outer(self.G[i,r],v[i])
                self.L[i] += L
                self.L_T[i] += L_T
                self.dstrain[i] += 0.5*(L + L_T)*P.dt
        self.strain += self.dstrain

        exec(P.S.law + '(self,P,G)')

        # self.stress += self.dstress # Cauchy stress rate
        for i in range(P.S.n): # FIXME - NEED TO REMOVE THIS LOOP
            self.stress[i] += self.dstress[i] - (np.dot(self.L[i],self.stress[i]) + np.dot(self.stress[i],self.L_T[i]))*P.dt # Lie stress rate

        if P.B.cyclic_lr: G.make_cyclic(P,G,['pressure', 'gammadot', 'dev_stress', 'yieldfunction', 'mu_s', 'mu', 'I'])

    def get_nodal_forces(self,P,G):
        """For every particle, get nodal forces. If cyclic boundaries are active, update the nodal cyclic forces.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        b = np.array([P.g*np.sin(P.theta),P.g*np.cos(P.theta),0.]) # body force
        bn = np.repeat(np.expand_dims(b,0),P.S.n,axis=0) # body force copied for each MP

        for r in range(4):
            for i in range(P.G.ns): G.fi[self.n[:,r]] += self.V*np.dot(self.stress[i],self.G[i,r]) # FIXME
            for j in range(2): G.fe[self.n[:,r],j] += self.N[:,r]*self.m*bn[:,j]

        if P.B.cyclic_lr: G.make_cyclic(P,G,['fi','fe'])

    def move_material_points(self,P,G):
        """For every particle, update the particle location.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """
        for r in range(4):
            for i in range(2):
                self.x[:,i] += P.dt*self.N[:,r]*G.q[self.n[:,r],i]    /G.m[self.n[:,r]]*G.valid_nodes[self.n[:,r]]
                self.v[:,i] += P.dt*self.N[:,r]*G.q_dot[self.n[:,r],i]/G.m[self.n[:,r]]*G.valid_nodes[self.n[:,r]]

        # HACK!!!!! THIS SHOULDN'T BE HERE!
        self.rho = 0.
        for r in range(4):
            self.rho += self.N[:,r]*G.m[self.n[:,r]]/G.V[self.n[:,r]]

    def get_nodal_gsd(self,P,G):
        """Measure the volume weighted average size at the node.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance.

        """

        for r in range(4):
            G.s_bar[self.valid_nodes[:,r]] += self.N[:,r]*self.m*self.size # mass weighted size
        if P.segregate:
            G.s_bar /= G.m

    def move_grainsize_on_grid(self,P,G):
        """Advect grainsize distribution across grid.

        :param P: A param.Param instance.

        """
        for r in range(4):
            self.phi[self.valid_nodes[:,r]] += self.N[self.valid_nodes[:,r],r]*G.dphi[self.valid_nodes[:,r]]
            self.pk[self.valid_nodes[:,r]]  += self.N[self.valid_nodes[:,r],r]*G.dpk[self.valid_nodes[:,r]]

        # re-normalise phi to account for `rounding erors`
        # if P.normalise_phi: # HACK: CHEATING!!!!!!!
        #     for p in range(P.phases):
        #         for i in range(P.S[p].n):
        #             self.S[p][i].phi[self.S[p][i].phi<0] = 0
        #             self.S[p][i].phi[self.S[p][i].phi>1] = 1
        #             self.S[p][i].phi /= sum(self.S[p][i].phi)
