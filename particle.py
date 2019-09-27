import sys
from numpy import *
from numpy.linalg import norm
from constit import *

class Particle():
    """A class defining an individual material point together with functions for mapping to/from the grid"""
    def __init__(self,x,y,P,G,p,phi=None):
        self.rho = P.S[p].rho # density (kg/m^3)
        self.x = array([x, y, 0.]) # position (m)
        try:
            self.v = P.S[p].initial_v
        except:
            self.v = array([0.,0.,0.]) # velocity (m/s)
        self.V = P.S[p].A*P.G.thickness # volume (m^3)
        self.m = self.rho*self.V # mass (kg)
        self.b = array([P.g*sin(P.theta),P.g*cos(P.theta),0.]) # body forces
        if P.pressure == 'lithostatic':
            K_0 = 1.
            if hasattr(P.G, 'top_gap'):
                # self.stress = cos(P.theta)*self.rho*P.max_g*(P.G.y_M - P.G.top_gap - y)*eye(3) # build in a pressure gradient
                sigma_yy = cos(P.theta)*self.rho*P.max_g*(P.G.y_M - P.G.top_gap - y) # build in a pressure gradient
                sigma_xy = sin(P.theta)*self.rho*P.max_g*(P.G.y_M - P.G.top_gap - y) # build in a pressure gradient
                sigma_xx = K_0*sigma_yy # build in a pressure gradient
                sigma_zz = K_0*sigma_yy # build in a pressure gradient
            else:
                # self.stress = cos(P.theta)*self.rho*P.max_g*(P.G.y_M - y)*eye(3) # build in a pressure gradient
                sigma_yy = cos(P.theta)*self.rho*P.max_g*(P.G.y_M - y) # build in a pressure gradient
                sigma_xy = sin(P.theta)*self.rho*P.max_g*(P.G.y_M - y) # build in a pressure gradient
                sigma_xx = K_0*sigma_yy # build in a pressure gradient
                sigma_zz = K_0*sigma_yy # build in a pressure gradient
            self.stress = array([[sigma_xx,sigma_xy,0],[sigma_xy,sigma_yy,0],[0,0,sigma_zz]])
            # self.strain = zeros((3,3)) # strain-free
            self.strain = eye(3)*trace(self.stress)/(9.*P.S[p].K) + (self.stress - eye(3)*trace(self.stress)/3.)/(2.*P.S[p].G)
        elif P.pressure == 'non-zero':
            self.strain = -1e-15*ones((3,3)) # mixed state compression
            self.stress = (P.S[p].K*trace(self.strain)*eye(3) +
                           2.*P.S[p].G*(self.strain - trace(self.strain)*eye(3)/3.))
        elif P.pressure == 'temperature':
            self.strain = -1e-5*eye(3) # pure compression everywhere
            self.stress = (P.S[p].K*trace(self.strain)*eye(3) +
                           2.*P.S[p].G*(self.strain - trace(self.strain)*eye(3)/3.))
        elif P.pressure == 'compression':
            self.stress = array([[-P.initial_pressure,0,0],
                                 [0,-P.initial_pressure,0],
                                 [0,0,-P.initial_pressure]])
            # self.strain = zeros((3,3)) # strain-free
            self.strain = eye(3)*trace(self.stress)/(9.*P.S[p].K) + (self.stress - eye(3)*trace(self.stress)/3.)/(2.*P.S[p].G)
        else:
            self.strain = zeros((3,3)) # corresponding strains
            self.stress = zeros((3,3))
        self.dstrain = zeros((3,3)) # change in strain
        self.dstress = zeros((3,3)) # change in strain
        self.de_ij = zeros((3,3))
        self.sigma_kk = trace(self.stress)
        self.dev_stress = self.stress - self.sigma_kk*eye(3)/3.
        self.gammadot = 0.
        self.pressure = self.sigma_kk/3.
        self.sigmav = self.stress[1,1]
        self.sigmah = self.stress[0,0]
        self.yieldfunction = 0.
        self.G = zeros((4,3)) # gradient of shape function
        self.N = zeros((4)) # basis function
        self.n_star = 0 # reference node
        self.pk = 0. # kinetic pressure
        if phi is not None:
            self.phi = array(phi)
        else:
            self.phi = array(P.S[p].phi)
        if P.initial_flow == 'steady':
            if P.S[0].law == 'viscous' or P.S[0].law == 'viscous_size':
                self.v = array([self.rho*P.max_g*sin(P.theta)*(P.G.y_M*y - y**2/2.)/P.S[p].mu_s,0.,0.]) # VISCOUS FLOW DOWN A SLOPE
            elif P.S[0].law == 'pouliquen' or P.S[0].law == 'pouliquen2D':
                s_bar = (P.G.s_m+P.G.s_M)/2.
                self.v = array([sqrt(abs(P.max_g)*P.G.s_bar_0)*(2./3.)*P.S[p].I_0*
                               (tan(abs(P.theta))-P.S[p].mu_0)/(P.S[p].mu_1-tan(abs(P.theta)))*
                               sqrt(P.S[p].packing*cos(abs(P.theta)))*
                               (P.G.y_M**1.5 - (P.G.y_M - y)**1.5)/P.G.s_bar_0**1.5,
                               0.,0.])
        elif P.initial_flow == 'steady_poiseulle': # https://link.springer.com/content/pdf/10.1007%2Fs10035-013-0447-3.pdf
            H = P.G.y_M - P.G.y_m
            y_star = y/H
            if y_star > 0.5: y_star = 0.5 - abs(y_star - 0.5) # y_star only defined on LHS
            K = abs(P.S[0].rho*P.max_g)
            Lambda = H/P.G.s_bar_0
            beta = sqrt(abs(K*H/P.initial_pressure))
            y_star_c = 0.5 - P.S[p].mu_0/(beta**2)
            if y_star < y_star_c:
                u_star = (Lambda*P.S[p].I_0)/(beta**3)*\
                        ((beta**2)*y_star + (P.S[p].mu_0 - P.S[p].mu_1)*log(-(2.*(beta**2)*y_star + 2*P.S[p].mu_1 - beta**2)/(beta**2 - 2*P.S[p].mu_1)))
            else:
                u_star = (Lambda*P.S[p].I_0)/(beta**3)*\
                         ((beta**2)/2. - P.S[p].mu_0 + (P.S[p].mu_0 - P.S[p].mu_1)*log((P.S[p].mu_1-P.S[p].mu_0)/(P.S[p].mu_1-(beta**2)/2.)))
            self.v = array([-u_star*sqrt(K*H/P.S[p].rho_s),0.,0.])

            stable = beta<sqrt(2*P.S[p].mu_1) # if True, this should work
            if not stable:
                print('WARNING: THIS SETUP WILL ACCELERATE FOREVER')
                print(sqrt(2*P.S[p].mu_0),beta,sqrt(2*P.S[p].mu_1))


        # elif P.initial_flow == 'steady_poiseulle': # TOTALLY FAKE
        #         H = P.G.y_M - P.G.y_m
        #         y_star = y/H
        #         if y_star > 0.5: y_star = 0.5 - abs(y_star - 0.5) # y_star only defined on LHS
        #         K = abs(P.S[0].rho*P.max_g)
        #         Lambda = H/P.G.s_bar_0
        #         beta = sqrt(abs(K*H/P.initial_pressure))
        #         y_star_c = 0.3
        #         # LINEAR
        #         # if y_star < y_star_c: u_star = -1.0/y_star_c*y_star
        #         # else:                 u_star = -1.0
        #         # PARABOLA
        #         # if y_star < y_star_c: u_star = 1.0*((y_star - 0.5)**2 - 1./4.)
        #         # else:                 u_star = 1.0*((y_star_c - 0.5)**2 - 1./4.)
        #         # CUBIC
        #         if y_star < y_star_c: u_star = 1.0*((y_star - 0.5)**4 - 0.5**4)
        #         else:                 u_star = 1.0*((y_star_c - 0.5)**4 - 0.5**4)
        #
        #         self.v = array([-u_star,0.,0.])

            # elif P.S[0].law == 'marks2012': self.v = array([self.rho*P.max_g*sin(P.theta)*(P.G.y_M*y - y**2/2.)/P.S[p].mu_s,0.,0.])
        # print(self.stress)
    def update_strain(self,P,G):
        """Calculate incremental strain at the particle level from the nodal velocity. Also update the total strain.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        self.dstrain = zeros((3,3))
        self.L = self.L_T = zeros((3,3))
        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                v = G.q[n]/G.m[n]
            else:
                v = array([0.,0.,0.])
            L = outer(v,self.G[r])
            L_T = outer(self.G[r],v)
            self.L += L
            self.L_T += L_T
            self.dstrain += 0.5*(L + L_T)*P.dt
        self.strain += self.dstrain

    def update_stress(self,P,G,p):
        """Calculate incremental stress at the particle level using the relevant constitutive model.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance
        :param p: Which material phase

        """
        exec(P.S[p].law + '(self,P,G,p)')
        # self.stress += self.dstress # Cauchy stress rate
        self.stress += self.dstress - (dot(self.L,self.stress) + dot(self.stress,self.L_T))*P.dt # Lie stress rate

    def get_reference_node(self,P,G):
        """Find the reference node for this material point.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        try:
            self.n_star = int(trunc((self.x[0] - P.G.x_m)/G.dx)
                            + trunc((self.x[1] - P.G.y_m)/G.dy)*P.G.nx)
        except ValueError:
            print('Something went wrong! Particle is out of the box :(')
            print(self.x)
            print(self.v)


    def get_basis_functions(self,P,G):
        """Get the basis functions for this material point.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        self.N, self.G = G.N(self.x-[G.X[self.n_star],G.Y[self.n_star],0.]) # linear basis function

    def get_nodal_mass_momentum(self,P,G):
        """Pass mass and momentum to nearby nodes.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            G.m[n] += self.N[r]*self.m
            G.q[n] += self.N[r]*self.m*self.v
            if P.segregate_grid:
                G.phim[n] += self.N[r]*self.m*self.phi
                G.pkm[n]  += self.N[r]*self.m*self.pk

    def add_to_boundary(self,P,G):
        """I have no idea what this does.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        mask = zeros_like(G.boundary_v)
        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            mask[n] = 1
        return mask

    def update_nodal_gsd(self,P,G):
        """Calculate the grainsize distribution at a node.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                G.s_bar[n] += self.N[r]*self.m*self.size # mass weighted size
#         for position,s in enumerate(P.GSD.s):
#             if s == self.s:
#                 pos = position
#         for r in range(4):
#             n = G.nearby_nodes(self.n_star,r,P)
#             G.gsd[n,pos] += self.N[r]*self.m/self.rho # total volume in size bin

    def map_gsd_to_mp(self,P,G):
        """Pass the nodal grainsize information back to the material point.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        self.s_bar = 0.
        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                self.s_bar += self.N[r]*G.s_bar[n]

    def map_phi_to_mp(self,P,G):
        """Pass the nodal grainsize information back to the material point.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """
        self.pk = 0.
        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                self.phi += self.N[r]*G.dphi[n]
                self.pk  += self.N[r]*G.pk[n]

    def move_material_points(self,P,G):
        """Update material point position and velocity using nearby nodal mass and momentum.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                self.x += P.dt*self.N[r]*G.q[n]/G.m[n]
                self.v += P.dt*self.N[r]*G.q_dot[n]/G.m[n]

        # HACK!!!!! THIS SHOULDN'T BE HERE!
        self.rho = 0.
        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            if G.m[n] > P.M_tol:
                self.rho += self.N[r]*G.m[n]/G.V[n]

    def move_rigid_points(self,P,G):
        """Update material point position and velocity for rigid particles with prescribed velocity.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        self.x += P.dt*P.intruder_v # Move at the defined rate
        self.v = P.intruder_v

    def recover_position(self,P,G):
        """Find out where this material point should be.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        x = zeros((3,))
        nodes = zeros((4))
        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            x += self.N[r]*array([G.X[n],G.Y[n],0.])
            nodes[r] = n
        return x, nodes

    def get_nodal_forces(self,P,G):
        """Calculate nearby nodal internal forces and body forces. Tractions are ignored for now.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        for r in range(4):
            n = G.nearby_nodes(self.n_star,r,P)
            G.fi[n] += self.V*dot(self.stress,self.G[r])
            G.fe[n] += self.N[r]*self.m*self.b
            # if there are tractions add something here

    def increment_gravity(self,P,G,p):
        """Calculate gravity components for this material point.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        if p == 0:
            self.b = array([P.g*sin(P.theta),P.g*cos(P.theta),0.]) # body forces

    def energy(self,P):
        """Calculate energies of this material point for plotting purposes.

        :param P: A param.Param instance.
        :param G: A grid.Grid instance

        """

        k = 0.5*self.m*sum(self.v**2) # kinetic
        h = (sum(self.x**2))**0.5*sin(P.theta+arctan(abs(self.x[1]/self.x[0]))) # height
        g = -self.m*P.g*h # gravitational potential
        s = 0.5*sum(sum(self.stress*self.strain))*self.V # elastic potential energy
        return [k,g,s]

#     def find_f(self,P,L):
#         """Calculate segregation scaling function.
#
#         .. warning::
#
#             Don't use this.
#
#         :param P: A param.Param instance.
#         :param G: A grid.Grid instance
#
#         """
#
#         self.f = P.GSD.s/(sum(P.GSD.s*self.phi))
#         self.rho_bar = sum(self.rho*self.phi)
#
#     def find_seg_velocity(self,P,L):
#         """Calculate segregation velocity.
#
#         .. warning::
#
#             Don't use this.
#
#         :param P: A param.Param instance.
#         :param G: A grid.Grid instance
#
#         """
#         self.u = zeros((len(self.neighbours),3))
#         for k in range(len(self.neighbours)):
#             j = int(self.neighbours[k])
#             n = (L.S[j].x - self.x)
#             # find seg velocity u
#             r = self.phi*self.rho*P.g
#             s = self.phi*dot(L.S[j].stress*L.S[j].f - self.stress*self.f,n)
#             m = dot(L.S[j].m*L.S[j].v - self.m*self.v,n) # material derivative of momentum
#             self.u[k] = self.gammadot/(self.rho_bar*self.phi*P.S.c)*(r - s - m) # seg velocity
#             print(self.u[k])
#
#     def find_flux(self,P,L):
#         """Calculate flux to move between material points.
#
#         .. warning::
#
#             Don't use this.
#
#         :param P: A param.Param instance.
#         :param G: A grid.Grid instance
#
#         """
#         self.dQ = zeros((len(self.neighbours)))
#         self.dm = zeros((len(self.neighbours)))
#         for k in range(len(self.neighbours)):
#             j = int(self.neighbours[k])
#             if j >= 0:
#                 n = (L.S[j].x - self.x)
#                 d = norm(n) # distance
#                 n_hat = n/d
#                 # use seg velocity as advection term below
#                 self.dQ[k] = -(P.S.D/d*(self.m - L.S[j].m) - dot(self.m*self.u[k] - L.S[j].m*L.S[j].u[0],n_hat)) # advection of seg velocity only
#                 self.dQ[k] = min(self.m/P.dt/2.,max(0.,self.dQ[k])) # min and max bounds on flux
#                 self.dm[k] = float(self.areas[k])*self.dQ[k]*P.dt
#
#     def move_flux(self,P,L):
#         """ For moving grainsize between material points.
#
#         .. warning::
#
#             Don't use this.
#
#         :param P: A param.Param instance.
#         :param G: A grid.Grid instance
#
#         """
#         for k in range(len(self.neighbours)):
#             j = int(self.neighbours[k])
#             if j >= 0:
#                 self.m += self.dm[k]
#                 L.S[j].m -= self.dm[k]
#                 #print self.rho
