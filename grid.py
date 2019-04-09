from numpy import *
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from integrators import increment_grainsize
from scipy.interpolate import RectBivariateSpline as interp2d
from astropy.convolution import convolve, Gaussian2DKernel, Box2DKernel
import sys

# from numba import jitclass          # import the decorator
# from numba import int32, float32, bool_    # import the types
#
# spec = [
#     ('x', float32[:]),          # an array field
#     ('y', float32[:]),          # an array field
#     ('dx', float32),          # an array field
#     ('dy', float32),          # an array field
#     ('X', float32[:]),          # an array field
#     ('Y', float32[:]),          # an array field
#     ('x_plot', float32[:]),          # an array field
#     ('y_plot', float32[:]),          # an array field
#     ('boundary_v', bool_[:]),          # an array field
#     ('boundary_h', bool_[:]),          # an array field
#     ('boundary_tot', bool_[:]),          # an array field
#     ('weighting', float32[:]),          # an array field
#     ('V', float32[:]),          # an array field
#     ('boundary_conveyor', bool_[:]),          # an array field
#     ('boundary_conveyor_triple', bool_[:]),          # an array field
#     ('fi', float32[:]),          # an array field
#     ('fe', float32[:]),          # an array field
#     ('m', float32[:]),          # an array field
#     ('q', float32[:]),          # an array field
#     ('q_dot', float32[:]),          # an array field
#     ('ext_f', float32[:]),          # an array field
#     ('gammadot', float32[:]),          # an array field
#     ('grad_gammadot', float32[:]),          # an array field
#     ('yieldfunction', float32[:]),          # an array field
#     ('pressure', float32[:]),          # an array field
#     ('sigmah', float32[:]),          # an array field
#     ('sigmav', float32[:]),          # an array field
#     ('dev_stress', float32[:]),          # an array field
#     ('dev_stress_dot', float32[:]),          # an array field
#     ('mu_s', float32[:]),          # an array field
#     ('mu', float32[:]),          # an array field
#     ('I', float32[:]),          # an array field
#     ('damping_force', float32[:]),          # an array field
#     ('s_bar', float32[:]),          # an array field
#     ('dphi', float32[:]),          # an array field
#     ('phim', float32[:]),          # an array field
#     ('phi', float32[:]),          # an array field
#     ('S', float32[:]),          # an array field
#     ('u_hat', float32[:]),          # an array field
#     ('v_hat', float32[:]),          # an array field
# ]
#
# @jitclass(spec)
class Grid():
    """
    This contains all of the methods which operate on the grid directly. The grid is assumed to be a regular lattice, numbered as:

    +----------+-----+-----+-----+-----------+
    |(ny-1)*nx | ... | ... | ... | (ny*nx-1) |
    +----------+-----+-----+-----+-----------+
    |...       | ... | ... | ... |   ...     |
    +----------+-----+-----+-----+-----------+
    |2*nx      | ... | ... | ... | (3*nx-1)  |
    +----------+-----+-----+-----+-----------+
    |nx        | ... | ... | ... | (2*nx-1)  |
    +----------+-----+-----+-----+-----------+
    |0         |  1  |  2  | ... |   nx-1    |
    +----------+-----+-----+-----+-----------+

    """
    def __init__(self,P):
        """
        Pass the parameter values to the grid on initialisation. Also set up boundaries and define the volume of every node.

        Parameters
        ----------
        P : Params class

        """
        # arrays to be called for actual locations of grid points
        self.x = linspace(P.G.x_m,P.G.x_M,P.G.nx)
        self.y = linspace(P.G.y_m,P.G.y_M,P.G.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.X = tile(self.x,P.G.ny)
        self.Y = repeat(self.y,P.G.nx)
        self.boundary(P)
        self.volume(P)

        # arrays to be used for plotting with pcolormesh - needs cell edges
        self.x_plot = hstack([P.G.x_m - self.dx/2.,(self.x[1:] + self.x[:-1])/2.,P.G.x_M + self.dx/2.])
        self.y_plot = hstack([P.G.y_m - self.dy/2.,(self.y[1:] + self.y[:-1])/2.,P.G.y_M + self.dy/2.])

    def boundary(self,P):
        """
        Build the boundary.

        This consists of three masks, the same shape as the grid, `boundary_v`, `boundary_h` and `boundary_tot`. These are boundaries that restrain flow in the left/right and up/down directions, respectively. `boundary_tot` is the sum of the two other boundaries.

        This respects the following arguments which may be present in the input file:
            * has_bottom
            * has_top
            * has_left
            * has_right
            * no_slip_bottom
            * wall
            * box
        """
        self.boundary_v = zeros((P.G.ny*P.G.nx),dtype=bool)
        self.boundary_h = zeros((P.G.ny*P.G.nx),dtype=bool)
        if P.B.has_bottom: self.boundary_h[:P.G.nx] = 1 # bottom
        if P.B.no_slip_bottom: self.boundary_v[:P.G.nx] = 1 # bottom
        if P.B.has_top: self.boundary_h[(P.G.ny-1)*P.G.nx:] = 1 # top
        if P.B.has_left: self.boundary_v[::P.G.nx] = 1 # left
        if P.B.has_right: self.boundary_v[P.G.nx-1::P.G.nx] = 1 # right
        if P.B.wall: self.boundary_v[(P.G.nx-1)//2:(P.G.ny*P.G.nx-1)//2:P.G.nx] = 1 # wall at centre
        if P.B.silo_left:
            self.boundary_v[1:(P.G.ny//2-3)*P.G.nx:P.G.nx] = 1 # wall at outlet
            self.boundary_v[(P.G.ny//2+3)*P.G.nx+1:-1:P.G.nx] = 1 # wall at outlet
        if P.B.silo_bottom:
            self.boundary_h[:P.G.nx//2-2] = 1 # wall at outlet
            self.boundary_h[P.G.nx//2+3:] = 1 # wall at outlet
        if P.B.box:
            l = 4 # min grid points in from left
            r = 8 # max grid points across from left
            b = 4 # min grid points up from bottom
            t = 8 # max grid points up from bottom
            self.boundary_v[b*P.G.ny + l:t*(P.G.ny+1):P.G.nx] = 1 # left wall
            self.boundary_v[b*P.G.ny + r:t*(P.G.ny+2):P.G.nx] = 1 # right wall
            self.boundary_h[b*P.G.ny + l:b*P.G.ny + r+1] = 1 # bottom wall
            self.boundary_h[t*P.G.ny + l:t*P.G.ny + r+1] = 1 # top wall
        if P.B.conveyor:
            self.boundary_conveyor = zeros((P.G.ny*P.G.nx),dtype=bool)
            self.boundary_conveyor_triple = zeros((P.G.ny*P.G.nx,3),dtype=bool)
            self.boundary_conveyor[2:P.G.nx-2] = 1
            self.boundary_conveyor_triple[:,0] = self.boundary_conveyor

        self.boundary_tot = self.boundary_v + self.boundary_h

    def volume(self,P): # get volume of each cell
        self.weighting = ones((P.G.ny*P.G.nx))
        self.weighting[:P.G.nx] /=2 # bottom
        self.weighting[(P.G.ny-1)*P.G.nx:] /= 2 # top
        if not P.B.cyclic_lr:
            self.weighting[::P.G.nx] /= 2 # left
            self.weighting[P.G.nx-1::P.G.nx] /= 2 # right
        self.V = self.weighting*self.dx*self.dy*P.G.thickness

    def wipe(self,P):
        self.fi = zeros((P.G.nx*P.G.ny,3)) # nodal internal force - {x,y,z}
        self.fe = zeros((P.G.nx*P.G.ny,3)) # nodal external force - {x,y,z}
        self.m = zeros((P.G.nx*P.G.ny)) # nodal mass
        self.q = zeros((P.G.nx*P.G.ny,3)) # nodal momentum
        self.q_dot = zeros((P.G.nx*P.G.ny,3)) # nodal change in momentum
        self.ext_f = zeros((P.G.nx*P.G.ny,3))

        self.gammadot = zeros((P.G.nx*P.G.ny)) # shear strain rate
        self.grad_gammadot = zeros((P.G.nx*P.G.ny,3)) # gradient of shear strain rate
        self.yieldfunction = zeros((P.G.nx*P.G.ny)) # yield function
        self.pressure = zeros((P.G.nx*P.G.ny)) # isotropic tension
        self.sigmah = zeros((P.G.nx*P.G.ny))
        self.sigmav = zeros((P.G.nx*P.G.ny))
        self.dev_stress = zeros((P.G.nx*P.G.ny)) # deviatoric stress norm
        self.dev_stress_dot = zeros((P.G.nx*P.G.ny)) # incremental deviatoric stress norm
        self.mu_s = zeros_like(self.m) # shear viscosity, for plotting with viscous_size
        self.mu = zeros_like(self.m)
        self.I = zeros_like(self.m)
        self.eta = zeros_like(self.m)

        self.damping_force = zeros((P.G.nx*P.G.ny,3)) # local non-viscous damping
        self.s_bar = zeros([P.G.nx*P.G.ny])
        self.dphi = zeros([P.G.nx*P.G.ny,P.G.ns])

        if P.segregate_grid:
            self.phim = zeros([P.G.nx*P.G.ny,P.G.ns])
            self.phi = zeros([P.G.nx*P.G.ny,P.G.ns])
            self.S = tile(P.G.s,[P.G.nx*P.G.ny,1])
            self.u_hat = zeros([P.G.nx*P.G.ny,P.G.ns])
            self.v_hat = zeros([P.G.nx*P.G.ny,P.G.ns])

    def nearby_nodes(self,n_star,r,P):
        if (r == 0) or (r == 1):
            return n_star + r
        else:
            return n_star + P.G.nx + 3-r

    def N(self,x):
        '''
        Calculate shape functions and their derivatives for Lagrange 4-node system
        '''
        N = zeros((4))
        G = zeros((4,3)) # {x,y,z}
        # Lagrange 4-node system
        N[0] = (self.dx-x[0])*(self.dy-x[1])/(self.dx*self.dy)
        N[1] = (x[0])*(self.dy-x[1])/(self.dx*self.dy)
        N[2] = (x[0])*(x[1])/(self.dx*self.dy)
        N[3] = (self.dx-x[0])*(x[1])/(self.dx*self.dy)
        G[0] = array([-(self.dy-x[1]),-(self.dx-x[0]),0.])/(self.dx*self.dy)
        G[1] = array([(self.dy-x[1]),-x[0],0.])/(self.dx*self.dy)
        G[2] = array([x[1],x[0],0.])/(self.dx*self.dy)
        G[3] = array([-x[1],(self.dx-x[0]),0.])/(self.dx*self.dy)
        return N, G

#     def apply_cyclic_BCs(self,P): # DOES NOT WORK !!!!
#         for param in [self.m, self.q]:#, self.s_bar]:#, self.q_dot]:
#             temp_store = param[::P.G.nx].copy() # left boundary
#             param[::P.G.nx] += param[P.G.nx-1::P.G.nx].copy() # add right to left boundary
#             param[P.G.nx-1::P.G.nx] += temp_store
# #         if P.B.cyclic_lr: self.s_bar /= self.m
#
#     def apply_cyclic_BCs_stress(self,P): # DOES NOT WORK !!!!
#         for param in [self.pressure, self.gammadot, self.dev_stress, self.yieldfunction]:
#             temp_store = param[::P.G.nx].copy() # left boundary
#             param[::P.G.nx] += param[P.G.nx-1::P.G.nx].copy() # add right to left boundary
#             param[P.G.nx-1::P.G.nx] += temp_store
#
#     def apply_cyclic_nodal_forces(self,P): # DOES NOT WORK !!!!
#         for param in [self.fi, self.fe]:
#             temp_store = param[::P.G.nx].copy() # left boundary
#             param[::P.G.nx] += param[P.G.nx-1::P.G.nx].copy() # add right to left boundary
#             param[P.G.nx-1::P.G.nx] += temp_store

    def make_cyclic(self,P,G,params): # DOES NOT WORK !!!!
        for label in params:
            param = getattr(G, label) # get the right attribute
            temp_store = param[::P.G.nx].copy() # left boundary
            param[::P.G.nx] += param[P.G.nx-1::P.G.nx].copy() # add right to left boundary
            param[P.G.nx-1::P.G.nx] += temp_store # add left to right boundary

    def BCs(self,P): # BCS directly affecting self.fe
        if P.B.vertical_force:
            self.ext_f[:P.G.nx,1] += P.q_v # bottom
            self.ext_f[(P.G.ny-1)*P.G.nx:,1] -= P.q_v # top
            self.fe[:,1] += 2.*self.ext_f[:,1]*self.m/P.S[0].rho/self.dy
        if P.B.horizontal_force:
            self.ext_f[::P.G.nx,0] += P.q_h # left
            self.ext_f[P.G.nx-1::P.G.nx,0] -= P.q_h # right
            self.fe[:,0] += 2.*self.ext_f[:,0]*self.m/P.S[0].rho/self.dx

        if P.mode == 'anisotropy' and P.t == 0:
            self.fe[P.G.nx*P.G.ny/2,2] = 1.

    def update_momentum(self,P):
        """Update the momentum at every nodal location, as calculated from the nodal change in moment, :math:`\dot q`. This also...

        :param P: A param.Param instance.

        """
        # if P.damping: self.damping_force = 0.8*abs(self.fe - self.fi)*sign(self.q_dot)
        self.q_dot = self.fe - self.fi #- self.damping_force
        if P.damping: self.q_dot *= 0.7

        self.q_dot[:,0] = self.q_dot[:,0]*(1.-self.boundary_v) # 0 at boundary
        self.q_dot[:,1] = self.q_dot[:,1]*(1.-self.boundary_h) # 0 at boundary
        if P.B.roughness:
                self.q_dot[:,0] = self.q_dot[:,0]*(1.-self.boundary_h) # top/bottom
                self.q_dot[:,1] = self.q_dot[:,1]*(1.-self.boundary_v) # sidewalls
        self.q += self.q_dot*P.dt

        if P.B.conveyor: self.q[self.boundary_conveyor_triple] = P.v_0*self.m[self.boundary_conveyor] # conveyor momentum

    def calculate_gammadot(self,P,G,smooth=False):
        """Calculate the bulk shear strain rate from the continuum measure of velocity.

        :param P: A param.Param instance.

        """
        u = (self.q[:,0]/self.m)#.reshape(P.G.ny,P.G.nx)
        v = (self.q[:,1]/self.m)#.reshape(P.G.ny,P.G.nx)
        # u[isnan(u)] = 0.
        # v[isnan(v)] = 0.
        # dudy,dudx = gradient(u,G.dy,G.dx)
        # dvdy,dvdx = gradient(v,G.dy,G.dx)
        gradu = self.calculate_gradient(P,G,u,smooth=P.smooth_gamma_dot)
        gradv = self.calculate_gradient(P,G,v,smooth=P.smooth_gamma_dot)
        dudy = gradu[:,1]
        dvdx = gradv[:,0]
        self.gammadot = (dudy + dvdx)#.flatten()
        # self.rate_of_shear = array([])

    def calculate_grad_gammadot(self,P,G):
        """Calculate the gradient of the absolute value of the shear strain rate using the built-in gradient method.

        :param P: A param.Param instance.

        """
        self.grad_gammadot = G.calculate_gradient(P,G,abs(G.gammadot),smooth=P.smooth_grad2)


    def calculate_gradient(self,P,G,Z,smooth=False):
        """Calculate the gradient of any property. Deals with grid points that have no mass (that should'nt contribute to the gradient).

        :param P: A param.Param instance.

         .. warning::

             This returns the gradient of the input field ... kind of ... sometimes
        """

        Z = ma.masked_where(G.m<P.M_tol,Z).reshape(P.G.ny,P.G.nx)
        dZdy,dZdx = gradient(Z,G.dy,G.dx)
        if smooth: # For details of astropy convolution process, see here: http://docs.astropy.org/en/stable/convolution/using.html
            # kernel = Box2DKernel(smooth) # smallest possible square kernel is 3
            kernel = Gaussian2DKernel(x_stddev=1,y_stddev=1)
            dZdy = convolve(dZdy, kernel, boundary='extend')
            dZdx = convolve(dZdx, kernel, boundary='extend')

        grad = array([dZdx.flatten(),
                      dZdy.flatten(),
                      zeros_like(G.m)]).T
        return grad

    def calculate_phi_increment(self,P):
        """Calculate the incremental change in phi across the grid.

        :param P: A param.Param instance.

         .. warning::

             This returns the gradient of the ABSOLUTE VALUE of the shear strain rate
        """
        self.dphi = increment_grainsize(P,self)
        self.dphi = nan_to_num(self.dphi)
