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

        # self.dx_plot = self.x_plot[1:] - self.x_plot[:-1]
        # self.dy_plot = self.y_plot[1:] - self.y_plot[:-1]
        self.dx_plot = hstack([self.dx/2.,self.dx*ones(P.G.nx-2),self.dx/2.])
        self.dy_plot = hstack([self.dy/2.,self.dy*ones(P.G.ny-2),self.dy/2.])

        # if P.B.cyclic_lr: self.DX = P.G.dx*ones([P.G.nx-1])
        # else: self.DX = hstack([P.G.dx/2.,P.G.dx*ones([P.G.nx-3]),P.G.dx/2.])
        # self.DY = hstack([P.G.dy/2.,P.G.dy*ones([P.G.ny-3]),P.G.dy/2.])

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
        if P.B.two_walls:
            self.boundary_v[(P.G.nx-1)//4:(P.G.ny*P.G.nx-1)//2:P.G.nx] = 1 # low wall at left
            self.boundary_v[(P.G.ny//2)*P.G.nx+(3*(P.G.nx-1))//4::P.G.nx] = 1 # upper wall at right
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

        self.pk = zeros([P.G.nx*P.G.ny])
        self.dpk = zeros([P.G.nx*P.G.ny])
        self.pkm = zeros([P.G.nx*P.G.ny])
        self.grad_pk = zeros((P.G.nx*P.G.ny,3))

        if P.segregate_grid:
            self.phim = zeros([P.G.nx*P.G.ny,P.G.ns])
            self.phi = zeros([P.G.nx*P.G.ny,P.G.ns])
            self.S = tile(P.G.s,[P.G.nx*P.G.ny,1])
            self.u_hat = zeros([P.G.nx*P.G.ny,P.G.ns])
            self.v_hat = zeros([P.G.nx*P.G.ny,P.G.ns])

    def nearby_nodes(self,n_star,r,P): # THIS VERSION WAS WORKING
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

    def make_cyclic(self,P,G,params):
        """
        WORKS (EXCEPT MAYBE WHEN nx < 3) !!!!
        """
        for label in params:
            param = getattr(G, label) # get the right attribute
            temp_store = param[::P.G.nx].copy() # left boundary
            param[::P.G.nx] += param[P.G.nx-1::P.G.nx].copy() # add right to left boundary
            param[P.G.nx-1::P.G.nx] += temp_store # add left to right boundary

    def BCs(self,P):
        """
        Boundary conditions are applied directly to the external force G.fe
        """
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

        # Impose orthogonal BCs
        if P.B.roughness:
            if P.B.wall_mu:
                bottom = self.boundary_h[:P.G.nx] = 1 # bottom
                top = self.boundary_h[(P.G.ny-1)*P.G.nx:] = 1 # top
                left = self.boundary_v[::P.G.nx] = 1 # left
                right = self.boundary_v[P.G.nx-1::P.G.nx] = 1 # right

                Ff_h = self.boundary_h*minimum(abs(self.q_dot[:,0]),P.B.wall_mu*abs(self.q_dot[:,1])) # min of tangential force and mu*normal force
                Ff_v = self.boundary_v*minimum(abs(self.q_dot[:,1]),P.B.wall_mu*abs(self.q_dot[:,0]))

                # JUST OPERATE WHEN NORMAL FORCE IS POINTING TOWARDS THE BOUNDARY
                Ff_v[left]   *= self.q_dot[:,0][left]   < 0
                Ff_v[right]  *= self.q_dot[:,0][right]  > 0
                Ff_h[top]    *= self.q_dot[:,1][top]    > 0
                Ff_h[bottom] *= self.q_dot[:,1][bottom] < 0

                self.q_dot[:,0] -= sign(self.q_dot[:,0])*Ff_h#*activated_wall # opposite in sign to the applied tangential force
                self.q_dot[:,1] -= sign(self.q_dot[:,1])*Ff_v#*activated_wall
            else:
                self.q_dot[:,0] = self.q_dot[:,0]*(1.-self.boundary_h) # top/bottom
                self.q_dot[:,1] = self.q_dot[:,1]*(1.-self.boundary_v) # sidewalls

        # Impose normal BCs
        self.q_dot[:,0] = self.q_dot[:,0]*(1.-self.boundary_v) # 0 at boundary
        self.q_dot[:,1] = self.q_dot[:,1]*(1.-self.boundary_h) # 0 at boundary

        self.q += self.q_dot*P.dt

        if P.B.conveyor: self.q[self.boundary_conveyor_triple] = P.v_0*self.m[self.boundary_conveyor] # conveyor momentum

    def calculate_gammadot(self,P,G,smooth=False):
        """Calculate the bulk shear strain rate from the continuum measure of velocity.

        :param P: A param.Param instance.

        """
        u = (self.q[:,0]/self.m)
        v = (self.q[:,1]/self.m)
        gradu = self.calculate_gradient(P,G,u,smooth=smooth)
        gradv = self.calculate_gradient(P,G,v,smooth=smooth)
        dudy = gradu[:,1]
        dvdx = gradv[:,0]
        self.gammadot = (dudy + dvdx)

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
        #     # kernel = Box2DKernel(smooth) # smallest possible square kernel is 3
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

    def update_pk(self,P,G):
        """Placeholder method until Ebrahim's model is finished.

        :param P: A param.Param instance.

        """
        decay_time = 0.1 # seconds


        # diffusivity = (length_scale**2)/(2.*decay_time) # definition of diffusivity?

        # I NEED TO IMPLEMENT BOUNDARY CONDITIONS FOR DIFFUSION PART
        #   1. At a boundary, reflection boundary?
        #   2. At a periodic boundary, use the other side
        #   3. At a free surface, it should be high!....
        # pk = ma.masked_where(G.m<P.M_tol,self.pk).reshape(P.G.ny,P.G.nx)
        # pk_pad_x = hstack([pk[:,0,newaxis], pk,  pk[:,-1, newaxis]])
        # pk_pad_y = vstack([pk[newaxis,0,:], pk,  pk[newaxis,-1,:]])
        # d2pk_dx2 = (roll(pk_pad_x,1,axis=1) - 2*pk_pad_x + roll(pk_pad_x,-1,axis=1))/P.G.dx**2
        # d2pk_dy2 = (roll(pk_pad_y,1,axis=0) - 2*pk_pad_y + roll(pk_pad_y,-1,axis=0))/P.G.dy**2
        # diff_term = (d2pk_dx2[:,1:-1] + d2pk_dy2[1:-1,:]).flatten()

        # if D > 0:
        # grad_pk  = self.calculate_gradient(P,G,self.pk,smooth=False)
        # grad2_Dpk_dx = self.calculate_gradient(P,G,diffusivity*grad_pk[:,0],smooth=False)[:,0]
        # grad2_Dpk_dy = self.calculate_gradient(P,G,diffusivity*grad_pk[:,1],smooth=False)[:,1]
        # diff_term = grad2_Dpk_dx + grad2_Dpk_dy

        # self.dpk = nan_to_num(P.l*self.s_bar*sqrt(abs(self.pressure/self.V))*abs(self.gammadot) - self.pk)/decay_time*P.dt #-

        # LATEST AND BEST WORKING EVOLUTION EQUATION:
        # max_gamma_dot_allowable = 100.0
        # sanitised_gamma_dot = minimum(abs(nan_to_num(self.gammadot)),max_gamma_dot_allowable)
        # self.dpk = (P.l*self.s_bar**2*sanitised_gamma_dot**2*self.m/self.I - self.pk)/decay_time*P.dt # p_k_steady = l*gamma_dot^2*d^2/I
        # self.grad_pk = self.calculate_gradient(P,G,self.pk.copy(),smooth=False)

        # self.I = maximum(minimum(self.I/self.m,1.0),1e-6)*self.m
        # self.I[G.m<P.M_tol] = nan
        self.pk = nan_to_num(-(self.pressure/self.m)*(self.I/self.m)) # tension positive!!
        self.grad_pk = self.calculate_gradient(P,G,self.pk,smooth=False)

        # # JUST USED FOR SEGREGATION MODEL - NOT ACTUALLY GRAD OF PK!!!!
        # grad_pk_mag = sqrt(self.grad_pk[:,0]**2 + self.grad_pk[:,1]**2)
        # grad_p = self.calculate_gradient(P,G,self.pressure.copy(),smooth=False)
        # grad_p_mag = sqrt(grad_p[:,0]**2 + grad_p[:,1]**2)
        # self.grad_pk[:,0] = -grad_pk_mag*grad_p[:,0]/grad_p_mag
        # self.grad_pk[:,1] = -grad_pk_mag*grad_p[:,1]/grad_p_mag
