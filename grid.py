from numpy import linspace,tile,repeat,zeros,array,sum,sign,ones,sqrt,abs
import matplotlib.pyplot as plt

class Grid():
    """
    GRID NUMBERING:
    
    (ny-1)*nx  ....     (ny*nx-1)
    ....       ....     ....
    2*nx -     ....    - (3*nx-1)
    nx   -     ....    - (2*nx-1)
    0    - 1 - 2 - ... - nx-1
    """
    def __init__(self,P):
        self.x = linspace(P.G.x_m,P.G.x_M,P.G.nx)
        self.y = linspace(P.G.y_m,P.G.y_M,P.G.ny)
        self.dx = self.x[1] - self.x[0] # grid spacing (m)
        self.dy = self.y[1] - self.y[0] # grid spacing (m)
        self.X = tile(self.x,P.G.ny)
        self.Y = repeat(self.y,P.G.nx)
        self.boundary(P)
        self.volume(P)
        
    def boundary(self,P):
        self.boundary_v = zeros((P.G.ny*P.G.nx))
        self.boundary_h = zeros((P.G.ny*P.G.nx))
        if P.B.has_bottom:
            self.boundary_h[:P.G.nx] = 1 # bottom
        if P.B.has_top:
            self.boundary_h[(P.G.ny-1)*P.G.nx:] = 1 # top
        if P.B.has_left:
            self.boundary_v[::P.G.nx] = 1 # left
        if P.B.has_right:
            self.boundary_v[P.G.nx-1::P.G.nx] = 1 # right
        if P.B.wall:
            self.boundary_v[10:(P.G.ny*P.G.nx/2.):P.G.nx] = 1. # wall
        self.boundary = sign(self.boundary_v + self.boundary_h)

    def volume(self,P): # get volume of each cell
        self.weighting = ones((P.G.ny*P.G.nx))
        self.weighting[:P.G.nx] /=2 # bottom
        self.weighting[(P.G.ny-1)*P.G.nx:] /= 2 # top
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
        self.yieldfunction = zeros((P.G.nx*P.G.ny)) # von Mises yield function
        self.pressure = zeros((P.G.nx*P.G.ny)) # normal stress
        self.sigmah = zeros((P.G.nx*P.G.ny)) # 
        self.sigmav = zeros((P.G.nx*P.G.ny)) # 
        self.dev_stress = zeros((P.G.nx*P.G.ny)) # 
        self.dev_stress_dot = zeros((P.G.nx*P.G.ny)) # 
        
        self.damping_force = zeros((P.G.nx*P.G.ny,3)) # local non-viscous damping
        try:
            self.gsd = zeros((P.G.nx*P.G.ny,P.GSD.ns))
            self.s_bar = zeros((P.G.nx*P.G.ny))
        except:
            pass
            
            
    def nearby_nodes(self,n_star,r,P):
        if (r == 0) or (r == 1):
            return n_star + r
        else:
            return n_star + P.G.nx + 3-r

    def N(self,x):
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

    def BCs(self,P): # BCS directly affecting self.fe
        if P.B.vertical_force:
            self.ext_f[:P.G.nx,1] += P.q*ones((P.G.nx)) # bottom
            self.ext_f[(P.G.ny-1)*P.G.nx:,1] += -P.q*ones((P.G.nx)) # top
            self.fe[:,1] += 2.*self.ext_f[:,1]*self.m/P.S.rho/self.dy

        if P.B.horizontal_force:
            self.ext_f[::P.G.nx,0] += P.q*ones((P.G.ny)) # left
            self.ext_f[P.G.nx-1::P.G.nx,0] += -P.q*ones((P.G.ny)) # right
            self.fe[:,0] += 2.*self.ext_f[:,0]*self.m/P.S.rho/self.dx
            
        if P.mode == 'anisotropy' and P.t == 0:
            self.fe[P.G.nx*P.G.ny/2,2] = 1.
            
    def update_momentum(self,P):
        if P.damping:
            self.damping_force = 0.8*abs(self.fe - self.fi)*sign(self.q_dot)
        self.q_dot = self.fe - self.fi - self.damping_force
        self.q_dot[:,0] = self.q_dot[:,0]*(1.-self.boundary_v) # 0 at boundary
        self.q_dot[:,1] = self.q_dot[:,1]*(1.-self.boundary_h) # 0 at boundary
        if P.B.roughness:
            self.q_dot[:,1] = self.q_dot[:,1]*(1.-self.boundary_v) # sidewalls
            self.q_dot[:,0] = self.q_dot[:,0]*(1.-self.boundary_h) # top/bottom
        self.q += self.q_dot*P.dt
