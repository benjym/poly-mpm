.. _tutorial:


***************
Tutorial
***************


Looking at the input file
=============================

Load up the file :file:`inputs/ball.py` in your favourite text editor. I like `atom <https://atom.io/>`_.

The Params class
=============================

This contains all of the parameters for controlling a single simulation. Some global information is stored here, such as timestepping and gravity values. Other information is stored in subclasses. A convenience method is also provided for updating gravity, or any other parameters, as many simulations require a slow ramp up in gravity, rather than applying a step function::

    class Params():
        def __init__(self,mode):
            self.dt = 1e-4 # timestep (s)
            self.savetime = 0.01
            self.t_f = 5.#100*self.dt # final time (s)
            self.nt = int(self.t_f/self.dt) # number of timesteps
            self.max_g = -10. # gravity (ms^-2)
            self.theta = 0.*pi/180. # slope angle (degrees)
            self.G = Grid_Params()
            self.B = Boundary_Params()
            self.O = Output_Params()
            self.S = Solid_Params()

        def update_forces(self):
            t_c = .5
            self.g=self.max_g

The Grid_Params class
=============================

This describes the :math:`x` and :math:`y` extents of the grid, as well as their resolution. This is pretty inflexible at the moment, and only allows for regular cartesian grids on a rectangular domain::

    class Grid_Params():
        def __init__(self):
            self.x_m = -1.0 # (m)
            self.x_M = 1.0 # (m)
            self.y_m = 0.0 # (m)
            self.y_M = 2.0 # (m)
            self.nx = 21 # number of grid edges in x direction
            self.ny = 21 # number of grid edges in y direction

The Boundary_Params class
=============================

Here all of the boundaries of the grid are defined. For this example, only the bottom of the grid has a boundary. It is also not rough (``self.roughness`` has not been set to  ``True``)::

    class Boundary_Params():
        def __init__(self):
            self.has_bottom = True

The Solid_Params class
=============================

This defines the properties of each material point. The lists ``self.X`` and ``self.Y`` define initial positions of each material point, up to ``self.n`` points in total. The constitutive model to use is defined by ``self.law`` --- this needs to be the same as the class in :file:`constit.py`. In this example, material points are placed onto a set of ``nr`` concentric circles with centre at ``c``, largest radius ``r`` and ``nphi`` material points around the largest circle.::

    class Solid_Params():
        def __init__(self):
            self.X = []
            self.Y = []
            self.n = 0

            self.rho = 1000. # density (kg/m^3)

    #         self.law = 'elastic'
            self.law = 'von_mises'
    #         self.law = 'dp'

            self.E = 1.e6 # elastic modulus (Pa)
            self.nu = 0.3 # poisson's ratio
            self.K = self.E/(3.*(1.-2.*self.nu)) # bulk modulus (Pa)
            self.G = self.E/(2.*(1.+self.nu)) # shear modulus (Pa)

            # for von_mises
            self.k = self.E/100.

            # for dp
            self.s = 5.
            self.beta = 10.
            self.mu = 0.5

            nr = 20 # particles in radial direction
            nphi = 50 # particles around circumference
            r = 0.3 # radius
            c = [0.,1.] # centre

            for i in linspace(0,r,nr):
                dnphi = around(nphi*i/r) # number in this ring
                for j in linspace(0,(1.-1./(dnphi))*2*pi,dnphi):
                    self.X.append(c[0]+i*sin(j))
                    self.Y.append(c[1]+i*cos(j))
                    self.n += 1
            self.A = pi*0.2**2/self.n # area (m^2)

The Output_Params class
=============================

This class defines what should be output from the simulation. Here, the continuum values and material point values are both show, using custom sizes for the respective figures.::

    class Output_Params():
        def __init__(self):
            self.continuum_fig_size = [10,6]
            self.mp_fig_size = [10,10]
            self.plot_continuum = True
            self.plot_material_points = True
