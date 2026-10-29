"""Regular two-dimensional grid; node numbering is x-fastest."""

import numpy as np
from transfers import divide, nodal_velocity
from integrators import segregation_update


class Grid:
    diagnostic_fields = (
        "pressure",
        "gammadot",
        "yieldfunction",
        "sigmah",
        "sigmav",
        "dev_stress",
        "dev_stress_dot",
        "mu_s",
        "mu",
        "I",
        "eta",
    )

    def __init__(self, P):
        if P.G.nx < 2 or P.G.ny < 2:
            raise ValueError("grid needs at least two nodes in each direction")
        self.x = np.linspace(P.G.x_m, P.G.x_M, P.G.nx)
        self.y = np.linspace(P.G.y_m, P.G.y_M, P.G.ny)
        self.dx, self.dy = self.x[1] - self.x[0], self.y[1] - self.y[0]
        if self.dx <= 0 or self.dy <= 0 or P.G.thickness <= 0:
            raise ValueError("grid extent and thickness must be positive")
        P.G.dx, P.G.dy = self.dx, self.dy
        self.X, self.Y = np.tile(self.x, P.G.ny), np.repeat(self.y, P.G.nx)
        self.boundary(P)
        self.volume(P)
        self.x_plot = np.r_[self.x[0] - self.dx / 2, (self.x[1:] + self.x[:-1]) / 2, self.x[-1] + self.dx / 2]
        self.y_plot = np.r_[self.y[0] - self.dy / 2, (self.y[1:] + self.y[:-1]) / 2, self.y[-1] + self.dy / 2]
        self.dx_plot = np.r_[self.dx / 2, np.full(P.G.nx - 2, self.dx), self.dx / 2]
        self.dy_plot = np.r_[self.dy / 2, np.full(P.G.ny - 2, self.dy), self.dy / 2]
        self.wipe(P)

    def boundary(self, P):
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
        self.boundary_v = np.zeros((P.G.ny * P.G.nx), dtype=bool)
        self.boundary_h = np.zeros((P.G.ny * P.G.nx), dtype=bool)
        if P.B.has_bottom:
            self.boundary_h[: P.G.nx] = 1  # bottom
        if P.B.no_slip_bottom:
            self.boundary_v[: P.G.nx] = 1  # bottom
        if P.B.has_top:
            self.boundary_h[(P.G.ny - 1) * P.G.nx :] = 1  # top
        if P.B.has_left:
            self.boundary_v[:: P.G.nx] = 1  # left
        if P.B.has_right:
            self.boundary_v[P.G.nx - 1 :: P.G.nx] = 1  # right
        if P.B.wall:
            self.boundary_v[(P.G.nx - 1) // 2 : (P.G.ny * P.G.nx - 1) // 2 : P.G.nx] = 1  # wall at centre
        if P.B.two_walls:
            self.boundary_v[(P.G.nx - 1) // 4 : (P.G.ny * P.G.nx - 1) // 2 : P.G.nx] = 1  # low wall at left
            self.boundary_v[(P.G.ny // 2) * P.G.nx + (3 * (P.G.nx - 1)) // 4 :: P.G.nx] = (
                1  # upper wall at right
            )
        if P.B.silo_left:
            self.boundary_v[1 : (P.G.ny // 2 - 3) * P.G.nx : P.G.nx] = 1  # wall at outlet
            self.boundary_v[(P.G.ny // 2 + 3) * P.G.nx + 1 : -1 : P.G.nx] = 1  # wall at outlet
        if P.B.silo_bottom:
            self.boundary_h[: P.G.nx // 2 - 2] = 1  # wall at outlet
            self.boundary_h[P.G.nx // 2 + 3 : P.G.nx] = 1  # wall at outlet
        if P.B.box:
            l = 4  # min grid points in from left
            r = 8  # max grid points across from left
            b = 4  # min grid points up from bottom
            t = 8  # max grid points up from bottom
            self.boundary_v[b * P.G.nx + l : t * P.G.nx + l + 1 : P.G.nx] = 1  # left wall
            self.boundary_v[b * P.G.nx + r : t * P.G.nx + r + 1 : P.G.nx] = 1  # right wall
            self.boundary_h[b * P.G.nx + l : b * P.G.nx + r + 1] = 1  # bottom wall
            self.boundary_h[t * P.G.nx + l : t * P.G.nx + r + 1] = 1  # top wall
        if P.B.conveyor:
            self.boundary_conveyor = np.zeros((P.G.ny * P.G.nx), dtype=bool)
            self.boundary_conveyor_triple = np.zeros((P.G.ny * P.G.nx, 3), dtype=bool)
            self.boundary_conveyor[2 : P.G.nx - 2] = 1
            self.boundary_conveyor_triple[:, 0] = self.boundary_conveyor

        self.boundary_tot = self.boundary_v + self.boundary_h

    def volume(self, P):  # get volume of each cell
        self.weighting = np.ones((P.G.ny * P.G.nx))
        self.weighting[: P.G.nx] /= 2  # bottom
        self.weighting[(P.G.ny - 1) * P.G.nx :] /= 2  # top
        if not P.B.cyclic_lr:
            self.weighting[:: P.G.nx] /= 2  # left
            self.weighting[P.G.nx - 1 :: P.G.nx] /= 2  # right
        self.V = self.weighting * self.dx * self.dy * P.G.thickness

    def wipe(self, P):
        count = P.G.nx * P.G.ny
        for name in ("fi", "fe", "q", "q_dot", "ext_f", "grad_gammadot", "grad_pk"):
            setattr(self, name, np.zeros((count, 3)))
        for name in ("m", "s_bar", "pk", "pkm", "dpk") + self.diagnostic_fields:
            setattr(self, name, np.zeros(count))
        for name in ("phi", "phim", "dphi"):
            setattr(self, name, np.zeros((count, P.G.ns)))
        self.valid_nodes = np.zeros(count, dtype=np.bool_)
        self.shear_rate = np.zeros(count)
        self.q_reference = np.zeros_like(self.q)
        self.constraint_impulse = np.zeros(3)
        self.constraint_angular_impulse = np.zeros(3)
        self.external_angular_impulse = np.zeros(3)
        self.force_scaling_impulse = np.zeros(3)
        self.force_scaling_angular_impulse = np.zeros(3)
        self.external_impulse = np.zeros(3)
        self.external_work = self.constraint_work = self.force_scaling_work = 0.0
        self.segregation_m = None
        self.segregation_phi = None
        self.segregation_particle_dphi = None

    def get_valid_nodes(self, P):
        self.valid_nodes = self.m > P.M_tol

    def make_cyclic(self, P, G, params):
        for name in params:
            value = getattr(G, name)
            combined = value[:: P.G.nx] + value[P.G.nx - 1 :: P.G.nx]
            value[:: P.G.nx] = combined
            value[P.G.nx - 1 :: P.G.nx] = combined

    def physical(self, name):
        """Convert a mass-weighted diagnostic to a physical nodal field once."""
        return divide(getattr(self, name), self.m)

    def BCs(self, P):
        # Preserve the legacy traction scaling and reference-material density.
        density = np.asarray(P.S[0].rho)
        if density.ndim and not density.size:
            if P.B.vertical_force or P.B.horizontal_force:
                raise ValueError("traction boundaries require a reference-material density")
            density = np.asarray(1.0)
        if density.ndim:
            if not np.all(density == density.flat[0]):
                if P.B.vertical_force or P.B.horizontal_force:
                    raise ValueError("traction boundaries require uniform reference-material density")
            density = density.flat[0] if density.size else 1.0
        if P.B.vertical_force:
            self.ext_f[: P.G.nx, 1] += P.q_v
            self.ext_f[-P.G.nx :, 1] -= P.q_v
            self.fe[:, 1] += 2 * self.ext_f[:, 1] * self.m / density / self.dy
        if P.B.horizontal_force:
            self.ext_f[:: P.G.nx, 0] += P.q_h
            self.ext_f[P.G.nx - 1 :: P.G.nx, 0] -= P.q_h
            self.fe[:, 0] += 2 * self.ext_f[:, 0] * self.m / density / self.dx
        if P.mode == "anisotropy" and P.t == 0:
            self.fe[P.G.nx * P.G.ny // 2, 2] = 1

    def wall_contact(self, P, velocity):
        """Wall nodes whose velocity does not move them away from their wall.

        With P.B.separable the domain edges are one-sided: a node moving away
        from its wall is left free. With wall_mu in the granular solver the
        domain edges are frictional contact rows in the solve and never fixed
        here. Interior walls stay two-sided.
        """
        h, v = self.boundary_h.copy(), self.boundary_v.copy()
        row, column = np.divmod(np.arange(len(self.m)), P.G.nx)
        if P.solver == "granular" and P.B.wall_mu:
            h &= (row != 0) & (row != P.G.ny - 1)
            v &= (column != 0) & (column != P.G.nx - 1)
        elif P.B.separable:
            h &= ~(((row == 0) & (velocity[:, 1] > 0)) | ((row == P.G.ny - 1) & (velocity[:, 1] < 0)))
            v &= ~(((column == 0) & (velocity[:, 0] > 0)) | ((column == P.G.nx - 1) & (velocity[:, 0] < 0)))
        return h, v

    def update_momentum(self, P):
        initial_velocity = divide(self.q, self.m[:, None])
        initial_momentum = self.q.copy()
        physical_force = self.fe - self.fi
        self.q_dot = physical_force.copy()
        if P.damping:
            self.q_dot *= 0.7
        self.contact_h, self.contact_v = self.wall_contact(P, divide(self.q + self.q_dot * P.dt, self.m[:, None]))
        if P.B.roughness:
            if P.B.wall_mu and P.solver != "granular":
                # Friction acts only where the normal force points into a wall.
                outward_h = np.zeros_like(self.m, dtype=np.bool_)
                outward_v = np.zeros_like(self.m, dtype=np.bool_)
                outward_h[: P.G.nx] = self.q_dot[: P.G.nx, 1] < 0
                outward_h[-P.G.nx :] = self.q_dot[-P.G.nx :, 1] > 0
                outward_v[:: P.G.nx] = self.q_dot[:: P.G.nx, 0] < 0
                outward_v[P.G.nx - 1 :: P.G.nx] = self.q_dot[P.G.nx - 1 :: P.G.nx, 0] > 0
                fh = (
                    self.boundary_h
                    * outward_h
                    * np.minimum(np.abs(self.q_dot[:, 0]), P.B.wall_mu * np.abs(self.q_dot[:, 1]))
                )
                fv = (
                    self.boundary_v
                    * outward_v
                    * np.minimum(np.abs(self.q_dot[:, 1]), P.B.wall_mu * np.abs(self.q_dot[:, 0]))
                )
                self.q_dot[:, 0] -= np.sign(self.q_dot[:, 0]) * fh
                self.q_dot[:, 1] -= np.sign(self.q_dot[:, 1]) * fv
            else:
                self.q_dot[self.contact_h, 0] = 0
                self.q_dot[self.contact_v, 1] = 0
        self.q_dot[self.contact_v, 0] = 0
        self.q_dot[self.contact_h, 1] = 0
        self.q += self.q_dot * P.dt
        if P.B.conveyor:
            self.q[self.boundary_conveyor, 0] = P.v_0 * self.m[self.boundary_conveyor]

        # All sums use unique physical nodes on periodic grids.
        unique = np.ones(len(self.m), dtype=bool)
        if P.B.cyclic_lr:
            unique[P.G.nx - 1 :: P.G.nx] = False
        midpoint = 0.5 * (initial_velocity + divide(self.q, self.m[:, None]))
        scaling = (-0.3 * physical_force * P.dt) if P.damping else np.zeros_like(self.q)
        reaction = self.q - initial_momentum - P.dt * physical_force - scaling
        coordinates = np.column_stack((self.X, self.Y, np.zeros_like(self.X)))
        self.external_angular_impulse = P.dt * np.cross(coordinates[unique], self.fe[unique]).sum(axis=0)
        self.constraint_angular_impulse += np.cross(coordinates[unique], reaction[unique]).sum(axis=0)
        self.force_scaling_impulse = scaling[unique].sum(axis=0)
        self.force_scaling_angular_impulse = np.cross(coordinates[unique], scaling[unique]).sum(axis=0)
        self.external_impulse = P.dt * self.fe[unique].sum(axis=0)
        self.constraint_impulse += reaction[unique].sum(axis=0)
        self.external_work = float(P.dt * np.sum(self.fe[unique] * midpoint[unique]))
        self.constraint_work += float(np.sum(reaction[unique] * midpoint[unique]))
        self.force_scaling_work = float(np.sum(scaling[unique] * midpoint[unique]))

    def calculate_gradient(self, P, G, Z, smooth=False):
        if not smooth:
            values = Z.reshape(P.G.ny, P.G.nx)
            valid = ~(G.m <= P.M_tol).reshape(P.G.ny, P.G.nx)
            dy, dx = np.gradient(values, G.dy, G.dx)
            # A central difference needs the two neighbouring values, even
            # when the centre itself is empty. Edge differences need both
            # endpoints. This matches the masked-array free-surface stencil.
            for derivative, axis in ((dy, 0), (dx, 1)):
                out = np.moveaxis(derivative, axis, 0)
                active = np.moveaxis(valid, axis, 0)
                out[1:-1] = np.where(active[:-2] & active[2:], out[1:-1], 0)
                out[0] = np.where(active[0] & active[1], out[0], 0)
                out[-1] = np.where(active[-2] & active[-1], out[-1], 0)
                # Masked division also masks nonfinite derivative results.
                out[~np.isfinite(out)] = 0
            return np.column_stack((dx.ravel(), dy.ravel(), np.zeros_like(G.m)))
        # Preserve masked free-surface differences; explicitly fill masked output
        # rather than exposing the arbitrary backing values of a masked array.
        values = np.ma.masked_where(G.m <= P.M_tol, Z).reshape(P.G.ny, P.G.nx)
        dy, dx = np.gradient(values, G.dy, G.dx)
        if smooth:
            from astropy.convolution import convolve, Gaussian2DKernel

            kernel = Gaussian2DKernel(x_stddev=1, y_stddev=1)
            dy = convolve(dy.filled(np.nan), kernel, boundary="extend")
            dx = convolve(dx.filled(np.nan), kernel, boundary="extend")
        return np.column_stack(
            (
                np.nan_to_num(np.ma.filled(dx, 0)).ravel(),
                np.nan_to_num(np.ma.filled(dy, 0)).ravel(),
                np.zeros_like(G.m),
            )
        )

    def calculate_gammadot(self, P, G, smooth=False):
        velocity = nodal_velocity(self, P)
        gu = self.calculate_gradient(P, G, velocity[:, 0], smooth)
        gv = self.calculate_gradient(P, G, velocity[:, 1], smooth)
        # Physical field, distinct from the deposited mass-weighted diagnostic.
        self.shear_rate = gu[:, 1] + gv[:, 0]

    def calculate_grad_gammadot(self, P, G):
        self.grad_gammadot = self.calculate_gradient(P, G, np.abs(self.shear_rate), P.smooth_grad2)

    def calculate_phi_increment(self, P, L=None):
        self.dphi, self.segregation_particle_dphi = segregation_update(
            P, self, None if L is None else L.materials
        )
        if not np.isfinite(self.dphi).all():
            raise FloatingPointError(f"non-finite segregation increment at timestep {P.tstep}")

    def update_pk(self, P, G):
        self.pk = -self.physical("pressure") * self.physical("I")
        self.grad_pk = self.calculate_gradient(P, G, self.pk)
