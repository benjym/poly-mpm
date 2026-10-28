"""Structure-of-arrays particle state and multiple-material coordination."""
from copy import copy
import numpy as np
from constit import EYE, get_law
from transfers import scatter, divide, gather, nodal_velocity


def field(value, shape, name):
    try:
        result = np.broadcast_to(np.asarray(value, dtype=float), shape).copy()
    except ValueError as exc:
        raise ValueError(f'{name} must broadcast to {shape}') from exc
    if not np.isfinite(result).all():
        raise ValueError(f'{name} must be finite')
    return result


class Particles:
    """One material's writable float64 particle arrays."""
    scalar_fields = ('rho', 'initial_rho', 'V', 'm', 'pressure', 'pk', 'sigma_kk',
                     'gammadot', 'yieldfunction', 'I', 'mu', 'eta', 'eta_limited', 'dp', 'de_kk')
    tensor_fields = ('strain', 'stress', 'dstrain', 'dstress', 'dev_stress', 'dev_dstress', 'de_ij')
    persistent_fields = ('x', 'v', 'x0', 'phi') + scalar_fields + tensor_fields

    def __init__(self, P, G, material=None, material_id=0):
        self.material = material if material is not None else P.S
        S = self.material
        self.material_id = material_id
        self.law = get_law(S.law)
        count = int(S.n)
        if count < 0 or count != S.n:
            raise ValueError('particle count must be a non-negative integer')
        for name in ('X', 'Y'):
            if np.shape(getattr(S, name)) != (count,):
                raise ValueError(f'{name} must have shape ({count},)')
        initial_v = getattr(S, 'initial_v', [0, 0, 0])
        if np.shape(initial_v) not in ((3,), (count, 3)):
            raise ValueError(f'initial_v must have shape (3,) or ({count}, 3)')
        composition = getattr(S, 'PHI', getattr(S, 'phi', np.ones(P.G.ns)/P.G.ns))
        if np.shape(composition) not in ((P.G.ns,), (count, P.G.ns)):
            raise ValueError(f'phi must have shape ({P.G.ns},) or ({count}, {P.G.ns})')
        self.x = np.zeros((count, 3))
        self.x[:, 0] = field(S.X, (count,), 'X')
        self.x[:, 1] = field(S.Y, (count,), 'Y')
        self.x0 = self.x.copy()
        self.v = field(initial_v, (count, 3), 'initial_v')
        self.phi = field(composition, (count, P.G.ns), 'phi')
        if np.any(self.phi < 0) or not np.allclose(self.phi.sum(axis=1), 1):
            raise ValueError('initial grain fractions must be nonnegative and sum to one')
        for name in self.scalar_fields:
            setattr(self, name, np.zeros(count))
        self.rho = field(S.rho, (count,), 'rho')
        self.initial_rho = self.rho.copy()
        self.V = field(S.A, (count,), 'A') * P.G.thickness
        self.m = self.rho * self.V
        if np.any(self.rho <= 0) or np.any(self.V <= 0):
            raise ValueError('initial density and volume must be positive')
        for name in self.tensor_fields:
            setattr(self, name, np.zeros((count, 3, 3)))
        self._scratch()
        self._initialise_stress(P)
        self._initialise_flow(P)
        self.check_finite(P)

    def __len__(self):
        return len(self.x)

    def _scratch(self):
        count = len(self)
        self.N = np.zeros((count, 4))
        self.G = np.zeros((count, 4, 3))
        self.n_star = np.zeros(count, dtype=np.intp)
        self.n = np.zeros((count, 4), dtype=np.intp)
        self.L = np.zeros((count, 3, 3))
        self.material.n = count

    @property
    def L_T(self):
        return self.L.swapaxes(1, 2)

    def _initialise_stress(self, P):
        S = self.material
        if P.pressure == 'lithostatic':
            depth = P.G.y_M - getattr(P.G, 'top_gap', 0) - self.x[:, 1]
            normal = np.cos(P.theta) * self.rho * P.max_g * depth
            shear = np.sin(P.theta) * self.rho * P.max_g * depth
            self.stress = normal[:, None, None] * EYE
            self.stress[:, 0, 1] = self.stress[:, 1, 0] = shear
        elif P.pressure == 'compression':
            self.stress[:] = -P.initial_pressure * EYE
        elif P.pressure in ('non-zero', 'temperature'):
            self.strain[:] = -1e-15 if P.pressure == 'non-zero' else -1e-5 * EYE
            tr = np.trace(self.strain, axis1=1, axis2=2)
            self.stress = S.K * tr[:, None, None] * EYE + 2 * S.G * (self.strain - tr[:, None, None] * EYE / 3)
        elif P.pressure is not None:
            raise ValueError(f'Unknown pressure initialisation {P.pressure!r}')
        if P.pressure in ('lithostatic', 'compression') and S.law != 'rigid':
            tr = np.trace(self.stress, axis1=1, axis2=2)
            self.strain = divide(tr, 9 * S.K)[:, None, None] * EYE
            if getattr(S, 'G', 0):
                self.strain += (self.stress - tr[:, None, None] * EYE / 3) / (2 * S.G)
        self.sigma_kk = np.trace(self.stress, axis1=1, axis2=2)
        self.pressure = self.sigma_kk / 3
        self.dev_stress = self.stress - self.pressure[:, None, None] * EYE

    def _initialise_flow(self, P):
        S = self.material
        y = self.x[:, 1]
        if P.initial_flow == 'steady':
            if S.law == 'viscous':
                self.v[:, 0] = self.rho * P.max_g * np.sin(P.theta) * (P.G.y_M * y - y*y/2) / S.mu_s
            elif S.law == 'pouliquen':
                size = P.G.s_bar_0
                a = (np.sqrt(abs(P.max_g) * size) * 2/3 * S.I_0
                     * (np.tan(abs(P.theta)) - S.mu_0) / (S.mu_1 - np.tan(abs(P.theta)))
                     * np.sqrt(S.packing * np.cos(abs(P.theta))))
                self.v[:, 0] = a * (P.G.y_M**1.5 - (P.G.y_M-y)**1.5) / size**1.5
                self.gammadot = 1.5 * a * np.sqrt(P.G.y_M-y) / size**1.5
                self.pk = P.l * size * np.sqrt(np.abs(self.pressure*self.rho)) * np.abs(self.gammadot)
            else:
                raise ValueError(f'steady initial flow is unsupported for {S.law}')
        elif P.initial_flow == 'steady_poiseulle':
            height = P.G.y_M - P.G.y_m
            ys = y / height
            ys = np.minimum(ys, 1-ys)
            k = abs(float(np.asarray(S.rho).flat[0]) * P.max_g)
            beta = np.sqrt(abs(k * height / P.initial_pressure))
            if beta >= np.sqrt(2*S.mu_1):
                raise ValueError('Poiseuille initialisation is outside its steady regime')
            yc = .5 - S.mu_0 / beta**2
            yp = np.minimum(ys, yc)
            u = (height/P.G.s_bar_0*S.I_0/beta**3
                 * (beta**2*yp + (S.mu_0-S.mu_1)
                    * np.log(-(2*beta**2*yp + 2*S.mu_1-beta**2)/(beta**2-2*S.mu_1))))
            self.v[:, 0] = -u * np.sqrt(k*height/S.rho_s)
        elif P.initial_flow not in (False, None):
            raise ValueError(f'Unknown initial flow {P.initial_flow!r}')
        if S.law == 'rigid':
            self.v[:] = self.prescribed_velocity(P)

    def prescribed_velocity(self, P):
        velocity = getattr(self.material, 'initial_v', getattr(P, 'intruder_v', [0, 0, 0]))
        # Per-particle prescriptions travel with their rows through removals.
        if np.ndim(velocity) == 2:
            return self.v.copy()
        return field(velocity, (len(self), 3), 'rigid velocity')

    def check_finite(self, P):
        for name in self.persistent_fields:
            value = getattr(self, name)
            bad = np.flatnonzero(~np.isfinite(value).reshape(len(self), -1).all(axis=1)) if len(self) else []
            if len(bad):
                raise FloatingPointError(f'{self.material.law}, material {self.material_id}, '
                                         f'timestep {P.tstep}, {name}, particles {bad[:10].tolist()}')

    def cyclic_lr(self, P, G=None):
        self.x[:, 0] = P.G.x_m + (self.x[:, 0]-P.G.x_m) % (P.G.x_M-P.G.x_m)

    def get_reference_nodes(self, P, G):
        if P.B.cyclic_lr:
            self.cyclic_lr(P)
        if not np.isfinite(self.x).all():
            self.check_finite(P)
        lo = np.array([P.G.x_m, P.G.y_m])
        hi = np.array([P.G.x_M, P.G.y_M])
        if np.any((self.x[:, :2] < lo) | (self.x[:, :2] > hi)):
            raise ValueError(f'material {self.material_id}: particle outside domain before cell lookup')
        cell = np.floor((self.x[:, :2]-lo) / [G.dx, G.dy]).astype(np.intp)
        cell = np.minimum(cell, [P.G.nx-2, P.G.ny-2])
        self.n_star = cell[:, 0] + cell[:, 1]*P.G.nx
        self.n = self.n_star[:, None] + [0, 1, P.G.nx+1, P.G.nx]

    def get_basis_functions(self, P, G):
        x = (self.x[:, 0]-G.X[self.n_star]) / G.dx
        y = (self.x[:, 1]-G.Y[self.n_star]) / G.dy
        self.N[:] = np.stack(((1-x)*(1-y), x*(1-y), x*y, (1-x)*y), axis=1)
        self.G[:, :, 0] = np.stack((-(1-y), 1-y, y, -y), axis=1) / G.dx
        self.G[:, :, 1] = np.stack((-(1-x), -x, x, 1-x), axis=1) / G.dy

    def get_nodal_mass_momentum(self, P, G):
        weight = self.N * self.m[:, None]
        scatter(G.m, self.n, weight)
        scatter(G.q, self.n, weight[:, :, None] * self.v[:, None, :])
        scatter(G.phim, self.n, weight[:, :, None] * self.phi[:, None, :])
        scatter(G.pkm, self.n, weight * self.pk[:, None])

    def update_stress_strain(self, P, G):
        if self.material.law == 'rigid':
            return
        velocity = nodal_velocity(G, P)[self.n]
        self.L = np.einsum('nri,nrj->nij', velocity, self.G)
        self.dstrain = .5 * (self.L+self.L_T) * P.dt
        self.strain += self.dstrain
        self.law(self, self.material, P)
        self.stress += self.dstress - (self.L @ self.stress + self.stress @ self.L_T)*P.dt
        self.check_finite(P)

    def get_nodal_forces(self, P, G):
        force = np.einsum('nij,nrj->nri', self.stress, self.G) * self.V[:, None, None]
        scatter(G.fi, self.n, force)
        body = np.array([P.g*np.sin(P.theta), P.g*np.cos(P.theta), 0])
        scatter(G.fe, self.n, self.N[:, :, None] * self.m[:, None, None] * body)

    def deposit_diagnostics(self, P, G):
        weight = self.N * self.m[:, None]
        diagnostics = dict(pressure=self.pressure, gammadot=self.gammadot,
                           sigmav=self.stress[:, 1, 1], sigmah=self.stress[:, 0, 0],
                           yieldfunction=self.yieldfunction, mu=self.mu, I=self.I,
                           eta=self.eta_limited if self.material.law == 'pouliquen' else self.eta,
                           dev_stress=np.linalg.norm(self.dev_stress, axis=(1, 2))/(np.sqrt(2) if self.material.law == 'pouliquen' else 1),
                           dev_stress_dot=np.linalg.norm(self.dev_dstress, axis=(1, 2))/P.dt)
        for name, value in diagnostics.items():
            scatter(getattr(G, name), self.n, weight*value[:, None])

    def move_material_points(self, P, G):
        if self.material.law == 'rigid':
            self.v[:] = self.prescribed_velocity(P)
            self.x += P.dt*self.v
        else:
            self.x += P.dt*gather(nodal_velocity(G, P), self.n, self.N, G.valid_nodes)
            self.v += P.dt*gather(nodal_velocity(G, P, 'q_dot'), self.n, self.N, G.valid_nodes)
            self.rho = gather(divide(G.m, G.V), self.n, self.N, G.valid_nodes)
        self.check_finite(P)

    def move_grainsize_on_grid(self, P, G, increment=None):
        if increment is None:
            increment = gather(G.dphi, self.n, self.N, G.valid_nodes)
        self.phi += increment
        self.pk += gather(G.dpk, self.n, self.N, G.valid_nodes)
        if P.normalise_phi:
            self.phi = np.clip(self.phi, 0, 1)
            total = self.phi.sum(axis=1, keepdims=True)
            if np.any(total == 0):
                raise FloatingPointError('Cannot normalise zero grain-fraction sum')
            self.phi /= total
        self.check_finite(P)

    def keep(self, mask):
        for name in self.persistent_fields:
            setattr(self, name, getattr(self, name)[mask])
        # Preserve per-particle injection specifications separately in coordinator.
        self._scratch()

    def append(self, batch):
        for name in self.persistent_fields:
            setattr(self, name, np.concatenate((getattr(self, name), getattr(batch, name))))
        self._scratch()


class ParticleSystem:
    """Coordinate material batches without individual particle objects."""
    def __init__(self, P, G):
        self.materials = [Particles(P, G, S, i) for i, S in enumerate(P.S)]
        self.templates = [copy(S) for S in P.S]

    def get_reference_nodes(self, P, G):
        for batch in self.materials:
            batch.get_reference_nodes(P, G)

    get_reference_node = get_reference_nodes

    def get_basis_functions(self, P, G):
        for batch in self.materials:
            batch.get_basis_functions(P, G)

    def get_nodal_mass_momentum(self, P, G):
        if P.segregate:
            G.segregation_m = np.zeros_like(G.m)
            G.segregation_phim = np.zeros_like(G.phim)
        for batch in self.materials:
            batch.get_nodal_mass_momentum(P, G)
            if P.segregate and batch.material.law != 'rigid':
                weight = batch.N*batch.m[:, None]
                scatter(G.segregation_m, batch.n, weight)
                scatter(G.segregation_phim, batch.n, weight[:, :, None]*batch.phi[:, None, :])
        if P.B.cyclic_lr:
            fields = ['m', 'q', 'phim', 'pkm']
            if P.segregate:
                fields += ['segregation_m', 'segregation_phim']
            G.make_cyclic(P, G, fields)
        if P.segregate:
            G.segregation_phi = divide(G.segregation_phim, G.segregation_m[:, None])
        G.get_valid_nodes(P)
        G.phi = divide(G.phim, G.m[:, None])
        G.s_bar = G.phi @ P.G.s
        G.pk = divide(G.pkm, G.m)
        G.q[G.boundary_v, 0] = 0
        G.q[G.boundary_h, 1] = 0
        if P.B.roughness:
            G.q[G.boundary_h, 0] = 0
            G.q[G.boundary_v, 1] = 0

    def update_stress_strain(self, P, G):
        for batch in self.materials:
            batch.update_stress_strain(P, G)
            batch.deposit_diagnostics(P, G)
        if P.B.cyclic_lr:
            G.make_cyclic(P, G, G.diagnostic_fields)

    def get_nodal_forces(self, P, G):
        for batch in self.materials:
            batch.get_nodal_forces(P, G)
        if P.B.cyclic_lr:
            G.make_cyclic(P, G, ['fi', 'fe'])

    def move_material_points(self, P, G):
        for batch in self.materials:
            batch.move_material_points(P, G)

    def move_grainsize_on_grid(self, P, G):
        for i, batch in enumerate(self.materials):
            if batch.material.law != 'rigid':
                increments = G.segregation_particle_dphi
                batch.move_grainsize_on_grid(P, G, None if increments is None else increments[i])

    def cyclic_lr(self, P, G=None):
        for batch in self.materials:
            batch.cyclic_lr(P)

    def process_boundaries(self, P, G):
        for batch, template in zip(self.materials, self.templates):
            keep = np.ones(len(batch), dtype=bool)
            if P.B.outlet_left:
                keep &= batch.x[:, 0] >= P.G.x_m + G.dx
            if P.B.outlet_bottom:
                keep &= batch.x[:, 1] >= P.G.y_m
            if not keep.all():
                batch.keep(keep)
            for inlet in ('right', 'top'):
                if not getattr(P.B, 'inlet_'+inlet) or not hasattr(template, 'inlet_rate'):
                    continue
                rate = template.inlet_rate
                if rate <= 0:
                    raise ValueError('inlet_rate must be positive')
                if not P.t or P.t % rate >= P.dt:
                    continue
                S = copy(template)
                if inlet == 'right':
                    S.Y = np.asarray(template.yp)
                    S.X = np.full(len(S.Y), P.G.x_M-template.gap[0])
                else:
                    S.X = np.array([(P.G.x_M-P.G.x_m)/2])
                    S.Y = np.array([P.G.y_M-template.gap[1]])
                S.n = len(S.X)
                for name in ('rho', 'A', 'initial_v', 'PHI', 'phi'):
                    value = getattr(S, name, None)
                    uniform_shape = {'initial_v': (3,), 'phi': (P.G.ns,), 'PHI': (P.G.ns,)}.get(name, ())
                    if value is not None and np.asarray(value).shape != uniform_shape:
                        raise ValueError(f'inlet {name} must be uniform per material')
                batch.append(Particles(P, G, S, batch.material_id))
            if P.B.cyclic_lr:
                batch.cyclic_lr(P)

    def recover_position(self, P, G):
        coordinates = np.column_stack((G.X, G.Y, np.zeros_like(G.X)))
        for batch in self.materials:
            np.testing.assert_allclose(np.sum(batch.N[:, :, None]*coordinates[batch.n], axis=1), batch.x, atol=1e-10)

    def measure_elastic_energy(self, P, G=None):
        kinetic = potential = strain = 0.
        direction = np.array([np.sin(P.theta), np.cos(P.theta), 0])
        for batch in self.materials:
            kinetic += .5*np.sum(batch.m*np.sum(batch.v**2, axis=1))
            potential -= P.g*np.sum(batch.m*(batch.x@direction))
            strain += .5*np.sum(batch.V*np.sum(batch.stress*batch.strain, axis=(1, 2)))
        return [P.t, kinetic, potential, strain]
