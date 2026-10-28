"""Independent scalar kernels derived from master's particle/mplist/constit code.

Fixes intentionally applied: L and transpose are distinct; repeated-node additions
remain scalar; spatial/grain axes cannot mix. No vector kernel is called here.
"""
from types import SimpleNamespace as NS
import numpy as np
from particles import Particles


def law(dstrain, stress, history, S, P, rho, phi, initial_rho):
    h = dict(history)
    eye = np.eye(3)
    tr = np.trace(dstrain)
    dev = dstrain-tr*eye/3
    if S.law == 'elastic':
        dstress = S.K*tr*eye+2*S.G*dev
        h['pressure'] = np.trace(stress)/3
        h['gammadot'] = np.linalg.norm(dev)/P.dt
    elif S.law == 'viscous':
        h['de_kk'], h['de_ij'] = tr/3, dev
        h['dp'] = S.K*tr/3
        h['pressure'] += h['dp']
        h['dev_stress'] = 2*S.mu_s*dev/P.dt
        dstress = h['pressure']*eye+h['dev_stress']-S.mu_v*tr/3*eye/P.dt-stress
        h['gammadot'] = np.linalg.norm(dev)/P.dt
    elif S.law == 'von_mises':
        work = abs(np.sum(h['dev_stress']*dev))
        norm = np.linalg.norm(h['dev_stress'])
        factor = work/(2*S.k**2)*(norm/(np.sqrt(2)*S.k))**(S.s-2) if norm else 0
        h['dev_dstress'] = 2*S.G*(dev-factor*h['dev_stress'])
        h['sigma_kk'] += 3*S.K*tr
        h['dev_stress'] = h['dev_stress']+h['dev_dstress']
        dstress = h['dev_stress']+h['sigma_kk']*eye/3-stress
        h['yieldfunction'] = norm/(np.sqrt(2)*S.k)-1
        h['pressure'] = np.trace(stress)/3
        h['gammadot'] = np.linalg.norm(dev)/P.dt
    elif S.law == 'pouliquen':
        h['de_kk'], h['de_ij'] = tr/3, dev
        h['gammadot'] = np.linalg.norm(2*dev/P.dt)
        pressure = abs(h['pressure'])
        h['I'] = h['gammadot']*np.dot(phi, P.G.s)*np.sqrt(S.rho_s/pressure) if pressure else 0.
        h['mu'] = S.mu_0+S.delta_mu*(h['I']/(S.I_0+h['I']) if h['I'] else (1 if h['gammadot'] and not pressure else 0))
        numerator = 2*np.sqrt(2)*h['mu']*pressure
        h['eta'] = numerator/h['gammadot'] if h['gammadot'] else (np.finfo(float).max if numerator else 0)
        h['eta_limited'] = min(h['eta'], S.eta_max)
        h['dev_stress'] = h['eta_limited']*dev/P.dt
        h['dp'] = S.K*tr/3
        h['pressure'] += h['dp']
        if h['pressure'] > 0 or rho < initial_rho*.8:
            h['dp'] -= h['pressure']
            h['pressure'] = 0.
        dstress = h['pressure']*eye+h['dev_stress']-stress
    else:
        dstress = np.zeros((3, 3))
    h['dstress'] = dstress
    return h


class ScalarParticles(Particles):
    def get_reference_nodes(self, P, G):
        if P.B.cyclic_lr:
            self.cyclic_lr(P)
        for i in range(len(self)):
            x, y = self.x[i, :2]
            if not (P.G.x_m <= x <= P.G.x_M and P.G.y_m <= y <= P.G.y_M):
                raise ValueError('outside domain')
            ix = min(int(np.floor((x-P.G.x_m)/G.dx)), P.G.nx-2)
            iy = min(int(np.floor((y-P.G.y_m)/G.dy)), P.G.ny-2)
            self.n_star[i] = ix+iy*P.G.nx
            self.n[i] = self.n_star[i]+np.array([0, 1, P.G.nx+1, P.G.nx])

    def get_basis_functions(self, P, G):
        for i in range(len(self)):
            x, y = self.x[i, :2]-[G.X[self.n_star[i]], G.Y[self.n_star[i]]]
            self.N[i] = [(G.dx-x)*(G.dy-y), x*(G.dy-y), x*y, (G.dx-x)*y]
            self.N[i] /= G.dx*G.dy
            self.G[i] = np.array([[-(G.dy-y), -(G.dx-x), 0],
                                  [G.dy-y, -x, 0], [y, x, 0], [-y, G.dx-x, 0]])/(G.dx*G.dy)

    def get_nodal_mass_momentum(self, P, G):
        for i in range(len(self)):
            for r, n in enumerate(self.n[i]):
                w = self.N[i, r]*self.m[i]
                G.m[n] += w
                G.q[n] += w*self.v[i]
                G.phim[n] += w*self.phi[i]
                G.pkm[n] += w*self.pk[i]

    def update_stress_strain(self, P, G):
        if self.material.law == 'rigid':
            return
        for i in range(len(self)):
            gradient = np.zeros((3, 3))
            for r, n in enumerate(self.n[i]):
                if G.m[n] > P.M_tol:
                    gradient += np.outer(G.q[n]/G.m[n], self.G[i, r])
            self.L[i] = gradient
            self.dstrain[i] = .5*(gradient+gradient.T)*P.dt
            self.strain[i] += self.dstrain[i]
            history = {key: getattr(self, key)[i].copy() for key in self.scalar_fields+self.tensor_fields}
            h = law(self.dstrain[i], self.stress[i], history, self.material, P,
                    self.rho[i], self.phi[i], self.initial_rho[i])
            for key, value in h.items():
                if key not in ('stress', 'strain', 'dstrain'):
                    getattr(self, key)[i] = value
            self.stress[i] += self.dstress[i]-(gradient@self.stress[i]+self.stress[i]@gradient.T)*P.dt

    def get_nodal_forces(self, P, G):
        body = np.array([P.g*np.sin(P.theta), P.g*np.cos(P.theta), 0])
        for i in range(len(self)):
            for r, n in enumerate(self.n[i]):
                G.fi[n] += self.V[i]*(self.stress[i]@self.G[i, r])
                G.fe[n] += self.N[i, r]*self.m[i]*body

    def deposit_diagnostics(self, P, G):
        for i in range(len(self)):
            fields = dict(pressure=self.pressure[i], gammadot=self.gammadot[i], sigmav=self.stress[i, 1, 1],
                          sigmah=self.stress[i, 0, 0], yieldfunction=self.yieldfunction[i],
                          mu=self.mu[i], I=self.I[i], eta=self.eta_limited[i] if self.material.law == 'pouliquen' else self.eta[i],
                          dev_stress=np.linalg.norm(self.dev_stress[i])/(np.sqrt(2) if self.material.law == 'pouliquen' else 1),
                          dev_stress_dot=np.linalg.norm(self.dev_dstress[i])/P.dt)
            for r, n in enumerate(self.n[i]):
                for key, value in fields.items():
                    getattr(G, key)[n] += self.N[i, r]*self.m[i]*value

    def move_material_points(self, P, G):
        if self.material.law == 'rigid':
            self.v[:] = self.prescribed_velocity(P)
            self.x += P.dt*self.v
            return
        for i in range(len(self)):
            self.rho[i] = 0
            for r, n in enumerate(self.n[i]):
                if G.m[n] > P.M_tol:
                    self.x[i] += P.dt*self.N[i, r]*G.q[n]/G.m[n]
                    self.v[i] += P.dt*self.N[i, r]*G.q_dot[n]/G.m[n]
                    self.rho[i] += self.N[i, r]*G.m[n]/G.V[n]
