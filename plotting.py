"""Array-aware plotting, saved fields, and legacy output-flag adapters."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from transfers import divide, nodal_velocity


def values(L, name):
    return np.concatenate([getattr(batch, name) for batch in L.materials], axis=0)


class Plotting:
    def savefig(self, P, name, dpi=150):
        path = Path(P.save_dir) / (name+'.png')
        path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(path, dpi=dpi, bbox_inches='tight')
        plt.close()

    def draw_grid(self, G):
        plt.plot(G.X, G.Y, '+', color='.7', markersize=3)
        plt.plot(G.X[G.boundary_h], G.Y[G.boundary_h], 'rx', markersize=3)
        plt.plot(G.X[G.boundary_v], G.Y[G.boundary_v], 'r+', markersize=3)

    draw_pretty_grid = draw_grid

    def draw_continuum(self, G, P):
        velocity = nodal_velocity(G, P)
        fields = (velocity[:, 0], velocity[:, 1], divide(G.m, G.V),
                  G.physical('pressure'), G.shear_rate, G.s_bar)
        names = ('u', 'v', 'density', 'pressure (tension positive)', 'shear rate', 'mean grain size')
        fig, axes = plt.subplots(2, 3, figsize=getattr(P.O, 'continuum_fig_size', (12, 7)), layout='constrained')
        for ax, name, data in zip(axes.flat, names, fields):
            masked = np.ma.masked_where(~G.valid_nodes, data).reshape(P.G.ny, P.G.nx)
            mesh = ax.pcolormesh(G.x_plot, G.y_plot, masked, shading='flat')
            ax.set(title=name, xlabel='x', ylabel='y', xlim=(P.G.x_m, P.G.x_M), ylim=(P.G.y_m, P.G.y_M))
            fig.colorbar(mesh, ax=ax)
        self.savefig(P, f'Continuum/{P.grid_save:05d}')
        P.grid_save += 1

    def draw_material_points(self, L, P, G, name=None):
        x, v, pressure = values(L, 'x'), values(L, 'v'), values(L, 'pressure')
        fig, axes = plt.subplots(1, 3, figsize=getattr(P.O, 'mp_fig_size', (12, 4)), layout='constrained')
        material = np.concatenate([np.full(len(b), i) for i, b in enumerate(L.materials)])
        for ax, title, data in zip(axes, ('material', 'speed', 'pressure'),
                                   (material, np.linalg.norm(v, axis=1), pressure)):
            points = ax.scatter(x[:, 0], x[:, 1], c=data, s=8)
            ax.set(title=title, xlim=(P.G.x_m, P.G.x_M), ylim=(P.G.y_m, P.G.y_M), aspect='equal')
            fig.colorbar(points, ax=ax)
        self.savefig(P, 'MP/'+(str(name) if name else f'{P.mp_save:05d}'))
        P.mp_save += 1

    def draw_gamma_dot(self, L, P, G):
        x = values(L, 'x')
        plt.figure()
        plt.scatter(x[:, 0], x[:, 1], c=values(L, 'gammadot'), s=8)
        plt.colorbar(label='strain rate')
        self.savefig(P, f'MP/gammadot_{P.mp_save:05d}')

    def draw_gsd_mp(self, L, P, G):
        x = values(L, 'x')
        phi = values(L, 'phi')
        fig, axes = plt.subplots(1, P.G.ns, squeeze=False, figsize=(4*P.G.ns, 4), layout='constrained')
        for j, ax in enumerate(axes.flat):
            points = ax.scatter(x[:, 0], x[:, 1], c=phi[:, j], vmin=0, vmax=1, s=8)
            ax.set(title=f'grain class {j}', xlim=(P.G.x_m, P.G.x_M), ylim=(P.G.y_m, P.G.y_M))
            fig.colorbar(points, ax=ax)
        self.savefig(P, f'GSD/MP_{P.mp_save:05d}')

    def draw_gsd_grid(self, L, P, G):
        fig, axes = plt.subplots(1, P.G.ns, squeeze=False, figsize=(4*P.G.ns, 4), layout='constrained')
        for j, ax in enumerate(axes.flat):
            mesh = ax.pcolormesh(G.x_plot, G.y_plot, G.phi[:, j].reshape(P.G.ny, P.G.nx), vmin=0, vmax=1)
            ax.set_title(f'grain class {j}')
            fig.colorbar(mesh, ax=ax)
        self.savefig(P, f'GSD/{P.grid_save:05d}')

    def _save(self, P, label, array, counter):
        directory = Path(P.save_dir)/'data'
        directory.mkdir(parents=True, exist_ok=True)
        index = getattr(P, counter, 0)
        np.save(directory/f'{label}_{index:05d}.npy', array)
        setattr(P, counter, index+1)

    def save_s_bar(self, L, P, G):
        self._save(P, 's_bar', G.s_bar.reshape(P.G.ny, P.G.nx), 's_bar_save')

    def save_density(self, L, P, G):
        self._save(P, 'density', divide(G.m, G.V).reshape(P.G.ny, P.G.nx), 'density_save')

    def save_u(self, L, P, G):
        velocity = nodal_velocity(G, P)
        index = getattr(P, 'u_save', 0)
        self._save(P, 'u', velocity[:, 0].reshape(P.G.ny, P.G.nx), 'u_save')
        P.v_save = index
        self._save(P, 'v', velocity[:, 1].reshape(P.G.ny, P.G.nx), 'v_save')

    def save_phi_MP(self, L, P, G):
        self._save(P, 'MP_phi', np.column_stack((values(L, 'x')[:, :2], values(L, 'phi'))), 'mp_phi_save')

    def measure_E(self, L, P, G):
        """Save array-valued macro strain and tangent ratios for each material."""
        directory = Path(P.save_dir)/'data'
        directory.mkdir(parents=True, exist_ok=True)
        for i, b in enumerate(L.materials):
            macro = divide(b.x0-b.x, [P.G.x_M-P.G.x_m, P.G.y_M-P.G.y_m, 1])
            np.savez(directory/f'stiffness_material_{i}.npz', macro_strain=macro,
                     strain=b.strain, stress=b.stress, incremental_modulus=divide(b.dstress, b.dstrain))

    def draw_energy(self, P):
        energy = np.asarray(P.O.energy)
        if not energy.size:
            return
        plt.figure()
        for i, name in enumerate(('kinetic', 'potential', 'strain'), start=1):
            plt.plot(energy[:, 0], energy[:, i], label=name)
        plt.legend()
        plt.xlabel('time')
        self.savefig(P, 'energy')
        directory = Path(P.save_dir)/'data'
        directory.mkdir(parents=True, exist_ok=True)
        np.save(directory/'energy.npy', energy)


def install_callbacks(P):
    plot = Plotting()
    O = P.O
    if O.measure_energy:
        O.energy = []  # append per step; dynamic timesteps need no preallocation

    def default_output(L, P, G):
        if O.plot_continuum:
            plot.draw_continuum(G, P)
        for flag, method in (('plot_material_points', 'draw_material_points'),
                             ('plot_gsd_mp', 'draw_gsd_mp'), ('plot_gsd_grid', 'draw_gsd_grid'),
                             ('save_s_bar', 'save_s_bar'), ('save_density', 'save_density'),
                             ('save_u', 'save_u'), ('save_phi_MP', 'save_phi_MP')):
            if getattr(O, flag):
                getattr(plot, method)(L, P, G)
        if O.plot_gsd_debug:
            plot.draw_gsd_grid(L, P, G)
        if getattr(O, 'plot_gamma_dot', False):
            plot.draw_gamma_dot(L, P, G)

    def each_step(L, P, G):
        if O.measure_energy:
            O.energy.append(L.measure_elastic_energy(P, G))

    def at_end(L, P, G):
        default_output(L, P, G)
        if O.measure_energy:
            plot.draw_energy(P)
        if O.measure_stiffness:
            plot.measure_E(L, P, G)

    if not hasattr(O, 'at_start'):
        O.at_start = lambda L, P, G: Path(P.save_dir).mkdir(parents=True, exist_ok=True)
    if not hasattr(O, 'after_every_timestep'):
        O.after_every_timestep = getattr(O, 'after_ever_timestep', each_step)
    O.after_ever_timestep = O.after_every_timestep
    if not hasattr(O, 'after_every_savetime'):
        O.after_every_savetime = default_output
    if not hasattr(O, 'at_end'):
        O.at_end = at_end
