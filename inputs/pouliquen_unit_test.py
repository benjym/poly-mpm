"""One cell of mu(I) material under boundary tractions.

Both tractions rise together to 100 Pa at t = 0.05 s; after that the horizontal traction keeps
rising while the vertical one falls, shearing the cell. The output records the constitutive
history of the first material point.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from plotting import Plotting

_MAX_Q = 1e2  # traction at the end of loading (Pa)
_T_C = 0.05  # end of loading (s)


class Params:
    def __init__(self, args):
        self.dt = 1e-4  # (s)
        self.savetime = 0.1  # (s)
        self.t_f = 0.1  # (s)
        self.output_dir = "im/pouliquen_unit_test"
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.S = Solid_Params(self.G)
        self.O = Output_Params()


class Grid_Params:
    def __init__(self):
        self.x_min, self.x_max = 0.0, 0.1  # (m)
        self.y_min, self.y_max = 0.0, 0.1  # (m)
        self.nx = 2
        self.ny = 2
        self.grain_sizes = [0.01]  # (m)


class Boundary_Params:
    def traction_vertical(self, t):
        return t * _MAX_Q / _T_C if t < _T_C else _MAX_Q - (t - _T_C) * _MAX_Q / _T_C

    def traction_horizontal(self, t):
        return t * _MAX_Q / _T_C if t < _T_C else _MAX_Q + (t - _T_C) * _MAX_Q / _T_C


class Solid_Params:
    def __init__(self, G):
        self.law = "pouliquen"
        self.rho = 2700.0  # (kg/m^3)
        self.phi_flow = 0.6  # packing fraction
        self.rho_s = self.rho / self.phi_flow  # solid density (kg/m^3)
        self.phi = np.ones(1)
        self.mu_s = np.tan(np.radians(20.9))
        self.mu_d = np.tan(np.radians(32.76))
        self.I_0 = 0.279
        E, nu = 1e8, 0.4  # elastic modulus (Pa) and Poisson's ratio
        self.K = E / (3.0 * (1.0 - 2.0 * nu))  # bulk modulus (Pa)
        self.eta_max = 1e12  # (Pa s)

        pts_per_cell = 3
        dx = (G.x_max - G.x_min) / (G.nx - 1)
        dy = (G.y_max - G.y_min) / (G.ny - 1)
        nx = (G.nx - 1) * pts_per_cell
        ny = (G.ny - 1) * pts_per_cell
        gap = np.array((dx, dy)) / (2 * pts_per_cell)
        xp = np.linspace(G.x_min + gap[0], G.x_max - gap[0], nx)
        yp = np.linspace(G.y_min + gap[1], G.y_max - gap[1], ny)
        self.X = np.tile(xp, ny)
        self.Y = np.repeat(yp, nx)
        self.n = len(self.X)
        self.A = (G.x_max - G.x_min) * (G.y_max - G.y_min) / self.n  # (m^2)


class Output_Params:
    """Array-valued constitutive history for the loading-path example."""

    def __init__(self):
        self._history = []

    def after_every_timestep(self, L, P, G):
        batch = L.materials[0]
        self._history.append(
            dict(
                t=P.t,
                **{
                    name: getattr(batch, name).copy()
                    for name in (
                        "mu",
                        "I",
                        "eta",
                        "eta_limited",
                        "pressure",
                        "dstrain",
                        "dev_stress",
                        "gammadot",
                    )
                },
            )
        )

    def at_end(self, L, P, G):
        if not self._history:
            return
        history = {name: np.array([row[name] for row in self._history]) for name in self._history[0]}
        directory = Path(P.save_dir) / "data"
        directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(directory / "constitutive_history.npz", **history)
        fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
        time = history["t"]
        for ax, name in zip(axes.flat, ("I", "pressure", "mu", "eta_limited", "gammadot")):
            ax.plot(time, history[name][:, 0])
            ax.set(xlabel="time", ylabel=name)
        tau = np.linalg.norm(history["dev_stress"][:, 0], axis=(1, 2)) / np.sqrt(2)
        axes.flat[-1].plot(time, tau)
        axes.flat[-1].set(xlabel="time", ylabel="deviatoric stress norm")
        Plotting().savefig(P, "0")
