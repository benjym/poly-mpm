import os
from numpy import *
import matplotlib.pyplot as plt


class Params:
    def __init__(self, mode):
        self.dt = 1e-4  # timestep (s)
        self.savetime = 0.1
        self.t_f = 0.1  # final time (s)
        self.nt = int(self.t_f / self.dt)  # number of timesteps
        self.max_q = 1e2
        self.max_g = 0.0
        self.pressure = "non-zero"
        self.theta = 0.0 * pi / 180.0  # slope angle (degrees)
        self.G = Grid_Params()
        self.B = Boundary_Params()
        self.S = Solid_Params(self, self.G)
        self.O = Output_Params(self.nt, self.S.n)
        self.supername = "im/pouliquen_unit_test"

    #         self.damping = True

    def update_forces(self):
        self.g = 0
        t_c = self.t_f / 2.0
        if self.t < t_c:
            self.q_h = self.t * self.max_q / t_c
            self.q_v = self.t * self.max_q / t_c
        else:
            self.q_h = (
                self.max_q + (self.t - t_c) * self.max_q / t_c
            )  # + 1e2*self.max_q*((self.t-t_c))#**3#/(self.t_f-t_c))**3
            self.q_v = self.max_q - (self.t - t_c) * self.max_q / t_c
            # self.q_v = 2*self.max_q - self.t/t_c*self.max_q


#         self.q_h = self.t*self.max_q/t_c
#         self.q_v = minimum(self.q_h,self.max_q)


class Grid_Params:
    def __init__(self):
        self.x_m = 0.0  # (m)
        self.x_M = 0.1  # (m)
        self.y_m = 0.0  # (m)
        self.y_M = 0.1  # (m)
        self.nx = 2
        self.ny = 2
        self.x = linspace(self.x_m, self.x_M, self.nx)
        self.y = linspace(self.y_m, self.y_M, self.ny)
        self.dx = self.x[1] - self.x[0]  # grid spacing (m)
        self.dy = self.y[1] - self.y[0]  # grid spacing (m)

        self.ns = 1
        self.s = [0.01]  # (m)


class Boundary_Params:
    def __init__(self):
        #         self.has_right = True
        #         self.has_left = True
        self.vertical_force = True
        self.horizontal_force = True


class Solid_Params:
    def __init__(self, P, G):
        self.X = []
        self.Y = []
        self.n = 0

        self.rho = 2700.0  # density (kg/m^3)
        self.packing = 0.6  # packing fraction
        self.rho_s = self.rho / self.packing  # solid density
        self.phi = ones([G.ns]) / float(G.ns)
        self.law = "pouliquen"
        self.mu_0 = tan(deg2rad(20.9))  # 0.3
        self.mu_1 = tan(deg2rad(32.76))
        self.delta_mu = self.mu_1 - self.mu_0
        self.I_0 = 0.279
        self.mu_v = 0
        self.E = 1e8
        self.nu = 0.4  # poissons ratio
        self.K = self.E / (3.0 * (1.0 - 2.0 * self.nu))
        self.G = self.E / (2.0 * (1.0 + self.nu))
        self.eta_max = 1e12  # *self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)
        # self.eta_max = 1e2*self.rho*sqrt(-P.max_g*(G.y_M-G.y_m)**3)

        self.pts_per_cell = 3
        self.x = (G.nx - 1) * self.pts_per_cell  # particles in x direction
        self.y = (G.ny - 1) * self.pts_per_cell  # particles in y direction
        gap = array((G.dx, G.dy)) / (2 * self.pts_per_cell)
        xp = linspace(G.x_m + gap[0], G.x_M - gap[0], self.x)
        yp = linspace(G.y_m + gap[1], G.y_M - gap[1], self.y)
        X = tile(xp, self.y)
        Y = repeat(yp, self.x)
        for i in range(self.x * self.y):
            self.X.append(X[i])
            self.Y.append(Y[i])
            self.n += 1
        self.A = (G.x_M - G.x_m) * (G.y_M - G.y_m) / self.n  # area (m^2)

        # Advective terms
        elastic_wave_speed = sqrt(self.K / self.rho)
        distance = minimum(G.dx, G.dy)
        critical_adv_time = distance / elastic_wave_speed
        # Diffusive terms
        l_0 = xp[1] - xp[0]  # initial distance between material points
        critical_diff_time = l_0**2 * self.rho / self.eta_max

        critical_time = minimum(critical_adv_time, critical_diff_time)
        CFL = critical_time / P.dt

        if CFL < 2:
            print("WARNING: STABILITY IS GUARANTEED TO BE POOR")
        print("CFL from elastic wave speed: " + str(critical_adv_time / P.dt))
        print("CFL from momentum diffusion: " + str(critical_diff_time / P.dt))
        print("Current CFL: " + str(CFL))


class Output_Params:
    """Array-valued constitutive history for the loading-path example."""

    def __init__(self, nt, n):
        self.history = []

    def after_every_timestep(self, L, P, G):
        batch = L.materials[0]
        self.history.append(
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
        from pathlib import Path
        from plotting import Plotting

        if not self.history:
            return
        history = {name: array([row[name] for row in self.history]) for name in self.history[0]}
        directory = Path(P.save_dir) / "data"
        directory.mkdir(parents=True, exist_ok=True)
        savez_compressed(directory / "constitutive_history.npz", **history)
        fig, axes = plt.subplots(2, 3, figsize=(12, 7), layout="constrained")
        time = history["t"]
        for ax, name in zip(axes.flat, ("I", "pressure", "mu", "eta_limited", "gammadot")):
            ax.plot(time, history[name][:, 0])
            ax.set(xlabel="time", ylabel=name)
        tau = linalg.norm(history["dev_stress"][:, 0], axis=(1, 2)) / sqrt(2)
        axes.flat[-1].plot(time, tau)
        axes.flat[-1].set(xlabel="time", ylabel="deviatoric stress norm")
        Plotting().savefig(P, "0")
