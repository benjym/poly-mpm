"""Run-time inventories and discrete work accounting; no hidden energy correction."""

import numpy as np


def inventory(P, L):
    species = np.zeros(P.G.ns)
    mass = kinetic = affine_kinetic = elastic = plastic = viscous = work = 0.0
    momentum = np.zeros(3)
    angular = np.zeros(3)
    minimum_J, maximum_J = np.inf, -np.inf
    minimum_rho, maximum_rho = np.inf, -np.inf
    for b in L.materials:
        mass += b.m.sum()
        species += np.sum(b.m[:, None] * b.phi, axis=0)
        momentum += np.sum(b.m[:, None] * b.v, axis=0)
        angular += np.sum(np.cross(b.x, b.m[:, None] * b.v), axis=0)
        kinetic += 0.5 * np.sum(b.m * np.sum(b.v**2, axis=1))
        if P.transfer == "apic" and b.material.law != "rigid":
            B = b.affine_momentum
            angular += np.sum(
                b.m[:, None]
                * np.stack(
                    (B[:, 2, 1] - B[:, 1, 2], B[:, 0, 2] - B[:, 2, 0], B[:, 1, 0] - B[:, 0, 1]), axis=1
                ),
                axis=0,
            )
            _, inverse = b.apic_matrices()
            affine_kinetic += 0.5 * np.sum(b.m * np.einsum("nij,njk,nik->n", B, inverse, B))
        elif P.transfer == "apic_gradient":
            # Current bilinear support covariance, not the stale pre-advection
            # stencil: this is the affine state carried into the next P2G.
            h = np.array([(P.G.x_M - P.G.x_m) / (P.G.nx - 1), (P.G.y_M - P.G.y_m) / (P.G.ny - 1)])
            coordinate = (b.x[:, :2] - [P.G.x_m, P.G.y_m]) / h
            fraction = coordinate - np.floor(coordinate)
            diagonal = h**2 * fraction * (1 - fraction)
            B = np.zeros_like(b.L)
            B[:, :, :2] = b.L[:, :, :2] * diagonal[:, None, :]
            angular += np.sum(
                b.m[:, None]
                * np.stack(
                    (B[:, 2, 1] - B[:, 1, 2], B[:, 0, 2] - B[:, 2, 0], B[:, 1, 0] - B[:, 0, 1]), axis=1
                ),
                axis=0,
            )
            affine_kinetic += 0.5 * np.sum(b.m * np.sum(b.L * B, axis=(1, 2)))
        elastic += b.elastic_energy.sum()
        plastic += b.plastic_work.sum()
        viscous += b.viscous_work.sum()
        work += b.material_work.sum()
        minimum_J = min(minimum_J, b.J.min(initial=np.inf))
        maximum_J = max(maximum_J, b.J.max(initial=-np.inf))
        minimum_rho = min(minimum_rho, b.rho.min(initial=np.inf))
        maximum_rho = max(maximum_rho, b.rho.max(initial=-np.inf))
    return dict(
        mass=float(mass),
        species_mass=species.tolist(),
        momentum=momentum.tolist(),
        angular_momentum=angular.tolist(),
        kinetic=float(kinetic),
        affine_kinetic=float(affine_kinetic),
        elastic_energy=float(elastic),
        plastic_work=float(plastic),
        viscous_work=float(viscous),
        material_work=float(work),
        minimum_J=float(minimum_J) if mass else None,
        maximum_J=float(maximum_J) if mass else None,
        minimum_density=float(minimum_rho) if mass else None,
        maximum_density=float(maximum_rho) if mass else None,
    )


class Diagnostics:
    """Compact history; boundary/inlet budgets and residuals remain explicit.

    External and constraint work use midpoint grid velocities. Residual includes
    transfer/time integration error and, for rate-form plastic/fluid laws, the
    difference between their quadratic elastic-energy estimate and exact energy.
    Periodic angular momentum requires an unwrapped domain, and is diagnostic
    only here. Prescribed rigid motion is externally controlled.
    """

    def __init__(self, P, L):
        self.initial = inventory(P, L)
        self.history = []
        self.last = None
        self.started = False
        self.external_work = self.constraint_work = self.force_scaling_work = 0.0
        self.boundary_mass = 0.0
        self.boundary_species = np.zeros(P.G.ns)
        self.constraint_impulse = np.zeros(3)
        self.external_impulse = np.zeros(3)
        self.boundary_momentum = np.zeros(3)
        self.boundary_angular = np.zeros(3)
        self.boundary_energy = 0.0
        self.external_angular_impulse = np.zeros(3)
        self.constraint_angular_impulse = np.zeros(3)
        self.force_scaling_impulse = np.zeros(3)
        self.force_scaling_angular_impulse = np.zeros(3)

    def begin(self, P, L):
        if not self.started:
            self.initial = inventory(P, L)
            self.started = True

    def finish(self, P, G, L, boundary_before):
        now = inventory(P, L)
        self.boundary_mass += now["mass"] - boundary_before["mass"]
        self.boundary_species += np.asarray(now["species_mass"]) - boundary_before["species_mass"]
        self.boundary_momentum += np.asarray(now["momentum"]) - boundary_before["momentum"]
        self.boundary_angular += np.asarray(now["angular_momentum"]) - boundary_before["angular_momentum"]
        self.boundary_energy += sum(
            now[key] - boundary_before[key]
            for key in ("kinetic", "affine_kinetic", "elastic_energy", "plastic_work", "viscous_work")
        )
        self.external_work += G.external_work
        self.constraint_work += G.constraint_work
        self.force_scaling_work += G.force_scaling_work
        self.constraint_impulse += G.constraint_impulse
        self.external_impulse += G.external_impulse
        self.external_angular_impulse += G.external_angular_impulse
        self.constraint_angular_impulse += G.constraint_angular_impulse
        self.force_scaling_impulse += G.force_scaling_impulse
        self.force_scaling_angular_impulse += G.force_scaling_angular_impulse
        supported = 0
        gradient_missing = 0
        for b in L.materials:
            if len(b):
                gradient_missing += int(np.count_nonzero(b.gradient_support_missing))
                coverage = np.sum(b.N * G.valid_nodes[b.n], axis=1)
                supported += int(np.count_nonzero(coverage < 1 - 1e-10))
        row = dict(
            time=float(P.t + P.dt),
            timestep=int(P.tstep + 1),
            mechanics=P.mechanics,
            interpolation=P.interpolation,
            quadrature=P.quadrature,
            transfer=P.transfer,
            maximum_stencil_entries=max((b.N.shape[1] for b in L.materials), default=0),
            **now,
            boundary_mass=float(self.boundary_mass),
            boundary_species=self.boundary_species.tolist(),
            boundary_momentum=self.boundary_momentum.tolist(),
            boundary_angular_momentum=self.boundary_angular.tolist(),
            boundary_energy=float(self.boundary_energy),
            mass_error=float(now["mass"] - self.initial["mass"] - self.boundary_mass),
            species_mass_error=(
                np.asarray(now["species_mass"]) - self.initial["species_mass"] - self.boundary_species
            ).tolist(),
            external_work=float(self.external_work),
            constraint_work=float(self.constraint_work),
            force_scaling_work=float(self.force_scaling_work),
            constraint_impulse=self.constraint_impulse.tolist(),
            external_impulse=self.external_impulse.tolist(),
            constraint_angular_impulse=self.constraint_angular_impulse.tolist(),
            external_angular_impulse=self.external_angular_impulse.tolist(),
            force_scaling_angular_impulse=self.force_scaling_angular_impulse.tolist(),
            timestep_within_estimate=bool(
                P.dt <= min(getattr(P, "timestep_limits", {}).values(), default=P.dt) * (1 + 1e-12)
            ),
            particles_with_incomplete_support=supported,
            particles_with_incomplete_gradient_support=gradient_missing,
            timestep_limits=getattr(P, "timestep_limits", {}),
        )
        if P.solver == "granular":
            row["granular_solver"] = G.granular_stats.copy()
        initial_energy = (
            self.initial["kinetic"] + self.initial["affine_kinetic"] + self.initial["elastic_energy"]
        )
        row["momentum_residual"] = (
            np.asarray(now["momentum"])
            - self.initial["momentum"]
            - self.external_impulse
            - self.constraint_impulse
            - self.force_scaling_impulse
            - self.boundary_momentum
        ).tolist()
        row["force_scaling_impulse"] = self.force_scaling_impulse.tolist()
        row["angular_momentum_residual"] = (
            np.asarray(now["angular_momentum"])
            - self.initial["angular_momentum"]
            - self.external_angular_impulse
            - self.constraint_angular_impulse
            - self.force_scaling_angular_impulse
            - self.boundary_angular
        ).tolist()
        row["angular_budget_has_unwrapped_coordinates"] = not P.B.cyclic_lr
        row["stored_energy_exact"] = P.mechanics == "updated" and all(
            b.material.law in ("elastic", "viscous", "bingham", "pouliquen", "linear_mu") for b in L.materials
        )
        row["energy_residual"] = float(
            now["kinetic"]
            + now["affine_kinetic"]
            + now["elastic_energy"]
            - initial_energy
            + now["plastic_work"]
            - self.initial["plastic_work"]
            + now["viscous_work"]
            - self.initial["viscous_work"]
            - self.external_work
            - self.constraint_work
            - self.force_scaling_work
            - self.boundary_energy
        )
        row["energy_budget_closed"] = not bool(
            self.boundary_mass
            or np.any(self.boundary_species)
            or any(
                getattr(P.B, flag) for flag in ("inlet_right", "inlet_top", "outlet_left", "outlet_bottom")
            )
            or any(b.material.law == "rigid" and len(b) for b in L.materials)
        )
        self.last = row
        if (P.tstep + 1) % P.diagnostics_stride == 0 or P.t + P.dt >= P.t_f:
            self.history.append(row)
        return row

    def complete(self):
        if self.last is not None and (
            not self.history or self.history[-1]["timestep"] != self.last["timestep"]
        ):
            self.history.append(self.last)

    def save(self, directory):
        from pathlib import Path
        import json

        path = Path(directory) / "data" / "diagnostics.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(dict(initial=self.initial, history=self.history), indent=2, allow_nan=False) + "\n"
        )
