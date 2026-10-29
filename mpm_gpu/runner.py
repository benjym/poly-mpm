"""Resident Metal timestep for one granular material (granular.step on the device).

Particle and grid state stay on the device between steps. Each step reads
the solver's state record between Newton/sweep batches and one record of
reduced scalars at the end, which supplies the diagnostics, work budgets and
timestep limits. NumPy arrays are written back only by snapshot(): at output
times, at the end of run() and after a standalone step(). Grid fields that
are deposits of particle diagnostics (granular.step's deposit_diagnostics,
grain fractions and shear rate) are formed from the snapshot with the CPU
deposit routines, using the start-of-step stencil exactly as the CPU step
does; they do not feed back into the mechanics.
"""

from weakref import WeakKeyDictionary

import numpy as np

from . import GPUError
from . import mechanics
from .primitives import MAX, MIN, SUM, ColumnReduction
from .solver import Solver
from .topology import join_positions, split_positions

MODULE = "runner.slang"
NODE_COLUMNS = (SUM,) * 8 + (MAX,)
ROW_COLUMNS = (SUM, MAX, SUM)
PACKING_COLUMNS = (MAX, SUM)
PARTICLE_COLUMNS = (SUM, SUM, SUM, SUM, SUM, SUM, SUM, MIN, MAX, MIN, MAX, MAX, MAX, SUM, SUM, SUM)
ROOT2 = np.sqrt(2.0)

_runtime = None
_runners = WeakKeyDictionary()


def runtime():
    """One Metal device per process; runners share it."""
    global _runtime
    if _runtime is None:
        from .runtime import MetalRuntime

        _runtime = MetalRuntime()
    return _runtime


def runner_for(P, G, L):
    """The cached runner for this particle system, rebuilt if the problem shape changed."""
    key = (P.G.nx, P.G.ny, bool(P.B.cyclic_lr), len(L.materials[0]), bool(P.B.wall_mu))
    runner = _runners.get(L)
    if runner is None or runner.key != key:
        runner = MetalRunner(P, G, L)
        _runners[L] = runner
    return runner


def check_callbacks(P):
    """Callbacks run between resident steps would see stale arrays; refuse them."""
    import initialise

    each = P.O.after_every_timestep
    if not (getattr(each, "__module__", "") == "plotting" and each.__name__ == "each_step"):
        raise ValueError(
            "Metal runs cannot call a custom after_every_timestep callback: particle arrays are "
            "resident on the device and only synchronised at output times"
        )
    if P.O.measure_energy:
        raise ValueError("Metal runs do not support O.measure_energy")
    controller = getattr(P.update_timestep, "_controller", P.update_timestep)
    if controller is not initialise.update_timestep:
        raise ValueError("Metal runs support only the built-in dynamic timestep controller")


class MetalRunner:
    def __init__(self, P, G, L):
        from granular import basis_volumes, wall_candidates

        self.rt = rt = runtime()
        self.key = (P.G.nx, P.G.ny, bool(P.B.cyclic_lr), len(L.materials[0]), bool(P.B.wall_mu))
        b = L.materials[0]
        self.count = n = len(b)
        if n == 0:
            raise ValueError("Metal runs require at least one particle")
        f32 = np.float32
        # Frictional walls: one static contact row per wall node (granular.wall_rows).
        empty = np.zeros(0, dtype=int)
        self.wall_nodes, self.wall_sides, blocks = (
            wall_candidates(P, G, 3) if P.B.wall_mu else (empty, empty, np.zeros((0, 3, 2)))
        )
        self.grid = grid = mechanics.device_grid(
            rt, P.G.nx, P.G.ny, G.dx, G.dy, periodic=P.B.cyclic_lr, wall_nodes=self.wall_nodes
        )
        T = grid.topology
        self.wall_blocks = rt.buffer(blocks if len(blocks) else np.zeros(6), dtype=f32)

        def vector(count, dtype=f32):
            return rt.buffer(count=count, dtype=dtype)

        node_row = np.zeros(grid.nodes, dtype=np.uint32)
        node_row[T.nodes] = np.arange(grid.stress_rows)
        node_row = node_row[T.canonical]
        self.node_row_host = node_row
        self.node_row = rt.buffer(node_row)
        volumes = np.concatenate((basis_volumes(P, G)[T.nodes], np.zeros(len(self.wall_nodes))))
        self.volumes = rt.buffer(volumes, dtype=f32)
        self.fixed = vector(grid.nodes, np.uint32)
        self.wall_changed = vector(1, np.uint32)
        # Particle state.
        self.particles = mechanics.Particles(
            n, *(vector(n, np.uint32) for _ in range(2)), vector(n), vector(n), vector(n)
        )
        self.mass, self.velocity, self.gradient = vector(n), vector(2 * n), vector(4 * n)
        self.deformation, self.jacobian, self.reference_volume = vector(9 * n), vector(n), vector(n)
        self.volume, self.density, self.shift = vector(n), vector(n), vector(n)
        self.stress, self.stress_change, self.strain = vector(4 * n), vector(4 * n), vector(4 * n)
        self.work, self.work_compensation = vector(n), vector(n)
        # Workspaces. Two operators alternate: the end-of-step binning is the
        # next step's start-of-step binning.
        self.operators = [mechanics.Operator(rt, grid, n), mechanics.Operator(rt, grid, n)]
        self.current = 0
        self.assembled = False
        self.preconditioner = mechanics.Preconditioner(rt, grid)
        self.precondition_status = vector(1, np.uint32)
        self.nodal = mechanics.GridState(rt, grid)
        self.gap, self.initial, self.previous = (
            vector(3 * grid.rows),
            vector(3 * grid.rows),
            vector(3 * grid.rows),
        )
        self.internal, self.momentum, self.rate = (
            vector(2 * grid.nodes),
            vector(2 * grid.nodes),
            vector(2 * grid.nodes),
        )
        self.status = vector(1, np.uint32)
        self.node_columns = vector(len(NODE_COLUMNS) * grid.nodes)
        self.row_columns = vector(len(ROW_COLUMNS) * grid.rows)
        self.packing_columns = vector(len(PACKING_COLUMNS) * grid.rows)
        self.particle_columns = vector(len(PARTICLE_COLUMNS) * n)
        self.reductions = [
            ColumnReduction(rt, grid.nodes, NODE_COLUMNS),
            ColumnReduction(rt, grid.rows, ROW_COLUMNS),
            ColumnReduction(rt, grid.rows, PACKING_COLUMNS),
            ColumnReduction(rt, n, PARTICLE_COLUMNS),
        ]
        self.segregation = None
        if P.segregate:
            from .segregation import Segregation

            self.segregation = Segregation(rt, grid, P, G, n)
            self.reductions += self.segregation.reductions
        self.record_size = sum(r.columns for r in self.reductions)
        self.record = vector(self.record_size)
        S = b.material
        self.solver = Solver(
            rt,
            grid,
            mu=S.mu,  # sqrt(2 / dimension) = 1 for the intrinsic 2D stress space
            dt=P.dt,
            tolerance=P.granular_tolerance,
            max_iterations=P.granular_max_iterations,
            wall_mu=P.B.wall_mu or 0.0,
        )
        self.last = None

    # ---- host <-> device ---------------------------------------------------------

    def upload(self, P, G, L):
        """Copy the host state to the device (start of run() or standalone step())."""
        from .mechanics import fixed_components

        rt, b, f32 = self.rt, L.materials[0], np.float32
        S = b.material
        lo, hi = np.array([P.G.x_m, P.G.y_m]), np.array([P.G.x_M, P.G.y_M])
        x = b.x[:, :2].copy()
        if P.B.cyclic_lr:
            x[:, 0] = P.G.x_m + (x[:, 0] - P.G.x_m) % (P.G.x_M - P.G.x_m)
        if not np.isfinite(b.x).all():
            b.check_finite(P)
        if np.any((x < lo) | (x > hi)):
            raise ValueError(f"material {b.material_id}: particle outside domain before cell lookup")
        if np.any(b.v[:, 2] != 0) or np.any(b.L[:, 2, :] != 0) or np.any(b.L[:, :, 2] != 0):
            raise ValueError("Metal runs require planar particle velocities and velocity gradients")
        cells, fractions = split_positions(x, lo, [G.dx, G.dy])
        p = self.particles
        rt.write(p.cell_x, cells[:, 0].view(np.uint32), np.uint32)
        rt.write(p.cell_y, cells[:, 1].view(np.uint32), np.uint32)
        rt.write(p.fraction_x, fractions[:, 0], f32)
        rt.write(p.fraction_y, fractions[:, 1], f32)
        rt.write(p.volume, b.m / S.rho_s, f32)
        rt.write(self.mass, b.m, f32)
        rt.write(self.velocity, b.v[:, :2], f32)
        rt.write(self.gradient, b.L[:, :2, :2], f32)
        rt.write(self.deformation, b.F, f32)
        rt.write(self.jacobian, b.J, f32)
        rt.write(self.reference_volume, b.V0, f32)
        rt.write(self.volume, b.V, f32)
        rt.write(self.density, b.rho, f32)
        rt.write(self.stress, _planar(b.stress), f32)
        rt.write(self.stress_change, _planar(b.dstress), f32)
        for buffer in (self.strain, self.work, self.work_compensation):
            rt.write(buffer, np.zeros(buffer.size // 4), f32)
        fixed = fixed_components(G.boundary_v, G.boundary_h, P.G.nx, P.B.separable, bool(P.B.wall_mu))
        rt.write(self.fixed, fixed, np.uint32)
        previous = getattr(G, "granular_contact", None)
        rows = self.grid.topology.nodes
        previous = np.zeros((len(rows), 3)) if previous is None else previous[rows]
        walls = getattr(G, "granular_wall_contact", None)
        walls = (
            np.zeros((len(self.wall_nodes), 3)) if walls is None else walls[self.wall_nodes, self.wall_sides]
        )
        rt.write(self.previous, np.vstack((previous, walls)), f32)
        # FP64 bases and constants; the device holds increments.
        self.base = dict(
            strain=b.strain.copy(),
            material_work=b.material_work.copy(),
            plastic_work=b.plastic_work.copy(),
        )
        self.constants = dict(
            mass=float(b.m.sum()),
            species=np.sum(b.m[:, None] * b.phi, axis=0),
            elastic=float(b.elastic_energy.sum()),
            viscous=float(b.viscous_work.sum()),
            material=float(b.material_work.sum()),
            plastic=float(b.plastic_work.sum()),
            solid_volume=float(np.sum(b.m / S.rho_s)),
            size=b.phi @ P.G.s,
        )
        self.speed = float(np.linalg.norm(b.v, axis=1).max(initial=0))
        if self.segregation is not None:
            self.segregation.upload(b, P)
        # Previous-step segregation maxima, as stability.limits reads G.grad_pk.
        self.segregation_rates = dict(
            gradient=float(np.linalg.norm(G.grad_pk, axis=1).max(initial=0)),
            component=float(np.max(np.abs(G.grad_pk), initial=0)),
            diffusivity=0.0,
        )
        if P.segregate and getattr(G, "shear_rate", None) is not None and len(G.shear_rate) == len(G.m):
            from integrators import _diffusivity

            self.segregation_rates["diffusivity"] = float(_diffusivity(G, P, G.phi).max(initial=0))
        self.assembled = False

    def snapshot(self, P, G, L):
        """Write the device state into the NumPy particle and grid arrays."""
        from interpolation import stencil
        from transfers import divide, scatter

        if self.last is None:
            return
        rt, b, f32 = self.rt, L.materials[0], np.float32
        S = b.material
        n, N = self.count, self.grid.nodes
        p = self.particles
        dt = self.last["dt"]
        origin, spacing = [P.G.x_m, P.G.y_m], [G.dx, G.dy]

        def read(buffer, shape=None, dtype=f32):
            value = rt.read(buffer, dtype).astype(np.float64 if dtype == f32 else dtype)
            return value if shape is None else value.reshape(shape)

        cells = np.column_stack(
            (read(p.cell_x, dtype=np.uint32).view(np.int32), read(p.cell_y, dtype=np.uint32).view(np.int32))
        )
        fractions = np.column_stack((read(p.fraction_x), read(p.fraction_y)))
        b.x[:, :2] = join_positions(cells, fractions, origin, spacing)
        b.x[:, 2] = 0
        b.v[:, :2] = read(self.velocity, (n, 2))
        b.v[:, 2] = 0
        L2 = read(self.gradient, (n, 2, 2))
        b.L = np.zeros((n, 3, 3))
        b.L[:, :2, :2] = L2
        b.F = read(self.deformation, (n, 3, 3))
        b.J = read(self.jacobian)
        b.V = read(self.volume)
        b.rho = read(self.density)
        b.stress = _tensor(read(self.stress, (n, 4)))
        b.dstress = _tensor(read(self.stress_change, (n, 4)))
        b.strain = self.base["strain"] + _tensor(read(self.strain, (n, 4)))
        b.dstrain = 0.5 * (b.L + b.L_T) * dt
        work = read(self.work) - read(self.work_compensation)
        b.material_work = self.base["material_work"] + work
        b.plastic_work = self.base["plastic_work"] + work
        size = self.constants["size"]
        if self.segregation is not None:
            b.phi, start = self.segregation.particles(read)
            size = start @ P.G.s
        b.pressure = np.trace(b.stress, axis1=1, axis2=2) / 3
        b.dp = np.trace(b.dstress, axis1=1, axis2=2) / 3
        b.sigma_kk = 3 * b.pressure
        eye = np.eye(3)
        b.dev_stress = b.stress - b.pressure[:, None, None] * eye
        b.dev_dstress = b.dstress - np.trace(b.dstress, axis1=1, axis2=2)[:, None, None] * eye / 3
        rate = 0.5 * (b.L + b.L_T)
        dev = rate.copy()
        dev[:, :2, :2] -= np.trace(rate, axis1=1, axis2=2)[:, None, None] * np.eye(2) / 2
        b.gammadot = ROOT2 * np.linalg.norm(dev, axis=(1, 2))
        phase_pressure = -np.trace(b.stress, axis1=1, axis2=2) / 3 * S.rho_s / b.rho
        b.I = np.divide(
            size * b.gammadot * np.sqrt(S.rho_s),
            np.sqrt(np.maximum(phase_pressure, 0)),
            out=np.zeros(n),
            where=phase_pressure > 0,
        )
        b.mu.fill(S.mu)
        b.yieldfunction = np.linalg.norm(b.dev_stress, axis=(1, 2)) / ROOT2 + S.mu * b.pressure
        # Start-of-step stencil of the last step, as the CPU leaves it.
        operator = self.operators[1 - self.current]
        keys = read(operator.keys, dtype=np.uint32).astype(np.int64)
        columns = P.G.nx - 1
        local = np.column_stack((read(operator.local_x), read(operator.local_y)))
        start = join_positions(np.column_stack((keys % columns, keys // columns)), local, origin, spacing)
        start = np.column_stack((start, np.zeros(n)))
        b.n, b.N, b.G, b.offsets = stencil(start, P, G)
        b.n_star = b.n[:, 0]
        b._apic_matrices = None
        # Grid.
        G.wipe(P)
        G.m = read(self.nodal.mass)
        G.get_valid_nodes(P)
        G.q_reference[:, :2] = read(self.nodal.reference, (N, 2))
        G.q[:, :2] = read(self.momentum, (N, 2))
        G.q_dot[:, :2] = read(self.rate, (N, 2))
        G.fi[:, :2] = read(self.internal, (N, 2))
        G.fe[:, :2] = G.m[:, None] * self.last["body"]
        contact = read(self.solver.contact, (self.grid.rows, 3))
        G.granular_contact = contact[self.node_row_host]
        if len(self.wall_nodes):
            from granular import WALL_SIDES

            G.granular_wall_contact = np.zeros((len(G.m), len(WALL_SIDES), 3))
            G.granular_wall_contact[self.wall_nodes, self.wall_sides] = contact[self.grid.stress_rows :]
        weight = b.N * b.m[:, None]
        if self.segregation is None:
            scatter(G.phim, b.n, weight[:, :, None] * b.phi[:, None, :])
        scatter(G.pkm, b.n, weight * b.pk[:, None])
        saved = P.dt
        P.dt = dt
        try:
            b.deposit_diagnostics(P, G)
        finally:
            P.dt = saved
        if P.B.cyclic_lr:
            G.make_cyclic(
                P, G, (["phim"] if self.segregation is None else []) + ["pkm", *G.diagnostic_fields]
            )
        G.phi = divide(G.phim, G.m[:, None])
        G.s_bar = G.phi @ P.G.s
        G.pk = divide(G.pkm, G.m)
        if self.segregation is not None:
            self.segregation.grid_fields(P, G, read)
        G.calculate_gammadot(P, G, P.smooth_gamma_dot)
        if self.segregation is not None:
            G.update_pk(P, G)
        self._budgets(G, self.last)
        G.granular_stats = dict(self.last["stats"])

    # ---- the timestep ---------------------------------------------------------

    def limits(self, P):
        """stability.limits for granular materials, from the device speed record."""
        h = min(P.G.dx, P.G.dy)
        result = {"advection": h / self.speed} if self.speed else {}
        for S in P.S:
            if hasattr(S, "critical_time"):
                value = float(S.critical_time(P))
                if not np.isfinite(value) or value <= 0:
                    raise ValueError("material critical_time must be positive and finite")
                result["material_0"] = value
        if P.segregate:
            rate = P.c * (P.G.s_M / P.G.s_m - 1) * self.segregation_rates["gradient"]
            if rate:
                result["segregation_advection"] = h / rate
            diffusion = self.segregation_rates["diffusivity"]
            if diffusion:
                result["segregation_diffusion"] = 1 / (2 * diffusion * (1 / P.G.dx**2 + 1 / P.G.dy**2))
        return {name: float(value * P.CFL) for name, value in result.items()}

    def update_timestep(self, P):
        """The built-in dynamic controller (initialise.update_timestep and limits)."""
        h = min(P.G.dx, P.G.dy)
        speed = self.last["nodal_speed"] if self.last else 0.0
        limits = [0.1 * h / speed] if speed > 0 else []
        if P.segregate:
            rate = P.c * (P.G.s_M / P.G.s_m - 1) * self.segregation_rates["component"]
            if rate > 0:
                limits.append(h / rate)
        for S in P.S:
            if hasattr(S, "critical_time"):
                limits.append(S.critical_time(P))
        candidate = P.CFL * min(limits) if limits else P.dt
        if not np.isfinite(candidate) or candidate <= 0:
            raise ValueError(f"invalid dynamic timestep {candidate}")
        P.dt = candidate
        restrictions = self.limits(P)
        if restrictions:
            P.dt = min(P.dt, min(restrictions.values()))

    def step(self, P, G, L):
        """Advance one timestep on the device; host arrays are not updated."""
        if not np.isfinite(P.dt) or P.dt <= 0:
            raise ValueError("dt must be finite and positive")
        P.timestep_limits = self.limits(P)
        if P.timestep_policy == "enforce" and P.dt > min(P.timestep_limits.values(), default=P.dt) * (
            1 + 1e-12
        ):
            raise ValueError("granular dt exceeds transport timestep estimate")
        if L.diagnostics is not None and not L.diagnostics.started:
            L.diagnostics.begin(P, L)
        P.update_forces()
        rt, grid, dt = self.rt, self.grid, float(P.dt)
        b = L.materials[0]
        S = b.material
        body = np.array([P.g * np.sin(P.theta), P.g * np.cos(P.theta)])
        n = self.count
        encoder = rt.device.create_command_encoder()
        operator = self.operators[self.current]
        if not self.assembled:
            mechanics.assemble(rt, encoder, grid, self.particles, operator)
        nodal = mechanics.p2g(
            rt,
            encoder,
            grid,
            operator,
            self.fixed,
            self.mass,
            self.velocity,
            self.gradient,
            dt=dt,
            m_tol=P.M_tol,
            body=body,
            affine=P.transfer == "apic_gradient",
            rough=P.B.roughness,
            state=self.nodal,
        )
        mechanics.wall_blocks(
            rt,
            encoder,
            grid,
            operator,
            self.wall_blocks,
            self.node_row,
            self.volumes,
            nodal.mass,
            P.M_tol,
            S.packing_max,
        )
        rt.dispatch(
            encoder,
            "packing_gap",
            grid.rows,
            module=MODULE,
            row_count=grid.rows,
            stress_rows=grid.stress_rows,
            packing_max=float(S.packing_max),
            dt=dt,
            occupied=operator.occupied,
            volumes=self.volumes,
            previous=self.previous,
            gap=self.gap,
            initial=self.initial,
        )
        solver = self.solver
        solver.dt = dt
        solver.tolerance = float(P.granular_tolerance)
        # One-sided walls: re-solve after clamping released wall nodes that the
        # solve drives into their wall (granular.step).
        for _ in range(grid.nodes if P.B.separable else 1):
            rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=self.precondition_status)
            pre = mechanics.preconditioner(
                rt,
                encoder,
                grid,
                operator,
                nodal.inverse_mass,
                dt,
                status=self.precondition_status,
                workspace=self.preconditioner,
            )
            result = solver.solve(
                operator, pre, nodal.trial, nodal.inverse_mass, self.gap, self.initial, encoder=encoder
            )
            if not P.B.separable or not result.converged:
                break
            encoder = rt.device.create_command_encoder()
            rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=self.wall_changed)
            mechanics.wall_recheck(
                rt, encoder, grid, self.fixed, nodal, solver.u, self.wall_changed, rough=P.B.roughness
            )
            rt.submit(encoder)
            if not rt.read(self.wall_changed, np.uint32)[0]:
                break
            # The re-solve starts from the attempt it replaces (same rows).
            encoder = rt.device.create_command_encoder()
            rt.dispatch(
                encoder,
                "copy_record",
                3 * grid.rows,
                module=MODULE,
                count=3 * grid.rows,
                offset=0,
                source=solver.contact,
                record=self.initial,
            )
        status = int(rt.read(operator.status, np.uint32)[0])
        if status & 1:
            raise ValueError(f"material {b.material_id}: particle outside domain before cell lookup")
        if status:
            raise GPUError(f"Metal operator assembly status {status}")
        if not result.converged:
            raise RuntimeError(
                f"granular velocity solve did not converge: residual={result.residual:.3g}, "
                f"sweeps={result.sweeps}, status={result.reasons}"
            )
        encoder = rt.device.create_command_encoder()
        rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=self.status)
        self._encode_finish(encoder, P, G, operator, body, dt)
        rt.submit(encoder)
        record = rt.read(self.record, np.float32).astype(np.float64)
        flags = int(rt.read(self.status, np.uint32)[0])
        if flags & 64:
            raise FloatingPointError(
                f"granular, material {b.material_id}, timestep {P.tstep}: invalid deformation determinant"
            )
        if flags & 2 or not np.isfinite(record).all():
            raise FloatingPointError(
                f"granular, material {b.material_id}, timestep {P.tstep}: nonfinite particle state"
            )
        self.current = 1 - self.current
        self.assembled = True
        self.last = self._record(P, record, result, body, dt)
        self.speed = self.last["speed"]
        self._budgets(G, self.last)
        G.granular_stats = dict(self.last["stats"])
        event = self._boundary_event(P, L)
        if event:
            # Inlets and outlets change the particle count: apply them to the
            # synchronised host state, as granular.step does, then rebuild.
            self.snapshot(P, G, L)
            L.process_boundaries(P, G)
        if L.diagnostics is not None:
            now = self.last["inventory"]
            before = dict(now, angular_momentum=self.last["angular_before"])
            if event:
                from diagnostics import inventory

                now = inventory(P, L)
            stats = G.granular_stats
            diagnostics_stats = dict(
                maximum_packing_after=self.last["packing_after"],
                integrated_overpacking_fraction=self.last["overpacking"],
                maximum_particle_packing=self.last["particle_packing"],
            )
            stats.update(diagnostics_stats)
            self.last["stats"].update(diagnostics_stats)
            L.diagnostics.finish(P, G, L, before, now=now, support=(self.last["incomplete_support"], 0))
        last = self.last
        if event:
            self.__init__(P, G, L)
            _runners[L] = self
            self.upload(P, G, L)
            self.last, self.speed = last, last["speed"]
        return last

    def _boundary_event(self, P, L):
        """Whether process_boundaries will add or remove particles this step."""
        if self.last["leaving"]:
            return True
        for side, template in (("right", L.templates[0]), ("top", L.templates[0])):
            rate = getattr(template, "inlet_rate", None)
            if getattr(P.B, "inlet_" + side) and rate and P.t and P.t % rate < P.dt:
                return True
        return False

    def _encode_finish(self, encoder, P, G, operator, body, dt):
        rt, grid, n = self.rt, self.grid, self.count
        b_transfer = P.transfer
        periodic = int(grid.periodic)
        self.solver.encode_finish(encoder)
        contact = self.solver.contact
        u = self.solver.u
        rt.dispatch(
            encoder,
            "grid_finish",
            grid.nodes,
            module=MODULE,
            node_count=grid.nodes,
            nx=grid.nx,
            periodic=periodic,
            dx=grid.dx,
            dy=grid.dy,
            x_m=float(P.G.x_m),
            y_m=float(P.G.y_m),
            dt=dt,
            m_tol=float(P.M_tol),
            body_x=float(body[0]),
            body_y=float(body[1]),
            transpose_rows=grid.transpose_rows,
            transpose_slots=grid.transpose_slots,
            gradients=operator.gradients,
            contact=contact,
            mass=self.nodal.mass,
            reference=self.nodal.reference,
            velocity=u,
            internal=self.internal,
            momentum=self.momentum,
            rate=self.rate,
            columns=self.node_columns,
        )
        rt.dispatch(
            encoder,
            "contact_rows",
            grid.rows,
            module=MODULE,
            row_count=grid.rows,
            neighbors=grid.neighbors,
            gradients=operator.gradients,
            contact=contact,
            velocity=u,
            occupied=operator.occupied,
            volumes=self.volumes,
            previous=self.previous,
            columns=self.row_columns,
        )
        mechanics.g2p(
            rt,
            encoder,
            grid,
            operator,
            self.particles,
            self.nodal,
            u,
            self.velocity,
            self.gradient,
            dt=dt,
            m_tol=P.M_tol,
            transfer=b_transfer,
            pic_fraction=P.pic_fraction,
            shift=self.shift,
        )
        mechanics.deform(
            rt,
            encoder,
            n,
            self.gradient,
            self.reference_volume,
            self.mass,
            self.deformation,
            self.jacobian,
            self.volume,
            self.density,
            self.status,
            dt=dt,
        )
        S = P.S[0]
        rt.dispatch(
            encoder,
            "particle_update",
            n,
            module=MODULE,
            count=n,
            nx=grid.nx,
            periodic=periodic,
            dt=dt,
            m_tol=float(P.M_tol),
            rho_s=float(S.rho_s),
            keys=operator.keys,
            local_x=operator.local_x,
            local_y=operator.local_y,
            mass=self.nodal.mass,
            contact=contact,
            node_row=self.node_row,
            gradient=self.gradient,
            density=self.density,
            volume=self.volume,
            stress=self.stress,
            stress_change=self.stress_change,
            strain=self.strain,
            work=self.work,
            work_compensation=self.work_compensation,
        )
        segregation = []
        if self.segregation is not None:
            segregation = self.segregation.encode(
                encoder,
                P,
                operator,
                self.nodal,
                u,
                self.mass,
                self.stress,
                self.gradient,
                self.density,
                self.status,
                dt,
            )
        rt.dispatch(
            encoder,
            "particle_inventory",
            n,
            module=MODULE,
            count=n,
            nx=grid.nx,
            periodic=periodic,
            dx=grid.dx,
            dy=grid.dy,
            x_m=float(P.G.x_m),
            y_m=float(P.G.y_m),
            m_tol=float(P.M_tol),
            rho_s=float(S.rho_s),
            affine=int(P.transfer == "apic_gradient"),
            outlets=int(bool(P.B.outlet_left)) + 2 * int(bool(P.B.outlet_bottom)),
            keys=operator.keys,
            local_x=operator.local_x,
            local_y=operator.local_y,
            nodal_mass=self.nodal.mass,
            mass=self.mass,
            velocity=self.velocity,
            gradient=self.gradient,
            cell_x=self.particles.cell_x,
            cell_y=self.particles.cell_y,
            fraction_x=self.particles.fraction_x,
            fraction_y=self.particles.fraction_y,
            shift=self.shift,
            jacobian=self.jacobian,
            density=self.density,
            work=self.work,
            work_compensation=self.work_compensation,
            columns=self.particle_columns,
        )
        # End-of-step binning: packing diagnostics now, start-of-step operator next step.
        after = mechanics.assemble(rt, encoder, grid, self.particles, self.operators[1 - self.current])
        rt.dispatch(
            encoder,
            "packing_rows",
            grid.rows,
            module=MODULE,
            row_count=grid.rows,
            packing_max=float(S.packing_max),
            occupied=after.occupied,
            volumes=self.volumes,
            columns=self.packing_columns,
        )
        offset = 0
        results = [
            reduction.encode(encoder, columns, self.status)
            for reduction, columns in zip(
                self.reductions[:4],
                (self.node_columns, self.row_columns, self.packing_columns, self.particle_columns),
            )
        ] + segregation
        for reduction, result in zip(self.reductions, results):
            rt.dispatch(
                encoder,
                "copy_record",
                reduction.columns,
                module=MODULE,
                count=reduction.columns,
                offset=offset,
                source=result,
                record=self.record,
            )
            offset += reduction.columns

    def _record(self, P, values, result, body, dt):
        node, row, packing, particle, seg = np.split(values, np.cumsum([9, 3, 2, 16]))
        c = self.constants
        species = c["species"].tolist()
        if self.segregation is not None:
            self.segregation_rates = dict(gradient=seg[0], component=seg[1], diffusivity=seg[2])
            species = seg[3:].tolist()
        work = particle[6]
        inventory = dict(
            mass=c["mass"],
            species_mass=species,
            momentum=[particle[0], particle[1], 0.0],
            angular_momentum=[0.0, 0.0, particle[2]],
            kinetic=particle[4],
            affine_kinetic=particle[5],
            elastic_energy=c["elastic"],
            plastic_work=c["plastic"] + work,
            viscous_work=c["viscous"],
            material_work=c["material"] + work,
            minimum_J=particle[7],
            maximum_J=particle[8],
            minimum_density=particle[9],
            maximum_density=particle[10],
        )
        solver = self.solver
        stats = dict(
            residual=result.residual,
            evaluations=0,
            sweeps=result.sweeps,
            converged=result.converged,
            method="velocity",
            backend="metal",
            operator_nnz=9 * self.grid.rows,
            operator_storage_bytes=self.grid.rows * 9 * 12,
            operator_format="device compact Q1 slots",
            assembled_delassus=False,
            residual_norm="conservative block Euclidean",
            colors=len(solver.colors),
            preconditioner_storage_bytes=self.grid.rows * 20,
            acceleration_steps=result.acceleration_steps,
            newton_steps=result.newton_steps,
            anderson_rejected=result.anderson_rejected,
            newton_attempts=result.newton_attempts,
            newton_rejected=result.newton_rejected,
            gmres_iterations=result.gmres_iterations,
            submissions=result.submissions,
            active_nodes=int(round(row[2])),
            stress_dimension=2,
            stress_basis="nodal",
            viscosity_trial_dissipation=0.0,
            contact_power=row[0],
            maximum_packing=row[1],
        )
        solid = c["solid_volume"]
        return dict(
            dt=dt,
            body=body,
            inventory=inventory,
            angular_before=[0.0, 0.0, particle[3]],
            speed=particle[11],
            nodal_speed=node[8],
            particle_packing=particle[12],
            incomplete_support=int(round(particle[13])),
            device_mass=particle[14],
            leaving=int(round(particle[15])),
            packing_after=packing[0],
            overpacking=packing[1] / solid if solid else 0.0,
            constraint_impulse=[node[0], node[1], 0.0],
            constraint_angular_impulse=[0.0, 0.0, node[2]],
            constraint_work=node[3],
            external_work=node[4],
            external_impulse=[node[5], node[6], 0.0],
            external_angular_impulse=[0.0, 0.0, node[7]],
            stats=stats,
        )

    @staticmethod
    def _budgets(G, record):
        G.constraint_impulse = np.array(record["constraint_impulse"])
        G.constraint_angular_impulse = np.array(record["constraint_angular_impulse"])
        G.constraint_work = float(record["constraint_work"])
        G.external_work = float(record["external_work"])
        G.external_impulse = np.array(record["external_impulse"])
        G.external_angular_impulse = np.array(record["external_angular_impulse"])
        G.force_scaling_work = 0.0
        G.force_scaling_impulse = np.zeros(3)
        G.force_scaling_angular_impulse = np.zeros(3)


def _planar(tensor):
    """(xx, yy, zz, xy) of planar symmetric 3x3 tensors."""
    return np.column_stack((tensor[:, 0, 0], tensor[:, 1, 1], tensor[:, 2, 2], tensor[:, 0, 1]))


def _tensor(planar):
    out = np.zeros((len(planar), 3, 3))
    out[:, 0, 0], out[:, 1, 1], out[:, 2, 2] = planar[:, 0], planar[:, 1], planar[:, 2]
    out[:, 0, 1] = out[:, 1, 0] = planar[:, 3]
    return out


def step(P, G, L):
    """Standalone step(): upload the host state, advance, write the result back."""
    runner = runner_for(P, G, L)
    runner.upload(P, G, L)
    runner.step(P, G, L)
    runner.snapshot(P, G, L)
