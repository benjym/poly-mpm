"""Resident compact Coulomb solver (granular_mixed.solve_velocity) on Metal.

All iterates, Anderson histories, TSQR/SVD factors, Krylov bases, Givens
rotations and acceptance decisions stay on the device. The host encodes
bounded batches (25 sweeps, or one Newton polish) and reads a 64-byte state
record between them. Every workspace is allocated once, so re-encoding a
batch allocates nothing.
"""

from dataclasses import dataclass

import numpy as np

from .primitives import Reduction

MODULE = "solver.slang"
HISTORY = 13
RESTART = 64
CYCLES = 4
FRACTIONS = (1.0, 0.5, 0.25, 0.125)
GMRES_RTOL = 1e-5
# Singular values below this fraction of the largest are dropped from the
# Anderson fit (the FP64 reference uses rcond=1e-10).
ANDERSON_RCOND = 1e-5
BATCH = 25
# Host chunking: stop encoding gated-off work once the device reports it done.
SWEEP_CHUNK = 5
SINGLE_ITERATION_ROWS = 1024
ARNOLDI_CHUNK = 16
KRYLOV_SIZE = 4420

STATE = (
    "flags",
    "result",
    "sweeps",
    "count",
    "head",
    "accelerated",
    "anderson_rejected",
    "newton_steps",
    "newton_rejected",
    "gmres_iterations",
    "krylov",
    "breakdowns",
    "newton_attempts",
    "gmres_converged",
)
REQUESTED = 256
STATUS = {
    2: "nonfinite residual or input",
    16: "invalid local metric",
    32: "local sliding solve unconverged",
    128: "compression without free velocity support",
    256: "sweep budget exhausted",
    1024: "nonfinite Krylov value (Newton candidate discarded)",
    2048: "nonfinite Anderson fit (candidate discarded)",
}
FATAL = 2 | 16 | 32 | 128 | 256


class SolverFailure(RuntimeError):
    def __init__(self, result):
        super().__init__(
            f"Metal granular velocity solve did not converge: residual={result.residual:.3g}, "
            f"sweeps={result.sweeps}, status={result.reasons}"
        )
        self.result = result


@dataclass
class SolveResult:
    converged: bool
    residual: float
    sweeps: int
    acceleration_steps: int
    anderson_rejected: int
    newton_attempts: int
    newton_steps: int
    newton_rejected: int
    gmres_iterations: int
    gmres_converged: int
    breakdowns: int
    status: int
    submissions: int
    magnitude: float
    norm: float
    scale: float

    @property
    def reasons(self):
        return [text for bit, text in STATUS.items() if self.status & bit]

    def require(self):
        if not self.converged:
            raise SolverFailure(self)
        return self


class Solver:
    """Workspaces sized by the grid; reusable across solves on that grid."""

    def __init__(self, rt, grid, *, mu, dt, tolerance=1e-4, max_iterations=200):
        if not np.isfinite(tolerance) or tolerance < 1e-5:
            raise ValueError("Metal FP32 requires a solver tolerance >= 1e-5")
        self.rt, self.grid = rt, grid
        self.mu, self.dt, self.tolerance = float(mu), float(dt), float(tolerance)
        self.budget = int(max_iterations) * 10
        R, N = grid.rows, grid.nodes
        n = self.size = 3 * R
        f32 = np.float32

        def vector(count):
            return rt.buffer(count=count, dtype=f32)

        colors = np.asarray(grid.topology.colors)
        order = np.argsort(colors, kind="stable")
        counts = np.bincount(colors)
        self.colors = [(int(o), int(c)) for o, c in zip(np.cumsum(counts) - counts, counts) if c]
        self.color_rows = rt.buffer(order.astype(np.uint32), label="color_rows")
        self.color_offsets = rt.buffer(np.cumsum([0] + [c for _, c in self.colors]).astype(np.uint32))
        # Small systems run each iteration as two single-workgroup kernels.
        self.single_iteration = grid.rows <= SINGLE_ITERATION_ROWS
        self.state = rt.buffer(count=16, label="solver_state")
        self.status = rt.buffer(count=1, label="solver_status")
        self.candidate_status = rt.buffer(count=1)
        self.scalars = vector(4)
        self.krylov = vector(KRYLOV_SIZE)
        self.z, self.candidate = vector(n), vector(n)
        self.u, self.candidate_u = vector(2 * N), vector(2 * N)
        self.contact = vector(n)
        self.row_values = [vector(R), vector(R)]
        self.rows = Reduction(rt, R, operation="max")
        self.compression = Reduction(rt, R, operation="max")
        self.starts, self.steps = vector(HISTORY * n), vector(HISTORY * n)
        self.coefficients = vector(HISTORY - 1)
        self.matrix = vector(HISTORY * n)
        # TSQR levels: each 256-row tile reduces to a 13-row R factor.
        self.tsqr = []
        rows = n
        while True:
            groups = (rows + 255) // 256
            self.tsqr.append((rows, groups, vector(HISTORY * HISTORY * groups)))
            if groups == 1:
                break
            rows = HISTORY * groups
        self.frozen, self.rhs, self.x, self.r = vector(6 * R), vector(n), vector(n), vector(n)
        self.w, self.du = vector(n), vector(2 * N)
        self.basis = vector((RESTART + 1) * n)
        self.chunks = (n + 255) // 256
        self.partials = vector((RESTART + 1) * self.chunks)
        self.first, self.second, self.squares = vector(RESTART + 1), vector(RESTART + 1), vector(1)
        # Arnoldi inner products: at most 64 chunks of 64 threads per vector.
        self.per_thread = -(-n // (64 * 64))
        self.groups = -(-n // (64 * self.per_thread))
        self.passes = [vector((RESTART + 1) * self.groups) for _ in range(3)]
        self.bindings = None

    # ---- encoding helpers ----------------------------------------------------

    def _dispatch(self, encoder, entry, threads, **bindings):
        self.rt.dispatch(encoder, entry, threads, module=MODULE, **bindings)

    def _single(self, encoder, entry, **bindings):
        self._dispatch(encoder, entry, 1, **bindings)

    def _fill(self, encoder, buffer, count, value=0):
        self.rt.dispatch(encoder, "fill_u32", count, count=count, value=value, output=buffer)

    def _refresh(self, encoder, gate, z, u):
        b = self.bindings
        self._dispatch(
            encoder,
            "velocity_refresh",
            self.grid.nodes,
            gate=gate,
            node_count=self.grid.nodes,
            dt=self.dt,
            state=self.state,
            scalars=self.scalars,
            transpose_rows=self.grid.transpose_rows,
            transpose_slots=self.grid.transpose_slots,
            gradients=b["gradients"],
            inverse_mass=b["inverse_mass"],
            trial=b["trial"],
            z=z,
            u=u,
        )

    def _residual(self, encoder, gate, z, u, status):
        b, R = self.bindings, self.grid.rows
        self._dispatch(
            encoder,
            "residual_rows",
            R,
            gate=gate,
            row_count=R,
            mu=self.mu,
            state=self.state,
            scalars=self.scalars,
            neighbors=self.grid.neighbors,
            gradients=b["gradients"],
            gap=b["gap"],
            z=z,
            u=u,
            residual=self.row_values[0],
        )
        return self.rows.encode(encoder, self.row_values[0], status)

    def _combine(self, encoder, gate, output, a, alpha=1.0, b=None, beta=0.0, size=None):
        size = self.size if size is None else size
        self._dispatch(
            encoder,
            "combine",
            size,
            gate=gate,
            size=size,
            alpha=float(alpha),
            beta=float(beta),
            state=self.state,
            a=a,
            b=a if b is None else b,
            output=output,
        )

    def _candidate(self, encoder, gate, newton):
        """Evaluate the candidate stress in self.candidate; commit it if it reduces the residual."""
        self._fill(encoder, self.candidate_status, 1)
        self._refresh(encoder, gate, self.candidate, self.candidate_u)
        result = self._residual(encoder, gate, self.candidate, self.candidate_u, self.candidate_status)
        self._single(
            encoder,
            "candidate_accept",
            gate=gate,
            newton=int(newton),
            residual=result,
            candidate_status=self.candidate_status,
            state=self.state,
            scalars=self.scalars,
        )
        self._combine(encoder, 4, self.z, self.candidate)
        self._combine(encoder, 4, self.u, self.candidate_u, size=2 * self.grid.nodes)
        self._single(
            encoder, "candidate_finish", tolerance=self.tolerance, state=self.state, scalars=self.scalars
        )

    def _dots(self, encoder, gate, vectors, include_self, w, output):
        total = vectors + include_self
        self._dispatch(
            encoder,
            "dot_partial",
            total * self.chunks,
            gate=gate,
            size=self.size,
            chunks=self.chunks,
            vectors=vectors,
            total=total,
            state=self.state,
            basis=self.basis,
            w=w,
            partials=self.partials,
        )
        self._dispatch(
            encoder,
            "dot_finish",
            total,
            gate=gate,
            chunks=self.chunks,
            total=total,
            state=self.state,
            partials=self.partials,
            output=output,
        )

    def _jv(self, encoder, gate, source, offset, output):
        b = self.bindings
        self._dispatch(
            encoder,
            "jv_force",
            self.grid.nodes,
            gate=gate,
            node_count=self.grid.nodes,
            dt=self.dt,
            offset=offset,
            state=self.state,
            scalars=self.scalars,
            transpose_rows=self.grid.transpose_rows,
            transpose_slots=self.grid.transpose_slots,
            gradients=b["gradients"],
            inverse_mass=b["inverse_mass"],
            source=source,
            du=self.du,
        )
        self._dispatch(
            encoder,
            "jv_rows",
            self.grid.rows,
            gate=gate,
            row_count=self.grid.rows,
            mu=self.mu,
            offset=offset,
            state=self.state,
            neighbors=self.grid.neighbors,
            gradients=b["gradients"],
            frozen=self.frozen,
            source=source,
            du=self.du,
            output=output,
        )

    # ---- batches ---------------------------------------------------------------

    def encode_setup(self, encoder, operator, preconditioner, trial, inverse_mass, gap, initial=None):
        """Scaling, warm start, first exact velocity and residual.

        trial and inverse_mass hold two components per canonical node, gap
        three per stress row; initial is a previous contact stress (not scaled).
        """
        self.bindings = dict(
            gradients=operator.gradients,
            trial=trial,
            inverse_mass=inverse_mass,
            gap=gap,
            diagonal=preconditioner.diagonal,
            rotation=preconditioner.rotation,
        )
        R = self.grid.rows
        for buffer, count in ((self.state, 16), (self.status, 1), (self.scalars, 4)):
            self._fill(encoder, buffer, count)
        magnitude, compression = self.row_values
        self._dispatch(
            encoder,
            "solver_rows_setup",
            R,
            row_count=R,
            neighbors=self.grid.neighbors,
            gradients=operator.gradients,
            trial=trial,
            gap=gap,
            magnitude=magnitude,
            compression=compression,
        )
        largest = self.rows.encode(encoder, magnitude, self.status)
        compressive = self.compression.encode(encoder, compression, self.status)
        self._single(
            encoder,
            "solver_setup",
            magnitude=largest,
            compression=compressive,
            norm=preconditioner.norm,
            state=self.state,
            status=self.status,
            scalars=self.scalars,
        )
        self._dispatch(
            encoder,
            "stress_initial",
            R,
            row_count=R,
            mu=self.mu,
            warm=int(initial is not None),
            state=self.state,
            scalars=self.scalars,
            initial=self.candidate if initial is None else initial,
            z=self.z,
        )
        self._refresh(encoder, 0, self.z, self.u)
        result = self._residual(encoder, 1, self.z, self.u, self.status)
        self._single(
            encoder,
            "solver_start",
            tolerance=self.tolerance,
            residual=result,
            state=self.state,
            status=self.status,
            scalars=self.scalars,
        )

    def encode_sweeps(self, encoder, iterations=BATCH, request=True):
        """Colored sweeps, exact refresh and safeguarded Anderson, then a Newton request."""
        b, R, n = self.bindings, self.grid.rows, self.size
        if self.single_iteration:
            for _ in range(iterations):
                self._iteration_single(encoder)
            if request:
                self._single(
                    encoder,
                    "newton_request",
                    tolerance=self.tolerance,
                    state=self.state,
                    scalars=self.scalars,
                )
            return
        for _ in range(iterations):
            self._single(encoder, "iteration_begin", budget=self.budget, state=self.state, status=self.status)
            self._dispatch(
                encoder, "history_begin", n, size=n, state=self.state, z=self.z, starts=self.starts
            )
            rows = dict(
                mu=self.mu,
                dt=self.dt,
                state=self.state,
                scalars=self.scalars,
                color_rows=self.color_rows,
                neighbors=self.grid.neighbors,
                gradients=b["gradients"],
                gap=b["gap"],
                diagonal=b["diagonal"],
                rotation=b["rotation"],
                inverse_mass=b["inverse_mass"],
                status=self.status,
                z=self.z,
                u=self.u,
            )
            for offset, count in self.colors:
                self._dispatch(encoder, "sweep_color", count, color_offset=offset, color_count=count, **rows)
            # Periodic refresh avoids accumulating incremental roundoff.
            self._refresh(encoder, 1, self.z, self.u)
            result = self._residual(encoder, 1, self.z, self.u, self.status)
            self._dispatch(
                encoder,
                "history_end",
                n,
                size=n,
                state=self.state,
                z=self.z,
                starts=self.starts,
                steps=self.steps,
            )
            self._single(
                encoder,
                "sweep_end",
                tolerance=self.tolerance,
                residual=result,
                state=self.state,
                status=self.status,
                scalars=self.scalars,
            )
            self._encode_anderson(encoder)
        if request:
            self._single(
                encoder, "newton_request", tolerance=self.tolerance, state=self.state, scalars=self.scalars
            )

    def _iteration_single(self, encoder):
        """One iteration as iteration_front, the TSQR/SVD dispatches and iteration_back."""
        b = self.bindings
        common = dict(
            mu=self.mu,
            dt=self.dt,
            node_count=self.grid.nodes,
            row_count=self.grid.rows,
            state=self.state,
            scalars=self.scalars,
            neighbors=self.grid.neighbors,
            gradients=b["gradients"],
            gap=b["gap"],
            inverse_mass=b["inverse_mass"],
            transpose_rows=self.grid.transpose_rows,
            transpose_slots=self.grid.transpose_slots,
            trial=b["trial"],
            z=self.z,
            u=self.u,
            starts=self.starts,
            steps=self.steps,
        )
        self._dispatch(
            encoder,
            "iteration_front",
            1024,
            budget=self.budget,
            tolerance=self.tolerance,
            colors=len(self.colors),
            status=self.status,
            color_offsets=self.color_offsets,
            color_rows=self.color_rows,
            diagonal=b["diagonal"],
            rotation=b["rotation"],
            matrix=self.matrix,
            **common,
        )
        self._encode_tsqr(encoder)
        self._dispatch(
            encoder,
            "iteration_back",
            1024,
            tolerance=self.tolerance,
            coefficients=self.coefficients,
            candidate=self.candidate,
            candidate_u=self.candidate_u,
            **common,
        )

    def _encode_tsqr(self, encoder):
        source = self.matrix
        for rows, groups, output in self.tsqr:
            self._dispatch(
                encoder,
                "tsqr_tile",
                256 * groups,
                gate=1 | 2,
                rows=rows,
                state=self.state,
                input=source,
                output=output,
            )
            source = output
        self._single(
            encoder,
            "anderson_solve",
            rcond=ANDERSON_RCOND,
            factor=source,
            state=self.state,
            status=self.status,
            coefficients=self.coefficients,
        )

    def _encode_anderson(self, encoder):
        n, gate = self.size, 1 | 2
        self._dispatch(
            encoder, "anderson_matrix", n, size=n, state=self.state, steps=self.steps, matrix=self.matrix
        )
        source = self.matrix
        for rows, groups, output in self.tsqr:
            self._dispatch(
                encoder,
                "tsqr_tile",
                256 * groups,
                gate=gate,
                rows=rows,
                state=self.state,
                input=source,
                output=output,
            )
            source = output
        self._single(
            encoder,
            "anderson_solve",
            rcond=ANDERSON_RCOND,
            factor=source,
            state=self.state,
            status=self.status,
            coefficients=self.coefficients,
        )
        self._dispatch(
            encoder,
            "anderson_candidate",
            self.grid.rows,
            row_count=self.grid.rows,
            mu=self.mu,
            state=self.state,
            coefficients=self.coefficients,
            starts=self.starts,
            steps=self.steps,
            z=self.z,
            candidate=self.candidate,
        )
        self._candidate(encoder, gate, newton=False)

    def encode_newton(self, encoder):
        """Matrix-free Newton: GMRES(64), four cycles, CGS2, then backtracking."""
        self._newton_rows(encoder)
        self.encode_krylov(encoder)
        self.encode_search(encoder)
        self._single(encoder, "newton_end", state=self.state)

    def _newton_chunked(self):
        """encode_newton as several submissions: Arnoldi steps go in chunks and
        encoding stops once the device has closed GMRES for the cycle."""
        rt = self.rt
        encoder = rt.device.create_command_encoder()
        self._newton_rows(encoder)
        self._fill(encoder, self.x, self.size)
        submissions = 0
        for cycle in range(CYCLES):
            self._cycle_start(encoder, cycle)
            for j in range(RESTART):
                self._arnoldi(encoder, j)
                if (j + 1) % ARNOLDI_CHUNK and j + 1 < RESTART:
                    continue
                rt.submit(encoder)
                submissions += 1
                flags = self.read_state()["flags"]
                encoder = rt.device.create_command_encoder()
                if not flags & 32:
                    break
            self._cycle_end(encoder)
            # Stop when Krylov is done (converged or breakdown) or Newton was abandoned.
            if flags & (128 | 8) != 8:
                break
        self.encode_search(encoder)
        self._single(encoder, "newton_end", state=self.state)
        rt.submit(encoder)
        return submissions + 1

    def _newton_rows(self, encoder):
        b, R = self.bindings, self.grid.rows
        self._single(encoder, "newton_begin", state=self.state)
        self._dispatch(
            encoder,
            "newton_rows",
            R,
            row_count=R,
            mu=self.mu,
            state=self.state,
            scalars=self.scalars,
            neighbors=self.grid.neighbors,
            gradients=b["gradients"],
            gap=b["gap"],
            z=self.z,
            u=self.u,
            frozen=self.frozen,
            rhs=self.rhs,
        )

    def encode_krylov(self, encoder):
        """Restarted GMRES for J x = -natural, from x = 0, at the frozen branch (one batch)."""
        self._fill(encoder, self.x, self.size)
        for cycle in range(CYCLES):
            self._cycle_start(encoder, cycle)
            for j in range(RESTART):
                self._arnoldi(encoder, j)
            self._cycle_end(encoder)

    def _cycle_start(self, encoder, cycle):
        newton = 1 | 8
        if cycle == 0:
            self._combine(encoder, newton, self.r, self.rhs)
        else:
            self._jv(encoder, newton, self.x, 0, self.w)
            self._combine(encoder, newton, self.r, self.rhs, 1.0, self.w, -1.0)
        self._dots(encoder, newton, 0, 1, self.r, self.squares)
        self._single(
            encoder,
            "gmres_start",
            cycle=cycle,
            rtol=GMRES_RTOL,
            squares=self.squares,
            state=self.state,
            status=self.status,
            krylov=self.krylov,
        )

    def _arnoldi(self, encoder, j):
        """One Arnoldi step in nine dispatches: J v_j (normalising v_j on the way),
        two classical Gram-Schmidt passes, the norm and the Givens update."""
        b, n, gmres = self.bindings, self.size, 1 | 8 | 32
        self._dispatch(
            encoder,
            "jv_force_normalise",
            max(n, self.grid.nodes),
            gate=gmres,
            node_count=self.grid.nodes,
            size=n,
            slot=j,
            dt=self.dt,
            state=self.state,
            scalars=self.scalars,
            transpose_rows=self.grid.transpose_rows,
            transpose_slots=self.grid.transpose_slots,
            gradients=b["gradients"],
            inverse_mass=b["inverse_mass"],
            krylov=self.krylov,
            source=self.r if j == 0 else self.w,
            basis=self.basis,
            du=self.du,
        )
        self._dispatch(
            encoder,
            "jv_rows",
            self.grid.rows,
            gate=gmres,
            row_count=self.grid.rows,
            mu=self.mu,
            offset=j * n,
            state=self.state,
            neighbors=self.grid.neighbors,
            gradients=b["gradients"],
            frozen=self.frozen,
            source=self.basis,
            du=self.du,
            output=self.w,
        )
        dots = dict(
            size=n, per_thread=self.per_thread, groups=self.groups, state=self.state, basis=self.basis
        )
        for h in (self.first, self.second):
            self._dispatch(
                encoder,
                "arnoldi_dots",
                64 * self.groups * (j + 1),
                gate=gmres,
                squared=0,
                w=self.w,
                partials=self.passes[0],
                **dots,
            )
            self._dispatch(
                encoder,
                "arnoldi_update",
                n,
                gate=gmres,
                size=n,
                vectors=j + 1,
                groups=self.groups,
                state=self.state,
                basis=self.basis,
                partials=self.passes[0],
                h=h,
                w=self.w,
            )
        self._dispatch(
            encoder,
            "arnoldi_dots",
            64 * self.groups,
            gate=gmres,
            squared=1,
            w=self.w,
            partials=self.passes[2],
            **dots,
        )
        self._single(
            encoder,
            "gmres_givens",
            step=j,
            groups=self.groups,
            first=self.first,
            second=self.second,
            squares=self.passes[2],
            state=self.state,
            status=self.status,
            krylov=self.krylov,
        )

    def _cycle_end(self, encoder):
        self._single(encoder, "gmres_solve", state=self.state, status=self.status, krylov=self.krylov)
        self._dispatch(
            encoder,
            "basis_update",
            self.size,
            size=self.size,
            state=self.state,
            krylov=self.krylov,
            basis=self.basis,
            x=self.x,
        )
        self._single(encoder, "gmres_end", state=self.state)

    def encode_search(self, encoder):
        """Backtracking on the unprojected Newton candidate, accepted on full-residual decrease."""
        search = 1 | 8 | 16
        for fraction in FRACTIONS:
            self._combine(encoder, search, self.candidate, self.z, 1.0, self.x, fraction)
            self._candidate(encoder, search, newton=True)

    def _normalise(self, encoder, gate, w, offset):
        self._dispatch(
            encoder,
            "normalise",
            self.size,
            gate=gate,
            size=self.size,
            offset=offset,
            state=self.state,
            krylov=self.krylov,
            w=w,
            basis=self.basis,
        )

    def encode_finish(self, encoder):
        self._dispatch(
            encoder,
            "final_contact",
            self.size,
            size=self.size,
            scalars=self.scalars,
            z=self.z,
            contact=self.contact,
        )

    # ---- host orchestration ----------------------------------------------------

    def read_state(self):
        return dict(zip(STATE, (int(v) for v in self.rt.read(self.state, np.uint32))))

    def solve(self, operator, preconditioner, trial, inverse_mass, gap, initial=None, *, encoder=None):
        """Run to convergence or explicit failure; returns a SolveResult.

        If `encoder` is given, the set-up is appended to it (for example after
        assembly and the preconditioner) and submitted here.
        """
        rt = self.rt
        if encoder is None:
            encoder = rt.device.create_command_encoder()
        self.encode_setup(encoder, operator, preconditioner, trial, inverse_mass, gap, initial)
        rt.submit(encoder)
        submissions = 1
        state = self.read_state()
        # Sweeps go in chunks within each 25-sweep Newton period; each chunk
        # advances the device sweep count or ends the solve.
        limit = 2 * self.budget + 4
        batches = 0
        while state["result"] == 0:
            if batches > limit:
                raise RuntimeError("Metal solver orchestration exceeded its submission bound")
            batches += 1
            if state["flags"] & REQUESTED:
                submissions += self._newton_chunked()
            else:
                encoder = rt.device.create_command_encoder()
                if state["sweeps"] >= self.budget:
                    self._single(
                        encoder, "iteration_begin", budget=self.budget, state=self.state, status=self.status
                    )
                else:
                    done = state["sweeps"] % BATCH
                    count = min(SWEEP_CHUNK, BATCH - done)
                    self.encode_sweeps(encoder, count, request=done + count == BATCH)
                rt.submit(encoder)
                submissions += 1
            state = self.read_state()
        encoder = rt.device.create_command_encoder()
        self.encode_finish(encoder)
        rt.submit(encoder)
        scalars = rt.read(self.scalars, np.float32)
        status = int(rt.read(self.status, np.uint32)[0])
        return SolveResult(
            converged=state["result"] == 1,
            residual=float(scalars[3]),
            sweeps=state["sweeps"],
            acceleration_steps=state["accelerated"],
            anderson_rejected=state["anderson_rejected"],
            newton_attempts=state["newton_attempts"],
            newton_steps=state["newton_steps"],
            newton_rejected=state["newton_rejected"],
            gmres_iterations=state["gmres_iterations"],
            gmres_converged=state["gmres_converged"],
            breakdowns=state["breakdowns"],
            status=status,
            submissions=submissions + 1,
            magnitude=float(scalars[0]),
            norm=float(scalars[1]),
            scale=float(scalars[2]),
        )
