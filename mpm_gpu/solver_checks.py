"""Solver hardware gate: the resident nonlinear solve against FP64 CPU references.

Converged device solutions are checked with an independent FP64 frozen natural
residual built from the CPU operator, not with the device's own estimate.
Failure paths (budget, nonfinite input, unsupported compression), rank-deficient
Anderson histories, GMRES breakdown and candidate rejection are exercised
explicitly.
"""

import numpy as np

from . import mechanics
from .fixtures import column
from .mechanics_checks import _cpu_grid_step, _inverse_mass, _require, _upload
from .solver import ANDERSON_RCOND, FATAL, HISTORY, Solver, SolverFailure

TOLERANCE = 1e-4
FROZEN_RESIDUAL = 2e-4  # gpu.md science gate
VELOCITY_TOLERANCE = 1e-3  # relative to max |u|, against the CPU at 1e-9
ANDERSON_TOLERANCE = 1e-4
JACOBIAN_TOLERANCE = 1e-5
NATURAL_TOLERANCE = 2e-5  # dimensionless, as the residual
PENDING = ()


class Case:
    """Device operator, preconditioner and solver inputs with FP64 references."""

    def __init__(
        self,
        rt,
        spacing,
        particles_per_cell,
        periodic=False,
        *,
        motion=False,
        max_iterations=200,
        tolerance=TOLERANCE,
    ):
        state = column(spacing, particles_per_cell, periodic=periodic, motion=motion)
        P, G, b = state["P"], state["G"], state["batch"]
        _, trial = _cpu_grid_step(state)
        self.P, self.B, self.active, self.gap = P, state["operator"], state["active"], state["gap"]
        self.trial, self.inverse_mass = trial, _inverse_mass(state)
        self.mu = b.material.mu  # sqrt(2 / dimension) = 1 in 2D
        self.dt = P.dt
        self.rt = rt
        self.grid, particles = _upload(rt, state)
        T = self.grid.topology
        row = np.full(self.grid.nodes, -1)
        row[T.nodes] = np.arange(self.grid.rows)
        self.rows = row[self.active]
        gap = np.zeros((self.grid.rows, 3))
        gap[self.rows] = self.gap
        self.buffers = dict(
            trial=rt.buffer(trial, dtype=np.float32),
            inverse_mass=rt.buffer(self.inverse_mass, dtype=np.float32),
            gap=rt.buffer(gap, dtype=np.float32),
        )
        self.encoder = rt.device.create_command_encoder()
        self.operator = mechanics.assemble(rt, self.encoder, self.grid, particles)
        self.preconditioner = mechanics.preconditioner(
            rt,
            self.encoder,
            self.grid,
            self.operator,
            self.buffers["inverse_mass"],
            self.dt,
            status=self.operator.status,
        )
        self.solver = Solver(
            rt, self.grid, mu=self.mu, dt=self.dt, tolerance=tolerance, max_iterations=max_iterations
        )
        # FP64 scaling exactly as in solve_velocity.
        q = (self.B @ trial.ravel()).reshape(-1, 3) + self.gap
        self.magnitude = max(float(np.max(-q[:, 0], initial=0)), float(np.max(np.abs(q[:, 1:]), initial=0)))
        from .mechanics_checks import _reference_metric

        _, self.norm, _ = _reference_metric(self.B, self.inverse_mass, self.dt)

    def solve(self, initial=None, **buffers):
        inputs = dict(self.buffers, **buffers)
        encoder, self.encoder = self.encoder, None
        return self.solver.solve(
            self.operator,
            self.preconditioner,
            inputs["trial"],
            inputs["inverse_mass"],
            inputs["gap"],
            initial,
            encoder=encoder,
        )

    def read_rows(self, buffer, width):
        return self.rt.read(buffer, np.float32).reshape(-1, width)[self.rows].astype(np.float64)

    def velocity(self, contact):
        return self.trial + self.dt * self.inverse_mass * self.B.transpose(contact.ravel()).reshape(-1, 2)

    def frozen(self, z, u):
        """FP64 natural residual and generalized Jacobian (solve_velocity's Newton block)."""
        B, mu = self.B, self.mu
        w = ((B @ u.ravel()).reshape(z.shape) + self.gap) / self.magnitude
        argument = z - w
        radius = mu * np.maximum(z[:, 0], 0)
        length = np.linalg.norm(argument[:, 1:], axis=1)
        sliding = length > radius
        unit = np.divide(
            argument[:, 1:], length[:, None], out=np.zeros_like(argument[:, 1:]), where=length[:, None] > 0
        )
        ratio = np.divide(radius, length, out=np.zeros(len(z)), where=length > 0)
        projected = argument.copy()
        projected[:, 0] = np.maximum(argument[:, 0], 0)
        projected[:, 1:] *= np.minimum(1, np.divide(radius, length, out=np.ones(len(z)), where=length > 0))[
            :, None
        ]
        natural = z - projected

        def jacobian(value):
            value = value.reshape(z.shape)
            du = self.dt / self.norm * self.inverse_mass * B.transpose(value.ravel()).reshape(-1, 2)
            dw = (B @ du.ravel()).reshape(z.shape)
            da = value - dw
            out = dw.copy()
            out[:, 0] = np.where(argument[:, 0] > 0, dw[:, 0], value[:, 0])
            tangent = ratio[:, None] * (da[:, 1:] - unit * np.sum(unit * da[:, 1:], axis=1)[:, None])
            tangent += (mu * (z[:, 0] >= 0) * value[:, 0])[:, None] * unit
            out[sliding, 1:] = value[sliding, 1:] - tangent[sliding]
            return out

        branches = np.column_stack((argument[:, 0] > 0, sliding, z[:, 0] >= 0))
        return natural, jacobian, branches

    def device_jacobian(self, frozen):
        """FP64 generalized Jacobian at the device's frozen branch (newton_rows layout)."""
        B, mu = self.B, self.mu
        open_, sliding, ratio = frozen[:, 0] > 0, frozen[:, 1] > 0, frozen[:, 2]
        unit, positive = frozen[:, 3:5], frozen[:, 5]

        def jacobian(value):
            value = value.reshape(-1, 3)
            du = self.dt / self.norm * self.inverse_mass * B.transpose(value.ravel()).reshape(-1, 2)
            dw = (B @ du.ravel()).reshape(value.shape)
            da = value - dw
            out = dw.copy()
            out[:, 0] = np.where(open_, dw[:, 0], value[:, 0])
            tangent = ratio[:, None] * (da[:, 1:] - unit * np.sum(unit * da[:, 1:], axis=1)[:, None])
            tangent += (mu * positive * value[:, 0])[:, None] * unit
            out[sliding, 1:] = value[sliding, 1:] - tangent[sliding]
            return out

        return jacobian

    def residual(self, contact):
        """Independent FP64 frozen natural residual of a device contact stress."""
        z = contact / (self.magnitude / self.norm)
        natural, _, _ = self.frozen(z, self.velocity(contact))
        return float(np.max(np.linalg.norm(natural, axis=1), initial=0))


def _summary(result):
    return dict(
        converged=result.converged,
        device_residual=result.residual,
        sweeps=result.sweeps,
        acceleration_steps=result.acceleration_steps,
        anderson_rejected=result.anderson_rejected,
        newton_attempts=result.newton_attempts,
        newton_steps=result.newton_steps,
        newton_rejected=result.newton_rejected,
        gmres_iterations=result.gmres_iterations,
        gmres_converged=result.gmres_converged,
        breakdowns=result.breakdowns,
        status_bits=result.status,
        submissions=result.submissions,
    )


def solve_check(rt, spacing, particles_per_cell, periodic, *, motion=False, single=None):
    """single forces (True) or excludes (False) the single-workgroup iteration path."""
    from granular_mixed import solve_velocity

    case = Case(rt, spacing, particles_per_cell, periodic, motion=motion)
    if single is not None:
        case.solver.single_iteration = single
    reference, stats = solve_velocity(
        case.B, case.trial.ravel(), case.inverse_mass.ravel(), case.gap, case.mu, case.dt, tolerance=1e-9
    )
    result = case.solve().require()
    _require(result.status & FATAL == 0, f"fatal status {result.reasons} on a converged solve")
    _require(result.residual <= TOLERANCE, f"device residual {result.residual:.3g}")
    contact = case.read_rows(case.solver.contact, 3)
    frozen = case.residual(contact)
    _require(frozen <= FROZEN_RESIDUAL, f"FP64 frozen natural residual {frozen:.3g}")
    u_reference = case.velocity(reference)
    u = case.rt.read(case.solver.u, np.float32).reshape(-1, 2).astype(np.float64)
    speed = np.abs(u_reference).max()
    velocity_error = float(np.abs(u - u_reference).max() / speed)
    _require(velocity_error <= VELOCITY_TOLERANCE, f"velocity differs from CPU by {velocity_error:.3g}")
    scale = np.abs(reference).max()
    contact_error = float(np.abs(contact - reference).max() / scale) if scale > 0 else 0.0
    # A warm start from the device solution needs no further sweeps; the
    # repeated solve reuses every workspace.
    allocated = rt.accounting.allocated_bytes
    case.encoder = rt.device.create_command_encoder()
    warm = case.solve(initial=case.solver.contact).require()
    _require(rt.accounting.allocated_bytes == allocated, "a repeated solve allocated device memory")
    _require(warm.sweeps <= 1, f"warm start took {warm.sweeps} sweeps")
    return dict(
        rows=case.grid.rows,
        active_rows=len(case.active),
        single_iteration=case.solver.single_iteration,
        colors=len(case.solver.colors),
        periodic=periodic,
        motion=motion,
        frozen_residual=frozen,
        velocity_error=velocity_error,
        contact_error=contact_error,
        cpu_sweeps=stats["sweeps"],
        cpu_newton_steps=stats.get("newton_steps", 0),
        warm_start_sweeps=warm.sweeps,
        **_summary(result),
    )


def _set_state(case, flags, error=None):
    state = np.zeros(16, dtype=np.uint32)
    state[0] = flags
    case.rt.write(case.solver.state, state, np.uint32)
    if error is not None:
        scalars = case.rt.read(case.solver.scalars, np.float32)
        scalars[3] = error
        case.rt.write(case.solver.scalars, scalars, np.float32)


def jacobian_check(rt, spacing, particles_per_cell, periodic, *, motion=True):
    """Device frozen branch, right-hand side and J v against the FP64 Newton block."""
    case = Case(rt, spacing, particles_per_cell, periodic, motion=motion)
    case.solve().require()
    solver = case.solver
    z = case.read_rows(solver.z, 3)
    u = rt.read(solver.u, np.float32).reshape(-1, 2).astype(np.float64)
    natural, jacobian, branches = case.frozen(z, u)
    v = np.zeros((case.grid.rows, 3))
    v[case.rows] = np.random.default_rng(5).normal(size=(len(case.rows), 3))
    source = rt.buffer(v, dtype=np.float32)
    _set_state(case, 1 | 8)
    encoder = rt.device.create_command_encoder()
    b = solver.bindings
    solver._dispatch(
        encoder,
        "newton_rows",
        case.grid.rows,
        row_count=case.grid.rows,
        mu=solver.mu,
        state=solver.state,
        scalars=solver.scalars,
        neighbors=case.grid.neighbors,
        gradients=b["gradients"],
        gap=b["gap"],
        z=solver.z,
        u=solver.u,
        frozen=solver.frozen,
        rhs=solver.rhs,
    )
    solver._jv(encoder, 1 | 8, source, 0, solver.w)
    rt.submit(encoder)
    frozen = case.read_rows(solver.frozen, 6)
    device_branches = np.column_stack((frozen[:, 0] > 0, frozen[:, 1] > 0, frozen[:, 5] > 0))
    agree = np.all(device_branches == branches, axis=1)
    _require(agree.mean() >= 0.99, f"frozen branches agree on only {agree.mean():.3f} of rows")
    rhs_error = float(np.abs(case.read_rows(solver.rhs, 3) + natural)[agree].max())
    _require(rhs_error <= NATURAL_TOLERANCE, f"natural residual differs by {rhs_error:.3g}")
    expected = jacobian(v[case.rows])
    actual = case.read_rows(solver.w, 3)
    jv_error = float(np.abs(actual - expected)[agree].max() / np.abs(expected).max())
    _require(jv_error <= JACOBIAN_TOLERANCE, f"J v differs by {jv_error:.3g}")
    return dict(
        rows=len(case.rows),
        sliding_rows=int(branches[:, 1].sum()),
        open_rows=int(branches[:, 0].sum()),
        branch_agreement=float(agree.mean()),
        natural_error=rhs_error,
        jv_error=jv_error,
    )


def anderson_check(rt, nx, ny, count, rank, *, seed=0):
    """Least-squares coefficients from device TSQR and Jacobi SVD against FP64 lstsq.

    The history has `count` states and differences of the given rank, so the
    remaining directions are dependent and must fall below the rank cutoff.
    """
    grid = mechanics.device_grid(rt, nx, ny, 1.0, 1.0)
    solver = Solver(rt, grid, mu=0.5, dt=1.0)
    n = solver.size
    rng = np.random.default_rng(seed)
    k = count - 1
    basis = rng.normal(size=(n, rank)) * np.logspace(0, -2, rank)
    df = basis @ rng.normal(size=(rank, k))
    f = rng.normal(size=n) * 0.1 + basis @ rng.normal(size=rank)
    # Steps F_i with F_{i+1} - F_i = df_i and F_newest = f, placed in ring order.
    head = (seed * 5) % HISTORY
    first = (head + HISTORY - count) % HISTORY
    F = np.zeros((count, n))
    F[-1] = f
    for i in range(k - 1, -1, -1):
        F[i] = F[i + 1] - df[:, i]
    steps = np.zeros((HISTORY, n))
    for i in range(count):
        steps[(first + i) % HISTORY] = F[i]
    F32 = steps.astype(np.float32)
    rt.write(solver.steps, F32, np.float32)
    state = np.zeros(16, dtype=np.uint32)
    state[0], state[3], state[4] = 1 | 2, count, head
    rt.write(solver.state, state, np.uint32)
    rt.write(solver.status, np.zeros(1), np.uint32)
    encoder = rt.device.create_command_encoder()
    source = solver.matrix
    solver._dispatch(
        encoder, "anderson_matrix", n, size=n, state=solver.state, steps=solver.steps, matrix=source
    )
    for rows, groups, output in solver.tsqr:
        solver._dispatch(
            encoder,
            "tsqr_tile",
            256 * groups,
            gate=3,
            rows=rows,
            state=solver.state,
            input=source,
            output=output,
        )
        source = output
    solver._single(
        encoder,
        "anderson_solve",
        rcond=ANDERSON_RCOND,
        factor=source,
        state=solver.state,
        status=solver.status,
        coefficients=solver.coefficients,
    )
    rt.submit(encoder)
    coefficients = rt.read(solver.coefficients, np.float32).astype(np.float64)
    _require(np.all(coefficients[k:] == 0), "coefficients beyond the history length")
    c = coefficients[:k]
    # FP64 reference on the FP32-rounded history the device actually saw.
    F64 = F32.astype(np.float64)[[(first + i) % HISTORY for i in range(count)]]
    D, target = np.diff(F64, axis=0).T, F64[-1]
    expected = np.linalg.lstsq(D, target, rcond=ANDERSON_RCOND)[0]
    # With full column rank, R_df and Q.T f are unique up to row signs.
    r_error = 0.0
    if rank == k:
        R = np.linalg.qr(np.column_stack((D, target)), mode="r")[:k]
        factor = rt.read(source, np.float32).reshape(13, 13)[:k, list(range(k)) + [12]].astype(np.float64)
        r_error = float(np.abs(np.abs(factor) - np.abs(R)).max() / np.abs(R).max())
    fit = np.linalg.norm(D @ c - target)
    best = np.linalg.norm(D @ expected - target)
    fit_error = float((fit - best) / np.linalg.norm(target))
    coefficient_error = float(np.linalg.norm(c - expected) / np.linalg.norm(expected))
    _require(r_error <= ANDERSON_TOLERANCE, f"TSQR factor differs by {r_error:.3g}")
    _require(fit_error <= ANDERSON_TOLERANCE, f"least-squares fit worse by {fit_error:.3g}")
    _require(coefficient_error <= 10 * ANDERSON_TOLERANCE, f"coefficients differ by {coefficient_error:.3g}")
    _require(int(rt.read(solver.status, np.uint32)[0]) == 0, "unexpected Anderson status")
    return dict(
        size=n,
        tsqr_levels=len(solver.tsqr),
        history=count,
        rank=rank,
        r_error=r_error,
        fit_error=fit_error,
        coefficient_error=coefficient_error,
    )


def gmres_check(rt, mode):
    """Krylov solve and line search at a prescribed frozen branch.

    identity: J = I, so GMRES converges in one step and x = -natural.
    singular: J = 0 (no free velocity, every row normal-open/tangent-sticking),
    so the first Hessenberg column vanishes; the breakdown leaves x = 0 and
    every line-search candidate is rejected.
    general: the branch of a converged sliding state with a random right-hand
    side; the linear residual is checked with the FP64 Jacobian.
    """
    case = Case(rt, 0.25, 2, False, motion=True)
    case.solve().require()
    solver, R = case.solver, case.grid.rows
    rng = np.random.default_rng(9)
    rhs = np.zeros((R, 3))
    rhs[case.rows] = rng.normal(scale=1e-3, size=(len(case.rows), 3))
    if mode == "general":
        z = case.read_rows(solver.z, 3)
        u = rt.read(solver.u, np.float32).reshape(-1, 2).astype(np.float64)
        _, jacobian, _ = case.frozen(z, u)
        _set_state(case, 1 | 8)
        encoder = rt.device.create_command_encoder()
        solver._dispatch(
            encoder,
            "newton_rows",
            R,
            row_count=R,
            mu=solver.mu,
            state=solver.state,
            scalars=solver.scalars,
            neighbors=case.grid.neighbors,
            gradients=solver.bindings["gradients"],
            gap=solver.bindings["gap"],
            z=solver.z,
            u=solver.u,
            frozen=solver.frozen,
            rhs=solver.rhs,
        )
        rt.submit(encoder)
    else:
        frozen = np.zeros((R, 6))
        if mode == "identity":
            frozen[:, 1] = 1  # sliding with zero ratio and z_n < 0: out = v
        else:
            frozen[:, 0] = 1  # normal open, tangent sticking: out = B du
            rt.write(solver.bindings["inverse_mass"], np.zeros(2 * case.grid.nodes), np.float32)
        rt.write(solver.frozen, frozen, np.float32)
    rt.write(solver.rhs, rhs, np.float32)
    # Search from the converged state, as if it were still above tolerance, so
    # that an unchanged candidate is rejected rather than ending the solve.
    error = float(rt.read(solver.scalars, np.float32)[3])
    solver.tolerance = 1e-5
    _require(error > solver.tolerance, f"residual {error:.3g} too small for the search test")
    _set_state(case, 1 | 8 | 16, error=error)
    rt.write(solver.status, np.zeros(1), np.uint32)
    encoder = rt.device.create_command_encoder()
    solver.encode_krylov(encoder)
    rt.submit(encoder)
    if mode == "singular":
        # Candidates are evaluated with the true inverse mass.
        rt.write(solver.bindings["inverse_mass"], case.inverse_mass, np.float32)
    encoder = rt.device.create_command_encoder()
    solver.encode_search(encoder)
    solver._single(encoder, "newton_end", state=solver.state)
    rt.submit(encoder)
    state = solver.read_state()
    x = case.read_rows(solver.x, 3)
    b = rhs[case.rows]
    result = dict(mode=mode, gmres_iterations=state["gmres_iterations"], breakdowns=state["breakdowns"])
    if mode == "identity":
        error = float(np.abs(x - b).max() / np.abs(b).max())
        _require(state["gmres_iterations"] == 1 and state["gmres_converged"] == 1, f"identity solve {state}")
        _require(error <= 1e-6, f"identity solution differs by {error:.3g}")
        result["solution_error"] = error
    elif mode == "singular":
        _require(
            state["breakdowns"] >= 1 and state["gmres_iterations"] == 0, f"no breakdown recorded {state}"
        )
        _require(np.all(x == 0), "breakdown produced a correction")
        _require(
            state["newton_rejected"] == 1 and state["newton_steps"] == 0,
            f"unchanged candidate accepted {state}",
        )
    else:
        linear = float(np.linalg.norm(jacobian(x) - b) / np.linalg.norm(b))
        _require(state["gmres_iterations"] >= 1, "no Krylov steps")
        _require(linear <= 1e-3, f"FP64 linear residual {linear:.3g}")
        result.update(linear_residual=linear, gmres_converged=state["gmres_converged"])
    _require(int(rt.read(solver.status, np.uint32)[0]) & 1024 == 0, "nonfinite Krylov value")
    return result


def newton_check(rt, spacing, particles_per_cell, periodic, *, motion=False, tolerance=TOLERANCE, accept):
    """One device Newton polish against FP64 GMRES and FP64 candidate residuals.

    Sweeps run until the device requests Newton. The device correction must
    match FP64 GMRES on the same frozen Jacobian, and its accept/reject
    decision must agree with FP64 residuals of the same candidates.
    """
    from scipy.sparse.linalg import LinearOperator, gmres

    case = Case(rt, spacing, particles_per_cell, periodic, motion=motion, tolerance=tolerance)
    solver = case.solver
    encoder, case.encoder = case.encoder, None
    b = case.buffers
    solver.encode_setup(encoder, case.operator, case.preconditioner, b["trial"], b["inverse_mass"], b["gap"])
    rt.submit(encoder)
    state = solver.read_state()
    while state["result"] == 0 and not state["flags"] & 256 and state["sweeps"] < solver.budget:
        encoder = rt.device.create_command_encoder()
        solver.encode_sweeps(encoder)
        rt.submit(encoder)
        state = solver.read_state()
    _require(state["flags"] & 256, f"no Newton request before the solve ended: {state}")
    scale = case.magnitude / case.norm
    z = case.read_rows(solver.z, 3)
    u = rt.read(solver.u, np.float32).reshape(-1, 2).astype(np.float64)
    natural, jacobian, _ = case.frozen(z, u)
    before = case.residual(z * scale)
    J = LinearOperator((natural.size,) * 2, matvec=lambda v: jacobian(v).ravel(), dtype=float)
    expected, _ = gmres(J, -natural.ravel(), rtol=1e-5, atol=0, restart=64, maxiter=4)
    encoder = rt.device.create_command_encoder()
    solver.encode_newton(encoder)
    rt.submit(encoder)
    after = solver.read_state()
    x = case.read_rows(solver.x, 3)
    correction_error = float(np.abs(x - expected.reshape(x.shape)).max() / np.abs(expected).max())
    linear_reference = float(np.linalg.norm(jacobian(expected) + natural) / np.linalg.norm(natural))
    # GMRES accuracy at the branch the device actually froze.
    frozen, rhs = case.read_rows(solver.frozen, 6), case.read_rows(solver.rhs, 3)
    linear = float(np.linalg.norm(case.device_jacobian(frozen)(x) - rhs) / np.linalg.norm(rhs))
    branches = np.column_stack((frozen[:, 0] > 0, frozen[:, 1] > 0, frozen[:, 5] > 0))
    agreement = float(np.mean(np.all(branches == case.frozen(z, u)[2], axis=1)))
    candidates = [case.residual((z + f * x) * scale) for f in (1.0, 0.5, 0.25, 0.125)]
    accepted = after["newton_steps"] == 1
    _require(after["newton_attempts"] == 1, f"Newton not attempted: {after}")
    # Corrections are compared through their linear residuals: rows close to a
    # branch threshold can be frozen differently in FP32 and FP64.
    _require(linear <= 1e-3, f"FP64 linear residual of the device correction {linear:.3g}")
    if accepted:
        committed = case.residual(case.read_rows(solver.z, 3) * scale)
        _require(
            committed < before, f"accepted candidate raises the FP64 residual {before:.3g} -> {committed:.3g}"
        )
        _require(after["count"] == 0, "Anderson history not cleared after a Newton step")
    else:
        _require(min(candidates) >= 0.99 * before, f"rejected a candidate that FP64 accepts: {candidates}")
    _require(
        accepted == accept, f"Newton {'rejected' if accept else 'accepted'}; FP64 candidates {candidates}"
    )
    return dict(
        sweeps=state["sweeps"],
        fp64_residual_before=before,
        fp64_candidate_residuals=candidates,
        accepted=accepted,
        correction_error=correction_error,
        linear_residual=linear,
        fp64_gmres_linear_residual=linear_reference,
        branch_agreement=agreement,
        gmres_iterations=after["gmres_iterations"],
        gmres_converged=after["gmres_converged"],
    )


def failure_check(rt, mode):
    """Explicit failures: exhausted budget, nonfinite input, compression without support."""
    if mode == "budget":
        case = Case(rt, 0.5, 4, True, max_iterations=2)
        expected = 256
        result = case.solve()
    elif mode == "nonfinite":
        case = Case(rt, 0.25, 2, False)
        trial = case.trial.copy()
        trial[np.flatnonzero(case.inverse_mass[:, 1])[0], 1] = np.nan
        expected = 2
        result = case.solve(trial=rt.buffer(trial, dtype=np.float32))
    else:
        case = Case(rt, 0.25, 2, False)
        zero = rt.buffer(np.zeros(2 * case.grid.nodes), dtype=np.float32)
        # The preconditioner was built with the true inverse mass; rebuild it.
        rt.submit(case.encoder)
        case.encoder = rt.device.create_command_encoder()
        case.preconditioner = mechanics.preconditioner(
            rt, case.encoder, case.grid, case.operator, zero, case.dt
        )
        expected = 128
        result = case.solve(inverse_mass=zero)
    _require(not result.converged, f"{mode}: reported convergence")
    _require(result.status & expected, f"{mode}: status {result.reasons}")
    try:
        result.require()
    except SolverFailure:
        pass
    else:
        raise AssertionError(f"{mode}: no exception from require()")
    return dict(mode=mode, sweeps=result.sweeps, status_bits=result.status, reasons=result.reasons)


def cases(rt):
    checks = [
        (
            f"solve_h{h}_ppc{n}_{'periodic' if periodic else 'walls'}{'_moving' if motion else ''}",
            lambda h=h, n=n, p=periodic, m=motion: solve_check(rt, h, n, p, motion=m),
        )
        for h, n, periodic, motion in (
            (0.25, 2, False, False),
            (0.25, 2, True, False),
            (0.1, 3, False, False),
            (0.1, 3, True, False),
            (0.5, 4, True, False),
            (1.0, 6, True, False),
            (0.1, 3, False, True),
            (0.1, 3, True, True),
            (0.05, 2, False, False),
            (0.05, 3, False, True),
        )
    ]
    # Fixtures up to 1,024 rows use the single-workgroup iteration; repeat some
    # with the multi-dispatch iteration used for larger grids.
    checks += [
        (
            f"solve_h{h}_ppc{n}_{'periodic' if periodic else 'walls'}{'_moving' if motion else ''}_multi_dispatch",
            lambda h=h, n=n, p=periodic, m=motion: solve_check(rt, h, n, p, motion=m, single=False),
        )
        for h, n, periodic, motion in ((0.1, 3, True, True), (0.05, 3, False, True), (0.5, 4, True, False))
    ]
    checks += [
        ("jacobian_walls_moving", lambda: jacobian_check(rt, 0.1, 3, False)),
        ("jacobian_periodic_moving", lambda: jacobian_check(rt, 0.1, 3, True)),
        ("anderson_full_rank", lambda: anderson_check(rt, 20, 20, 13, 12)),
        ("anderson_rank_deficient", lambda: anderson_check(rt, 20, 20, 13, 4, seed=1)),
        ("anderson_short_history", lambda: anderson_check(rt, 7, 5, 4, 3, seed=2)),
        ("anderson_three_levels", lambda: anderson_check(rt, 190, 190, 13, 6, seed=3)),
        ("newton_accepted", lambda: newton_check(rt, 0.05, 2, False, tolerance=1e-5, accept=True)),
        ("newton_rejected", lambda: newton_check(rt, 0.5, 4, True, accept=False)),
        ("gmres_identity", lambda: gmres_check(rt, "identity")),
        ("gmres_singular_breakdown", lambda: gmres_check(rt, "singular")),
        ("gmres_general", lambda: gmres_check(rt, "general")),
        ("failure_budget", lambda: failure_check(rt, "budget")),
        ("failure_nonfinite", lambda: failure_check(rt, "nonfinite")),
        ("failure_unsupported", lambda: failure_check(rt, "unsupported")),
    ]
    return checks
