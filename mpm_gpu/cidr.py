"""Critical-state CIDR step (cidr.py, granular_solver="cidr") resident on Metal, in FP32.

The algorithm is the CPU one: the method of multipliers with semismooth Newton
subproblems on the Moreau envelope of psi, whose gradient comes from an exact
proximal step per stress row. Two changes suit the device:

- Each Newton system is solved matrix-free by Jacobi-preconditioned conjugate
  gradients (relative residual `CG_RTOL`) instead of a sparse factorisation.
- The line search uses the sign of the directional derivative, which FP32
  resolves, instead of objective values, whose decrease near convergence it does
  not: along the Newton direction the objective is convex, so a step with a
  nonpositive derivative decreases it. A full step that overshoots slightly (CG
  is inexact) is replaced by the secant step between the derivatives at 0 and 1.

CG scalars stay on the device; the host reads a flag every `CG_CHUNK`
iterations and one or two scalars per Newton iteration.
"""

from dataclasses import dataclass

import numpy as np

from .primitives import Reduction

MODULE = "cidr.slang"
CG_CHUNK = 16
# Grids up to this many nodes run CG in one 1024-thread workgroup, `steps` iterations per
# dispatch; larger ones dispatch per operation.
BLOCK_NODES = 1 << 11
CG_RTOL = 1e-3
CG_LIMIT = 4000
TABLE = 1024
# A decrement this small relative to the kinetic scale is FP32 roundoff of the gradient.
DECREMENT = 1e-12
# A Newton step smaller than this, relative to the trial velocity in the mass norm, ends Newton.
STEP = 1e-6
# Default envelope smoothing in units of dt / rho: iterative linear solves favour a smoother
# envelope (fewer, longer Newton steps) than the CPU's 0.075; the converged step does not
# depend on it.
GAMMA = 1.0
DEBUG = False


@dataclass
class CidrResult:
    converged: bool
    outer: int
    newton: int
    cg: int
    jammed: int
    active: int
    warm: bool
    status: int  # bit 2: a local proximal Newton hit its limit; bit 4: CG breakdown
    submissions: int


class CidrSolver:
    """Resident CIDR solve on the compact operator; `u` and `contact` hold the result."""

    def __init__(self, rt, grid, material, *, dt, gamma, tolerance, max_outer=200, max_newton=200):
        from cidr import CriticalStateCIDR

        self.rt, self.grid = rt, grid
        S = material
        self.dt, self.gamma, self.tolerance = float(dt), float(gamma), float(tolerance)
        self.max_outer, self.max_newton = max_outer, max_newton
        model = CriticalStateCIDR(S.mu_s, S.mu_d, S.I_0, 2 * S.grain_diameter * np.sqrt(S.rho_s))
        lo, hi = np.log(1e-12), np.log(1e3)
        x = np.linspace(lo, hi, TABLE)
        I = np.exp(x)
        alpha = model.alpha(I)
        slope = (2 * model.mu(I) + 2 * I * model.dmu(I) - alpha) / 2.5  # d alpha / d ln I
        self.law = dict(
            mu_s=float(S.mu_s),
            mu_d=float(S.mu_d),
            i0=float(S.I_0),
            k=float(model.k),
            a0=float(model.alpha0),
            table_lo=float(lo),
            table_step=float(x[1] - x[0]),
            table_count=TABLE,
        )
        f32 = np.float32
        self.table = rt.buffer(np.column_stack((alpha, slope)), dtype=f32)
        Delta = S.phi_max - S.phi_min
        self.modulus = dict(
            phi_min=float(S.phi_min),
            phi_max=float(S.phi_max),
            modulus_scale=float(S.cidr_p_ref / (2 * S.mu_s * (Delta / S.cidr_deficit - 1))),
        )
        R, N = grid.rows, 2 * grid.nodes

        def vector(count, dtype=f32):
            return rt.buffer(count=count, dtype=dtype)

        self.bulk, self.phi, self.C, self.lam_abs = vector(R), vector(R), vector(R), vector(R)
        self.lam, self.grad, self.force, self.contact = (
            vector(3 * R),
            vector(3 * R),
            vector(3 * R),
            vector(3 * R),
        )
        self.hessian, self.rows = vector(6 * R), vector(3 * R)
        self.jammed = vector(R, np.uint32)
        self.pressure = vector(R)
        self.u, self.u_trial, self.u_previous = vector(N), vector(N), vector(N)
        self.gradient, self.block = vector(N), vector(3 * grid.nodes)
        self.groups = (grid.nodes + 255) // 256
        self.partials = vector(2 * self.groups)
        self.x, self.r, self.z, self.p, self.q = vector(N), vector(N), vector(N), vector(N), vector(N)
        self.work, self.work2, self.work3, self.row_work = vector(N), vector(N), vector(N), vector(R)
        self.scalars = vector(8)
        self.packed = vector(2)
        self.zero = rt.buffer(np.zeros(N), dtype=f32)
        self.flags = vector(2, np.uint32)
        self.cg_iterations = vector(1, np.uint32)
        self.cg_estimate = CG_CHUNK  # iterations of the previous CG solve: when to first check
        self.fresh = False  # whether the envelope buffers hold the envelope at self.u
        self.status = vector(1, np.uint32)
        self.node_sums = [Reduction(rt, N, operation="sum") for _ in range(3)]
        self.row_sum = Reduction(rt, R, operation="sum")
        self.row_max = Reduction(rt, R, operation="max")

    def encode_finish(self, encoder):
        """The contact is final after solve(); nothing to encode."""

    # ---- encoding helpers --------------------------------------------------------------

    def _dispatch(self, encoder, entry, threads, **bindings):
        self.rt.dispatch(encoder, entry, threads, module=MODULE, **bindings)

    def _sum(self, encoder, index, source):
        return self.node_sums[index].encode(encoder, source, self.status)

    def _envelope(self, encoder, operator, velocity, gamma):
        g = self.grid
        self._dispatch(
            encoder,
            "cidr_envelope",
            g.rows,
            row_count=g.rows,
            gamma=float(gamma),
            **self.law,
            table=self.table,
            neighbors=g.neighbors,
            gradients=operator.gradients,
            velocity=velocity,
            occupied=operator.occupied,
            bulk=self.bulk,
            modulus=self.C,
            lam=self.lam,
            grad=self.grad,
            force=self.force,
            hessian=self.hessian,
            jammed=self.jammed,
            pressure=self.pressure,
            flags=self.flags,
        )

    def _read(self, buffer, dtype=np.float32):
        return self.rt.read(buffer, dtype)

    def _submit(self, encoder):
        self.rt.submit(encoder)
        self.submissions += 1

    # ---- the step --------------------------------------------------------------------

    def solve(self, operator, nodal, initial, volume, encoder, *, cold_gamma, body):
        """Solve the step from the trial velocity; `initial` holds the previous contact.

        cold_gamma: the starting gamma of a cold start (no previous stress), or None to
        decide from the multiplier: a start with any nonzero multiplier is warm.
        """
        rt, g = self.rt, self.grid
        R, N = g.rows, 2 * g.nodes
        self.submissions = 0
        bins = operator.bins
        self._dispatch(
            encoder,
            "cidr_bulk_rows",
            g.rows,
            row_count=g.stress_rows,
            nx=g.nx,
            ny=g.ny,
            periodic=int(g.periodic),
            row_nodes=g.row_nodes,
            cell_offsets=bins.offsets,
            cell_counts=bins.counts,
            permutation=bins.permutation,
            local_x=operator.local_x,
            local_y=operator.local_y,
            volume=volume,
            bulk=self.bulk,
        )
        self._dispatch(
            encoder,
            "cidr_rows_setup",
            R,
            row_count=R,
            **self.modulus,
            occupied=operator.occupied,
            bulk=self.bulk,
            initial=initial,
            phi=self.phi,
            modulus=self.C,
            lam=self.lam,
            lam_abs=self.lam_abs,
        )
        rt.dispatch(encoder, "fill_u32", 2, count=2, value=0, output=self.flags)
        rt.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=self.status)
        self._dispatch(encoder, "cidr_copy", N, count=N, source=nodal.trial, output=self.u)
        largest = self.row_max.encode(encoder, self.lam_abs, self.status)
        self._dispatch(
            encoder,
            "cidr_change",
            N,
            count=N,
            dt=self.dt,
            inverse_mass=nodal.inverse_mass,
            velocity=self.u,
            previous=self.u,
            trial=nodal.trial,
            work=self.work,
            work2=self.work2,
            work3=self.work3,
        )
        scale_sum = self._sum(encoder, 2, self.work3)
        self._submit(encoder)
        warm = float(self._read(largest)[0]) > 0
        scale = float(self._read(scale_sum)[0])
        gamma = self.gamma if warm or cold_gamma is None else max(cold_gamma, self.gamma)
        if warm:
            # Start from the transferred velocity before the body force (granular.cidr_contact).
            encoder = rt.device.create_command_encoder()
            self._dispatch(
                encoder,
                "cidr_unforced",
                g.nodes,
                node_count=g.nodes,
                dt=self.dt,
                body_x=float(body[0]),
                body_y=float(body[1]),
                inverse_mass=nodal.inverse_mass,
                velocity=self.u,
            )
            self._submit(encoder)
        newton_total = cg_total = 0
        converged = False
        outer = 0
        if scale > 0:
            for outer in range(1, self.max_outer + 1):
                encoder = rt.device.create_command_encoder()
                self._dispatch(encoder, "cidr_copy", N, count=N, source=self.u, output=self.u_previous)
                self._submit(encoder)
                its, cgs = self._newton(operator, nodal, gamma, scale)
                newton_total += its
                cg_total += cgs
                # Multiplier update lam = grad psi_gamma(D(u) + gamma lam), and the velocity change.
                encoder = rt.device.create_command_encoder()
                if not self.fresh:
                    self._envelope(encoder, operator, self.u, gamma)
                self.fresh = False  # lam changes, so the envelope at u does too
                self._dispatch(encoder, "cidr_copy", 3 * R, count=3 * R, source=self.grad, output=self.lam)
                self._dispatch(
                    encoder,
                    "cidr_change",
                    N,
                    count=N,
                    dt=self.dt,
                    inverse_mass=nodal.inverse_mass,
                    velocity=self.u,
                    previous=self.u_previous,
                    trial=nodal.trial,
                    work=self.work,
                    work2=self.work2,
                    work3=self.work3,
                )
                sums = [self._sum(encoder, i, w) for i, w in enumerate((self.work, self.work2, self.work3))]
                self._submit(encoder)
                change, size, trial = (float(self._read(s)[0]) for s in sums)
                if gamma <= self.gamma and change <= self.tolerance**2 * max(size, trial):
                    converged = True
                    break
                gamma = max(gamma / np.sqrt(10), self.gamma)
        else:
            converged = True
        encoder = rt.device.create_command_encoder()
        self._dispatch(
            encoder,
            "cidr_contact",
            R,
            row_count=R,
            phi=self.phi,
            lam=self.lam,
            contact=self.contact,
            jammed=self.jammed,
            work=self.row_work,
        )
        jammed = self.row_sum.encode(encoder, self.row_work, self.status)
        self._submit(encoder)
        local = int(self._read(self.flags, np.uint32)[1])
        if int(self._read(self.status, np.uint32)[0]) & 2:
            raise FloatingPointError("Metal CIDR step: nonfinite reduction")
        occupied = self._read(operator.occupied)
        return CidrResult(
            converged=converged,
            outer=outer,
            newton=newton_total,
            cg=cg_total,
            jammed=int(round(float(self._read(jammed)[0]))),
            active=int(np.count_nonzero(occupied[: g.stress_rows] > 0)),
            warm=warm,
            status=local,
            submissions=self.submissions,
        )

    def _newton(self, operator, nodal, gamma, scale):
        """Semismooth Newton on the envelope objective at fixed gamma, from self.u.

        On return the envelope buffers hold the envelope at self.u when self.fresh is set.
        """
        rt, g = self.rt, self.grid
        R, N = g.rows, 2 * g.nodes
        common = dict(
            transpose_rows=g.transpose_rows,
            transpose_slots=g.transpose_slots,
            gradients=operator.gradients,
            inverse_mass=nodal.inverse_mass,
        )
        cg_total = 0
        for iteration in range(1, self.max_newton + 1):
            encoder = rt.device.create_command_encoder()
            if not self.fresh:
                self._envelope(encoder, operator, self.u, gamma)
            self.fresh = True
            self._dispatch(
                encoder,
                "cidr_gradient",
                g.nodes,
                node_count=g.nodes,
                dt=self.dt,
                **common,
                velocity=self.u,
                trial=nodal.trial,
                force=self.force,
                hessian=self.hessian,
                gradient=self.gradient,
                block=self.block,
                partials=self.partials,
            )
            self._dispatch(
                encoder,
                "cidr_cg_init",
                1,
                rtol=CG_RTOL,
                groups=self.groups,
                partials=self.partials,
                scalars=self.scalars,
                flags=self.flags,
                iterations=self.cg_iterations,
            )
            self._dispatch(
                encoder,
                "cidr_cg_start",
                g.nodes,
                node_count=g.nodes,
                gradient=self.gradient,
                block=self.block,
                x=self.x,
                r=self.r,
                z=self.z,
                p=self.p,
            )
            done, encoded = False, 0
            while not done:
                # Encode up to the expected count before the first check, then in chunks.
                chunk = max(CG_CHUNK, self.cg_estimate - encoded) if encoded == 0 else CG_CHUNK
                if g.nodes <= BLOCK_NODES:
                    self._dispatch(
                        encoder,
                        "cidr_cg_block",
                        1024,
                        row_count=R,
                        node_count=g.nodes,
                        dt=self.dt,
                        steps=chunk,
                        limit=CG_LIMIT,
                        neighbors=g.neighbors,
                        **common,
                        hessian=self.hessian,
                        block=self.block,
                        rows=self.rows,
                        x=self.x,
                        r=self.r,
                        z=self.z,
                        p=self.p,
                        q=self.q,
                        scalars=self.scalars,
                        flags=self.flags,
                        iterations=self.cg_iterations,
                    )
                else:
                    for _ in range(chunk):
                        self._cg_iteration(encoder, operator, common)
                encoded += chunk
                self._submit(encoder)
                done = bool(self._read(self.flags, np.uint32)[0])
                encoder = rt.device.create_command_encoder()
            count = int(self._read(self.cg_iterations, np.uint32)[0])
            self.cg_estimate = max(CG_CHUNK, int(0.9 * count))
            cg_total += count
            # Decrement -G . delta and the step's size |delta|^2_(M/dt), read together.
            self._dispatch(encoder, "cidr_products", N, count=N, a=self.gradient, b=self.x, work=self.work)
            slope_sum = self._sum(encoder, 0, self.work)
            self._dispatch(
                encoder,
                "cidr_change",
                N,
                count=N,
                dt=self.dt,
                inverse_mass=nodal.inverse_mass,
                velocity=self.x,
                previous=self.zero,
                trial=self.zero,
                work=self.work2,
                work2=self.work3,
                work3=self.work3,
            )
            size_sum = self._sum(encoder, 1, self.work2)
            self._dispatch(encoder, "cidr_pack", 1, a=slope_sum, b=size_sum, output=self.packed)
            self._submit(encoder)
            slope, size = (float(v) for v in self._read(self.packed))
            # Stop on the step's size, which does not depend on the stiffness 1 / gamma of
            # jammed rows; the decrement test is a backstop at FP32 roundoff.
            if size <= STEP**2 * scale or not -slope > DECREMENT * scale:
                return iteration - 1, cg_total
            step, at_step = self._line_search(operator, nodal, gamma, slope)
            if DEBUG:
                print(f"    newton {iteration}: decrement/scale {-slope / scale:.2e} step {step:.4f}")
            if step == 0.0:
                self.fresh = False
                return iteration - 1, cg_total  # no representable decrease
            encoder = rt.device.create_command_encoder()
            if step != at_step:
                self._dispatch(
                    encoder,
                    "cidr_axpy",
                    N,
                    count=N,
                    t=step,
                    base=self.u,
                    direction=self.x,
                    output=self.u_trial,
                )
            self._dispatch(encoder, "cidr_copy", N, count=N, source=self.u_trial, output=self.u)
            self._submit(encoder)
            # The envelope was last evaluated at u + at_step delta.
            self.fresh = step == at_step
        raise RuntimeError(f"Metal CIDR Newton did not converge at gamma {gamma:.2e}: decrement {-slope:.2e}")

    def _cg_iteration(self, encoder, operator, common):
        g = self.grid
        self._dispatch(
            encoder,
            "cidr_rows_product",
            g.rows,
            row_count=g.rows,
            gate=1,
            neighbors=g.neighbors,
            gradients=operator.gradients,
            source=self.p,
            hessian=self.hessian,
            output=self.rows,
            flags=self.flags,
        )
        self._dispatch(
            encoder,
            "cidr_nodes_product",
            g.nodes,
            node_count=g.nodes,
            dt=self.dt,
            **common,
            source=self.p,
            rows=self.rows,
            output=self.q,
            partials=self.partials,
            flags=self.flags,
        )
        self._dispatch(
            encoder,
            "cidr_cg_alpha",
            1,
            groups=self.groups,
            partials=self.partials,
            scalars=self.scalars,
            flags=self.flags,
        )
        self._dispatch(
            encoder,
            "cidr_cg_update",
            g.nodes,
            node_count=g.nodes,
            block=self.block,
            q=self.q,
            x=self.x,
            r=self.r,
            z=self.z,
            p=self.p,
            scalars=self.scalars,
            partials=self.partials,
            flags=self.flags,
        )
        self._dispatch(
            encoder,
            "cidr_cg_beta",
            1,
            limit=CG_LIMIT,
            groups=self.groups,
            partials=self.partials,
            scalars=self.scalars,
            flags=self.flags,
            iterations=self.cg_iterations,
        )
        self._dispatch(
            encoder,
            "cidr_cg_direction",
            2 * g.nodes,
            count=2 * g.nodes,
            z=self.z,
            p=self.p,
            scalars=self.scalars,
            flags=self.flags,
        )

    def _line_search(self, operator, nodal, gamma, slope):
        """Step along delta with a nonpositive derivative: t = 1, else regula falsi (Illinois).

        phi'(0) = slope < 0; the objective is convex along delta, so phi' increases. Returns
        the step and the last step at which the envelope was evaluated.
        """
        at_one = self._derivative(operator, nodal, gamma, 1.0)
        # A full step that overshoots the minimiser along delta only slightly (CG is inexact)
        # still decreases a near-quadratic objective by about half the predicted decrease.
        if at_one <= -0.1 * slope:
            return 1.0, 1.0
        lo, f_lo, hi, f_hi = 0.0, slope, 1.0, at_one
        last, side = 1.0, 0
        for _ in range(8):
            t = (lo * f_hi - hi * f_lo) / (f_hi - f_lo)
            f = self._derivative(operator, nodal, gamma, t)
            last = t
            if f <= 0:
                if f > 0.1 * slope:  # close enough to the minimiser along delta
                    return t, last
                lo, f_lo = t, f
                if side == -1:
                    f_hi *= 0.5
                side = -1
            else:
                hi, f_hi = t, f
                if side == 1:
                    f_lo *= 0.5
                side = 1
        return lo, last

    def _derivative(self, operator, nodal, gamma, t):
        """d/dt of the envelope objective at u + t delta; leaves u + t delta in u_trial."""
        rt, g = self.rt, self.grid
        R, N = g.rows, 2 * g.nodes
        encoder = rt.device.create_command_encoder()
        self._dispatch(
            encoder, "cidr_axpy", N, count=N, t=float(t), base=self.u, direction=self.x, output=self.u_trial
        )
        self._envelope(encoder, operator, self.u_trial, gamma)
        self._dispatch(
            encoder,
            "cidr_slope_nodes",
            N,
            count=N,
            dt=self.dt,
            inverse_mass=nodal.inverse_mass,
            velocity=self.u_trial,
            trial=nodal.trial,
            direction=self.x,
            work=self.work,
        )
        nodes = self._sum(encoder, 0, self.work)
        self._dispatch(
            encoder,
            "cidr_slope_rows",
            R,
            row_count=R,
            neighbors=g.neighbors,
            gradients=operator.gradients,
            direction=self.x,
            force=self.force,
            work=self.row_work,
        )
        rows = self.row_sum.encode(encoder, self.row_work, self.status)
        self._dispatch(encoder, "cidr_pack", 1, a=nodes, b=rows, output=self.packed)
        self._submit(encoder)
        return float(self._read(self.packed).astype(np.float64).sum())
