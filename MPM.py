#!/usr/bin/env python3
"""Array MPM timestep and command-line entry point."""

import sys
from pathlib import Path
from rich.progress import BarColumn, Progress, TextColumn, TimeElapsedColumn, TimeRemainingColumn
import numpy as np
import checkpoint
import initialise


def step(P, G, L, *, _limit_cache=None):
    """Advance one complete timestep, without output or changing the clock."""
    if getattr(P, "backend", "cpu") == "metal":
        from mpm_gpu.config import validate
        from mpm_gpu.runner import step as metal_step

        validate(P)
        if _limit_cache is not None:
            _limit_cache.clear()
        return metal_step(P, G, L)
    if P.solver == "granular":
        from granular import step as granular_step

        if _limit_cache is not None:
            _limit_cache.clear()
        return granular_step(P, G, L)
    if not np.isfinite(P.dt) or P.dt <= 0:
        raise ValueError("dt must be finite and positive")
    from stability import limits

    if P.mechanics == "updated":
        P.timestep_limits = limits(P, G, L, cache=_limit_cache)
        safe = min(P.timestep_limits.values(), default=P.dt)
        if P.timestep_policy == "enforce" and P.dt > safe * (1 + 1e-12):
            raise ValueError(
                f"dt={P.dt:g} exceeds estimated stable dt={safe:g}; use dynamic stepping or reduce dt"
            )
    if _limit_cache is not None:
        # No cached estimates may survive a state update or output callback.
        _limit_cache.clear()
    if L.diagnostics is not None:
        L.diagnostics.begin(P, L)
    G.wipe(P)
    P.update_forces()
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    if P.O.check_positions:
        L.recover_position(P, G)
    L.get_nodal_mass_momentum(P, G)
    L.update_stress_strain(P, G)
    L.get_nodal_forces(P, G)
    G.BCs(P)
    G.update_momentum(P)
    G.calculate_gammadot(P, G, P.smooth_gamma_dot)
    if P.segregate:
        G.update_pk(P, G)
        G.calculate_phi_increment(P, L)
        L.move_grainsize_on_grid(P, G)
    L.move_material_points(P, G)
    if L.diagnostics is not None:
        from diagnostics import inventory

        boundary_before = inventory(P, L)
    L.process_boundaries(P, G)
    if L.diagnostics is not None:
        L.diagnostics.finish(P, G, L, boundary_before)


def run(P, G, L, *, max_steps=None, output=True, progress=False):
    """Run prepared state. max_steps permits bounded smoke tests and previews."""
    engine = None
    if getattr(P, "backend", "cpu") == "metal":
        # Resident state: arrays are synchronised only at output times and on return.
        from mpm_gpu.config import validate
        from mpm_gpu.runner import check_callbacks, runner_for

        validate(P)
        if output:
            check_callbacks(P)
        engine = runner_for(P, G, L)
        engine.upload(P, G, L)
    with Progress(
        TextColumn("{task.description}"),
        BarColumn(),
        TextColumn("{task.percentage:>5.1f}%"),
        TextColumn("t = {task.fields[time]:.6g} s"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        disable=not progress,
    ) as display:
        # Rich treats a zero total as indeterminate; an empty run is complete.
        task = display.add_task("Simulating", total=P.t_f or 1, completed=P.t if P.t_f else 1, time=P.t)
        if output:
            P.O.at_start(L, P, G)
        count = 0
        # Kept on P so that a continued or restarted run keeps the save schedule.
        if getattr(P, "next_save", None) is None or P.next_save <= P.t:
            P.next_save = P.t + P.savetime
        limit_cache = {}
        while P.t < P.t_f and (max_steps is None or count < max_steps):
            if progress:
                display.update(task, description="Simulating")
            nominal_dt = P.dt
            remaining = P.t_f - P.t
            if remaining <= 16 * np.finfo(float).eps * max(abs(P.t_f), nominal_dt):
                P.t = P.t_f
                break
            P.dt = min(P.dt, remaining)
            if P.time_stepping == "dynamic" and P.mechanics == "updated":
                from stability import limits

                restrictions = engine.limits(P) if engine else limits(P, G, L, cache=limit_cache)
                if restrictions:
                    P.dt = min(P.dt, min(restrictions.values()))
            used_dt = P.dt
            if engine:
                engine.step(P, G, L)
            else:
                step(P, G, L, _limit_cache=limit_cache)
            P.t += used_dt
            P.tstep += 1
            count += 1
            if progress:
                display.update(task, completed=P.t, time=P.t)
            if output and not engine:
                # check_callbacks admits only the built-in no-op for Metal.
                P.O.after_every_timestep(L, P, G)
            if P.t >= P.next_save:
                # The schedule advances without output too, so it is restart state.
                P.next_save += (np.floor((P.t - P.next_save) / P.savetime) + 1) * P.savetime
                if output:
                    if progress:
                        display.update(task, description="Saving output")
                    if engine:
                        engine.snapshot(P, G, L)
                    P.O.after_every_savetime(L, P, G)
                    if P.O.save_checkpoint:
                        checkpoint.save(Path(P.save_dir) / "data" / "checkpoint.pkl", P, G, L)
            P.dt = nominal_dt
            if P.time_stepping == "dynamic" and P.t < P.t_f:
                if engine:
                    engine.update_timestep(P)
                elif getattr(P.update_timestep, "_uses_limit_cache", False):
                    P.update_timestep(P, G, _limit_cache=limit_cache)
                else:
                    P.update_timestep(P, G)
        if engine:
            engine.snapshot(P, G, L)
        if L.diagnostics is not None:
            L.diagnostics.complete()
        if output:
            if progress:
                display.update(task, description="Finalizing", completed=P.t if P.t_f else 1, time=P.t)
            P.O.at_end(L, P, G)
            if L.diagnostics is not None and P.O.save_diagnostics:
                L.diagnostics.save(P.save_dir)
            if P.O.save_checkpoint:
                checkpoint.save(Path(P.save_dir) / "data" / "checkpoint.pkl", P, G, L)
        if progress:
            display.update(
                task,
                description="Complete" if P.t >= P.t_f else "Stopped",
                completed=P.t if P.t_f else 1,
                time=P.t,
            )
    return P, G, L


def main(params=None):
    try:
        state = initialise.get_parameters(sys.argv[1:] if params is None else params)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    run(*state, progress=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
