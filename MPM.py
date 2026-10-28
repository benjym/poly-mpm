#!/usr/bin/env python3
"""Array MPM timestep and command-line entry point."""
import sys
import numpy as np
import initialise


def step(P, G, L):
    """Advance one complete timestep, without output or changing the clock."""
    if not np.isfinite(P.dt) or P.dt <= 0:
        raise ValueError('dt must be finite and positive')
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
    L.process_boundaries(P, G)


def run(P, G, L, *, max_steps=None, output=True, progress=False):
    """Run prepared state. max_steps permits bounded smoke tests and previews."""
    if output:
        P.O.at_start(L, P, G)
    count = 0
    next_save = P.t+P.savetime
    while P.t < P.t_f and (max_steps is None or count < max_steps):
        nominal_dt = P.dt
        P.dt = min(P.dt, P.t_f-P.t)
        used_dt = P.dt
        step(P, G, L)
        P.t += used_dt
        P.tstep += 1
        count += 1
        if output:
            P.O.after_every_timestep(L, P, G)
            if P.t >= next_save:
                P.O.after_every_savetime(L, P, G)
                next_save += (np.floor((P.t-next_save)/P.savetime)+1)*P.savetime
        if progress:
            print(f'{P.t/P.t_f:.1%} complete, t = {P.t:.6g}', end='\r')
        P.dt = nominal_dt
        if P.time_stepping == 'dynamic':
            P.update_timestep(P, G)
    if output:
        P.O.at_end(L, P, G)
    return P, G, L


def main(params=None):
    state = initialise.get_parameters(sys.argv[1:] if params is None else params)
    run(*state, progress=True)
    print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
