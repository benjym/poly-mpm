#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

The main module. This contains the loop which will be iterated over every timestep.

"""

import sys
import matplotlib
matplotlib.use('Agg')
# import matplotlib.pyplot as plt
from numpy import zeros, sum, seterr, mean
from plotting import Plotting
import initialise
import MPM
seterr(all='ignore')

def main(params):
    """This is the main loop which is repeated every timestep. Currently, this follows the Update Strain Last paradigm (does it?!?).

    :param mode: The name of the input file to use.
    :type mode: str
    :returns: int -- The return code.

    """

    P,G,L = initialise.get_parameters(params)
    plot = Plotting()
#     if P.O.plot_material_points: plot.draw_material_points(L,P,G,'initial')

    while P.t <= P.t_f:# Time march
        G.wipe(P) # Discard previous grid
        P.update_forces() # update time dependent gravity
        L.update_forces(P,G) # pass body forces to material points
        L.get_reference_node(P,G) # Find node down and left
        L.get_basis_functions(P,G) # Make basis functions
        if P.O.check_positions: L.recover_position(P,G) # check basis functions
        L.get_nodal_mass_momentum(P,G) # Initialise from grid state
        if P.B.cyclic_lr:
#             G.apply_cyclic_BCs(P)
            G.make_cyclic(P,G,['m','q'])
        L.update_stress_strain(P,G) # Update stress and strain
        L.get_nodal_forces(P,G) # Compute internal and external forces
        G.BCs(P) # Add external forces from BCs
        G.update_momentum(P) # Compute rate of momentum and update nodes
        G.calculate_gammadot(P,G)
        if P.segregate_grid:
            G.calculate_grad_gammadot(P,G)
            G.calculate_phi_increment(P)
            L.move_grainsize_on_grid(P,G)

        L.move_material_points(P,G) # Update material points (position and velocity)
        # Move/Remove any particles that have left the grid
        if P.B.outlet_left: L.outlet_left(P,G)
        if P.B.outlet_bottom: L.outlet_bottom(P,G)
        if P.B.inlet_right: L.inlet_right(P,G)
        if P.B.cyclic_lr: L.cyclic_lr(P,G)

        # Output related things
        if P.O.measure_energy: P.O.energy[P.tstep] = L.measure_elastic_energy(P,G) # measure energy
        print('{0:.4f}'.format(P.t*100./P.t_f) + '% Complete, t = ' +
              '{0:.4f}'.format(P.t) + ', g = ' + str(P.g), end='\r')

        if P.t%P.savetime < P.dt:
            if P.O.plot_gsd_mp: plot.draw_gsd_mp(L,P,G)
            if P.O.plot_gsd_grid: plot.draw_gsd_grid(L,P,G)
#             plot.draw_voronoi(P,G)
            if P.O.plot_continuum: plot.draw_continuum(G,P)
            if P.O.plot_material_points: plot.draw_material_points(L,P,G)
            if P.mode == 'anisotropy': plot.draw_gamma_dot(L,P,G)
            if P.O.measure_energy: P.O.measure_E(P,L)
            if P.O.save_u: plot.save_u(L,P,G)
            if P.O.save_s_bar: plot.save_s_bar(L,P,G)
        if P.mode == 'dp_unit_test' or P.mode == 'dp_rate_unit_test': P.O.store_p_q(P,G,L,P.tstep)
        if P.mode == 'pouliquen_unit_test': P.O.store_mu(P,G,L,P.tstep)

        # Increment time
        P.t += P.dt
        P.tstep += 1
#         for p in range(P.phases):
#             if (P.S[p].law is 'von_mises' or P.S[p].law is 'dp') and not P.has_yielded:
#             if P.S[p].law is 'von_mises' and not P.has_yielded:
#                 L.update_timestep(P) # Check for yielding and reduce timestep

    # Final things to do
    if P.O.plot_material_points: plot.draw_material_points(L,P,G,'final')
    if P.O.measure_stiffness: P.O.measure_E(P,L)
    if P.O.measure_energy: plot.draw_energy(P)
    if P.O.plot_gsd_mp: plot.draw_gsd_mp(L,P,G)
    if P.O.plot_gsd_grid: plot.draw_gsd_grid(L,P,G)
    if P.O.save_u: plot.save_u(L,P,G)
    if P.O.save_s_bar: plot.save_s_bar(L,P,G)
    if P.mode == 'dp_unit_test' or P.mode == 'dp_rate_unit_test': P.O.draw_p_q(P,G,L,plot,P.tstep)
    if P.mode == 'pouliquen_unit_test': P.O.draw_mu(P,G,L,plot,P.tstep)
    print('')
    return 0

if __name__ == '__main__':
    params = sys.argv[1:]
    status = main(params)
    sys.exit(status)
