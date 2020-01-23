#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

The main module. This contains the loop which will be iterated over every timestep.

"""

import sys
import matplotlib
matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
from plotting import Plotting
import initialise
import MPM
np.seterr(all='ignore')

def main(params):
    """This is the main loop which is repeated every timestep. Currently, this follows the Update Strain Last paradigm (does it?!?).

    :param mode: The name of the input file to use.
    :type mode: str
    :returns: int -- The return code.

    """

    P,G,L = initialise.get_parameters(params)
    plot = Plotting()
    P.O.at_start()

    while P.t <= P.t_f:# Time march
        G.wipe(P) # Discard previous grid
        P.update_forces() # update time dependent gravity
        L.get_reference_nodes(P,G) # Find nearby nodes for each material point
        L.get_basis_functions(P,G) # Make basis functions
        L.get_nodal_mass_momentum(P,G) # Initialise from grid state
        G.get_valid_nodes(P) # get nodes with mass near them
        if P.B.cyclic_lr: G.make_cyclic(P,G,['m','q'])
        L.update_stress_strain(P,G) # Update stress and strain
        L.get_nodal_forces(P,G) # Compute internal and external forces
        G.BCs(P) # Add external forces from BCs
        G.update_momentum(P) # Compute rate of momentum and update nodes
        G.calculate_gammadot(P,G)
        if P.segregate:
            G.update_pk(P,G)
            # if P.B.cyclic_lr: G.make_cyclic(P,G,['phi','pk','s_bar'])
            G.calculate_phi_increment(P)
            L.move_grainsize_on_grid(P,G)
            G.make_cyclic(P,G,['eta'])
        L.move_material_points(P,G) # Update material points (position and velocity)
        if P.B.cyclic_lr: L.cyclic_lr(P,G)

        print('{0:.4f}'.format(P.t*100./P.t_f) + '% Complete, t = ' +
              '{0:.4f}'.format(P.t) + ', g = ' + str(P.g), end='\r')

        if (P.t%P.savetime < P.dt) and (P.t != P.dt): # ignore first timestep
            P.O.after_every_savetime(L,P,G)
        P.O.after_ever_timestep(L,P,G)

        P.t += P.dt # Increment time
        P.tstep += 1
        if P.time_stepping == 'dynamic': P.update_timestep(P,G)
    P.O.at_end(L,P,G)
    print('')
    return 0

if __name__ == '__main__':
    params = sys.argv[1:]
    status = main(params)
    sys.exit(status)
