"""Column collapse with critical-state CIDR, the default granular solver.

The geometry of granular_collapse (a 0.5 m by 1 m column released against a wall)
with glass beads: mu(I) of Jop et al., 3 mm grains. The column (phi = 0.6) sits
on the cap at rest at twice the pressure under its own weight, so it is denser
than the cap for its load. See docs/granular-cidr-design.md.
"""

import numpy as np

from inputs.granular_collapse import Params as CollapseParams


class Params(CollapseParams):
    def __init__(self, args):
        super().__init__(args)
        del self.granular_solver  # the default, "cidr"
        S = self.S
        S.packing_max = 0.605
        S.mu = np.tan(np.radians(20.9))
        S.mu_d = np.tan(np.radians(32.76))
        S.I0 = 0.279
        S.grain_diameter = 0.003
        S.phi_min = 0.05
        S.cidr_deficit = 0.005
        S.cidr_p_ref = 2 * S.rho * abs(self.max_g) * 1.0
