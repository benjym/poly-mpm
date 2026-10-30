"""Square drum of bi_square_drum with the experimental rigid-grain granular solver.

The geometry, rotation, grain sizes and grid segregation follow bi_square_drum.
The Pouliquen elastoplastic material is replaced by rigid grains with a
constant friction coefficient mu = tan(20.9 deg), the quasi-static friction
mu_s of bi_square_drum, because the Metal backend runs only the rigid-grain
step. The rate dependence of mu(I) and the elasticity are therefore absent.

    python MPM.py granular_square_drum 11 2 metal
"""

from types import SimpleNamespace as NS
import numpy as np
from inputs.bi_square_drum import Params as SquareDrumParams


class Params(SquareDrumParams):
    def __init__(self, args):
        super().__init__(args)
        self.backend = args[3]
        self.solver = "granular"
        self.interpolation = "linear"
        self.transfer = "flip"
        self.granular_solver = "velocity"
        self.granular_tolerance = 1e-4
        # The solve stalls between 1e-4 and 3e-4 at the default 2,000 sweeps once
        # the bed starts to avalanche (near 32 deg tilt); 20,000 sweeps converge
        # over the first second at ny = 11.
        self.granular_max_iterations = 2000
        self.dt = 1e-3
        self.initial_stress = None
        # One-sided walls: grains next to a wall leave it when gravity turns
        # away. Two-sided walls hold the mass mapped to wall nodes, which leaves
        # a row of grains stuck inside the first cell along each wall.
        self.B.one_sided = True
        self.output_dir = f"im/granular_drum/ny_{self.G.ny}/Fr_{self._Fr}/"
        legacy = self.S[0]
        self.S = NS(
            law="granular",
            n=legacy.n,
            X=np.asarray(legacy.X),
            Y=np.asarray(legacy.Y),
            A=legacy.A,
            rho=legacy.rho,
            rho_s=legacy.rho_s,
            phi_max=legacy.phi_flow,
            mu_s=legacy.mu_s,
            critical_time=self.critical_time,
        )

    def critical_time(self, P):
        """Time for a grain at rest to fall one grid spacing under gravity.

        Rigid grains have no elastic wave speed, so without this bound the
        dynamic controller (CFL * this time, alongside the advection limit) has
        no limit while the bed is at rest.
        """
        return float(np.sqrt(min(P.G.dx, P.G.dy) / abs(P.gravity)))
