"""Exploratory collapse with the coupled rigid-grain prototype."""

import numpy as np
from inputs.granular_column import Params as ColumnParams


class Params(ColumnParams):
    def __init__(self, args):
        super().__init__(args)
        self.G.nx = 9
        self.G.x_M = 2.0
        self.B.roughness = True
        self.t_f = 0.3
        self.savetime = 0.05
        keep = self.S.X < 0.5
        self.S.X = self.S.X[keep]
        self.S.Y = self.S.Y[keep]
        self.S.n = int(np.count_nonzero(keep))
