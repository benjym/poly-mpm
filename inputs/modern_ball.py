"""Elastic-ball example using updated mechanics, quadratic splines, and APIC."""

from inputs.ball import Params as BallParams


class Params(BallParams):
    def __init__(self, args):
        super().__init__(args)
        self.mechanics = "updated"
        self.interpolation = "quadratic"
        self.transfer = "apic"
        self.time_stepping = "dynamic"
        self.CFL = 0.15
        # The inherited controller supplies a timestep; central wave/viscous
        # restrictions additionally constrain it before the first and later steps.
