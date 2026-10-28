"""Updated-Lagrangian kinematics with a rotation-preserving exponential map."""

import numpy as np
from scipy.linalg import expm


def advance(batch, P):
    increment = expm(batch.L * P.dt)
    candidate = increment @ batch.F
    determinant = np.linalg.det(candidate)
    bad = np.flatnonzero((determinant <= 0) | ~np.isfinite(determinant))
    if len(bad):
        raise FloatingPointError(
            f"{batch.material.law}, material {batch.material_id}, timestep {P.tstep}: "
            f"invalid deformation determinant, particles {bad[:10].tolist()}"
        )
    batch.F = candidate
    batch.J = determinant
    batch.V = batch.V0 * batch.J
    batch.rho = batch.m / batch.V
