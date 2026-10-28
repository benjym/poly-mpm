"""Array constitutive kernels. Tension is positive; equations follow the legacy solver.

Experimental scalar implementations are retained in legacy/constit_scalar.py.
Kernels update particle history and dstress; the caller applies the Lie correction.
"""
import numpy as np

EYE = np.eye(3)
SUPPORTED = ('rigid', 'elastic', 'von_mises', 'viscous', 'pouliquen')


def invariants(a):
    trace = np.trace(a, axis1=-2, axis2=-1)
    return trace, a - trace[:, None, None] * EYE / 3


def rigid(mp, material, params):
    mp.dstress.fill(0)


def elastic(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    mp.dstress = material.K * trace[:, None, None] * EYE + 2 * material.G * dev
    mp.gammadot = np.linalg.norm(dev, axis=(1, 2)) / params.dt
    mp.pressure = np.trace(mp.stress, axis1=1, axis2=2) / 3


def viscous(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    mp.de_kk = trace / 3
    mp.de_ij = dev
    mp.dp = material.K * mp.de_kk
    mp.pressure += mp.dp
    mp.dev_stress = 2 * material.mu_s * dev / params.dt
    volumetric = material.mu_v * mp.de_kk[:, None, None] * EYE / params.dt
    mp.dstress = mp.pressure[:, None, None] * EYE + mp.dev_stress - volumetric - mp.stress
    mp.gammadot = np.linalg.norm(dev, axis=(1, 2)) / params.dt


def von_mises(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    work = np.abs(np.sum(mp.dev_stress * dev, axis=(1, 2)))
    norm = np.linalg.norm(mp.dev_stress, axis=(1, 2))
    # For s < 2 the legacy expression is singular at zero stress. Its product
    # with the zero work/stress is defined here by its zero-increment limit.
    factor = np.zeros_like(norm)
    active = norm > 0
    factor[active] = (work[active] / (2 * material.k**2)
                      * (norm[active] / (np.sqrt(2) * material.k))**(material.s - 2))
    mp.dev_dstress = 2 * material.G * (dev - factor[:, None, None] * mp.dev_stress)
    mp.sigma_kk += 3 * material.K * trace
    mp.dev_stress += mp.dev_dstress
    mp.dstress = mp.dev_stress + mp.sigma_kk[:, None, None] * EYE / 3 - mp.stress
    mp.yieldfunction = norm / (np.sqrt(2) * material.k) - 1
    mp.gammadot = np.linalg.norm(dev, axis=(1, 2)) / params.dt
    mp.pressure = np.trace(mp.stress, axis1=1, axis2=2) / 3


def pouliquen(mp, material, params):
    trace, dev = invariants(mp.dstrain)
    mp.de_kk = trace / 3
    mp.de_ij = dev
    mp.gammadot = np.linalg.norm(2 * dev / params.dt, axis=(1, 2))
    size = mp.phi @ params.G.s
    pressure = np.abs(mp.pressure)
    mp.I = np.zeros_like(pressure)
    positive = pressure > 0
    mp.I[positive] = mp.gammadot[positive] * size[positive] * np.sqrt(material.rho_s / pressure[positive])
    # I is undefined at zero pressure: store zero for that diagnostic and
    # evaluate friction using its infinite-I limit without a huge sentinel.
    ratio = np.divide(mp.I, material.I_0 + mp.I, out=np.zeros_like(mp.I), where=mp.I > 0)
    ratio[(~positive) & (mp.gammadot > 0)] = 1
    mp.mu = material.mu_0 + material.delta_mu * ratio
    numerator = 2 * np.sqrt(2) * mp.mu * pressure
    mp.eta = np.divide(numerator, mp.gammadot, out=np.zeros_like(pressure), where=mp.gammadot > 0)
    mp.eta[(mp.gammadot == 0) & (numerator > 0)] = np.finfo(float).max
    mp.eta_limited = np.minimum(mp.eta, material.eta_max)
    mp.dev_stress = mp.eta_limited[:, None, None] * dev / params.dt
    mp.dp = material.K * mp.de_kk
    mp.pressure += mp.dp
    cutoff = (mp.pressure > 0) | (mp.rho < mp.initial_rho * .8)
    mp.dp[cutoff] -= mp.pressure[cutoff]
    mp.pressure[cutoff] = 0
    mp.dstress = mp.pressure[:, None, None] * EYE + mp.dev_stress - mp.stress


REGISTRY = {name: globals()[name] for name in SUPPORTED}


def get_law(name):
    try:
        return REGISTRY[name]
    except KeyError:
        raise ValueError(f"Unsupported material law {name!r}. Supported: {', '.join(SUPPORTED)}. "
                         "Experimental laws have not been converted or validated.") from None
