"""Legacy KT grain-size integrator with explicit spatial/grain axes.

ax=0 denotes x; ax=1 denotes y. The KT reconstruction and flux equations
are unchanged. Padding is applied to complete spatial rows, never flat buffers.
"""
import numpy as np
from transfers import divide


def div0(a, b):
    return divide(a, b)


def minmod(a, b):
    return .5*(np.sign(a)+np.sign(b))*np.minimum(np.abs(a), np.abs(b))


def BC(c, v):
    v[1] = -divide(c[2]*v[2], c[1])
    v[0] = -divide(c[1]*v[1], c[0])
    v[-2] = -divide(c[-3]*v[-3], c[-2])
    v[-1] = -divide(c[-2]*v[-2], c[-1])
    return v


def _oriented(values, P, ax):
    result = values.reshape((P.G.ny, P.G.nx)+values.shape[1:])
    if ax == 0:
        result = result.swapaxes(0, 1)
        if P.B.cyclic_lr:
            result = result[:-1]  # final grid endpoint duplicates the first
    return result


def _pad(values, periodic):
    padding = [(2, 2)]+[(0, 0)]*(values.ndim-1)
    return np.pad(values, padding, mode='wrap' if periodic else 'edge')


def pad_KT(u, P, ax):
    return _pad(_oriented(u, P, ax), ax == 0 and P.B.cyclic_lr), (P.G.dx if ax == 0 else P.G.dy)


def KT2_flux(phi, G, P, ax):
    periodic = ax == 0 and P.B.cyclic_lr
    inverse_mean = np.sum(phi / P.G.s, axis=-1, keepdims=True)
    fc = divide(1., inverse_mean*P.G.s)-1
    gradient = _pad(_oriented(G.grad_pk[:, ax], P, ax), periodic)[..., None]
    # Legacy variable diffusion: 10*d^3*|gamma|/sqrt(I). At I=0,
    # the old integrator discarded nonfinite flux; retain its zero-flux limit.
    diffusivity = divide(10*G.s_bar**3*np.abs(G.shear_rate), np.sqrt(np.maximum(G.physical('I'), 0)))
    diffusion = _pad(_oriented(diffusivity, P, ax), periodic)[..., None]
    spacing = P.G.dx if ax == 0 else P.G.dy
    derivative = np.gradient(phi, spacing, axis=0)
    boundary = _pad(_oriented(G.boundary_tot | (G.m <= P.M_tol), P, ax), periodic)
    flux = P.c*fc*gradient-diffusion*derivative
    flux[boundary] = 0
    return flux


def KT(P, G, ax):
    C, dx = pad_KT(G.phi, P, ax)
    V = KT2_flux(C, G, P, ax)
    if not (ax == 0 and P.B.cyclic_lr):
        V = BC(C, V)
    cx = minmod((C-np.roll(C, 1, axis=0))/dx, (np.roll(C, -1, axis=0)-C)/dx)
    cpr = np.roll(C, 1, axis=0)-np.roll(cx, 1, axis=0)*dx/2
    cpl = C+cx*dx/2
    cmr = C-cx*dx/2
    cml = np.roll(C, -1, axis=0)+np.roll(cx, -1, axis=0)*dx/2
    vp = (V+np.roll(V, 1, axis=0))/2
    vm = (V+np.roll(V, -1, axis=0))/2
    ap = np.maximum(np.abs(V), np.abs(np.roll(V, 1, axis=0)))
    am = np.maximum(np.abs(V), np.abs(np.roll(V, -1, axis=0)))
    gc = P.dt*(vp*(cpr+cpl)-vm*(cmr+cml)+ap*(cpr-cpl)-am*(cmr-cml))/(2*dx)
    result = gc[2:-2]
    if ax == 0:
        if P.B.cyclic_lr:
            result = np.concatenate((result, result[:1]), axis=0)
        result = result.swapaxes(0, 1)
    return result.reshape(-1, P.G.ns)


def increment_grainsize(P, G):
    return KT(P, G, 0)+KT(P, G, 1)
