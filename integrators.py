#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

A module containing functions for solving partial differential equations on a grid. So far only used for segregation equations.

"""

from numpy import *

def div0( a, b ):
    """This function replaces nans with zeros when dividing by zero.

    :param a: array - Numerator
    :param b: array - Demoniator
    :returns: array - a/b, with infs and nans replaced by 0

    """
    with errstate(divide='ignore', invalid='ignore'):
        c = true_divide( a, b )
        c[ ~ isfinite( c )] = 0  # -inf inf NaN
    return c

def minmod(a,b):
    """Minmod flux limiter.

    :param a: (array)
    :param b: (array)
    :returns: Minmod limited function of a and b

    """
    return 0.5*(sign(a) + sign(b)) * minimum(abs(a),abs(b)) # minmod

def maxmod(a,b):
    """Maxmod flux limiter.

    :param a: (array)
    :param b: (array)
    :returns: Maxmod limited function of a and b

    """
    return 0.5*(sign(a) + sign(b)) * maximum(abs(a),abs(b))

def superbee(a,b):
    """Superbee flux limiter.

    :param a: (array)
    :param b: (array)
    :returns: Superbee limited function of a and b

    """
    return maxmod(minmod(a,2.*b),minmod(2.*a,b))

### NT SCHEME ###

# def pad_NT(P,u,ax):
#     if ax == 0:
#         dx = P.G.dx
#         u = u.reshape(P.G.ny,P.G.nx,P.G.ns).transpose([1,0,2]).reshape(-1,P.G.ns)
#         U = vstack([u[:P.G.ny], u[:P.G.ny], u, u[-P.G.ny:], u[-P.G.ny:]])
#         u = u.reshape(P.G.nx,P.G.ny,P.G.ns)
#         U = U.reshape(P.G.nx+4,P.G.ny,P.G.ns)
#     elif ax == 1:
#         dx = P.G.dy
#         U = vstack([u[:P.G.nx], u[:P.G.nx], u, u[-P.G.nx:], u[-P.G.nx:]])
#         u = u.reshape(P.G.ny,P.G.nx,P.G.ns)
#         U = U.reshape(P.G.ny+4,P.G.nx,P.G.ns)
#     return u,U,dx
#
# def flux(phi,G,P,ax):
#     S_bar = zeros_like(phi)
#     S = zeros_like(phi)
#     S_1_bar = zeros_like(phi)
#
#     boundary = G.boundary_tot.astype(bool) # true boundaries
#     boundary[G.m < P.M_tol] = True # empty cells
#     boundary[:P.G.nx] = True # bottom
#     boundary[(P.G.ny-1)*P.G.nx:] = True # top
#
#     for i in range(P.G.ns): S[:,:,i] = P.G.s[i]
# #     for i in range(P.G.ns): S_bar[:,:,i] = sum(phi*S,axis=2)
#     for i in range(P.G.ns): S_1_bar[:,:,i] = sum(phi*(1./S),axis=2) # INVERSE OF HARMONIC MEAN!
#
#     if ax == 0:
#         pad = (phi.shape[0] - P.G.nx)//2
#         g = G.grad_pk[:,ax].reshape(P.G.ny,P.G.nx).T.flatten() # P.G.nx*P.G.ny
#         boundary = boundary.reshape(P.G.ny,P.G.nx).T.flatten()
#         for i in range(pad):
#             g = hstack([g[:P.G.ny], g,  g[-P.G.ny:]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
#             boundary = hstack([boundary[:P.G.ny], boundary,  boundary[-P.G.ny:]])
#         g = tile(g,[P.G.ns,1]).T.reshape(P.G.nx+2*pad,P.G.ny,P.G.ns)
#         boundary = tile(boundary,[P.G.ns,1]).T.reshape(P.G.nx+2*pad,P.G.ny,P.G.ns)
#         # flux_diff = P.D*(roll(phi,-1,axis=0) - roll(phi,1,axis=0))/(2.*P.G.dx)
#
#     elif ax == 1:
#         pad = (phi.shape[0] - P.G.ny)//2
#         g = G.grad_pk[:,ax] # P.G.ny,P.G.nx
#         g = tile(g,[P.G.ns,1]).T # P.G.ny*P.G.nx,P.G.ns
#         g = g.reshape(-1,P.G.ns) # P.G.ny*P.G.nx,P.G.ns
#         for i in range(pad):
#             g = vstack([g[:P.G.nx], g,  g[-P.G.nx:]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
#             boundary = hstack([boundary[:P.G.nx], boundary,  boundary[-P.G.nx:]])
#         g = g.reshape(P.G.ny+2*pad,P.G.nx,P.G.ns)
#         boundary = tile(boundary,[P.G.ns,1]).T.reshape(P.G.ny+2*pad,P.G.nx,P.G.ns)
#         # flux_diff = P.D*(roll(phi,-1,axis=0) - roll(phi,1,axis=0))/(2.*P.G.dy)
#
#
#     f_c = 1./(S_1_bar*S) - 1.
#     flux = P.c*f_c*g #- flux_diff
#     flux[boundary] = 0
#     return flux#, boundary
#
# def NT(P,G,ax): # ax=0 is horizontal, ax=1 is vertical
#     u,U,dx = pad_NT(P,G.phi,ax)
#
#     U_prime = limiter(U[2:]-U[1:-1],U[1:-1]-U[:-2]) # each point centred at x+1/2
#     F = flux(U,G,P,ax)
#     F_prime = limiter(F[2:]-F[1:-1],F[1:-1]-F[:-2]) # same here is a bit to the right
#
#     u_staggered = U[1:-1] + P.dt*(F_prime)/(2.*dx) # again a bit to the right
#     F_staggered = flux(u_staggered,G,P,ax)
#
#     u_next_staggered = (0.5*(U[1:-2] + U[2:-1]) +
#                         0.125*(U_prime[:-1] - U_prime[1:]) +
#                         P.dt*(F_staggered[1:] - F_staggered[:-1])/dx
#                         )
#
#     u_next = (u_next_staggered[1:] + u_next_staggered[:-1])/2.
#
#     if   ax == 0: return (u_next - u).transpose([1,0,2]).reshape(-1,P.G.ns)
#     elif ax == 1: return (u_next - u).reshape(-1,P.G.ns)

def BC(c,v):
    v[1] = -div0(c[2]*v[2],c[1])
    v[0] = -div0(c[1]*v[1],c[0])
    v[-2] = -div0(c[-3]*v[-3],c[-2])
    v[-1] = -div0(c[-2]*v[-2],c[-1])
    return v


### KT SCHEME ###

def pad_KT(u,P,ax):
    if ax == 0:
        dx = P.G.dx
        u = u.reshape(P.G.ny,P.G.nx,P.G.ns).transpose([1,0,2]).reshape(-1,P.G.ns)
        U = vstack([u[:P.G.ny], u[:P.G.ny], u, u[-P.G.ny:], u[-P.G.ny:]])
        u = u.reshape(P.G.nx,P.G.ny,P.G.ns)
        U = U.reshape(P.G.nx+4,P.G.ny,P.G.ns)
    elif ax == 1:
        dx = P.G.dy
        U = vstack([u[:P.G.nx], u[:P.G.nx], u, u[-P.G.nx:], u[-P.G.nx:]])
        u = u.reshape(P.G.ny,P.G.nx,P.G.ns)
        U = U.reshape(P.G.ny+4,P.G.nx,P.G.ns)
    return U,dx

def KT(P,G,ax):
    C,dx = pad_KT(G.phi,P,ax)
#     V = BC(C,flux(C,G,P,ax))
    V = BC(C,KT2_flux(C,G,P,ax))

    # KT, Eq 4.2
    cx = limiter((C - roll(C,1,axis=0))/dx, (roll(C,-1,axis=0) - C)/dx)
    cpr = roll(C,1,axis=0)  - roll(cx,1,axis=0)*dx/2. # c^+_{j+1/2}
    cpl = C          + cx*dx/2. # c^+_{j-1/2}
    cmr = C          - cx*dx/2. # c^-_{j+1/2}
    cml = roll(C,-1,axis=0) + roll(cx,-1,axis=0)*dx/2. # c^-_{j-1/2}
    vp = (V + roll(V, 1,axis=0))/2.
    vm = (V + roll(V,-1,axis=0))/2.
    ap = maximum(abs(V),abs(roll(V, 1,axis=0)))
    am = maximum(abs(V),abs(roll(V,-1,axis=0)))
    gc = P.dt*(vp*(cpr + cpl) - vm*(cmr + cml) + ap*(cpr - cpl) - am*(cmr - cml) )/(2.*dx)

    if   ax == 0: return gc[2:-2].transpose([1,0,2]).reshape(-1,P.G.ns)
    elif ax == 1: return gc[2:-2].reshape(-1,P.G.ns)

# def KT_flux(C,G,P,ax):
#     if ax == 0:
#         dx = P.G.dx
#         V = zeros_like(C)
#     elif ax == 1:
#         dx = P.G.dy
#         V = ones_like(C)
#         V[:,:,1] = -1.
#     dCdx,tt,tt1 = gradient(C,dx)
#     return V*(1-C) - P.D*dCdx

def KT2_flux(phi,G,P,ax):
    S = zeros_like(phi)
    S_1_bar = zeros_like(phi)

    boundary = G.boundary_tot.astype(bool) # true boundaries
    boundary[G.m < P.M_tol] = True # empty cells
    # boundary[:P.G.nx] = True # bottom
    # boundary[(P.G.ny-1)*P.G.nx:] = True # top

    for i in range(P.G.ns): S[:,:,i] = P.G.s[i]
    for i in range(P.G.ns): S_1_bar[:,:,i] = sum(phi*(1./S),axis=2) # INVERSE OF HARMONIC MEAN!

    D = 1e-1 #P.D
    # D = P.l*(G.s_bar**2.)*abs(G.gammadot)/sqrt(G.I/G.m) # from Pierre, D = l*gamma_dot*d^2/sqrt(I), l \approx 10

    if ax == 0: # x direction
        pad = (phi.shape[0] - P.G.nx)//2
        g = G.grad_pk[:,ax].reshape(P.G.ny,P.G.nx).T.flatten()
        # D = D.reshape(P.G.ny,P.G.nx).T.flatten()
        boundary = boundary.reshape(P.G.ny,P.G.nx).T.flatten()
        for i in range(pad):
            if P.B.cyclic_lr:
                g = hstack([g[-P.G.ny:], g, g[:P.G.ny]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
                # D = hstack([D[-P.G.ny:], D, D[:P.G.ny]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
            else:
                g = hstack([g[:P.G.ny],  g, g[-P.G.ny:]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
                # D = hstack([D[:P.G.ny],  D, D[-P.G.ny:]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
        # boundary = hstack([boundary[:P.G.ny], boundary,  boundary[-P.G.ny:]])
        boundary = hstack([boundary[:boundary.shape[0]//2], zeros([P.G.ny*pad*2],dtype=bool), boundary[boundary.shape[0]//2:]]) # keep just the edges as boundaries
        g = tile(g,[P.G.ns,1]).T.reshape(P.G.nx+2*pad,P.G.ny,P.G.ns)
        # D = tile(D,[P.G.ns,1]).T.reshape(P.G.nx+2*pad,P.G.ny,P.G.ns)
        boundary = tile(boundary,[P.G.ns,1]).T.reshape(P.G.nx+2*pad,P.G.ny,P.G.ns)
        dCdx,tt,tt1 = gradient(phi,P.G.dx)
        # print(g.shape,boundary.shape,D.shape)

    elif ax == 1: # y direction
        pad = (phi.shape[0] - P.G.ny)//2
        g = G.grad_pk[:,ax] # P.G.ny*P.G.nx
        # print(P.G.nx,P.G.ny,P.G.nx*P.G.ny,g.shape)
        g = tile(g,[P.G.ns,1]).T # P.G.ny*P.G.nx,P.G.ns
        g = g.reshape(-1,P.G.ns) # P.G.ny*P.G.nx,P.G.ns
        # print(g.shape)
        # D = tile(D,[P.G.ns,1]).T.reshape(-1,P.G.ns)
        for i in range(pad):
            g = vstack([g[:P.G.nx], g,  g[-P.G.nx:]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
            # D = vstack([D[:P.G.nx], D,  D[-P.G.nx:]]) # (P.G.ny+2*pad)*P.G.nx,P.G.ns
        boundary = hstack([boundary[:boundary.shape[0]//2], zeros([P.G.nx*pad*2],dtype=bool), boundary[boundary.shape[0]//2:]]) # keep just the edges as boundaries
        g = g.reshape(P.G.ny+2*pad,P.G.nx,P.G.ns)
        # D = D.reshape(P.G.ny+2*pad,P.G.nx,P.G.ns)
        boundary = tile(boundary,[P.G.ns,1]).T.reshape(P.G.ny+2*pad,P.G.nx,P.G.ns)
        dCdx,tt,tt1 = gradient(phi,P.G.dy)
        # print(g.shape,boundary.shape,D.shape)

    f_c = 1./(S_1_bar*S) - 1. # NOTE: FLIPPED TO MAKE COMPRESSION POSITIVE - JFM PAPER HAS TENSION POSITIVE

    flux = P.c*f_c*g - D*dCdx
    flux[boundary] = 0

    return flux#, boundary

def limiter(a,b):
#    return minmod(a,b)
    return superbee(a,b)

def increment_grainsize(P,G):
    dphi = KT(P,G,0) + KT(P,G,1)
    # dphi = NT(P,G,0) + NT(P,G,1)
    dphi = normalise_phi_increment(P,G,dphi) # THIS IS A HACK!!!!!!!
    return dphi

def apply_fixed_BC(P,G): # NEVER USE THIS! JUST FOR TESTING!
    phi = G.phi.reshape(P.G.ny,P.G.nx,P.G.ns)
    phi[ :4,:, 0] = 1
    phi[ :4,:,-1] = 0
    phi[-5:,:, 0] = 0
    phi[-5:,:,-1] = 1
    return phi.reshape(P.G.ny*P.G.nx,P.G.ns)


def normalise_phi_increment(P,G,dphi): # NEVER USE THIS! JUST FOR TESTING!
    d = dphi.copy()
    for i in range(P.G.ns):
        phi = G.phi[:,i] + dphi[:,i]
        phi_tot = sum(sum(phi))
        phi[phi<0] = 0.0
        phi[phi>1] = 1.0
        phi *= phi_tot/sum(sum(phi)) # end up with same amount
        d[:,i] = phi - G.phi[:,i]
    return d

if __name__ == "__main__":
    import initialise
    from numpy import random, maximum, ones
    from plotting import Plotting
    import matplotlib.pyplot as plt

    plot = Plotting()
    P,G,L = initialise.get_parameters(['bi_seg_test','23','2','51'])
    G.wipe(P)
    L.get_reference_node(P,G) # Find node down and left
    L.get_basis_functions(P,G) # Make basis functions
    L.get_nodal_mass_momentum(P,G) # Initialise from grid state
    P.O.plot_gsd_debug = True
    plot.draw_gsd_grid(L,P,G)

    # G.grad_pk = -100.*ones([P.G.ny*P.G.nx,3])
    # G.grad_pk[:,0] = 0.
    G.gammadot = 1*ones([P.G.ny*P.G.nx])
    G.I_0 = 1*ones([P.G.ny*P.G.nx])
    G.I = 1*ones([P.G.ny*P.G.nx])
    G.m = 1*ones([P.G.ny*P.G.nx])
    G.pressure = G.pk*1000000. # run faster, also set up a gradient field

    G.I = G.I_0*G.s_bar
    G.pk = G.pressure*G.I
    G.grad_pk = G.calculate_gradient(P,G,G.pk,smooth=False)
    G.grad_pk[:] = G.grad_pk[21]

    while P.t <= P.t_f:
        # G.phi = apply_fixed_BC(P,G)
        G.phi += increment_grainsize(P,G)

        # G.phi += KT(P,G,1) # just in vertical direction
        G.s_bar = G.phi[:,0]*P.G.s[0] + G.phi[:,1]*P.G.s[1]

        # G.I = G.I_0*G.s_bar
        # G.pk = G.pressure*G.I
        # G.grad_pk = G.calculate_gradient(P,G,G.pk,smooth=False)

        P.t += P.dt
        P.tstep += 1
        if P.tstep%100 == 0:
            P.grid_save += 1
            plot.draw_gsd_grid(L,P,G)
            print(' t = ' + str(P.t), end='\r')
