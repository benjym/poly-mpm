"""Conservative central-upwind segregation on the MPM nodal control volumes.

ax=0 denotes x and ax=1 denotes y. Face fluxes transport species mass;
concentration increments are obtained by dividing by participating nodal mass.
"""
import numpy as np
from transfers import divide


def div0(a, b):
    return divide(a, b)


def minmod(a, b):
    return .5*(np.sign(a)+np.sign(b))*np.minimum(np.abs(a), np.abs(b))


def _mixture(G):
    # A standalone grid test may not have assembled particle participation yet.
    if getattr(G, 'segregation_m', None) is None:
        return G.m, G.phi
    return G.segregation_m, G.segregation_phi


def _oriented(values, P, ax):
    result = values.reshape((P.G.ny, P.G.nx)+values.shape[1:])
    if ax == 0:
        result = result.swapaxes(0, 1)
        if P.B.cyclic_lr:
            result = result[:-1]  # final endpoint duplicates the first
    return result


def _pad(values, periodic):
    padding = [(2, 2)]+[(0, 0)]*(values.ndim-1)
    return np.pad(values, padding, mode='wrap' if periodic else 'edge')


def pad_KT(u, P, ax):
    return _pad(_oriented(u, P, ax), ax == 0 and P.B.cyclic_lr), (P.G.dx if ax == 0 else P.G.dy)


def _advection(phi, forcing, sizes):
    """Physical relative flux and a common bound on its Jacobian's row norm."""
    harmonic_inverse = np.sum(phi/sizes, axis=-1, keepdims=True)
    velocity = forcing*(divide(1., harmonic_inverse*sizes)-1)
    flux = phi*velocity
    rank = divide(forcing*phi, sizes*harmonic_inverse**2)
    inverse_size = 1/sizes
    row_norm = np.abs(velocity-rank*inverse_size)
    row_norm += np.abs(rank)*(inverse_size.sum()-inverse_size)
    return flux, np.max(row_norm, axis=-1, keepdims=True)


def _diffusivity(G, P, phi):
    size = phi @ P.G.s
    # Keep the existing constitutive diffusion coefficient and its I=0 limit.
    return divide(10*size**3*np.abs(G.shear_rate), np.sqrt(np.maximum(G.physical('I'), 0)))


def KT2_flux(phi, G, P, ax):
    """Physical concentration flux, for inspection on a padded spatial array.

    Diffusion is a flux itself, not a velocity to be multiplied by phi again.
    The timestep uses reconstructed face states via _face_transfers below.
    """
    periodic = ax == 0 and P.B.cyclic_lr
    mass, composition = _mixture(G)
    forcing = _pad(_oriented(P.c*G.grad_pk[:, ax], P, ax), periodic)[..., None]
    diffusion = _pad(_oriented(_diffusivity(G, P, composition), P, ax), periodic)[..., None]
    spacing = P.G.dx if ax == 0 else P.G.dy
    advective, _ = _advection(phi, forcing, P.G.s)
    flux = advective-diffusion*np.gradient(phi, spacing, axis=0)
    active = _pad(_oriented(mass > P.M_tol, P, ax), periodic)
    flux[~active] = 0
    return flux


def _slopes(C, active, dx):
    """One limiter for the entire mixture preserves sum_k d(phi_k)/dx = 0."""
    left = (C-np.roll(C, 1, axis=0))/dx
    right = (np.roll(C, -1, axis=0)-C)/dx
    centre = .5*(left+right)
    bound = np.where(left*right > 0, np.minimum(np.abs(left), np.abs(right)), 0.)
    ratio = np.divide(bound, np.abs(centre), out=np.ones_like(centre), where=centre != 0)
    theta = np.minimum(1., np.min(ratio, axis=-1, keepdims=True))
    stencil = active & np.roll(active, 1, axis=0) & np.roll(active, -1, axis=0)
    return centre*theta*stencil[..., None]


def _face_transfers(P, G, ax):
    periodic = ax == 0 and P.B.cyclic_lr
    mass, composition = _mixture(G)
    C, dx = pad_KT(composition, P, ax)
    nodal_mass = _oriented(mass, P, ax)
    active = _pad(nodal_mass > P.M_tol, periodic)
    slopes = _slopes(C, active, dx)
    # N+1 interfaces surround N unique control volumes, including closed edges
    # or the duplicated periodic seam flux.
    left = C[1:-2]+.5*dx*slopes[1:-2]
    right = C[2:-1]-.5*dx*slopes[2:-1]
    forcing = _pad(_oriented(P.c*G.grad_pk[:, ax], P, ax), periodic)
    face_forcing = .5*(forcing[1:-2]+forcing[2:-1])[..., None]
    fl, al = _advection(left, face_forcing, P.G.s)
    fr, ar = _advection(right, face_forcing, P.G.s)
    # A scalar speed shared by all species makes the numerical diffusion sum
    # vanish, since both reconstructed mixtures sum to one.
    flux = .5*(fl+fr)-.5*np.maximum(al, ar)*(right-left)
    diffusion = _pad(_oriented(_diffusivity(G, P, composition), P, ax), periodic)
    face_diffusion = .5*(diffusion[1:-2]+diffusion[2:-1])[..., None]
    flux -= face_diffusion*(C[2:-1]-C[1:-2])/dx
    density = _pad(_oriented(divide(mass, G.V), P, ax), periodic)
    rho_l, rho_r = density[1:-2], density[2:-1]
    face_density = divide(2*rho_l*rho_r, rho_l+rho_r)
    flux *= face_density[..., None]
    # Occupied outer nodes own half control volumes: close only the exterior
    # face, not the face between the boundary node and the interior.
    wall = _oriented(G.boundary_v if ax == 0 else G.boundary_h, P, ax).copy()
    if not periodic:
        wall[0] = wall[-1] = False
    wall = _pad(wall, periodic)
    blocked = (~active[1:-2] | ~active[2:-1] | wall[1:-2] | wall[2:-1])
    flux[blocked] = 0
    if periodic:
        flux[-1] = flux[0]
    else:
        flux[0] = flux[-1] = 0
    if ax == 0:
        area = G.dy_plot*P.G.thickness
    else:
        widths = np.full(P.G.nx, G.dx) if P.B.cyclic_lr else G.dx_plot
        area = widths*P.G.thickness
    transfer = P.dt*area[None, :, None]*flux
    indices = _pad(_oriented(np.arange(len(G.m)), P, ax), periodic)
    left_nodes, right_nodes = indices[1:-2], indices[2:-1]
    selection = slice(None, -1) if periodic else slice(1, -1)
    left_nodes = left_nodes[selection]
    right_nodes = right_nodes[selection]
    transfer = transfer[selection]
    if ax == 1 and P.B.cyclic_lr:
        left_nodes, right_nodes = left_nodes[:, :-1], right_nodes[:, :-1]
        transfer = transfer[:, :-1]
    left_nodes, right_nodes = left_nodes.ravel(), right_nodes.ravel()
    transfer = transfer.reshape(-1, P.G.ns)
    distinct = left_nodes != right_nodes
    return left_nodes[distinct], right_nodes[distinct], transfer[distinct]


def _concentration_increment(P, G, left, right, transfer):
    mass, _ = _mixture(G)
    change = np.zeros_like(G.phi)
    for k in range(P.G.ns):
        change[:, k] = (np.bincount(right, weights=transfer[:, k], minlength=len(mass))
                        - np.bincount(left, weights=transfer[:, k], minlength=len(mass)))
    if P.B.cyclic_lr:
        change[P.G.nx-1::P.G.nx] = change[::P.G.nx]
    return divide(change, mass[:, None])


def _particle_exchanges(P, G, left, right, transfer, particles):
    """Conservative, inventory-aware return transfer with one limiter per face.

    Donate each outgoing species in proportion to its actual particle inventory
    on that face's source node. Exchange equal incoming mass on the same particle,
    so its total mass is unchanged. No PIC smoothing or fraction clipping occurs.
    """
    from scipy.sparse import coo_array
    batches = [b for b in particles if b.material.law != 'rigid' and len(b)]
    if not batches or not len(transfer):
        return np.zeros_like(transfer), [np.zeros_like(b.phi) for b in particles]
    mass, composition = _mixture(G)
    inventory = mass[:, None]*composition
    phi = np.concatenate([b.phi for b in batches])
    nodes = np.concatenate([b.n for b in batches]).copy()
    weights = np.concatenate([b.N for b in batches]).copy()
    if P.B.cyclic_lr:
        nodes[nodes % P.G.nx == P.G.nx-1] -= P.G.nx-1
    weights[mass[nodes] <= P.M_tol] = 0
    rows = np.repeat(np.arange(len(phi)), 4)
    interpolation = coo_array((weights.ravel(), (rows, nodes.ravel())),
                              shape=(len(phi), len(mass))).tocsr()

    def outgoing(positive, negative):
        result = np.zeros_like(inventory)
        for k in range(P.G.ns):
            result[:, k] = (np.bincount(left, weights=positive[:, k], minlength=len(mass))
                            + np.bincount(right, weights=negative[:, k], minlength=len(mass)))
        return result

    positive, negative = np.maximum(transfer, 0), np.maximum(-transfer, 0)
    demand = phi*(interpolation @ divide(outgoing(positive, negative), inventory))
    ratios = np.divide(np.maximum(phi, 0), demand, out=np.ones_like(phi), where=demand > 0)
    node_ratio = np.ones_like(inventory)
    for k in range(P.G.ns):
        present = (weights.ravel() > 0) & (phi[rows, k] > 0)
        np.minimum.at(node_ratio[:, k], nodes.ravel()[present], ratios[rows[present], k])
    node_ratio[inventory <= 0] = 0
    theta = np.minimum(np.min(np.where(positive > 0, node_ratio[left], 1), axis=1),
                       np.min(np.where(negative > 0, node_ratio[right], 1), axis=1))
    theta[theta < 1] *= 1-32*np.finfo(float).eps
    transfer = transfer*theta[:, None]
    positive, negative = np.maximum(transfer, 0), np.maximum(-transfer, 0)
    # Physical and numerical mixture fluxes sum to zero analytically. Balance
    # their positive/negative totals at roundoff before assigning an exchange.
    total_positive = positive.sum(axis=1, keepdims=True)
    total_negative = negative.sum(axis=1, keepdims=True)
    exchange = np.minimum(total_positive, total_negative)
    positive *= divide(exchange, total_positive)
    negative *= divide(exchange, total_negative)
    transfer = positive-negative
    out = phi*(interpolation @ divide(outgoing(positive, negative), inventory))
    incoming = np.zeros_like(phi)
    donor_left = divide(positive, inventory[left])
    donor_right = divide(negative, inventory[right])
    share_left = divide(negative, negative.sum(axis=1, keepdims=True))
    share_right = divide(positive, positive.sum(axis=1, keepdims=True))
    # Class loops keep memory O(particles*classes), not O(particles*classes^2).
    for target in range(P.G.ns):
        coefficients = np.zeros_like(inventory)
        for donor in range(P.G.ns):
            coefficients[:, donor] = (
                np.bincount(left, weights=donor_left[:, donor]*share_left[:, target], minlength=len(mass))
                + np.bincount(right, weights=donor_right[:, donor]*share_right[:, target], minlength=len(mass)))
        incoming[:, target] = np.sum(phi*(interpolation @ coefficients), axis=1)
    delta = incoming-out
    increments = []
    start = 0
    for batch in particles:
        if batch.material.law == 'rigid' or not len(batch):
            increments.append(np.zeros_like(batch.phi))
        else:
            increments.append(delta[start:start+len(batch)])
            start += len(batch)
    return transfer, increments


def KT(P, G, ax):
    return _concentration_increment(P, G, *_face_transfers(P, G, ax))


def segregation_update(P, G, particles=None):
    """Return the conservative grid increment and optional particle exchanges."""
    x = _face_transfers(P, G, 0)
    y = _face_transfers(P, G, 1)
    left, right, transfer = (np.concatenate((a, b), axis=0) for a, b in zip(x, y))
    increments = None
    if particles is not None:
        transfer, increments = _particle_exchanges(P, G, left, right, transfer, particles)
    return _concentration_increment(P, G, left, right, transfer), increments


def increment_grainsize(P, G):
    return segregation_update(P, G)[0]
