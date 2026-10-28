"""Independent FP64 CPU reference states for the mechanics hardware gate.

The states mirror benchmarks/granular_prototype.configuration, which is not
packaged. Particles are perturbed off the lattice so every cell is crossed.
"""

import numpy as np


def column(
    spacing,
    particles_per_cell,
    *,
    periodic=False,
    perturbation=0.3,
    seed=7,
    edges=False,
    transfer="flip",
    pic_fraction=0.0,
    motion=False,
):
    """Prepared granular column with the CPU compact operator assembled in FP64."""
    from granular import assemble
    from initialise import prepare
    from inputs.granular_column import Params

    P = Params([])
    P.G.nx = round((P.G.x_M - P.G.x_m) / spacing) + 1
    P.G.ny = round((P.G.y_M - P.G.y_m) / spacing) + 1
    s = spacing / particles_per_cell
    # Edge cases fill the domain so that every node near an edge carries mass.
    width, height = (P.G.x_M - P.G.x_m, P.G.y_M - P.G.y_m) if edges else (1.0, 1.0)
    x, y = np.meshgrid(P.G.x_m + np.arange(s / 2, width, s), P.G.y_m + np.arange(s / 2, height, s))
    P.S.X, P.S.Y = x.ravel(), y.ravel()
    P.S.n, P.S.A = x.size, s**2
    P.dt = 0.005
    P.granular_solver = "velocity"
    P.granular_max_nodes = 10**6  # CPU prototype guard; the fixtures exceed it
    P.transfer, P.pic_fraction = transfer, pic_fraction
    if periodic:
        P.B.has_left = P.B.has_right = False
        P.B.cyclic_lr = True
    P, G, L = prepare(P)
    b = L.materials[0]
    b.x[:, :2] += np.random.default_rng(seed).uniform(-perturbation, perturbation, (len(b), 2)) * s
    if edges:
        # Exercise the closed upper edge and, without periodicity, the right edge.
        b.x[np.argmax(b.x[:, 1]), 1] = P.G.y_M
        if not periodic:
            b.x[np.argmax(b.x[:, 0]), 0] = P.G.x_M
    if motion:
        # Particle velocities and a retained planar velocity gradient (apic_gradient).
        rng = np.random.default_rng(seed + 1)
        b.v[:, :2] = rng.normal(scale=0.1, size=(len(b), 2))
        b.L = np.zeros((len(b), 3, 3))
        b.L[:, :2, :2] = rng.normal(scale=0.5, size=(len(b), 2, 2))
    G.wipe(P)
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    B, free, active, occupied, gap = assemble(P, G, b, compact=True)
    return dict(P=P, G=G, batch=b, operator=B, active=active, free=free, occupied=occupied, gap=gap)
