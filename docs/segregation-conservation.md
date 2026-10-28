# Segregation transport

The segregation operator transports grain fractions while preserving each
mobile species mass and each particle's total mass. Rigid materials do not
participate. Enable it with `segregate=True` and supply grain diameters,
composition, and the model coefficients described in [the tutorial](tutorial.rst).

For species `k`, the conserved inventory is

```
Q[k] = sum_p m[p] phi[p,k] = sum_i M[i] phi_grid[i,k]
```

Periodic endpoint duplicates count once. Because nodal masses can differ,
conserving an unweighted sum of concentration increments is insufficient. Both
weighted species conservation and local mixture closure `sum_k phi[k]=1` are
required.

## Grid flux and mass exchanges

The relative concentration flux uses segregation and diffusion
coefficients:

    F[k] = c*phi[k]*(1/(s[k]*sum_j(phi[j]/s[j])) - 1)*grad(pk)
           - D*grad(phi[k])

A scalar limiter is shared across the whole reconstructed mixture. Both sides of
a face use one common numerical speed, bounded by the advective Jacobian row
norm. Physical diffusion uses the cell-centre concentration difference separately
from advection and numerical stabilisation. These choices preserve zero total
relative mixture flux without selecting a special dependent grain class.

Each face exchanges a species-mass vector `T[k]`, rather than evolving an
unweighted concentration. Flux density uses the harmonic mean of adjacent mobile
mass densities, `M[i]/V[i]`, and the face's actual area. Boundary control volumes
have half widths. Closed boundaries block exterior faces; empty nodes and internal
walls block transport through their normal faces. Periodic faces are evaluated
once, then copied to the display endpoint.

The same `T[k]` is removed from one node and added to its neighbour. Thus each
species is globally conservative. Concentration increments are obtained by
dividing the resulting nodal mass increments by the participating nodal mass.
Rigid particles are excluded from both the transport inventory and the return
transfer.

## Bounded return to particles

The return transfer uses actual species inventories so depleted classes cannot
supply grains they do not own, and diffusion can mix pure layers.

For a face carrying class k from node i, its donor weights are

    donor[p,i,k] = m[p]*N[p,i]*phi[p,k] / Q[i,k]

where `Q[i,k]` is the deposited mobile species mass. A donor particle loses only
grains it owns and receives equal mass from the counter-flowing classes on that
face. Each particle's total mass stays fixed, and the opposite side receives the
same species masses. For uniform composition this reduces to an incremental concentration gather. When physical and numerical fluxes are zero, particle
composition is unchanged: there is no unconditional PIC smoothing.

Outgoing demand is budgeted over all faces in both directions before updating.
One factor per face limits every class together, so limiting preserves both
species conservation and zero mixture flux. A small floating-point balance of
oppositely signed face totals precedes the return transfer. Fractions are never
clipped or renormalised by this algorithm.

The grid API is `G.calculate_phi_increment(P, L)` for the coupled update. It
stages both `G.dphi` and inventory-aware particle increments. The main timestep
and benchmark driver pass `L`. Calling it without `L`, or calling
`increment_grainsize(P, G)`, returns the standalone grid update; its caller is
responsible for its explicit timestep stability.

## Verification

Checks include nonuniform nodal/particle masses, local fraction sums, coupled
species inventories, rigid particles, free surfaces, closed and periodic faces,
class-order permutation, both spatial directions, pure-layer mixing, and
transport steps far beyond the unbounded explicit update. A node-aligned
sinusoidal mixture reproduces the discrete heat-equation step for constant
diffusivity. The slope examples also run for 1,000 steps with normalisation off.

| Grain classes | Fraction-sum error, 100 steps | Fraction-sum error, 1,000 steps | Max species mass error, 1,000 steps (kg) |
|---:|---:|---:|---:|
| 2 | 6.66e-16 | 3.11e-15 | 1.14e-13 |
| 3 | 4.44e-16 | 2.11e-15 | 5.68e-14 |
| 8 | 4.44e-16 | 1.55e-15 | 2.84e-14 |

Measurements and environment information are in
`benchmarks/segregation_drift.json`. Reproduce them with `python -m
benchmarks.segregation_drift --output /tmp/poly-mpm-segregation.json`.
The [benchmark guide](https://github.com/benjym/poly-mpm/blob/main/benchmarks/README.md) distinguishes report configurations.

These checks establish conservation, bounds, and specific diffusion behaviour;
they do not validate the segregation constitutive law against experiments. The
diffusion coefficient has a special zero-inertial-number convention. Limiting can reduce transport when particle inventories cannot supply
a proposed exchange, so quantitative predictions still require mesh/timestep
convergence checks.

The central-upwind framework and its reconstruction/flux distinction are
described by [Kurganov and Lin (2007)](https://math.umd.edu/~tadmor/centpack/publications/files/Kurganov%20Lin%20CiCP2006-centpack.pdf).
