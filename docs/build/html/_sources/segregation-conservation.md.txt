# Source and correction of segregation drift

## Where the drift came from

There were three discretisation errors, rather than accumulated floating-point
roundoff. Particle-to-grid deposition already conserved species mass.

**The wrong global quantity was conserved.** The previous KT update telescoped
`sum_i delta_phi[i,k]` to zero. The total particle inventory is instead
`Q[k] = sum_p m[p]*phi[p,k]`. The former incremental gather obeys

    delta_Q[k] = sum_i M[i]*delta_phi[i,k]

because `M[i] = sum_p m[p]*N[p,i]`. Therefore an unweighted concentration sum is
insufficient when nodal masses differ. In the slope fixture, boundary nodes carry
13.5 kg and interior nodes 27 kg. On the very first bidisperse timestep,
particle-to-grid mass error is zero, but the old update changes the two species
by approximately +6.78e-7 kg and -6.78e-7 kg. The error accumulates after that.
Periodic duplicate endpoints must also count only once in this identity.

**Diffusion was multiplied by concentration twice.** `KT2_flux` returned
`v[k] - D*grad(phi[k])` and the KT routine treated that entire quantity as a
velocity, multiplying it by the reconstructed `phi[k]`. The resulting diffusion
flux was `-D*phi[k]*grad(phi[k])`, rather than `-D*grad(phi[k])`. Summing the former
over species is generally nonzero even when the fractions sum to one.

**The mixture was discretised as independent scalar equations.** Each class got
a different numerical diffusion speed and an independently limited slope.
Those operations do not preserve `sum_k delta_phi[i,k] = 0`. Physical segregation
velocities do sum to zero when weighted by composition; that identity must also
hold for the numerical face flux, including its stabilisation terms.

The old conservation test used uniform `M[i] = 1` and checked only global sums.
It could not detect the first error and did not check local mixture closure.

## Corrected grid update

The relative concentration flux retains the intended segregation and diffusion
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

The original FLIP-style concentration increment can make a particle's depleted
grain class negative, even if the nodal mixture remains positive. A limiter based
only on that gather was too restrictive: it could freeze diffusion across pure
layers. The corrected return transfer uses the actual species inventories.

For a face carrying class k from node i, its donor weights are

    donor[p,i,k] = m[p]*N[p,i]*phi[p,k] / Q[i,k]

where `Q[i,k]` is the deposited mobile species mass. A donor particle loses only
grains it owns and receives equal mass from the counter-flowing classes on that
face. Each particle's total mass stays fixed, and the opposite side receives the
same species masses. For uniform composition this reduces to the original
incremental gather. When physical and numerical fluxes are zero, particle
composition is unchanged: there is no unconditional PIC smoothing.

Outgoing demand is budgeted over all faces in both directions before updating.
One factor per face limits every class together, so limiting preserves both
species conservation and zero mixture flux. A small floating-point balance of
oppositely signed face totals precedes the return transfer. Fractions are never
clipped or renormalised by this algorithm.

The grid API is now `G.calculate_phi_increment(P, L)` for the coupled update. It
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

| Grain classes | Old fraction-sum error, 100 steps | Corrected, 100 steps | Corrected, 1,000 steps | Max species mass change, 1,000 steps (kg) |
|---:|---:|---:|---:|---:|
| 2 | 0.000157 | 6.66e-16 | 3.11e-15 | 1.14e-13 |
| 3 | 0.000349 | 4.44e-16 | 2.11e-15 | 5.68e-14 |
| 8 | 0.000673 | 4.44e-16 | 1.55e-15 | 2.84e-14 |

Raw measurements are in `benchmarks/segregation_drift_before.json` and
`benchmarks/segregation_drift.json`. Regenerate the corrected report with
`python -m benchmarks.segregation_drift`. The small performance snapshot is
separate from the historical conversion benchmarks.

These checks establish conservation, bounds, and specific diffusion behaviour;
they do not validate the segregation constitutive law against experiments. The
existing diffusion coefficient's zero-inertial-number convention remains
unchanged. Limiting can reduce transport when particle inventories cannot supply
a proposed exchange, so quantitative predictions still require mesh/timestep
convergence checks.

The central-upwind framework and its reconstruction/flux distinction are
described by [Kurganov and Lin (2007)](https://math.umd.edu/~tadmor/centpack/publications/files/Kurganov%20Lin%20CiCP2006-centpack.pdf).
The defects and corrections above were established directly from this code and
the independent conservation/physics checks, rather than attributed to the paper.
