# Source audit and citing-paper review

This audit compares implementation choices before changing the constitutive
model. It does not reproduce the published 3D collapse experiment. The original
FLIP/intrinsic-2D/zero-viscosity configuration remains the default.

## Forward literature search

This is a focused review of direct author follow-ups, rather than an exhaustive
citation census. Only public papers and public code were used. The findings
change the order of the next checks, rather than establishing a new model.

**Rousseau, G., Métivet, T., Rousseau, H., Daviet, G., and
Bertails-Descoubes, F. (2023).** *Revisiting the role of friction coefficients
in granular collapses: confrontation of 3-D non-smooth simulations with
experiments*. Journal of Fluid Mechanics, 975, A14.
[doi:10.1017/jfm.2023.835](https://doi.org/10.1017/jfm.2023.835).

This directly validates a descendant of the 2016 method against collapse
experiments. Crucially, Section 2.3 explicitly acknowledges interpolation-induced
constraint violations and volume loss; parameters were adjusted to keep total
relative volume losses below 3%. This is not a bound on maximum local packing.
Section 4.1 uses roughly 216 randomly placed particles per cell and a timestep
of 8.3e-4 s. Appendix A retains APIC and lumped mass. The paper compares free
surfaces and static/flowing boundaries, and notes early-onset discrepancies.
These observations support a matched numerical audit, not the assumption that
the source method has exact density preservation. The paper does not specify
the volume estimator sufficiently here to equate it with our diagnostics.

**Daviet, G. (2026).** *Mixed Material Point Methods for Stiff Elastoplasticity*.
ACM Transactions on Graphics, 45(4), Article 151.
[doi:10.1145/3811345](https://doi.org/10.1145/3811345);
[author paper and code](https://research.nvidia.com/labs/prl/mixed_mpm/).

This is a direct extension of the earlier mixed formulation. Section 4.2 explores
discontinuous stress spaces, predominantly Q1–P0 and Q1–dP1. Section 5 develops
a different local flow-rule solver and GPU iteration. Section 7.1 explicitly
states that its trilinear mixed elements do not satisfy the inf-sup condition,
and reports slow asymptotic convergence of colored Gauss–Seidel/Jacobi.
Advection still has a CFL restriction. Nonassociated plasticity is handled with
a bipotential; it should not be described as an ordinary convex minimization
merely from the abstract. This is a strong candidate for a later solver and
space comparison, rather than evidence that sparse storage alone fixes our
packing error or guarantees a particular speedup.

## Original implementation inspected

The original [Sand6 1.0 archive](https://gdaviet.fr/files/mpm/sand6-1.0.tar.gz)
is available from the [author's code page](https://gdaviet.fr/code/sand6/).
The archive SHA256 and hashes of relevant source files are recorded in
[the frozen fixture](../tests/fixtures/sand6_frozen.json).

The code confirms the gradient-carried affine transfer in
`src/simu/DynParticles.cc`, the small-viscosity trial solve and mass-only contact
increment in `src/simu/PhaseSolver.cc`, and positive packing-gap coefficients.
Particle splitting and merging are compiled in; our comparison does not include
that population control. `scenes/collapse.conf` sets `enforceMaxFrac` to zero,
although the code also contains an optional geometric projection. We therefore
do not add a geometric correction as part of reproducing that configuration.

### Executed frozen-cloud oracle

The wrapper [sand6_oracle.cpp](../benchmarks/sand6_oracle.cpp) calls the original
`FormBuilder::addTauDu` and `addDuDv`, rather than reimplementing their formulas.
It was compiled against the archive's Eigen/bogus headers and original
`FormBuilder.cc`. The only compatibility edit was in bundled Eigen:
`Transpositions.h`, `matrix.derived(), trt.derived()` became
`matrix.derived(), trt`; original assembly code was unchanged.

```sh
c++ -std=c++11 -O1 -w -I "$SAND6/src" -I "$SAND6/vendor/eigen3" \
    -I "$SAND6/vendor" benchmarks/sand6_oracle.cpp \
    "$SAND6/src/simu/FormBuilder.cc" -Wl,-dead_strip -o /tmp/sand6-oracle
```

This command records the macOS build used for the audit; other linkers may use
different dead-code elimination flags. `$SAND6` is the extracted archive folder.
The input is a particle count followed by unit-cell `x y z grain_volume` rows.
The stored fixture uses all eight combinations of `x,y,z` in `{0.25,0.75}`,
with grain volume `0.0375` each. The wrapper prints B (48x24), then the D:D
matrix (24x24). Summing tied in-plane velocities and stress nodes across the
two z layers reduces these to the stored 2D matrices.

The original tensor mapping is orthonormal for `(tau:sigma)/2`, while ours uses
`tau:sigma`. Consequently, source B equals `sqrt(2)` times our plane-strain B;
the stress coefficients scale reciprocally. Our `2*eta*D:D` stiffness also
matches. Regression tests compare the executed-source fixture, lumped masses,
affine integrals, friction normalization, and packing-gap factors. This is one
frozen cloud and a constrained reduction, not a general implementation proof.

## Controlled implementation comparison

```sh
uv run --locked python -m benchmarks.granular_source_audit --output /tmp/granular-source-audit.json
```

The committed [report](../benchmarks/granular_source_audit.json) fingerprints the
implementation and fixture. Cases use the same domain, walls, initial mass and
particle cloud, dt=0.005 s, and final time 0.3 s. A denser-cloud case changes
only particle sampling. Each is one bounded timing, not a performance claim.

| Transfer | Stress | Viscosity | Particles/cell/axis | Peak packing / target | Maximum integrated overpacking |
| --- | --- | ---: | ---: | ---: | ---: |
| FLIP | Intrinsic 2D | 0 | 2 | 1.03537 | 0.818% |
| Source gradient | Intrinsic 2D | 0 | 2 | 1.05747 | 1.363% |
| Source gradient | Intrinsic 2D | 0.001 | 2 | 1.05747 | 1.363% |
| Source gradient | Plane strain | 0.001 | 2 | 1.05705 | 1.380% |
| Source gradient | Intrinsic 2D | 0.001 | 6 | 1.01949 | 0.259% |

Integrated overpacking is explicitly `sum(max(phi-phi_max,0)*nodal_volume)`
divided by total fixed grain volume. It measures how much grain volume exceeds
nodal capacity. It is neither the paper's unspecified global volume-loss
estimator nor the deformation-based `sum(V)` diagnostic. A low value does not
establish agreement with the paper's 3% criterion.

Source-style transfer changes the trajectory but does not remove the coarse-grid
overshoot. The tested small viscosity has negligible impact on this diagnostic.
Increasing particle sampling reduces integrated excess here. The deformation
volume changes are retained separately in JSON; F/J is not the source packing
variable and does not feed the contact operator.

The plane-strain solve initially stalled at t=0.21 s in our natural-map least
squares. An equivalent second-order-cone Fischer–Burmeister merit map with the
nonassociated Coulomb shift resolves this case. Its result must still satisfy
the original natural-map tolerance; solver termination alone is never accepted.
Tests cover sticking, sliding, opening and rejection of associated dilation.
This was a solver limitation in our implementation, not a demonstration that
the underlying granular model lacks a solution.

## Audit options and remaining scope

- `transfer='apic_gradient'` is restricted to this granular branch. It carries
  the previous grid velocity gradient directly, with PIC particle velocities.
  Existing moment-based quadratic APIC is unchanged. Affine energy/angular
  diagnostics include the current bilinear support moment.
- `granular_stress='plane_strain'` retains the 3D trace and deviator with planar
  motion and independent sigma_zz. It requires the hybrid solver. Cone friction
  is `sqrt(2/3)*mu`, keeping physical `tau_eq/p=mu`. Intrinsic 2D remains default.
- `granular_viscosity` defaults to zero. A positive value solves the implicit
  Newtonian trial stage, then uses the mass-only contact increment, matching
  the source split. Viscous work is recorded; splitting errors remain in the
  energy residual.

Before calling this a faithful trajectory reproduction, match source geometry,
sidewall/contact treatment, sampling/randomization and split/merge behavior;
identify its actual total-volume estimator and reproduce an author-code run.
The newer mixed spaces are an informed next experiment after those checks,
rather than a correction justified solely by our current density plots.
