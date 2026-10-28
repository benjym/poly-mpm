# Solver implementation checklist

This file tracks implementation and evidence, rather than feature suggestions.
The current authorised tranche is stages 1–4, in order. Later stages remain
backlog; marking a capability complete requires its acceptance checks to pass.

## Current tranche

- [x] 1. Verification and diagnostics
  - [x] Track mass, species inventory, momentum, angular momentum, kinetic energy,
        material work/dissipation, external work, and constraint impulses.
  - [x] Report support loss, invalid deformation, density/volume ranges, and
        timestep restrictions; do not silently hide invalid state.
  - [x] Analytical/convergence checks and reproducible benchmark report.
- [x] 2. Consistent deformation, volume, density, and stress
  - [x] Store reference volume, deformation gradient, determinant, and histories;
        update volume/density from deformation with constant particle mass.
  - [x] Objective stress integration and physical fluid bulk/shear conventions;
        explicit compatibility mode for historical scalar comparisons.
  - [x] Rotation, affine compression, stress/work, inlet/outlet, and species tests.
- [x] 3. Interchangeable interpolation and quadratic B-splines
  - [x] One stencil interface for weights, gradients, nodes, and physical offsets.
  - [x] Linear and quadratic implementations, including physical boundaries and
        periodic seams; no clipping/renormalisation that breaks completeness.
  - [x] Generalise every transfer, including segregation, to stencil size.
  - [x] Partition, affine reproduction, gradients, cell crossing, and convergence.
- [x] 4. APIC
  - [x] Matched particle/grid affine transfers with explicit affine momentum.
  - [x] PIC and FLIP controls; conservative rank-deficient-support handling.
  - [x] Translation, affine reproduction, transfer angular momentum/energy,
        force-driven evolution, boundary and segregation integration checks.
- [x] Refresh user docs, example audit, and rendered manual; record migration and
      remaining numerical limitations; run the complete regression suite.

- [x] Black formatting, CI formatting/tests/docs/package checks, and Codecov reporting.

## Acceptance evidence

The completed tranche is checked by the full regression suite, the independent
`tests/test_modern_solver.py` checks, and the source-fingerprinted
`benchmarks/verification.json` report. The report includes complete `modern_ball`
and `dp_rate_unit_test` runs, standing-wave refinement, remapping/dissipation
comparisons, and caveats about non-monotone errors. The wheel and strict HTML
manual are also checked. These are implementation checks, not experimental
validation. Rate-form plasticity and local granular rheology retain the limits
described in the updated-methods guide.

## Deferred implementation backlog

- [ ] Confirm PyPI project ownership, choose a release version, test the distributions
      on TestPyPI, and configure Trusted Publishing before the first upload.
- [ ] Establish coverage thresholds after reviewing the first Codecov baseline.

- [ ] Broader hydrostatic/translated-stress/large-deformation and application
      benchmarks, separating remapping, spatial, and temporal errors.
- [ ] Benchmark-driven assessment of GIMP/CPDI or higher-order interpolation,
      if quadratic boundary/interior accuracy is insufficient.

- [ ] Locking mitigation: F-bar / mixed pressure formulation, with incompressible benchmarks.
- [ ] Curved rigid boundaries, contact reactions, friction and separate body fields.
- [ ] Implicit/quasistatic solver and constitutive consistent tangents.
- [ ] Finite-strain multiplicative plasticity and calibrated granular constitutive extensions.
- [ ] Nonlocal granular rheology / size effects, with independent validation data.
- [ ] Checkpoint/restart including controller and constitutive state.
- [ ] Conservative particle adaptivity and/or grid refinement.
- [ ] Three-dimensional mechanics, transport, boundaries, and verification.
- [ ] Pore-pressure / solid–fluid coupling if required by the intended applications.
- [ ] Profile and accelerate verified kernels; retain a numerical reference backend.
- [ ] Experimental validation and application-specific parameter calibration.
