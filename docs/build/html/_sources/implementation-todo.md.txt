# Solver implementation checklist

This checklist distinguishes implemented capabilities with verification evidence
from outstanding implementation and validation work. A checked capability has
passed its stated tests; it is not a guarantee for every application.

## Implemented capabilities

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
        explicit compatibility configuration.
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
- [x] User guides, input example matrix, rendered manual and regression checks.

- [x] Black formatting, CI formatting/tests/docs/package checks, and Codecov reporting.
- [x] `pypi.yml` Trusted Publishing on version-tag pushes and manual dispatch,
      with tag/version consistency checks and no GitHub environment.
- [x] Production PyPI release `1.0.0`, verified with a clean index installation,
      console entry point and bounded APIC simulation.
- [x] Hydrostatic equilibrium, translated self-equilibrated stress, and nonlinear
      manufactured dynamics with analytical references and regression checks.
- [x] Separate dynamic grid/timestep studies, particle-quadrature studies, frozen
      exact-state force checks, fixed-cloud remaps, and prescribed affine patches.
- [x] Present-tense user documentation with measured verification limits.
- [x] Optional finite-domain integration for linear/quadratic bases, including
      F-convected rectangles, exact knot-split integration, periodic seams and
      explicit rejection of inadmissible wall geometry/support.
- [x] Shared integration weights/gradients for force, strain, APIC, diagnostics and
      segregation; conservation, power pairing and geometry/inlet/outlet checks.
- [x] Controlled point/domain cell-phase and grid-refinement comparisons, plus a
      runnable `domain_column` example and integration-mode diagnostics.

## Acceptance evidence

The regression suite includes `tests/test_modern_solver.py` and
`tests/test_broader_verification.py`. Source-fingerprinted reports in
`benchmarks/verification.json`, `benchmarks/broader_verification.json` and
`benchmarks/quadrature_verification.json` cover
complete example runs, standing waves, equilibrium, finite strain, and operator
controls. The [verification guide](verification.md) explains measured errors and
non-monotone levels. Distribution builds and strict HTML documentation checks
cover packaging and user instructions. These checks verify implementation;
experimental validation and application calibration remain separate tasks.

## Deferred implementation backlog

- [ ] Establish coverage thresholds after reviewing the first Codecov baseline.

- [ ] Control remaining gaps/overlaps of independently convected domains, and
      verify convergence for broader nonuniform deformation and material interfaces.
- [ ] Well-balanced hydrostatics with quadratic interpolation on general particle clouds.
- [ ] Higher-order position integration consistent with deformation/stress evolution;
      isolate timestep accuracy from the number of remaps.
- [ ] Broader application benchmarks and experimental validation.
- [ ] Assess shared-corner/topology tracking or alternative quadrature against the
      remaining geometric integration defects.
- [ ] Profile and accelerate rotated/sheared polygon integration while retaining
      the exact reference path and its conservation checks.

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
