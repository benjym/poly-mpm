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
- [x] Opt-in compact velocity-based rigid-grain iteration, scalar local flow
      solves and cached metrics, with coupled reference and physical checks.
- [x] Experimental Q1–P0 stress space, independent assembly checks and a
      controlled sampling/grid/packing study that retains failed gates.

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

### Implemented from the backlog

- [x] Checkpoint/restart including controller and constitutive state:
      `O.save_checkpoint` and `checkpoint.save`/`load`. A restarted CPU run is
      bitwise identical to an uninterrupted one for elastic, quadratic APIC,
      segregation, finite-domain and rigid-grain runs (`tests/test_checkpoint.py`).
- [x] Profile and accelerate rotated/sheared polygon integration while retaining
      the exact reference path and its conservation checks: batched clipping is
      about 8 times faster and agrees with the per-domain reference to 4e-15.
- [x] Profile and accelerate verified kernels, first pass: the per-step finite
      check and the elastic timestep estimate (closed-form planar singular
      values, LAPACK kept for nonplanar F) cut 10–32% per step in the profiled
      examples. The remaining time is spread over `einsum` products and the
      stress update with no single hotspot.
- [x] Metal: silo walls, inlets and outlets (host boundary events, checked
      against the CPU), and the plotting cost in `gpu.md`.
- [x] Locking mitigation: opt-in F-bar volumetric averaging (`fbar`) with an
      incompressible cantilever benchmark; at nu = 0.499 the tip deflection
      ratio rises from 0.148 to 0.941 of beam theory. A mixed pressure
      formulation was not added.

### Actionable without further input

None remaining; every item below needs a decision first.

### Needs your input (second pass)

Each item below depends on a decision, data or a scope choice that the code
cannot settle. The question to answer is given for each.

- [ ] Establish coverage thresholds after reviewing the first Codecov baseline.
      *Question:* which thresholds (overall and patch) should CI enforce?
- [ ] Resolve Q1–P0 iterative convergence and flow-accuracy limits, or evaluate
      dP1 stresses, before selecting a compact discretization for GPU porting.
      *Question:* tied to the modelling decision for refined connected flows in
      `gpu-handoff.md` (fixed stress spacing, a rheology change or uncertified
      solves); which direction?
- [ ] Repeated connected 100k-particle Metal benchmark (`gpu.md`). Blocked by the
      same decision: connected refined domains do not converge on either backend.
- [ ] Metal `increment` stopping and `granular_viscosity`. Both are options
      within the same decision: increment stopping is the "uncertified solves"
      route and the viscosity a regularisation of the rheology. *Question:* is
      either the chosen direction, so that it should run on Metal?
- [ ] Control remaining gaps/overlaps of independently convected domains, and
      verify convergence for broader nonuniform deformation and material interfaces.
      *Question:* which approach (shared-corner topology tracking, CPDI-style
      shared corners or an alternative quadrature)? This also settles the
      shared-corner assessment item below.
- [ ] Assess shared-corner/topology tracking or alternative quadrature against the
      remaining geometric integration defects. *Question:* as above.
- [ ] Well-balanced hydrostatics with quadratic interpolation on general particle
      clouds. *Question:* which well-balancing approach (hydrostatic
      reconstruction, a reference-state subtraction or quadrature changes), and
      which particle clouds must it cover?
- [ ] Higher-order position integration consistent with deformation/stress
      evolution; isolate timestep accuracy from the number of remaps.
      *Question:* which integrator family (midpoint/RK2 positions with matching
      F updates, or a convected-grid scheme), and is a change to the default
      update acceptable?
- [ ] Mixed pressure formulation beyond F-bar, and F-bar checks for plasticity,
      large-strain rotation and mixed-material interfaces. *Question:* which
      materials and benchmarks matter, and is F-bar enough for them?
- [ ] Curved rigid boundaries, contact reactions, friction and separate body
      fields. Frictional and one-sided walls on the domain edges now exist.
      *Question:* which geometries and body interactions are needed?
- [ ] Implicit/quasistatic solver and constitutive consistent tangents.
      *Question:* is this in scope, and for which materials?
- [ ] Finite-strain multiplicative plasticity and calibrated granular
      constitutive extensions. *Question:* which yield/flow laws, and which data
      to calibrate against?
- [ ] Nonlocal granular rheology / size effects, with independent validation
      data. *Question:* which nonlocal model and which validation data?
- [ ] Conservative particle adaptivity and/or grid refinement.
      *Question:* splitting/merging, adaptive grids or both, and for which flows?
- [ ] Three-dimensional mechanics, transport, boundaries, and verification.
      *Question:* is 3D in scope?
- [ ] Pore-pressure / solid–fluid coupling if required by the intended
      applications. *Question:* do the applications need it?
- [ ] Broader application benchmarks and experimental validation.
      *Question:* which applications and data sets?
- [ ] Experimental validation and application-specific parameter calibration.
      *Question:* as above.
