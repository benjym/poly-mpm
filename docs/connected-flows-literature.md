# Refined connected flows: literature review

The rigid-grain solver converges on coarse grids and on disconnected domains,
but on refined connected domains every solver tried stalls at natural residuals
of 1e-2 to 1e-3 ([granular-mixed](granular-mixed.md)). Coloured Gauss–Seidel
with Anderson acceleration, exact Newton, the assembled least-squares reference
and an interior-point fixed point on De Saxcé's shift all fail alike. The
frictionless and associated (convex) versions of the same problems solve
readily, a dilatant flow rule lowers the residual without removing the
growth with resolution, and holding the stress spacing at 0.05 m while the
velocity grid refines converges. This review asks whether the difficulty lies in
the solver, the discretisation or the continuum model, and how others handle it.

The model in question is a rigid-plastic granular continuum: rigid grains with a
hard packing limit, constant Coulomb friction `mu`, a non-dilatant
(non-associated) flow rule and no rate dependence, written as a cone
complementarity problem in nodal stresses on bilinear velocity and nodal (Q1)
or cell (P0) stress spaces.

Papers listed under each heading were read in full or in the sections relevant
here. Papers seen only as abstracts or metadata are listed separately at the
end and are not used as evidence.

## Well-posedness of the continuum model

Barker, Schaeffer, Bohorquez and Gray (2015) linearise the incompressible
`mu(I)` equations in 2D and find short-wave growth rates proportional to the
square of the wavenumber in an ill-posed region of parameter space. With
constant friction, the Coulomb case, the equations are ill-posed for every
deformation state; this recovers Schaeffer (1987), which was not read here. For
the `mu(I)` curve of Jop et al. the equations are well-posed only over an
intermediate band of inertial number, about two decades wide, and ill-posed as
`I -> 0` and at large `I`. The authors list column collapse, silos and rotating
drums as problematic because they contain quasi-static zones. In ill-posed
cases, grid refinement resolves faster-growing modes and numerical results do
not converge; coarse grids or numerical viscosity suppress the growth, which the
authors describe as an unsatisfactory regularisation.

Barker and Gray (2017) show that for incompressible flow, well-posedness as
`I -> 0` requires `mu(0) = 0`. A finite static friction `mu_s > 0` is therefore
incompatible with a well-posed incompressible local model. Their partial
regularisation replaces the low-`I` branch by a curve that decays
logarithmically to zero, with one new parameter. Decelerating chute flows with
the Jop et al. curve show localised strain-rate patterns that depend on
resolution; with the regularised curve, results at 256² and 512² agree. The
cost is the loss of a true yield stress: static material creeps slowly. Very
large `I` remains ill-posed.

Barker, Schaeffer, Shearer and Gray (2017) add compressibility: a yield
condition `|tau| = Y(p, phi, I)` and a dilatant flow rule
`div u = 2 f(p, phi, I) |D|`. The system is linearly well-posed in 2D for all
states if `Y` and `f` are linked by a specific relation and satisfy two
monotonicity inequalities. Without rate dependence the relation reduces to the
critical-state flow rule derived from normality, which is well-posed except
possibly along two wave directions. Well-posedness therefore needs a
near-associated coupling between the pressure dependence of the yield surface
and dilatancy, plus rate dependence; an arbitrary dilatancy function does not
qualify. The proof is 2D; 3D is conjectured.

Heyman, Delannay, Tabuteau and Valance (2017) reach a similar result with
`phi(I) = phi_max - Delta phi I`, a rate-dependent equation of state and a bulk
friction `mu_b`. The problem is well-posed if `mu_b > 1 - 7 mu / 6`. In the
incompressible limit the ill-posedness of Barker et al. (2015) returns. With
constant `mu`, compressible Coulomb shear flow is unstable but well-posed,
provided `mu_b` exceeds the bound; a larger `mu_b` places the most unstable mode
at a finite wavenumber, which sets a required grid size.

Goddard and Lee (2017) relate the Hadamard instability to loss of ellipticity of
the inertia-free equations, the classical material instability behind
localisation. A dissipative hyperstress (a gradient term) restores ellipticity
and gives a finite shear-band width, at the cost of a new length scale and
higher-order boundary conditions.

Li and Henann (2019) show that the nonlocal granular fluidity model suppresses
the Hadamard instability in all conditions, leaving instability only at finite
wavelength with bounded growth. In the quasi-static limit the shortest unstable
wavelength depends only on `mu_s`; for `mu_s = 0.3819` and nonlocal amplitude
`A = 0.48` it is about 24.5 grain diameters. Plane-strain compression with the
nonlocal model gives the same shear-band width at mesh sizes 2d, d, d/2 and
d/4, while the local model gives bands that sharpen with the mesh. A tensorial
inertial-number gradient model does not remove the instability, so adding a
gradient does not by itself guarantee well-posedness.

### Consequences for the solver

Our bulk law is the constant-friction, incompressible-when-packed,
rate-independent case that Barker et al. (2015) identify as ill-posed for every
deformation. If the discrete step inherits this, refining the grid admits modes
the step cannot control, and no algorithm should be expected to reach a fixed
tolerance at every resolution. This is an inference: none of these papers
analyses the discrete cone complementarity problem, and the stall has not been
shown to come from ill-posedness. The rest of our evidence fits it:

- The frictionless and associated problems are convex minimisations, and the
  rate-independent limit of the well-posed compressible model of Barker et al.
  (2017) is normality, so these versions sit on the well-posed side.
- A fixed coarser stress spacing acts as a wavenumber cutoff, the numerical
  regularisation that Barker et al. (2015) describe as resolution-dependent.
- A dilatant flow rule `w_n = beta |w_t|` helps only partly. Barker et al. (2017)
  and Heyman et al. (2017) need a specific yield–dilatancy coupling and rate
  dependence or bulk dissipation; dilatancy alone leaves ill-posed directions.

## Solution algorithms

Gräser and Sander (2019) analyse truncated nonsmooth Newton multigrid (TNNMG):
a nonlinear block Gauss–Seidel smoother followed by a Newton correction on the
inactive subspace, solved inexactly with one multigrid V-cycle and a line
search, with no penalty parameters. Global convergence is proved for a smooth
energy plus convex, block-separable nonsmooth terms. Sander (2017) applies it to
primal von Mises plasticity with kinematic hardening on one to six uniform
refinement levels: iterations to reach 1e-9 stayed mostly bounded independent of
the mesh size, became erratic while the plastic front crossed the domain, and
the method was about 40 times faster than predictor–corrector with a direct
solver on the finest grid. Hardening makes the problem strictly convex; perfect
and non-associated plasticity are not treated. Kornhuber et al. (2008) report
monotone multigrid for frictionless two-body contact with mostly 20 to 25
V-cycles per time step and warn that plain multigrid on linearised contact can
fail.

Bleyer (2022) reviews conic programming in nonsmooth mechanics. Interior-point
solves (Mosek) of elastoplastic steps take 10 to 18 iterations for about 250,000
variables and are described as robust to load step and mesh refinement, with one
mesh tabulated. The review states that non-associated behaviour has no
potential, can only be written as a mixed complementarity problem, and is
treated in the literature by sequences of associated conic programs that need
further work and have no convergence guarantee. This matches our velocity-space
study, where each convex subproblem took 23 to 46 interior-point iterations up to
9,963 stress unknowns and the outer non-associated fixed point stalled.

Treskatis, Moyers-González and Price (2016) solve Bingham flow on the stress
dual with an accelerated proximal gradient method (FISTA), showing that the
augmented Lagrangian method ALG2 is an unaccelerated special case. Rates are
`O(1/sqrt(k))` for ALG2 and `O(1/k)` with acceleration; the step size comes from
the plastic viscosity, independent of the mesh. A rigid-plastic law has no
viscosity, so this smoothness is absent and the remaining inertial term grows
like `h^-2`; dual gradient methods would then need iterations that grow with
refinement. That last point is an inference. Ionescu, Mangeney, Bouchut and
Roche (2015) solve granular column collapse with an augmented Lagrangian method
and state that a minimum viscosity of 0.1 Pa s is necessary for it; Coulomb wall
friction is also regularised. They report no mesh study.

Song, Fan, Ascher and Pai (2026) split non-associated Coulomb friction into a
strongly convex cone quadratic program with the De Saxcé coupling frozen and an
outer update of that coupling (forward–backward–forward instead of Picard). The
coupling is Lipschitz but not monotone, a direct consequence of
non-associativity, so convergence is shown only empirically. A proximal term
makes the inner problem uniquely solvable when the Delassus operator is singular,
which addresses static indeterminacy. They report that on globally coupled
contact graphs, such as a masonry arch with force chains through the whole
structure, the explicit correction loses contractivity and load-bearing steps
need many more outer iterations, while weakly coupled problems need few. This is
the behaviour of our De Saxcé fixed point on connected domains.

No paper read shows mesh-independent convergence for a non-associated,
rate-independent Coulomb or Drucker–Prager problem with indeterminate stress.
The methods with that evidence (TNNMG, monotone multigrid, interior point) need
a convex energy or an associated law.

## Practice in granular continuum codes

Daviet and Bertails-Descoubes (2016) introduced the semi-implicit MPM closest to
ours: rigid Drucker–Prager with an exact cone and a packing complementarity, no
elasticity, a small Newtonian viscosity (1 mPa s), and optionally `mu(I)`
evaluated explicitly per particle. Velocity and stress both use trilinear
functions on the same grid nodes, the constraint is imposed at nodes, and a
projected Gauss–Seidel solver is warm-started to a tolerance of 1e-3 times
`rho g L`, which they report is satisfied in most cases, typically in 25 to 250
iterations. Grids reach about 121,000 active nodes; there is no refinement
study. Rousseau, Métivet, Rousseau, Daviet and Bertails-Descoubes (2023) use the
same solver with constant `mu`, unregularised, against experimental column
collapses at up to 250 by 150 cells in 2D. The yield constraint holds only at
nodes; they report small but noticeable violations as volume loss and tuned
parameters to keep it below 3%. They do not report a mesh-convergence study.

In his thesis, Daviet (2016) solves Bingham and Drucker–Prager flows with finite
elements and shows convergence under uniform refinement on Poiseuille flow for
several velocity–stress pairs; P2–P1d (quadratic velocity, discontinuous linear
stress) gives a good ratio of accuracy to cost and P4–P1d the most accurate
results. The thesis notes that existence and uniqueness guarantees are weak and
the stress is not unique for a given velocity, which did not cause problems in
its tests because flowing regions set the boundary conditions for rigid ones. A
1D compatibility argument rules out several pairs for equilibrium under gravity
when the packing fraction lives in the stress space, and the MPM avoids this by
carrying it on particles. Ionescu et al. (2015) use Taylor–Hood velocity and
pressure with a discontinuous linear stress deviator and state that the spaces
cannot be chosen independently.

The other codes read keep finite elastic stiffness. Dunatunga and Kamrin (2015)
use explicit MPM with a hypoelastic-plastic `mu(I)` law, no plastic dilation and
a stress-free disconnected state, at grid spacings of 0.0125 and 0.00625 m, with
no grid-convergence study; they note high-frequency stress noise and smooth
stresses only for plotting. Klár et al. (2016) model sand with Hencky elasticity
and a Drucker–Prager return map, using a Young's modulus about a hundred times
below their stated physical value. Gaume et al. (2018) model snow with Cam-clay
elastoplasticity and scale the softening with the grid spacing so that
dissipated energy is independent of the mesh. Finite stiffness makes stress a
function of elastic strain, which removes the static indeterminacy; none of
these papers shows grid convergence for constant friction.

Henann and Kamrin (2016) implement the nonlocal fluidity model in implicit
finite elements, with fluidity as an extra nodal field, F-bar elements and
consistent tangents giving quadratic Newton convergence. Their results converge
better than linearly and are mesh-independent at element sizes of a third, a
half and one grain diameter, while the local model under identical conditions is
mesh-dependent.

Zhao and Choo (2020) show that equal-order low-order MPM violates the inf–sup
condition in the undrained limit, producing checkerboard pressure, and
stabilise it by penalising the difference between the pressure and its cell
average (a polynomial pressure projection). Too little stabilisation leaves the
checkerboard and too much over-smooths.

## Options

Every published route to grid-independent results gives up at least one of
incompressibility, rate independence or locality, or keeps elasticity. The
current model keeps all four, which the well-posedness results above suggest is
not possible. The options below are ordered from least to most change to the
model.

| Option | Evidence | What it changes | Solver consequence |
| --- | --- | --- | --- |
| Accept the regularisation and stop on a loose criterion | Current practice of Daviet and Bertails-Descoubes (2016), Rousseau et al. (2023) and Daviet (2026): tolerance 1e-3 `rho g L` or a plastic-increment stop, moderate grids, no refinement study | Nothing in the model; results depend on resolution | Works now (`granular_stopping = "increment"`, coarse stress spacing); not convergent |
| Fixed coarser stress spacing | Our study converges at 0.05 m stress spacing up to 204,800 particles; Daviet's thesis favours low-order stress with higher-order velocity | A numerical length for the stress field | Converges; Barker et al. (2015) call this kind of cutoff resolution-dependent |
| Stable discretisation (higher-order velocity, discontinuous stress, or projection stabilisation) | Daviet thesis (P2–P1d), Ionescu et al. (2015), Zhao and Choo (2020) | Discretisation only | Addresses indeterminacy and checkerboarding, not continuum ill-posedness |
| Partial `mu(I)` regularisation | Barker and Gray (2017): grid-independent at 256² and 512² | Rate dependence; no static yield stress (slow creep) | Friction depends on the unknown strain rate; still ill-posed at large `I` |
| Compressible critical-state rheology with near-normal dilatancy and rate dependence | Barker et al. (2017), Heyman et al. (2017): linearly well-posed in 2D | Replaces the hard packing limit by `phi(p, I)`; new parameters (`phi(I)`, `mu_b` or the yield–dilatancy link) | Near-associated coupling, close to the convex problems our solvers already scale on; no published nonlinear grid-convergence study |
| Nonlocal granular fluidity | Li and Henann (2019), Henann and Kamrin (2016): mesh-independent | Extra fluidity field and equation, one parameter, a grain-scale length | The cooperativity length (about 25 grain diameters) is below our grid spacing for millimetre grains, so on our grids it may act as a local model |
| Finite elastic stiffness | Dunatunga and Kamrin (2015), Klár et al. (2016), Gaume et al. (2018) | Drops rigid grains; small timestep or implicit solve | Removes static indeterminacy; no grid-convergence evidence for constant `mu` |

## Recommendation

The literature points at the model rather than the solver. A further search for
a better algorithm for the current rate-independent, non-associated,
incompressible-when-packed law is unlikely to succeed: the well-posedness
results predict non-convergent behaviour under refinement, the closest
published code accepts loose, uncertified solves at moderate resolution, and
the algorithm papers that scale need convexity.

The option with the best fit to this code is the compressible critical-state
rheology with near-normal dilatancy and rate dependence (Barker et al. 2017;
Heyman et al. 2017). It keeps rigid grains and locality, it is proved linearly
well-posed in 2D, and its near-associated structure moves the per-step problem
towards the convex problems on which interior-point solves already showed
nearly mesh-independent iteration counts in our velocity-space study. Its costs
are the hard packing limit, which becomes a pressure- and rate-dependent packing
fraction, two or three new parameters, and the absence of a published nonlinear
grid-convergence study, which this code would have to provide.

Before committing, one experiment can separate model from solver at low cost:
repeat the frozen first-step refinement study (`granular_stress_scaling`) with
the associated flow rule `beta = mu` and a small rate-dependent bulk viscosity
in the volumetric constraint, which together approximate the well-posed
critical-state limit, and solve with the convex interior-point path. If the
iteration counts stay bounded from 0.05 to 0.0125 m where the non-associated
law stalls, the model is the obstacle and the compressible rheology is worth
implementing. If they do not, the discretisation (stress indeterminacy, the
Q1–Q1 pair) needs attention first, along the lines of Daviet's thesis and Zhao
and Choo (2020).

This experiment has been run ([regularisation study](granular-mixed.md#regularised-laws-under-refinement)).
On the frozen first step, the associated law converges in one interior-point
solve on grids from 0.05 to 0.00625 m (25 to 41 iterations) and its velocity
field converges at about first order, while the non-associated law fails beyond
0.05 m with residuals that grow under refinement. A shear viscosity alone does
not restore convergence. The model, not the solver, is the obstacle; the next
step is a critical-state dilatancy law and a test of an evolving flow.

For production runs in the meantime, the fixed coarser stress spacing is the
only option already shown to converge in this code. It should be reported as a
numerical regularisation with a stated length.

## Published models

The review and the regularisation study point to a compressible, critical-state
rheology. Three complete published models of that kind were read.

**Inertial CIDR (Schaeffer et al. 2019).** The yield and dilatancy functions are
`Y = mu(Psi(phi)) (I / Psi(phi)) p` and
`f = (mu(Psi(phi)) / 4) (Psi(phi) / I - I / Psi(phi))`, where `Psi` inverts the
steady relation `phi = Phi(I)`. The model reduces to `mu(I)` and `Phi(I)` in
isochoric flow, adds no parameters beyond those two curves, and is proved
linearly well-posed for all deformations and packing fractions. Pressure follows
from a local equation of state rather than a constraint. One-dimensional shear
cell computations converge under refinement and agree with transient DEM, while
the same computations with `mu(I), Phi(I)` blow up at a time that shrinks as the
grid refines. Appendix B shows that the stresses are derivatives of a convex
dissipation potential, so an implicit time step is a convex minimisation, the
structure on which the associated problem converged in the regularisation
study. The model covers only the inertial regime below the critical packing
`phi_c`: pressure is zero without shear unless the material is compacting,
there is no static yield stress, and pressure diverges as `phi -> phi_c`. The
authors leave the transition to the quasi-static regime, where stresses scale
with grain stiffness, to future work.

**Critical-state CIDR (Barker et al. 2017; Fannon, Moyles and Fowler 2019).** The
yield and dilatancy functions are `Y = alpha(I) p - p^2 / C(phi)` and
`f = beta(I) - 2 p / C(phi)`, with `alpha` and `beta` derived from `mu(I)` by
integrals that make the model isochoric at `|tau| = mu(I) p`, and
`C(phi) = (P / b) psi / (1 - psi)`, `psi = (phi - phi_min) / (phi_max - phi_min)`.
It is proved linearly well-posed in 2D and has a static yield stress, which
Schaeffer et al. (2019) name as the route outside the inertial regime. Fannon et
al. use it with the `mu(I)` parameters of Jop et al. for glass beads
(`mu_1 = tan 20.9 deg`, `mu_2 = tan 32.76 deg`, `I_0 = 0.279`, the values in
`bi_square_drum`) and find that its linear stability agrees with chute
experiments better than the incompressible model. Its packing function `C` is
described by those authors as speculative: `b` and the pressure scale `P` are
free, and the scale `rho_s g d` used by Barker et al. is not frame-indifferent
(Schaeffer et al. 2019 suggest a stiffness scale `k / d` instead). No
two-dimensional transient computation with it has been published, and whether
its stresses derive from a convex dissipation potential was not established in
the papers read.

**Critical-state `mu(I)` elastoplasticity (Blatny, Gray and Gaume 2024).**
Modified Cam-clay with Hencky elasticity, whose critical-state line depends on
`I` so that flow recovers `mu(I)`, implemented in 2D and 3D MPM and compared
with granular-collapse experiments. The collapse front converges under grid
refinement, but shear bands sharpen as the grid refines; the authors attribute
this to strain localisation and name gradient or nonlocal terms as remedies. It
needs a Young's modulus, so it replaces the rigid-grain formulation with
explicit or implicit elastoplasticity.

No published model combines rigid grains, a static yield stress, a proof of
well-posedness and validated two-dimensional transient computations.

## Sources read

- Barker, T., Schaeffer, D. G., Bohorquez, P. and Gray, J. M. N. T. (2015).
  Well-posed and ill-posed behaviour of the μ(I)-rheology for granular flow.
  *J. Fluid Mech.* 779, 794–818. https://doi.org/10.1017/jfm.2015.412
- Barker, T. and Gray, J. M. N. T. (2017). Partial regularisation of the
  incompressible μ(I)-rheology for granular flow. *J. Fluid Mech.* 828, 5–32.
  https://doi.org/10.1017/jfm.2017.428
- Barker, T., Schaeffer, D. G., Shearer, M. and Gray, J. M. N. T. (2017).
  Well-posed continuum equations for granular flow with compressibility and
  μ(I)-rheology. *Proc. R. Soc. A* 473, 20160846.
  https://doi.org/10.1098/rspa.2016.0846
- Heyman, J., Delannay, R., Tabuteau, H. and Valance, A. (2017). Compressibility
  regularizes the μ(I)-rheology for dense granular flows. *J. Fluid Mech.* 830,
  553–568. https://doi.org/10.1017/jfm.2017.612 (read as arXiv:1609.01502)
- Goddard, J. D. and Lee, J. (2017). On the stability of the μ(I) rheology for
  granular flow. *J. Fluid Mech.* 833, 302–331.
  https://doi.org/10.1017/jfm.2017.651 (read as arXiv:1710.08139)
- Li, S. and Henann, D. L. (2019). Material stability and instability in
  non-local continuum models for dense granular materials. *J. Fluid Mech.* 871,
  799–830. https://doi.org/10.1017/jfm.2019.311
- Henann, D. L. and Kamrin, K. (2016). A finite element implementation of the
  nonlocal granular rheology. *Int. J. Numer. Methods Eng.* 108, 273–302.
  https://doi.org/10.1002/nme.5213
- Henann, D. L. and Kamrin, K. (2013). A predictive, size-dependent continuum
  model for dense granular flows. *PNAS* 110, 6730–6735.
  https://doi.org/10.1073/pnas.1219153110 (read through a summary of the PMC
  full text; not used as evidence above)
- Dunatunga, S. and Kamrin, K. (2015). Continuum modelling and simulation of
  granular flows through their many phases. *J. Fluid Mech.* 779, 483–513.
  https://doi.org/10.1017/jfm.2015.383
- Daviet, G. and Bertails-Descoubes, F. (2016). A semi-implicit material point
  method for the continuum simulation of granular materials. *ACM Trans. Graph.*
  35(4), 102. https://doi.org/10.1145/2897824.2925877
- Daviet, G. (2016). *Modèles et algorithmes pour la simulation du contact
  frottant dans les matériaux complexes*. PhD thesis, Université Grenoble Alpes.
  https://gdaviet.fr/files/these-gdaviet-archivage.pdf (chapters 6 and 7)
- Rousseau, H., Métivet, T., Rousseau, R., Daviet, G. and Bertails-Descoubes, F.
  (2023). Revisiting the role of friction coefficients in granular collapses:
  confrontation of 3-D non-smooth simulations with experiments. *J. Fluid Mech.*
  975, A14. https://doi.org/10.1017/jfm.2023.835 (read as arXiv:2211.04865)
- Ionescu, I. R., Mangeney, A., Bouchut, F. and Roche, O. (2015). Viscoplastic
  modeling of granular column collapse with pressure-dependent rheology.
  *J. Non-Newton. Fluid Mech.* 219, 1–18.
  https://doi.org/10.1016/j.jnnfm.2015.02.006
- Klár, G., Gast, T., Pradhana, A., Fu, C., Schroeder, C., Jiang, C. and Teran, J.
  (2016). Drucker–Prager elastoplasticity for sand animation. *ACM Trans. Graph.*
  35(4). https://math.ucdavis.edu/~jteran/papers/KGPSJT16.pdf
- Gaume, J., Gast, T., Teran, J., van Herwijnen, A. and Jiang, C. (2018). Dynamic
  anticrack propagation in snow. *Nat. Commun.* 9, 3047.
  https://doi.org/10.1038/s41467-018-05181-w
- Schaeffer, D. G., Barker, T., Tsuji, D., Gremaud, P., Shearer, M. and Gray,
  J. M. N. T. (2019). Constitutive relations for compressible granular flow in
  the inertial regime. *J. Fluid Mech.* 874, 926–951.
  https://doi.org/10.1017/jfm.2019.476
- Fannon, J. S., Moyles, I. R. and Fowler, A. C. (2019). Application of the
  compressible I-dependent rheology to chute and shear flow instabilities.
  *J. Fluid Mech.* (read as the accepted manuscript, sections 1–3).
- Blatny, L., Gray, J. M. N. T. and Gaume, J. (2024). A critical state
  μ(I)-rheology model for cohesive granular flows. *J. Fluid Mech.* 997, A67.
  https://doi.org/10.1017/jfm.2024.643 (abstract, sections 1–2 and appendix C)
- Zhao, Y. and Choo, J. (2020). Stabilized material point methods for coupled
  large deformation and fluid flow in porous materials. *Comput. Methods Appl.
  Mech. Eng.* 362, 112742. (read as arXiv:1905.00671)
- Treskatis, T., Moyers-González, M. A. and Price, C. J. (2016). An accelerated
  dual proximal gradient method for applications in viscoplasticity.
  *J. Non-Newton. Fluid Mech.* 238, 115–130.
  https://doi.org/10.1016/j.jnnfm.2016.09.004 (read as arXiv:1509.05084)
- Gräser, C. and Sander, O. (2019). Truncated nonsmooth Newton multigrid methods
  for block-separable minimization problems. *IMA J. Numer. Anal.* 39(1),
  454–481. https://doi.org/10.1093/imanum/dry073 (read as arXiv:1709.04992)
- Sander, O. (2017). Solving primal plasticity increment problems in the time of
  a single predictor–corrector iteration. arXiv:1707.03733.
- Kornhuber, R., Krause, R., Sander, O., Deuflhard, P. and Ertel, S. (2008). A
  monotone multigrid solver for two body contact problems in biomechanics.
  *Comput. Visual. Sci.* (read as preprint; volume not confirmed)
- Bleyer, J. (2022). Applications of conic programming in non-smooth mechanics.
  arXiv:2202.00928.
- Song, Fan, Ascher and Pai (2026). A splitting architecture for exact reduced
  Coulomb friction. *Comput. Graph. Forum* 45(8). (DOI not confirmed)

Seen only as abstracts or metadata, and not used as evidence: Schaeffer (1987);
Kamrin and Koval (2012); Bleyer (2018), *Comput. Methods Appl. Mech. Eng.* 330;
Krabbenhoft, Lyamin, Huang and Vicente da Silva (2012), *Comput. Geotech.* 43;
Acary, Cadoux, Lemaréchal and Malick (2011), *ZAMM* 91; Daviet and
Bertails-Descoubes (2016), *J. Non-Newton. Fluid Mech.* 234; Lagrée, Staron and
Popinet (2011), *J. Fluid Mech.* 686; Chauchat and Médale (2014), *J. Comput.
Phys.* 256; Coombs et al. (2018), *Comput. Methods Appl. Mech. Eng.* 333.
