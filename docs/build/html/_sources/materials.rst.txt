Choosing a material
===================

Use one of the exact registry names below as ``P.S.law`` (or on each entry in
``P.S`` when it is a list). Every material also needs the particle geometry,
initial density, and represented area described in :doc:`tutorial`. All tensor
components use SI units and tension-positive Cauchy stress.

.. list-table::
   :header-rows: 1
   :widths: 16 30 54

   * - Name
     - Required model parameters
     - Meaning / scope
   * - ``rigid``
     - Optional ``initial_v``
     - Prescribed velocity; no deformable stress update or new contact model.
   * - ``elastic``
     - ``K``, ``G``
     - Legacy incremental isotropic elasticity, with bulk/shear moduli in Pa.
   * - ``viscous``
     - ``K``, ``mu_s``, ``mu_v``
     - Legacy viscous law; viscosities in Pa s. Retains inherited volumetric
       normalisation/sign and stress-rate conventions; not identical to the
       Newtonian limit of the corrected Bingham law.
   * - ``von_mises``
     - ``K``, ``G``, ``k``, ``s``
     - Legacy von Mises H² incremental law. ``k`` is a stress scale in Pa and
       ``s`` its dimensionless exponent; this is not a standard radial return map.
   * - ``pouliquen``
     - ``K``, ``rho_s``, ``mu_0``, ``delta_mu``, ``I_0``, ``eta_max``
     - Retained local granular friction law, with legacy invariant conventions,
       viscosity cap, and density cutoff. Also supply ``mu_1`` and ``packing``
       for the example's steady-flow initialisation; ``mu_1=mu_0+delta_mu``.
   * - ``bingham``
     - ``K``, ``mu_s``, ``mu_0``, ``tau_0``; optional ``mu_v``
     - Corrected bi-viscosity approximation to Bingham flow. Small creep below
       yield is permitted; ``gamma_c`` is derived and need not be supplied.
   * - ``dp``
     - ``K``, ``G``, ``mu``, ``beta``
     - Cohesionless Drucker-Prager plasticity. ``mu`` is the ``q/p`` yield slope,
       not a friction angle; ``beta`` controls dilation.
   * - ``dp_rate``
     - As ``dp``, plus ``t_star``
     - Linear Perzyna regularisation. ``t_star`` is the shear-overstress relaxation
       time in seconds for zero dilation; zero gives the ``dp`` update.
   * - ``linear_mu``
     - ``K``, ``rho_s``, ``mu_0``, ``b``, ``eta_max``
     - Dense flowing ``mu(I)=mu_0+b I`` law with a viscosity cap. ``b=0`` gives
       constant friction. It is not a dilute/collisional granular model.

``K`` and ``G`` are Pa; ``rho_s`` is solid-grain density in kg/m³; grain diameters
``P.G.s`` are metres; ``eta_max`` and viscosities are Pa s; ``tau_0`` is Pa.
Friction coefficients, ``I_0``, ``b``, and ``beta`` are dimensionless. For flow laws,
``phi @ P.G.s`` supplies the local mean diameter where used.

The restored-law :doc:`constitutive audit <constitutive-models>` gives equations,
allowed parameter ranges, sign conventions, integration, and analytical checks.
The original five laws retain earlier equations and scalar-comparison tests;
those comparisons are not experimental validation or an audit of every physical
assumption. Select and calibrate a law for the intended regime rather than
switching names and retaining unrelated parameters.

Initialisation can require additional parameters: confined stress and the legacy
``non-zero``/``temperature`` modes may use elastic moduli; steady-flow options are
only implemented for particular laws. Start with a working input using the same
initialisation and inspect its setup before changing it. New algebraic Bingham
and linear-friction fluids do not need an elastic shear modulus.

``elastic_time``, ``ken_kamrin``, and ``ken_simple`` are retired with explicit
reasons in the audit. Their archived source is not imported. ``viscous_size`` and
the separate legacy ``P.F`` compressible-fluid path remain unimplemented.
