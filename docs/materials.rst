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
     - Prescribed motion with no deformable stress evolution.
   * - ``elastic``
     - ``K``, ``G``
     - Compressible neo-Hookean elasticity, with bulk/shear moduli in Pa. Legacy mode
       retains incremental isotropic elasticity.
   * - ``viscous``
     - ``K``, ``eta``, ``eta_v``
     - Newtonian fluid; viscosities in Pa s. The zero-yield Bingham limit gives the same stress law.
   * - ``von_mises``
     - ``K``, ``G``, ``yield_stress``
     - Perfect J2 return map; ``yield_stress`` is the von Mises equivalent stress at
       yield, in Pa. ``s`` only controls the legacy H² law.
   * - ``pouliquen``
     - ``K``, ``rho_s``, ``mu_s``, ``mu_d``, ``I_0``, ``eta_max``
     - Local granular friction law ``mu(I)=mu_s+(mu_d-mu_s) I/(I_0+I)`` with a viscosity cap.
       Uses engineering shear rate and current physical density. Also supply the
       packing fraction ``phi_flow`` for the example's steady-flow initialisation.
   * - ``bingham``
     - ``K``, ``eta``, ``eta_0``, ``yield_stress``; optional ``eta_v``
     - Bi-viscosity approximation to Bingham flow: plastic viscosity ``eta``,
       unyielded viscosity ``eta_0``. Small creep below
       yield is permitted; ``gamma_c`` is derived and need not be supplied.
   * - ``dp``
     - ``K``, ``G``, ``mu_s``, ``beta``
     - Cohesionless Drucker-Prager plasticity. ``mu_s`` is the ``q/p`` yield slope,
       not a friction angle; ``beta`` controls dilation.
   * - ``dp_rate``
     - As ``dp``, plus ``t_star``
     - Linear Perzyna regularisation. ``t_star`` is the shear-overstress relaxation
       time in seconds for zero dilation; zero gives the ``dp`` update.
   * - ``linear_mu``
     - ``K``, ``rho_s``, ``mu_s``, ``b``, ``eta_max``
     - Dense flowing ``mu(I)=mu_s+b I`` law with a viscosity cap. ``b=0`` gives
       constant friction. It is not a dilute/collisional granular model.
   * - ``granular``
     - ``rho_s``, ``mu_s``, ``phi_max`` (at most 1); optional ``critical_time(P)``
     - Rigid grains with constant Coulomb friction and a maximum packing
       fraction, solved as a contact problem. Requires ``solver='granular'`` and
       is the only material in the run. See :doc:`granular-prototype`.

``K`` and ``G`` are Pa; ``rho_s`` is solid-grain density in kg/m³; grain diameters
``P.G.s`` are metres; ``eta_max`` and viscosities are Pa s; ``yield_stress`` is Pa.
Friction coefficients, ``I_0``, ``b``, and ``beta`` are dimensionless. For flow laws,
``phi @ P.G.s`` supplies the local mean diameter where used.

The :doc:`constitutive reference <constitutive-models>` gives equations,
allowed parameter ranges, sign conventions, integration, and analytical checks.
See :doc:`numerical-foundations` for deformation, transfer, and diagnostic conventions.
These checks are not experimental validation. Select and calibrate a law for the intended regime rather than
switching names and retaining unrelated parameters.

Initialisation can require additional parameters: confined stress and the legacy
``non-zero``/``temperature`` modes may use elastic moduli; steady-flow options are
only implemented for particular laws. Start with a working input using the same
initialisation and inspect its setup before changing it. Algebraic Bingham
and linear-friction fluids do not need an elastic shear modulus.

``elastic_time``, ``ken_kamrin``, and ``ken_simple`` are retired with explicit
reasons in the constitutive reference. Their archived source is not imported. ``viscous_size`` and
the separate legacy ``P.F`` compressible-fluid path remain unimplemented.
