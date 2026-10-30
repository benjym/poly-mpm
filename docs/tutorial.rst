Create a simulation
===================

Inputs are Python modules under ``inputs/``. Each exposes ``Params(args)``;
``args[0]`` is the example name and subsequent entries are its CLI arguments.
The loader imports ``inputs.<name>`` and passes the complete argument list to
``Params``. It does not load arbitrary filesystem paths or JSON input files.
Declare positional arguments in ``inputs.ARGUMENTS`` (in ``inputs/__init__.py``),
including their converters and help text. The loader validates these before
constructing ``Params`` and exposes them through ``python MPM.py <name> --help``.
Inputs without a declaration accept no positional arguments.

A complete minimal input
------------------------

From the checkout root, copy the supplied example and run it::

    cp docs/examples/minimal.py inputs/my_sim.py
    uv run --locked python MPM.py my_sim

On Windows PowerShell use ``Copy-Item docs/examples/minimal.py inputs/my_sim.py``.
The input is intentionally short: 0.01 s of a falling elastic block. It writes
plots and particle arrays under ``im/my_sim/``. Its duration is too short for the
block to reach the floor; this checks setup and output rather than impact physics.

.. literalinclude:: examples/minimal.py
   :language: python

.. note::

   This is a setup example, not a general stable-timestep prescription. If you
   change grid spacing, stiffness, density, viscosity, or loading, reassess ``dt``
   and check timestep and spatial convergence.

Simulation parameters
---------------------

Every name an input sets must be a setting listed in :doc:`features`;
``initialise.prepare`` rejects any other name and suggests the nearest setting
for a misspelling. Helpers an input keeps on its objects must start with an
underscore; otherwise make them local variables.

``dt``, ``t_f`` and ``savetime`` are seconds. ``gravity`` is acceleration in
m/s²; negative values give downward gravity at ``gravity_angle=0``.
``gravity_angle`` is in radians, and the body acceleration is
``g * [sin(angle), cos(angle), 0]`` with ``g`` the gravity times
``gravity_ramp(t)``. ``gravity_ramp`` and ``gravity_angle`` may be functions of
time, evaluated at the start of each step; for example
``self.gravity_ramp = lambda t: 1 - np.exp(-3 * t**2 / 0.5**2)`` raises gravity
over about 0.5 s, and ``self.gravity_angle = lambda t: omega * t`` rotates it.
Example-specific arguments belong in the input.

``output_dir`` sets the output folder; the default is
``im/<example>/<first-material-law>/``. The run loop uses physical time and
stops at ``t_f``. A bounded run with ``max_steps`` stops earlier and still calls
its end-of-run output callback.

``time_stepping='dynamic'`` enables the timestep controller, with ``cfl``
and optional material ``critical_time(P)`` functions. It is not a universal
stability guarantee: the central controller estimates elastic, viscous,
advection, and segregation restrictions, but nonlinear and model-specific
stability still requires verification. Start with conservative values and
verify the behaviour of the chosen input.

Grid and boundaries
-------------------

``P.G`` describes a regular rectangular Cartesian grid. ``nx`` and ``ny`` count
**nodes**, not cells; spacing is computed from the extents and node counts.
The extents are ``x_min``, ``x_max``, ``y_min`` and ``y_max``. ``thickness``
defaults to 1 m. Coordinates are metres. Grain-size segregation also needs
``grain_sizes``: the representative positive diameter of each size class in
metres. The default is one class of diameter 0.001 m.

``P.B`` describes the boundaries. ``walls`` lists the domain edges with walls,
for example ``{"bottom", "left", "right"}``; each restrains the normal grid
velocity on its edge. ``wall_contact = "no_slip"`` also restrains the
tangential velocity on every wall and ``"no_slip_bottom"`` only on the bottom.
``periodic = True`` makes the left and right boundaries periodic. Particle positions must lie within the domain before cell lookup,
except for permitted periodic wrapping and the explicit inlet/outlet paths.

These settings impose grid constraints. They do not define curved contact surfaces
or independent contact between material bodies. Boundary tractions
(``traction_vertical``, ``traction_horizontal``) retain the
legacy reference-density scaling; see :doc:`array-conversion` before using them
for quantitative loading.

Materials and particle geometry
--------------------------------

``P.S`` can be one material namespace or a list of them. Each material needs:

* ``law``: an exact supported registry name, with its model parameters.
* ``n``, ``X``, ``Y``: particle count and initial coordinates, each coordinate
  array of shape ``(n,)``. Python lists are accepted and copied into NumPy arrays.
* ``rho``: initial density in kg/m³, positive and finite.
* ``A``: represented area **per particle**, in m². Volume is ``A*thickness`` and
  mass is ``rho*A*thickness``. Density and area can be scalar or length ``n``.
* Optional ``initial_v``: velocity in m/s, shape ``(3,)`` or ``(n, 3)``.
* For mixtures, ``phi`` of shape ``(ns,)`` or ``(n, ns)``, with ``ns`` the
  number of grain sizes; rows must be nonnegative and sum to one.

A material may set only the fields its law reads (listed in :doc:`features`);
any other field is rejected.

Two-by-two particles per grid cell is a convenient starting arrangement, not an
accuracy guarantee. Choose represented areas consistently with the occupied
geometry, then vary particle count and grid spacing independently in convergence
checks. The physical grains represented by size classes are distinct from the
numerical material points.

For elastic material, supply bulk modulus ``K`` and shear modulus ``G`` in Pa;
``E`` and ``nu`` are convenient local variables in the input, from which it
computes ``K`` and ``G``. The :doc:`materials` page lists parameters. ``elastic``
uses neo-Hookean elasticity, while ``von_mises`` uses a J2 return map; see
:doc:`numerical-foundations`.

To initialise confined stress, set ``initial_stress='compression'`` and a positive
``initial_pressure`` in Pa. This produces negative normal Cauchy stress. The
``lithostatic`` mode uses gravity, depth and slope. The particle ``pressure``
array is tension positive.

Grain-size segregation
----------------------

Start with ``inputs/bi_seg_test.py`` for a working configuration. Set
``segregate=True``, ``G.grain_sizes``, the composition ``phi`` and the
segregation coefficients ``segregation_rate`` and ``kinetic_pressure_scale``.
The current implementation transports conservative
species-mass exchanges; arbitrary clipping/renormalisation is not needed for the
validated cases. See :doc:`segregation-conservation` for equations, coefficient
roles, diagnostics, and limitations.

Inspecting state and customising output
---------------------------------------

``initialise.get_parameters`` returns ``(P, G, L)``. Access material arrays through
``L.materials[0].x``, ``.v``, ``.stress``, and ``.phi``. There are no scalar particle
objects at ``L.S[p][i]``. ``P.S`` is normalised to a material list.

Use ``MPM.run(P, G, L, max_steps=..., output=...)`` to advance the simulation and
clock. ``MPM.step`` is a lower-level single-step operator: it does not advance
``P.t`` or ``P.tstep`` and does not invoke output callbacks.

See :doc:`outputs` for figures, fields, file formats, and the callback interface.
For finite-domain integration, set ``quadrature='domain'`` and supply each
material's reference ``domain_size=(width_x,width_y)`` during construction. Its
product must match ``A``. The :doc:`quadrature` guide and ``domain_column`` input
show a complete setup; this mode requires a source checkout.

The :doc:`array-conversion` reference gives full array shapes, defaults, and
current numerical limitations, including fixed particle mass/volume and the
reconstructed-density diagnostic.
