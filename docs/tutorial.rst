Create a simulation
===================

Inputs are Python modules under ``inputs/``. Each exposes ``Params(args)``;
``args[0]`` is the example name and subsequent entries are its CLI arguments.
The loader imports ``inputs.<name>`` and passes the complete argument list to
``Params``. It does not load arbitrary filesystem paths or JSON input files.

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

``dt``, ``t_f``, and ``savetime`` are seconds. ``max_g`` is acceleration in m/s²;
negative values give downward gravity at ``theta=0``. Internal ``theta`` is in
radians. Body acceleration is ``g * [sin(theta), cos(theta), 0]``. ``g`` defaults
to ``max_g``; an optional zero-argument ``update_forces`` method may update it
using ``self.t``. Ramped loads and example-specific arguments belong in the input.

``save_dir`` sets the output folder directly. If absent, ``supername`` is used;
otherwise the default is ``im/<example>/<first-material-law>/``. ``nt`` need not
be supplied: the run loop uses physical time and stops at ``t_f``. A bounded run
with ``max_steps`` stops earlier and still calls its end-of-run output callback.

``time_stepping='dynamic'`` enables the existing timestep controller, with ``CFL``
and optional material ``critical_time(P)`` functions. It is not a universal
stability guarantee: the current controller does not automatically derive every
elastic, viscous, and transport restriction. Start with conservative values and
verify the behaviour of the chosen input.

Grid and boundaries
-------------------

``P.G`` describes a regular rectangular Cartesian grid. ``nx`` and ``ny`` count
**nodes**, not cells; spacing is computed from the extents and node counts.
``thickness`` defaults to 1 m. Coordinates are metres. Grain-size segregation
also needs ``ns`` and ``s``: the number of size classes and their representative
positive diameters in metres. The default is one class of diameter 0.001 m.

``P.B`` contains boundary flags. Missing flags default to false. ``has_bottom``
restrains normal motion at the bottom; add ``no_slip_bottom`` to restrain
horizontal motion there. ``has_top``, ``has_left``, and ``has_right`` similarly
restrain normal grid velocity. ``cyclic_lr`` makes the left/right boundary
periodic. Particle positions must lie within the domain before cell lookup,
except for permitted periodic wrapping and the explicit inlet/outlet paths.

These flags impose grid constraints. They do not define curved contact surfaces
or independent contact between material bodies. Traction flags retain the
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
* For mixtures, ``phi`` of shape ``(ns,)`` or ``(n, ns)``; rows must be
  nonnegative and sum to one. ``PHI`` is an accepted higher-priority alias.

Two-by-two particles per grid cell is a convenient starting arrangement, not an
accuracy guarantee. Choose represented areas consistently with the occupied
geometry, then vary particle count and grid spacing independently in convergence
checks. The physical grains represented by size classes are distinct from the
numerical material points.

For elastic material, supply bulk modulus ``K`` and shear modulus ``G`` in Pa;
``E`` and ``nu`` are convenient input variables, but the loader does not derive
``K`` or ``G`` from them. The :doc:`constitutive-models` page describes the restored
models and their parameter restrictions. For the original five laws, consult the
chosen input and ``constit.py`` as well: they retain legacy conventions.

To initialise confined stress, set ``pressure='compression'`` and a positive
``initial_pressure`` in Pa. This produces negative normal Cauchy stress. The
``lithostatic`` mode uses gravity, depth, and slope. Do not confuse these input
mode strings with the particle ``pressure`` array, which is tension positive.

Grain-size segregation
----------------------

Start with ``inputs/bi_seg_test.py`` for a working configuration. Set
``segregate=True`` (or the compatible alias ``segregate_grid=True``), grain sizes,
composition, and the model's segregation coefficients ``c`` and ``l``.
``segregate_mp`` is unsupported. The current implementation transports conservative
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

See :doc:`outputs` for flags, file formats, and the callback interface.
The :doc:`array-conversion` reference gives full array shapes, defaults, and
current numerical limitations, including fixed particle mass/volume and the
reconstructed-density diagnostic.
