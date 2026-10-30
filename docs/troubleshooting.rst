Troubleshooting and limitations
===============================

Installation and command errors
-------------------------------

``uv`` cannot find the pinned Python
    Update uv using the installation method you originally used, then retry
    ``uv sync --extra dev --locked``. If your platform has no compatible
    interpreter or dependency wheel, inspect that package's error rather than
    silently changing numerical dependencies. Python 3.14.10 is the tested pin;
    an arbitrary later minor version is not necessarily supported by every wheel.

``ModuleNotFoundError`` or ``poly-mpm: command not found``
    From the checkout root, run ``uv sync --extra dev --locked`` and use
    ``uv run --locked python MPM.py ball``. For a pip environment, activate it
    first. ``python -m pip`` may be absent from a uv-created environment; use
    ``uv sync`` to manage that environment.

The example cannot be imported
    Use its module name, such as ``ball``, rather than ``inputs/ball.py``. Put new
    modules under ``inputs/`` and expose ``Params(args)``. Run the documented
    commands from the checkout root.

An example reports a usage error
    Run ``uv run --locked python MPM.py --help`` to list examples, or
    ``uv run --locked python MPM.py <example> --help`` for its argument order,
    types, and choices. Missing, invalid, and extra arguments are rejected before
    setup. The Python loader raises ``ValueError`` for these errors. Check the
    :doc:`example matrix <array-conversion>` for runnable commands.

    For the square drum, use
    ``uv run --locked python MPM.py bi_square_drum 11 2``. The arguments are the
    number of grid nodes per side (at least 3) and the finite, positive Froude
    number. For ``bi_dam_break``, the arrangement must be ``top``, ``bot``, or
    ``mix``. Grain size class counts must be at least two. ``silo`` accepts only
    its grid node count: use ``silo 11``; its two size classes are fixed, and the
    formerly documented trailing ``2`` was ignored.

``unknown setting``
    Each name an input sets on ``P``, ``P.G``, ``P.B``, ``P.O`` or a material must
    be a setting listed in :doc:`features`; the error names the closest setting
    when the name looks like a misspelling. Keep helper values as local variables
    or give them names that start with an underscore. The error ``... is required``
    names a setting with no default, and an input that defines ``update_forces`` is
    rejected: use ``gravity_ramp``, ``gravity_angle`` and ``B.traction_vertical``
    or ``B.traction_horizontal``, which may be functions of time.

``Unsupported material law``
    Check the exact registry name and :doc:`constitutive-models`. The retired
    ``elastic_time``, ``ken_kamrin``, and ``ken_simple`` implementations are not
    available. ``viscous_size`` and a separate compressible-fluid phase are
    unimplemented.

Output and runtime
------------------

No plot window appears
    Expected: figures are saved as PNGs. Check ``output_dir``, ``O.figures``,
    ``O.fields`` and whether the input replaces a default callback. ``output=False`` disables
    all output, including custom histories.

No files after a short preview
    The preview may end before ``savetime``, and an input's custom callback may
    save only at save times. Set ``P.savetime=P.dt`` for a few steps and choose a
    separate directory, as in :doc:`getting_started`.

A run appears to take forever
    Check ``dt``, ``t_f``, and plot frequency in the input. Some examples
    request millions of steps. Start with ``run(..., max_steps=100, output=False)``.
    The CLI has no generic duration or maximum-step switch. Ctrl-C interrupts a
    run and does not guarantee end callbacks run; set ``O.save_checkpoint`` to
    keep a restart file from the last save time (see :doc:`outputs`).

A Matplotlib cache warning appears
    Matplotlib needs a writable cache/config directory. If your home directory
    is read-only, set ``MPLCONFIGDIR`` to a writable directory before starting
    Python; for example on macOS/Linux,
    ``MPLCONFIGDIR=/tmp/poly-mpm-mpl uv run --locked python MPM.py ball``.

Numerical errors
----------------

A particle leaves the domain
    The solver fails before invalid cell lookup. Check domain extents, inlet/
    outlet or periodic configuration, velocity and loading, then timestep. Moving
    beyond an unconstrained grid edge is not an implicit free-surface boundary
    extension or an automatically expanding computational domain.

A model reports nonfinite state
    The error identifies the law, material, timestep, field, and particle rows.
    Check parameter units, initial conditions, model domain, and timestep. Reduce
    the timestep as a diagnostic, then perform convergence checks. Finite output
    alone is not proof of an accurate or physically valid solution.

Grain fractions do not sum to one
    Initial rows must be nonnegative and sum to one. Use ``phi`` of shape
    ``(ns,)`` or ``(n, ns)``. The segregation tests conserve species mass
    without clipping or renormalisation. Investigate inventory and boundary fluxes
    rather than enabling ``segregation_normalise`` to hide drift.

What the solver currently does
------------------------------

This is a research solver with a regular 2D grid and explicit mechanics sharing
one nodal momentum field. It evolves deformation, volume and density and
supports linear/quadratic interpolation plus FLIP/PIC/APIC transfers. APIC requires
a quadratic basis and planar motion. There is no GIMP, independent multi-body
contact, general curved boundary or implicit/quasistatic solve.
Rate-form plasticity is not a multiplicative finite-strain plasticity model.

See :doc:`numerical-foundations`, :doc:`array-conversion`, :doc:`constitutive-models`,
and :doc:`segregation-conservation` for exact scope. Smoke-tested examples, equation
checks, conservation tests, experimental validation, and convergence studies are
different levels of evidence. Choose the evidence required for your application.
