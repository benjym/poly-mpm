.. _getting_started:

Getting started
***************

Installation
============

Python 3.12 is the tested baseline. Install in a virtual environment from a checkout::

    python3.12 -m venv .venv
    source .venv/bin/activate
    python -m pip install -e '.[dev]'

Running examples
================

Each input module supplies one simulation's parameters::

    python MPM.py ball
    python MPM.py bi_seg_test 22 2 11

Default output is ``im/<example>/<first-material-law>/``. An input may override
this with ``supername`` or ``save_dir``. Original examples may request long runs.

Supported laws are rigid, elastic, viscous, von Mises, and Pouliquen. Unsupported
experimental laws fail at initialisation. See ``docs/array-conversion.md`` for
the complete example matrix, numerical corrections, and remaining limitations.

Array state and tests
=====================

``initialise.get_parameters(args)`` returns ``(P, G, L)``. ``L.materials`` contains
one array collection per material; ``P.S`` is a material list. All output callbacks
receive ``(L, P, G)``. Run bounded simulations with::

    from initialise import get_parameters
    from MPM import run
    run(*get_parameters(['ball']), max_steps=100, output=False)

Run the test suite and the benchmark from the repository root::

    python -m pytest -q -W error::RuntimeWarning
    python -m benchmarks.audit_examples
    python -m benchmarks.compare
