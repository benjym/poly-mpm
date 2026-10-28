Getting started
===============

Requirements and installation
-----------------------------

Use a checkout of this repository, Git, and a current version of
`uv <https://docs.astral.sh/uv/getting-started/installation/>`_. Python 3.14 is the
tested baseline; ``.python-version`` pins 3.14.10. uv can download the pinned
interpreter if it is absent. Initial setup requires network access and disk
space for Python, NumPy, SciPy, Astropy, and Matplotlib.

From a terminal::

    git clone https://github.com/benjym/poly-mpm.git
    cd poly-mpm
    uv sync --extra dev --locked
    uv run --locked python -c "import numpy; from constit import SUPPORTED; print(SUPPORTED)"

``uv sync`` installs the project as well as its dependencies. ``--extra dev``
includes pytest; omit it if you only want to run simulations. Run the commands
below from the checkout root. The snippets use ``uv run`` so activation is not
required. An existing checkout should use the branch containing these docs;
older branches may still contain the scalar implementation.

Alternatively, with Python 3.14 already installed::

    python3.14 -m venv .venv
    source .venv/bin/activate
    python -m pip install -e '.[dev]'

On Windows PowerShell, activate with ``.venv\Scripts\Activate.ps1`` instead.
For this pip workflow, replace ``uv run --locked python`` below with ``python``.
The pip install resolves compatible dependency versions; ``uv.lock`` records the
exact tested environment. Updating Python alone does not modernise an old checkout.

First complete simulation
-------------------------

Run the small elastic ball example::

    uv run --locked python MPM.py ball

It simulates one second with 500 particles and 1,000 timesteps. The terminal shows
progress and ends with a newline. Matplotlib runs without a GUI; figures are saved
to disk, not opened in a window. Expect ``Continuum/*.png`` and ``MP/*.png``
under ``im/ball/elastic/``. Wall-clock duration depends on your machine and plot
rendering. Repeating a run in the same directory can overwrite its files.

The installed command is equivalent::

    uv run --locked poly-mpm ball

A bounded preview
-----------------

Some examples request millions of steps. The CLI takes an example name and its
positional arguments; it has no ``--max-steps`` or ``--output`` options. Use the
Python API for a short preview and a separate output directory::

    uv run --locked python - <<'PY'
    from initialise import get_parameters
    from MPM import run

    P, G, L = get_parameters(['ball'])
    P.save_dir = 'im/ball-preview/'
    P.savetime = P.dt
    run(P, G, L, max_steps=3)
    print(P.t, P.tstep, P.save_dir)
    PY

This shell heredoc works in macOS/Linux shells. On Windows, put its Python body
in ``preview.py`` and run ``uv run --locked python preview.py``. Set
``output=False`` in ``run`` to suppress all output callbacks, including custom
history recording. In that mode there will be no saved plots or arrays.

Choosing an example
-------------------

* ``ball``: a small first run with elastic material and saved plots.
* ``dp_unit_test`` / ``dp_rate_unit_test``: Drucker-Prager loading histories;
  full runs have 50,000 / 5,000 steps respectively.
* ``bi_seg_test 22 2 11``: segregation at 22 degrees, two grain-size classes,
  and 11 vertical grid nodes. Replace ``2`` with ``8`` for eight classes. The
  input converts its angle argument to its internal signed-radian convention.
  This is a long research example; start with a bounded run.

The :doc:`complete example matrix <array-conversion>` gives all required
arguments and identifies unsupported inputs. Examples with arguments must receive
them in the listed order; there is no general argument parser for their settings.
Bounded smoke checks establish compatibility, not stability over every original
run duration or physical calibration.

Supported material names are ``rigid``, ``elastic``, ``viscous``, ``von_mises``,
``pouliquen``, ``bingham``, ``dp``, ``dp_rate``, and ``linear_mu``. See
:doc:`constitutive-models` before selecting parameters. Retired and missing laws
fail at initialisation; they are not silently replaced.

Checking the installation
-------------------------

With the dev extra installed::

    uv run --locked --extra dev python -m pytest -q -W error::RuntimeWarning
    uv run --locked python -m benchmarks.audit_examples --output /tmp/poly-mpm-examples.json

The tests and benchmark modules require the source checkout. The audit performs
five bounded timesteps per supported example; its output file records each
example's status. On Windows, choose a writable output path instead of ``/tmp``.
