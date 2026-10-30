Output and Python API
=====================

Where files go
--------------

Output paths are relative to the process's working directory. ``P.output_dir``
sets the output directory; the default is ``im/<example>/<first-material-law>/``. Choose a different directory for each run
if you want to preserve results. Existing filenames may be overwritten.

Matplotlib uses the noninteractive Agg backend, so no plot window is expected.
Open saved PNGs in an image viewer. Simulation output is ignored by Git.

Figures and fields
------------------

``P.O.figures`` and ``P.O.fields`` list the figures and arrays to write, for
example ``self.figures = {"continuum", "material_points"}``. Both default to
empty. The other output settings are listed in :doc:`features`. The default
adapter writes the requested figures and arrays at each ``savetime`` and at
the end. A custom callback replaces the corresponding default stage, so some
examples intentionally produce different output.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Setting
     - Result under the output directory
   * - figure ``continuum``
     - ``Continuum/*.png``: velocity, density, pressure, shear rate, mean size.
   * - figure ``material_points``
     - ``MP/*.png``: material identity, speed, pressure.
   * - figure ``gsd_material_points`` / ``gsd_grid``
     - ``GSD/*.png``: grain fractions on particles or grid.
   * - figure ``shear_rate``
     - ``MP/gammadot_*.png``: shear rate on the particles.
   * - field ``u``
     - ``data/u_*.npy`` and ``data/v_*.npy``, each ``(ny, nx)``.
   * - field ``density`` / ``s_bar``
     - ``data/density_*.npy`` / ``data/s_bar_*.npy``, each ``(ny, nx)``.
   * - field ``phi_material_points``
     - ``data/MP_phi_*.npy``, columns ``x, y, phi_0, ..., phi_(ns-1)``.
   * - ``measure_energy``
     - Default callback records time, kinetic, potential, and the legacy strain-
       energy estimate in ``data/energy.npy``; ``energy.png`` plots them.
   * - ``measure_stiffness``
     - Default end callback saves ``data/stiffness_material_*.npz``.

``O.save_diagnostics`` defaults to true and saves the
separate verification ledger to ``data/diagnostics.json``. See
:doc:`numerical-foundations` for work budgets, sampling, and limitations.

The legacy strain-energy estimate is not a validated stored-energy or dissipation
budget for every nonlinear material. Likewise, the saved grid density is a
reconstruction from nodal mass, not the physical particle density evolved through ``F``.

Read a saved particle array with NumPy::

    import numpy as np
    from pathlib import Path

    path = next(Path('im/my_sim/data').glob('MP_phi_*.npy'))
    particles = np.load(path)
    xy = particles[:, :2]
    phi = particles[:, 2:]
    print(xy.shape, phi.shape, phi.sum(axis=1))

The DP loading examples save ``data/pq_history.npz`` and ``pq_history.png``.
Here ``p`` is compression positive and ``q=sqrt(3 J2)``. Other saved ``pressure``
fields use the solver's tension-positive convention. Bingham's history pressure
excludes the instantaneous bulk-viscous contribution to mean Cauchy stress.

Callbacks
---------

All callbacks receive ``(L, P, G)``:

* ``at_start``: once before the run loop, including a bounded run.
* ``after_every_timestep``: after state and clock advancement.
* ``after_every_savetime``: when a step crosses the next requested save time.
* ``at_end``: once when the run ends, including a ``max_steps`` stop.

A callback may be set on ``P.O`` as a function or defined as a method of the
input's output class. A custom callback replaces that stage's default adapter; it
does not automatically add default plotting or energy recording. Explicitly call
``plotting.Plotting`` methods if both behaviours are needed.

For example, define this method on the input's output class::

    def after_every_savetime(self, L, P, G):
        from pathlib import Path
        import numpy as np
        from plotting import Plotting

        Plotting().draw_material_points(L, P, G)
        np.save(Path(P.save_dir)/f'stress_{P.tstep:06d}.npy',
                L.materials[0].stress)

Setting ``run(..., output=False)`` suppresses every callback. Large output
frequencies can dominate runtime and disk usage; a small timestep does not imply
that every step needs a plot.

Inspecting live state
---------------------

``L.materials`` contains one batch per material. Position and velocity arrays have
shape ``(n, 3)``, stress and strain ``(n, 3, 3)``, and composition ``(n, ns)``.
Geometry is two dimensional; the three-component/tensor representation is not a
3D solver. ``P.S`` holds material parameters, while ``L.materials`` holds evolving
particle state.

Grid constitutive diagnostics such as ``G.pressure`` are mass-weighted sums.
Use ``G.physical('pressure')`` to obtain the physical field. ``G.phi``,
``G.s_bar``, and ``G.shear_rate`` are already physical fields; do not divide them
by mass again. Periodic endpoint columns duplicate one physical grid node.

``run`` advances ``P.t`` and ``P.tstep``. ``step`` intentionally does not.

Checkpoint and restart
----------------------

``O.save_checkpoint = True`` writes ``data/checkpoint.pkl`` at each save time
and at the end of the run, replacing the previous file. ``checkpoint.save`` and
``checkpoint.load`` do the same from Python. A checkpoint holds the evolving
state of ``P`` (clock, timestep, save schedule), the grid, every particle batch
and the diagnostics. Inputs, materials, callbacks and controllers are rebuilt
from the same input, so restart with the same example and arguments::

    import checkpoint
    from initialise import get_parameters
    from MPM import run

    P, G, L = get_parameters(['bi_seg_test', '22', '2', '11'])
    checkpoint.load(P.save_dir + 'data/checkpoint.pkl', P, G, L)
    run(P, G, L)

On the CPU a restarted run is bitwise identical to an uninterrupted one, checked
for elastic, quadratic APIC, segregation, finite-domain and rigid-grain runs with
one-sided and frictional walls. A Metal restart re-uploads the host state; in the
frictional-wall column it reproduced the uninterrupted run exactly. The file is a
Python pickle: load only checkpoints you wrote.

Diagnostics also record the selected ``interpolation``, ``quadrature`` and
``transfer``, plus ``maximum_stencil_entries``. Finite-domain
support varies with particle geometry; see :doc:`quadrature`.
