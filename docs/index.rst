poly-mpm user guide
===================

poly-mpm is a small two-dimensional material point method solver for solid and
viscous/granular flow, including polydisperse grain-size segregation. It uses
NumPy particle arrays and a regular Cartesian grid with selectable linear or quadratic interpolation.

Start with :doc:`getting_started`, then use :doc:`tutorial` to create an input.
:doc:`features` lists every setting with its default, the solvers and backends
that support it and the page that describes it.
The numerical reports explain what has been checked and what remains limited;
the examples are research fixtures, not a catalogue of calibrated materials.

.. toctree::
   :maxdepth: 2
   :caption: Using the solver

   getting_started
   tutorial
   features
   materials
   outputs
   troubleshooting

.. toctree::
   :maxdepth: 2
   :caption: Methods and reference

   numerical-foundations
   granular-prototype
   granular-source-audit
   granular-mixed
   connected-flows-literature
   gpu
   gpu-handoff
   quadrature
   verification
   implementation-todo
   releasing
   constitutive-models
   array-conversion
   segregation-conservation
   code
   README

* :ref:`genindex`
* :ref:`modindex`
