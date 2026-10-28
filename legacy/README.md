# Archived scalar implementations

`constit_scalar.py` and `integrators.py` preserve the unfinished `np-array`
branch's constitutive and segregation source at commit `b4e2f4f`. They are
reference material, are not imported by the solver, and are not included in the
installed package. Their old validation claims are historical, not current
support guarantees.

The active registry supports rigid, elastic, viscous, von Mises, and Pouliquen
materials. Experimental laws must be ported and tested before registry admission.
The independent, corrected scalar comparison kernels live in `tests/reference.py`.
