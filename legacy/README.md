# Archived scalar implementations

`constit_scalar.py` and `integrators.py` preserve the unfinished `np-array`
branch's constitutive and segregation source at commit `b4e2f4f`. They are
reference material, are not imported by the solver, and are not included in the
installed package. Their old validation claims are historical, not current
support guarantees.

The active registry now also supports Bingham, Drucker–Prager, rate-dependent
Drucker–Prager, and linear friction. Their kernels are physically corrected
implementations, not parity copies of this archive. `elastic_time`, `ken_kamrin`,
and `ken_simple` are explicitly retired; they are retained here only to preserve
provenance. See the [constitutive audit](../docs/constitutive-models.md) for reasons,
replacement options, equations, and parameter conventions.
The independent, corrected scalar comparison kernels for the original five laws
live in `tests/reference.py`.
