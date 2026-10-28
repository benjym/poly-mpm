"""Opt-in GPU support. Importing this package never imports SlangPy.

The Metal runtime and primitives are experimental. The resident MPM runner is
not enabled until the hardware validation gates documented in docs/gpu.md pass.
"""

SLANGPY_VERSION = "0.43.1"


class GPUError(RuntimeError):
    """An explicit GPU failure; never triggers a CPU solver fallback."""


class GPUUnavailable(GPUError):
    """The required platform, dependency, compiler or device is unavailable."""


class ValidationPending(GPUError):
    """A dependent implementation stage has not passed its hardware gate."""
