"""Real Metal tests. Opt in with POLY_MPM_TEST_METAL=1 on Apple silicon.

An opted-in but unavailable device FAILS; it is never silently skipped/faked.
"""

import os

import pytest

from mpm_gpu.runtime import MetalRuntime
from mpm_gpu.validate import validate_runtime

pytestmark = pytest.mark.skipif(
    os.environ.get("POLY_MPM_TEST_METAL") != "1",
    reason="requires explicit Apple silicon Metal hardware validation",
)


def test_runtime_gate_on_actual_metal_device():
    with MetalRuntime(debug=True) as runtime:
        report = validate_runtime(runtime)
    failures = [row for row in report["cases"] if row["status"] != "passed"]
    assert not failures, failures


def test_mechanics_operator_on_actual_metal_device():
    from mpm_gpu.validate import validate_mechanics

    with MetalRuntime(debug=True) as runtime:
        report = validate_mechanics(runtime)
    failures = [row for row in report["cases"] if row["status"] != "passed"]
    assert not failures, failures
