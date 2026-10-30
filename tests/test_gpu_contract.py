"""CPU-verifiable GPU contracts. These tests are not evidence of Metal execution."""

import json
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose
import pytest

from inputs.granular_column import Params
from initialise import prepare
from mpm_gpu import GPUError, GPUUnavailable
from mpm_gpu.config import defaults
from mpm_gpu.primitives import MAX_COUNT, check_count
from mpm_gpu.runtime import check_platform
from mpm_gpu.topology import q1_topology, split_positions, join_positions
from mpm_gpu.validate import cpu_checks, main


def metal_params():
    P = Params([])
    P.backend = "metal"
    return P


def test_cpu_import_does_not_load_native_gpu_dependency():
    # Test in a fresh interpreter; other tests may have imported modules already.
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import mpm_gpu, mpm_gpu.runtime, mpm_gpu.validate; import sys; assert 'slangpy' not in sys.modules",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_cpu_backend_and_numerical_defaults_are_preserved():
    P, _, _ = prepare(Params([]))
    assert P.backend == "cpu"
    assert P.granular_tolerance == 1e-9
    # The example names its Coulomb solver; without it the default is CIDR, whose
    # material fields this Coulomb column lacks.
    P = Params([])
    del P.granular_solver
    with pytest.raises(ValueError, match="cidr"):
        prepare(P)


def test_metal_defaults_and_explicit_values():
    P = SimpleNamespace(backend="metal")
    defaults(P)
    assert P.granular_tolerance == 1e-4
    assert P.granular_cidr_tolerance == 1e-4
    P = SimpleNamespace(backend="metal", granular_tolerance=2e-5, granular_cidr_tolerance=3e-5)
    defaults(P)
    assert P.granular_tolerance == 2e-5
    assert P.granular_cidr_tolerance == 3e-5


@pytest.mark.parametrize("backend", ["cuda", "auto", None, "Metal"])
def test_unknown_backend_is_rejected(backend):
    P = Params([])
    P.backend = backend
    with pytest.raises(ValueError, match="backend must"):
        prepare(P)


def test_metal_opt_in_cannot_fall_through_to_cpu(monkeypatch):
    """A Metal step either runs on the device or fails; it never calls the CPU step."""
    import granular
    from MPM import step

    def refuse(*args, **kwargs):
        raise AssertionError("the CPU granular step ran for a Metal backend")

    monkeypatch.setattr(granular, "step", refuse)
    P, G, L = prepare(metal_params())
    try:
        step(P, G, L)
    except GPUError:
        pass  # no usable Metal device on this host


def test_switching_an_already_prepared_state_is_validated():
    from MPM import run, step

    P, G, L = prepare(Params([]))
    original = L.materials[0].x.copy()
    P.backend = "metal"
    for function in (step, run):
        # The CPU defaults (hybrid solver, 1e-9 tolerance) are not Metal configurations.
        with pytest.raises(ValueError, match="Metal"):
            function(P, G, L)
        assert_allclose(L.materials[0].x, original, rtol=0, atol=0)


@pytest.mark.parametrize("tolerance", [1e-9, 9.9e-6, float("nan"), float("inf"), 0])
def test_metal_precision_contract_never_silently_loosens_tolerance(tolerance):
    P = metal_params()
    P.granular_tolerance = tolerance
    with pytest.raises(ValueError, match="tolerance"):
        prepare(P)
    if np.isfinite(tolerance):
        assert P.granular_tolerance == tolerance


@pytest.mark.parametrize("tolerance", [1e-5, 2e-5, 1e-4])
def test_supported_tolerance_reaches_hardware_gate(tolerance):
    P = metal_params()
    P.granular_tolerance = tolerance
    prepare(P)
    assert P.granular_tolerance == tolerance


@pytest.mark.parametrize(
    "name,value",
    [
        ("granular_solver", "hybrid"),
        ("granular_stress_basis", "P0"),
        ("granular_stress", "plane_strain"),
        ("granular_viscosity", 0.01),
    ],
)
def test_unsupported_metal_features_fail_before_native_initialization(name, value):
    P = metal_params()
    setattr(P, name, value)
    with pytest.raises(ValueError, match="Metal"):
        prepare(P)


def test_metal_segregation_is_accepted_with_supported_settings():
    P = metal_params()
    P.G.grain_sizes = [0.001, 0.003]
    P.S.phi = np.tile([0.5, 0.5], (P.S.n, 1))
    P.segregate, P.segregation_rate = True, 1e-2
    prepare(P)
    P = metal_params()
    P.segregate, P.smooth_shear_rate = True, True
    with pytest.raises(ValueError, match="smooth_shear_rate"):
        prepare(P)


@pytest.mark.parametrize(
    "name,value",
    [("inlets", {"top"}), ("inlets", {"right"}), ("outlets", {"bottom"}), ("outlets", {"left"})]
    + [("obstacles", {"silo_bottom"})],
)
def test_particle_count_changing_boundaries_accepted(name, value):
    # The runner applies them on the host (runner_checks.boundary_events_check).
    P = metal_params()
    setattr(P.B, name, value)
    prepare(P)


@pytest.mark.parametrize("name", ["conveyor_speed", "traction_vertical", "traction_horizontal"])
def test_unsupported_boundaries_rejected(name):
    P = metal_params()
    setattr(P.B, name, 1.0)
    with pytest.raises(ValueError):
        prepare(P)


@pytest.mark.parametrize(
    "nx,ny,periodic", [(2, 2, False), (2, 4, True), (3, 3, True), (4, 6, True), (9, 7, False), (9, 7, True)]
)
def test_topology_matches_independent_q1_cell_incidence(nx, ny, periodic):
    T = q1_topology(nx, ny, periodic=periodic)
    adjacency = {int(node): set() for node in T.nodes}
    canonical = T.canonical.reshape(ny, nx)
    for y in range(ny - 1):
        for x in range(nx - 1):
            corners = set(map(int, canonical[y : y + 2, x : x + 2].ravel()))
            for node in corners:
                adjacency[node].update(corners)
    for row, node in enumerate(T.nodes):
        support = T.neighbors[row]
        assert set(support[support >= 0]) == adjacency[int(node)]
    # Transpose really is the adjoint of the padded sparse gather, including seams.
    rng = np.random.default_rng(82)
    gradients = rng.normal(size=T.neighbors.shape)
    gradients[T.neighbors < 0] = 0
    velocity = rng.normal(size=nx * ny)
    stress = rng.normal(size=len(T.nodes))
    forward = np.sum(gradients * velocity[np.maximum(T.neighbors, 0)], axis=1)
    transpose = np.zeros(nx * ny)
    for node in range(nx * ny):
        rows, slots = T.transpose_rows[node], T.transpose_slots[node]
        good = rows >= 0
        transpose[node] = np.sum(gradients[rows[good], slots[good]] * stress[rows[good]])
    assert_allclose(stress @ forward, velocity @ transpose, rtol=1e-14, atol=1e-13)
    for color in np.unique(T.colors):
        support = T.neighbors[T.colors == color].ravel()
        support = support[support >= 0]
        assert len(support) == len(np.unique(support))


def test_local_coordinates_do_not_lose_subcell_motion_on_translated_grid():
    origin = np.array([2**26, -(2**27)], dtype=float)
    spacing = np.array([0.25, 0.5])
    local = np.array([[0.125, 0.75], [0.5, -0.25], [1.0, 1.125]])
    positions = origin + local
    cells, fraction = split_positions(positions, origin, spacing)
    assert_allclose(join_positions(cells, fraction, origin, spacing), positions, rtol=0, atol=0)
    small = np.array([0.0001, -0.0001]) / spacing
    advanced = join_positions(cells, fraction + small.astype(np.float32), origin, spacing)
    assert_allclose(advanced - positions, [0.0001, -0.0001] * np.ones((3, 1)), rtol=0, atol=4e-8)
    # Ordinary world-space float32 storage would have lost all of this motion.
    assert np.all(positions.astype(np.float32) == (positions + small * spacing).astype(np.float32))


def test_fraction_rounding_carries_to_next_integer_cell():
    cells, fraction = split_positions([[1 - 1e-9, -1e-9]], [0, 0], [1, 1])
    assert np.all((fraction >= 0) & (fraction < 1))
    assert cells.tolist() == [[1, 0]]


@pytest.mark.parametrize(
    "points,origin,spacing",
    [([[float("nan"), 0]], [0, 0], [1, 1]), ([[0, 0]], [0, 0], [0, 1]), ([[2**32, 0]], [0, 0], [1, 1])],
)
def test_coordinate_conversion_rejects_invalid_or_overflowing_state(points, origin, spacing):
    with pytest.raises(ValueError):
        split_positions(points, origin, spacing)


@pytest.mark.parametrize("count", [0, -1, 0.5, True, 2**32])
def test_bad_dispatch_counts_are_rejected(count):
    with pytest.raises(ValueError):
        check_count(count)


def test_scan_maximum_avoids_ceil_group_integer_overflow():
    assert check_count(MAX_COUNT) == MAX_COUNT
    assert MAX_COUNT + 255 < 2**32


def test_cpu_validation_report_is_not_a_gpu_pass(tmp_path):
    output = tmp_path / "cpu.json"
    assert main(["--cpu-only", "--stage", "all", "--output", str(output)]) == 0
    report = json.loads(output.read_text())
    assert report["status"] == "cpu_checks_passed_gpu_unvalidated"
    assert report["cpu_only"] is True
    assert report["stages"] == {}
    assert "shaders/primitives.slang" in report["source_sha256"]


def test_unavailable_gpu_failure_still_writes_report(tmp_path, monkeypatch):
    import mpm_gpu.validate as module

    def unavailable(**kwargs):
        raise GPUUnavailable("no test GPU")

    monkeypatch.setattr(module, "MetalRuntime", unavailable)
    output = tmp_path / "failure.json"
    assert main(["--output", str(output)]) == 1
    report = json.loads(output.read_text())
    assert report["status"] == "unavailable"
    assert "no test GPU" in report["error"]
    assert "traceback" in report


def test_platform_gate_does_not_try_importing_native_code_on_intel(monkeypatch):
    import mpm_gpu.runtime as runtime

    monkeypatch.setattr(runtime, "environment", lambda: dict(system="Darwin", machine="x86_64", macos="26.0"))
    with pytest.raises(GPUUnavailable, match="native Apple silicon"):
        check_platform()


def test_cpu_self_checks():
    assert cpu_checks()["status"] == "passed"


def test_platform_gate_detects_xcode_forwarder_without_real_toolchain(monkeypatch):
    import mpm_gpu.runtime as runtime

    monkeypatch.setattr(runtime, "environment", lambda: dict(system="Darwin", machine="arm64", macos="26.0"))
    calls = iter(
        [
            SimpleNamespace(returncode=0, stdout="/Xcode/usr/bin/metal\n", stderr=""),
            SimpleNamespace(returncode=1, stdout="", stderr="missing Metal Toolchain"),
        ]
    )
    monkeypatch.setattr(runtime.subprocess, "run", lambda *args, **kwargs: next(calls))
    with pytest.raises(GPUUnavailable, match="present but cannot run"):
        check_platform()
