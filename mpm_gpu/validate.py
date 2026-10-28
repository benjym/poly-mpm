"""Hardware gates for the experimental Metal port; failures are saved as JSON."""

import argparse
import hashlib
from importlib import resources
import json
from pathlib import Path
from time import perf_counter
import traceback

import numpy as np

from . import GPUError
from .primitives import bins, exclusive_scan, reduce
from .runtime import MetalRuntime, environment
from .topology import q1_topology, split_positions, join_positions

STAGES = ("runtime", "mechanics", "solver", "runner")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def fingerprint():
    root = resources.files("mpm_gpu")
    names = [
        "__init__.py",
        "config.py",
        "runtime.py",
        "primitives.py",
        "topology.py",
        "mechanics.py",
        "mechanics_checks.py",
        "fixtures.py",
        "solver.py",
        "solver_checks.py",
        "runner.py",
        "runner_checks.py",
        "validate.py",
        "shaders/primitives.slang",
        "shaders/mechanics.slang",
        "shaders/solver.slang",
        "shaders/runner.slang",
    ]
    return {name: hashlib.sha256(root.joinpath(name).read_bytes()).hexdigest() for name in names}


def cpu_checks():
    """Packaging and static-reference checks, explicitly not Metal validation."""
    shader = resources.files("mpm_gpu.shaders").joinpath("primitives.slang").read_text()
    require("void scan_tiles" in shader and "void guarded_marker" in shader, "shader resources missing")
    mechanics = resources.files("mpm_gpu.shaders").joinpath("mechanics.slang").read_text()
    require(
        "void assemble_gradients" in mechanics and "void apply_bt" in mechanics, "mechanics shader missing"
    )
    solver = resources.files("mpm_gpu.shaders").joinpath("solver.slang").read_text()
    require(
        all(
            f"void {entry}" in solver
            for entry in ("sweep_color", "tsqr_tile", "anderson_solve", "gmres_givens")
        ),
        "solver shader missing",
    )
    runner = resources.files("mpm_gpu.shaders").joinpath("runner.slang").read_text()
    require("void grid_finish" in runner and "void particle_update" in runner, "runner shader missing")
    from types import SimpleNamespace
    from .config import defaults

    P = SimpleNamespace(backend="metal")
    defaults(P)
    require(P.granular_solver == "velocity" and P.granular_tolerance == 1e-4, "Metal defaults")
    P = SimpleNamespace(backend="metal", granular_tolerance=3e-5, granular_solver="hybrid")
    defaults(P)
    require(P.granular_tolerance == 3e-5 and P.granular_solver == "hybrid", "explicit input was overridden")
    cases = []
    for nx, ny, periodic in ((2, 2, False), (2, 4, True), (3, 5, True), (9, 7, True), (9, 7, False)):
        T = q1_topology(nx, ny, periodic=periodic)
        for color in np.unique(T.colors):
            support = T.neighbors[T.colors == color].ravel()
            support = support[support >= 0]
            require(len(np.unique(support)) == len(support), "coloring has a write conflict")
        for row, neighbors in enumerate(T.neighbors):
            valid = neighbors[neighbors >= 0]
            require(len(np.unique(valid)) == len(valid), "periodic alias duplicated")
            for slot, node in enumerate(neighbors):
                if node >= 0:
                    require(
                        np.any((T.transpose_rows[node] == row) & (T.transpose_slots[node] == slot)),
                        "transpose incidence missing",
                    )
        cases.append(dict(nx=nx, ny=ny, periodic=periodic, colors=int(T.colors.max()) + 1))
    points = np.array([[1e8 + 0.125, -2e8 + 0.75], [1e8 + 0.499999999, -2e8 - 0.01]])
    origin = [1e8, -2e8]
    cells, fraction = split_positions(points, origin, [0.25, 0.5])
    reconstructed = join_positions(cells, fraction, origin, [0.25, 0.5])
    require(np.allclose(reconstructed, points, rtol=0, atol=2e-8), "split coordinate reconstruction")
    require(np.all((fraction >= 0) & (fraction < 1)), "invalid local coordinates")
    return dict(status="passed", evidence="CPU packaging/configuration/topology checks only", topology=cases)


def _case(report, name, function):
    started = perf_counter()
    try:
        details = function() or {}
        result = dict(name=name, status="passed", seconds=perf_counter() - started, **details)
    except Exception as exc:
        result = dict(
            name=name,
            status="failed",
            seconds=perf_counter() - started,
            error=f"{type(exc).__name__}: {exc}",
            traceback=traceback.format_exc(),
        )
    report.append(result)
    print(f"{name}: {result['status']}", flush=True)
    return result


def abi_check(rt):
    integers = np.array([0, 1, 2**31, 2**32 - 1, 0x12345678], dtype=np.uint32)
    reals = np.array([-2.0, -0.0, 0.5, 16.0, 2**-20], dtype=np.float32)
    a, b = rt.buffer(integers), rt.buffer(reals, dtype=np.float32)
    c, d = rt.buffer(count=len(integers)), rt.buffer(count=len(reals), dtype=np.float32)
    encoder = rt.device.create_command_encoder()
    rt.dispatch(
        encoder,
        "abi_probe",
        len(reals),
        count=len(reals),
        integers=a,
        reals=b,
        integer_output=c,
        real_output=d,
    )
    rt.submit(encoder)
    require(np.array_equal(rt.read(c, np.uint32), integers ^ np.uint32(0xA5A5A5A5)), "uint32 ABI mismatch")
    require(np.array_equal(rt.read(d, np.float32), 2 * reals + np.float32(0.25)), "float32 ABI mismatch")
    return dict(elements=len(reals))


def scan_check(rt, size):
    values = np.random.default_rng(size).integers(0, 11, size, dtype=np.uint32)
    source = rt.buffer(values)
    encoder = rt.device.create_command_encoder()
    output, scratch = exclusive_scan(rt, encoder, source, size)
    rt.submit(encoder)
    expected = np.cumsum(values, dtype=np.uint64) - values
    require(np.array_equal(rt.read(output, np.uint32), expected), f"exclusive scan mismatch at length {size}")
    return dict(elements=size, scratch_buffers=len(scratch))


def bins_check(rt, size, bin_count, mode):
    rng = np.random.default_rng(size + bin_count)
    keys = rng.integers(0, bin_count, size, dtype=np.uint32)
    if mode == "one_bin":
        keys[:] = bin_count - 1
    elif mode == "sparse":
        keys[:] = np.where(np.arange(size) % 2, 0, bin_count - 1)
    elif mode == "invalid":
        keys[-1] = bin_count
    device_keys = rt.buffer(keys)
    encoder = rt.device.create_command_encoder()
    result = bins(rt, encoder, device_keys, size, bin_count)
    rt.submit(encoder)
    status = int(rt.read(result.status, np.uint32)[0])
    valid = keys < bin_count
    require(status == (0 if np.all(valid) else 1), "invalid-key status mismatch")
    counts = rt.read(result.counts, np.uint32)
    offsets = rt.read(result.offsets, np.uint32)
    permutation = rt.read(result.permutation, np.uint32)
    expected = np.bincount(keys[valid].astype(int), minlength=bin_count)
    require(np.array_equal(counts, expected), "integer histogram mismatch")
    require(np.array_equal(offsets, np.cumsum(expected) - expected), "bin offsets mismatch")
    require(np.array_equal(rt.read(result.cursors, np.uint32), counts), "scatter counts mismatch")
    used = permutation[: int(counts.sum(dtype=np.uint64))]
    require(np.array_equal(np.sort(used), np.flatnonzero(valid)), "scatter lost/duplicated particle indices")
    require(
        np.array_equal(keys[used], np.repeat(np.arange(bin_count), counts.astype(int))),
        "scatter bin membership mismatch",
    )
    return dict(elements=size, bin_count=bin_count, distribution=mode)


def reduction_check(rt, operation, *, bad=None):
    # Exact integers permit an exact FP32 sum oracle below 2**24.
    values = np.random.default_rng(99).integers(-10, 11, 100003).astype(np.float32)
    if bad is not None:
        values[11] = bad
    source = rt.buffer(values, dtype=np.float32)
    encoder = rt.device.create_command_encoder()
    output, status, scratch = reduce(rt, encoder, source, len(values), operation=operation)
    rt.submit(encoder)
    flags = int(rt.read(status, np.uint32)[0])
    if bad is not None:
        require(flags & 2, "nonfinite reduction was not reported")
    else:
        require(flags == 0, "unexpected nonfinite reduction status")
        actual = rt.read(output, np.float32)[0]
        expected = getattr(np, operation)(values.astype(np.float64))
        require(actual == expected, f"{operation} reduction mismatch")
    return dict(elements=len(values), operation=operation, expected_nonfinite=bad is not None)


def guarded_check(rt, active):
    n = 1003
    control = rt.buffer(np.array([active], dtype=np.uint32))
    output = rt.buffer(count=n)
    encoder = rt.device.create_command_encoder()
    rt.dispatch(encoder, "fill_u32", n, count=n, value=7, output=output)
    rt.dispatch(encoder, "guarded_marker", n, count=n, control=control, value=42, output=output)
    rt.submit(encoder)
    expected = np.where(np.arange(n) < active, 42, 7)
    require(np.array_equal(rt.read(output, np.uint32), expected), "device-guarded dispatch mismatch")
    return dict(active=active, elements=n)


def barrier_check(rt):
    n = 100003
    buffer = rt.buffer(count=n)
    encoder = rt.device.create_command_encoder()
    for value in range(8):
        rt.dispatch(encoder, "fill_u32", n, count=n, value=value, output=buffer)
    started = perf_counter()
    rt.submit(encoder)
    rt.wait()
    elapsed = perf_counter() - started
    require(np.all(rt.read(buffer, np.uint32) == 7), "dispatch/barrier sequence mismatch")
    return dict(dispatches=8, elements=n, wall_seconds=elapsed)


def indirect_probe(rt):
    """Capability record. An unsupported dispatch may silently do nothing."""
    n = 1003
    control = rt.buffer(np.array([1], dtype=np.uint32))
    args = rt.buffer(count=3, indirect=True)
    output = rt.buffer(count=n)
    encoder = rt.device.create_command_encoder()
    rt.dispatch(encoder, "fill_u32", n, count=n, value=7, output=output)
    rt.dispatch(encoder, "make_dispatch_args", 1, count=n, control=control, args=args)
    rt.dispatch_indirect(encoder, "indirect_marker", args, count=n, value=42, output=output)
    rt.submit(encoder)
    values = rt.read(output, np.uint32)
    return dict(supported=bool(np.all(values == 42)), elements_written=int(np.sum(values == 42)))


def timestamp_probe(rt):
    pool = rt.device.create_query_pool(type=rt.spy.QueryType.timestamp, count=2)
    buffer = rt.buffer(count=1024)
    encoder = rt.device.create_command_encoder()
    encoder.write_timestamp(pool, 0)
    rt.dispatch(encoder, "fill_u32", 1024, count=1024, value=1, output=buffer)
    encoder.write_timestamp(pool, 1)
    rt.submit(encoder)
    rt.wait()
    frequency = rt.device.info.timestamp_frequency
    try:
        ticks = rt.timestamps(pool)
    except Exception as exc:
        return dict(supported=False, timestamp_frequency=frequency, error=f"{type(exc).__name__}: {exc}")
    return dict(supported=bool(frequency > 0 and ticks[1] > ticks[0]), timestamp_frequency=frequency)


def capabilities(rt):
    """Backend features recorded but not required; mechanics must not use unsupported ones."""
    result = {}
    for name, probe in (("indirect_dispatch", indirect_probe), ("gpu_timestamps", timestamp_probe)):
        try:
            result[name] = probe(rt)
        except Exception as exc:
            result[name] = dict(supported=False, error=f"{type(exc).__name__}: {exc}")
    return result


def primitive_benchmark(rt, size):
    """Paired primitive timings; not an MPM speedup or solver benchmark."""
    bin_count = max(1, size // 4)
    keys = np.random.default_rng(size).integers(0, bin_count, size, dtype=np.uint32)
    source = rt.buffer(keys)
    result = None

    def gpu():
        # Metal command buffers are single-use: encode (and allocate) each sample
        # outside the timed region. SlangPy 0.43.1 Metal has no timestamp
        # queries, so host wall time covers submission, execution and wait.
        nonlocal result
        started = perf_counter()
        encoder = rt.device.create_command_encoder()
        result = bins(rt, encoder, source, size, bin_count)
        command = encoder.finish()
        encode = perf_counter() - started
        started = perf_counter()
        rt.device.submit_command_buffer(command)
        rt.accounting.submissions += 1
        rt.wait()
        return dict(wall_seconds=perf_counter() - started, encode_seconds=encode)

    def cpu():
        started = perf_counter()
        counts = np.bincount(keys.astype(int), minlength=bin_count)
        offsets = np.cumsum(counts) - counts
        permutation = np.argsort(keys, kind="stable")
        require(len(offsets) == bin_count and len(permutation) == size, "CPU binning shape")
        return dict(wall_seconds=perf_counter() - started)

    gpu()
    cpu()
    samples = []
    for pair in range(5):
        if pair % 2:
            cpu_result, gpu_result = cpu(), gpu()
        else:
            gpu_result, cpu_result = gpu(), cpu()
        samples.append(
            dict(pair=pair, order="cpu-gpu" if pair % 2 else "gpu-cpu", cpu=cpu_result, gpu=gpu_result)
        )
    require(int(rt.read(result.status, np.uint32)[0]) == 0, "benchmark binning failed")
    require(
        np.array_equal(rt.read(result.counts, np.uint32), np.bincount(keys.astype(int), minlength=bin_count)),
        "benchmark output mismatch",
    )
    return dict(
        elements=size,
        bin_count=bin_count,
        samples=samples,
        scope="Resident integer binning only. NumPy stable sort and GPU unordered scatter differ within bins. Not an MPM speedup.",
        warmup_pairs=1,
        setup_and_upload_in_timing=False,
        allocation_bytes_total=rt.accounting.allocated_bytes,
    )


def validate_runtime(rt, *, benchmark=False):
    results = []
    checks = [("scalar_abi", lambda: abi_check(rt))]
    checks += [(f"scan_{n}", lambda n=n: scan_check(rt, n)) for n in (1, 255, 256, 257, 65537, 100003)]
    checks += [
        (f"bins_{mode}", lambda mode=mode: bins_check(rt, 10003, 513, mode))
        for mode in ("random", "one_bin", "sparse", "invalid")
    ]
    checks += [(f"reduction_{op}", lambda op=op: reduction_check(rt, op)) for op in ("sum", "min", "max")]
    checks += [
        ("nonfinite_nan", lambda: reduction_check(rt, "sum", bad=np.nan)),
        ("nonfinite_inf", lambda: reduction_check(rt, "max", bad=np.inf)),
    ]
    checks += [
        ("guarded_all", lambda: guarded_check(rt, 1003)),
        ("guarded_none", lambda: guarded_check(rt, 0)),
        ("guarded_partial", lambda: guarded_check(rt, 517)),
        ("barriers", lambda: barrier_check(rt)),
    ]
    for name, check in checks:
        _case(results, name, check)
    if benchmark and all(r["status"] == "passed" for r in results):
        for size in (10000, 100000):
            _case(results, f"binning_benchmark_{size}", lambda size=size: primitive_benchmark(rt, size))
    return dict(
        status="passed" if all(r["status"] == "passed" for r in results) else "failed",
        cases=results,
        capabilities=capabilities(rt),
    )


def _validate_cases(module, rt):
    results = []
    for name, check in module.cases(rt):
        _case(results, name, check)
    if any(r["status"] != "passed" for r in results):
        status = "failed"
    else:
        status = "partial" if module.PENDING else "passed"
    return dict(status=status, cases=results, pending=list(module.PENDING))


def validate_mechanics(rt):
    from . import mechanics_checks

    return _validate_cases(mechanics_checks, rt)


def validate_solver(rt):
    from . import solver_checks

    return _validate_cases(solver_checks, rt)


def validate_runner(rt):
    from . import runner, runner_checks

    # MPM.run and MPM.step reach the device through the runner's shared runtime.
    previous = runner._runtime
    runner._runtime = rt
    try:
        return _validate_cases(runner_checks, rt)
    finally:
        runner._runtime = previous
        runner._runners.clear()


def main(argv=None):
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--stage", choices=(*STAGES, "all"), default="all")
    parser.add_argument(
        "--cpu-only",
        action="store_true",
        help="Check packaging/configuration without importing SlangPy; never certifies Metal",
    )
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--debug", action="store_true", help="Enable graphics API validation layers")
    parser.add_argument("--output", type=Path, default=Path("metal-validation.json"))
    args = parser.parse_args(argv)
    report = dict(
        schema_version=1,
        environment=environment(),
        source_sha256={},
        requested_stage=args.stage,
        cpu_only=args.cpu_only,
        status="running",
        stages={},
    )
    code = 1
    try:
        report["source_sha256"] = fingerprint()
        report["cpu_checks"] = cpu_checks()
        if args.cpu_only:
            report["status"] = "cpu_checks_passed_gpu_unvalidated"
            code = 0
        else:
            with MetalRuntime(debug=args.debug) as rt:
                report["environment"] = rt.environment
                report["startup_seconds"] = rt.startup_seconds
                report["stages"]["runtime"] = validate_runtime(rt, benchmark=args.benchmark)
                # A stage runs only after every preceding stage has passed.
                last = len(STAGES) if args.stage == "all" else STAGES.index(args.stage) + 1
                requested = STAGES[1:last]
                previous = report["stages"]["runtime"]["status"]
                for stage, validate in (
                    ("mechanics", validate_mechanics),
                    ("solver", validate_solver),
                    ("runner", validate_runner),
                ):
                    if stage in requested and previous == "passed":
                        report["stages"][stage] = validate(rt)
                        previous = report["stages"][stage]["status"]
                report["accounting"] = rt.metrics()
            statuses = [stage["status"] for stage in report["stages"].values()]
            if "failed" in statuses:
                report["status"] = "failed"
            elif requested:
                for stage in requested:
                    report["stages"].setdefault(
                        stage,
                        dict(
                            status="not_implemented",
                            reason="Dependent on the preceding Apple silicon hardware validation gate; see docs/gpu.md",
                        ),
                    )
                if all(report["stages"][stage]["status"] == "passed" for stage in requested):
                    report["status"] = "passed"
                    code = 0
                else:
                    report["status"] = "runtime_passed_later_stages_pending"
                    code = 3
            else:
                report["status"] = "passed"
                code = 0
    except Exception as exc:
        report["status"] = "unavailable" if isinstance(exc, GPUError) else "failed"
        report["error"] = f"{type(exc).__name__}: {exc}"
        report["traceback"] = traceback.format_exc()
    finally:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        print(f"{report['status']}: {args.output}", flush=True)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
