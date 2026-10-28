"""Small SlangPy Metal adapter with explicit synchronization/accounting."""

from contextlib import ExitStack
from dataclasses import dataclass, asdict
from importlib import import_module, metadata, resources
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter

import numpy as np

from . import GPUError, GPUUnavailable, SLANGPY_MINIMUM


def environment():
    return dict(
        python=sys.version,
        system=platform.system(),
        machine=platform.machine(),
        macos=platform.mac_ver()[0],
        platform=platform.platform(),
        slangpy_minimum=".".join(map(str, SLANGPY_MINIMUM)),
    )


def check_platform():
    info = environment()
    if info["system"] != "Darwin" or info["machine"] != "arm64":
        raise GPUUnavailable(
            "Metal jobs require native Apple silicon Python on macOS 26+. "
            "This host can run python -m mpm_gpu.validate --cpu-only."
        )
    if not info["macos"] or int(info["macos"].split(".")[0]) < 26:
        raise GPUUnavailable("The pinned SlangPy wheel requires macOS 26 or later")
    try:
        result = subprocess.run(
            ["xcrun", "--find", "metal"], capture_output=True, text=True, timeout=15, check=False
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise GPUUnavailable("Cannot locate the Metal compiler with xcrun --find metal") from exc
    if result.returncode or not result.stdout.strip():
        raise GPUUnavailable(
            "Metal compiler unavailable. Install Xcode and its Metal Toolchain component, "
            "then check xcrun --find metal. " + result.stderr.strip()
        )
    info["metal_compiler"] = result.stdout.strip()
    # Xcode can ship a forwarding executable even when its toolchain component
    # is absent. Finding that executable alone is not a compiler check.
    try:
        version = subprocess.run(
            ["xcrun", "metal", "--version"], capture_output=True, text=True, timeout=15, check=False
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise GPUUnavailable("Cannot execute the Metal compiler") from exc
    if version.returncode:
        raise GPUUnavailable(
            "The Metal compiler executable is present but cannot run. "
            "Install Xcode's Metal Toolchain component. " + version.stderr.strip()
        )
    info["metal_compiler_version"] = version.stdout.strip()
    return info


@dataclass
class Accounting:
    uploaded_bytes: int = 0
    readback_bytes: int = 0
    readback_seconds: float = 0.0
    wait_seconds: float = 0.0
    allocated_bytes: int = 0
    submissions: int = 0
    dispatches: int = 0
    compile_seconds: float = 0.0


class MetalRuntime:
    """Own device, package-resource lifetime and compiled kernels.

    Buffers are scalar, tightly packed uint32/float32 arrays. Vector layout and
    hidden NumPy synchronization are deliberately avoided at this boundary.
    """

    def __init__(self, *, debug=False):
        self.environment = check_platform()
        try:
            installed = metadata.version("slangpy")
        except metadata.PackageNotFoundError as exc:
            raise GPUUnavailable("Install the optional GPU dependency: pip install 'poly-mpm[gpu]'") from exc
        if tuple(int(part) for part in installed.split(".")[:3]) < SLANGPY_MINIMUM:
            raise GPUUnavailable(f"slangpy {installed} is older than the validated minimum")
        self._resources = ExitStack()
        self.accounting = Accounting()
        self._kernels = {}
        self._closed = False
        try:
            self.spy = import_module("slangpy")
            self.shader_dir = Path(
                self._resources.enter_context(resources.as_file(resources.files("mpm_gpu.shaders")))
            )
            started = perf_counter()
            self.device = self.spy.Device(
                type=self.spy.DeviceType.metal,
                enable_debug_layers=debug,
                compiler_options={"include_paths": [str(self.shader_dir)]},
            )
            self.startup_seconds = perf_counter() - started
            self.environment.update(
                slangpy=installed,
                adapter=self.device.info.adapter_name,
                api=self.device.info.api_name,
                timestamp_frequency=self.device.info.timestamp_frequency,
            )
            if self.device.info.type != self.spy.DeviceType.metal:
                raise GPUUnavailable("SlangPy did not create a Metal device")
        except Exception as exc:
            self._resources.close()
            raise GPUUnavailable(f"Metal initialization failed: {exc}") from exc

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        if not self._closed:
            self.device.wait()
            self._resources.close()
            self._closed = True

    def kernel(self, entry, module="primitives.slang"):
        key = module, entry
        if key not in self._kernels:
            started = perf_counter()
            try:
                program = self.device.load_program(module_name=module, entry_point_names=[entry])
                self._kernels[key] = self.device.create_compute_kernel(program=program)
            except Exception as exc:
                raise GPUError(f"Metal shader compilation failed for {module}:{entry}: {exc}") from exc
            self.accounting.compile_seconds += perf_counter() - started
        return self._kernels[key]

    def buffer(self, data=None, *, count=None, dtype=np.uint32, indirect=False, label=""):
        dtype = np.dtype(dtype)
        if dtype not in (np.dtype("uint32"), np.dtype("float32")):
            raise ValueError("GPU buffers require uint32 or float32 scalar storage")
        if data is not None:
            data = np.ascontiguousarray(data, dtype=dtype).ravel()
            count = data.size
        if count is None or int(count) != count or count < 1:
            raise ValueError("GPU buffers require a positive element count")
        usage = (
            self.spy.BufferUsage.shader_resource
            | self.spy.BufferUsage.unordered_access
            | self.spy.BufferUsage.copy_source
            | self.spy.BufferUsage.copy_destination
        )
        if indirect:
            usage |= self.spy.BufferUsage.indirect_argument
        buffer = self.device.create_buffer(
            element_count=int(count), struct_size=4, usage=usage, data=data, label=label
        )
        self.accounting.allocated_bytes += int(count) * 4
        if data is not None:
            self.accounting.uploaded_bytes += data.nbytes
        return buffer

    def dispatch(self, encoder, entry, threads, *, module="primitives.slang", **bindings):
        if threads < 1:
            raise ValueError("dispatch count must be positive")
        self.kernel(entry, module).dispatch(
            thread_count=[int(threads), 1, 1], command_encoder=encoder, **bindings
        )
        self.accounting.dispatches += 1
        # Dispatch boundaries are global synchronization points. No spin locks.
        encoder.global_barrier()

    def dispatch_indirect(self, encoder, entry, args, *, module="primitives.slang", **bindings):
        kernel = self.kernel(entry, module)
        with encoder.begin_compute_pass() as compute:
            root = compute.bind_pipeline(kernel.pipeline)
            cursor = self.spy.ShaderCursor(root)
            # Entry-point parameters occupy the entry-point scope in Slang.
            target = cursor.find_entry_point(0)
            for name, value in bindings.items():
                target[name] = value
            compute.dispatch_compute_indirect(self.spy.BufferOffsetPair(args, 0))
        self.accounting.dispatches += 1
        encoder.global_barrier()

    def submit(self, encoder):
        self.device.submit_command_buffer(encoder.finish())
        self.accounting.submissions += 1

    def wait(self):
        started = perf_counter()
        self.device.wait()
        self.accounting.wait_seconds += perf_counter() - started

    def read(self, buffer, dtype):
        started = perf_counter()
        # Buffer.to_numpy synchronizes its readback. Count the actual byte copy.
        raw = np.asarray(buffer.to_numpy())
        result = np.frombuffer(raw.tobytes(), dtype=dtype).copy()
        self.accounting.readback_bytes += raw.nbytes
        self.accounting.readback_seconds += perf_counter() - started
        return result

    def timestamps(self, pool):
        started = perf_counter()
        ticks = pool.get_results(0, 2)
        self.accounting.readback_bytes += 16
        self.accounting.readback_seconds += perf_counter() - started
        return ticks

    def metrics(self):
        return asdict(self.accounting)
