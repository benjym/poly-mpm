"""Resident scan, binning and reductions. No intermediate GPU readbacks."""

from dataclasses import dataclass

import numpy as np

TILE = 256
MAX_COUNT = 2**32 - TILE


def check_count(count, *, name="count"):
    if isinstance(count, bool) or int(count) != count or not 1 <= count <= MAX_COUNT:
        raise ValueError(f"{name} must be an integer in [1, {MAX_COUNT}]")
    return int(count)


@dataclass
class Bins:
    counts: object
    offsets: object
    cursors: object
    permutation: object
    status: object
    scratch: list


def exclusive_scan(runtime, encoder, source, count):
    """Hierarchical uint32 scan. Caller guarantees total fits uint32.

    Returns output and retained temporary resources, since a command buffer may
    outlive this call. At most ceil(log_256(count)) levels; zero host transfers.
    """
    count = check_count(count)
    groups = (count + TILE - 1) // TILE
    output = runtime.buffer(count=count)
    totals = runtime.buffer(count=groups)
    retained = [output, totals]
    runtime.dispatch(
        encoder, "scan_tiles", groups * TILE, count=count, input=source, output=output, totals=totals
    )
    if groups > 1:
        offsets, scratch = exclusive_scan(runtime, encoder, totals, groups)
        retained.extend(scratch)
        runtime.dispatch(encoder, "scan_offsets", count, count=count, offsets=offsets, output=output)
    return output, retained


def bins(runtime, encoder, keys, count, bin_count):
    """Integer histogram -> scan -> scatter, with explicit invalid-key status.

    Ordering within a bin is unspecified; each particle index appears once.
    Only integer atomics are used. Consumers gather by bin on the device.
    """
    count = check_count(count)
    bin_count = check_count(bin_count, name="bin_count")
    counts = runtime.buffer(count=bin_count)
    cursors = runtime.buffer(count=bin_count)
    status = runtime.buffer(count=1)
    permutation = runtime.buffer(count=count)
    for buffer, size in ((counts, bin_count), (cursors, bin_count), (status, 1)):
        runtime.dispatch(encoder, "fill_u32", size, count=size, value=0, output=buffer)
    # Make invalid entries recognizable instead of reading uninitialized memory.
    runtime.dispatch(encoder, "fill_u32", count, count=count, value=2**32 - 1, output=permutation)
    runtime.dispatch(
        encoder, "histogram", count, count=count, bin_count=bin_count, keys=keys, counts=counts, status=status
    )
    offsets, scratch = exclusive_scan(runtime, encoder, counts, bin_count)
    runtime.dispatch(
        encoder,
        "scatter_bins",
        count,
        count=count,
        bin_count=bin_count,
        keys=keys,
        offsets=offsets,
        cursors=cursors,
        permutation=permutation,
    )
    return Bins(counts, offsets, cursors, permutation, status, scratch)


def reduce(runtime, encoder, source, count, *, operation="sum", status=None):
    """Hierarchical FP32 tree reduction, including nonfinite/overflow status."""
    if operation not in ("sum", "max", "min"):
        raise ValueError("operation must be sum, max or min")
    count = check_count(count)
    retained = []
    if status is None:
        status = runtime.buffer(count=1)
        retained.append(status)
        runtime.dispatch(encoder, "fill_u32", 1, count=1, value=0, output=status)
    while True:
        groups = (count + TILE - 1) // TILE
        result = runtime.buffer(count=groups, dtype=np.float32)
        retained.append(result)
        runtime.dispatch(
            encoder,
            "reduce_tiles",
            groups * TILE,
            count=count,
            operation=("sum", "max", "min").index(operation),
            input=source,
            output=result,
            status=status,
        )
        if groups == 1:
            return result, status, retained
        source, count = result, groups


class Reduction:
    """Preallocated hierarchical FP32 reduction, so re-encoding allocates nothing.

    The result is element 0 of `result`. Nonfinite inputs set status bit 2 in
    whichever status buffer is passed to encode().
    """

    def __init__(self, runtime, count, *, operation="max"):
        if operation not in ("sum", "max", "min"):
            raise ValueError("operation must be sum, max or min")
        self.runtime, self.count = runtime, check_count(count)
        self.operation = ("sum", "max", "min").index(operation)
        self.levels = []
        while True:
            groups = (count + TILE - 1) // TILE
            self.levels.append((count, runtime.buffer(count=groups, dtype=np.float32)))
            if groups == 1:
                break
            count = groups
        self.result = self.levels[-1][1]

    def encode(self, encoder, source, status):
        for count, output in self.levels:
            groups = (count + TILE - 1) // TILE
            self.runtime.dispatch(
                encoder,
                "reduce_tiles",
                groups * TILE,
                count=count,
                operation=self.operation,
                input=source,
                output=output,
                status=status,
            )
            source = output
        return self.result


class Scan:
    """Preallocated hierarchical exclusive uint32 scan (see exclusive_scan)."""

    def __init__(self, runtime, count):
        self.runtime, self.count = runtime, check_count(count)
        groups = (self.count + TILE - 1) // TILE
        self.output = runtime.buffer(count=self.count)
        self.totals = runtime.buffer(count=groups)
        self.inner = Scan(runtime, groups) if groups > 1 else None

    def encode(self, encoder, source):
        groups = (self.count + TILE - 1) // TILE
        self.runtime.dispatch(
            encoder,
            "scan_tiles",
            groups * TILE,
            count=self.count,
            input=source,
            output=self.output,
            totals=self.totals,
        )
        if self.inner is not None:
            offsets = self.inner.encode(encoder, self.totals)
            self.runtime.dispatch(
                encoder, "scan_offsets", self.count, count=self.count, offsets=offsets, output=self.output
            )
        return self.output


class Binner:
    """Preallocated histogram -> scan -> scatter with the same contract as bins()."""

    def __init__(self, runtime, count, bin_count):
        self.runtime = runtime
        self.count = check_count(count)
        self.bin_count = check_count(bin_count, name="bin_count")
        self.counts = runtime.buffer(count=self.bin_count)
        self.cursors = runtime.buffer(count=self.bin_count)
        self.status = runtime.buffer(count=1)
        self.permutation = runtime.buffer(count=self.count)
        self.scan = Scan(runtime, self.bin_count)
        self.offsets = self.scan.output

    def encode(self, encoder, keys):
        rt, n, k = self.runtime, self.count, self.bin_count
        for buffer, size in ((self.counts, k), (self.cursors, k), (self.status, 1)):
            rt.dispatch(encoder, "fill_u32", size, count=size, value=0, output=buffer)
        rt.dispatch(encoder, "fill_u32", n, count=n, value=2**32 - 1, output=self.permutation)
        rt.dispatch(
            encoder, "histogram", n, count=n, bin_count=k, keys=keys, counts=self.counts, status=self.status
        )
        self.scan.encode(encoder, self.counts)
        rt.dispatch(
            encoder,
            "scatter_bins",
            n,
            count=n,
            bin_count=k,
            keys=keys,
            offsets=self.offsets,
            cursors=self.cursors,
            permutation=self.permutation,
        )
        return self


SUM, MAX, MIN = 0, 1, 2


class ColumnReduction:
    """Preallocated reduction of several columns at once, each with its own operation.

    Input is column-major: column c occupies elements [c * count, (c + 1) * count).
    The results are elements 0..columns-1 of `result`.
    """

    def __init__(self, runtime, count, operations):
        self.runtime, self.count = runtime, check_count(count)
        self.columns = len(operations)
        self.operations = runtime.buffer(np.asarray(operations, dtype=np.uint32))
        self.levels = []
        while True:
            groups = (count + TILE - 1) // TILE
            self.levels.append((count, groups, runtime.buffer(count=groups * self.columns, dtype=np.float32)))
            if groups == 1:
                break
            count = groups
        self.result = self.levels[-1][2]

    def encode(self, encoder, source, status):
        for count, groups, output in self.levels:
            self.runtime.dispatch(
                encoder,
                "reduce_columns",
                groups * TILE * self.columns,
                count=count,
                groups=groups,
                operations=self.operations,
                input=source,
                output=output,
                status=status,
            )
            source = output
        return self.result
