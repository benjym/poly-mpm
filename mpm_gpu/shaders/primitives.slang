// All storage is scalar (4-byte stride). No FP64 or floating-point atomics.
// Every cross-workgroup dependency uses another dispatch and an API barrier.

[shader("compute")]
[numthreads(256, 1, 1)]
void fill_u32(uint3 tid : SV_DispatchThreadID,
              uniform uint count, uniform uint value, RWStructuredBuffer<uint> output)
{
    if (tid.x < count) output[tid.x] = value;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void histogram(uint3 tid : SV_DispatchThreadID, uniform uint count,
               uniform uint bin_count, StructuredBuffer<uint> keys,
               RWStructuredBuffer<uint> counts, RWStructuredBuffer<uint> status)
{
    if (tid.x >= count) return;
    uint key = keys[tid.x];
    uint ignored;
    if (key < bin_count) InterlockedAdd(counts[key], 1, ignored);
    else InterlockedOr(status[0], 1, ignored);
}

groupshared uint scan_values[256];

[shader("compute")]
[numthreads(256, 1, 1)]
void scan_tiles(uint3 tid : SV_DispatchThreadID, uint3 lane : SV_GroupThreadID,
                uint3 group : SV_GroupID, uniform uint count,
                StructuredBuffer<uint> input, RWStructuredBuffer<uint> output,
                RWStructuredBuffer<uint> totals)
{
    uint original = tid.x < count ? input[tid.x] : 0;
    scan_values[lane.x] = original;
    GroupMemoryBarrierWithGroupSync();
    // Hillis-Steele inclusive scan, with separate read/write barriers.
    for (uint offset = 1; offset < 256; offset *= 2)
    {
        uint addend = lane.x >= offset ? scan_values[lane.x - offset] : 0;
        GroupMemoryBarrierWithGroupSync();
        scan_values[lane.x] += addend;
        GroupMemoryBarrierWithGroupSync();
    }
    if (tid.x < count) output[tid.x] = scan_values[lane.x] - original;
    if (lane.x == 255) totals[group.x] = scan_values[255];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void scan_offsets(uint3 tid : SV_DispatchThreadID, uniform uint count,
                  StructuredBuffer<uint> offsets, RWStructuredBuffer<uint> output)
{
    if (tid.x < count) output[tid.x] += offsets[tid.x / 256];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void scatter_bins(uint3 tid : SV_DispatchThreadID, uniform uint count,
                  uniform uint bin_count, StructuredBuffer<uint> keys,
                  StructuredBuffer<uint> offsets, RWStructuredBuffer<uint> cursors,
                  RWStructuredBuffer<uint> permutation)
{
    if (tid.x >= count || keys[tid.x] >= bin_count) return;
    uint slot;
    InterlockedAdd(cursors[keys[tid.x]], 1, slot);
    permutation[offsets[keys[tid.x]] + slot] = tid.x;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void sort_bins(uint3 tid : SV_DispatchThreadID, uniform uint bin_count,
               StructuredBuffer<uint> offsets, StructuredBuffer<uint> counts,
               RWStructuredBuffer<uint> permutation)
{
    // The atomic scatter leaves each bin in scheduling order. Sorting every
    // bin's indices makes the order ascending, as a stable sort by key gives,
    // so floating-point gathers over a bin are reproducible run to run.
    // Insertion sort: bins hold a few particles; cost is quadratic in a bin.
    uint b = tid.x;
    if (b >= bin_count) return;
    uint begin = offsets[b], end = begin + counts[b];
    for (uint i = begin + 1; i < end; i++)
    {
        uint value = permutation[i];
        uint j = i;
        while (j > begin && permutation[j - 1] > value)
        {
            permutation[j] = permutation[j - 1];
            j--;
        }
        permutation[j] = value;
    }
}

groupshared float reduce_values[256];

[shader("compute")]
[numthreads(256, 1, 1)]
void reduce_tiles(uint3 tid : SV_DispatchThreadID, uint3 lane : SV_GroupThreadID,
                  uint3 group : SV_GroupID, uniform uint count, uniform uint operation,
                  StructuredBuffer<float> input, RWStructuredBuffer<float> output,
                  RWStructuredBuffer<uint> status)
{
    // op 0: sum; 1: max; 2: min. Identity is finite to simplify Metal lowering.
    float identity = operation == 0 ? 0.0 : (operation == 1 ? -3.402823466e38 : 3.402823466e38);
    float value = tid.x < count ? input[tid.x] : identity;
    if (tid.x < count && !isfinite(value))
    {
        uint ignored;
        InterlockedOr(status[0], 2, ignored);
        value = identity;
    }
    reduce_values[lane.x] = value;
    GroupMemoryBarrierWithGroupSync();
    for (uint offset = 128; offset > 0; offset /= 2)
    {
        if (lane.x < offset)
        {
            float a = reduce_values[lane.x], b = reduce_values[lane.x + offset];
            reduce_values[lane.x] = operation == 0 ? a + b : (operation == 1 ? max(a,b) : min(a,b));
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (lane.x == 0)
    {
        output[group.x] = reduce_values[0];
        if (!isfinite(reduce_values[0]))
        {
            uint ignored;
            InterlockedOr(status[0], 2, ignored);
        }
    }
}

[shader("compute")]
[numthreads(1, 1, 1)]
void make_dispatch_args(uniform uint count, StructuredBuffer<uint> control,
                        RWStructuredBuffer<uint> args)
{
    // A device control record decides whether optional work is dispatched.
    args[0] = control[0] != 0 ? (count + 255) / 256 : 0;
    args[1] = 1;
    args[2] = 1;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void indirect_marker(uint3 tid : SV_DispatchThreadID, uniform uint count,
                     uniform uint value, RWStructuredBuffer<uint> output)
{
    if (tid.x < count) output[tid.x] = value;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void guarded_marker(uint3 tid : SV_DispatchThreadID, uniform uint count,
                    StructuredBuffer<uint> control, uniform uint value,
                    RWStructuredBuffer<uint> output)
{
    // SlangPy 0.43.1 Metal has no indirect compute dispatch. The host launches
    // an upper bound; a device record sets the active extent.
    if (tid.x < min(count, control[0])) output[tid.x] = value;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void abi_probe(uint3 tid : SV_DispatchThreadID, uniform uint count,
               StructuredBuffer<uint> integers, StructuredBuffer<float> reals,
               RWStructuredBuffer<uint> integer_output, RWStructuredBuffer<float> real_output)
{
    if (tid.x < count)
    {
        integer_output[tid.x] = integers[tid.x] ^ 0xa5a5a5a5u;
        real_output[tid.x] = 2.0 * reals[tid.x] + 0.25;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void reduce_columns(uint3 lane : SV_GroupThreadID, uint3 group : SV_GroupID, uniform uint count,
                    uniform uint groups, StructuredBuffer<uint> operations,
                    StructuredBuffer<float> input, RWStructuredBuffer<float> output,
                    RWStructuredBuffer<uint> status)
{
    // Workgroup g reduces tile g % groups of column g / groups. Column-major
    // input and output, so each level feeds the next unchanged.
    uint column = group.x / groups, tile = group.x % groups;
    uint operation = operations[column];
    uint i = tile * 256 + lane.x;
    float identity = operation == 0 ? 0.0 : (operation == 1 ? -3.402823466e38 : 3.402823466e38);
    float value = i < count ? input[column * count + i] : identity;
    if (i < count && !isfinite(value))
    {
        uint ignored;
        InterlockedOr(status[0], 2, ignored);
        value = identity;
    }
    reduce_values[lane.x] = value;
    GroupMemoryBarrierWithGroupSync();
    for (uint offset = 128; offset > 0; offset /= 2)
    {
        if (lane.x < offset)
        {
            float a = reduce_values[lane.x], b = reduce_values[lane.x + offset];
            reduce_values[lane.x] = operation == 0 ? a + b : (operation == 1 ? max(a, b) : min(a, b));
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (lane.x == 0) output[column * groups + tile] = reduce_values[0];
}
