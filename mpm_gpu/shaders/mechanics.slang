// Compact Q1 nodal-stress operator: one integrated gradient 2-vector per
// stress-node/velocity-node pair. Row-owned assembly and B.u, node-owned
// B.T z. No atomics in floating point; every write has a single owner.

static const uint INVALID = 0xffffffffu;
static const float HALF_ROOT2 = 0.70710678118654752;

uint canonical(uint column, uint row, uint nx, uint periodic)
{
    return (periodic != 0 && column == nx - 1 ? 0 : column) + nx * row;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void particle_cells(uint3 tid : SV_DispatchThreadID, uniform uint count,
                    uniform uint nx, uniform uint ny, uniform uint periodic,
                    StructuredBuffer<uint> cell_x, StructuredBuffer<uint> cell_y,
                    StructuredBuffer<float> fraction_x, StructuredBuffer<float> fraction_y,
                    RWStructuredBuffer<uint> keys, RWStructuredBuffer<float> local_x,
                    RWStructuredBuffer<float> local_y, RWStructuredBuffer<uint> status)
{
    if (tid.x >= count) return;
    int columns = int(nx) - 1;
    int rows = int(ny) - 1;
    int cx = asint(cell_x[tid.x]);
    int cy = asint(cell_y[tid.x]);
    float fx = fraction_x[tid.x];
    float fy = fraction_y[tid.x];
    if (periodic != 0) cx = ((cx % columns) + columns) % columns;
    // The closed upper edge belongs to the last cell, as in interpolation.stencil.
    else if (cx == columns && fx == 0.0) { cx = columns - 1; fx = 1.0; }
    if (cy == rows && fy == 0.0) { cy = rows - 1; fy = 1.0; }
    bool valid = cx >= 0 && cx < columns && cy >= 0 && cy < rows
                 && fx >= 0.0 && fx <= 1.0 && fy >= 0.0 && fy <= 1.0;
    if (!valid)
    {
        uint ignored;
        InterlockedOr(status[0], 1, ignored);
        keys[tid.x] = uint(columns * rows);
        return;
    }
    keys[tid.x] = uint(cx + columns * cy);
    local_x[tid.x] = fx;
    local_y[tid.x] = fy;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void assemble_gradients(uint3 tid : SV_DispatchThreadID, uniform uint row_count,
                        uniform uint nx, uniform uint ny, uniform uint periodic,
                        uniform float dx, uniform float dy,
                        StructuredBuffer<uint> row_nodes, StructuredBuffer<uint> neighbors,
                        StructuredBuffer<uint> cell_offsets, StructuredBuffer<uint> cell_counts,
                        StructuredBuffer<uint> permutation, StructuredBuffer<float> local_x,
                        StructuredBuffer<float> local_y, StructuredBuffer<float> volume,
                        RWStructuredBuffer<float> gradients, RWStructuredBuffer<uint> status)
{
    uint a = tid.x;
    if (a >= row_count) return;
    uint node = row_nodes[a];
    int x = int(node % nx);
    int y = int(node / nx);
    int columns = int(nx) - 1;
    float2 sum[9];
    for (uint s = 0; s < 9; s++) sum[s] = float2(0.0, 0.0);
    int previous = -1;
    for (int ox = -1; ox <= 0; ox++)
    {
        int cx = x + ox;
        if (periodic != 0) cx = ((cx % columns) + columns) % columns;
        else if (cx < 0 || cx >= columns) continue;
        // A one-cell periodic strip sees the same column from both sides.
        if (cx == previous) continue;
        previous = cx;
        for (int oy = -1; oy <= 0; oy++)
        {
            int cy = y + oy;
            if (cy < 0 || cy >= int(ny) - 1) continue;
            // Corner order matches interpolation.stencil: (0,0), (1,0), (1,1), (0,1).
            uint corners[4] = {
                canonical(uint(cx), uint(cy), nx, periodic),
                canonical(uint(cx + 1), uint(cy), nx, periodic),
                canonical(uint(cx + 1), uint(cy + 1), nx, periodic),
                canonical(uint(cx), uint(cy + 1), nx, periodic),
            };
            uint slots[4];
            float stress[4];
            for (uint c = 0; c < 4; c++)
            {
                stress[c] = corners[c] == node ? 1.0 : 0.0;
                slots[c] = INVALID;
                for (uint s = 0; s < 9; s++)
                    if (neighbors[9 * a + s] == corners[c]) slots[c] = s;
                if (slots[c] == INVALID)
                {
                    uint ignored;
                    InterlockedOr(status[0], 4, ignored);
                    return;
                }
            }
            uint cell = uint(cx + columns * cy);
            uint begin = cell_offsets[cell];
            uint end = begin + cell_counts[cell];
            for (uint k = begin; k < end; k++)
            {
                uint p = permutation[k];
                float fx = local_x[p];
                float fy = local_y[p];
                float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
                float weight = 0.0;
                for (uint c = 0; c < 4; c++) weight += stress[c] * w[c];
                if (weight == 0.0) continue;
                float scale = volume[p] * weight;
                float2 g[4] = {
                    float2(-(1 - fy) / dx, -(1 - fx) / dy),
                    float2((1 - fy) / dx, -fx / dy),
                    float2(fy / dx, fx / dy),
                    float2(-fy / dx, (1 - fx) / dy),
                };
                for (uint c = 0; c < 4; c++) sum[slots[c]] += scale * g[c];
            }
        }
    }
    for (uint s = 0; s < 9; s++)
    {
        gradients[2 * (9 * a + s)] = sum[s].x;
        gradients[2 * (9 * a + s) + 1] = sum[s].y;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void apply_b(uint3 tid : SV_DispatchThreadID, uniform uint row_count,
             StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
             StructuredBuffer<float> velocity, RWStructuredBuffer<float> output)
{
    // Intrinsic Frobenius-orthonormal basis: (I, diag(1,-1), offdiag)/sqrt(2).
    uint a = tid.x;
    if (a >= row_count) return;
    float3 w = float3(0.0, 0.0, 0.0);
    for (uint s = 0; s < 9; s++)
    {
        uint j = neighbors[9 * a + s];
        if (j == INVALID) continue;
        float gx = gradients[2 * (9 * a + s)];
        float gy = gradients[2 * (9 * a + s) + 1];
        float ux = velocity[2 * j];
        float uy = velocity[2 * j + 1];
        w += float3(gx * ux + gy * uy, gx * ux - gy * uy, gy * ux + gx * uy);
    }
    w *= HALF_ROOT2;
    output[3 * a] = w.x;
    output[3 * a + 1] = w.y;
    output[3 * a + 2] = w.z;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void apply_bt(uint3 tid : SV_DispatchThreadID, uniform uint node_count,
              StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
              StructuredBuffer<float> gradients, StructuredBuffer<float> stress,
              RWStructuredBuffer<float> output)
{
    uint j = tid.x;
    if (j >= node_count) return;
    float2 f = float2(0.0, 0.0);
    for (uint i = 0; i < 9; i++)
    {
        uint a = transpose_rows[9 * j + i];
        if (a == INVALID) break;
        uint s = transpose_slots[9 * j + i];
        float gx = gradients[2 * (9 * a + s)];
        float gy = gradients[2 * (9 * a + s) + 1];
        float z0 = stress[3 * a];
        float z1 = stress[3 * a + 1];
        float z2 = stress[3 * a + 2];
        f += float2((z0 + z1) * gx + z2 * gy, (z0 - z1) * gy + z2 * gx);
    }
    f *= HALF_ROOT2;
    output[2 * j] = f.x;
    output[2 * j + 1] = f.y;
}
