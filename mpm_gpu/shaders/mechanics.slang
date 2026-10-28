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

// Per-slot 3x2 block of the intrinsic basis applied to a gradient: rows are
// stress components k, columns velocity components c, B[k][c] = basis_kcd g_d.
float3x2 block(float gx, float gy)
{
    return HALF_ROOT2 * float3x2(gx, gy, gx, -gy, gy, gx);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void absolute_force(uint3 tid : SV_DispatchThreadID, uniform uint node_count,
                    StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
                    StructuredBuffer<float> gradients, RWStructuredBuffer<float> output)
{
    // |B.T| 1 per velocity component, for the conservative norm bound.
    uint j = tid.x;
    if (j >= node_count) return;
    float2 f = float2(0.0, 0.0);
    for (uint i = 0; i < 9; i++)
    {
        uint a = transpose_rows[9 * j + i];
        if (a == INVALID) break;
        uint s = transpose_slots[9 * j + i];
        float3x2 B = block(gradients[2 * (9 * a + s)], gradients[2 * (9 * a + s) + 1]);
        for (uint k = 0; k < 3; k++) f += abs(B[k]);
    }
    output[2 * j] = f.x;
    output[2 * j + 1] = f.y;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void local_metric(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform float dt,
                  StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                  StructuredBuffer<float> inverse_mass, StructuredBuffer<float> absolute,
                  RWStructuredBuffer<float> metric, RWStructuredBuffer<float> bound)
{
    // metric = dt sum_s B_s M^-1 B_s^T (upper triangle 00,01,02,11,12,22);
    // bound = max_k dt (|B| M^-1 |B.T| 1)_k >= ||B M^-1 B.T||_inf.
    uint a = tid.x;
    if (a >= row_count) return;
    float3x3 W = float3x3(0, 0, 0, 0, 0, 0, 0, 0, 0);
    float3 b = float3(0.0, 0.0, 0.0);
    for (uint s = 0; s < 9; s++)
    {
        uint j = neighbors[9 * a + s];
        if (j == INVALID) continue;
        float3x2 B = block(gradients[2 * (9 * a + s)], gradients[2 * (9 * a + s) + 1]);
        float2 m = float2(inverse_mass[2 * j], inverse_mass[2 * j + 1]);
        float2 weighted = m * float2(absolute[2 * j], absolute[2 * j + 1]);
        for (uint k = 0; k < 3; k++)
        {
            b[k] += dot(abs(B[k]), weighted);
            for (uint l = 0; l < 3; l++) W[k][l] += dot(B[k] * m, B[l]);
        }
    }
    W *= dt;
    metric[6 * a] = W[0][0];
    metric[6 * a + 1] = W[0][1];
    metric[6 * a + 2] = W[0][2];
    metric[6 * a + 3] = W[1][1];
    metric[6 * a + 4] = W[1][2];
    metric[6 * a + 5] = W[2][2];
    bound[a] = dt * max(b.x, max(b.y, b.z));
}

[shader("compute")]
[numthreads(256, 1, 1)]
void metric_finish(uint3 tid : SV_DispatchThreadID, uniform uint row_count,
                   StructuredBuffer<float> norm, RWStructuredBuffer<float> metric,
                   RWStructuredBuffer<float> diagonal, RWStructuredBuffer<float> rotation,
                   RWStructuredBuffer<uint> status)
{
    // Normalise by the global bound, add the proximal regulariser and cache
    // the analytic tangential eigensystem. A preconditioner, not compliance.
    uint a = tid.x;
    if (a >= row_count) return;
    float scale = norm[0];
    if (!(scale > 0.0) || !isfinite(scale))
    {
        if (a == 0)
        {
            uint ignored;
            InterlockedOr(status[0], 8, ignored);
        }
        scale = 1.0;
    }
    float W[6];
    for (uint i = 0; i < 6; i++) W[i] = metric[6 * a + i] / scale;
    for (uint i = 0; i < 6; i++) metric[6 * a + i] = W[i];
    float regularizer = 1e-6 * max(W[0] + W[3] + W[5], 1e-12);
    // Symmetric Jacobi rotation of [[t11, t12], [t12, t22]] (Golub and Van
    // Loan, algorithm 8.5.1), stored as columns (c, s), (-s, c). Metal's
    // atan2(0, 0) is not usable for isotropic blocks.
    float t11 = W[3], t12 = W[4], t22 = W[5];
    float c = 1.0, s = 0.0, first = t11, second = t22;
    if (t12 != 0.0)
    {
        float tau = (t22 - t11) / (2.0 * t12);
        float t = (tau >= 0.0 ? 1.0 : -1.0) / (abs(tau) + sqrt(1.0 + tau * tau));
        c = rsqrt(1.0 + t * t);
        s = -t * c;
        first = t11 - t * t12;
        second = t22 + t * t12;
    }
    diagonal[3 * a] = W[0] + regularizer;
    diagonal[3 * a + 1] = first + regularizer;
    diagonal[3 * a + 2] = second + regularizer;
    rotation[2 * a] = c;
    rotation[2 * a + 1] = s;
}

// Exact local solve for a positive diagonal metric (granular_mixed.diagonal_coulomb).
// Status bit 16: nonpositive/nonfinite metric. Bit 32: sliding solve unconverged.
float3 coulomb_local(float3 d, float3 q, float mu, inout uint flags)
{
    if (!(d.x > 0.0 && d.y > 0.0 && d.z > 0.0) || !all(isfinite(d)) || !all(isfinite(q)))
    {
        flags |= 16;
        return float3(0.0, 0.0, 0.0);
    }
    float normal = max(-q.x / d.x, 0.0);
    float radius = mu * normal;
    float2 tangent = -q.yz / d.yz;
    // Opening (normal = 0) and zero friction leave no tangential reaction.
    if (!(radius > 0.0)) return float3(normal, 0.0, 0.0);
    if (length(tangent) <= radius) return float3(normal, tangent);
    // Sliding: ||b / (d_t + alpha)|| = 1, with b scaled by the yield radius.
    float2 dt = d.yz;
    float2 b = q.yz / radius;
    float size = length(b);
    float lower = max(0.0, size - max(dt.x, dt.y));
    float upper = max(lower, size - min(dt.x, dt.y));
    float alpha = lower;
    float2 value = b / (dt + alpha);
    bool converged = false;
    for (uint i = 0; i < 64; i++)
    {
        float2 denominator = dt + alpha;
        value = b / denominator;
        float magnitude = length(value);
        float error = magnitude - 1.0;
        if (abs(error) <= 4e-7) { converged = true; break; }
        if (error > 0.0) lower = alpha;
        else upper = alpha;
        // FP32 stagnation: the bracket has collapsed to adjacent values.
        if (upper - lower <= 2.4e-7 * max(upper, 1e-30)) { converged = true; break; }
        float derivative = -dot(value, value / denominator) / magnitude;
        float newton = alpha - error / derivative;
        float next = (newton >= lower && newton <= upper) ? newton : 0.5 * (lower + upper);
        if (next == alpha) { converged = true; break; }
        alpha = next;
    }
    if (!converged) flags |= 32;
    return float3(normal, -radius * value);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void local_coulomb_test(uint3 tid : SV_DispatchThreadID, uniform uint count,
                        StructuredBuffer<float> diagonal, StructuredBuffer<float> q,
                        StructuredBuffer<float> mu, RWStructuredBuffer<float> output,
                        RWStructuredBuffer<uint> status)
{
    uint i = tid.x;
    if (i >= count) return;
    uint flags = 0;
    float3 r = coulomb_local(float3(diagonal[3 * i], diagonal[3 * i + 1], diagonal[3 * i + 2]),
                             float3(q[3 * i], q[3 * i + 1], q[3 * i + 2]), mu[i], flags);
    output[3 * i] = r.x;
    output[3 * i + 1] = r.y;
    output[3 * i + 2] = r.z;
    status[i] = flags;
}
