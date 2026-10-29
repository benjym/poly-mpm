// Compact Q1 nodal-stress operator: one integrated gradient 2-vector per
// stress-node/velocity-node pair. Row-owned assembly and B.u, node-owned
// B.T z. No atomics in floating point; every write has a single owner.

static const uint INVALID = 0xffffffffu;
static const float HALF_ROOT2 = 0.70710678118654752;
static const float FLT_EPSILON = 1.1920928955078125e-7;

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
                        RWStructuredBuffer<float> gradients, RWStructuredBuffer<float> occupied,
                        RWStructuredBuffer<uint> status)
{
    uint a = tid.x;
    if (a >= row_count) return;
    uint node = row_nodes[a];
    int x = int(node % nx);
    int y = int(node / nx);
    int columns = int(nx) - 1;
    float2 sum[9];
    for (uint s = 0; s < 9; s++) sum[s] = float2(0.0, 0.0);
    float occupancy = 0.0;
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
                occupancy += scale;
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
    occupied[a] = occupancy;
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

// ---- Particle-grid transfers (granular.step ordering) ----------------------

static const uint FIXED_X = 1;
static const uint FIXED_Y = 2;
// fixed_components: bits 1/2 an x/y wall; one-sided edges leave it while the
// node moves in +x (4), -x (8), +y (16) or -y (32).
static const uint LEAVES_POSITIVE_X = 4;
static const uint LEAVES_NEGATIVE_X = 8;
static const uint LEAVES_POSITIVE_Y = 16;
static const uint LEAVES_NEGATIVE_Y = 32;

// grid.wall_contact: the walls (bit 1 x, bit 2 y) a node at velocity v presses.
uint wall_contact(uint walls, float2 v)
{
    uint contact = walls & (FIXED_X | FIXED_Y);
    if (((walls & LEAVES_POSITIVE_X) != 0 && v.x > 0.0) || ((walls & LEAVES_NEGATIVE_X) != 0 && v.x < 0.0))
        contact &= ~FIXED_X;
    if (((walls & LEAVES_POSITIVE_Y) != 0 && v.y > 0.0) || ((walls & LEAVES_NEGATIVE_Y) != 0 && v.y < 0.0))
        contact &= ~FIXED_Y;
    return contact;
}

// Fixed components of the walls in contact; a rough wall also fixes the tangent.
uint fixed_of(uint contact, uint rough)
{
    return rough != 0 && contact != 0 ? FIXED_X | FIXED_Y : contact;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void p2g_nodes(uint3 tid : SV_DispatchThreadID, uniform uint node_count,
               uniform uint nx, uniform uint ny, uniform uint periodic,
               uniform float dx, uniform float dy, uniform float dt, uniform float m_tol,
               uniform float body_x, uniform float body_y, uniform uint affine, uniform uint rough,
               StructuredBuffer<uint> fixed, StructuredBuffer<uint> cell_offsets,
               StructuredBuffer<uint> cell_counts, StructuredBuffer<uint> permutation,
               StructuredBuffer<float> local_x, StructuredBuffer<float> local_y,
               StructuredBuffer<float> particle_mass, StructuredBuffer<float> particle_velocity,
               StructuredBuffer<float> particle_gradient,
               RWStructuredBuffer<float> mass, RWStructuredBuffer<float> reference,
               RWStructuredBuffer<float> momentum, RWStructuredBuffer<float> trial,
               RWStructuredBuffer<float> inverse_mass, RWStructuredBuffer<uint> contact)
{
    // Node-owned gather of mass and momentum over the four adjacent cells.
    // An alias column gathers as its canonical node, so the periodic merge
    // and copy-back of grid.make_cyclic need no separate pass.
    uint j = tid.x;
    if (j >= node_count) return;
    uint column = j % nx, row = j / nx;
    uint node = canonical(column, row, nx, periodic);
    int x = int(node % nx), y = int(node / nx);
    int columns = int(nx) - 1;
    float m = 0.0;
    float2 q = float2(0.0, 0.0);
    int previous = -1;
    for (int ox = -1; ox <= 0; ox++)
    {
        int cx = x + ox;
        if (periodic != 0) cx = ((cx % columns) + columns) % columns;
        else if (cx < 0 || cx >= columns) continue;
        if (cx == previous) continue;
        previous = cx;
        for (int oy = -1; oy <= 0; oy++)
        {
            int cy = y + oy;
            if (cy < 0 || cy >= int(ny) - 1) continue;
            uint corners[4] = {
                canonical(uint(cx), uint(cy), nx, periodic),
                canonical(uint(cx + 1), uint(cy), nx, periodic),
                canonical(uint(cx + 1), uint(cy + 1), nx, periodic),
                canonical(uint(cx), uint(cy + 1), nx, periodic),
            };
            float2 corner_offset[4] = { float2(0, 0), float2(1, 0), float2(1, 1), float2(0, 1) };
            uint cell = uint(cx + columns * cy);
            uint begin = cell_offsets[cell];
            uint end = begin + cell_counts[cell];
            for (uint k = begin; k < end; k++)
            {
                uint p = permutation[k];
                float fx = local_x[p], fy = local_y[p];
                float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
                float2 v = float2(particle_velocity[2 * p], particle_velocity[2 * p + 1]);
                float2x2 L = float2x2(particle_gradient[4 * p], particle_gradient[4 * p + 1],
                                      particle_gradient[4 * p + 2], particle_gradient[4 * p + 3]);
                for (uint c = 0; c < 4; c++)
                {
                    if (corners[c] != node) continue;
                    float weight = w[c] * particle_mass[p];
                    float2 velocity = v;
                    // apic_gradient: v_p + L_p (x_i - x_p), with L retained between steps.
                    if (affine != 0)
                        velocity += mul(L, (corner_offset[c] - float2(fx, fy)) * float2(dx, dy));
                    m += weight;
                    q += weight * velocity;
                }
            }
        }
    }
    bool valid = m > m_tol;
    reference[2 * j] = q.x;
    reference[2 * j + 1] = q.y;
    // Walls in contact at the gathered velocity mask the momentum
    // (particles.get_nodal_mass_momentum); walls in contact after dt * m g mask
    // the force and fix the trial velocity (grid.update_momentum).
    uint walls = fixed[j];
    uint gathered = fixed_of(wall_contact(walls, m != 0.0 ? q / m : float2(0.0, 0.0)), rough);
    if ((gathered & FIXED_X) != 0) q.x = 0.0;
    if ((gathered & FIXED_Y) != 0) q.y = 0.0;
    float2 force = m * float2(body_x, body_y);
    uint touching = wall_contact(walls, m != 0.0 ? (q + dt * force) / m : float2(0.0, 0.0));
    contact[j] = touching;
    uint mask = fixed_of(touching, rough);
    if ((mask & FIXED_X) != 0) { q.x = 0.0; force.x = 0.0; }
    if ((mask & FIXED_Y) != 0) { q.y = 0.0; force.y = 0.0; }
    q += dt * force;
    mass[j] = m;
    momentum[2 * j] = q.x;
    momentum[2 * j + 1] = q.y;
    float2 velocity = valid ? q / m : float2(0.0, 0.0);
    trial[2 * j] = velocity.x;
    trial[2 * j + 1] = velocity.y;
    // Free components: resolved, unique (not a periodic alias) and not fixed.
    bool unique = node == j;
    inverse_mass[2 * j] = valid && unique && (mask & FIXED_X) == 0 ? 1.0 / m : 0.0;
    inverse_mass[2 * j + 1] = valid && unique && (mask & FIXED_Y) == 0 ? 1.0 / m : 0.0;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void wall_recheck(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform uint rough,
                  StructuredBuffer<uint> fixed, RWStructuredBuffer<float> velocity,
                  RWStructuredBuffer<uint> contact, RWStructuredBuffer<float> trial,
                  RWStructuredBuffer<float> inverse_mass, RWStructuredBuffer<uint> changed)
{
    // granular.step: a released wall node that the solve drives into its wall
    // is clamped again before the step is re-solved.
    uint j = tid.x;
    if (j >= node_count) return;
    uint before = contact[j];
    uint after = before | wall_contact(fixed[j], float2(velocity[2 * j], velocity[2 * j + 1]));
    if (after == before) return;
    contact[j] = after;
    uint mask = fixed_of(after, rough);
    if ((mask & FIXED_X) != 0) { trial[2 * j] = 0.0; inverse_mass[2 * j] = 0.0; }
    if ((mask & FIXED_Y) != 0) { trial[2 * j + 1] = 0.0; inverse_mass[2 * j + 1] = 0.0; }
    uint ignored;
    InterlockedOr(changed[0], 1, ignored);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void g2p_particles(uint3 tid : SV_DispatchThreadID, uniform uint count,
                   uniform uint nx, uniform uint periodic, uniform float dx, uniform float dy,
                   uniform float dt, uniform float m_tol, uniform uint transfer,
                   uniform float pic_fraction, StructuredBuffer<uint> keys,
                   StructuredBuffer<float> local_x, StructuredBuffer<float> local_y,
                   StructuredBuffer<float> mass, StructuredBuffer<float> reference,
                   StructuredBuffer<float> velocity,
                   RWStructuredBuffer<uint> cell_x, RWStructuredBuffer<uint> cell_y,
                   RWStructuredBuffer<float> fraction_x, RWStructuredBuffer<float> fraction_y,
                   RWStructuredBuffer<float> particle_velocity,
                   RWStructuredBuffer<float> particle_gradient, RWStructuredBuffer<float> shift)
{
    // transfer: 0 FLIP with PIC blend, 1 PIC, 2 apic_gradient. Uses the
    // start-of-step stencil, as particles.move_material_points does.
    uint p = tid.x;
    if (p >= count) return;
    uint columns = nx - 1;
    uint cx = keys[p] % columns, cy = keys[p] / columns;
    float fx = local_x[p], fy = local_y[p];
    uint corners[4] = {
        canonical(cx, cy, nx, periodic),
        canonical(cx + 1, cy, nx, periodic),
        canonical(cx + 1, cy + 1, nx, periodic),
        canonical(cx, cy + 1, nx, periodic),
    };
    float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
    float2 g[4] = {
        float2(-(1 - fy) / dx, -(1 - fx) / dy),
        float2((1 - fy) / dx, -fx / dy),
        float2(fy / dx, fx / dy),
        float2(-fy / dx, (1 - fx) / dy),
    };
    float2 pic = float2(0.0, 0.0), change = float2(0.0, 0.0);
    float2x2 L = float2x2(0, 0, 0, 0);
    for (uint c = 0; c < 4; c++)
    {
        uint j = corners[c];
        float m = mass[j];
        if (!(m > m_tol)) continue;
        float2 v = float2(velocity[2 * j], velocity[2 * j + 1]);
        float2 start = float2(reference[2 * j], reference[2 * j + 1]) / m;
        pic += w[c] * v;
        change += w[c] * (v - start);
        L += float2x2(v.x * g[c].x, v.x * g[c].y, v.y * g[c].x, v.y * g[c].y);
    }
    float2 old = float2(particle_velocity[2 * p], particle_velocity[2 * p + 1]);
    float2 next = transfer == 0 ? (1 - pic_fraction) * (old + change) + pic_fraction * pic : pic;
    particle_velocity[2 * p] = next.x;
    particle_velocity[2 * p + 1] = next.y;
    particle_gradient[4 * p] = L[0][0];
    particle_gradient[4 * p + 1] = L[0][1];
    particle_gradient[4 * p + 2] = L[1][0];
    particle_gradient[4 * p + 3] = L[1][1];
    // x += dt * pic in split coordinates; carry whole cells out of the fraction.
    float2 local = float2(fx, fy) + dt * pic / float2(dx, dy);
    float2 carry = floor(local);
    local -= carry;
    // Rounding can leave exactly one after subtraction of the carry.
    if (local.x >= 1.0) { local.x = 0.0; carry.x += 1.0; }
    if (local.y >= 1.0) { local.y = 0.0; carry.y += 1.0; }
    int unwrapped = int(cx) + int(carry.x);
    int nextx = periodic != 0 ? ((unwrapped % int(columns)) + int(columns)) % int(columns) : unwrapped;
    shift[p] = float(unwrapped - nextx);
    cell_x[p] = asuint(nextx);
    cell_y[p] = asuint(int(cy) + int(carry.y));
    fraction_x[p] = local.x;
    fraction_y[p] = local.y;
}

// exp(A) for a real 2x2 A: exp(s) [ch(d) I + sh(d)/d (A - s I)], s = tr/2,
// d^2 = -det(A - s I), with series for small |d| and either real or
// imaginary d. Applied after scaling by 2^-k; then squared k times.
float2x2 exponential(float2x2 A)
{
    float size = max(abs(A[0][0]) + abs(A[1][0]), abs(A[0][1]) + abs(A[1][1]));
    int k = size > 0.5 ? int(ceil(log2(size / 0.5))) : 0;
    k = min(k, 60);
    A *= exp2(-float(k));
    float s = 0.5 * (A[0][0] + A[1][1]);
    float2x2 B = A - float2x2(s, 0, 0, s);
    float delta = -(B[0][0] * B[1][1] - B[0][1] * B[1][0]);
    float even, odd;
    if (abs(delta) < 1e-3)
    {
        even = 1.0 + delta * (0.5 + delta / 24.0);
        odd = 1.0 + delta * (1.0 / 6.0 + delta / 120.0);
    }
    else if (delta > 0.0)
    {
        float d = sqrt(delta);
        even = cosh(d);
        odd = sinh(d) / d;
    }
    else
    {
        float d = sqrt(-delta);
        even = cos(d);
        odd = sin(d) / d;
    }
    float2x2 E = exp(s) * (float2x2(even, 0, 0, even) + odd * B);
    for (int i = 0; i < k; i++) E = mul(E, E);
    return E;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void deform_particles(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform float dt,
                      StructuredBuffer<float> particle_gradient, StructuredBuffer<float> reference_volume,
                      StructuredBuffer<float> particle_mass, RWStructuredBuffer<float> deformation,
                      RWStructuredBuffer<float> jacobian, RWStructuredBuffer<float> volume,
                      RWStructuredBuffer<float> density, RWStructuredBuffer<uint> status)
{
    // F <- blockdiag(exp(L dt), 1) F on the full 3x3 F (kinematics.advance).
    // J <- J exp(tr(L dt)) equals det(F) exactly; an FP32 3x3 determinant of a
    // strongly stretched F cancels catastrophically.
    uint p = tid.x;
    if (p >= count) return;
    float2x2 L = float2x2(particle_gradient[4 * p], particle_gradient[4 * p + 1],
                          particle_gradient[4 * p + 2], particle_gradient[4 * p + 3]);
    float2x2 E = exponential(dt * L);
    float3x3 F;
    for (uint i = 0; i < 9; i++) F[i / 3][i % 3] = deformation[9 * p + i];
    float3x3 increment = float3x3(E[0][0], E[0][1], 0, E[1][0], E[1][1], 0, 0, 0, 1);
    F = mul(increment, F);
    float J = jacobian[p] * exp(dt * (L[0][0] + L[1][1]));
    bool finite = true;
    for (uint i = 0; i < 9; i++) finite = finite && isfinite(F[i / 3][i % 3]);
    if (!(J > 0.0) || !isfinite(J) || !finite)
    {
        uint ignored;
        InterlockedOr(status[0], 64, ignored);
        return;
    }
    for (uint i = 0; i < 9; i++) deformation[9 * p + i] = F[i / 3][i % 3];
    jacobian[p] = J;
    volume[p] = reference_volume[p] * J;
    density[p] = particle_mass[p] / volume[p];
}
