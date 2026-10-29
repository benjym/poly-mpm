// End-of-step kernels for the resident granular timestep (granular.step).
// Column outputs are column-major [column][element] for ColumnReduction.

import mechanics;

[shader("compute")]
[numthreads(256, 1, 1)]
void packing_gap(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform uint stress_rows,
                 uniform float packing_max, uniform float dt, StructuredBuffer<float> occupied, StructuredBuffer<float> volumes,
                 StructuredBuffer<float> previous, RWStructuredBuffer<float> gap,
                 RWStructuredBuffer<float> initial)
{
    // granular.assemble: gap_n = max(clearance, 0) / (sqrt(2) dt), with
    // clearances at summation roundoff set to zero (here FP32 roundoff).
    // Rows without grain volume carry no constraint and start from zero; wall
    // rows have no gap and start from the previous traction.
    uint a = tid.x;
    if (a >= row_count) return;
    float occupancy = occupied[a];
    float capacity = packing_max * volumes[a];
    float clearance = capacity - occupancy;
    if (abs(clearance) <= 64.0 * FLT_EPSILON * max(capacity, occupancy)) clearance = 0.0;
    bool active = occupancy > 0.0;
    gap[3 * a] = active ? max(clearance, 0.0) / (sqrt(2.0) * dt) : 0.0;
    gap[3 * a + 1] = 0.0;
    gap[3 * a + 2] = 0.0;
    bool warm = active || a >= stress_rows;
    for (uint k = 0; k < 3; k++) initial[3 * a + k] = warm ? previous[3 * a + k] : 0.0;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void grid_finish(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform uint nx,
                 uniform uint periodic, uniform float dx, uniform float dy, uniform float x_m,
                 uniform float y_m, uniform float dt, uniform float m_tol, uniform float body_x,
                 uniform float body_y, StructuredBuffer<uint> transpose_rows,
                 StructuredBuffer<uint> transpose_slots, StructuredBuffer<float> gradients,
                 StructuredBuffer<float> contact, StructuredBuffer<float> mass,
                 StructuredBuffer<float> reference, RWStructuredBuffer<float> velocity,
                 RWStructuredBuffer<float> internal, RWStructuredBuffer<float> momentum,
                 RWStructuredBuffer<float> rate, RWStructuredBuffer<float> columns)
{
    // Final nodal velocity on periodic aliases, internal force fi = -B.T s
    // over stress rows (wall traction is a constraint), q = m u, q_dot and the
    // per-node constraint/external budget terms.
    uint j = tid.x;
    if (j >= node_count) return;
    uint column = j % nx, row = j / nx;
    uint c = canonical(column, row, nx, periodic);
    float2 u = float2(velocity[2 * c], velocity[2 * c + 1]);
    if (c != j)
    {
        velocity[2 * j] = u.x;
        velocity[2 * j + 1] = u.y;
    }
    float2 force = float2(0.0, 0.0);
    for (uint i = 0; i < 9; i++)
    {
        uint a = transpose_rows[9 * c + i];
        if (a == INVALID) break;
        uint s = transpose_slots[9 * c + i];
        if (s == WALL_SLOT) continue;
        force += mul(float3(contact[3 * a], contact[3 * a + 1], contact[3 * a + 2]), row_block(a, s, gradients));
    }
    float m = mass[j];
    float2 start = float2(reference[2 * j], reference[2 * j + 1]);
    float2 q = m * u;
    float2 fe = m * float2(body_x, body_y);
    internal[2 * j] = -force.x;
    internal[2 * j + 1] = -force.y;
    momentum[2 * j] = q.x;
    momentum[2 * j + 1] = q.y;
    rate[2 * j] = (q.x - start.x) / dt;
    rate[2 * j + 1] = (q.y - start.y) / dt;
    // Sums over unique physical nodes only.
    float weight = c == j ? 1.0 : 0.0;
    float2 reaction = q - start - dt * (fe + force);
    float2 midpoint = 0.5 * ((m > m_tol ? start / m : float2(0.0, 0.0)) + u);
    float2 x = float2(x_m + float(column) * dx, y_m + float(row) * dy);
    float values[8] = {
        reaction.x, reaction.y, x.x * reaction.y - x.y * reaction.x, dot(reaction, midpoint),
        dt * dot(fe, midpoint), dt * fe.x, dt * fe.y, dt * (x.x * fe.y - x.y * fe.x),
    };
    for (uint k = 0; k < 8; k++) columns[k * node_count + j] = weight * values[k];
    // Largest nodal velocity component, for initialise.update_timestep.
    columns[8 * node_count + j] = m > m_tol ? max(abs(u.x), abs(u.y)) : 0.0;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void contact_rows(uint3 tid : SV_DispatchThreadID, uniform uint row_count,
                  StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                  StructuredBuffer<float> contact, StructuredBuffer<float> velocity,
                  StructuredBuffer<float> occupied, StructuredBuffer<float> volumes,
                  RWStructuredBuffer<float> previous, RWStructuredBuffer<float> columns)
{
    // Contact power s . B u, nodal packing, active-row count; keep the
    // contact for the next warm start.
    uint a = tid.x;
    if (a >= row_count) return;
    float3 w = float3(0.0, 0.0, 0.0);
    for (uint s = 0; s < 9; s++)
    {
        uint j = neighbors[9 * a + s];
        if (j == INVALID) continue;
        w += mul(row_block(a, row_slot(a, s, neighbors), gradients), float2(velocity[2 * j], velocity[2 * j + 1]));
    }
    float3 z = float3(contact[3 * a], contact[3 * a + 1], contact[3 * a + 2]);
    for (uint k = 0; k < 3; k++) previous[3 * a + k] = z[k];
    float volume = volumes[a];
    columns[a] = dot(z, w);
    columns[row_count + a] = volume > 0.0 ? occupied[a] / volume : 0.0;
    columns[2 * row_count + a] = occupied[a] > 0.0 ? 1.0 : 0.0;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void packing_rows(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform float packing_max,
                  StructuredBuffer<float> occupied, StructuredBuffer<float> volumes,
                  RWStructuredBuffer<float> columns)
{
    // End-of-step nodal packing from current positions (granular.packing_fraction).
    uint a = tid.x;
    if (a >= row_count) return;
    float volume = volumes[a];
    float packing = volume > 0.0 ? occupied[a] / volume : 0.0;
    columns[a] = packing;
    columns[row_count + a] = max(packing - packing_max, 0.0) * volume;
}

void kahan(RWStructuredBuffer<float> sum, RWStructuredBuffer<float> compensation, uint p, float value)
{
    float y = value - compensation[p];
    float t = sum[p] + y;
    compensation[p] = (t - sum[p]) - y;
    sum[p] = t;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void particle_update(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform uint nx,
                     uniform uint periodic, uniform float dt, uniform float m_tol, uniform float rho_s,
                     StructuredBuffer<uint> keys, StructuredBuffer<float> local_x,
                     StructuredBuffer<float> local_y, StructuredBuffer<float> mass,
                     StructuredBuffer<float> contact, StructuredBuffer<uint> node_row,
                     StructuredBuffer<float> gradient, StructuredBuffer<float> density,
                     StructuredBuffer<float> volume, RWStructuredBuffer<float> stress,
                     RWStructuredBuffer<float> stress_change, RWStructuredBuffer<float> strain,
                     RWStructuredBuffer<float> work, RWStructuredBuffer<float> work_compensation)
{
    // Stress from the gathered nodal contact (start-of-step stencil, valid
    // nodes only), scaled by rho / rho_s; strain and work increments.
    // Planar symmetric tensors are stored as (xx, yy, zz, xy).
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
    float3 coefficient = float3(0.0, 0.0, 0.0);
    for (uint c = 0; c < 4; c++)
    {
        if (!(mass[corners[c]] > m_tol)) continue;
        uint a = node_row[corners[c]];
        coefficient += w[c] * float3(contact[3 * a], contact[3 * a + 1], contact[3 * a + 2]);
    }
    float scale = -HALF_ROOT2 * density[p] / rho_s;
    float xx = scale * (coefficient.x + coefficient.y);
    float yy = scale * (coefficient.x - coefficient.y);
    float xy = scale * coefficient.z;
    float4 next = float4(xx, yy, 0.5 * (xx + yy), xy);
    float4 old = float4(stress[4 * p], stress[4 * p + 1], stress[4 * p + 2], stress[4 * p + 3]);
    float L00 = gradient[4 * p], L01 = gradient[4 * p + 1];
    float L10 = gradient[4 * p + 2], L11 = gradient[4 * p + 3];
    float4 increment = dt * float4(L00, L11, 0.0, 0.5 * (L01 + L10));
    for (uint k = 0; k < 4; k++)
    {
        stress[4 * p + k] = next[k];
        stress_change[4 * p + k] = next[k] - old[k];
        strain[4 * p + k] += increment[k];
    }
    float power = next.x * increment.x + next.y * increment.y + 2.0 * next.w * increment.w;
    kahan(work, work_compensation, p, volume[p] * power);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void particle_inventory(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform uint nx,
                        uniform uint periodic, uniform float dx, uniform float dy, uniform float x_m,
                        uniform float y_m, uniform float m_tol, uniform float rho_s,
                        uniform uint affine, uniform uint outlets, StructuredBuffer<uint> keys,
                        StructuredBuffer<float> local_x, StructuredBuffer<float> local_y,
                        StructuredBuffer<float> nodal_mass, StructuredBuffer<float> mass,
                        StructuredBuffer<float> velocity, StructuredBuffer<float> gradient,
                        StructuredBuffer<uint> cell_x, StructuredBuffer<uint> cell_y,
                        StructuredBuffer<float> fraction_x, StructuredBuffer<float> fraction_y,
                        StructuredBuffer<float> shift, StructuredBuffer<float> jacobian,
                        StructuredBuffer<float> density, StructuredBuffer<float> work,
                        StructuredBuffer<float> work_compensation, RWStructuredBuffer<float> columns)
{
    // diagnostics.inventory after advection, plus the angular momentum before
    // the periodic wrap (process_boundaries), the CFL speed, support counts and
    // particles that process_boundaries removes (outlets bit 1 left, bit 2 bottom).
    uint p = tid.x;
    if (p >= count) return;
    float m = mass[p];
    float2 v = float2(velocity[2 * p], velocity[2 * p + 1]);
    float2 f = float2(fraction_x[p], fraction_y[p]);
    float2 x = float2(x_m + (float(asint(cell_x[p])) + f.x) * dx, y_m + (float(asint(cell_y[p])) + f.y) * dy);
    float before = x.x + shift[p] * dx;
    float L00 = gradient[4 * p], L01 = gradient[4 * p + 1];
    float L10 = gradient[4 * p + 2], L11 = gradient[4 * p + 3];
    float spin = 0.0, affine_kinetic = 0.0;
    if (affine != 0)
    {
        float2 d = float2(dx * dx, dy * dy) * f * (1.0 - f);
        spin = m * (L10 * d.x - L01 * d.y);
        affine_kinetic = 0.5 * m * (L00 * L00 * d.x + L01 * L01 * d.y + L10 * L10 * d.x + L11 * L11 * d.y);
    }
    // Start-of-step stencil coverage by valid nodes (Diagnostics.finish).
    uint columns_x = nx - 1;
    uint cx = keys[p] % columns_x, cy = keys[p] / columns_x;
    float gx = local_x[p], gy = local_y[p];
    uint corners[4] = {
        canonical(cx, cy, nx, periodic),
        canonical(cx + 1, cy, nx, periodic),
        canonical(cx + 1, cy + 1, nx, periodic),
        canonical(cx, cy + 1, nx, periodic),
    };
    float w[4] = { (1 - gx) * (1 - gy), gx * (1 - gy), gx * gy, (1 - gx) * gy };
    float coverage = 0.0;
    for (uint c = 0; c < 4; c++) coverage += nodal_mass[corners[c]] > m_tol ? w[c] : 0.0;
    float total = work[p] - work_compensation[p];
    bool leaving = ((outlets & 1) != 0 && x.x < x_m + dx) || ((outlets & 2) != 0 && x.y < y_m);
    float values[16] = {
        m * v.x, m * v.y, m * (x.x * v.y - x.y * v.x) + spin, m * (before * v.y - x.y * v.x) + spin,
        0.5 * m * dot(v, v), affine_kinetic, total, jacobian[p], jacobian[p], density[p], density[p],
        length(v), density[p] / rho_s, coverage < 1.0 - 1e-6 ? 1.0 : 0.0, m, leaving ? 1.0 : 0.0,
    };
    for (uint k = 0; k < 16; k++) columns[k * count + p] = values[k];
}

[shader("compute")]
[numthreads(64, 1, 1)]
void copy_record(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform uint offset,
                 StructuredBuffer<float> source, RWStructuredBuffer<float> record)
{
    if (tid.x < count) record[offset + tid.x] = source[tid.x];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void occupancy_rows(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform uint nx,
                    uniform uint ny, uniform uint periodic, StructuredBuffer<uint> row_nodes,
                    StructuredBuffer<uint> cell_offsets, StructuredBuffer<uint> cell_counts,
                    StructuredBuffer<uint> permutation, StructuredBuffer<float> local_x,
                    StructuredBuffer<float> local_y, StructuredBuffer<float> volume,
                    RWStructuredBuffer<float> occupied)
{
    // Occupied grain volume per stress row, without the gradient assembly.
    uint a = tid.x;
    if (a >= row_count) return;
    uint node = row_nodes[a];
    int x = int(node % nx), y = int(node / nx);
    int columns = int(nx) - 1;
    float occupancy = 0.0;
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
            uint cell = uint(cx + columns * cy);
            uint begin = cell_offsets[cell];
            uint end = begin + cell_counts[cell];
            for (uint k = begin; k < end; k++)
            {
                uint p = permutation[k];
                float fx = local_x[p], fy = local_y[p];
                float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
                for (uint c = 0; c < 4; c++)
                    if (corners[c] == node) occupancy += volume[p] * w[c];
            }
        }
    }
    occupied[a] = occupancy;
}
