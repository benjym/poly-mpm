// Grid segregation for the resident granular step (integrators.segregation_update).
// Node- and face-owned gathers only; species are the innermost index ([item][k]).
// Static face lists, padding and wall closures come from the host, built with
// the same helpers as the CPU scheme.

import mechanics;

static const uint NS_MAX = 8;
static const float FLT_EPS = 1.1920928955078125e-7;

float safe_divide(float a, float b)
{
    return b != 0.0 ? a / b : 0.0;
}

// ---- Nodal fractions and diagnostic deposits ---------------------------------

[shader("compute")]
[numthreads(64, 1, 1)]
void seg_nodes(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform uint nx, uniform uint ny,
               uniform uint periodic, uniform uint ns, uniform float rho_s, uniform float root_rho_s,
               StructuredBuffer<uint> cell_offsets, StructuredBuffer<uint> cell_counts,
               StructuredBuffer<uint> permutation, StructuredBuffer<float> local_x,
               StructuredBuffer<float> local_y, StructuredBuffer<float> particle_mass,
               StructuredBuffer<float> phi, StructuredBuffer<float> sizes,
               StructuredBuffer<float> stress, StructuredBuffer<float> gradient,
               StructuredBuffer<float> density, RWStructuredBuffer<float> phim,
               RWStructuredBuffer<float> pressure, RWStructuredBuffer<float> inertial)
{
    // sum N m phi_k, sum N m p and sum N m I over the start-of-step stencil, with
    // the particle pressure and inertial number of this step (deposit_diagnostics).
    uint j = tid.x;
    if (j >= node_count) return;
    uint node = canonical(j % nx, j / nx, nx, periodic);
    int x = int(node % nx), y = int(node / nx);
    int columns = int(nx) - 1;
    float sums[NS_MAX];
    for (uint k = 0; k < NS_MAX; k++) sums[k] = 0.0;
    float p_sum = 0.0, i_sum = 0.0;
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
            uint begin = cell_offsets[cell], end = begin + cell_counts[cell];
            for (uint e = begin; e < end; e++)
            {
                uint p = permutation[e];
                float fx = local_x[p], fy = local_y[p];
                float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
                float weight = 0.0;
                for (uint c = 0; c < 4; c++) weight += corners[c] == node ? w[c] : 0.0;
                if (weight == 0.0) continue;
                float m = weight * particle_mass[p];
                float size = 0.0;
                for (uint k = 0; k < ns; k++)
                {
                    float value = phi[ns * p + k];
                    sums[k] += m * value;
                    size += value * sizes[k];
                }
                float xx = stress[4 * p], yy = stress[4 * p + 1], zz = stress[4 * p + 2];
                float trace = xx + yy + zz;
                p_sum += m * trace / 3.0;
                // gammadot = sqrt(2) |dev(sym L)| in 2D; phase pressure from rho_s / rho.
                float L00 = gradient[4 * p], L01 = gradient[4 * p + 1];
                float L10 = gradient[4 * p + 2], L11 = gradient[4 * p + 3];
                float half = 0.5 * (L00 - L11), shear = 0.5 * (L01 + L10);
                float gammadot = 1.4142135623730951 * sqrt(2.0 * (half * half + shear * shear));
                float phase = -trace / 3.0 * rho_s / density[p];
                i_sum += phase > 0.0 ? m * size * gammadot * root_rho_s / sqrt(phase) : 0.0;
            }
        }
    }
    for (uint k = 0; k < ns; k++) phim[ns * j + k] = sums[k];
    pressure[j] = p_sum;
    inertial[j] = i_sum;
}

// Masked central difference of grid.calculate_gradient (non-smooth branch) on
// the full node array, along x (axis 0) or y (axis 1).
float masked_derivative(uint j, uint axis, uint nx, uint ny, float spacing, StructuredBuffer<float> mass,
                        float m_tol, float values[3], bool valid[3])
{
    // values/valid: [previous, centre, next] along the axis (as available).
    uint position = axis == 0 ? j % nx : j / nx;
    uint count = axis == 0 ? nx : ny;
    float result;
    if (position == 0) result = valid[1] && valid[2] ? (values[2] - values[1]) / spacing : 0.0;
    else if (position == count - 1) result = valid[0] && valid[1] ? (values[1] - values[0]) / spacing : 0.0;
    else result = valid[0] && valid[2] ? (values[2] - values[0]) / (2.0 * spacing) : 0.0;
    return isfinite(result) ? result : 0.0;
}

float nodal_value(uint field, uint j, StructuredBuffer<float> mass, float m_tol, StructuredBuffer<float> velocity,
                  StructuredBuffer<float> pressure, StructuredBuffer<float> inertial)
{
    // field 0/1: nodal velocity components (zero below m_tol); 2: pk = -p_phys I_phys.
    float m = mass[j];
    if (field < 2) return m > m_tol ? velocity[2 * j + field] : 0.0;
    return -safe_divide(pressure[j], m) * safe_divide(inertial[j], m);
}

float axis_derivative(uint field, uint j, uint axis, uint nx, uint ny, float spacing, float m_tol,
                      StructuredBuffer<float> mass, StructuredBuffer<float> velocity,
                      StructuredBuffer<float> pressure, StructuredBuffer<float> inertial)
{
    uint position = axis == 0 ? j % nx : j / nx;
    uint count = axis == 0 ? nx : ny;
    uint stride = axis == 0 ? 1 : nx;
    float values[3] = { 0.0, 0.0, 0.0 };
    bool valid[3] = { false, false, false };
    for (int o = -1; o <= 1; o++)
    {
        int q = int(position) + o;
        if (q < 0 || q >= int(count)) continue;
        uint n = uint(int(j) + o * int(stride));
        values[o + 1] = nodal_value(field, n, mass, m_tol, velocity, pressure, inertial);
        valid[o + 1] = !(mass[n] <= m_tol);
    }
    return masked_derivative(j, axis, nx, ny, spacing, mass, m_tol, values, valid);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void seg_fields(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform uint nx, uniform uint ny,
                uniform uint ns, uniform float dx, uniform float dy, uniform float m_tol,
                StructuredBuffer<float> mass, StructuredBuffer<float> velocity,
                StructuredBuffer<float> pressure, StructuredBuffer<float> inertial,
                StructuredBuffer<float> phim, StructuredBuffer<float> sizes,
                RWStructuredBuffer<float> phi_nodes, RWStructuredBuffer<float> grad_pk,
                RWStructuredBuffer<float> diffusivity, RWStructuredBuffer<float> columns)
{
    // phi = phim / m; shear rate du/dy + dv/dx; grad pk; diffusivity
    // 10 s^3 |shear| / sqrt(max(I_phys, 0)); and the stability maxima.
    uint j = tid.x;
    if (j >= node_count) return;
    float m = mass[j];
    float size = 0.0;
    for (uint k = 0; k < ns; k++)
    {
        float value = safe_divide(phim[ns * j + k], m);
        phi_nodes[ns * j + k] = value;
        size += value * sizes[k];
    }
    float shear = axis_derivative(0, j, 1, nx, ny, dy, m_tol, mass, velocity, pressure, inertial)
                  + axis_derivative(1, j, 0, nx, ny, dx, m_tol, mass, velocity, pressure, inertial);
    float gx = axis_derivative(2, j, 0, nx, ny, dx, m_tol, mass, velocity, pressure, inertial);
    float gy = axis_derivative(2, j, 1, nx, ny, dy, m_tol, mass, velocity, pressure, inertial);
    grad_pk[2 * j] = gx;
    grad_pk[2 * j + 1] = gy;
    float inertial_phys = safe_divide(inertial[j], m);
    float d = safe_divide(10.0 * size * size * size * abs(shear), sqrt(max(inertial_phys, 0.0)));
    diffusivity[j] = d;
    columns[j] = length(float2(gx, gy));
    columns[node_count + j] = max(abs(gx), abs(gy));
    columns[2 * node_count + j] = d;
}

// ---- Central-upwind face transfers --------------------------------------------

void advection(float C[NS_MAX], float forcing, uint ns, StructuredBuffer<float> sizes, out float flux[NS_MAX],
               out float speed)
{
    // integrators._advection: relative flux and a bound on its Jacobian's row norm.
    float harmonic = 0.0, inverse_sum = 0.0;
    for (uint k = 0; k < ns; k++)
    {
        harmonic += C[k] / sizes[k];
        inverse_sum += 1.0 / sizes[k];
    }
    speed = 0.0;
    for (uint k = 0; k < NS_MAX; k++) flux[k] = 0.0;
    for (uint k = 0; k < ns; k++)
    {
        float velocity = forcing * (safe_divide(1.0, harmonic * sizes[k]) - 1.0);
        flux[k] = C[k] * velocity;
        float rank = safe_divide(forcing * C[k], sizes[k] * harmonic * harmonic);
        float inverse = 1.0 / sizes[k];
        float row = abs(velocity - rank * inverse) + abs(rank) * (inverse_sum - inverse);
        speed = max(speed, row);
    }
}

void limited_slope(uint before, uint centre, uint after, uint ns, float spacing, float m_tol,
                   StructuredBuffer<float> mass, RWStructuredBuffer<float> phi_nodes, out float slope[NS_MAX])
{
    // integrators._slopes: one limiter for the whole mixture.
    float theta = 1.0;
    float centred[NS_MAX];
    for (uint k = 0; k < NS_MAX; k++) { slope[k] = 0.0; centred[k] = 0.0; }
    for (uint k = 0; k < ns; k++)
    {
        float c = phi_nodes[ns * centre + k];
        float left = (c - phi_nodes[ns * before + k]) / spacing;
        float right = (phi_nodes[ns * after + k] - c) / spacing;
        float middle = 0.5 * (left + right);
        float bound = left * right > 0.0 ? min(abs(left), abs(right)) : 0.0;
        float ratio = middle != 0.0 ? bound / abs(middle) : 1.0;
        theta = min(theta, ratio);
        centred[k] = middle;
    }
    bool stencil = mass[before] > m_tol && mass[centre] > m_tol && mass[after] > m_tol;
    for (uint k = 0; k < ns; k++) slope[k] = stencil ? centred[k] * theta : 0.0;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void seg_faces(uint3 tid : SV_DispatchThreadID, uniform uint face_count, uniform uint ns, uniform float dt,
               uniform float c_seg, uniform float m_tol, StructuredBuffer<uint> faces,
               StructuredBuffer<float> face_geometry, StructuredBuffer<float> mass,
               StructuredBuffer<float> nodal_volume, StructuredBuffer<float> grad_pk,
               StructuredBuffer<float> diffusivity, StructuredBuffer<float> sizes,
               RWStructuredBuffer<float> phi_nodes, RWStructuredBuffer<float> transfer)
{
    // faces: [a, b, before_a, after_b, axis, wall] per face; geometry: [spacing, area].
    uint f = tid.x;
    if (f >= face_count) return;
    uint a = faces[6 * f], b = faces[6 * f + 1], ap = faces[6 * f + 2], bn = faces[6 * f + 3];
    uint axis = faces[6 * f + 4];
    bool wall = faces[6 * f + 5] != 0;
    float spacing = face_geometry[2 * f], area = face_geometry[2 * f + 1];
    float sa[NS_MAX], sb[NS_MAX];
    limited_slope(ap, a, b, ns, spacing, m_tol, mass, phi_nodes, sa);
    limited_slope(a, b, bn, ns, spacing, m_tol, mass, phi_nodes, sb);
    float left[NS_MAX], right[NS_MAX];
    for (uint k = 0; k < NS_MAX; k++) { left[k] = 0.0; right[k] = 0.0; }
    for (uint k = 0; k < ns; k++)
    {
        left[k] = phi_nodes[ns * a + k] + 0.5 * spacing * sa[k];
        right[k] = phi_nodes[ns * b + k] - 0.5 * spacing * sb[k];
    }
    float forcing = 0.5 * (c_seg * grad_pk[2 * a + axis] + c_seg * grad_pk[2 * b + axis]);
    float fl[NS_MAX], fr[NS_MAX], al, ar;
    advection(left, forcing, ns, sizes, fl, al);
    advection(right, forcing, ns, sizes, fr, ar);
    float speed = max(al, ar);
    float face_diffusion = 0.5 * (diffusivity[a] + diffusivity[b]);
    float rho_a = safe_divide(mass[a], nodal_volume[a]), rho_b = safe_divide(mass[b], nodal_volume[b]);
    float face_density = safe_divide(2.0 * rho_a * rho_b, rho_a + rho_b);
    bool blocked = !(mass[a] > m_tol) || !(mass[b] > m_tol) || wall;
    for (uint k = 0; k < ns; k++)
    {
        float flux = 0.5 * (fl[k] + fr[k]) - 0.5 * speed * (right[k] - left[k]);
        flux -= face_diffusion * (phi_nodes[ns * b + k] - phi_nodes[ns * a + k]) / spacing;
        transfer[ns * f + k] = blocked ? 0.0 : dt * area * face_density * flux;
    }
}

// ---- Inventory-aware particle exchange (integrators._particle_exchanges) -------

[shader("compute")]
[numthreads(256, 1, 1)]
void seg_outgoing(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform uint ns,
                  StructuredBuffer<uint> node_faces, StructuredBuffer<float> mass,
                  RWStructuredBuffer<float> phi_nodes, StructuredBuffer<float> transfer,
                  RWStructuredBuffer<float> ratio)
{
    // outgoing / inventory per node: positive parts where the node is a left
    // endpoint, negative parts where it is a right endpoint.
    // node_faces: up to 4 signed entries per node, (face + 1) or -(face + 1) packed as uint.
    uint j = tid.x;
    if (j >= node_count) return;
    for (uint k = 0; k < ns; k++)
    {
        float out = 0.0;
        for (uint s = 0; s < 4; s++)
        {
            uint entry = node_faces[4 * j + s];
            if (entry == 0) continue;
            int signed_entry = asint(entry);
            uint f = uint(abs(signed_entry)) - 1;
            float t = transfer[ns * f + k];
            out += signed_entry > 0 ? max(t, 0.0) : max(-t, 0.0);
        }
        ratio[ns * j + k] = safe_divide(out, mass[j] * phi_nodes[ns * j + k]);
    }
}

void corners_of(uint p, uint nx, uint periodic, StructuredBuffer<uint> keys, StructuredBuffer<float> local_x,
                StructuredBuffer<float> local_y, StructuredBuffer<float> mass, float m_tol, out uint corners[4],
                out float w[4])
{
    uint columns = nx - 1;
    uint cx = keys[p] % columns, cy = keys[p] / columns;
    float fx = local_x[p], fy = local_y[p];
    corners[0] = canonical(cx, cy, nx, periodic);
    corners[1] = canonical(cx + 1, cy, nx, periodic);
    corners[2] = canonical(cx + 1, cy + 1, nx, periodic);
    corners[3] = canonical(cx, cy + 1, nx, periodic);
    float raw[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
    for (uint c = 0; c < 4; c++) w[c] = mass[corners[c]] > m_tol ? raw[c] : 0.0;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void seg_demand(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform uint nx, uniform uint periodic,
                uniform uint ns, uniform float m_tol, StructuredBuffer<uint> keys,
                StructuredBuffer<float> local_x, StructuredBuffer<float> local_y, StructuredBuffer<float> mass,
                StructuredBuffer<float> phi, StructuredBuffer<float> ratio, RWStructuredBuffer<float> limit)
{
    // limit = max(phi, 0) / demand, demand = phi * interpolated outgoing ratio.
    uint p = tid.x;
    if (p >= count) return;
    uint corners[4];
    float w[4];
    corners_of(p, nx, periodic, keys, local_x, local_y, mass, m_tol, corners, w);
    for (uint k = 0; k < ns; k++)
    {
        float interpolated = 0.0;
        for (uint c = 0; c < 4; c++) interpolated += w[c] * ratio[ns * corners[c] + k];
        float value = phi[ns * p + k];
        float demand = value * interpolated;
        limit[ns * p + k] = demand > 0.0 ? max(value, 0.0) / demand : 1.0;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void seg_node_limit(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform uint nx, uniform uint ny,
                    uniform uint periodic, uniform uint ns, uniform float m_tol,
                    StructuredBuffer<uint> cell_offsets, StructuredBuffer<uint> cell_counts,
                    StructuredBuffer<uint> permutation, StructuredBuffer<float> local_x,
                    StructuredBuffer<float> local_y, StructuredBuffer<float> mass, StructuredBuffer<float> phi,
                    RWStructuredBuffer<float> phi_nodes, StructuredBuffer<float> limit,
                    RWStructuredBuffer<float> node_limit)
{
    // Smallest particle limit among particles that hold species k and weight this node.
    uint j = tid.x;
    if (j >= node_count) return;
    uint node = canonical(j % nx, j / nx, nx, periodic);
    float smallest[NS_MAX];
    for (uint k = 0; k < NS_MAX; k++) smallest[k] = 1.0;
    if (mass[node] > m_tol)
    {
        int x = int(node % nx), y = int(node / nx);
        int columns = int(nx) - 1;
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
                uint begin = cell_offsets[cell], end = begin + cell_counts[cell];
                for (uint e = begin; e < end; e++)
                {
                    uint p = permutation[e];
                    float fx = local_x[p], fy = local_y[p];
                    float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
                    bool touches = false;
                    for (uint c = 0; c < 4; c++) touches = touches || (corners[c] == node && w[c] > 0.0);
                    if (!touches) continue;
                    for (uint k = 0; k < ns; k++)
                        if (phi[ns * p + k] > 0.0) smallest[k] = min(smallest[k], limit[ns * p + k]);
                }
            }
        }
    }
    for (uint k = 0; k < ns; k++)
        node_limit[ns * j + k] = mass[j] * phi_nodes[ns * j + k] > 0.0 ? smallest[k] : 0.0;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void seg_limit_faces(uint3 tid : SV_DispatchThreadID, uniform uint face_count, uniform uint ns,
                     StructuredBuffer<uint> faces, StructuredBuffer<float> node_limit,
                     RWStructuredBuffer<float> transfer)
{
    // One limiter per face, then balance positive and negative totals.
    uint f = tid.x;
    if (f >= face_count) return;
    uint a = faces[6 * f], b = faces[6 * f + 1];
    float theta = 1.0, positive = 0.0, negative = 0.0;
    for (uint k = 0; k < ns; k++)
    {
        float t = transfer[ns * f + k];
        if (t > 0.0) theta = min(theta, node_limit[ns * a + k]);
        if (t < 0.0) theta = min(theta, node_limit[ns * b + k]);
    }
    if (theta < 1.0) theta *= 1.0 - 32.0 * FLT_EPS;
    for (uint k = 0; k < ns; k++)
    {
        float t = transfer[ns * f + k] * theta;
        positive += max(t, 0.0);
        negative += max(-t, 0.0);
    }
    float exchange = min(positive, negative);
    float up = safe_divide(exchange, positive), down = safe_divide(exchange, negative);
    for (uint k = 0; k < ns; k++)
    {
        float t = transfer[ns * f + k] * theta;
        transfer[ns * f + k] = max(t, 0.0) * up - max(-t, 0.0) * down;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void seg_coefficients(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform uint nx,
                      uniform uint periodic, uniform uint ns, StructuredBuffer<uint> node_faces,
                      StructuredBuffer<float> mass, RWStructuredBuffer<float> phi_nodes,
                      StructuredBuffer<float> transfer, RWStructuredBuffer<float> ratio,
                      RWStructuredBuffer<float> coefficients, RWStructuredBuffer<float> change)
{
    // Outgoing ratio with the limited transfers, exchange coefficients
    // [target][donor] and the grid concentration change (for G.dphi).
    uint j = tid.x;
    if (j >= node_count) return;
    uint node = canonical(j % nx, j / nx, nx, periodic);
    float inventory[NS_MAX];
    for (uint k = 0; k < NS_MAX; k++) inventory[k] = k < ns ? mass[node] * phi_nodes[ns * node + k] : 0.0;
    float out[NS_MAX], net[NS_MAX];
    float coefficient[NS_MAX * NS_MAX];
    for (uint k = 0; k < NS_MAX; k++) { out[k] = 0.0; net[k] = 0.0; }
    for (uint k = 0; k < NS_MAX * NS_MAX; k++) coefficient[k] = 0.0;
    for (uint s = 0; s < 4; s++)
    {
        uint entry = node_faces[4 * node + s];
        if (entry == 0) continue;
        int signed_entry = asint(entry);
        bool is_left = signed_entry > 0;
        uint f = uint(abs(signed_entry)) - 1;
        float positive_sum = 0.0, negative_sum = 0.0;
        for (uint k = 0; k < ns; k++)
        {
            float t = transfer[ns * f + k];
            positive_sum += max(t, 0.0);
            negative_sum += max(-t, 0.0);
        }
        for (uint donor = 0; donor < ns; donor++)
        {
            float t = transfer[ns * f + donor];
            float given = is_left ? max(t, 0.0) : max(-t, 0.0);
            out[donor] += given;
            net[donor] += is_left ? -t : t;
            float share_scale = safe_divide(given, inventory[donor]);
            for (uint target = 0; target < ns; target++)
            {
                float u = transfer[ns * f + target];
                float share = is_left ? safe_divide(max(-u, 0.0), negative_sum) : safe_divide(max(u, 0.0), positive_sum);
                coefficient[NS_MAX * target + donor] += share_scale * share;
            }
        }
    }
    for (uint k = 0; k < ns; k++)
    {
        ratio[ns * j + k] = safe_divide(out[k], inventory[k]);
        change[ns * j + k] = safe_divide(net[k], mass[node]);
        for (uint d = 0; d < ns; d++) coefficients[ns * ns * j + ns * k + d] = coefficient[NS_MAX * k + d];
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void seg_particles(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform uint nx, uniform uint periodic,
                   uniform uint ns, uniform float m_tol, uniform uint normalise, StructuredBuffer<uint> keys,
                   StructuredBuffer<float> local_x, StructuredBuffer<float> local_y, StructuredBuffer<float> mass,
                   StructuredBuffer<float> ratio, StructuredBuffer<float> coefficients,
                   RWStructuredBuffer<float> phi, RWStructuredBuffer<uint> status)
{
    // phi += incoming - outgoing; each particle's mass is unchanged.
    uint p = tid.x;
    if (p >= count) return;
    uint corners[4];
    float w[4];
    corners_of(p, nx, periodic, keys, local_x, local_y, mass, m_tol, corners, w);
    float old[NS_MAX], next[NS_MAX];
    for (uint k = 0; k < NS_MAX; k++) { old[k] = k < ns ? phi[ns * p + k] : 0.0; next[k] = old[k]; }
    for (uint target = 0; target < ns; target++)
    {
        float out = 0.0, incoming = 0.0;
        for (uint c = 0; c < 4; c++)
        {
            out += w[c] * ratio[ns * corners[c] + target];
            for (uint donor = 0; donor < ns; donor++)
                incoming += old[donor] * w[c] * coefficients[ns * ns * corners[c] + ns * target + donor];
        }
        next[target] = old[target] + incoming - old[target] * out;
    }
    if (normalise != 0)
    {
        float total = 0.0;
        for (uint k = 0; k < ns; k++)
        {
            next[k] = clamp(next[k], 0.0, 1.0);
            total += next[k];
        }
        if (total == 0.0)
        {
            uint ignored;
            InterlockedOr(status[0], 4, ignored);
            return;
        }
        for (uint k = 0; k < ns; k++) next[k] /= total;
    }
    bool finite = true;
    for (uint k = 0; k < ns; k++) finite = finite && isfinite(next[k]);
    if (!finite)
    {
        uint ignored;
        InterlockedOr(status[0], 2, ignored);
        return;
    }
    for (uint k = 0; k < ns; k++) phi[ns * p + k] = next[k];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void seg_copy(uint3 tid : SV_DispatchThreadID, uniform uint count, StructuredBuffer<float> source,
              RWStructuredBuffer<float> target)
{
    if (tid.x < count) target[tid.x] = source[tid.x];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void seg_species(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform uint ns,
                 StructuredBuffer<float> mass, StructuredBuffer<float> phi, RWStructuredBuffer<float> columns)
{
    // Column k holds m phi_k per particle, for the species inventory.
    uint p = tid.x;
    if (p >= count) return;
    for (uint k = 0; k < ns; k++) columns[k * count + p] = mass[p] * phi[ns * p + k];
}
