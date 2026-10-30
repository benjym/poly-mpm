// Critical-state CIDR step (cidr.py) on Metal, in FP32.
//
// Rows are the canonical stress rows of the compact operator; a row is active
// iff it holds grain volume. Per row: the strain rate D = B u / Vs, the exact
// proximal step of psi (jam test and local Newton), the Moreau-envelope
// gradient and generalised Hessian. Per velocity component: the Newton
// gradient, a Jacobi diagonal and the matrix-free product with
// K = M/dt + B^T (W / Vs^2) H B. CG scalars stay on the device.

import mechanics;

static const float ROOT2 = 1.41421356237;
static const float PI = 3.14159265359;
static const uint SAMPLES = 128;

// Scalars (float buffer `scalars`).
static const uint RZ = 0;        // r . z
static const uint PQ = 1;        // p . K p
static const uint RR = 2;        // r . r
static const uint BB = 3;        // |rhs|^2
static const uint ALPHA = 4;
static const uint BETA = 5;
static const uint RZ_NEW = 6;
static const uint TOL2 = 7;      // (rtol |rhs|)^2

// Flags (uint buffer `flags`).
static const uint CG_DONE = 0;
static const uint STATUS = 1;    // bit 1: local Newton unconverged, bit 2: CG breakdown

struct Law
{
    float mu_s, mu_d, i0, k, a0;
    float lo, step;              // alpha table: ln I = lo + step * i
    uint count;
};

float mu_of(Law m, float I) { return m.mu_s + (m.mu_d - m.mu_s) * I / (m.i0 + I); }
float dmu_of(Law m, float I) { return (m.mu_d - m.mu_s) * m.i0 / ((m.i0 + I) * (m.i0 + I)); }

// alpha(I) by cubic Hermite in ln I; the slopes come from the ODE alpha + 5/2 I alpha' = 2 mu + 2 I mu'.
float alpha_of(Law m, StructuredBuffer<float> table, float I)
{
    float x = clamp((log(max(I, 1e-30)) - m.lo) / m.step, 0.0, float(m.count - 1) - 1e-3);
    uint i = uint(x);
    float t = x - float(i);
    float y0 = table[2 * i], y1 = table[2 * i + 2];
    float d0 = table[2 * i + 1] * m.step, d1 = table[2 * i + 3] * m.step;
    float t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * y0 + (t3 - 2 * t2 + t) * d0 + (-2 * t3 + 3 * t2) * y1 + (t3 - t2) * d1;
}

struct Functions { float I, a, da, b, db; };

Functions functions(Law m, StructuredBuffer<float> table, float s, float p)
{
    Functions f;
    f.I = m.k * s / sqrt(p);
    f.a = alpha_of(m, table, f.I);
    float mu = mu_of(m, f.I), dmu = dmu_of(m, f.I);
    f.da = f.I > 1e-10 ? (2 * mu + 2 * f.I * dmu - f.a) / (2.5 * f.I) : 0.0;
    f.b = 2 * (f.a - mu);
    f.db = 2 * (f.da - dmu);
    return f;
}

struct Psi { float psi, ps, pd, pss, psd, pdd, p; };

float gp_of(Law m, StructuredBuffer<float> table, float s, float p, float C)
{
    return 2 * s * (functions(m, table, s, p).b - 2 * p / C);
}

// psi(s, delta) and its derivatives at packing modulus C (CriticalStateCIDR.evaluate_C).
Psi evaluate(Law m, StructuredBuffer<float> table, float s, float delta, float C)
{
    Psi out;
    float cap = m.a0 * C;
    float tiny = 1e-30 * cap;
    bool at_cap = gp_of(m, table, s, cap, C) - delta >= 0;
    if (!at_cap && gp_of(m, table, s, tiny, C) - delta <= 0)
    {
        out.psi = out.ps = out.pd = out.pss = out.psd = out.pdd = out.p = 0;
        return out;
    }
    float p = cap;
    if (!at_cap)
    {
        float lo = tiny, hi = cap;
        for (uint i = 0; i < 48; i++)
        {
            p = 0.5 * (lo + hi);
            if (gp_of(m, table, s, p, C) - delta > 0) lo = p; else hi = p;
        }
        p = 0.5 * (lo + hi);
    }
    Functions f = functions(m, table, s, p);
    float mu = mu_of(m, f.I);
    float g = 2 * p * s * ((5 * f.a - 4 * mu) / 3 - p / C);
    float Y = f.a * p - p * p / C;
    float g_ss = 2 * f.da * p * f.I / max(s, 1e-30);
    out.psi = g - p * delta;
    out.ps = 2 * Y;
    out.pd = -p;
    out.p = p;
    if (at_cap)
    {
        out.pss = g_ss;
        out.psd = 0;
        out.pdd = 0;
    }
    else
    {
        float g_sp = 2 * (f.a - 2 * p / C) - f.da * f.I;
        float g_pp = 2 * s * (-2 / C - f.I * f.db / (2 * p));
        out.pss = g_ss - g_sp * g_sp / g_pp;
        out.psd = g_sp / g_pp;
        out.pdd = -1 / g_pp;
    }
    return out;
}

// sup_{0 <= p <= a0 C} [2 s Y0(p) - p delta], the limit of psi(t s, t delta) / t.
float static_support(Law m, float s, float delta, float C)
{
    float p = s > 0 ? C * (2 * s * m.a0 - delta) / (4 * s) : (delta < 0 ? 1e30 : 0.0);
    p = clamp(p, 0.0, m.a0 * C);
    return 2 * s * (m.a0 * p - p * p / C) - p * delta;
}

float support_value(Law m, float th, float f0, float fv, float C)
{
    return fv * cos(th) + f0 * sin(th) - static_support(m, cos(th) / ROOT2, ROOT2 * sin(th), C);
}

// F = psi(sigma / sqrt 2, sqrt 2 d0) + r/2 [(sigma - nv)^2 + (d0 - q0)^2] with gradient and Hessian.
struct Local { float F; float2 g; float2x2 H; Psi psi; };

Local local_terms(Law m, StructuredBuffer<float> table, float sigma, float d0, float nv, float q0, float r, float C)
{
    Local out;
    out.psi = evaluate(m, table, sigma / ROOT2, ROOT2 * d0, C);
    Psi e = out.psi;
    out.F = e.psi + 0.5 * r * ((sigma - nv) * (sigma - nv) + (d0 - q0) * (d0 - q0));
    out.g = float2(e.ps / ROOT2 + r * (sigma - nv), ROOT2 * e.pd + r * (d0 - q0));
    out.H = float2x2(0.5 * e.pss + r, e.psd, e.psd, 2 * e.pdd + r);
    return out;
}

// Exact proximal step argmin psi(d) + r/2 |d - q|^2 and the envelope's gradient and Hessian.
struct Envelope { float3 grad; float3x3 H; uint jammed; uint failed; float p; };

Envelope envelope(Law m, StructuredBuffer<float> table, float3 q, float r, float C)
{
    Envelope out;
    out.failed = 0;
    float q0 = q.x;
    float2 qv = float2(q.y, q.z);
    float nv = length(qv);
    float f0 = r * q0, fv = r * nv;
    // Jam test: max over unit e of f.e - h(e), on samples then by ternary search.
    uint best = 0;
    float bestv = -1e38;
    for (uint i = 0; i < SAMPLES; i++)
    {
        float th = -0.5 * PI + PI * float(i) / float(SAMPLES - 1);
        float v = support_value(m, th, f0, fv, C);
        if (v > bestv) { bestv = v; best = i; }
    }
    float stepth = PI / float(SAMPLES - 1);
    float lo = max(-0.5 * PI, -0.5 * PI + stepth * (float(best) - 1));
    float hi = min(0.5 * PI, -0.5 * PI + stepth * (float(best) + 1));
    for (uint i = 0; i < 32; i++)
    {
        float m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
        if (support_value(m, m1, f0, fv, C) < support_value(m, m2, f0, fv, C)) lo = m1; else hi = m2;
    }
    float theta = 0.5 * (lo + hi);
    float margin = support_value(m, theta, f0, fv, C);
    float scale = r * length(q);
    float3x3 identity = float3x3(1, 0, 0, 0, 1, 0, 0, 0, 1);
    if (margin <= 1e-6 * max(scale, 1e-30))
    {
        out.jammed = 1;
        out.grad = r * q;
        out.H = r * identity;
        out.p = 0;
        return out;
    }
    out.jammed = 0;
    float P = m.a0 * C;
    if (nv <= 1e-7 * max(abs(q0), 1e-30))
    {
        // Pure volumetric: soft threshold, psi(0, delta) = P max(-delta, 0).
        float d0 = q0 >= 0 ? q0 : min(q0 + ROOT2 * P / r, 0.0);
        out.grad = float3(r * (q0 - d0), r * q.y, r * q.z);
        out.H = r * identity;
        if (d0 == q0) out.H[0][0] = 0;
        out.p = d0 < 0 ? P : 0;
        return out;
    }
    // Damped Newton in (sigma, d0) from (|q_dev|, q0) or from just outside the cone.
    float sigma = nv, d0 = q0;
    float exit_ = margin / r;
    float se = exit_ * cos(theta), de = exit_ * sin(theta);
    if (se > 0 && local_terms(m, table, se, de, nv, q0, r, C).F < local_terms(m, table, nv, q0, nv, q0, r, C).F)
    {
        sigma = se;
        d0 = de;
    }
    Local L = local_terms(m, table, sigma, d0, nv, q0, r, C);
    bool converged = false;
    for (uint k = 0; k < 60; k++)
    {
        float det = L.H[0][0] * L.H[1][1] - L.H[0][1] * L.H[1][0];
        float2 step = -float2(L.H[1][1] * L.g.x - L.H[0][1] * L.g.y, L.H[0][0] * L.g.y - L.H[1][0] * L.g.x) / det;
        float slope = dot(L.g, step);
        // Stop on the gradient, which FP32 resolves to its relative precision against the
        // stress scale r |q|; the decrement it cannot resolve near the minimum.
        if (!(length(L.g) > 1e-6 * r * length(q)) || !(-slope > 0))
        {
            converged = true;
            break;
        }
        float t = 1.0;
        bool moved = false;
        for (uint j = 0; j < 40; j++)
        {
            float st = sigma + t * step.x;
            if (st > 0)
            {
                Local T = local_terms(m, table, st, d0 + t * step.y, nv, q0, r, C);
                // Sufficient decrease, or (below FP32 resolution of F) a smaller gradient.
                if (T.F <= L.F + 1e-4 * t * slope || length(T.g) < (1 - 1e-4 * t) * length(L.g))
                {
                    sigma = st;
                    d0 += t * step.y;
                    L = T;
                    moved = true;
                    break;
                }
            }
            t *= 0.5;
        }
        if (!moved)
        {
            converged = true;   // no decrease representable in FP32
            break;
        }
    }
    out.failed = converged ? 0 : 1;
    out.p = L.psi.p;
    // Gradient: grad psi(d), which equals r (q - d) at the solution without its cancellation.
    float2 n = qv / nv;
    Psi e = L.psi;
    out.grad = float3(ROOT2 * e.pd, e.ps / ROOT2 * n.x, e.ps / ROOT2 * n.y);
    // Hessian of psi at d in the components, then r H (H + r I)^-1 without cancellation.
    float3x3 Hp;
    float transverse = e.ps / (ROOT2 * sigma);
    float2x2 outer = float2x2(n.x * n.x, n.x * n.y, n.y * n.x, n.y * n.y);
    float2x2 dev = 0.5 * e.pss * outer + transverse * (float2x2(1, 0, 0, 1) - outer);
    Hp = float3x3(2 * e.pdd, e.psd * n.x, e.psd * n.y,
                  e.psd * n.x, dev[0][0], dev[0][1],
                  e.psd * n.y, dev[1][0], dev[1][1]);
    float3x3 A = Hp + r * identity;
    float3x3 Ainv = inverse3(A);
    float3x3 M = r * mul(Hp, Ainv);
    out.H = 0.5 * (M + transpose(M));
    return out;
}

float3x3 inverse3(float3x3 a)
{
    float3x3 c;
    c[0][0] = a[1][1] * a[2][2] - a[1][2] * a[2][1];
    c[0][1] = a[0][2] * a[2][1] - a[0][1] * a[2][2];
    c[0][2] = a[0][1] * a[1][2] - a[0][2] * a[1][1];
    c[1][0] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
    c[1][1] = a[0][0] * a[2][2] - a[0][2] * a[2][0];
    c[1][2] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
    c[2][0] = a[1][0] * a[2][1] - a[1][1] * a[2][0];
    c[2][1] = a[0][1] * a[2][0] - a[0][0] * a[2][1];
    c[2][2] = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    float det = a[0][0] * c[0][0] + a[0][1] * c[1][0] + a[0][2] * c[2][0];
    return c / det;
}

float3 row_strain(uint a, StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                  StructuredBuffer<float> velocity)
{
    float3 w = float3(0.0, 0.0, 0.0);
    for (uint s = 0; s < 9; s++)
    {
        uint j = neighbors[9 * a + s];
        if (j == INVALID) continue;
        w += mul(row_block(a, row_slot(a, s, neighbors), gradients), float2(velocity[2 * j], velocity[2 * j + 1]));
    }
    return w;
}

float3 row_strain_rw(uint a, StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                     RWStructuredBuffer<float> velocity)
{
    float3 w = float3(0.0, 0.0, 0.0);
    for (uint s = 0; s < 9; s++)
    {
        uint j = neighbors[9 * a + s];
        if (j == INVALID) continue;
        w += mul(row_block(a, row_slot(a, s, neighbors), gradients), float2(velocity[2 * j], velocity[2 * j + 1]));
    }
    return w;
}

// ---- set-up -------------------------------------------------------------------------------

[shader("compute")]
[numthreads(64, 1, 1)]
void cidr_bulk_rows(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform uint nx,
                    uniform uint ny, uniform uint periodic, StructuredBuffer<uint> row_nodes,
                    StructuredBuffer<uint> cell_offsets, StructuredBuffer<uint> cell_counts,
                    StructuredBuffer<uint> permutation, StructuredBuffer<float> local_x,
                    StructuredBuffer<float> local_y, StructuredBuffer<float> volume,
                    RWStructuredBuffer<float> bulk)
{
    // Bulk volume sum_p V_p N_a(x_p) per stress row, the gather of assemble_gradients.
    uint a = tid.x;
    if (a >= row_count) return;
    uint node = row_nodes[a];
    int x = int(node % nx);
    int y = int(node / nx);
    int columns = int(nx) - 1;
    float total = 0.0;
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
                float fx = local_x[p];
                float fy = local_y[p];
                float w[4] = { (1 - fx) * (1 - fy), fx * (1 - fy), fx * fy, (1 - fx) * fy };
                float weight = 0.0;
                for (uint c = 0; c < 4; c++) weight += corners[c] == node ? w[c] : 0.0;
                total += volume[p] * weight;
            }
        }
    }
    bulk[a] = total;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_rows_setup(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform float phi_min,
                     uniform float phi_max, uniform float modulus_scale, StructuredBuffer<float> occupied,
                     StructuredBuffer<float> bulk, StructuredBuffer<float> initial,
                     RWStructuredBuffer<float> phi, RWStructuredBuffer<float> modulus,
                     RWStructuredBuffer<float> lam, RWStructuredBuffer<float> lam_abs)
{
    // Packing and C (cidr.packing_modulus) per row; the multiplier starts from the previous
    // contact: lam = -contact phi.
    uint a = tid.x;
    if (a >= row_count) return;
    float solid = occupied[a];
    bool active = solid > 0.0 && bulk[a] > 0.0;
    float f = active ? solid / bulk[a] : 0.0;
    phi[a] = f;
    float q = clamp((f - phi_min) / (phi_max - phi_min), 1e-7, 1 - 1e-7);
    modulus[a] = active ? modulus_scale * q / (1 - q) : 0.0;
    float m = 0.0;
    for (uint k = 0; k < 3; k++)
    {
        float v = active ? -initial[3 * a + k] * f : 0.0;
        lam[3 * a + k] = v;
        m = max(m, abs(v));
    }
    lam_abs[a] = m;
}

// ---- per-row envelope -----------------------------------------------------------------------

[shader("compute")]
[numthreads(64, 1, 1)]
void cidr_envelope(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform float gamma,
                   uniform float mu_s, uniform float mu_d, uniform float i0, uniform float k,
                   uniform float a0, uniform float table_lo, uniform float table_step, uniform uint table_count,
                   StructuredBuffer<float> table, StructuredBuffer<uint> neighbors,
                   StructuredBuffer<float> gradients, StructuredBuffer<float> velocity,
                   StructuredBuffer<float> occupied, StructuredBuffer<float> bulk,
                   StructuredBuffer<float> modulus, StructuredBuffer<float> lam,
                   RWStructuredBuffer<float> grad, RWStructuredBuffer<float> force,
                   RWStructuredBuffer<float> hessian, RWStructuredBuffer<uint> jammed,
                   RWStructuredBuffer<float> pressure, RWStructuredBuffer<uint> flags)
{
    // x = B u / Vs + gamma lam; envelope gradient g (the new multiplier), force g Vb / Vs
    // for B^T, and the Hessian block scaled by Vb / Vs^2 for the Newton matrix.
    uint a = tid.x;
    if (a >= row_count) return;
    float solid = occupied[a];
    float C = modulus[a];
    for (uint i = 0; i < 3; i++) { grad[3 * a + i] = 0; force[3 * a + i] = 0; }
    for (uint i = 0; i < 6; i++) hessian[6 * a + i] = 0;
    jammed[a] = 0;
    pressure[a] = 0;
    if (!(solid > 0.0) || !(C > 0.0)) return;
    Law m = { mu_s, mu_d, i0, k, a0, table_lo, table_step, table_count };
    float3 x = row_strain(a, neighbors, gradients, velocity) / solid
               + gamma * float3(lam[3 * a], lam[3 * a + 1], lam[3 * a + 2]);
    Envelope e = envelope(m, table, x, 1.0 / gamma, C);
    float ratio = bulk[a] / solid;
    float weight = ratio / solid;
    for (uint i = 0; i < 3; i++)
    {
        grad[3 * a + i] = e.grad[i];
        force[3 * a + i] = ratio * e.grad[i];
    }
    hessian[6 * a + 0] = weight * e.H[0][0];
    hessian[6 * a + 1] = weight * e.H[0][1];
    hessian[6 * a + 2] = weight * e.H[0][2];
    hessian[6 * a + 3] = weight * e.H[1][1];
    hessian[6 * a + 4] = weight * e.H[1][2];
    hessian[6 * a + 5] = weight * e.H[2][2];
    jammed[a] = e.jammed;
    pressure[a] = e.p;
    if (e.failed != 0)
    {
        uint ignored;
        InterlockedOr(flags[STATUS], 2, ignored);
    }
}

float3x3 row_hessian(RWStructuredBuffer<float> h, uint a)
{
    return float3x3(h[6 * a], h[6 * a + 1], h[6 * a + 2],
                    h[6 * a + 1], h[6 * a + 3], h[6 * a + 4],
                    h[6 * a + 2], h[6 * a + 4], h[6 * a + 5]);
}

float2 node_force(uint j, StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
                  StructuredBuffer<float> gradients, RWStructuredBuffer<float> rows)
{
    float2 f = float2(0.0, 0.0);
    for (uint i = 0; i < 9; i++)
    {
        uint a = transpose_rows[9 * j + i];
        if (a == INVALID) break;
        f += mul(float3(rows[3 * a], rows[3 * a + 1], rows[3 * a + 2]), row_block(a, transpose_slots[9 * j + i], gradients));
    }
    return f;
}

float mass_dt(StructuredBuffer<float> inverse_mass, uint i, float dt)
{
    float w = inverse_mass[i];
    return w > 0.0 ? 1.0 / (w * dt) : 0.0;
}

// ---- Newton gradient, preconditioner and products -----------------------------------------

groupshared float sum_a[256];
groupshared float sum_b[256];

// Group sums of two values; thread 0 of each group writes them to partials[2 g], [2 g + 1].
void group_sums(uint lane, uint group, float a, float b, RWStructuredBuffer<float> partials)
{
    sum_a[lane] = a;
    sum_b[lane] = b;
    GroupMemoryBarrierWithGroupSync();
    for (uint offset = 128; offset > 0; offset /= 2)
    {
        if (lane < offset)
        {
            sum_a[lane] += sum_a[lane + offset];
            sum_b[lane] += sum_b[lane + offset];
        }
        GroupMemoryBarrierWithGroupSync();
    }
    if (lane == 0)
    {
        partials[2 * group] = sum_a[0];
        partials[2 * group + 1] = sum_b[0];
    }
}

// The 2x2 block preconditioner of node j: (d00, d01, d11), identity on fixed components.
float2 precondition(RWStructuredBuffer<float> block, uint j, float2 r)
{
    float a = block[3 * j], b = block[3 * j + 1], c = block[3 * j + 2];
    float det = a * c - b * b;
    return float2(c * r.x - b * r.y, a * r.y - b * r.x) / det;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_gradient(uint3 tid : SV_DispatchThreadID, uint3 lane : SV_GroupThreadID, uint3 group : SV_GroupID,
                   uniform uint node_count, uniform float dt,
                   StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
                   StructuredBuffer<float> gradients, StructuredBuffer<float> inverse_mass,
                   StructuredBuffer<float> velocity, StructuredBuffer<float> trial,
                   RWStructuredBuffer<float> force, RWStructuredBuffer<float> hessian,
                   RWStructuredBuffer<float> gradient, RWStructuredBuffer<float> block,
                   RWStructuredBuffer<float> partials)
{
    // G = M/dt (u - u*) + B^T force on free components, the 2x2 block of K per node, and the
    // CG start: x = 0, r = -G, z = P^-1 r, p = z; partials of (r . z, G . G).
    uint j = tid.x;
    float rz = 0.0, gg = 0.0;
    if (j < node_count)
    {
        float2 f = node_force(j, transpose_rows, transpose_slots, gradients, force);
        float d00 = 0.0, d01 = 0.0, d11 = 0.0;
        for (uint i = 0; i < 9; i++)
        {
            uint a = transpose_rows[9 * j + i];
            if (a == INVALID) break;
            float3x2 Bs = row_block(a, transpose_slots[9 * j + i], gradients);
            float3x3 H = row_hessian(hessian, a);
            float3 c0 = float3(Bs[0][0], Bs[1][0], Bs[2][0]);
            float3 c1 = float3(Bs[0][1], Bs[1][1], Bs[2][1]);
            float3 h0 = mul(H, c0), h1 = mul(H, c1);
            d00 += dot(c0, h0);
            d01 += dot(c0, h1);
            d11 += dot(c1, h1);
        }
        bool free0 = inverse_mass[2 * j] > 0.0, free1 = inverse_mass[2 * j + 1] > 0.0;
        float2 g = float2(free0 ? mass_dt(inverse_mass, 2 * j, dt) * (velocity[2 * j] - trial[2 * j]) + f.x : 0.0,
                          free1 ? mass_dt(inverse_mass, 2 * j + 1, dt) * (velocity[2 * j + 1] - trial[2 * j + 1]) + f.y : 0.0);
        block[3 * j] = free0 ? mass_dt(inverse_mass, 2 * j, dt) + d00 : 1.0;
        block[3 * j + 1] = free0 && free1 ? d01 : 0.0;
        block[3 * j + 2] = free1 ? mass_dt(inverse_mass, 2 * j + 1, dt) + d11 : 1.0;
        gradient[2 * j] = g.x;
        gradient[2 * j + 1] = g.y;
        float2 r = -g;
        float2 z = precondition(block, j, r);
        rz = dot(r, z);
        gg = dot(g, g);
    }
    group_sums(lane.x, group.x, rz, gg, partials);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_cg_start(uint3 tid : SV_DispatchThreadID, uniform uint node_count, StructuredBuffer<float> gradient,
                   RWStructuredBuffer<float> block, RWStructuredBuffer<float> x, RWStructuredBuffer<float> r,
                   RWStructuredBuffer<float> z, RWStructuredBuffer<float> p)
{
    uint j = tid.x;
    if (j >= node_count) return;
    float2 rj = -float2(gradient[2 * j], gradient[2 * j + 1]);
    float2 zj = precondition(block, j, rj);
    for (uint c = 0; c < 2; c++)
    {
        x[2 * j + c] = 0;
        r[2 * j + c] = rj[c];
        z[2 * j + c] = zj[c];
        p[2 * j + c] = zj[c];
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_rows_product(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform uint gate,
                       StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                       StructuredBuffer<float> source, RWStructuredBuffer<float> hessian,
                       RWStructuredBuffer<float> output, RWStructuredBuffer<uint> flags)
{
    // output = (Vb / Vs^2) H B source per row (gated off once CG has converged).
    uint a = tid.x;
    if (a >= row_count || (gate != 0 && flags[CG_DONE] != 0)) return;
    float3 w = mul(row_hessian(hessian, a), row_strain(a, neighbors, gradients, source));
    for (uint k = 0; k < 3; k++) output[3 * a + k] = w[k];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_nodes_product(uint3 tid : SV_DispatchThreadID, uint3 lane : SV_GroupThreadID, uint3 group : SV_GroupID,
                        uniform uint node_count, uniform float dt,
                        StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
                        StructuredBuffer<float> gradients, StructuredBuffer<float> inverse_mass,
                        StructuredBuffer<float> source, RWStructuredBuffer<float> rows,
                        RWStructuredBuffer<float> output, RWStructuredBuffer<float> partials,
                        RWStructuredBuffer<uint> flags)
{
    // output = K source = M/dt source + B^T rows on free components; partials of source . output.
    if (flags[CG_DONE] != 0) return;
    uint j = tid.x;
    float pq = 0.0;
    if (j < node_count)
    {
        float2 f = node_force(j, transpose_rows, transpose_slots, gradients, rows);
        for (uint c = 0; c < 2; c++)
        {
            uint i = 2 * j + c;
            bool free = inverse_mass[i] > 0.0;
            float v = free ? mass_dt(inverse_mass, i, dt) * source[i] + f[c] : 0.0;
            output[i] = v;
            pq += source[i] * v;
        }
    }
    group_sums(lane.x, group.x, pq, 0.0, partials);
}

// ---- conjugate gradients -------------------------------------------------------------------

[shader("compute")]
[numthreads(1, 1, 1)]
void cidr_cg_init(uint3 tid : SV_DispatchThreadID, uniform float rtol, uniform uint groups,
                  RWStructuredBuffer<float> partials, RWStructuredBuffer<float> scalars,
                  RWStructuredBuffer<uint> flags, RWStructuredBuffer<uint> iterations)
{
    float rz = 0.0, bb = 0.0;
    for (uint g = 0; g < groups; g++) { rz += partials[2 * g]; bb += partials[2 * g + 1]; }
    scalars[RZ] = rz;
    scalars[BB] = bb;
    scalars[TOL2] = rtol * rtol * bb;
    flags[CG_DONE] = bb > 0.0 ? 0 : 1;
    iterations[0] = 0;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void cidr_cg_alpha(uint3 tid : SV_DispatchThreadID, uniform uint groups, RWStructuredBuffer<float> partials,
                   RWStructuredBuffer<float> scalars, RWStructuredBuffer<uint> flags)
{
    if (flags[CG_DONE] != 0) return;
    float d = 0.0;
    for (uint g = 0; g < groups; g++) d += partials[2 * g];
    if (!(d > 0.0))
    {
        flags[CG_DONE] = 1;
        uint ignored;
        InterlockedOr(flags[STATUS], 4, ignored);
        scalars[ALPHA] = 0;
        return;
    }
    scalars[ALPHA] = scalars[RZ] / d;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_cg_update(uint3 tid : SV_DispatchThreadID, uint3 lane : SV_GroupThreadID, uint3 group : SV_GroupID,
                    uniform uint node_count, RWStructuredBuffer<float> block,
                    StructuredBuffer<float> q, RWStructuredBuffer<float> x, RWStructuredBuffer<float> r,
                    RWStructuredBuffer<float> z, StructuredBuffer<float> p, RWStructuredBuffer<float> scalars,
                    RWStructuredBuffer<float> partials, RWStructuredBuffer<uint> flags)
{
    // x += alpha p, r -= alpha q, z = P^-1 r; partials of (r . z, r . r).
    if (flags[CG_DONE] != 0) return;
    uint j = tid.x;
    float rz = 0.0, rr = 0.0;
    if (j < node_count)
    {
        float alpha = scalars[ALPHA];
        float2 rj;
        for (uint c = 0; c < 2; c++)
        {
            uint i = 2 * j + c;
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
            rj[c] = r[i];
        }
        float2 zj = precondition(block, j, rj);
        z[2 * j] = zj.x;
        z[2 * j + 1] = zj.y;
        rz = dot(rj, zj);
        rr = dot(rj, rj);
    }
    group_sums(lane.x, group.x, rz, rr, partials);
}

[shader("compute")]
[numthreads(1, 1, 1)]
void cidr_cg_beta(uint3 tid : SV_DispatchThreadID, uniform uint limit, uniform uint groups,
                  RWStructuredBuffer<float> partials, RWStructuredBuffer<float> scalars,
                  RWStructuredBuffer<uint> flags, RWStructuredBuffer<uint> iterations)
{
    if (flags[CG_DONE] != 0) return;
    float rz = 0.0, rr = 0.0;
    for (uint g = 0; g < groups; g++) { rz += partials[2 * g]; rr += partials[2 * g + 1]; }
    iterations[0] += 1;
    if (rr <= scalars[TOL2] || iterations[0] >= limit)
    {
        flags[CG_DONE] = 1;
        return;
    }
    scalars[BETA] = rz / scalars[RZ];
    scalars[RZ] = rz;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_cg_direction(uint3 tid : SV_DispatchThreadID, uniform uint count, StructuredBuffer<float> z,
                       RWStructuredBuffer<float> p, RWStructuredBuffer<float> scalars,
                       RWStructuredBuffer<uint> flags)
{
    uint i = tid.x;
    if (i >= count || flags[CG_DONE] != 0) return;
    p[i] = z[i] + scalars[BETA] * p[i];
}

// ---- conjugate gradients in one workgroup ---------------------------------------------------

static const uint BLOCK = 1024;
groupshared float block_a[1024];
groupshared float block_b[1024];

float2 block_sums(uint lane, float a, float b)
{
    block_a[lane] = a;
    block_b[lane] = b;
    GroupMemoryBarrierWithGroupSync();
    for (uint offset = BLOCK / 2; offset > 0; offset /= 2)
    {
        if (lane < offset)
        {
            block_a[lane] += block_a[lane + offset];
            block_b[lane] += block_b[lane + offset];
        }
        GroupMemoryBarrierWithGroupSync();
    }
    float2 total = float2(block_a[0], block_b[0]);
    GroupMemoryBarrierWithGroupSync();
    return total;
}

[shader("compute")]
[numthreads(1024, 1, 1)]
void cidr_cg_block(uint3 lane3 : SV_GroupThreadID, uniform uint row_count, uniform uint node_count,
                   uniform float dt, uniform uint steps, uniform uint limit,
                   StructuredBuffer<uint> neighbors, StructuredBuffer<uint> transpose_rows,
                   StructuredBuffer<uint> transpose_slots, StructuredBuffer<float> gradients,
                   StructuredBuffer<float> inverse_mass, RWStructuredBuffer<float> hessian,
                   RWStructuredBuffer<float> block, RWStructuredBuffer<float> rows,
                   RWStructuredBuffer<float> x, RWStructuredBuffer<float> r, RWStructuredBuffer<float> z,
                   RWStructuredBuffer<float> p, RWStructuredBuffer<float> q,
                   RWStructuredBuffer<float> scalars, RWStructuredBuffer<uint> flags,
                   RWStructuredBuffer<uint> iterations)
{
    // `steps` CG iterations in one workgroup: barriers replace dispatch boundaries. The same
    // arithmetic as cidr_rows_product, cidr_nodes_product, cidr_cg_update and cidr_cg_direction.
    uint lane = lane3.x;
    for (uint step = 0; step < steps; step++)
    {
        AllMemoryBarrierWithGroupSync();
        if (flags[CG_DONE] != 0) return;
        for (uint a = lane; a < row_count; a += BLOCK)
        {
            float3 w = mul(row_hessian(hessian, a), row_strain_rw(a, neighbors, gradients, p));
            for (uint k = 0; k < 3; k++) rows[3 * a + k] = w[k];
        }
        AllMemoryBarrierWithGroupSync();
        float pq = 0.0;
        for (uint j = lane; j < node_count; j += BLOCK)
        {
            float2 f = node_force(j, transpose_rows, transpose_slots, gradients, rows);
            for (uint c = 0; c < 2; c++)
            {
                uint i = 2 * j + c;
                float v = inverse_mass[i] > 0.0 ? mass_dt(inverse_mass, i, dt) * p[i] + f[c] : 0.0;
                q[i] = v;
                pq += p[i] * v;
            }
        }
        float d = block_sums(lane, pq, 0.0).x;
        if (!(d > 0.0))
        {
            if (lane == 0)
            {
                flags[CG_DONE] = 1;
                uint ignored;
                InterlockedOr(flags[STATUS], 4, ignored);
            }
            return;
        }
        float alpha = scalars[RZ] / d;
        AllMemoryBarrierWithGroupSync();
        float rz = 0.0, rr = 0.0;
        for (uint j = lane; j < node_count; j += BLOCK)
        {
            float2 rj;
            for (uint c = 0; c < 2; c++)
            {
                uint i = 2 * j + c;
                x[i] += alpha * p[i];
                r[i] -= alpha * q[i];
                rj[c] = r[i];
            }
            float2 zj = precondition(block, j, rj);
            z[2 * j] = zj.x;
            z[2 * j + 1] = zj.y;
            rz += dot(rj, zj);
            rr += dot(rj, rj);
        }
        float2 sums = block_sums(lane, rz, rr);
        uint count = iterations[0] + 1;
        bool done = sums.y <= scalars[TOL2] || count >= limit;
        float beta = sums.x / scalars[RZ];
        AllMemoryBarrierWithGroupSync();
        if (lane == 0)
        {
            iterations[0] = count;
            scalars[RZ] = sums.x;
            if (done) flags[CG_DONE] = 1;
        }
        if (done) return;
        for (uint i = lane; i < 2 * node_count; i += BLOCK) p[i] = z[i] + beta * p[i];
    }
}

// ---- line search and outer step ------------------------------------------------------------

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_axpy(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform float t,
               StructuredBuffer<float> base, StructuredBuffer<float> direction, RWStructuredBuffer<float> output)
{
    uint i = tid.x;
    if (i >= count) return;
    output[i] = base[i] + t * direction[i];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_products(uint3 tid : SV_DispatchThreadID, uniform uint count, StructuredBuffer<float> a,
                   StructuredBuffer<float> b, RWStructuredBuffer<float> work)
{
    uint i = tid.x;
    if (i >= count) return;
    work[i] = a[i] * b[i];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_slope_nodes(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform float dt,
                      StructuredBuffer<float> inverse_mass, StructuredBuffer<float> velocity,
                      StructuredBuffer<float> trial, StructuredBuffer<float> direction,
                      RWStructuredBuffer<float> work)
{
    // M/dt (u - u*) . delta per component, for the line-search derivative.
    uint i = tid.x;
    if (i >= count) return;
    work[i] = mass_dt(inverse_mass, i, dt) * (velocity[i] - trial[i]) * direction[i];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_slope_rows(uint3 tid : SV_DispatchThreadID, uniform uint row_count, StructuredBuffer<uint> neighbors,
                     StructuredBuffer<float> gradients, StructuredBuffer<float> direction,
                     RWStructuredBuffer<float> force, RWStructuredBuffer<float> work)
{
    // force . B delta per row (force = g Vb / Vs), the stress part of the derivative.
    uint a = tid.x;
    if (a >= row_count) return;
    float3 w = row_strain(a, neighbors, gradients, direction);
    work[a] = force[3 * a] * w.x + force[3 * a + 1] * w.y + force[3 * a + 2] * w.z;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_change(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform float dt,
                 StructuredBuffer<float> inverse_mass, StructuredBuffer<float> velocity,
                 StructuredBuffer<float> previous, StructuredBuffer<float> trial,
                 RWStructuredBuffer<float> work, RWStructuredBuffer<float> work2, RWStructuredBuffer<float> work3)
{
    // Mass-weighted |u - u_prev|^2, |u|^2 and |u*|^2 per component.
    uint i = tid.x;
    if (i >= count) return;
    float m = mass_dt(inverse_mass, i, dt);
    float d = velocity[i] - previous[i];
    work[i] = m * d * d;
    work2[i] = m * velocity[i] * velocity[i];
    work3[i] = m * trial[i] * trial[i];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_contact(uint3 tid : SV_DispatchThreadID, uniform uint row_count, StructuredBuffer<float> phi,
                  StructuredBuffer<float> lam, RWStructuredBuffer<float> contact,
                  StructuredBuffer<uint> jammed, RWStructuredBuffer<float> work)
{
    // contact = -lam / phi, the intrinsic stress the other solvers return; work = jammed.
    uint a = tid.x;
    if (a >= row_count) return;
    float f = phi[a];
    for (uint k = 0; k < 3; k++) contact[3 * a + k] = f > 0.0 ? -lam[3 * a + k] / f : 0.0;
    work[a] = float(jammed[a]);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_unforced(uint3 tid : SV_DispatchThreadID, uniform uint node_count, uniform float dt,
                   uniform float body_x, uniform float body_y, StructuredBuffer<float> inverse_mass,
                   RWStructuredBuffer<float> velocity)
{
    // Warm start: the transferred velocity before this step's body force, on free components.
    uint j = tid.x;
    if (j >= node_count) return;
    if (inverse_mass[2 * j] > 0.0) velocity[2 * j] -= dt * body_x;
    if (inverse_mass[2 * j + 1] > 0.0) velocity[2 * j + 1] -= dt * body_y;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void cidr_pack(uint3 tid : SV_DispatchThreadID, StructuredBuffer<float> a, StructuredBuffer<float> b,
               RWStructuredBuffer<float> output)
{
    // Two reduction results into one buffer, so the host reads them with one synchronisation.
    output[0] = a[0];
    output[1] = b[0];
}

[shader("compute")]
[numthreads(256, 1, 1)]
void cidr_copy(uint3 tid : SV_DispatchThreadID, uniform uint count, StructuredBuffer<float> source,
               RWStructuredBuffer<float> output)
{
    uint i = tid.x;
    if (i >= count) return;
    output[i] = source[i];
}

// Test entry: the envelope at given rows of strain rates.
[shader("compute")]
[numthreads(64, 1, 1)]
void cidr_envelope_test(uint3 tid : SV_DispatchThreadID, uniform uint count, uniform float r,
                        uniform float mu_s, uniform float mu_d, uniform float i0, uniform float k,
                        uniform float a0, uniform float table_lo, uniform float table_step,
                        uniform uint table_count, StructuredBuffer<float> table,
                        StructuredBuffer<float> q, StructuredBuffer<float> modulus,
                        RWStructuredBuffer<float> grad, RWStructuredBuffer<float> hessian,
                        RWStructuredBuffer<uint> jammed)
{
    uint a = tid.x;
    if (a >= count) return;
    Law m = { mu_s, mu_d, i0, k, a0, table_lo, table_step, table_count };
    Envelope e = envelope(m, table, float3(q[3 * a], q[3 * a + 1], q[3 * a + 2]), r, modulus[a]);
    for (uint i = 0; i < 3; i++) grad[3 * a + i] = e.grad[i];
    for (uint i = 0; i < 3; i++)
        for (uint j = 0; j < 3; j++) hessian[9 * a + 3 * i + j] = e.H[i][j];
    jammed[a] = e.jammed;
}
