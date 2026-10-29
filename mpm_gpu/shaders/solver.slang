// Resident compact Coulomb solver (granular_mixed.solve_velocity) for linear
// Q1 nodal stresses. Every iterate, history, small factorisation and control
// decision stays on the device. The host encodes bounded batches; each kernel
// reads a device flag word and returns early when its gate is closed. Single
// thread control kernels own every scalar decision, in separate dispatches.

import mechanics;

// Device state record (uint).
static const uint FLAGS = 0;
static const uint RESULT = 1;            // 0 running, 1 converged, 2 failed
static const uint SWEEPS = 2;
static const uint COUNT = 3;             // Anderson history states held
static const uint HEAD = 4;              // next history slot
static const uint ACCELERATED = 5;
static const uint ANDERSON_REJECTED = 6;
static const uint NEWTON_STEPS = 7;
static const uint NEWTON_REJECTED = 8;
static const uint GMRES_ITERATIONS = 9;
static const uint KRYLOV = 10;           // Arnoldi steps in the current cycle
static const uint BREAKDOWNS = 11;
static const uint NEWTON_ATTEMPTS = 12;
static const uint GMRES_CONVERGED = 13;

// Flag bits.
static const uint RUNNING = 1;
static const uint ANDERSON = 2;
static const uint ACCEPTED = 4;
static const uint NEWTON = 8;
static const uint SEARCH = 16;
static const uint GMRES = 32;
static const uint UPDATE = 64;
static const uint KRYLOV_DONE = 128;
static const uint REQUESTED = 256;

// Scalar record (float).
static const uint MAGNITUDE = 0;
static const uint NORM = 1;
static const uint SCALE = 2;
static const uint ERROR = 3;

// Small Krylov record (float): Givens-reduced Hessenberg and right-hand side.
static const uint RESTART = 64;
static const uint KG = 0;
static const uint KCS = 65;
static const uint KSN = 129;
static const uint KH = 193;              // column-major, 65 rows
static const uint KY = 4353;
static const uint KINV = 4417;
static const uint KTARGET = 4418;

// Status bits. Fatal: 2 nonfinite residual or input, 16/32 local solve,
// 128 compression without free velocity support, 256 sweep budget exhausted.
// Informational: 1024 nonfinite Krylov value, 2048 nonfinite Anderson fit.
static const uint FATAL = 2 | 16 | 32 | 128 | 256;

static const uint HISTORY = 13;          // 12 differences

bool open(StructuredBuffer<uint> state, uint gate)
{
    return (state[FLAGS] & gate) == gate;
}

bool open(RWStructuredBuffer<uint> state, uint gate)
{
    return (state[FLAGS] & gate) == gate;
}

float3 load3(RWStructuredBuffer<float> v, uint i)
{
    return float3(v[3 * i], v[3 * i + 1], v[3 * i + 2]);
}

void store3(RWStructuredBuffer<float> v, uint i, float3 x)
{
    v[3 * i] = x.x;
    v[3 * i + 1] = x.y;
    v[3 * i + 2] = x.z;
}

float3 strain(uint a, StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
              RWStructuredBuffer<float> velocity)
{
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
    return HALF_ROOT2 * w;
}

float2 transpose_force(uint j, StructuredBuffer<uint> transpose_rows,
                       StructuredBuffer<uint> transpose_slots, StructuredBuffer<float> gradients,
                       RWStructuredBuffer<float> stress, uint offset)
{
    float2 f = float2(0.0, 0.0);
    for (uint i = 0; i < 9; i++)
    {
        uint a = transpose_rows[9 * j + i];
        if (a == INVALID) break;
        uint s = transpose_slots[9 * j + i];
        float gx = gradients[2 * (9 * a + s)];
        float gy = gradients[2 * (9 * a + s) + 1];
        float z0 = stress[offset + 3 * a];
        float z1 = stress[offset + 3 * a + 1];
        float z2 = stress[offset + 3 * a + 2];
        f += float2((z0 + z1) * gx + z2 * gy, (z0 - z1) * gy + z2 * gx);
    }
    return HALF_ROOT2 * f;
}

// Normal reaction nonnegative, tangential reaction inside the Coulomb disk.
float3 feasible(float3 v, float mu)
{
    v.x = max(v.x, 0.0);
    float size = length(v.yz);
    v.yz *= size > 0.0 ? min(1.0, mu * v.x / size) : 1.0;
    return v;
}

// Euclidean natural residual of one row, with the nonassociated radius mu z_n.
float natural(float3 z, float3 w, float mu)
{
    float3 argument = z - w;
    float radius = mu * max(z.x, 0.0);
    float size = length(argument.yz);
    float3 projected = float3(max(argument.x, 0.0),
                              argument.yz * (size > 0.0 ? min(1.0, radius / size) : 1.0));
    return length(z - projected);
}

void fail(RWStructuredBuffer<uint> state)
{
    state[RESULT] = 2;
    state[FLAGS] = 0;
}

void converge(RWStructuredBuffer<uint> state)
{
    state[RESULT] = 1;
    state[FLAGS] = 0;
}

// ---- Set-up ----------------------------------------------------------------

[shader("compute")]
[numthreads(256, 1, 1)]
void solver_rows_setup(uint3 tid : SV_DispatchThreadID, uniform uint row_count,
                       StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                       RWStructuredBuffer<float> trial, StructuredBuffer<float> gap,
                       RWStructuredBuffer<float> magnitude, RWStructuredBuffer<float> compression)
{
    // q = B u_trial + gap. Positive normal gaps do not enter the scaling.
    uint a = tid.x;
    if (a >= row_count) return;
    float3 q = strain(a, neighbors, gradients, trial)
               + float3(gap[3 * a], gap[3 * a + 1], gap[3 * a + 2]);
    magnitude[a] = max(max(-q.x, 0.0), max(abs(q.y), abs(q.z)));
    compression[a] = -q.x;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void solver_setup(StructuredBuffer<float> magnitude, StructuredBuffer<float> compression,
                  StructuredBuffer<float> norm, RWStructuredBuffer<uint> state,
                  RWStructuredBuffer<uint> status, RWStructuredBuffer<float> scalars)
{
    float m = magnitude[0], c = compression[0], n = norm[0];
    scalars[MAGNITUDE] = m;
    scalars[NORM] = n;
    scalars[SCALE] = 0.0;
    scalars[ERROR] = 0.0;
    uint ignored;
    if ((status[0] & 2) != 0 || !isfinite(m) || !isfinite(n))
    {
        InterlockedOr(status[0], 2, ignored);
        fail(state);
    }
    else if (!(n > 0.0) && c > 0.0)
    {
        InterlockedOr(status[0], 128, ignored);
        fail(state);
    }
    else if (m == 0.0 || !(n > 0.0)) converge(state);
    else
    {
        scalars[SCALE] = m / n;
        state[FLAGS] = RUNNING;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void stress_initial(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform float mu,
                    uniform uint warm, StructuredBuffer<uint> state, StructuredBuffer<float> scalars,
                    StructuredBuffer<float> initial, RWStructuredBuffer<float> z)
{
    // Begin with an admissible warm start in units of the solve scale.
    uint a = tid.x;
    if (a >= row_count) return;
    float3 value = float3(0.0, 0.0, 0.0);
    if (warm != 0 && open(state, RUNNING))
        value = feasible(float3(initial[3 * a], initial[3 * a + 1], initial[3 * a + 2])
                         / scalars[SCALE], mu);
    store3(z, a, value);
}

// ---- Exact velocity refresh and natural residual ---------------------------

[shader("compute")]
[numthreads(256, 1, 1)]
void velocity_refresh(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint node_count,
                      uniform float dt, StructuredBuffer<uint> state, StructuredBuffer<float> scalars,
                      StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
                      StructuredBuffer<float> gradients, StructuredBuffer<float> inverse_mass,
                      RWStructuredBuffer<float> trial, RWStructuredBuffer<float> z,
                      RWStructuredBuffer<float> u)
{
    // u = u_trial + dt scale M^-1 B.T z, recomputed rather than accumulated.
    uint j = tid.x;
    if (j >= node_count || !open(state, gate)) return;
    float2 f = transpose_force(j, transpose_rows, transpose_slots, gradients, z, 0);
    float factor = dt * scalars[SCALE];
    u[2 * j] = trial[2 * j] + factor * inverse_mass[2 * j] * f.x;
    u[2 * j + 1] = trial[2 * j + 1] + factor * inverse_mass[2 * j + 1] * f.y;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void residual_rows(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint row_count,
                   uniform float mu, StructuredBuffer<uint> state, StructuredBuffer<float> scalars,
                   StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                   StructuredBuffer<float> gap, RWStructuredBuffer<float> z,
                   RWStructuredBuffer<float> u, RWStructuredBuffer<float> residual)
{
    uint a = tid.x;
    if (a >= row_count || !open(state, gate)) return;
    float m = scalars[MAGNITUDE];
    float inverse = m > 0.0 ? 1.0 / m : 0.0;
    float3 w = (strain(a, neighbors, gradients, u) + float3(gap[3 * a], gap[3 * a + 1], gap[3 * a + 2]))
               * inverse;
    residual[a] = natural(load3(z, a), w, mu);
}

[shader("compute")]
[numthreads(1, 1, 1)]
void solver_start(uniform float tolerance, StructuredBuffer<float> residual,
                  RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status,
                  RWStructuredBuffer<float> scalars)
{
    if ((state[FLAGS] & RUNNING) == 0) return;
    float error = residual[0];
    scalars[ERROR] = error;
    if ((status[0] & FATAL) != 0) fail(state);
    else if (error <= tolerance) converge(state);
}

// ---- Colored Gauss-Seidel sweep --------------------------------------------

void begin_logic(uint budget, RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status)
{
    if ((state[FLAGS] & RUNNING) == 0) return;
    uint ignored;
    if ((status[0] & FATAL) != 0) fail(state);
    else if (state[SWEEPS] >= budget)
    {
        InterlockedOr(status[0], 256, ignored);
        fail(state);
    }
}

[shader("compute")]
[numthreads(1, 1, 1)]
void iteration_begin(uniform uint budget, RWStructuredBuffer<uint> state,
                     RWStructuredBuffer<uint> status)
{
    begin_logic(budget, state, status);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void history_begin(uint3 tid : SV_DispatchThreadID, uniform uint size,
                   StructuredBuffer<uint> state, RWStructuredBuffer<float> z,
                   RWStructuredBuffer<float> starts)
{
    uint e = tid.x;
    if (e >= size || !open(state, RUNNING)) return;
    starts[state[HEAD] * size + e] = z[e];
}

void sweep_row(uint a, float mu, float dt, RWStructuredBuffer<float> scalars,
               StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
               StructuredBuffer<float> gap, StructuredBuffer<float> diagonal,
               StructuredBuffer<float> rotation, StructuredBuffer<float> inverse_mass,
               RWStructuredBuffer<uint> status, RWStructuredBuffer<float> z, RWStructuredBuffer<float> u)
{
    // Rows of one color share no velocity node (static Q1 coloring), so the
    // incremental velocity update of each row has a single writer.
    float3 w = (strain(a, neighbors, gradients, u) + float3(gap[3 * a], gap[3 * a + 1], gap[3 * a + 2]))
               / scalars[MAGNITUDE];
    float c = rotation[2 * a], s = rotation[2 * a + 1];
    float3 previous = load3(z, a);
    float3 old = float3(previous.x, c * previous.y + s * previous.z, -s * previous.y + c * previous.z);
    float3 rotated = float3(w.x, c * w.y + s * w.z, -s * w.y + c * w.z);
    float3 d = float3(diagonal[3 * a], diagonal[3 * a + 1], diagonal[3 * a + 2]);
    uint flags = 0;
    float3 local = coulomb_local(d, rotated - d * old, mu, flags);
    if (flags != 0)
    {
        uint ignored;
        InterlockedOr(status[0], flags, ignored);
    }
    float3 next = float3(local.x, c * local.y - s * local.z, s * local.y + c * local.z);
    store3(z, a, next);
    float3 delta = next - previous;
    float factor = dt * scalars[SCALE];
    for (uint k = 0; k < 9; k++)
    {
        uint j = neighbors[9 * a + k];
        if (j == INVALID) continue;
        float3x2 B = block(gradients[2 * (9 * a + k)], gradients[2 * (9 * a + k) + 1]);
        float2 f = delta.x * B[0] + delta.y * B[1] + delta.z * B[2];
        u[2 * j] += factor * inverse_mass[2 * j] * f.x;
        u[2 * j + 1] += factor * inverse_mass[2 * j + 1] * f.y;
    }
}

[shader("compute")]
[numthreads(64, 1, 1)]
void sweep_color(uint3 tid : SV_DispatchThreadID, uniform uint color_offset,
                 uniform uint color_count, uniform float mu, uniform float dt,
                 StructuredBuffer<uint> state, RWStructuredBuffer<float> scalars,
                 StructuredBuffer<uint> color_rows, StructuredBuffer<uint> neighbors,
                 StructuredBuffer<float> gradients, StructuredBuffer<float> gap,
                 StructuredBuffer<float> diagonal, StructuredBuffer<float> rotation,
                 StructuredBuffer<float> inverse_mass, RWStructuredBuffer<uint> status,
                 RWStructuredBuffer<float> z, RWStructuredBuffer<float> u)
{
    uint t = tid.x;
    if (t >= color_count || !open(state, RUNNING)) return;
    sweep_row(color_rows[color_offset + t], mu, dt, scalars, neighbors, gradients, gap, diagonal,
              rotation, inverse_mass, status, z, u);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void history_end(uint3 tid : SV_DispatchThreadID, uniform uint size,
                 StructuredBuffer<uint> state, RWStructuredBuffer<float> z,
                 RWStructuredBuffer<float> starts, RWStructuredBuffer<float> steps)
{
    uint e = tid.x;
    if (e >= size || !open(state, RUNNING)) return;
    uint slot = state[HEAD] * size + e;
    steps[slot] = z[e] - starts[slot];
}

void sweep_end_logic(float tolerance, float error, RWStructuredBuffer<uint> state,
                     RWStructuredBuffer<uint> status, RWStructuredBuffer<float> scalars)
{
    if ((state[FLAGS] & RUNNING) == 0) return;
    state[SWEEPS] += 1;
    state[HEAD] = (state[HEAD] + 1) % HISTORY;
    state[COUNT] = min(state[COUNT] + 1, HISTORY);
    scalars[ERROR] = error;
    if ((status[0] & FATAL) != 0) fail(state);
    else if (error <= tolerance) converge(state);
    else if (state[COUNT] > 1) state[FLAGS] |= ANDERSON;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void sweep_end(uniform float tolerance, StructuredBuffer<float> residual,
               RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status,
               RWStructuredBuffer<float> scalars)
{
    sweep_end_logic(tolerance, residual[0], state, status, scalars);
}

// ---- Anderson acceleration: TSQR and one-sided Jacobi SVD -------------------

[shader("compute")]
[numthreads(256, 1, 1)]
void anderson_matrix(uint3 tid : SV_DispatchThreadID, uniform uint size,
                     StructuredBuffer<uint> state, RWStructuredBuffer<float> steps,
                     RWStructuredBuffer<float> matrix)
{
    // Row e of [df_0 .. df_11 | f_newest], oldest difference first; unused
    // columns are zero and fall below the rank cutoff.
    uint e = tid.x;
    if (e >= size || !open(state, RUNNING | ANDERSON)) return;
    uint count = state[COUNT];
    uint first = (state[HEAD] + HISTORY - count) % HISTORY;
    for (uint c = 0; c < HISTORY - 1; c++)
    {
        float value = 0.0;
        if (c + 1 < count)
        {
            uint i0 = (first + c) % HISTORY, i1 = (first + c + 1) % HISTORY;
            value = steps[i1 * size + e] - steps[i0 * size + e];
        }
        matrix[HISTORY * e + c] = value;
    }
    matrix[HISTORY * e + HISTORY - 1] = steps[((state[HEAD] + HISTORY - 1) % HISTORY) * size + e];
}

groupshared float tile[256][13];
groupshared float partial[256][13];
groupshared float squares[256];

[shader("compute")]
[numthreads(256, 1, 1)]
void tsqr_tile(uint3 lane : SV_GroupThreadID, uint3 group : SV_GroupID, uniform uint gate,
               uniform uint rows, StructuredBuffer<uint> state, StructuredBuffer<float> input,
               RWStructuredBuffer<float> output)
{
    // Householder R factor of one 256-row tile of a 13-column matrix. Stacked
    // tile factors are reduced by the same kernel until one tile remains.
    if (!open(state, gate)) return;  // uniform across the group
    uint i = lane.x;
    uint base = group.x * 256;
    uint m = min(256u, rows - base);
    for (uint c = 0; c < 13; c++) tile[i][c] = i < m ? input[13 * (base + i) + c] : 0.0;
    GroupMemoryBarrierWithGroupSync();
    for (uint c = 0; c < 13; c++)
    {
        float x = tile[i][c];
        squares[i] = i > c ? x * x : 0.0;
        GroupMemoryBarrierWithGroupSync();
        for (uint offset = 128; offset > 0; offset /= 2)
        {
            if (i < offset) squares[i] += squares[i + offset];
            GroupMemoryBarrierWithGroupSync();
        }
        float sigma = squares[0];
        float alpha = tile[c][c];
        GroupMemoryBarrierWithGroupSync();
        if (sigma == 0.0) continue;  // uniform: the column is already reduced
        float beta = alpha >= 0.0 ? -sqrt(alpha * alpha + sigma) : sqrt(alpha * alpha + sigma);
        float head = alpha - beta;
        float v = i == c ? head : (i > c ? x : 0.0);
        float vv = head * head + sigma;
        for (uint d = 0; d < 13; d++) partial[i][d] = d > c ? v * tile[i][d] : 0.0;
        GroupMemoryBarrierWithGroupSync();
        for (uint offset = 128; offset > 0; offset /= 2)
        {
            if (i < offset)
                for (uint d = c + 1; d < 13; d++) partial[i][d] += partial[i + offset][d];
            GroupMemoryBarrierWithGroupSync();
        }
        for (uint d = c + 1; d < 13; d++) tile[i][d] -= 2.0 * partial[0][d] / vv * v;
        if (i == c) tile[i][c] = beta;
        else if (i > c) tile[i][c] = 0.0;
        GroupMemoryBarrierWithGroupSync();
    }
    if (i < 13)
        for (uint c = 0; c < 13; c++)
            output[13 * (13 * group.x + i) + c] = c >= i ? tile[i][c] : 0.0;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void anderson_solve(uniform float rcond, StructuredBuffer<float> factor,
                    RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status,
                    RWStructuredBuffer<float> coefficients)
{
    // min ||df c - f|| from R = [R_df, Q.T f]: one-sided Jacobi SVD of R_df,
    // then the minimum-norm solution over singular values above the cutoff.
    if (!open(state, RUNNING | ANDERSON)) return;
    uint k = state[COUNT] - 1;
    float G[12][12];
    float V[12][12];
    float b[12];
    for (uint r = 0; r < 12; r++)
    {
        b[r] = factor[13 * r + 12];
        for (uint c = 0; c < 12; c++)
        {
            G[r][c] = factor[13 * r + c];
            V[r][c] = r == c ? 1.0 : 0.0;
        }
    }
    for (uint sweep = 0; sweep < 30; sweep++)
    {
        bool rotated = false;
        for (uint p = 0; p + 1 < k; p++)
            for (uint q = p + 1; q < k; q++)
            {
                float alpha = 0.0, beta = 0.0, gamma = 0.0;
                for (uint r = 0; r < 12; r++)
                {
                    alpha += G[r][p] * G[r][p];
                    beta += G[r][q] * G[r][q];
                    gamma += G[r][p] * G[r][q];
                }
                if (gamma == 0.0 || abs(gamma) <= 1e-7 * sqrt(alpha * beta)) continue;
                rotated = true;
                float zeta = (beta - alpha) / (2.0 * gamma);
                float t = (zeta >= 0.0 ? 1.0 : -1.0) / (abs(zeta) + sqrt(1.0 + zeta * zeta));
                float cs = rsqrt(1.0 + t * t), sn = cs * t;
                for (uint r = 0; r < 12; r++)
                {
                    float gp = G[r][p], gq = G[r][q];
                    G[r][p] = cs * gp - sn * gq;
                    G[r][q] = sn * gp + cs * gq;
                    float vp = V[r][p], vq = V[r][q];
                    V[r][p] = cs * vp - sn * vq;
                    V[r][q] = sn * vp + cs * vq;
                }
            }
        if (!rotated) break;
    }
    float sigma[12];
    float largest = 0.0;
    for (uint c = 0; c < 12; c++)
    {
        float s = 0.0;
        for (uint r = 0; r < 12; r++) s += G[r][c] * G[r][c];
        sigma[c] = c < k ? sqrt(s) : 0.0;
        largest = max(largest, sigma[c]);
    }
    float x[12];
    for (uint r = 0; r < 12; r++) x[r] = 0.0;
    for (uint c = 0; c < k; c++)
    {
        if (!(sigma[c] > rcond * largest)) continue;
        float projection = 0.0;
        for (uint r = 0; r < 12; r++) projection += G[r][c] * b[r];
        projection /= sigma[c] * sigma[c];
        for (uint r = 0; r < 12; r++) x[r] += projection * V[r][c];
    }
    bool finite = true;
    for (uint r = 0; r < 12; r++)
    {
        finite = finite && isfinite(x[r]);
        coefficients[r] = x[r];
    }
    if (!finite)
    {
        uint ignored;
        InterlockedOr(status[0], 2048, ignored);
        state[FLAGS] &= ~ANDERSON;
        state[ANDERSON_REJECTED] += 1;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void anderson_candidate(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform float mu,
                        StructuredBuffer<uint> state, StructuredBuffer<float> coefficients,
                        RWStructuredBuffer<float> starts, RWStructuredBuffer<float> steps,
                        RWStructuredBuffer<float> z, RWStructuredBuffer<float> candidate)
{
    // z - (dx + df) c, projected to feasibility.
    uint a = tid.x;
    if (a >= row_count || !open(state, RUNNING | ANDERSON)) return;
    uint size = 3 * row_count;
    uint count = state[COUNT];
    uint first = (state[HEAD] + HISTORY - count) % HISTORY;
    float3 value = load3(z, a);
    for (uint c = 0; c + 1 < count; c++)
    {
        uint i0 = ((first + c) % HISTORY) * size + 3 * a;
        uint i1 = ((first + c + 1) % HISTORY) * size + 3 * a;
        float3 dx = float3(starts[i1], starts[i1 + 1], starts[i1 + 2])
                    - float3(starts[i0], starts[i0 + 1], starts[i0 + 2]);
        float3 df = float3(steps[i1], steps[i1 + 1], steps[i1 + 2])
                    - float3(steps[i0], steps[i0 + 1], steps[i0 + 2]);
        value -= coefficients[c] * (dx + df);
    }
    store3(candidate, a, feasible(value, mu));
}

// ---- Candidate acceptance (Anderson and Newton) -----------------------------

void accept_logic(uint gate, uint newton, float error, bool finite, RWStructuredBuffer<uint> state,
                  RWStructuredBuffer<float> scalars)
{
    // Accept only a finite candidate that reduces the full natural residual.
    if ((state[FLAGS] & gate) != gate) return;
    if (finite && isfinite(error) && error < scalars[ERROR])
    {
        scalars[ERROR] = error;
        state[FLAGS] |= ACCEPTED;
        if (newton != 0)
        {
            state[NEWTON_STEPS] += 1;
            state[COUNT] = 0;
            state[FLAGS] &= ~SEARCH;
        }
        else state[ACCELERATED] += 1;
    }
    else if (newton == 0) state[ANDERSON_REJECTED] += 1;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void candidate_accept(uniform uint gate, uniform uint newton, StructuredBuffer<float> residual,
                      StructuredBuffer<uint> candidate_status, RWStructuredBuffer<uint> state,
                      RWStructuredBuffer<float> scalars)
{
    accept_logic(gate, newton, residual[0], (candidate_status[0] & 2) == 0, state, scalars);
}

void finish_logic(float tolerance, RWStructuredBuffer<uint> state, RWStructuredBuffer<float> scalars)
{
    state[FLAGS] &= ~(ACCEPTED | ANDERSON);
    if ((state[FLAGS] & RUNNING) != 0 && scalars[ERROR] <= tolerance) converge(state);
}

[shader("compute")]
[numthreads(1, 1, 1)]
void candidate_finish(uniform float tolerance, RWStructuredBuffer<uint> state,
                      RWStructuredBuffer<float> scalars)
{
    finish_logic(tolerance, state, scalars);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void combine(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint size,
             uniform float alpha, uniform float beta, StructuredBuffer<uint> state,
             RWStructuredBuffer<float> a, RWStructuredBuffer<float> b,
             RWStructuredBuffer<float> output)
{
    // output = alpha a + beta b; beta = 0 ignores b entirely.
    uint e = tid.x;
    if (e >= size || !open(state, gate)) return;
    output[e] = beta == 0.0 ? alpha * a[e] : alpha * a[e] + beta * b[e];
}

// ---- Newton: frozen generalized Jacobian and restarted GMRES ----------------

[shader("compute")]
[numthreads(1, 1, 1)]
void newton_request(uniform float tolerance, RWStructuredBuffer<uint> state,
                    StructuredBuffer<float> scalars)
{
    // Polish nearly converged but degenerate stick/slip transitions.
    float error = scalars[ERROR];
    if ((state[FLAGS] & RUNNING) != 0 && tolerance < error && error < 1e-3)
        state[FLAGS] |= REQUESTED;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void newton_begin(RWStructuredBuffer<uint> state)
{
    uint flags = state[FLAGS];
    if ((flags & (RUNNING | REQUESTED)) != (RUNNING | REQUESTED)) return;
    state[FLAGS] = (flags & ~(REQUESTED | KRYLOV_DONE)) | NEWTON | SEARCH;
    state[NEWTON_ATTEMPTS] += 1;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void newton_rows(uint3 tid : SV_DispatchThreadID, uniform uint row_count, uniform float mu,
                 StructuredBuffer<uint> state, StructuredBuffer<float> scalars,
                 StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                 StructuredBuffer<float> gap, RWStructuredBuffer<float> z,
                 RWStructuredBuffer<float> u, RWStructuredBuffer<float> frozen,
                 RWStructuredBuffer<float> rhs)
{
    // Frozen projection branch per row and right-hand side -natural residual.
    uint a = tid.x;
    if (a >= row_count || !open(state, RUNNING | NEWTON)) return;
    float3 w = (strain(a, neighbors, gradients, u) + float3(gap[3 * a], gap[3 * a + 1], gap[3 * a + 2]))
               / scalars[MAGNITUDE];
    float3 value = load3(z, a);
    float3 argument = value - w;
    float radius = mu * max(value.x, 0.0);
    float size = length(argument.yz);
    float2 unit = size > 0.0 ? argument.yz / size : float2(0.0, 0.0);
    float ratio = size > 0.0 ? radius / size : 0.0;
    float3 projected = float3(max(argument.x, 0.0),
                              argument.yz * (size > 0.0 ? min(1.0, radius / size) : 1.0));
    frozen[6 * a] = argument.x > 0.0 ? 1.0 : 0.0;
    frozen[6 * a + 1] = size > radius ? 1.0 : 0.0;
    frozen[6 * a + 2] = ratio;
    frozen[6 * a + 3] = unit.x;
    frozen[6 * a + 4] = unit.y;
    frozen[6 * a + 5] = value.x >= 0.0 ? 1.0 : 0.0;
    store3(rhs, a, projected - value);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void jv_force(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint node_count,
              uniform float dt, uniform uint offset, StructuredBuffer<uint> state,
              StructuredBuffer<float> scalars, StructuredBuffer<uint> transpose_rows,
              StructuredBuffer<uint> transpose_slots, StructuredBuffer<float> gradients,
              StructuredBuffer<float> inverse_mass, RWStructuredBuffer<float> source,
              RWStructuredBuffer<float> du)
{
    // du = dt / norm M^-1 B.T v.
    uint j = tid.x;
    if (j >= node_count || !open(state, gate)) return;
    float2 f = transpose_force(j, transpose_rows, transpose_slots, gradients, source, offset);
    float factor = dt / scalars[NORM];
    du[2 * j] = factor * inverse_mass[2 * j] * f.x;
    du[2 * j + 1] = factor * inverse_mass[2 * j + 1] * f.y;
}

float3 jv_row(uint a, float mu, StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
              RWStructuredBuffer<float> frozen, RWStructuredBuffer<float> source, uint offset,
              RWStructuredBuffer<float> du)
{
    // Generalized Jacobian of the natural residual at the frozen branch.
    float3 dw = strain(a, neighbors, gradients, du);
    float3 v = float3(source[offset + 3 * a], source[offset + 3 * a + 1], source[offset + 3 * a + 2]);
    float3 out = dw;
    out.x = frozen[6 * a] > 0.0 ? dw.x : v.x;
    if (frozen[6 * a + 1] > 0.0)
    {
        float2 unit = float2(frozen[6 * a + 3], frozen[6 * a + 4]);
        float2 da = v.yz - dw.yz;
        float2 tangent = frozen[6 * a + 2] * (da - unit * dot(unit, da))
                         + mu * frozen[6 * a + 5] * v.x * unit;
        out.yz = v.yz - tangent;
    }
    return out;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void jv_rows(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint row_count,
             uniform float mu, uniform uint offset, StructuredBuffer<uint> state,
             StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
             RWStructuredBuffer<float> frozen, RWStructuredBuffer<float> source,
             RWStructuredBuffer<float> du, RWStructuredBuffer<float> output)
{
    // Generalized Jacobian of the natural residual at the frozen branch.
    uint a = tid.x;
    if (a >= row_count || !open(state, gate)) return;
    float3 out = jv_row(a, mu, neighbors, gradients, frozen, source, offset, du);
    store3(output, a, out);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void dot_partial(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint size,
                 uniform uint chunks, uniform uint vectors, uniform uint total,
                 StructuredBuffer<uint> state, RWStructuredBuffer<float> basis,
                 RWStructuredBuffer<float> w, RWStructuredBuffer<float> partials)
{
    // Chunked inner products V_i . w for i < vectors; index `vectors` is w . w.
    uint t = tid.x;
    if (t >= total * chunks || !open(state, gate)) return;
    uint i = t / chunks, chunk = t % chunks;
    uint end = min(size, 256 * (chunk + 1));
    float sum = 0.0;
    for (uint e = 256 * chunk; e < end; e++)
        sum += (i < vectors ? basis[i * size + e] : w[e]) * w[e];
    partials[t] = sum;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void dot_finish(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint chunks,
                uniform uint total, StructuredBuffer<uint> state,
                RWStructuredBuffer<float> partials, RWStructuredBuffer<float> output)
{
    uint i = tid.x;
    if (i >= total || !open(state, gate)) return;
    float sum = 0.0;
    for (uint chunk = 0; chunk < chunks; chunk++) sum += partials[i * chunks + chunk];
    output[i] = sum;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void orthogonalise(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint size,
                   uniform uint vectors, StructuredBuffer<uint> state,
                   RWStructuredBuffer<float> basis, RWStructuredBuffer<float> h,
                   RWStructuredBuffer<float> w)
{
    // One classical Gram-Schmidt pass; called twice per Arnoldi step.
    uint e = tid.x;
    if (e >= size || !open(state, gate)) return;
    float value = w[e];
    for (uint i = 0; i < vectors; i++) value -= h[i] * basis[i * size + e];
    w[e] = value;
}

[shader("compute")]
[numthreads(256, 1, 1)]
void normalise(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint size,
               uniform uint offset, StructuredBuffer<uint> state, StructuredBuffer<float> krylov,
               RWStructuredBuffer<float> w, RWStructuredBuffer<float> basis)
{
    uint e = tid.x;
    if (e >= size || !open(state, gate)) return;
    basis[offset + e] = w[e] * krylov[KINV];
}

[shader("compute")]
[numthreads(1, 1, 1)]
void gmres_start(uniform uint cycle, uniform float rtol, StructuredBuffer<float> squares,
                 RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status,
                 RWStructuredBuffer<float> krylov)
{
    uint flags = state[FLAGS];
    if ((flags & (RUNNING | NEWTON)) != (RUNNING | NEWTON) || (flags & KRYLOV_DONE) != 0) return;
    float beta = sqrt(squares[0]);
    if (!isfinite(beta))
    {
        uint ignored;
        InterlockedOr(status[0], 1024, ignored);
        state[FLAGS] = flags | KRYLOV_DONE;
        return;
    }
    if (cycle == 0) krylov[KTARGET] = rtol * beta;
    if (!(beta > krylov[KTARGET]))
    {
        state[FLAGS] = flags | KRYLOV_DONE;
        return;
    }
    for (uint i = 0; i <= RESTART; i++) krylov[KG + i] = 0.0;
    krylov[KG] = beta;
    krylov[KINV] = 1.0 / beta;
    state[KRYLOV] = 0;
    state[FLAGS] = flags | GMRES | UPDATE;
}

void givens_column(uint j, float column[RESTART + 1], float norm2, RWStructuredBuffer<uint> state,
                   RWStructuredBuffer<uint> status, RWStructuredBuffer<float> krylov)
{
    // Reduce Hessenberg column j by the stored Givens rotations; then the new
    // rotation, the residual estimate and the termination decision.
    float total = 0.0;
    for (uint i = 0; i <= j; i++) total += column[i] * column[i];
    float next = sqrt(norm2);
    total = sqrt(total + next * next);
    for (uint i = 0; i < j; i++)
    {
        float c = krylov[KCS + i], s = krylov[KSN + i];
        float upper = c * column[i] + s * column[i + 1];
        column[i + 1] = -s * column[i] + c * column[i + 1];
        column[i] = upper;
    }
    float scale = max(abs(column[j]), next);
    float r = scale > 0.0 ? scale * length(float2(column[j], next) / scale) : 0.0;
    bool finite = isfinite(r) && isfinite(total);
    for (uint i = 0; i < j; i++) finite = finite && isfinite(column[i]);
    if (!finite || r == 0.0)
    {
        // Nonfinite or singular column: keep the previous steps only.
        uint ignored;
        if (!finite) InterlockedOr(status[0], 1024, ignored);
        else state[BREAKDOWNS] += 1;
        state[FLAGS] = (state[FLAGS] & ~GMRES) | KRYLOV_DONE;
        return;
    }
    float c = column[j] / r, s = next / r;
    krylov[KCS + j] = c;
    krylov[KSN + j] = s;
    column[j] = r;
    for (uint i = 0; i <= j; i++) krylov[KH + i + (RESTART + 1) * j] = column[i];
    float g = krylov[KG + j];
    krylov[KG + j] = c * g;
    krylov[KG + j + 1] = -s * g;
    krylov[KINV] = next > 0.0 ? 1.0 / next : 0.0;
    state[KRYLOV] = j + 1;
    state[GMRES_ITERATIONS] += 1;
    if (abs(krylov[KG + j + 1]) <= krylov[KTARGET])
    {
        state[FLAGS] = (state[FLAGS] & ~GMRES) | KRYLOV_DONE;
        state[GMRES_CONVERGED] += 1;
    }
    else if (next <= 1e-6 * total)
    {
        // Invariant Krylov subspace in FP32: restarting would repeat it.
        state[BREAKDOWNS] += 1;
        state[FLAGS] = (state[FLAGS] & ~GMRES) | KRYLOV_DONE;
    }
    else if (j + 1 == RESTART) state[FLAGS] &= ~GMRES;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void gmres_givens(uniform uint step, uniform uint groups, StructuredBuffer<float> first,
                  StructuredBuffer<float> second, StructuredBuffer<float> squares,
                  RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status,
                  RWStructuredBuffer<float> krylov)
{
    // Hessenberg column `step` from both Gram-Schmidt passes.
    if (!open(state, RUNNING | NEWTON | GMRES)) return;
    float column[RESTART + 1] = {};
    for (uint i = 0; i <= RESTART; i++) column[i] = i <= step ? first[i] + second[i] : 0.0;
    float norm2 = 0.0;
    for (uint g = 0; g < groups; g++) norm2 += squares[g];
    givens_column(step, column, norm2, state, status, krylov);
}

[shader("compute")]
[numthreads(1, 1, 1)]
void gmres_solve(RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status,
                 RWStructuredBuffer<float> krylov)
{
    // Back substitution with the reduced Hessenberg factor.
    if (!open(state, RUNNING | NEWTON | UPDATE)) return;
    int k = int(state[KRYLOV]);
    bool finite = true;
    for (int i = k - 1; i >= 0; i--)
    {
        float value = krylov[KG + i];
        for (int l = i + 1; l < k; l++) value -= krylov[KH + i + (RESTART + 1) * l] * krylov[KY + l];
        value /= krylov[KH + i + (RESTART + 1) * i];
        finite = finite && isfinite(value);
        krylov[KY + i] = value;
    }
    if (!finite)
    {
        uint ignored;
        InterlockedOr(status[0], 1024, ignored);
        state[KRYLOV] = 0;
        state[FLAGS] |= KRYLOV_DONE;
    }
}

[shader("compute")]
[numthreads(256, 1, 1)]
void basis_update(uint3 tid : SV_DispatchThreadID, uniform uint size,
                  StructuredBuffer<uint> state, StructuredBuffer<float> krylov,
                  RWStructuredBuffer<float> basis, RWStructuredBuffer<float> x)
{
    uint e = tid.x;
    if (e >= size || !open(state, RUNNING | NEWTON | UPDATE)) return;
    float value = x[e];
    for (uint i = 0; i < state[KRYLOV]; i++) value += krylov[KY + i] * basis[i * size + e];
    x[e] = value;
}

[shader("compute")]
[numthreads(1, 1, 1)]
void gmres_end(RWStructuredBuffer<uint> state)
{
    state[FLAGS] &= ~(GMRES | UPDATE);
}

[shader("compute")]
[numthreads(1, 1, 1)]
void newton_end(RWStructuredBuffer<uint> state)
{
    uint flags = state[FLAGS];
    if ((flags & NEWTON) == 0) return;
    if ((flags & SEARCH) != 0) state[NEWTON_REJECTED] += 1;
    state[FLAGS] = flags & ~(NEWTON | SEARCH | GMRES | UPDATE | KRYLOV_DONE);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void final_contact(uint3 tid : SV_DispatchThreadID, uniform uint size,
                   StructuredBuffer<float> scalars, RWStructuredBuffer<float> z,
                   RWStructuredBuffer<float> contact)
{
    uint e = tid.x;
    if (e >= size) return;
    contact[e] = z[e] * scalars[SCALE];
}

// ---- Arnoldi step ------------------------------------------------------------
// Inner products use one 64-thread workgroup per (vector, chunk); chunk g holds
// elements (g * per_thread + m) * 64 + t, m < per_thread, and per_thread is
// chosen so that there are at most 64 chunks. Updates reduce their own
// coefficients from those partials. All kernels stream the basis once.

groupshared float wave_totals[2];
groupshared float coefficient_shared[RESTART + 1];

[shader("compute")]
[numthreads(256, 1, 1)]
void jv_force_normalise(uint3 tid : SV_DispatchThreadID, uniform uint gate, uniform uint node_count,
                        uniform uint size, uniform uint slot, uniform float dt,
                        StructuredBuffer<uint> state, StructuredBuffer<float> scalars,
                        StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
                        StructuredBuffer<float> gradients, StructuredBuffer<float> inverse_mass,
                        StructuredBuffer<float> krylov, RWStructuredBuffer<float> source,
                        RWStructuredBuffer<float> basis, RWStructuredBuffer<float> du)
{
    // v_slot = source / h (stored), and du = dt / norm M^-1 B.T v_slot.
    uint t = tid.x;
    if (!open(state, gate)) return;
    float scale = krylov[KINV];
    if (t < size) basis[slot * size + t] = source[t] * scale;
    if (t >= node_count) return;
    float2 f = transpose_force(t, transpose_rows, transpose_slots, gradients, source, 0) * scale;
    float factor = dt / scalars[NORM];
    du[2 * t] = factor * inverse_mass[2 * t] * f.x;
    du[2 * t + 1] = factor * inverse_mass[2 * t + 1] * f.y;
}

[shader("compute")]
[numthreads(64, 1, 1)]
void arnoldi_dots(uint3 lane : SV_GroupThreadID, uint3 group : SV_GroupID, uniform uint gate,
                  uniform uint size, uniform uint per_thread, uniform uint groups, uniform uint squared,
                  StructuredBuffer<uint> state, RWStructuredBuffer<float> basis,
                  RWStructuredBuffer<float> w, RWStructuredBuffer<float> partials)
{
    // partials[i * groups + g] = chunk g of V_i . w, or of w . w when squared != 0.
    if (!open(state, gate)) return;
    uint t = lane.x, i = group.x / groups, g = group.x % groups;
    float sum = 0.0;
    for (uint m = 0; m < per_thread; m++)
    {
        uint e = (g * per_thread + m) * 64 + t;
        if (e >= size) break;
        float value = w[e];
        sum += (squared != 0 ? value : basis[i * size + e]) * value;
    }
    float total = WaveActiveSum(sum);
    if (WaveGetLaneIndex() == 0) wave_totals[t / WaveGetLaneCount()] = total;
    GroupMemoryBarrierWithGroupSync();
    if (t == 0) partials[i * groups + g] = wave_totals[0] + (WaveGetLaneCount() < 64 ? wave_totals[1] : 0.0);
}

[shader("compute")]
[numthreads(256, 1, 1)]
void arnoldi_update(uint3 tid : SV_DispatchThreadID, uint3 lane : SV_GroupThreadID, uniform uint gate,
                    uniform uint size, uniform uint vectors, uniform uint groups,
                    StructuredBuffer<uint> state, RWStructuredBuffer<float> basis,
                    StructuredBuffer<float> partials, RWStructuredBuffer<float> h,
                    RWStructuredBuffer<float> w)
{
    // One classical Gram-Schmidt pass, w -= sum_i h_i V_i, with h reduced from
    // the chunk partials in every workgroup.
    if (!open(state, gate)) return;
    uint t = lane.x, e = tid.x;
    for (uint i = t; i < vectors; i += 256)
    {
        float sum = 0.0;
        for (uint k = 0; k < groups; k++) sum += partials[i * groups + k];
        coefficient_shared[i] = sum;
        if (e == i) h[i] = sum;
    }
    GroupMemoryBarrierWithGroupSync();
    if (e >= size) return;
    float value = w[e];
    for (uint i = 0; i < vectors; i++) value -= coefficient_shared[i] * basis[i * size + e];
    w[e] = value;
}

// ---- One solver iteration in two single-workgroup kernels ---------------------
// For small systems (up to about 1,000 stress rows on the M4) an iteration is
// faster as two 1,024-thread workgroups than as about 33 dispatches. iteration_front runs the budget check, history, colored sweep,
// refresh, residual, sweep_end and the Anderson matrix; the TSQR and SVD
// dispatches follow; iteration_back evaluates and accepts the candidate. The
// control decisions use the same functions as the multi-dispatch path.

groupshared float iteration_wave[32];
groupshared uint iteration_flag;

float workgroup_max(uint t, float value)
{
    // Largest value over the 1,024 threads; leaves a barrier behind it.
    float wave = WaveActiveMax(value);
    if (WaveGetLaneIndex() == 0) iteration_wave[t / WaveGetLaneCount()] = wave;
    GroupMemoryBarrierWithGroupSync();
    float result = -3.402823466e38;
    for (uint k = 0; k < 1024 / WaveGetLaneCount(); k++) result = max(result, iteration_wave[k]);
    GroupMemoryBarrierWithGroupSync();
    return result;
}

void refresh_nodes(uint t, uint node_count, float dt, RWStructuredBuffer<float> scalars,
                   StructuredBuffer<uint> transpose_rows, StructuredBuffer<uint> transpose_slots,
                   StructuredBuffer<float> gradients, StructuredBuffer<float> inverse_mass,
                   RWStructuredBuffer<float> trial, RWStructuredBuffer<float> z, RWStructuredBuffer<float> u)
{
    float factor = dt * scalars[SCALE];
    for (uint j = t; j < node_count; j += 1024)
    {
        float2 f = transpose_force(j, transpose_rows, transpose_slots, gradients, z, 0);
        u[2 * j] = trial[2 * j] + factor * inverse_mass[2 * j] * f.x;
        u[2 * j + 1] = trial[2 * j + 1] + factor * inverse_mass[2 * j + 1] * f.y;
    }
}

float residual_max(uint t, uint row_count, float mu, RWStructuredBuffer<float> scalars,
                   StructuredBuffer<uint> neighbors, StructuredBuffer<float> gradients,
                   StructuredBuffer<float> gap, RWStructuredBuffer<float> z, RWStructuredBuffer<float> u,
                   inout bool finite)
{
    float m = scalars[MAGNITUDE];
    float inverse = m > 0.0 ? 1.0 / m : 0.0;
    float largest = -3.402823466e38;
    bool ok = true;
    for (uint a = t; a < row_count; a += 1024)
    {
        float3 w = (strain(a, neighbors, gradients, u) + float3(gap[3 * a], gap[3 * a + 1], gap[3 * a + 2]))
                   * inverse;
        float value = natural(load3(z, a), w, mu);
        if (isfinite(value)) largest = max(largest, value);
        else ok = false;
    }
    if (t == 0) iteration_flag = 0;
    GroupMemoryBarrierWithGroupSync();
    if (!ok) InterlockedOr(iteration_flag, 1);
    float result = workgroup_max(t, largest);
    finite = iteration_flag == 0;
    return result;
}

[shader("compute")]
[numthreads(1024, 1, 1)]
void iteration_front(uint3 lane : SV_GroupThreadID, uniform uint budget, uniform float tolerance,
                     uniform float mu, uniform float dt, uniform uint colors, uniform uint node_count,
                     uniform uint row_count, RWStructuredBuffer<uint> state, RWStructuredBuffer<uint> status,
                     RWStructuredBuffer<float> scalars, StructuredBuffer<uint> color_offsets,
                     StructuredBuffer<uint> color_rows, StructuredBuffer<uint> neighbors,
                     StructuredBuffer<float> gradients, StructuredBuffer<float> gap,
                     StructuredBuffer<float> diagonal, StructuredBuffer<float> rotation,
                     StructuredBuffer<float> inverse_mass, StructuredBuffer<uint> transpose_rows,
                     StructuredBuffer<uint> transpose_slots, RWStructuredBuffer<float> trial,
                     RWStructuredBuffer<float> z, RWStructuredBuffer<float> u,
                     RWStructuredBuffer<float> starts, RWStructuredBuffer<float> steps,
                     RWStructuredBuffer<float> matrix)
{
    uint t = lane.x, size = 3 * row_count;
    if (t == 0) begin_logic(budget, state, status);
    DeviceMemoryBarrierWithGroupSync();
    if (!open(state, RUNNING)) return;  // uniform: every thread reads after the barrier
    uint head = state[HEAD];
    for (uint e = t; e < size; e += 1024) starts[head * size + e] = z[e];
    DeviceMemoryBarrierWithGroupSync();
    for (uint c = 0; c < colors; c++)
    {
        for (uint k = color_offsets[c] + t; k < color_offsets[c + 1]; k += 1024)
            sweep_row(color_rows[k], mu, dt, scalars, neighbors, gradients, gap, diagonal, rotation,
                      inverse_mass, status, z, u);
        DeviceMemoryBarrierWithGroupSync();
    }
    // Periodic refresh avoids accumulating incremental roundoff.
    refresh_nodes(t, node_count, dt, scalars, transpose_rows, transpose_slots, gradients, inverse_mass,
                  trial, z, u);
    DeviceMemoryBarrierWithGroupSync();
    bool finite = true;
    float error = residual_max(t, row_count, mu, scalars, neighbors, gradients, gap, z, u, finite);
    for (uint e = t; e < size; e += 1024) steps[head * size + e] = z[e] - starts[head * size + e];
    DeviceMemoryBarrierWithGroupSync();
    if (t == 0)
    {
        uint ignored;
        if (!finite) InterlockedOr(status[0], 2, ignored);
        sweep_end_logic(tolerance, error, state, status, scalars);
    }
    DeviceMemoryBarrierWithGroupSync();
    if (!open(state, RUNNING | ANDERSON)) return;
    uint count = state[COUNT];
    uint first = (state[HEAD] + HISTORY - count) % HISTORY;
    uint newest = (state[HEAD] + HISTORY - 1) % HISTORY;
    for (uint e = t; e < size; e += 1024)
    {
        for (uint c = 0; c < HISTORY - 1; c++)
        {
            float value = 0.0;
            if (c + 1 < count)
            {
                uint i0 = (first + c) % HISTORY, i1 = (first + c + 1) % HISTORY;
                value = steps[i1 * size + e] - steps[i0 * size + e];
            }
            matrix[HISTORY * e + c] = value;
        }
        matrix[HISTORY * e + HISTORY - 1] = steps[newest * size + e];
    }
}

[shader("compute")]
[numthreads(1024, 1, 1)]
void iteration_back(uint3 lane : SV_GroupThreadID, uniform float tolerance, uniform float mu,
                    uniform float dt, uniform uint node_count, uniform uint row_count,
                    RWStructuredBuffer<uint> state, RWStructuredBuffer<float> scalars,
                    StructuredBuffer<float> coefficients, StructuredBuffer<uint> neighbors,
                    StructuredBuffer<float> gradients, StructuredBuffer<float> gap,
                    StructuredBuffer<float> inverse_mass, StructuredBuffer<uint> transpose_rows,
                    StructuredBuffer<uint> transpose_slots, RWStructuredBuffer<float> trial,
                    RWStructuredBuffer<float> starts, RWStructuredBuffer<float> steps,
                    RWStructuredBuffer<float> z, RWStructuredBuffer<float> u,
                    RWStructuredBuffer<float> candidate, RWStructuredBuffer<float> candidate_u)
{
    uint t = lane.x, size = 3 * row_count;
    if (!open(state, RUNNING | ANDERSON))
    {
        if (t == 0) finish_logic(tolerance, state, scalars);
        return;
    }
    uint count = state[COUNT];
    uint first = (state[HEAD] + HISTORY - count) % HISTORY;
    for (uint a = t; a < row_count; a += 1024)
    {
        // z - (dx + df) c, projected to feasibility.
        float3 value = load3(z, a);
        for (uint c = 0; c + 1 < count; c++)
        {
            uint i0 = ((first + c) % HISTORY) * size + 3 * a;
            uint i1 = ((first + c + 1) % HISTORY) * size + 3 * a;
            float3 dx = float3(starts[i1], starts[i1 + 1], starts[i1 + 2])
                        - float3(starts[i0], starts[i0 + 1], starts[i0 + 2]);
            float3 df = float3(steps[i1], steps[i1 + 1], steps[i1 + 2])
                        - float3(steps[i0], steps[i0 + 1], steps[i0 + 2]);
            value -= coefficients[c] * (dx + df);
        }
        store3(candidate, a, feasible(value, mu));
    }
    DeviceMemoryBarrierWithGroupSync();
    refresh_nodes(t, node_count, dt, scalars, transpose_rows, transpose_slots, gradients, inverse_mass,
                  trial, candidate, candidate_u);
    DeviceMemoryBarrierWithGroupSync();
    bool finite = true;
    float error = residual_max(t, row_count, mu, scalars, neighbors, gradients, gap, candidate, candidate_u,
                               finite);
    if (t == 0) accept_logic(RUNNING | ANDERSON, 0, error, finite, state, scalars);
    DeviceMemoryBarrierWithGroupSync();
    if ((state[FLAGS] & ACCEPTED) != 0)
    {
        for (uint e = t; e < size; e += 1024) z[e] = candidate[e];
        for (uint e = t; e < 2 * node_count; e += 1024) u[e] = candidate_u[e];
    }
    DeviceMemoryBarrierWithGroupSync();
    if (t == 0) finish_logic(tolerance, state, scalars);
}
