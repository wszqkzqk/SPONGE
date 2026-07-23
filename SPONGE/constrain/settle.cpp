#include "settle.h"

#include "velocity_projection.h"

static __global__ void remember_triangle_BA_CA(
    const int num_triangle_local, const CONSTRAIN_TRIANGLE* triangles,
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    VECTOR* last_triangle_BA, VECTOR* last_triangle_CA)
{
    CONSTRAIN_TRIANGLE triangle;
#ifdef USE_GPU
    int triangle_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (triangle_i < num_triangle_local)
#else
#pragma omp parallel for private(triangle)
    for (int triangle_i = 0; triangle_i < num_triangle_local; triangle_i++)
#endif
    {
        triangle = triangles[triangle_i];
        last_triangle_BA[triangle_i] = Get_Periodic_Displacement(
            crd[triangle.atom_B], crd[triangle.atom_A], cell, rcell);
        last_triangle_CA[triangle_i] = Get_Periodic_Displacement(
            crd[triangle.atom_C], crd[triangle.atom_A], cell, rcell);
    }
}

static __global__ void remember_pair_AB(const int num_task_local,
                                        const CONSTRAIN_PAIR* pairs,
                                        const VECTOR* crd, const LTMatrix3 cell,
                                        const LTMatrix3 rcell,
                                        VECTOR* last_pair_AB)
{
    CONSTRAIN_PAIR pair;
#ifdef USE_GPU
    int pair_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (pair_i < num_task_local)
#else
#pragma omp parallel for private(pair)
    for (int pair_i = 0; pair_i < num_task_local; pair_i++)
#endif
    {
        pair = pairs[pair_i];
        last_pair_AB[pair_i] = Get_Periodic_Displacement(
            crd[pair.atom_j_serial], crd[pair.atom_i_serial], cell, rcell);
    }
}

// 对几何信息进行转化
// 输入：rAB、rAC、rBC：三角形三边长
// 输入：mA、mB、mC：ABC三个的质量
// 输出：ra rb rc rd re：位置参数，当刚体三角形质心放置于原点时
// A点放置于(0,ra,0)，B点放置于(rc, rb, 0)，C点放置于(rd, re, 0)
static __device__ __host__ void Get_Rabcde_From_SSS(
    float rAB, float rAC, float rBC, float mA, float mB, float mC, float& ra,
    float& rb, float& rc, float& rd, float& re)
{
    float mTotal = mA + mB + mC;
    float Ax = 0;
    float Ay = 0;
    float Bx = -rAB;
    float By = 0;
    float costemp = (rBC * rBC - rAC * rAC - rAB * rAB) / (2 * rAC * rAB);
    float Cx = rAC * costemp;
    float sintemp = sqrtf(1.0f - costemp * costemp);
    float Cy = rAC * sintemp;

    float Ox = (Bx * mB + Cx * mC) / mTotal;
    float Oy = Cy * mC / mTotal;

    Ax -= Ox;
    Ay -= Oy;
    Bx -= Ox;
    By -= Oy;
    Cx -= Ox;
    Cy -= Oy;

    costemp = 1.0f / sqrtf(1.0f + Ax * Ax / Ay / Ay);
    sintemp = costemp * Ax / Ay;

    ra = Ax * sintemp + Ay * costemp;

    rc = Bx * costemp - By * sintemp;
    rb = Bx * sintemp + By * costemp;
    rd = Cx * costemp - Cy * sintemp;
    re = Cx * sintemp + Cy * costemp;

    if (ra < 0)
    {
        ra *= -1;
        rb *= -1;
        re *= -1;
    }
}

// 核心部分
//  部分参考了Shuichi & Peter: SETTLE: An Analytical Version of the SHAKE and
//  RATTLE Algorithm for Rigid Water Models A B C 三个点，O 质心
// 输入：rB0 上一步的B原子坐标（A为原点）；rC0 上一步的C原子坐标（A为原点）
// rA1 这一步的A原子坐标（质心为原点） rB1 这一步的B原子坐标（质心为原点）；rC1
// 这一步的C原子坐标（质心为原点） ra rb rc rd
// re：位置参数：当刚体三角形质心放置于原点，A点放置于(0,ra,0)，B点放置于(rc,
// rb, 0)，C点放置于(rd, re, 0)
//  mA、mB、mC：ABC三个的质量 dt:步长
//  half_exp_gamma_plus_half, exp_gamma: 同simple_constrain
// 输出：rA3 这一步限制后的A原子坐标（质心为原点） rB3
// 这一步限制后的B原子坐标（质心为 原点） rC3
// 这一步限制后的C原子坐标（质心为原点）vA vB vC 约束后的速度（原位替换） virial
// virial_vector 约束后的维里（原位替换）
static __device__ void SETTLE_DO_TRIANGLE(
    VECTOR rB0, VECTOR rC0, VECTOR rA1, VECTOR rB1, VECTOR rC1, float ra,
    float rb, float rc, float rd, float re, float mA, float mB, float mC,
    float dt, float half_exp_gamma_plus_half, float exp_gamma, VECTOR& rA3,
    VECTOR& rB3, VECTOR& rC3, VECTOR& vA, VECTOR& vB, VECTOR& vC,
    LTMatrix3& virial_tensor, int triangle_i)
{
    // 第0步：构建新坐标系
    // z轴垂直于上一步的BA和BC。 VECTOR ^ VECTOR 是外积
    VECTOR base_vector_z = rB0 ^ rC0;
    // x轴垂直于z轴和这一步的AO
    VECTOR base_vector_x = rA1 ^ base_vector_z;
    // y轴垂直于z轴和x轴
    VECTOR base_vector_y = base_vector_z ^ base_vector_x;
    // 归一化
    base_vector_x =
        rnorm3df(base_vector_x.x, base_vector_x.y, base_vector_x.z) *
        base_vector_x;
    base_vector_y =
        rnorm3df(base_vector_y.x, base_vector_y.y, base_vector_y.z) *
        base_vector_y;
    base_vector_z =
        rnorm3df(base_vector_z.x, base_vector_z.y, base_vector_z.z) *
        base_vector_z;

    // 第1步：投影至新坐标系
    //      rA0d = {0, 0, 0};
    VECTOR rB0d = {base_vector_x * rB0, base_vector_y * rB0, 0};
    VECTOR rC0d = {base_vector_x * rC0, base_vector_y * rC0, 0};
    VECTOR rA1d = {0, 0, base_vector_z * rA1};
    VECTOR rB1d = {base_vector_x * rB1, base_vector_y * rB1,
                   base_vector_z * rB1};
    VECTOR rC1d = {base_vector_x * rC1, base_vector_y * rC1,
                   base_vector_z * rC1};

    // 第2步：绕base_vector_y旋转psi，绕base_vector_x旋转phi得到rX2d
    float sinphi = rA1d.z / ra;
    float cosphi = sqrtf(1.0f - sinphi * sinphi);
    float sinpsi =
        (rB1d.z - rC1d.z - (rb - re) * sinphi) / ((rd - rc) * cosphi);
    float cospsi = sqrtf(1.0f - sinpsi * sinpsi);

    VECTOR rA2d = {0.0f, ra * cosphi, rA1d.z};
    VECTOR rB2d = {rc * cospsi, rb * cosphi + rc * sinpsi * sinphi, rB1d.z};
    VECTOR rC2d = {rd * cospsi, re * cosphi + rd * sinpsi * sinphi, rC1d.z};

    // 第3步：计算辅助变量 alpha、beta、gamma
    float alpha =
        rB2d.x * rB0d.x + rC2d.x * rC0d.x + rB2d.y * rB0d.y + rC2d.y * rC0d.y;
    float beta =
        -rB2d.x * rB0d.y - rC2d.x * rC0d.y + rB2d.y * rB0d.x + rC2d.y * rC0d.x;
    float gamma =
        rB1d.y * rB0d.x - rB1d.x * rB0d.y + rC1d.y * rC0d.x - rC1d.x * rC0d.y;

    // 第4步：绕base_vector_z旋转theta
    float temp = alpha * alpha + beta * beta;
    float sintheta =
        (alpha * gamma - beta * sqrtf(temp - gamma * gamma)) / temp;
    float costheta = sqrt(1.0f - sintheta * sintheta);
    VECTOR rA3d = {-rA2d.y * sintheta, rA2d.y * costheta, rA2d.z};
    VECTOR rB3d = {rB2d.x * costheta - rB2d.y * sintheta,
                   rB2d.x * sintheta + rB2d.y * costheta, rB2d.z};
    VECTOR rC3d = {rC2d.x * costheta - rC2d.y * sintheta,
                   rC2d.x * sintheta + rC2d.y * costheta, rC2d.z};

    // 第5步：投影回去
    rA3 = {rA3d.x * base_vector_x.x + rA3d.y * base_vector_y.x +
               rA3d.z * base_vector_z.x,
           rA3d.x * base_vector_x.y + rA3d.y * base_vector_y.y +
               rA3d.z * base_vector_z.y,
           rA3d.x * base_vector_x.z + rA3d.y * base_vector_y.z +
               rA3d.z * base_vector_z.z};

    rB3 = {rB3d.x * base_vector_x.x + rB3d.y * base_vector_y.x +
               rB3d.z * base_vector_z.x,
           rB3d.x * base_vector_x.y + rB3d.y * base_vector_y.y +
               rB3d.z * base_vector_z.y,
           rB3d.x * base_vector_x.z + rB3d.y * base_vector_y.z +
               rB3d.z * base_vector_z.z};

    rC3 = {rC3d.x * base_vector_x.x + rC3d.y * base_vector_y.x +
               rC3d.z * base_vector_z.x,
           rC3d.x * base_vector_x.y + rC3d.y * base_vector_y.y +
               rC3d.z * base_vector_z.y,
           rC3d.x * base_vector_x.z + rC3d.y * base_vector_y.z +
               rC3d.z * base_vector_z.z};

    // 第6步：计算约束造成的速度变化和维里变化
    // 节约寄存器，把不用的rX1d拿来当delta vX用
    temp = exp_gamma / dt / half_exp_gamma_plus_half;
    rA1d = temp * (rA3 - rA1);
    rB1d = temp * (rB3 - rB1);
    rC1d = temp * (rC3 - rC1);

    vA = vA + rA1d;
    vB = vB + rB1d;
    vC = vC + rC1d;
    // 节约寄存器，把不用的rX0d拿来当FX用
    temp = 1.0f / dt / dt / half_exp_gamma_plus_half;
    // rA0d = temp * mA * (rA3 - rA1);
    rB0d = temp * mB * (rB3 - rB1);
    rC0d = temp * mC * (rC3 - rC1);

    virial_tensor = Get_Virial_From_Force_Dis(rB0d, rB0) +
                    Get_Virial_From_Force_Dis(rC0d, rC0);
}

static __global__ void settle_triangle(
    int num_task_local, CONSTRAIN_TRIANGLE* triangles, const float* d_mass,
    VECTOR* crd, LTMatrix3 cell, LTMatrix3 rcell, VECTOR* last_triangle_BA,
    VECTOR* last_triangle_CA, float dt, float exp_gamma,
    float half_exp_gamma_plus_half, VECTOR* vel, LTMatrix3* virial_tensor)
{
    CONSTRAIN_TRIANGLE triangle;
    VECTOR rO;
    VECTOR rA, rB, rC;
    float mA, mB, mC;
#ifdef USE_GPU
    int triangle_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (triangle_i < num_task_local)
#else
#pragma omp parallel for private(triangle, rO, rA, rB, rC, mA, mB, mC)
    for (int triangle_i = 0; triangle_i < num_task_local; triangle_i++)
#endif
    {
        triangle = triangles[triangle_i];
        rA = crd[triangle.atom_A];
        rB = Get_Periodic_Displacement(crd[triangle.atom_B], rA, cell, rcell);
        rC = Get_Periodic_Displacement(crd[triangle.atom_C], rA, cell, rcell);
        mA = d_mass[triangle.atom_A];
        mB = d_mass[triangle.atom_B];
        mC = d_mass[triangle.atom_C];

        rO = 1.0f / (mA + mB + mC) * (mB * rB + mC * rC) + rA;
        rA = rA - rO;
        rB = rB + rA;
        rC = rC + rA;

        SETTLE_DO_TRIANGLE(
            last_triangle_BA[triangle_i], last_triangle_CA[triangle_i], rA, rB,
            rC, triangle.ra, triangle.rb, triangle.rc, triangle.rd, triangle.re,
            mA, mB, mC, dt, half_exp_gamma_plus_half, exp_gamma, rA, rB, rC,
            vel[triangle.atom_A], vel[triangle.atom_B], vel[triangle.atom_C],
            virial_tensor[triangle_i], triangle_i);

        crd[triangle.atom_A] = rA + rO;
        crd[triangle.atom_B] = rB + rO;
        crd[triangle.atom_C] = rC + rO;
    }
}

static __global__ void settle_pair(int num_task_local, CONSTRAIN_PAIR* pairs,
                                   const float* d_mass, VECTOR* crd,
                                   LTMatrix3 cell, LTMatrix3 rcell,
                                   VECTOR* last_pair_AB, float dt,
                                   float exp_gamma,
                                   float half_exp_gamma_plus_half, VECTOR* vel,
                                   LTMatrix3* virial_tensor)
{
    CONSTRAIN_PAIR pair;
    VECTOR r1, r2, kr2;
    float mA, mB, r0r0, r1r1, r1r2, r2r2, k;
#ifdef USE_GPU
    int pair_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (pair_i < num_task_local)
#else
#pragma omp parallel for private(pair, r1, r2, kr2, mA, mB, r0r0, r1r1, r1r2, \
                                     r2r2, k)
    for (int pair_i = 0; pair_i < num_task_local; pair_i++)
#endif
    {
        pair = pairs[pair_i];

        r1 = Get_Periodic_Displacement(crd[pair.atom_j_serial],
                                       crd[pair.atom_i_serial], cell, rcell);
        r2 = last_pair_AB[pair_i];
        mA = d_mass[pair.atom_i_serial];
        mB = d_mass[pair.atom_j_serial];

        r0r0 = pair.constant_r * pair.constant_r;
        r1r1 = r1 * r1;
        r1r2 = r1 * r2;
        r2r2 = r2 * r2;

        k = (sqrt(r1r2 * r1r2 - r1r1 * r2r2 + r2r2 * r0r0) - r1r2) / r2r2;
        kr2 = k * r2;

        r1 = -mB * pair.constrain_k * kr2;
        kr2 = mA * pair.constrain_k * kr2;

        crd[pair.atom_i_serial] = crd[pair.atom_i_serial] + r1;
        crd[pair.atom_j_serial] = crd[pair.atom_j_serial] + kr2;

        k = exp_gamma / dt / half_exp_gamma_plus_half;
        vel[pair.atom_i_serial] = vel[pair.atom_i_serial] + k * r1;
        vel[pair.atom_j_serial] = vel[pair.atom_j_serial] + k * kr2;

        r1 = k * mB / dt / exp_gamma * kr2;
        virial_tensor[pair_i] = Get_Virial_From_Force_Dis(r1, r2);
    }
}

struct SETTLE_PAIR_NUMERIC_DIAGNOSTIC
{
    int index_valid;
    int ownership_valid;
    int atom_i;
    int atom_j;
    int global_i;
    int global_j;
    float constant_r;
    float constrain_k;
    float mass_i;
    float mass_j;
    VECTOR crd_i;
    VECTOR crd_j;
    VECTOR vel_i;
    VECTOR vel_j;
    VECTOR frc_i;
    VECTOR frc_j;
    VECTOR r1;
    VECTOR r2;
    float r1r1;
    float r1r2;
    float r2r2;
    float radicand;
    float root;
    float k;
};

static __global__ void diagnose_settle_pairs(
    const int pair_numbers, const int owned_atom_numbers,
    const int active_atom_numbers, const CONSTRAIN_PAIR* pairs,
    const int* atom_local, const VECTOR* frc, const float* mass,
    const VECTOR* crd, const VECTOR* vel, const LTMatrix3 cell,
    const LTMatrix3 rcell, const VECTOR* last_pair_AB,
    SETTLE_PAIR_NUMERIC_DIAGNOSTIC* diagnostics)
{
#ifdef USE_GPU
    int pair_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (pair_i < pair_numbers)
#else
#pragma omp parallel for
    for (int pair_i = 0; pair_i < pair_numbers; pair_i++)
#endif
    {
        SETTLE_PAIR_NUMERIC_DIAGNOSTIC result = {};
        const CONSTRAIN_PAIR pair = pairs[pair_i];
        result.atom_i = pair.atom_i_serial;
        result.atom_j = pair.atom_j_serial;
        result.constant_r = pair.constant_r;
        result.constrain_k = pair.constrain_k;
        result.index_valid =
            result.atom_i >= 0 && result.atom_i < active_atom_numbers &&
            result.atom_j >= 0 && result.atom_j < active_atom_numbers;
        if (result.index_valid)
        {
            result.global_i = atom_local[result.atom_i];
            result.global_j = atom_local[result.atom_j];
            result.ownership_valid = result.atom_i < owned_atom_numbers &&
                                     result.atom_j < owned_atom_numbers;
        }
        if (result.ownership_valid)
        {
            result.mass_i = mass[result.atom_i];
            result.mass_j = mass[result.atom_j];
            result.crd_i = crd[result.atom_i];
            result.crd_j = crd[result.atom_j];
            result.vel_i = vel[result.atom_i];
            result.vel_j = vel[result.atom_j];
            result.frc_i = frc[result.atom_i];
            result.frc_j = frc[result.atom_j];
            result.r1 = Get_Periodic_Displacement(result.crd_j, result.crd_i,
                                                  cell, rcell);
            result.r2 = last_pair_AB[pair_i];
            const float r0r0 = pair.constant_r * pair.constant_r;
            result.r1r1 = result.r1 * result.r1;
            result.r1r2 = result.r1 * result.r2;
            result.r2r2 = result.r2 * result.r2;
            result.radicand = result.r1r2 * result.r1r2 -
                              result.r1r1 * result.r2r2 + result.r2r2 * r0r0;
            result.root = sqrt(result.radicand);
            result.k = (result.root - result.r1r2) / result.r2r2;
        }
        diagnostics[pair_i] = result;
    }
}

static bool Settle_Diagnostic_Float_Is_Finite(const float& value)
{
    return Float_Memory_Is_Finite(&value);
}

static bool Settle_Diagnostic_Vector_Is_Finite(const VECTOR& value)
{
    return Settle_Diagnostic_Float_Is_Finite(value.x) &&
           Settle_Diagnostic_Float_Is_Finite(value.y) &&
           Settle_Diagnostic_Float_Is_Finite(value.z);
}

static bool Settle_Pair_Diagnostic_Is_Valid(
    const SETTLE_PAIR_NUMERIC_DIAGNOSTIC& value)
{
    return value.index_valid && value.ownership_valid &&
           Settle_Diagnostic_Float_Is_Finite(value.constant_r) &&
           value.constant_r > 0.0f &&
           Settle_Diagnostic_Float_Is_Finite(value.constrain_k) &&
           Settle_Diagnostic_Float_Is_Finite(value.mass_i) &&
           Settle_Diagnostic_Float_Is_Finite(value.mass_j) &&
           value.mass_i + value.mass_j > 0.0f &&
           Settle_Diagnostic_Vector_Is_Finite(value.crd_i) &&
           Settle_Diagnostic_Vector_Is_Finite(value.crd_j) &&
           Settle_Diagnostic_Vector_Is_Finite(value.vel_i) &&
           Settle_Diagnostic_Vector_Is_Finite(value.vel_j) &&
           Settle_Diagnostic_Vector_Is_Finite(value.frc_i) &&
           Settle_Diagnostic_Vector_Is_Finite(value.frc_j) &&
           Settle_Diagnostic_Vector_Is_Finite(value.r1) &&
           Settle_Diagnostic_Vector_Is_Finite(value.r2) &&
           Settle_Diagnostic_Float_Is_Finite(value.r1r1) &&
           Settle_Diagnostic_Float_Is_Finite(value.r1r2) &&
           Settle_Diagnostic_Float_Is_Finite(value.r2r2) && value.r2r2 > 0.0f &&
           Settle_Diagnostic_Float_Is_Finite(value.radicand) &&
           value.radicand >= 0.0f &&
           Settle_Diagnostic_Float_Is_Finite(value.root) &&
           Settle_Diagnostic_Float_Is_Finite(value.k);
}

static unsigned int Settle_Diagnostic_Float_Bits(const float value)
{
    static_assert(sizeof(unsigned int) == sizeof(float),
                  "SPONGE requires 32-bit IEEE-754 floats");
    unsigned int bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static bool Synchronize_Settle_Diagnostic_Kernel(CONTROLLER* controller,
                                                 const char* function,
                                                 const char* kernel_name)
{
#ifdef GPU_ARCH_NAME
    const deviceError_t launch_error = deviceGetLastError();
    if (launch_error != 0)
    {
        const std::string reason = std::string("Reason:\n\tthe ") +
                                   kernel_name + " kernel failed to launch\n";
        controller->Throw_Device_Error(launch_error, function, reason.c_str());
        return false;
    }
    const deviceError_t synchronization_error = hostDeviceSynchronize();
    if (synchronization_error != 0)
    {
        const std::string reason = std::string("Reason:\n\tthe ") +
                                   kernel_name +
                                   " kernel failed during execution\n";
        controller->Throw_Device_Error(synchronization_error, function,
                                       reason.c_str());
        return false;
    }
#endif
    return true;
}

static void Preflight_Settle_Pairs_Once(
    CONTROLLER* controller, const int pair_numbers,
    const int owned_atom_numbers, const int active_atom_numbers,
    const CONSTRAIN_PAIR* pairs, const int* atom_local, const VECTOR* frc,
    const float* mass, const VECTOR* crd, const VECTOR* vel,
    const LTMatrix3 cell, const LTMatrix3 rcell, const VECTOR* last_pair_AB,
    const float dt, const float v_factor, const float x_factor,
    std::vector<SETTLE_PAIR_NUMERIC_DIAGNOSTIC>* preflight)
{
    preflight->clear();
    if (pair_numbers <= 0)
    {
        controller->printf(
            "SETTLE_PAIR_DIAG stage=preflight status=no_pairs\n");
        return;
    }

    SETTLE_PAIR_NUMERIC_DIAGNOSTIC* d_diagnostics = NULL;
    Device_Malloc_Safely((void**)&d_diagnostics,
                         sizeof(SETTLE_PAIR_NUMERIC_DIAGNOSTIC) * pair_numbers);
    constexpr int diagnostic_threads = 128;
    Launch_Device_Kernel(
        diagnose_settle_pairs,
        (pair_numbers + diagnostic_threads - 1) / diagnostic_threads,
        diagnostic_threads, 0, NULL, pair_numbers, owned_atom_numbers,
        active_atom_numbers, pairs, atom_local, frc, mass, crd, vel, cell,
        rcell, last_pair_AB, d_diagnostics);
    if (!Synchronize_Settle_Diagnostic_Kernel(
            controller, "Preflight_Settle_Pairs_Once", "SETTLE pair preflight"))
    {
        deviceFree(d_diagnostics);
        return;
    }
    std::vector<SETTLE_PAIR_NUMERIC_DIAGNOSTIC> diagnostics(pair_numbers);
    deviceMemcpy(diagnostics.data(), d_diagnostics,
                 sizeof(SETTLE_PAIR_NUMERIC_DIAGNOSTIC) * pair_numbers,
                 deviceMemcpyDeviceToHost);
    deviceFree(d_diagnostics);

    int first_invalid = -1;
    int invalid_count = 0;
    int minimum_index = -1;
    float minimum_radicand = FLT_MAX;
    for (int pair_i = 0; pair_i < pair_numbers; pair_i++)
    {
        const SETTLE_PAIR_NUMERIC_DIAGNOSTIC& value = diagnostics[pair_i];
        if (!Settle_Pair_Diagnostic_Is_Valid(value))
        {
            if (first_invalid < 0) first_invalid = pair_i;
            invalid_count++;
        }
        if (Settle_Diagnostic_Float_Is_Finite(value.radicand) &&
            value.radicand < minimum_radicand)
        {
            minimum_radicand = value.radicand;
            minimum_index = pair_i;
        }
    }

    if (first_invalid < 0)
    {
        const SETTLE_PAIR_NUMERIC_DIAGNOSTIC& minimum =
            diagnostics[minimum_index];
        controller->printf(
            "SETTLE_PAIR_DIAG stage=preflight status=finite pairs=%d "
            "min_pair=%d "
            "local/global=(%d,%d)/(%d,%d) r0=%.9g radicand=%.9g\n",
            pair_numbers, minimum_index, minimum.atom_i, minimum.atom_j,
            minimum.global_i, minimum.global_j, minimum.constant_r,
            minimum.radicand);
        preflight->swap(diagnostics);
        return;
    }

    const SETTLE_PAIR_NUMERIC_DIAGNOSTIC& bad = diagnostics[first_invalid];
    const char* failure = !bad.index_valid ? "an out-of-range local mapping"
                          : !bad.ownership_valid
                              ? "a cross-domain ghost mapping"
                              : "an invalid analytic update";
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown, "Preflight_Settle_Pairs_Once",
        "Reason:\n\tthe first SETTLE pair preflight found %s "
        "(pair=%d, invalid_pairs=%d/%d, index/ownership_valid=(%d,%d), "
        "owned/active_atoms=(%d,%d), "
        "local=(%d,%d), global=(%d,%d), r0=%.9g, constrain_k=%.9g, "
        "mass=(%.9g,%.9g), dt/v_factor/x_factor=(%.9g,%.9g,%.9g), "
        "r1=(%.9g,%.9g,%.9g), r2=(%.9g,%.9g,%.9g), "
        "r1r1/r1r2/r2r2=(%.9g,%.9g,%.9g), "
        "radicand/root/k=(%.9g,%.9g,%.9g), radicand_bits=0x%08x, "
        "crd_i=(%.9g,%.9g,%.9g), crd_j=(%.9g,%.9g,%.9g), "
        "vel_i=(%.9g,%.9g,%.9g), vel_j=(%.9g,%.9g,%.9g), "
        "frc_i=(%.9g,%.9g,%.9g), frc_j=(%.9g,%.9g,%.9g))\n",
        failure, first_invalid, invalid_count, pair_numbers, bad.index_valid,
        bad.ownership_valid, owned_atom_numbers, active_atom_numbers,
        bad.atom_i, bad.atom_j, bad.global_i, bad.global_j, bad.constant_r,
        bad.constrain_k, bad.mass_i, bad.mass_j, dt, v_factor, x_factor,
        bad.r1.x, bad.r1.y, bad.r1.z, bad.r2.x, bad.r2.y, bad.r2.z, bad.r1r1,
        bad.r1r2, bad.r2r2, bad.radicand, bad.root, bad.k,
        Settle_Diagnostic_Float_Bits(bad.radicand), bad.crd_i.x, bad.crd_i.y,
        bad.crd_i.z, bad.crd_j.x, bad.crd_j.y, bad.crd_j.z, bad.vel_i.x,
        bad.vel_i.y, bad.vel_i.z, bad.vel_j.x, bad.vel_j.y, bad.vel_j.z,
        bad.frc_i.x, bad.frc_i.y, bad.frc_i.z, bad.frc_j.x, bad.frc_j.y,
        bad.frc_j.z);
}

static void Diagnose_Actual_Settle_Pair_Output_Once(
    CONTROLLER* controller, const int owned_atom_numbers, const int* atom_local,
    const VECTOR* crd, const VECTOR* vel, const float dt, const float v_factor,
    const float x_factor,
    const std::vector<SETTLE_PAIR_NUMERIC_DIAGNOSTIC>& preflight)
{
    if (!Synchronize_Settle_Diagnostic_Kernel(
            controller, "Diagnose_Actual_Settle_Pair_Output_Once",
            "actual SETTLE pair"))
        return;

    if (owned_atom_numbers <= 0)
    {
        controller->printf(
            "SETTLE_PAIR_DIAG stage=post_pair status=no_owned_atoms\n");
        return;
    }

    const std::size_t atom_count = static_cast<std::size_t>(owned_atom_numbers);
    std::vector<VECTOR> coordinates(atom_count);
    std::vector<VECTOR> velocities(atom_count);
    std::vector<int> global_atoms(atom_count);
    deviceMemcpy(coordinates.data(), crd, atom_count * sizeof(VECTOR),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(velocities.data(), vel, atom_count * sizeof(VECTOR),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(global_atoms.data(), atom_local, atom_count * sizeof(int),
                 deviceMemcpyDeviceToHost);

    int first_invalid = -1;
    int invalid_atoms = 0;
    for (int atom_i = 0; atom_i < owned_atom_numbers; atom_i++)
    {
        if (!Settle_Diagnostic_Vector_Is_Finite(coordinates[atom_i]) ||
            !Settle_Diagnostic_Vector_Is_Finite(velocities[atom_i]))
        {
            if (first_invalid < 0) first_invalid = atom_i;
            invalid_atoms++;
        }
    }
    if (first_invalid < 0)
    {
        controller->printf(
            "SETTLE_PAIR_DIAG stage=post_pair status=finite atoms=%d "
            "pairs=%zu\n",
            owned_atom_numbers, preflight.size());
        return;
    }

    int writer_pair = -1;
    int writer_matches = 0;
    for (std::size_t pair_i = 0; pair_i < preflight.size(); pair_i++)
    {
        const SETTLE_PAIR_NUMERIC_DIAGNOSTIC& value = preflight[pair_i];
        if (value.atom_i == first_invalid || value.atom_j == first_invalid)
        {
            if (writer_pair < 0) writer_pair = static_cast<int>(pair_i);
            writer_matches++;
        }
    }

    SETTLE_PAIR_NUMERIC_DIAGNOSTIC writer = {};
    VECTOR post_crd_i = {};
    VECTOR post_crd_j = {};
    VECTOR post_vel_i = {};
    VECTOR post_vel_j = {};
    if (writer_pair >= 0)
    {
        writer = preflight[writer_pair];
        post_crd_i = coordinates[writer.atom_i];
        post_crd_j = coordinates[writer.atom_j];
        post_vel_i = velocities[writer.atom_i];
        post_vel_j = velocities[writer.atom_j];
    }

    const VECTOR& bad_crd = coordinates[first_invalid];
    const VECTOR& bad_vel = velocities[first_invalid];
    const char* bad_field = !Settle_Diagnostic_Vector_Is_Finite(bad_crd)
                                ? "coordinate"
                                : "velocity";
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown,
        "Diagnose_Actual_Settle_Pair_Output_Once",
        "Reason:\n\tthe actual SETTLE pair kernel first produced a "
        "non-finite %s (invalid_atoms=%d/%d, local/global atom=%d/%d, "
        "writer_pair/matches=%d/%d, writer local/global=(%d,%d)/(%d,%d), "
        "r0=%.9g, constrain_k=%.9g, mass=(%.9g,%.9g), "
        "dt/v_factor/x_factor=(%.9g,%.9g,%.9g), "
        "pre_crd_i=(%.9g,%.9g,%.9g), pre_crd_j=(%.9g,%.9g,%.9g), "
        "pre_vel_i=(%.9g,%.9g,%.9g), pre_vel_j=(%.9g,%.9g,%.9g), "
        "pre_frc_i=(%.9g,%.9g,%.9g), pre_frc_j=(%.9g,%.9g,%.9g), "
        "pre_r1=(%.9g,%.9g,%.9g), pre_r2=(%.9g,%.9g,%.9g), "
        "pre_r1r1/r1r2/r2r2=(%.9g,%.9g,%.9g), "
        "pre_radicand/root/k=(%.9g,%.9g,%.9g), "
        "radicand_bits=0x%08x, "
        "bad_post_crd=(%.9g,%.9g,%.9g;0x%08x,0x%08x,0x%08x), "
        "bad_post_vel=(%.9g,%.9g,%.9g;0x%08x,0x%08x,0x%08x), "
        "writer_post_crd_i=(%.9g,%.9g,%.9g), "
        "writer_post_crd_j=(%.9g,%.9g,%.9g), "
        "writer_post_vel_i=(%.9g,%.9g,%.9g), "
        "writer_post_vel_j=(%.9g,%.9g,%.9g))\n",
        bad_field, invalid_atoms, owned_atom_numbers, first_invalid,
        global_atoms[first_invalid], writer_pair, writer_matches, writer.atom_i,
        writer.atom_j, writer.global_i, writer.global_j, writer.constant_r,
        writer.constrain_k, writer.mass_i, writer.mass_j, dt, v_factor,
        x_factor, writer.crd_i.x, writer.crd_i.y, writer.crd_i.z,
        writer.crd_j.x, writer.crd_j.y, writer.crd_j.z, writer.vel_i.x,
        writer.vel_i.y, writer.vel_i.z, writer.vel_j.x, writer.vel_j.y,
        writer.vel_j.z, writer.frc_i.x, writer.frc_i.y, writer.frc_i.z,
        writer.frc_j.x, writer.frc_j.y, writer.frc_j.z, writer.r1.x,
        writer.r1.y, writer.r1.z, writer.r2.x, writer.r2.y, writer.r2.z,
        writer.r1r1, writer.r1r2, writer.r2r2, writer.radicand, writer.root,
        writer.k, Settle_Diagnostic_Float_Bits(writer.radicand), bad_crd.x,
        bad_crd.y, bad_crd.z, Settle_Diagnostic_Float_Bits(bad_crd.x),
        Settle_Diagnostic_Float_Bits(bad_crd.y),
        Settle_Diagnostic_Float_Bits(bad_crd.z), bad_vel.x, bad_vel.y,
        bad_vel.z, Settle_Diagnostic_Float_Bits(bad_vel.x),
        Settle_Diagnostic_Float_Bits(bad_vel.y),
        Settle_Diagnostic_Float_Bits(bad_vel.z), post_crd_i.x, post_crd_i.y,
        post_crd_i.z, post_crd_j.x, post_crd_j.y, post_crd_j.z, post_vel_i.x,
        post_vel_i.y, post_vel_i.z, post_vel_j.x, post_vel_j.y, post_vel_j.z);
}

static __global__ void Sum_Virial_Tensor_To_Stress(int N,
                                                   LTMatrix3* virial_tensor,
                                                   LTMatrix3* stress,
                                                   const LTMatrix3 rcell)
{
    LTMatrix3 virial = {0, 0, 0, 0, 0, 0};
    float factor = rcell.a11 * rcell.a22 * rcell.a33;
#ifdef USE_GPU
    int tid = blockDim.x * blockDim.y * blockIdx.x + blockDim.y * threadIdx.x +
              threadIdx.y;
    if (tid < N)
    {
        virial = virial + virial_tensor[tid];
    }
#else
    float v11 = 0.0f, v21 = 0.0f, v22 = 0.0f;
    float v31 = 0.0f, v32 = 0.0f, v33 = 0.0f;
#pragma omp parallel for reduction(+ : v11, v21, v22, v31, v32, v33)
    for (int tid = 0; tid < N; tid++)
    {
        v11 += virial_tensor[tid].a11;
        v21 += virial_tensor[tid].a21;
        v22 += virial_tensor[tid].a22;
        v31 += virial_tensor[tid].a31;
        v32 += virial_tensor[tid].a32;
        v33 += virial_tensor[tid].a33;
    }
    virial = {v11, v21, v22, v31, v32, v33};
#endif
    virial = factor * virial;
    Warp_Sum_To(stress, virial, warpSize);
}

void SETTLE::Initial(CONTROLLER* controller, CONSTRAIN* constrain,
                     float* h_mass, const char* module_name)
{
    if (module_name == NULL)
    {
        strcpy(this->module_name, "settle");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }
    if (controller->Command_Exist("settle_disable") &&
        controller->Get_Bool("settle_disable", "Main_Initial"))
    {
        return;
    }
    if (constrain->constrain_pair_numbers > 0)
    {
        this->constrain = constrain;
        controller->printf("START INITIALIZING SETTLE:\n");
        // 遍历搜出constrain里的三角形
        int* linker_numbers = NULL;
        int* linker_atoms = NULL;
        float* link_r = NULL;
        Malloc_Safely((void**)&linker_numbers,
                      sizeof(int) * constrain->atom_numbers);
        Malloc_Safely((void**)&linker_atoms,
                      2 * sizeof(int) * constrain->atom_numbers);
        Malloc_Safely((void**)&link_r,
                      3 * sizeof(float) * constrain->atom_numbers);
        for (int i = 0; i < constrain->atom_numbers; i++)
        {
            linker_numbers[i] = 0;
        }
        int atom_i, atom_j;
        CONSTRAIN_PAIR pair;
        for (int i = 0; i < constrain->constrain_pair_numbers; i++)
        {
            pair = constrain->h_constrain_pair[i];
            atom_i = pair.atom_i_serial;
            atom_j = pair.atom_j_serial;

            if (linker_numbers[atom_i] < 2 && linker_numbers[atom_j] < 2)
            {
                linker_atoms[2 * atom_i + linker_numbers[atom_i]] = atom_j;
                linker_atoms[2 * atom_j + linker_numbers[atom_j]] = atom_i;
                link_r[3 * atom_i + linker_numbers[atom_i]] = pair.constant_r;
                link_r[3 * atom_j + linker_numbers[atom_j]] = pair.constant_r;
                linker_numbers[atom_i]++;
                linker_numbers[atom_j]++;
            }
            else
            {
                linker_numbers[atom_i] = 3;
                linker_numbers[atom_j] = 3;
            }
        }
        triangle_numbers = 0;
        pair_numbers = 0;
        for (int i = 0; i < constrain->atom_numbers; i++)
        {
            if (linker_numbers[i] == 2)
            {
                atom_i = linker_atoms[2 * i];
                atom_j = linker_atoms[2 * i + 1];
                if (linker_numbers[atom_i] == 2 &&
                    linker_numbers[atom_j] == 2 &&
                    ((linker_atoms[2 * atom_i] == i &&
                      linker_atoms[2 * atom_i + 1] == atom_j) ||
                     (linker_atoms[2 * atom_i + 1] == i &&
                      linker_atoms[2 * atom_i] == atom_j)))
                {
                    triangle_numbers++;
                    linker_numbers[atom_i] = -2;
                    linker_numbers[atom_j] = -2;
                    if (linker_atoms[2 * atom_i + 1] == atom_j)
                    {
                        link_r[3 * i + 2] = link_r[3 * atom_i + 1];
                    }
                    else
                    {
                        link_r[3 * i + 2] = link_r[3 * atom_i];
                    }
                }
                else
                {
                    linker_numbers[i] = 3;
                    linker_numbers[atom_i] = 3;
                    linker_numbers[atom_j] = 3;
                }
            }
            else if (linker_numbers[i] == 1)
            {
                atom_i = linker_atoms[2 * i];
                if (linker_numbers[atom_i] == 1)
                {
                    pair_numbers++;
                    linker_numbers[atom_i] = -1;
                }
                else
                {
                    linker_numbers[i] = 3;
                    linker_numbers[atom_i] = 3;
                }
            }
        }
        controller->printf("    rigid triangle numbers is %d\n",
                           triangle_numbers);
        controller->printf("    rigid pair numbers is %d\n", pair_numbers);
        if (triangle_numbers > 0 || pair_numbers > 0)
        {
            Malloc_Safely((void**)&h_triangles,
                          sizeof(CONSTRAIN_TRIANGLE) * triangle_numbers);
            Malloc_Safely((void**)&h_pairs,
                          sizeof(CONSTRAIN_PAIR) * pair_numbers);

            Device_Malloc_Safely((void**)&last_triangle_BA,
                                 sizeof(VECTOR) * triangle_numbers);
            Device_Malloc_Safely((void**)&last_triangle_CA,
                                 sizeof(VECTOR) * triangle_numbers);
            Device_Malloc_Safely((void**)&last_pair_AB,
                                 sizeof(VECTOR) * pair_numbers);
            Device_Malloc_Safely(
                (void**)&virial_tensor,
                sizeof(LTMatrix3) * (triangle_numbers + pair_numbers));
            int triangle_i = 0;
            int pair_i = 0;
            for (int i = 0; i < constrain->atom_numbers; i++)
            {
                if (linker_numbers[i] == 2)
                {
                    linker_numbers[i] = -2;
                    atom_i = linker_atoms[2 * i];
                    atom_j = linker_atoms[2 * i + 1];
                    h_triangles[triangle_i].atom_A = i;
                    h_triangles[triangle_i].atom_B = atom_i;
                    h_triangles[triangle_i].atom_C = atom_j;
                    Get_Rabcde_From_SSS(
                        link_r[3 * i], link_r[3 * i + 1], link_r[3 * i + 2],
                        h_mass[i], h_mass[atom_i], h_mass[atom_j],
                        h_triangles[triangle_i].ra, h_triangles[triangle_i].rb,
                        h_triangles[triangle_i].rc, h_triangles[triangle_i].rd,
                        h_triangles[triangle_i].re);
                    triangle_i++;
                }
                if (linker_numbers[i] == 1)
                {
                    atom_j = linker_atoms[2 * i];
                    linker_numbers[i] = -1;
                    h_pairs[pair_i].atom_i_serial = i;
                    h_pairs[pair_i].atom_j_serial = atom_j;
                    h_pairs[pair_i].constant_r = link_r[3 * i];
                    h_pairs[pair_i].constrain_k =
                        1.0f / (h_mass[i] + h_mass[atom_j]);
                    pair_i++;
                }
            }

            Device_Malloc_And_Copy_Safely(
                (void**)&d_triangles, h_triangles,
                sizeof(CONSTRAIN_TRIANGLE) * triangle_numbers);
            Device_Malloc_And_Copy_Safely(
                (void**)&d_pairs, h_pairs,
                sizeof(CONSTRAIN_PAIR) * pair_numbers);

            Device_Malloc_Safely((void**)&d_triangles_local,
                                 sizeof(CONSTRAIN_TRIANGLE) * triangle_numbers);
            Device_Malloc_Safely((void**)&d_pairs_local,
                                 sizeof(CONSTRAIN_PAIR) * pair_numbers);
            Device_Malloc_Safely((void**)&d_num_triangle_local, sizeof(int));
            Device_Malloc_Safely((void**)&d_num_pair_local, sizeof(int));
            Device_Malloc_Safely((void**)&d_delta_vel_local,
                                 sizeof(VECTOR) * constrain->atom_numbers);

            // 原来的重塑
            int new_constrain_pair_numbers = constrain->constrain_pair_numbers -
                                             3 * triangle_numbers -
                                             pair_numbers;
            int new_pair_i = 0;

            CONSTRAIN_PAIR* new_h_constrain_pair = NULL;
            Malloc_Safely((void**)&new_h_constrain_pair,
                          sizeof(CONSTRAIN_PAIR) * new_constrain_pair_numbers);

            for (int i = 0; i < constrain->constrain_pair_numbers; i++)
            {
                pair = constrain->h_constrain_pair[i];
                atom_i = pair.atom_i_serial;
                ;
                if (linker_numbers[atom_i] > 0)
                {
                    new_h_constrain_pair[new_pair_i] = pair;
                    new_pair_i++;
                }
            }
            constrain->constrain_pair_numbers = new_constrain_pair_numbers;

            Free_Host_And_Device_Pointer((void**)&constrain->h_constrain_pair,
                                         (void**)&constrain->d_constrain_pair);

            constrain->h_constrain_pair = new_h_constrain_pair;
            Device_Malloc_And_Copy_Safely(
                (void**)&constrain->d_constrain_pair,
                constrain->h_constrain_pair,
                sizeof(CONSTRAIN_PAIR) * new_constrain_pair_numbers);

            controller->printf(
                "    remaining simple constrain pair numbers is %d\n",
                new_pair_i);
            for (int i = 0; i < constrain->constrain_pair_numbers; i++)
            {
                pair = constrain->h_constrain_pair[i];
                atom_i = pair.atom_i_serial;
            }
            free(linker_numbers);
            free(linker_atoms);
            free(link_r);
            is_initialized = 1;
            controller->printf("END INITIALIZING SETTLE\n\n");
        }
        else
        {
            controller->printf("SETTLE IS NOT INITIALIZED\n\n");
        }
    }
    else
    {
        controller->printf("SETTLE IS NOT INITIALIZED\n\n");
    }
}

void SETTLE::Remember_Last_Coordinates(const VECTOR* crd, const LTMatrix3 cell,
                                       const LTMatrix3 rcell)
{
    if (!is_initialized) return;

    Launch_Device_Kernel(remember_pair_AB,
                         (num_pair_local + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, num_pair_local,
                         d_pairs_local, crd, cell, rcell, last_pair_AB);

    Launch_Device_Kernel(
        remember_triangle_BA_CA,
        (num_triangle_local + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, num_triangle_local,
        d_triangles_local, crd, cell, rcell, last_triangle_BA,
        last_triangle_CA);
}

static __global__ void get_local_device(int triangle_numbers, int pair_numbers,
                                        const CONSTRAIN_TRIANGLE* d_triangles,
                                        CONSTRAIN_TRIANGLE* d_triangles_local,
                                        const CONSTRAIN_PAIR* d_pairs,
                                        CONSTRAIN_PAIR* d_pairs_local,
                                        const int* atom_local_id,
                                        const char* atom_local_label,
                                        int* d_num_triangle_local,
                                        int* d_num_pair_local)
{
#ifdef USE_GPU
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid != 0) return;
#endif
    d_num_triangle_local[0] = 0;
    for (int i = 0; i < triangle_numbers; i++)
    {
        int atom_a = d_triangles[i].atom_A;
        int atom_b = d_triangles[i].atom_B;
        int atom_c = d_triangles[i].atom_C;
        if (atom_local_label[atom_a])
        {
            d_triangles_local[d_num_triangle_local[0]] = d_triangles[i];
            d_triangles_local[d_num_triangle_local[0]].atom_A =
                atom_local_id[atom_a];
            d_triangles_local[d_num_triangle_local[0]].atom_B =
                atom_local_id[atom_b];
            d_triangles_local[d_num_triangle_local[0]].atom_C =
                atom_local_id[atom_c];
            d_num_triangle_local[0] += 1;
        }
    }
    d_num_pair_local[0] = 0;
    for (int i = 0; i < pair_numbers; i++)
    {
        int atom_a = d_pairs[i].atom_i_serial;
        int atom_b = d_pairs[i].atom_j_serial;
        if (atom_local_label[atom_a])
        {
            d_pairs_local[d_num_pair_local[0]] = d_pairs[i];
            d_pairs_local[d_num_pair_local[0]].atom_i_serial =
                atom_local_id[atom_a];
            d_pairs_local[d_num_pair_local[0]].atom_j_serial =
                atom_local_id[atom_b];
            d_num_pair_local[0] += 1;
        }
    }
}

void SETTLE::Get_Local(const int* atom_local_id, const char* atom_local_label,
                       const int local_atom_numbers)
{
    if (!is_initialized) return;
    this->local_atom_numbers = local_atom_numbers;
    num_triangle_local = 0;
    num_pair_local = 0;
    Launch_Device_Kernel(get_local_device, 1, 1, 0, NULL, triangle_numbers,
                         pair_numbers, d_triangles, d_triangles_local, d_pairs,
                         d_pairs_local, atom_local_id, atom_local_label,
                         d_num_triangle_local, d_num_pair_local);
    deviceMemcpy(&num_triangle_local, d_num_triangle_local, sizeof(int),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&num_pair_local, d_num_pair_local, sizeof(int),
                 deviceMemcpyDeviceToHost);
}

void SETTLE::Do_SETTLE(CONTROLLER* controller, const int* atom_local,
                       const int active_atom_numbers, const VECTOR* frc,
                       const float* d_mass, VECTOR* crd, const LTMatrix3 cell,
                       const LTMatrix3 rcell, VECTOR* vel,
                       const int need_pressure, LTMatrix3* d_stress)
{
    if (!is_initialized) return;
    const bool diagnose_pair_state = diagnose_pair_state_once;
    std::vector<SETTLE_PAIR_NUMERIC_DIAGNOSTIC> pair_preflight;
    if (diagnose_pair_state)
    {
        diagnose_pair_state_once = false;
        Preflight_Settle_Pairs_Once(
            controller, num_pair_local, local_atom_numbers, active_atom_numbers,
            d_pairs_local, atom_local, frc, d_mass, crd, vel, cell, rcell,
            last_pair_AB, constrain->dt, constrain->v_factor,
            constrain->x_factor, &pair_preflight);
    }
    Launch_Device_Kernel(settle_pair,
                         (num_pair_local + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, num_pair_local,
                         d_pairs_local, d_mass, crd, cell, rcell, last_pair_AB,
                         constrain->dt, constrain->v_factor,
                         constrain->x_factor, vel,
                         virial_tensor + triangle_numbers);

    if (diagnose_pair_state)
    {
        Diagnose_Actual_Settle_Pair_Output_Once(
            controller, local_atom_numbers, atom_local, crd, vel, constrain->dt,
            constrain->v_factor, constrain->x_factor, pair_preflight);
    }

    Launch_Device_Kernel(
        settle_triangle,
        (num_triangle_local + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, num_triangle_local,
        d_triangles_local, d_mass, crd, cell, rcell, last_triangle_BA,
        last_triangle_CA, constrain->dt, constrain->v_factor,
        constrain->x_factor, vel, virial_tensor);

    if (need_pressure)
    {
        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        Launch_Device_Kernel(Sum_Virial_Tensor_To_Stress,
                             (num_triangle_local + num_pair_local +
                              CONTROLLER::device_max_thread - 1) /
                                 CONTROLLER::device_max_thread,
                             blockSize, 0, NULL, num_triangle_local,
                             virial_tensor, d_stress, rcell);
        Launch_Device_Kernel(Sum_Virial_Tensor_To_Stress,
                             (num_triangle_local + num_pair_local +
                              CONTROLLER::device_max_thread - 1) /
                                 CONTROLLER::device_max_thread,
                             blockSize, 0, NULL, num_pair_local,
                             virial_tensor + triangle_numbers, d_stress, rcell);
    }
}

static __device__ __host__ __forceinline__ bool
compute_velocity_constraint_correction_settle(
    const int atom_i, const int atom_j, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const float* mass_inverse, const VECTOR* vel,
    VECTOR* correction_i, VECTOR* correction_j, const float relative_tolerance,
    bool* constraint_violated)
{
    if (constraint_violated != NULL) *constraint_violated = false;
    float mass_i_inverse = mass_inverse[atom_i];
    float mass_j_inverse = mass_inverse[atom_j];
    if (mass_i_inverse == 0.0f && mass_j_inverse == 0.0f) return false;

    VECTOR dr =
        Get_Periodic_Displacement(crd[atom_i], crd[atom_j], cell, rcell);
    float dr2 = dr * dr;
    if (dr2 < 1e-12f)
    {
        if (constraint_violated != NULL) *constraint_violated = true;
        return false;
    }

    const VECTOR velocity_i = vel[atom_i];
    const VECTOR velocity_j = vel[atom_j];
    VECTOR v_diff = velocity_i - velocity_j;
    float denom = (mass_i_inverse + mass_j_inverse) * dr2;
    if (denom < 1e-20f)
    {
        if (constraint_violated != NULL) *constraint_violated = true;
        return false;
    }

    const float projection = dr * v_diff;
    if (constraint_violated != NULL)
    {
        const float tolerance = Velocity_Constraint_Residual_Tolerance(
            dr2, velocity_i, velocity_j, v_diff, relative_tolerance);
        *constraint_violated = fabsf(projection) > tolerance;
    }
    float lambda = projection / denom;
    correction_i[0] = (-mass_i_inverse * lambda) * dr;
    correction_j[0] = (mass_j_inverse * lambda) * dr;
    return true;
}

static __global__ void project_velocity_to_settle_pairs(
    const int pair_numbers, const CONSTRAIN_PAIR* pairs, const VECTOR* crd,
    const LTMatrix3 cell, const LTMatrix3 rcell, const float* mass_inverse,
    const VECTOR* vel, VECTOR* delta_vel, const float relative_tolerance,
    int* violation)
{
#ifdef USE_GPU
    int pair_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (pair_i < pair_numbers)
#else
#pragma omp parallel for
    for (int pair_i = 0; pair_i < pair_numbers; pair_i++)
#endif
    {
        CONSTRAIN_PAIR cp = pairs[pair_i];
        VECTOR correction_i, correction_j;
        bool constraint_violated = false;
        if (compute_velocity_constraint_correction_settle(
                cp.atom_i_serial, cp.atom_j_serial, crd, cell, rcell,
                mass_inverse, vel, &correction_i, &correction_j,
                relative_tolerance,
                violation != NULL ? &constraint_violated : NULL))
        {
            atomicAdd(&delta_vel[cp.atom_i_serial].x, correction_i.x);
            atomicAdd(&delta_vel[cp.atom_i_serial].y, correction_i.y);
            atomicAdd(&delta_vel[cp.atom_i_serial].z, correction_i.z);
            atomicAdd(&delta_vel[cp.atom_j_serial].x, correction_j.x);
            atomicAdd(&delta_vel[cp.atom_j_serial].y, correction_j.y);
            atomicAdd(&delta_vel[cp.atom_j_serial].z, correction_j.z);
        }
        if (constraint_violated && violation != NULL) atomicExch(violation, 1);
    }
}

static __global__ void project_velocity_to_settle_triangles(
    const int triangle_numbers, const CONSTRAIN_TRIANGLE* triangles,
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    const float* mass_inverse, const VECTOR* vel, VECTOR* delta_vel,
    const float relative_tolerance, int* violation)
{
#ifdef USE_GPU
    int tri_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (tri_i < triangle_numbers)
#else
#pragma omp parallel for
    for (int tri_i = 0; tri_i < triangle_numbers; tri_i++)
#endif
    {
        CONSTRAIN_TRIANGLE tri = triangles[tri_i];
        VECTOR correction_i, correction_j;
        bool constraint_violated = false;
        if (compute_velocity_constraint_correction_settle(
                tri.atom_A, tri.atom_B, crd, cell, rcell, mass_inverse, vel,
                &correction_i, &correction_j, relative_tolerance,
                violation != NULL ? &constraint_violated : NULL))
        {
            atomicAdd(&delta_vel[tri.atom_A].x, correction_i.x);
            atomicAdd(&delta_vel[tri.atom_A].y, correction_i.y);
            atomicAdd(&delta_vel[tri.atom_A].z, correction_i.z);
            atomicAdd(&delta_vel[tri.atom_B].x, correction_j.x);
            atomicAdd(&delta_vel[tri.atom_B].y, correction_j.y);
            atomicAdd(&delta_vel[tri.atom_B].z, correction_j.z);
        }
        if (constraint_violated && violation != NULL) atomicExch(violation, 1);
        if (compute_velocity_constraint_correction_settle(
                tri.atom_A, tri.atom_C, crd, cell, rcell, mass_inverse, vel,
                &correction_i, &correction_j, relative_tolerance,
                violation != NULL ? &constraint_violated : NULL))
        {
            atomicAdd(&delta_vel[tri.atom_A].x, correction_i.x);
            atomicAdd(&delta_vel[tri.atom_A].y, correction_i.y);
            atomicAdd(&delta_vel[tri.atom_A].z, correction_i.z);
            atomicAdd(&delta_vel[tri.atom_C].x, correction_j.x);
            atomicAdd(&delta_vel[tri.atom_C].y, correction_j.y);
            atomicAdd(&delta_vel[tri.atom_C].z, correction_j.z);
        }
        if (constraint_violated && violation != NULL) atomicExch(violation, 1);
        if (compute_velocity_constraint_correction_settle(
                tri.atom_B, tri.atom_C, crd, cell, rcell, mass_inverse, vel,
                &correction_i, &correction_j, relative_tolerance,
                violation != NULL ? &constraint_violated : NULL))
        {
            atomicAdd(&delta_vel[tri.atom_B].x, correction_i.x);
            atomicAdd(&delta_vel[tri.atom_B].y, correction_i.y);
            atomicAdd(&delta_vel[tri.atom_B].z, correction_i.z);
            atomicAdd(&delta_vel[tri.atom_C].x, correction_j.x);
            atomicAdd(&delta_vel[tri.atom_C].y, correction_j.y);
            atomicAdd(&delta_vel[tri.atom_C].z, correction_j.z);
        }
        if (constraint_violated && violation != NULL) atomicExch(violation, 1);
    }
}

static __global__ void apply_settle_velocity_correction(
    const int local_atom_numbers, VECTOR* vel, VECTOR* crd,
    const VECTOR* delta_vel, const float velocity_factor,
    const float coordinate_factor)
{
#ifdef USE_GPU
    int atom_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (atom_i < local_atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < local_atom_numbers; atom_i++)
#endif
    {
        VECTOR delta = velocity_factor * delta_vel[atom_i];
        vel[atom_i] = vel[atom_i] + delta;
        if (coordinate_factor != 0.0f)
        {
            crd[atom_i] = crd[atom_i] + coordinate_factor * delta;
        }
    }
}

bool SETTLE::Project_Velocity_To_Constraint_Manifold(VECTOR* vel, VECTOR* crd,
                                                     const float* mass_inverse,
                                                     const LTMatrix3 cell,
                                                     const LTMatrix3 rcell,
                                                     bool update_coordinates)
{
    if (!is_initialized || local_atom_numbers <= 0) return true;
    bool has_settle_constraints = num_triangle_local > 0 || num_pair_local > 0;
    if (!has_settle_constraints) return true;

    constexpr int legacy_projection_iterations = 8;
    constexpr int maximum_velocity_only_iterations = 512;
    constexpr float relative_tolerance = 1.0e-5f;
    const int projection_iterations = update_coordinates
                                          ? legacy_projection_iterations
                                          : maximum_velocity_only_iterations;
    // SETTLE extracts disjoint pairs and triangles.  A pair can be projected
    // directly; a triangle has local constraint degree two, so damping its
    // parallel Jacobi update by one half is sufficient and does not couple
    // convergence to unrelated constraints left for SHAKE.
    const int settle_constraint_degree = num_triangle_local > 0 ? 2 : 1;
    const float velocity_factor =
        update_coordinates
            ? 1.0f
            : 1.0f / static_cast<float>(settle_constraint_degree);
    int* d_violation = NULL;
    if (!update_coordinates &&
        !Device_Malloc_Safely((void**)&d_violation, sizeof(int)))
        return false;
    bool converged = update_coordinates;
    for (int iter = 0; iter < projection_iterations; ++iter)
    {
        if (!update_coordinates) deviceMemset(d_violation, 0, sizeof(int));
        deviceMemset(d_delta_vel_local, 0, sizeof(VECTOR) * local_atom_numbers);
        if (num_pair_local > 0 && d_pairs_local != NULL)
        {
            Launch_Device_Kernel(
                project_velocity_to_settle_pairs,
                (num_pair_local + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, num_pair_local,
                d_pairs_local, crd, cell, rcell, mass_inverse, vel,
                d_delta_vel_local, relative_tolerance, d_violation);
        }
        if (num_triangle_local > 0 && d_triangles_local != NULL)
        {
            Launch_Device_Kernel(
                project_velocity_to_settle_triangles,
                (num_triangle_local + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, num_triangle_local,
                d_triangles_local, crd, cell, rcell, mass_inverse, vel,
                d_delta_vel_local, relative_tolerance, d_violation);
        }
        if (!update_coordinates)
        {
            int violation = 0;
            deviceMemcpy(&violation, d_violation, sizeof(int),
                         deviceMemcpyDeviceToHost);
            if (violation == 0)
            {
                converged = true;
                break;
            }
        }
        Launch_Device_Kernel(
            apply_settle_velocity_correction,
            (local_atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, local_atom_numbers, vel,
            crd, d_delta_vel_local, velocity_factor,
            update_coordinates ? 0.5f * constrain->dt : 0.0f);
    }
    if (d_violation != NULL) deviceFree(d_violation);
    return converged;
}

void SETTLE::update_ug_connectivity(CONECT* connectivity)
{
    if (!is_initialized) return;
    for (int i = 0; i < pair_numbers; i++)
    {
        CONSTRAIN_PAIR p = h_pairs[i];
        (*connectivity)[p.atom_i_serial].insert(p.atom_j_serial);
        (*connectivity)[p.atom_j_serial].insert(p.atom_i_serial);
    }
    for (int i = 0; i < triangle_numbers; i++)
    {
        CONSTRAIN_TRIANGLE t = h_triangles[i];
        (*connectivity)[t.atom_A].insert(t.atom_B);
        (*connectivity)[t.atom_A].insert(t.atom_C);
        (*connectivity)[t.atom_B].insert(t.atom_A);
        (*connectivity)[t.atom_B].insert(t.atom_C);
        (*connectivity)[t.atom_C].insert(t.atom_A);
        (*connectivity)[t.atom_C].insert(t.atom_B);
    }
}
