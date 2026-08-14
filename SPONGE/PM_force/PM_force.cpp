#include "PM_force.h"

#include <cerrno>
#include <cfloat>
#include <climits>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <new>
#include <string>

#include "../Lennard_Jones_force/pair_activity.h"
#include "pme_excluded_pair.h"
/*
    2025-10-14 SPONGE Particle Mesh算法
    目前支持单进程Particle-Mesh-Ewald 与 PMC-IZ
   计算，在头文件预留了MPI-FFT多进程接口
    单进程下，PM进程与PP进程设置为同一个进程；多进程下PM独享一个进程

    未来可能的改进：使用PSWF作为分裂核与插值核，改进Particle
   Mesh算法，需要修改：
    - 插值核 与 修正系数计算。可以引用NUFFT相关代码
    - 近程项、长程修正项等
*/

// constants
#define PI 3.1415926f
#define INVSQRTPI 0.56418958835977f
#define TWO_DIVIDED_BY_SQRT_PI 1.1283791670218446f
static __device__ float PME_Ma[4] = {1.0 / 6.0, -0.5, 0.5, -1.0 / 6.0};
static __device__ float PME_Mb[4] = {0, 0.5, -1, 0.5};
static __device__ float PME_Mc[4] = {0, 0.5, 0, -0.5};
static __device__ float PME_Md[4] = {0, 1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0};
static __device__ float PME_dMa[4] = {0.5, -1.5, 1.5, -0.5};
static __device__ float PME_dMb[4] = {0, 1, -2, 1};
static __device__ float PME_dMc[4] = {0, 0.5, 0, -0.5};

// 计算B样条插值的递归函数 (compact-window function, or SI kernel)
static float M_(float u, int n)
{
    if (n == 2)
    {
        if (u > 2 || u < 0) return 0;
        return 1 - abs(u - 1);
    }
    else
        return u / (n - 1) * M_(u, n - 1) +
               (n - u) / (n - 1) * M_(u - 1, n - 1);
}

// 修正B样条插值的递归函数 (influence function)
static float getb(int k, int NFFT, int B_order)
{
    FFT_COMPLEX tempc, tempc2, res;
    float tempf;
    REAL(tempc2) = 0;
    IMAGINARY(tempc2) = 0;

    REAL(tempc) = 0;
    IMAGINARY(tempc) = 2 * (B_order - 1) * PI * k / NFFT;
    res = expc(tempc);

    for (int kk = 0; kk < (B_order - 1); kk++)
    {
        REAL(tempc) = 0;
        IMAGINARY(tempc) = 2.0f * PI * k / NFFT * kk;
        tempc = expc(tempc);
        tempf = M_(kk + 1, B_order);
        REAL(tempc2) += tempf * REAL(tempc);
        IMAGINARY(tempc2) += tempf * IMAGINARY(tempc);
    }
    res = divc(res, tempc2);
    return REAL(res) * REAL(res) + IMAGINARY(res) * IMAGINARY(res);
}

// PMC_IZ Method
static __global__ void Build_PMC_IZ_C(const int PME_Nfft, int fftx, int ffty,
                                      int fftz,
                                      float box_length_inverse_x_square,
                                      float box_length_inverse_y_square,
                                      float grid_length_of_z, float beta,
                                      float scalor, FFT_COMPLEX* C)
{
    SIMPLE_DEVICE_FOR(tid, PME_Nfft)
    {
        int ffta = (fftx / 2 + 1);
        int grid_x = tid % ffta;
        int grid_y = (tid % (ffta * ffty)) / ffta;
        int grid_z = tid / ffty / ffta;
        if (grid_x >= fftx / 2)
        {
            grid_x = fftx - grid_x;
        }
        if (grid_y >= ffty / 2)
        {
            grid_y = ffty - grid_y;
        }
        if (grid_z >= fftz / 2)
        {
            grid_z = fftz - grid_z;
        }
        float z = grid_length_of_z * grid_z;
        float A = 2.0f * CONSTANT_Pi *
                  sqrtf(grid_x * grid_x * box_length_inverse_x_square +
                        grid_y * grid_y * box_length_inverse_y_square);
        float AB = A / beta / 2.0f;
        float zb2 = z * beta;
        float AB_minus_zb2 = AB - zb2;
        float AB_plus_zb2 = AB + zb2;
        float temp_f =
            expf(-A * z) * (erfcf(AB_minus_zb2) +
                            expf(2.0f * A * z - AB_plus_zb2 * AB_plus_zb2) *
                                erfcxf(AB_plus_zb2));
        temp_f = temp_f / A;
        if (grid_x == 0 && grid_y == 0)
        {
            temp_f =
                2.0f / sqrtf(CONSTANT_Pi) / beta * (1.0f - expf(-zb2 * zb2)) -
                2.0f * z * erff(zb2);
        }
        REAL(C[tid]) = scalor * temp_f;
        IMAGINARY(C[tid]) = 0.;
    }
}

static __global__ void Build_PMC_IZ_BC_Final(const int Nfft, int fftx, int ffty,
                                             int fftz, const FFT_COMPLEX* C,
                                             const FFT_COMPLEX* B, float* BC)
{
    SIMPLE_DEVICE_FOR(tid, Nfft)
    {
        int fftc = fftz / 2 + 1;
        int ffta = fftx / 2 + 1;
        int zi = tid % fftc;
        int yi = (tid / fftc) % ffty;
        int xi = tid / fftc / ffty;
        if (xi >= fftx / 2)
        {
            xi = fftx - xi;
        }
        int ti = zi * ffta * ffty + yi * ffta + xi;
        float b = REAL(B[ti]);
        BC[tid] = REAL(C[ti]) / (b * b);
    }
}

static void Build_PMC_IZ_BC(CONTROLLER* controller, int fftx, int ffty,
                            int fftz, int PME_Nfft, int PME_Nall, int PME_Nin,
                            float box_length_inverse_x_square,
                            float box_length_inverse_y_square,
                            float grid_length_of_z, float beta, float scalor,
                            float** BC)
{
    Device_Malloc_Safely((void**)BC, sizeof(float) * PME_Nfft);
    FFT_SIZE_t n2d[2] = {ffty, fftx};
    FFT_RESULT result;
    FFT_HANDLE plan_2d_many_c2r, plan_3d_temp_r2c;

    result = SPONGE_FFT_WRAPPER::Make_FFT_Plan(&plan_2d_many_c2r, fftz, 2, n2d,
                                               FFT_C2R);
    if (result != FFT_SUCCESS)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMallocFailed, "Build_PMC_IZ_Efficient_Potential",
            "Reason:\n\tFail to create the batched 2D FFT plan");
    }
    FFT_SIZE_t n3d[3] = {fftz, ffty, fftx};

    result = SPONGE_FFT_WRAPPER::Make_FFT_Plan(&plan_3d_temp_r2c, 1, 3, n3d,
                                               FFT_R2C);
    if (result != FFT_SUCCESS)
    {
        controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                       "Build_PMC_IZ_Efficient_Potential",
                                       "Reason:\n\tFail to create the "
                                       "temporary 3D Real to Complex FFT plan");
    }
    FFT_COMPLEX *B, *C;
    float *d_FB, *h_FB, *FC;
    int temp_Nfft = (fftx / 2 + 1) * ffty * fftz;
    Device_Malloc_Safely((void**)&B, sizeof(FFT_COMPLEX) * temp_Nfft);
    Device_Malloc_Safely((void**)&C, sizeof(FFT_COMPLEX) * temp_Nfft);
    Device_Malloc_Safely((void**)&FC, sizeof(float) * PME_Nall);
    Malloc_Safely((void**)&h_FB, sizeof(float) * PME_Nall);

    for (int i = 0; i < PME_Nall; i = i + 1)
    {
        h_FB[i] = 0.;
    }
    float temp_b_spline[3] = {1. / 6., 2. / 3., 1. / 6.};
    for (int k = -1; k <= 1; k = k + 1)
    {
        for (int j = -1; j <= 1; j = j + 1)
        {
            for (int i = -1; i <= 1; i = i + 1)
            {
                float weight = temp_b_spline[k + 1] * temp_b_spline[j + 1] *
                               temp_b_spline[i + 1];
                int kk, jj, ii;
                if (k < 0)
                {
                    kk = k + fftz;
                }
                else
                {
                    kk = k;
                }
                if (j < 0)
                {
                    jj = j + ffty;
                }
                else
                {
                    jj = j;
                }
                if (i < 0)
                {
                    ii = i + fftx;
                }
                else
                {
                    ii = i;
                }
                h_FB[ii + jj * fftx + kk * fftx * ffty] = weight;
            }
        }
    }
    Device_Malloc_And_Copy_Safely((void**)&d_FB, h_FB,
                                  sizeof(float) * PME_Nall);
    SPONGE_FFT_WRAPPER::R2C(plan_3d_temp_r2c, d_FB, B);
    Launch_Device_Kernel(Build_PMC_IZ_C,
                         (temp_Nfft + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, temp_Nfft,
                         fftx, ffty, fftz, box_length_inverse_x_square,
                         box_length_inverse_y_square, grid_length_of_z, beta,
                         scalor, C);

    SPONGE_FFT_WRAPPER::C2R(plan_2d_many_c2r, C, FC);
    SPONGE_FFT_WRAPPER::R2C(plan_3d_temp_r2c, FC, C);
    Launch_Device_Kernel(Build_PMC_IZ_BC_Final,
                         (PME_Nfft + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, PME_Nfft, fftx,
                         ffty, fftz, C, B, BC[0]);

    Free_Single_Device_Pointer((void**)&FC);
    Free_Single_Device_Pointer((void**)&C);
    Free_Single_Device_Pointer((void**)&B);
    Free_Host_And_Device_Pointer((void**)&h_FB, (void**)&d_FB);
    SPONGE_FFT_WRAPPER::Destroy_FFT_Plan(&plan_2d_many_c2r);
    SPONGE_FFT_WRAPPER::Destroy_FFT_Plan(&plan_3d_temp_r2c);
}

static void Set_PME_Validation_Error(std::string* error,
                                     const std::string& message)
{
    if (error != nullptr) *error = message;
}

// These parsers deliberately accept only ordinary base-10 configuration
// tokens.  strtod/strtoll alone accept leading whitespace, hexadecimal
// floats, NaN/Inf, and prefixes followed by garbage, none of which is a valid
// PME control value.
static bool Is_Strict_PME_Decimal(const char* token)
{
    if (token == nullptr || token[0] == '\0') return false;
    const char* cursor = token;
    if (*cursor == '+' || *cursor == '-') ++cursor;

    bool has_digit = false;
    while (*cursor >= '0' && *cursor <= '9')
    {
        has_digit = true;
        ++cursor;
    }
    if (*cursor == '.')
    {
        ++cursor;
        while (*cursor >= '0' && *cursor <= '9')
        {
            has_digit = true;
            ++cursor;
        }
    }
    if (!has_digit) return false;

    if (*cursor == 'e' || *cursor == 'E')
    {
        ++cursor;
        if (*cursor == '+' || *cursor == '-') ++cursor;
        const char* exponent_start = cursor;
        while (*cursor >= '0' && *cursor <= '9') ++cursor;
        if (cursor == exponent_start) return false;
    }
    return *cursor == '\0';
}

static bool Try_Parse_PME_Integer(const char* token, int* value,
                                  std::string* error)
{
    if (error != nullptr) error->clear();
    if (token == nullptr || token[0] == '\0')
    {
        Set_PME_Validation_Error(error, "the value is empty");
        return false;
    }
    const char* cursor = token;
    if (*cursor == '+' || *cursor == '-') ++cursor;
    const char* digit_start = cursor;
    while (*cursor >= '0' && *cursor <= '9') ++cursor;
    if (cursor == digit_start || *cursor != '\0')
    {
        Set_PME_Validation_Error(error,
                                 "expected one complete base-10 integer token");
        return false;
    }

    errno = 0;
    char* end = nullptr;
    const long long parsed = std::strtoll(token, &end, 10);
    if (errno == ERANGE || end == token || *end != '\0' || parsed < INT_MIN ||
        parsed > INT_MAX)
    {
        Set_PME_Validation_Error(error,
                                 "the integer is outside the supported int "
                                 "range");
        return false;
    }
    *value = static_cast<int>(parsed);
    return true;
}

static bool Try_Parse_PME_Finite_Double(const char* token, double* value,
                                        std::string* error)
{
    if (error != nullptr) error->clear();
    if (!Is_Strict_PME_Decimal(token))
    {
        Set_PME_Validation_Error(
            error, "expected one complete decimal floating-point token");
        return false;
    }

    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(token, &end);
    if (errno == ERANGE || end == token || *end != '\0' ||
        !Double_Memory_Is_Finite(&parsed))
    {
        Set_PME_Validation_Error(
            error, "the floating-point value must be finite and in range");
        return false;
    }
    *value = parsed;
    return true;
}

static bool Try_Parse_PME_Positive_Normal_Float(const char* token, float* value,
                                                std::string* error)
{
    double parsed = 0.0;
    if (!Try_Parse_PME_Finite_Double(token, &parsed, error)) return false;
    if (!(parsed > 0.0) || parsed > FLT_MAX)
    {
        Set_PME_Validation_Error(
            error, "the value is outside the positive float range");
        return false;
    }
    const float checked = static_cast<float>(parsed);
    if (!Float_Memory_Is_Normal(&checked) || !(checked > 0.0f))
    {
        Set_PME_Validation_Error(
            error,
            "the value is not representable as a positive normal "
            "float");
        return false;
    }
    *value = checked;
    return true;
}

static bool Try_Validate_PME_Process_Count(int mpi_size, int cc_mpi_size,
                                           int pm_mpi_size, std::string* error)
{
    if (error != nullptr) error->clear();
    if (mpi_size <= 0)
    {
        Set_PME_Validation_Error(error,
                                 "MPI_size must be a positive rank count");
        return false;
    }
    if (cc_mpi_size < 0 || cc_mpi_size >= mpi_size)
    {
        Set_PME_Validation_Error(
            error,
            "CC_MPI_size must be nonnegative and smaller than "
            "MPI_size");
        return false;
    }
    if (pm_mpi_size < 0)
    {
        Set_PME_Validation_Error(error, "PM.MPI_size must be nonnegative");
        return false;
    }
    if (mpi_size == 1)
    {
        if (cc_mpi_size != 0 || pm_mpi_size > 1)
        {
            Set_PME_Validation_Error(
                error,
                "a one-rank run permits only PM.MPI_size 0 or 1 and "
                "no dedicated CC rank");
            return false;
        }
        return true;
    }

    const long long available_for_pp_and_pm =
        static_cast<long long>(mpi_size) - cc_mpi_size;
    if (pm_mpi_size >= available_for_pp_and_pm)
    {
        Set_PME_Validation_Error(
            error,
            "PM.MPI_size must leave at least one PP rank after "
            "dedicated PM and CC ranks are assigned");
        return false;
    }
    return true;
}

static bool Try_Validate_PME_Implemented_Process_Roles(int cc_mpi_size,
                                                       std::string* error)
{
    if (error != nullptr) error->clear();
    if (cc_mpi_size != 0)
    {
        Set_PME_Validation_Error(
            error,
            "PME does not yet implement a distinct CC-rank role; "
            "CC_MPI_size must be 0");
        return false;
    }
    return true;
}

static bool PME_Split_Product_Equals(INT_VECTOR split, int rank_count)
{
    if (split.int_x <= 0 || split.int_y <= 0 || split.int_z <= 0 ||
        rank_count <= 0)
        return false;
    const long long xy = static_cast<long long>(split.int_x) * split.int_y;
    // Once xy is bounded by the int rank count, multiplying by the final int
    // dimension cannot overflow long long.
    return xy <= rank_count && xy * split.int_z == rank_count;
}

// Finds the same 2/3/5/7-smooth, four-aligned dimensions as the legacy
// Get_Fft_Patameter search, but enumerates a finite set.  This avoids both the
// float-to-int overflow and the unbounded `tempi += 4` loop on hostile input.
static bool Try_Get_PME_Fft_Parameter(double mesh_length, int* dimension,
                                      std::string* error)
{
    if (error != nullptr) error->clear();
    if (!Double_Memory_Is_Finite(&mesh_length) || !(mesh_length > 0.0))
    {
        Set_PME_Validation_Error(
            error,
            "box length divided by grid_spacing must be finite and "
            "positive");
        return false;
    }
    const double shifted = mesh_length + 3.0;
    if (!Double_Memory_Is_Finite(&shifted) ||
        shifted > static_cast<double>(INT_MAX))
    {
        Set_PME_Validation_Error(
            error,
            "the automatic FFT dimension exceeds the supported int "
            "range");
        return false;
    }

    long long candidate = static_cast<long long>(std::ceil(shifted));
    candidate = candidate >> 2 << 2;
    if (candidate < 4) candidate = 4;
    if (candidate >= 60 && candidate <= 68)
        candidate = 64;
    else if (candidate >= 120 && candidate <= 136)
        candidate = 128;
    else if (candidate >= 240 && candidate <= 272)
        candidate = 256;
    else if (candidate >= 480 && candidate <= 544)
        candidate = 512;
    else if (candidate >= 960 && candidate <= 1088)
        candidate = 1024;

    long long best = static_cast<long long>(INT_MAX) + 1;
    for (long long factor2 = 4; factor2 <= INT_MAX;)
    {
        for (long long factor23 = factor2; factor23 <= INT_MAX;)
        {
            for (long long factor235 = factor23; factor235 <= INT_MAX;)
            {
                for (long long smooth = factor235; smooth <= INT_MAX;)
                {
                    if (smooth >= candidate && smooth < best) best = smooth;
                    if (smooth > INT_MAX / 7) break;
                    smooth *= 7;
                }
                if (factor235 > INT_MAX / 5) break;
                factor235 *= 5;
            }
            if (factor23 > INT_MAX / 3) break;
            factor23 *= 3;
        }
        if (factor2 > INT_MAX / 2) break;
        factor2 *= 2;
    }

    if (best > INT_MAX)
    {
        Set_PME_Validation_Error(
            error, "no supported 2/3/5/7-smooth FFT dimension fits in int");
        return false;
    }
    *dimension = static_cast<int>(best);
    return true;
}

struct PME_Grid_Shape
{
    int all = 0;
    int input_plane = 0;
    int transformed = 0;
    int pmc_transformed = 0;
};

static bool PME_Allocation_Count_Fits(int count, std::size_t element_size)
{
    return count >= 0 &&
           static_cast<std::size_t>(count) <=
               std::numeric_limits<std::size_t>::max() / element_size;
}

static bool PME_MPI_Byte_Count_Fits(int count, std::size_t element_size)
{
    return count >= 0 && element_size > 0 &&
           static_cast<std::size_t>(count) <=
               static_cast<std::size_t>(INT_MAX) / element_size;
}

static int PME_MPI_Byte_Count(int count, std::size_t element_size)
{
    return static_cast<int>(static_cast<std::size_t>(count) * element_size);
}

// Point-to-point PME messages use fixed protocol tags.  Source and destination
// ranks already disambiguate independent PP transfers, so encoding a rank in a
// tag needlessly makes otherwise-valid layouts depend on MPI_TAG_UB.
static constexpr int PME_MPI_TAG_DOMAIN_MIN_CORNER = 0;
static constexpr int PME_MPI_TAG_DOMAIN_MAX_CORNER = 1;
static constexpr int PME_MPI_TAG_DOMAIN_PP_COUNT = 2;
static constexpr int PME_MPI_TAG_DOMAIN_PP_RANKS = 3;
static constexpr int PME_MPI_TAG_DOMAIN_SPLIT = 4;
static constexpr int PME_MPI_TAG_DOMAIN_PM_ASSIGNMENT = 5;
static constexpr int PME_MPI_TAG_ATOM_COUNT = 6;
static constexpr int PME_MPI_TAG_COORDINATES = 7;
static constexpr int PME_MPI_TAG_CHARGES = 8;
static constexpr int PME_MPI_TAG_ATOM_IDS = 9;
static constexpr int PME_MPI_TAG_FORCES = 10;

static bool Try_Build_PME_Atom_Count_Prefix(const std::vector<int>& counts,
                                            int global_atom_count,
                                            std::vector<int>* prefixes,
                                            std::string* error)
{
    if (error != nullptr) error->clear();
    if (prefixes == nullptr)
    {
        Set_PME_Validation_Error(error, "the atom-count prefix output is null");
        return false;
    }
    if (global_atom_count <= 0 || counts.empty() ||
        counts.size() > static_cast<std::size_t>(INT_MAX))
    {
        Set_PME_Validation_Error(
            error,
            "the global atom count and PP atom-count list must be positive "
            "and fit in int");
        return false;
    }

    try
    {
        std::vector<int> checked_prefixes(counts.size(), 0);
        int prefix = 0;
        for (std::size_t i = 0; i < counts.size(); ++i)
        {
            const int count = counts[i];
            if (count < 0 || count > global_atom_count - prefix ||
                !PME_MPI_Byte_Count_Fits(count, sizeof(VECTOR)) ||
                !PME_MPI_Byte_Count_Fits(count, sizeof(float)) ||
                !PME_MPI_Byte_Count_Fits(count, sizeof(int)))
            {
                Set_PME_Validation_Error(
                    error, "PP atom-count entry " + std::to_string(i) +
                               " is negative or exceeds the remaining global "
                               "atom/MPI byte capacity");
                return false;
            }
            checked_prefixes[i] = prefix;
            prefix += count;
        }
        if (prefix != global_atom_count)
        {
            Set_PME_Validation_Error(
                error, "the PP atom counts sum to " + std::to_string(prefix) +
                           ", but the PME global atom count is " +
                           std::to_string(global_atom_count));
            return false;
        }
        *prefixes = std::move(checked_prefixes);
    }
    catch (const std::bad_alloc&)
    {
        Set_PME_Validation_Error(
            error, "host allocation failed while validating PP atom counts");
        return false;
    }
    return true;
}

static bool Try_Build_PME_Grid_Shape(int fftx, int ffty, int fftz,
                                     PME_Grid_Shape* shape, std::string* error)
{
    if (error != nullptr) error->clear();
    // Fourth-order B-spline interpolation requires at least four mesh points
    // per periodic direction.
    if (fftx < 4 || ffty < 4 || fftz < 4)
    {
        Set_PME_Validation_Error(
            error, "fftx, ffty, and fftz must each be at least 4");
        return false;
    }

    const long long xy = static_cast<long long>(fftx) * ffty;
    const long long yz = static_cast<long long>(ffty) * fftz;
    if (xy > INT_MAX || yz > INT_MAX)
    {
        Set_PME_Validation_Error(
            error, "an intermediate FFT dimension product exceeds INT_MAX");
        return false;
    }
    const long long all = xy * fftz;
    const long long transformed = xy * (static_cast<long long>(fftz) / 2 + 1);
    const long long pmc_xy = (static_cast<long long>(fftx) / 2 + 1) * ffty;
    if (all > INT_MAX || transformed > INT_MAX || pmc_xy > INT_MAX ||
        pmc_xy * fftz > INT_MAX)
    {
        Set_PME_Validation_Error(
            error,
            "the FFT grid product exceeds the supported INT_MAX "
            "element count");
        return false;
    }

    PME_Grid_Shape checked;
    checked.all = static_cast<int>(all);
    checked.input_plane = static_cast<int>(yz);
    checked.transformed = static_cast<int>(transformed);
    checked.pmc_transformed = static_cast<int>(pmc_xy * fftz);
    if (!PME_Allocation_Count_Fits(checked.all, sizeof(float)) ||
        !PME_Allocation_Count_Fits(checked.transformed, sizeof(FFT_COMPLEX)) ||
        !PME_Allocation_Count_Fits(checked.transformed, sizeof(LTMatrix3)) ||
        !PME_Allocation_Count_Fits(checked.pmc_transformed,
                                   sizeof(FFT_COMPLEX)))
    {
        Set_PME_Validation_Error(
            error, "an FFT allocation byte size exceeds SIZE_MAX");
        return false;
    }
    *shape = checked;
    return true;
}

// Solving erfc(beta * cutoff) == tolerance * cutoff in the dimensionless
// variable beta*cutoff avoids multiplying an ever-growing beta by cutoff.
// Both loops have hard bounds, and every result is checked before conversion
// back to the float used by the force kernels.
static bool Try_Compute_PME_Beta(double cutoff, double tolerance, float* beta,
                                 std::string* error)
{
    if (error != nullptr) error->clear();
    if (!Double_Memory_Is_Finite(&cutoff) || !(cutoff > 0.0))
    {
        Set_PME_Validation_Error(error, "cutoff must be finite and positive");
        return false;
    }
    if (!Double_Memory_Is_Finite(&tolerance) || !(tolerance > 0.0))
    {
        Set_PME_Validation_Error(
            error, "Direct_Tolerance must be finite and positive");
        return false;
    }
    const double target = tolerance * cutoff;
    if (!Double_Memory_Is_Finite(&target) || !(target > 0.0) || !(target < 1.0))
    {
        Set_PME_Validation_Error(
            error,
            "Direct_Tolerance * cutoff must be strictly between 0 "
            "and 1 so that a positive Ewald beta exists");
        return false;
    }

    double low = 0.0;
    double high = 1.0;
    bool bracketed = false;
    for (int iteration = 0; iteration < 64; ++iteration)
    {
        const double value = std::erfc(high);
        if (!Double_Memory_Is_Finite(&value))
        {
            Set_PME_Validation_Error(
                error,
                "erfc returned a non-finite value while bracketing "
                "the Ewald beta");
            return false;
        }
        if (value <= target)
        {
            bracketed = true;
            break;
        }
        high *= 2.0;
        if (!Double_Memory_Is_Finite(&high)) break;
    }
    if (!bracketed)
    {
        Set_PME_Validation_Error(
            error, "failed to bracket the Ewald beta within 64 expansions");
        return false;
    }

    for (int iteration = 0; iteration < 96; ++iteration)
    {
        const double middle = low + 0.5 * (high - low);
        if (std::erfc(middle) > target)
            low = middle;
        else
            high = middle;
    }
    const double solved = (low + 0.5 * (high - low)) / cutoff;
    if (!Double_Memory_Is_Finite(&solved) || !(solved > 0.0) ||
        solved > FLT_MAX)
    {
        Set_PME_Validation_Error(
            error, "the solved Ewald beta is not a finite positive float");
        return false;
    }
    const float checked = static_cast<float>(solved);
    if (!Float_Memory_Is_Normal(&checked) || !(checked > 0.0f))
    {
        Set_PME_Validation_Error(
            error, "the solved Ewald beta is not a positive normal float");
        return false;
    }
    *beta = checked;
    return true;
}

// ene += factor * charge_sum^2
static __global__ void device_add(float* ene, float factor, float* charge_sum)
{
    ene[0] += factor * charge_sum[0] * charge_sum[0];
}

static __global__ void charge_square_kernel(int element_number,
                                            const float* charge,
                                            float* charge_square)
{
    SIMPLE_DEVICE_FOR(i, element_number)
    {
        float q = charge[i];
        charge_square[i] = q * q;
    }
}

//--------Particle Mesh Ewald Method----------

void Particle_Mesh::Initial(CONTROLLER* controller, int atom_numbers,
                            LTMatrix3 cell, LTMatrix3 rcell, VECTOR box_length,
                            float cutoff,
                            int no_direct_interaction_virtual_atom_numbers,
                            const char* module_name)
{
    if (module_name == NULL)
    {
        strcpy(this->module_name, "PM");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }

    controller->printf("START INITIALIZING PME:\n");
    this->cutoff = cutoff;

    controller->printf("    PME backend library: %s\n", FFT_LIBRARY_NAME);

    std::string validation_error;
    int parsed_integer = 0;

    if (!Float_Memory_Is_Normal(&cutoff) || !(cutoff > 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tPME cutoff must be a finite positive normal float; "
            "got %.9g\n",
            cutoff);
        return;
    }

    tolerance = 0.00001f;
    if (controller->Command_Exist(this->module_name, "Direct_Tolerance"))
    {
        const char* token =
            controller->Command(this->module_name, "Direct_Tolerance");
        if (!Try_Parse_PME_Positive_Normal_Float(token, &tolerance,
                                                 &validation_error))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tinvalid %s.Direct_Tolerance value '%s': %s; the "
                "value must be a finite positive normal float\n",
                this->module_name, token, validation_error.c_str());
            return;
        }
    }

    if (CONTROLLER::PP_MPI_size == 1)
    {
        exclude_factor = 1.0f;
    }
    else
    {
        exclude_factor = 0.5f;
    }

    fftx = -1;
    ffty = -1;
    fftz = -1;
    if (controller->Command_Exist(this->module_name, "fftx"))
    {
        const char* token = controller->Command(this->module_name, "fftx");
        if (!Try_Parse_PME_Integer(token, &parsed_integer, &validation_error) ||
            parsed_integer < 4)
        {
            if (validation_error.empty())
                validation_error = "the dimension is smaller than 4";
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tinvalid %s.fftx value '%s': %s; fourth-order "
                "PME requires an integer dimension of at least 4\n",
                this->module_name, token, validation_error.c_str());
            return;
        }
        fftx = parsed_integer;
    }
    if (controller->Command_Exist(this->module_name, "ffty"))
    {
        const char* token = controller->Command(this->module_name, "ffty");
        if (!Try_Parse_PME_Integer(token, &parsed_integer, &validation_error) ||
            parsed_integer < 4)
        {
            if (validation_error.empty())
                validation_error = "the dimension is smaller than 4";
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tinvalid %s.ffty value '%s': %s; fourth-order "
                "PME requires an integer dimension of at least 4\n",
                this->module_name, token, validation_error.c_str());
            return;
        }
        ffty = parsed_integer;
    }
    if (controller->Command_Exist(this->module_name, "fftz"))
    {
        const char* token = controller->Command(this->module_name, "fftz");
        if (!Try_Parse_PME_Integer(token, &parsed_integer, &validation_error) ||
            parsed_integer < 4)
        {
            if (validation_error.empty())
                validation_error = "the dimension is smaller than 4";
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tinvalid %s.fftz value '%s': %s; fourth-order "
                "PME requires an integer dimension of at least 4\n",
                this->module_name, token, validation_error.c_str());
            return;
        }
        fftz = parsed_integer;
    }

    PM_MPI_size = 0;
    if (controller->Command_Exist(this->module_name, "MPI_size"))
    {
        const char* token = controller->Command(this->module_name, "MPI_size");
        if (!Try_Parse_PME_Integer(token, &PM_MPI_size, &validation_error))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tinvalid %s.MPI_size value '%s': %s\n",
                this->module_name, token, validation_error.c_str());
            return;
        }
    }
    else
    {
        PM_MPI_size = controller->PM_MPI_size;
    }
    if (!Try_Validate_PME_Process_Count(controller->MPI_size,
                                        controller->CC_MPI_size, PM_MPI_size,
                                        &validation_error) ||
        controller->MPI_rank < 0 ||
        controller->MPI_rank >= controller->MPI_size)
    {
        if (controller->MPI_rank < 0 ||
            controller->MPI_rank >= controller->MPI_size)
            validation_error = "MPI_rank is outside [0, MPI_size)";
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tinvalid PME process layout (MPI_size=%d, MPI_rank=%d, "
            "CC_MPI_size=%d, PM.MPI_size=%d): %s\n",
            controller->MPI_size, controller->MPI_rank, controller->CC_MPI_size,
            PM_MPI_size, validation_error.c_str());
        return;
    }

    if (!Try_Validate_PME_Implemented_Process_Roles(controller->CC_MPI_size,
                                                    &validation_error))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorNotImplemented, "Particle_Mesh::Initial",
            "Reason:\n\t%s, but got %d. Refusing this layout prevents CC "
            "ranks from being misclassified as PM ranks.\n",
            validation_error.c_str(), controller->CC_MPI_size);
        return;
    }

    // 2025-10-14: temporary disable multi-process PME
    if (PM_MPI_size > 1)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorNotImplemented, "Particle_Mesh::Initial",
            "Reason:\n\tMulti-process PME is not supported yet; use "
            "PM.MPI_size = 0 or 1.\n");
        return;
    }

    pm_pp_corres.assign(static_cast<std::size_t>(PM_MPI_size), {});
    pm_pp_num.assign(static_cast<std::size_t>(PM_MPI_size), 0);
    min_corner_set.assign(static_cast<std::size_t>(PM_MPI_size), VECTOR{});
    max_corner_set.assign(static_cast<std::size_t>(PM_MPI_size), VECTOR{});
    pm_corres_pp_num = 0;
    pm_corres_pp_rank_set.clear();
    pm_corres_pp_atom_number.clear();
    pm_corres_pp_atom_number_prefix.clear();
    pp_corres_pm_rank = -1;
    reported_pp_atom_numbers = -1;
    atom_mapping_is_valid = false;
    pm_dom_dec_split_num = {0, 0, 0};
    neighbor_num.fill(0);
    for (std::vector<int>& neighbors : neighbor_dir) neighbors.clear();

    if (!PM_MPI_size)
    {
        controller->printf("    PM reciprocal ranks disabled\n");
    }

    if (atom_numbers <= 0 || atom_numbers > INT_MAX / 64 ||
        no_direct_interaction_virtual_atom_numbers < 0 ||
        no_direct_interaction_virtual_atom_numbers > INT_MAX - atom_numbers ||
        !PME_Allocation_Count_Fits(atom_numbers, 64 * sizeof(int)) ||
        !PME_Allocation_Count_Fits(
            atom_numbers + no_direct_interaction_virtual_atom_numbers,
            sizeof(VECTOR)))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "Particle_Mesh::Initial",
            "Reason:\n\tinvalid PME atom counts: atom_numbers=%d and "
            "no_direct_interaction_virtual_atom_numbers=%d; counts must be "
            "positive/nonnegative and all PME index and allocation products "
            "must fit\n",
            atom_numbers, no_direct_interaction_virtual_atom_numbers);
        return;
    }
    this->atom_numbers = atom_numbers;
    this->max_atom_numbers = atom_numbers;

    const float box_extents[3] = {box_length.x, box_length.y, box_length.z};
    for (int axis = 0; axis < 3; ++axis)
    {
        if (!Float_Memory_Is_Normal(&box_extents[axis]) ||
            !(box_extents[axis] > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tPME box extent %c must be a finite positive "
                "normal float; got %.9g\n",
                'x' + axis, box_extents[axis]);
            return;
        }
    }
    const double volume_double = static_cast<double>(cell.a11) * cell.a22 *
                                 static_cast<double>(cell.a33);
    if (!Double_Memory_Is_Finite(&volume_double) || !(volume_double > 0.0) ||
        volume_double > FLT_MAX)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tPME cell volume must be finite, positive, and "
            "representable as float; got %.17g\n",
            volume_double);
        return;
    }
    const float volume = static_cast<float>(volume_double);
    if (!Float_Memory_Is_Normal(&volume))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tPME cell volume is not a positive normal float\n");
        return;
    }

    float grid_spacing = 1.0f;
    if (controller->Command_Exist(this->module_name, "grid_spacing"))
    {
        const char* token =
            controller->Command(this->module_name, "grid_spacing");
        if (!Try_Parse_PME_Positive_Normal_Float(token, &grid_spacing,
                                                 &validation_error))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tinvalid %s.grid_spacing value '%s': %s; the "
                "value must be a finite positive normal float\n",
                this->module_name, token, validation_error.c_str());
            return;
        }
    }
    controller->printf("    grid_spacing: %f Angstrom\n", grid_spacing);
    int* fft_dimensions[3] = {&fftx, &ffty, &fftz};
    for (int axis = 0; axis < 3; ++axis)
    {
        if (*fft_dimensions[axis] >= 0) continue;
        const double mesh_length =
            static_cast<double>(box_extents[axis]) / grid_spacing;
        if (!Try_Get_PME_Fft_Parameter(mesh_length, fft_dimensions[axis],
                                       &validation_error))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tcannot derive fft%c from box extent %.9g and "
                "grid_spacing %.9g: %s\n",
                'x' + axis, box_extents[axis], grid_spacing,
                validation_error.c_str());
            return;
        }
    }

    controller->printf("    fftx: %d\n", fftx);
    controller->printf("    ffty: %d\n", ffty);
    controller->printf("    fftz: %d\n", fftz);

    PME_Grid_Shape grid_shape;
    if (!Try_Build_PME_Grid_Shape(fftx, ffty, fftz, &grid_shape,
                                  &validation_error))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "Particle_Mesh::Initial",
            "Reason:\n\tinvalid PME FFT grid (%d, %d, %d): %s\n", fftx, ffty,
            fftz, validation_error.c_str());
        return;
    }
    PME_Nall = grid_shape.all;
    PME_Nin = grid_shape.input_plane;
    PME_Nfft = grid_shape.transformed;

    if (!Try_Compute_PME_Beta(cutoff, tolerance, &beta, &validation_error))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tcannot determine PME beta from cutoff %.9g and "
            "Direct_Tolerance %.9g: %s\n",
            cutoff, tolerance, validation_error.c_str());
        return;
    }
    controller->printf("    beta: %f\n", beta);

    const double neutralizing_factor_double =
        -0.5 * static_cast<double>(CONSTANT_Pi) /
        (static_cast<double>(beta) * beta * volume);
    if (!Double_Memory_Is_Finite(&neutralizing_factor_double) ||
        neutralizing_factor_double < -FLT_MAX ||
        neutralizing_factor_double > FLT_MAX)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tthe PME neutralizing factor is not a finite float\n");
        return;
    }
    neutralizing_factor = static_cast<float>(neutralizing_factor_double);
    if (!Float_Memory_Is_Normal(&neutralizing_factor))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tthe PME neutralizing factor underflows float\n");
        return;
    }

    int i, kx, ky, kz, index;
    FFT_RESULT errP1, errP2;
    update_interval = 1;
    if (controller->Command_Exist("PME", "update_interval"))
    {
        const char* token = controller->Command("PME", "update_interval");
        if (!Try_Parse_PME_Integer(token, &update_interval,
                                   &validation_error) ||
            update_interval <= 0)
        {
            if (validation_error.empty())
                validation_error = "the value is not positive";
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
                "Reason:\n\tinvalid PME.update_interval value '%s': %s; the "
                "value must be a positive integer\n",
                token, validation_error.c_str());
            return;
        }
    }

    Device_Malloc_Safely((void**)&num_ghost_dir_id,
                         sizeof(int) * max_atom_numbers * 6);
    Device_Malloc_Safely((void**)&charge_sum, sizeof(float));
    Device_Malloc_Safely((void**)&charge_square, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&force_backup, sizeof(VECTOR) * atom_numbers);
    deviceMemset(force_backup, 0, sizeof(VECTOR) * atom_numbers);
    Device_Malloc_Safely((void**)&PME_uxyz,
                         sizeof(UNSIGNED_INT_VECTOR) * atom_numbers);
    Device_Malloc_Safely((void**)&PME_frxyz, sizeof(VECTOR) * atom_numbers);
    Reset_List((int*)PME_uxyz, 1 << 30, 3 * atom_numbers);

    Device_Malloc_Safely((void**)&PME_Q, sizeof(float) * PME_Nall);
    Device_Malloc_Safely((void**)&PME_FQ, sizeof(FFT_COMPLEX) * PME_Nfft);
    Device_Malloc_Safely((void**)&PME_FBCFQ, sizeof(float) * PME_Nall);

    Device_Malloc_Safely((void**)&PME_atom_near,
                         sizeof(int) * 64 * atom_numbers);
    deviceMemset(PME_atom_near, 0, sizeof(int) * 64 * atom_numbers);

    FFT_SIZE_t n3d[3] = {fftx, ffty, fftz};
    errP1 =
        SPONGE_FFT_WRAPPER::Make_FFT_Plan(&PME_plan_r2c, 1, 3, n3d, FFT_R2C);
    errP2 =
        SPONGE_FFT_WRAPPER::Make_FFT_Plan(&PME_plan_c2r, 1, 3, n3d, FFT_C2R);
    if (errP1 != FFT_SUCCESS || errP2 != FFT_SUCCESS)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Initial",
            "Reason:\n\tError occurs when create fft plan of PME");
    }

#ifdef USE_CUDA
    // 单进程路径启用倒易链双流重叠：FFT plan 绑定到独立的非阻塞流，
    // 多进程（PM/PP 分离）路径不绑流，FFT 仍随默认流串行执行
    if (controller->MPI_size == 1)
    {
        deviceStreamCreateWithFlags(&pme_stream, deviceStreamNonBlocking);
        deviceEventCreateWithFlags(&pme_event_start, deviceEventDisableTiming);
        deviceEventCreateWithFlags(&pme_event_done, deviceEventDisableTiming);
        SPONGE_FFT_WRAPPER::Set_FFT_Stream(PME_plan_r2c, pme_stream);
        SPONGE_FFT_WRAPPER::Set_FFT_Stream(PME_plan_c2r, pme_stream);
        pme_overlap_enabled = true;
    }
#endif

    Device_Malloc_And_Copy_Safely((void**)&d_reciprocal_ene, &reciprocal_ene,
                                  sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_self_ene, &self_ene,
                                  sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_direct_ene, &direct_ene,
                                  sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_correction_ene, &correction_ene,
                                  sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_ee_ene, &ee_ene, sizeof(float));
    Device_Malloc_Safely((void**)&d_direct_atom_energy,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_correction_atom_energy,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&atom_id_l_g, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&atom_id_g_l, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&atom_id_g_l_candidate,
                         sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&atom_id_validation_error, 3 * sizeof(int));
    Device_Malloc_Safely(
        (void**)&g_crd,
        sizeof(VECTOR) *
            (atom_numbers + no_direct_interaction_virtual_atom_numbers));
    Device_Malloc_Safely(
        (void**)&g_frc,
        sizeof(VECTOR) *
            (atom_numbers + no_direct_interaction_virtual_atom_numbers));
    deviceMemset(atom_id_l_g, 0, sizeof(int) * atom_numbers);
    deviceMemset(atom_id_g_l, 0, sizeof(int) * atom_numbers);
    deviceMemset(atom_id_g_l_candidate, -1, sizeof(int) * atom_numbers);
    deviceMemset(atom_id_validation_error, 0, 3 * sizeof(int));
    deviceMemset(g_crd, 0,
                 sizeof(VECTOR) * (atom_numbers +
                                   no_direct_interaction_virtual_atom_numbers));
    deviceMemset(g_frc, 0,
                 sizeof(VECTOR) * (atom_numbers +
                                   no_direct_interaction_virtual_atom_numbers));
    deviceMemset(d_direct_atom_energy, 0, sizeof(float) * atom_numbers);
    deviceMemset(d_correction_atom_energy, 0, sizeof(float) * atom_numbers);

    calculate_reciprocal_part = true;
    if (controller->Command_Exist("PME", "calculate_reciprocal_part"))
    {
        calculate_reciprocal_part = controller->Get_Bool(
            "PME", "calculate_reciprocal_part", "Particle_Mesh::Initial");
    }
    calculate_excluded_part = true;
    if (controller->Command_Exist("PME", "calculate_excluded_part"))
    {
        calculate_excluded_part = controller->Get_Bool(
            "PME", "calculate_excluded_part", "Particle_Mesh::Initial");
    }
    bool use_pmc_iz = false;
    if (controller->Command_Exist("PME", "replaced_by_PMC_IZ"))
    {
        use_pmc_iz = controller->Get_Bool("PME", "replaced_by_PMC_IZ",
                                          "Particle_Mesh::Initial");
    }

    // 计算B-Spline修正系数 * 泊松算子因子， 用于倒空间乘法
    if (calculate_reciprocal_part)
    {
        if (use_pmc_iz)
        {
            controller->printf("    PMC-IZ will be used instead of PME\n");
            if (controller->Command_Choice("mode", "npt"))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorConflictingCommand, "Particle_Mesh::Initial",
                    "Reason:\n\tPMC-IZ can not be used in NPT mode");
            }
            Build_PMC_IZ_BC(
                controller, fftx, ffty, fftz, PME_Nfft, PME_Nall, PME_Nin,
                1.0f / box_length.x / box_length.x,
                1.0f / box_length.y / box_length.y, box_length.z / fftz, beta,
                CONSTANT_Pi / PME_Nall / box_length.x / box_length.y, &PME_BC);
        }
        else
        {
            float *B1 = NULL, *B2 = NULL, *B3 = NULL, *h_PME_BC = NULL,
                  *h_PME_BC0 = NULL;
            LTMatrix3* h_PME_virial_BC = NULL;
            B1 = (float*)malloc(sizeof(float) * fftx);
            ;
            B2 = (float*)malloc(sizeof(float) * ffty);
            B3 = (float*)malloc(sizeof(float) * fftz);
            h_PME_BC0 = (float*)malloc(sizeof(float) * PME_Nfft);
            h_PME_BC = (float*)malloc(sizeof(float) * PME_Nfft);
            h_PME_virial_BC = (LTMatrix3*)malloc(sizeof(LTMatrix3) * PME_Nfft);
            if (B1 == NULL || B2 == NULL || B3 == NULL || h_PME_BC0 == NULL ||
                h_PME_BC == NULL)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorMallocFailed, "Particle_Mesh::Initial",
                    "Reason:\n\tError occurs when malloc PME_BC of PME");
            }
            for (kx = 0; kx < fftx; kx++)
            {
                B1[kx] = getb(kx, fftx, 4);
            }

            for (ky = 0; ky < ffty; ky++)
            {
                B2[ky] = getb(ky, ffty, 4);
            }

            for (kz = 0; kz < fftz; kz++)
            {
                B3[kz] = getb(kz, fftz, 4);
            }

            float kxrp, kyrp, kzrp;
            float mprefactor = PI * PI / beta / beta;
            float msq;
            VECTOR m;
            for (kx = 0; kx < fftx; kx++)
            {
                kxrp = kx;
                if (kx > fftx / 2) kxrp = kx - fftx;
                for (ky = 0; ky < ffty; ky++)
                {
                    kyrp = ky;
                    if (ky > ffty / 2) kyrp = ky - ffty;
                    for (kz = 0; kz <= fftz / 2; kz++)
                    {
                        kzrp = kz;
                        m = {kxrp, kyrp, kzrp};
                        m = MultiplyTranspose(m, rcell);
                        msq = m * m;

                        index = kx * ffty * (fftz / 2 + 1) +
                                ky * (fftz / 2 + 1) + kz;

                        if (kx + ky + kz == 0)
                        {
                            h_PME_BC[index] = 0;
                            h_PME_virial_BC[index] = {0, 0, 0, 0, 0, 0};
                        }
                        else
                        {
                            h_PME_BC[index] = (float)1.0 / PI / msq *
                                              exp(-mprefactor * msq) / volume;
                            h_PME_virial_BC[index].a11 =
                                1 -
                                2 / msq * (1 + mprefactor * msq) * m.x * m.x;
                            h_PME_virial_BC[index].a21 =
                                0 -
                                2 / msq * (1 + mprefactor * msq) * m.y * m.x;
                            h_PME_virial_BC[index].a22 =
                                1 -
                                2 / msq * (1 + mprefactor * msq) * m.y * m.y;
                            h_PME_virial_BC[index].a31 =
                                0 -
                                2 / msq * (1 + mprefactor * msq) * m.z * m.x;
                            h_PME_virial_BC[index].a32 =
                                0 -
                                2 / msq * (1 + mprefactor * msq) * m.z * m.y;
                            h_PME_virial_BC[index].a33 =
                                1 -
                                2 / msq * (1 + mprefactor * msq) * m.z * m.z;
                        }
                        h_PME_BC0[index] = B1[kx] * B2[ky] * B3[kz];
                        h_PME_BC[index] *= h_PME_BC0[index];
                        h_PME_virial_BC[index] =
                            0.5f * h_PME_BC[index] * h_PME_virial_BC[index];
                    }
                }
            }

            Device_Malloc_Safely((void**)&PME_BC, sizeof(float) * PME_Nfft);
            Device_Malloc_Safely((void**)&PME_BC0, sizeof(float) * PME_Nfft);
            Device_Malloc_Safely((void**)&PME_Virial_BC,
                                 sizeof(LTMatrix3) * PME_Nfft);
            deviceMemcpy(PME_BC, h_PME_BC, sizeof(float) * PME_Nfft,
                         deviceMemcpyHostToDevice);
            deviceMemcpy(PME_BC0, h_PME_BC0, sizeof(float) * PME_Nfft,
                         deviceMemcpyHostToDevice);
            deviceMemcpy(PME_Virial_BC, h_PME_virial_BC,
                         sizeof(LTMatrix3) * PME_Nfft,
                         deviceMemcpyHostToDevice);
            free(B1);
            free(B2);
            free(B3);
            free(h_PME_BC0);
            free(h_PME_BC);
            free(h_PME_virial_BC);
        }
    }
    is_initialized = 1;
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial(this->module_name, "%.2f");
        if (controller->Command_Exist(this->module_name, "print_detail"))
        {
            print_detail = controller->Get_Bool(
                this->module_name, "print_detail", "Particle_Mesh::Initial");
            if (print_detail)
            {
                controller->Step_Print_Initial("PM_direct", "%.2f");
                controller->Step_Print_Initial("PM_reciprocal", "%.2f");
                controller->Step_Print_Initial("PM_self", "%.2f");
                controller->Step_Print_Initial("PM_correction", "%.2f");
            }
        }
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
    }
    controller->printf("END INITIALIZING PME\n\n");
}

void Particle_Mesh::Clear()
{
    if (is_initialized)
    {
        is_initialized = 0;
        Free_Single_Device_Pointer((void**)&PME_uxyz);
        Free_Single_Device_Pointer((void**)&PME_frxyz);
        Free_Single_Device_Pointer((void**)&PME_Q);
        Free_Single_Device_Pointer((void**)&PME_FQ);
        Free_Single_Device_Pointer((void**)&PME_FBCFQ);
        Free_Single_Device_Pointer((void**)&PME_BC);
        Free_Single_Device_Pointer((void**)&PME_Virial_BC);
        Free_Single_Device_Pointer((void**)&PME_BC0);
        Free_Single_Device_Pointer((void**)&charge_sum);
        Free_Single_Device_Pointer((void**)&charge_square);
        Free_Single_Device_Pointer((void**)&num_ghost_dir_id);

        Free_Single_Device_Pointer((void**)&atom_id_l_g);
        Free_Single_Device_Pointer((void**)&atom_id_g_l);
        Free_Single_Device_Pointer((void**)&atom_id_g_l_candidate);
        Free_Single_Device_Pointer((void**)&atom_id_validation_error);
        Free_Single_Device_Pointer((void**)&g_crd);
        Free_Single_Device_Pointer((void**)&g_frc);

        // Free_Single_Device_Pointer((void**)&MPI_PME_Q);
        // Free_Single_Device_Pointer((void**)&MPI_PME_FQ);
        // Free_Single_Device_Pointer((void**)&MPI_PME_FBCFQ);

        Free_Single_Device_Pointer((void**)&PME_atom_near);
        Free_Single_Device_Pointer((void**)&force_backup);

        SPONGE_FFT_WRAPPER::Destroy_FFT_Plan(&PME_plan_r2c);
        SPONGE_FFT_WRAPPER::Destroy_FFT_Plan(&PME_plan_c2r);

#ifdef USE_CUDA
        if (pme_overlap_enabled)
        {
            deviceEventDestroy(pme_event_start);
            deviceEventDestroy(pme_event_done);
            deviceStreamDestroy(pme_stream);
            pme_overlap_enabled = false;
        }
#endif
        pme_async_chain_active = false;
        pme_async_scheduled = false;

        Free_Host_And_Device_Pointer(NULL, (void**)&d_reciprocal_ene);
        Free_Host_And_Device_Pointer(NULL, (void**)&d_self_ene);
        Free_Host_And_Device_Pointer(NULL, (void**)&d_direct_ene);
        Free_Host_And_Device_Pointer(NULL, (void**)&d_correction_ene);
        Free_Host_And_Device_Pointer(NULL, (void**)&d_ee_ene);
        Free_Single_Device_Pointer((void**)&d_direct_atom_energy);
        Free_Single_Device_Pointer((void**)&d_correction_atom_energy);
    }
    pm_pp_corres.clear();
    pm_pp_num.clear();
    pm_corres_pp_rank_set.clear();
    pm_corres_pp_atom_number.clear();
    pm_corres_pp_atom_number_prefix.clear();
    min_corner_set.clear();
    max_corner_set.clear();
    for (std::vector<int>& neighbors : neighbor_dir) neighbors.clear();
    neighbor_num.fill(0);
    pm_corres_pp_num = 0;
    pp_corres_pm_rank = -1;
    reported_pp_atom_numbers = -1;
    atom_mapping_is_valid = false;
    pm_dom_dec_split_num = {0, 0, 0};
}

// 计算每个原子所在的网格点以及其周围64个网格点的索引
__global__ void PME_Atom_Near(const VECTOR* crd, int* PME_atom_near,
                              const int PME_Nin, const LTMatrix3 cell,
                              const LTMatrix3 rcell, const int atom_numbers,
                              const int fftx, const int ffty, const int fftz,
                              UNSIGNED_INT_VECTOR* PME_uxyz, VECTOR* PME_frxyz,
                              VECTOR* force_backup)
{
    SIMPLE_DEVICE_FOR(atom, atom_numbers)
    {
        force_backup[atom] = {0.0f, 0.0f, 0.0f};
        UNSIGNED_INT_VECTOR* temp_uxyz = &PME_uxyz[atom];
        VECTOR frac_crd = crd[atom] * rcell;
        frac_crd = frac_crd - floorf(frac_crd);
        if (!isfinite(frac_crd.x) || !isfinite(frac_crd.y) ||
            !isfinite(frac_crd.z))
        {
            frac_crd = {0.0f, 0.0f, 0.0f};
        }
        int k, tempux, tempuy, tempuz;
        frac_crd.x *= fftx;
        tempux = (int)frac_crd.x;
        tempux = tempux < 0 ? 0 : (tempux < fftx ? tempux : fftx - 1);
        PME_frxyz[atom].x = frac_crd.x - tempux;
        PME_frxyz[atom].x = PME_frxyz[atom].x - floorf(PME_frxyz[atom].x);
        frac_crd.y *= ffty;
        tempuy = (int)frac_crd.y;
        tempuy = tempuy < 0 ? 0 : (tempuy < ffty ? tempuy : ffty - 1);
        PME_frxyz[atom].y = frac_crd.y - tempuy;
        PME_frxyz[atom].y = PME_frxyz[atom].y - floorf(PME_frxyz[atom].y);
        frac_crd.z *= fftz;
        tempuz = (int)frac_crd.z;
        tempuz = tempuz < 0 ? 0 : (tempuz < fftz ? tempuz : fftz - 1);
        PME_frxyz[atom].z = frac_crd.z - tempuz;
        PME_frxyz[atom].z = PME_frxyz[atom].z - floorf(PME_frxyz[atom].z);
        if (tempux != (*temp_uxyz).uint_x || tempuy != (*temp_uxyz).uint_y ||
            tempuz != (*temp_uxyz).uint_z)
        {
            (*temp_uxyz).uint_x = tempux;
            (*temp_uxyz).uint_y = tempuy;
            (*temp_uxyz).uint_z = tempuz;
            int* temp_near = PME_atom_near + atom * 64;
            int kx, ky, kz;
            for (k = 0; k < 64; k++)
            {
                kx = k / 16;
                ky = (k - 16 * kx) / 4;
                kz = k % 4;

                kx = tempux - kx;

                if (kx < 0) kx += fftx;
                if (kx >= fftx) kx -= fftx;
                ky = tempuy - ky;
                if (ky < 0) ky += ffty;
                if (ky >= ffty) ky -= ffty;
                kz = tempuz - kz;
                if (kz < 0) kz += fftz;
                if (kz >= fftz) kz -= fftz;
                temp_near[k] = kx * PME_Nin + ky * fftz + kz;
            }
        }
    }
}

// 将原子电荷分配到其周围的64个网格点上
__global__ void PME_Q_Spread(int* PME_atom_near, const float* charge,
                             const VECTOR* PME_frxyz, float* PME_Q,
                             const int atom_numbers, const int PME_Nall)
{
    SIMPLE_DEVICE_FOR(atom, atom_numbers)
    {
        int k;
        float tempf, tempQ, tempf2;
        int* temp_near = PME_atom_near + atom * 64;
        VECTOR temp_frxyz = PME_frxyz[atom];
        float tempcharge = charge[atom];

        unsigned int kx;
#ifdef USE_GPU
        for (k = threadIdx.y; k < 64; k = k + blockDim.y)
#else
        for (k = 0; k < 64; k++)
#endif
        {
            kx = k / 16;
            tempf = temp_frxyz.x;
            tempf2 = tempf * tempf;
            tempf = PME_Ma[kx] * tempf * tempf2 + PME_Mb[kx] * tempf2 +
                    PME_Mc[kx] * tempf + PME_Md[kx];

            tempQ = tempcharge * tempf;

            kx = (k - kx * 16) / 4;
            tempf = temp_frxyz.y;
            tempf2 = tempf * tempf;
            tempf = PME_Ma[kx] * tempf * tempf2 + PME_Mb[kx] * tempf2 +
                    PME_Mc[kx] * tempf + PME_Md[kx];

            tempQ = tempQ * tempf;

            kx = k % 4;
            tempf = temp_frxyz.z;
            tempf2 = tempf * tempf;
            tempf = PME_Ma[kx] * tempf * tempf2 + PME_Mb[kx] * tempf2 +
                    PME_Mc[kx] * tempf + PME_Md[kx];
            tempQ = tempQ * tempf;

            int near_index = temp_near[k];
            if ((unsigned int)near_index < (unsigned int)PME_Nall)
            {
                atomicAdd(&PME_Q[near_index], tempQ);
            }
        }
    }
}

// 对FFT后的电荷密度进行修正
__global__ void PME_BCFQ(FFT_COMPLEX* PME_FQ, float* PME_BC, int PME_Nfft)
{
    SIMPLE_DEVICE_FOR(index, PME_Nfft)
    {
        float tempf = PME_BC[index];
        FFT_COMPLEX tempc = PME_FQ[index];
        REAL(PME_FQ[index]) = REAL(tempc) * tempf;
        IMAGINARY(PME_FQ[index]) = IMAGINARY(tempc) * tempf;
    }
}

// 计算每个原子受力
static __global__ void PME_Final(int* PME_atom_near, const float* charge,
                                 const float* PME_Q, VECTOR* force,
                                 const VECTOR* PME_frxyz, const LTMatrix3 rcell,
                                 const int fftx, const int ffty, const int fftz,
                                 const int atom_numbers, const int PME_Nall)
{
#ifdef GPU_ARCH_NAME
    int atom = blockDim.y * blockIdx.x + threadIdx.y;
    if (atom < atom_numbers)
#else
#pragma omp parallel for
    for (int atom = 0; atom < atom_numbers; atom++)
#endif
    {
        int k, kx;
        float tempdx, tempdy, tempdz, tempx, tempy, tempz, tempdQf;
        VECTOR tempdQ;
        float tempf, tempf2;
        float temp_charge = charge[atom];
        int* temp_near = PME_atom_near + atom * 64;
        VECTOR temp_frxyz = PME_frxyz[atom];
        VECTOR tempnv = {0, 0, 0};
#ifdef USE_GPU
        for (k = threadIdx.x; k < 64; k = k + blockDim.x)
#else
        for (k = 0; k < 64; k++)
#endif
        {
            int near_index = temp_near[k];
            if ((unsigned int)near_index >= (unsigned int)PME_Nall)
            {
                continue;
            }
            tempdQf = -PME_Q[near_index] * temp_charge;

            kx = k / 16;
            tempf = temp_frxyz.x;
            tempf2 = tempf * tempf;
            tempx = PME_Ma[kx] * tempf * tempf2 + PME_Mb[kx] * tempf2 +
                    PME_Mc[kx] * tempf + PME_Md[kx];
            tempdx = PME_dMa[kx] * tempf2 + PME_dMb[kx] * tempf + PME_dMc[kx];

            kx = (k - kx * 16) / 4;
            tempf = temp_frxyz.y;
            tempf2 = tempf * tempf;
            tempy = PME_Ma[kx] * tempf * tempf2 + PME_Mb[kx] * tempf2 +
                    PME_Mc[kx] * tempf + PME_Md[kx];
            tempdy = PME_dMa[kx] * tempf2 + PME_dMb[kx] * tempf + PME_dMc[kx];

            kx = k % 4;
            tempf = temp_frxyz.z;
            tempf2 = tempf * tempf;
            tempz = PME_Ma[kx] * tempf * tempf2 + PME_Mb[kx] * tempf2 +
                    PME_Mc[kx] * tempf + PME_Md[kx];
            tempdz = PME_dMa[kx] * tempf2 + PME_dMb[kx] * tempf + PME_dMc[kx];

            tempdQ.x = tempdx * tempy * tempz * fftx;
            tempdQ.y = tempdy * tempx * tempz * ffty;
            tempdQ.z = tempdz * tempx * tempy * fftz;
            tempdQ = tempdQf * MultiplyTranspose(tempdQ, rcell);
            tempnv = tempnv + tempdQ;
        }
        Warp_Sum_To(force + atom, tempnv, 8);
    }
}

// sum += list1 * list2
__global__ void PME_Energy_Product(const int element_number, const float* list1,
                                   const float* list2, float* sum)
{
#ifdef USE_GPU
    if (threadIdx.x == 0)
    {
        sum[0] = 0.;
    }
    __syncthreads();
#else
    sum[0] = 0;
#endif
    float lin = 0.0;
#ifdef USE_GPU
    for (int i = threadIdx.x; i < element_number; i = i + blockDim.x)
#else
#pragma omp parallel for reduction(+ : lin)
    for (int i = 0; i < element_number; i++)
#endif
    {
        lin = lin + list1[i] * list2[i];
    }
    atomicAdd(sum, lin);
}

static __global__ void PME_Excluded_Force_With_Atom_Energy_Correction(
    const int atom_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const float* charge, const float pme_beta,
    const int* excluded_list_start, const int* excluded_list,
    const int* excluded_atom_numbers, VECTOR* frc, float* atom_ene,
    float* this_ene, LTMatrix3* atom_virial)
{
    SIMPLE_DEVICE_FOR(atom_i, atom_numbers)
    {
        int excluded_numbers = excluded_atom_numbers[atom_i];
        if (excluded_numbers > 0)
        {
            int list_start = excluded_list_start[atom_i];
            int list_end = list_start + excluded_numbers;
            int atom_j;

            float charge_i = charge[atom_i];
            float charge_j;
            VECTOR r1 = crd[atom_i], r2;
            VECTOR dr;
            float dr2;

            float frc_abs = 0.;
            VECTOR frc_lin;
            VECTOR frc_record = {0., 0., 0.};
            LTMatrix3 virial_record = {0, 0, 0, 0, 0, 0};
            float ene_lin = 0.;

            for (int i = list_start; i < list_end; i = i + 1)
            {
                atom_j = excluded_list[i];
                r2 = crd[atom_j];
                charge_j = charge[atom_j];
                if (!PairwiseInteraction::Coulomb_Is_Active(charge_i, charge_j))
                {
                    continue;
                }

                dr = Get_Periodic_Displacement(r2, r1, cell, rcell);
                dr2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;

                const PME_Excluded_Radial_Kernels radial =
                    Get_PME_Excluded_Radial_Kernels(dr2, pme_beta);
                const float charge_product = charge_i * charge_j;
                frc_abs = charge_product * radial.force;
                frc_lin = frc_abs * dr;
                ene_lin -= charge_product * radial.energy;
                frc_record = frc_record + frc_lin;
                atomicAdd(frc + atom_j, -frc_lin);
                virial_record =
                    virial_record - Get_Virial_From_Force_Dis(frc_lin, dr);
            }  // atom_j cycle
            atomicAdd(frc + atom_i, frc_record);
            atomicAdd(atom_virial + atom_i, virial_record);
            atomicAdd(atom_ene + atom_i, ene_lin);
            this_ene[atom_i] = ene_lin;
        }  // if need excluded
    }
}

void Particle_Mesh::PME_Excluded_Force_With_Atom_Energy(
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    const float* charge, const int* excluded_list_start,
    const int* excluded_list, const int* excluded_atom_numbers, VECTOR* frc,
    int need_energy, float* atom_ene, LTMatrix3* atom_virial)
{
    if (is_initialized && calculate_excluded_part)
    {
        if (need_energy)
            deviceMemset(d_correction_atom_energy, 0,
                         sizeof(float) * atom_numbers);
        if (CONTROLLER::MPI_rank != 0) return;
        Launch_Device_Kernel(
            PME_Excluded_Force_With_Atom_Energy_Correction,
            (atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers, crd, cell,
            rcell, charge, beta, excluded_list_start, excluded_list,
            excluded_atom_numbers, frc, atom_ene, d_correction_atom_energy,
            atom_virial);
    }
}

static __global__ void PME_Add_Energy_To_Potential(float* d_ene,
                                                   float* d_self_ene,
                                                   float* d_reciprocal_ene)
{
    d_ene[0] += d_self_ene[0] + d_reciprocal_ene[0];
}

static __global__ void device_add_force(const int atom_numbers,
                                        float update_interval, VECTOR* force,
                                        const VECTOR* force_backup)
{
    SIMPLE_DEVICE_FOR(atom_i, atom_numbers)
    {
        force[atom_i] = force[atom_i] + update_interval * force_backup[atom_i];
    }
}

static __global__ void PME_Sum_Virial(const int nfft,
                                      const LTMatrix3* virial_BC,
                                      const FFT_COMPLEX* FQ, LTMatrix3* virial,
                                      int fftz)
{
    LTMatrix3 vir = {0, 0, 0, 0, 0, 0};
#ifdef USE_GPU
    int tid = blockDim.x * blockIdx.x * blockDim.y + threadIdx.x * blockDim.y +
              threadIdx.y;
    for (int index = tid; index < nfft;
         index += blockDim.x * blockDim.y * gridDim.x)
    {
        int fftc = fftz / 2 + 1;
        int nz = index % fftc;
        float factor = (nz == 0 || nz == fftc - 1) ? 0.5f : 1.0f;
        FFT_COMPLEX FQ0 = FQ[index];
        LTMatrix3 vir0 =
            factor * (REAL(FQ0) * REAL(FQ0) + IMAGINARY(FQ0) * IMAGINARY(FQ0)) *
            virial_BC[index];
        vir = vir - vir0;
    }
#else
    float v11 = 0.0f, v21 = 0.0f, v22 = 0.0f;
    float v31 = 0.0f, v32 = 0.0f, v33 = 0.0f;
#pragma omp parallel for reduction(+ : v11, v21, v22, v31, v32, v33)
    for (int index = 0; index < nfft; index++)
    {
        int fftc = fftz / 2 + 1;
        int nz = index % fftc;
        float factor = (nz == 0 || nz == fftc - 1) ? 0.5f : 1.0f;
        FFT_COMPLEX FQ0 = FQ[index];
        LTMatrix3 vir0 =
            factor * (REAL(FQ0) * REAL(FQ0) + IMAGINARY(FQ0) * IMAGINARY(FQ0)) *
            virial_BC[index];
        v11 -= vir0.a11;
        v21 -= vir0.a21;
        v22 -= vir0.a22;
        v31 -= vir0.a31;
        v32 -= vir0.a32;
        v33 -= vir0.a33;
    }
    vir = {v11, v21, v22, v31, v32, v33};
#endif
    Warp_Sum_To(virial, vir, warpSize);
}

void Particle_Mesh::PME_Reciprocal_Chain(const VECTOR* crd,
                                         const LTMatrix3 cell,
                                         const LTMatrix3 rcell,
                                         const float* charge, int need_virial,
                                         LTMatrix3* d_virial,
                                         deviceStream_t stream)
{
    // 计算插值索引
    deviceMemsetAsync(PME_Q, 0, sizeof(float) * PME_Nall, stream);
    Launch_Device_Kernel(
        PME_Atom_Near,
        (atom_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, stream, crd, PME_atom_near, PME_Nin,
        cell, rcell, atom_numbers, fftx, ffty, fftz, PME_uxyz, PME_frxyz,
        force_backup);

    dim3 blockSize = {CONTROLLER::device_max_thread / 64, 64};

    // 电荷Bspline插值
    Launch_Device_Kernel(PME_Q_Spread,
                         (atom_numbers + blockSize.x - 1) / blockSize.x,
                         blockSize, 0, stream, PME_atom_near, charge,
                         PME_frxyz, PME_Q, atom_numbers, PME_Nall);

    // do FFT
    SPONGE_FFT_WRAPPER::R2C(PME_plan_r2c, PME_Q, PME_FQ);

    // 修正Bspline插值
    blockSize = {CONTROLLER::device_warp,
                 CONTROLLER::device_max_thread / CONTROLLER::device_warp};
    if (need_virial)
        Launch_Device_Kernel(
            PME_Sum_Virial,
            (PME_Nfft + 4 * CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            blockSize, 0, stream, PME_Nfft, PME_Virial_BC, PME_FQ, d_virial,
            fftz);

    Launch_Device_Kernel(PME_BCFQ,
                         (PME_Nfft + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, stream, PME_FQ,
                         PME_BC, PME_Nfft);

    // do inverse FFT
    SPONGE_FFT_WRAPPER::C2R(PME_plan_c2r, PME_FQ, PME_FBCFQ);

    // 计算势能和力
    blockSize = {8, CONTROLLER::device_max_thread / 8};
    // kernel 内 atom 索引为 blockDim.y * blockIdx.x + threadIdx.y，
    // grid 按 blockSize.y 划分
    Launch_Device_Kernel(PME_Final,
                         (atom_numbers + blockSize.y - 1) / blockSize.y,
                         blockSize, 0, stream, PME_atom_near, charge,
                         PME_FBCFQ, force_backup, PME_frxyz, rcell, fftx, ffty,
                         fftz, atom_numbers, PME_Nall);
}

void Particle_Mesh::PME_Reciprocal_Energy_Tail(const float* charge,
                                               float* d_potential)
{
    Launch_Device_Kernel(PME_Energy_Product, 1, CONTROLLER::device_max_thread,
                         0, NULL, PME_Nall, PME_Q, PME_FBCFQ,
                         d_reciprocal_ene);
    Scale_List(d_reciprocal_ene, 0.5f, 1);

    Launch_Device_Kernel(charge_square_kernel,
                         (atom_numbers + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
                         charge, charge_square);
    Sum_Of_List(charge_square, d_self_ene, atom_numbers);

    Scale_List(d_self_ene, -beta / sqrt(PI), 1);

    Sum_Of_List(charge, charge_sum, atom_numbers);

    Launch_Device_Kernel(device_add, 1, 1, 0, NULL, d_self_ene,
                         neutralizing_factor, charge_sum);

    Launch_Device_Kernel(PME_Add_Energy_To_Potential, 1, 1, 0, NULL,
                         d_potential, d_self_ene, d_reciprocal_ene);
}

void Particle_Mesh::PME_Reciprocal_Force_With_Energy_And_Virial(
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    const float* charge, VECTOR* force, int need_virial, int need_energy,
    LTMatrix3* d_virial, float* d_potential, int step, bool exact_state)
{
    if (is_initialized && calculate_reciprocal_part)
    {
        if (need_energy)
        {
            deviceMemset(d_reciprocal_ene, 0, sizeof(float));
            deviceMemset(d_self_ene, 0, sizeof(float));
        }
        const bool scheduled_force_update = step % update_interval == 0;
        // MTS force impulses remain on their physical schedule.  Energy,
        // virial, and transactional evaluations nevertheless require a mesh
        // built from the current coordinates and box; using the previous
        // mesh would compare two different states in an MC acceptance test.
        if (scheduled_force_update || need_energy || need_virial || exact_state)
        {
            PME_Reciprocal_Chain(crd, cell, rcell, charge, need_virial,
                                 d_virial, NULL);
            if (scheduled_force_update)
            {
                Launch_Device_Kernel(
                    device_add_force,
                    (atom_numbers + CONTROLLER::device_max_thread - 1) /
                        CONTROLLER::device_max_thread,
                    CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
                    update_interval, force, force_backup);
            }
        }
        if (need_energy)
        {
            PME_Reciprocal_Energy_Tail(charge, d_potential);
        }
    }
}

void Particle_Mesh::PME_Reciprocal_Force_Async_Start(
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    const float* charge, int need_virial, int need_energy, LTMatrix3* d_virial,
    int step, bool exact_state)
{
    pme_async_chain_active = false;
    pme_async_scheduled = false;
    if (!is_initialized || !calculate_reciprocal_part) return;
    const bool scheduled_force_update = step % update_interval == 0;
    // 触发条件与串行版本完全一致：MTS 调度步，或需要能量/维里/精确态时
    if (!(scheduled_force_update || need_energy || need_virial || exact_state))
        return;
    deviceStream_t stream = NULL;
#ifdef USE_CUDA
    if (pme_overlap_enabled)
    {
        // 在 legacy 默认流上记录 crd 定型事件，pme_stream 等它，保证倒易链
        // 读到的是当期坐标（含此前的 ghost 更新、邻居表更新等全部默认流工作）
        deviceEventRecord(pme_event_start, NULL);
        deviceStreamWaitEvent(pme_stream, pme_event_start, 0);
        stream = pme_stream;
    }
#endif
    PME_Reciprocal_Chain(crd, cell, rcell, charge, need_virial, d_virial,
                         stream);
#ifdef USE_CUDA
    if (pme_overlap_enabled)
    {
        deviceEventRecord(pme_event_done, pme_stream);
    }
#endif
    pme_async_chain_active = true;
    pme_async_scheduled = scheduled_force_update;
}

void Particle_Mesh::PME_Reciprocal_Force_Async_Join(const float* charge,
                                                    VECTOR* force,
                                                    int need_energy,
                                                    float* d_potential)
{
    if (!is_initialized || !calculate_reciprocal_part) return;
    if (pme_async_chain_active)
    {
#ifdef USE_CUDA
        if (pme_overlap_enabled)
        {
            // 默认流等待 pme_stream 上的倒易链结束，之后 device_add_force
            // 与能量尾部在默认流上按原串行顺序执行，保证逐位一致
            deviceStreamWaitEvent(NULL, pme_event_done, 0);
        }
#endif
        if (pme_async_scheduled)
        {
            Launch_Device_Kernel(
                device_add_force,
                (atom_numbers + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
                update_interval, force, force_backup);
        }
        pme_async_chain_active = false;
        pme_async_scheduled = false;
    }
    if (need_energy)
    {
        deviceMemset(d_reciprocal_ene, 0, sizeof(float));
        deviceMemset(d_self_ene, 0, sizeof(float));
        PME_Reciprocal_Energy_Tail(charge, d_potential);
    }
}

// 计算PME的位移势能和Virial张量
static __global__ void up_box_bc(int fftx, int ffty, int fftz, float* PME_BC,
                                 float* PME_BC0, LTMatrix3* PME_virial_BC,
                                 float mprefactor, LTMatrix3 rcell,
                                 float volume)
{
    float kxrp, kyrp, kzrp;
    int ky, kz, index;
    float msq;
    VECTOR m;
    LTMatrix3 virial_bc_local;
    float bc_local;
#ifdef USE_GPU
    for (int kx = blockIdx.x * blockDim.x + threadIdx.x; kx < fftx;
         kx += blockDim.x * gridDim.x)
#else
#pragma omp parallel for firstprivate(kxrp, kyrp, kzrp, ky, kz, index, msq, m, \
                                          virial_bc_local, bc_local)
    for (int kx = 0; kx < fftx; kx++)
#endif
    {
        kxrp = kx;
        if (kx > fftx / 2) kxrp = kx - fftx;
#ifdef USE_GPU
        for (ky = blockIdx.y * blockDim.y + threadIdx.y; ky < ffty;
             ky += blockDim.y * gridDim.y)
#else
        for (ky = 0; ky < ffty; ky++)
#endif
        {
            kyrp = ky;
            if (ky > ffty / 2) kyrp = ky - ffty;
#ifdef USE_GPU
            for (kz = threadIdx.z; kz <= fftz / 2; kz += blockDim.z)
#else
            for (kz = 0; kz <= fftz / 2; kz++)
#endif
            {
                kzrp = kz;
                m = {kxrp, kyrp, kzrp};
                m = MultiplyTranspose(m, rcell);
                msq = m * m;

                index = kx * ffty * (fftz / 2 + 1) + ky * (fftz / 2 + 1) + kz;

                if (kx + ky + kz == 0)
                {
                    PME_BC[index] = 0;
                    PME_virial_BC[index] = {0, 0, 0, 0, 0, 0};
                }
                else
                {
                    bc_local = (float)1.0 / PI / msq * exp(mprefactor * msq) /
                               volume * PME_BC0[index];
                    virial_bc_local.a11 =
                        1 - 2 / msq * (1 + mprefactor * msq) * m.x * m.x;
                    virial_bc_local.a21 =
                        0 - 2 / msq * (1 + mprefactor * msq) * m.y * m.x;
                    virial_bc_local.a22 =
                        1 - 2 / msq * (1 + mprefactor * msq) * m.y * m.y;
                    virial_bc_local.a31 =
                        0 - 2 / msq * (1 + mprefactor * msq) * m.z * m.x;
                    virial_bc_local.a32 =
                        0 - 2 / msq * (1 + mprefactor * msq) * m.z * m.y;
                    virial_bc_local.a33 =
                        1 - 2 / msq * (1 + mprefactor * msq) * m.z * m.z;
                    PME_virial_BC[index] = 0.5f * bc_local * virial_bc_local;
                    PME_BC[index] = bc_local;
                }
            }
        }
    }
}

static void Scale_Positions_Device(const LTMatrix3 g, VECTOR* crd, float dt)
{
    VECTOR r_dash;
    r_dash.x = crd[0].x +
               dt * (crd[0].x * g.a11 + crd[0].y * g.a21 + crd[0].z * g.a31);
    r_dash.y = crd[0].y + dt * (crd[0].y * g.a22 + crd[0].z * g.a32);
    r_dash.z = crd[0].z + dt * crd[0].z * g.a33;
    crd[0] = r_dash;
}

void Particle_Mesh::Update_Box(LTMatrix3 cell, LTMatrix3 rcell, LTMatrix3 g,
                               float dt)
{
    float volume = cell.a11 * cell.a22 * cell.a33;
    neutralizing_factor = -0.5 * CONSTANT_Pi / (beta * beta * volume);
    float mprefactor = PI * PI / -beta / beta;
    dim3 blockSize = {8, 8, CONTROLLER::device_max_thread / 64};
    dim3 gridSize = {64, 64};
    Launch_Device_Kernel(up_box_bc, gridSize, blockSize, 0, NULL, fftx, ffty,
                         fftz, PME_BC, PME_BC0, PME_Virial_BC, mprefactor,
                         rcell, volume);
    Scale_Positions_Device(g, &min_corner, dt);
    Scale_Positions_Device(g, &max_corner, dt);
}

//-------domain-decomposition and communication----------------

struct PME_Domain_Layout
{
    INT_VECTOR pm_split = {0, 0, 0};
    std::vector<std::vector<int>> pm_to_pp;
    std::vector<VECTOR> minimum_corners;
    std::vector<VECTOR> maximum_corners;
};

// 找出n的所有因子，存入factor_set，并按从大到小排序
static std::vector<int> Find_PME_Factors(int n)
{
    std::vector<int> factors;
    for (int factor = 1; factor <= n / factor; ++factor)
    {
        if (n % factor != 0) continue;
        factors.push_back(factor);
        if (factor != n / factor) factors.push_back(n / factor);
    }
    std::sort(factors.begin(), factors.end(), std::greater<int>());
    return factors;
}

static bool Try_Build_PME_Domain_Layout(int pm_size, int pp_size,
                                        INT_VECTOR pp_split, VECTOR box_length,
                                        PME_Domain_Layout* layout,
                                        std::string* error)
{
    if (error != nullptr) error->clear();
    if (pm_size <= 0 || pp_size <= 0)
    {
        Set_PME_Validation_Error(error,
                                 "PM and PP rank counts must both be positive");
        return false;
    }
    if (pp_size % pm_size != 0)
    {
        Set_PME_Validation_Error(
            error,
            "the PP rank count must be divisible by the PM rank "
            "count");
        return false;
    }
    // Rank lists are sent as MPI_BYTE with an int byte count.
    if (pp_size > INT_MAX / static_cast<int>(sizeof(int)))
    {
        Set_PME_Validation_Error(
            error,
            "the PP rank list byte count exceeds MPI's int count "
            "range");
        return false;
    }
    if (pp_split.int_x <= 0 || pp_split.int_y <= 0 || pp_split.int_z <= 0)
    {
        Set_PME_Validation_Error(
            error, "all PP domain split dimensions must be positive");
        return false;
    }
    const long long pp_xy =
        static_cast<long long>(pp_split.int_x) * pp_split.int_y;
    if (pp_xy > pp_size || pp_xy * pp_split.int_z != pp_size)
    {
        Set_PME_Validation_Error(
            error,
            "the PP domain split product must equal the PP rank "
            "count without overflow");
        return false;
    }
    const float extents[3] = {box_length.x, box_length.y, box_length.z};
    for (float extent : extents)
    {
        if (!Float_Memory_Is_Normal(&extent) || !(extent > 0.0f))
        {
            Set_PME_Validation_Error(
                error,
                "all PME box extents must be finite positive normal "
                "floats");
            return false;
        }
    }

    const std::vector<int> factors_x = Find_PME_Factors(pp_split.int_x);
    const std::vector<int> factors_y = Find_PME_Factors(pp_split.int_y);
    INT_VECTOR pm_split = {0, 0, 0};
    bool found_split = false;
    for (int split_x : factors_x)
    {
        if (pm_size % split_x != 0) continue;
        const int after_x = pm_size / split_x;
        for (int split_y : factors_y)
        {
            if (after_x % split_y != 0) continue;
            const int split_z = after_x / split_y;
            if (pp_split.int_z % split_z != 0) continue;
            pm_split = {split_x, split_y, split_z};
            found_split = true;
            break;
        }
        if (found_split) break;
    }
    if (!found_split)
    {
        Set_PME_Validation_Error(
            error,
            "no PM split compatible with the PP domain split was "
            "found");
        return false;
    }
    const long long pm_product = static_cast<long long>(pm_split.int_x) *
                                 pm_split.int_y * pm_split.int_z;
    if (pm_product != pm_size)
    {
        Set_PME_Validation_Error(
            error,
            "the selected PM domain split product is inconsistent "
            "with the PM rank count");
        return false;
    }

    try
    {
        PME_Domain_Layout checked;
        checked.pm_split = pm_split;
        checked.pm_to_pp.resize(static_cast<std::size_t>(pm_size));
        checked.minimum_corners.resize(static_cast<std::size_t>(pm_size));
        checked.maximum_corners.resize(static_cast<std::size_t>(pm_size));
        const int pp_per_pm = pp_size / pm_size;
        for (std::vector<int>& ranks : checked.pm_to_pp)
            ranks.reserve(static_cast<std::size_t>(pp_per_pm));

        const int pm_xy = pm_split.int_x * pm_split.int_y;
        for (int k = 0; k < pm_split.int_z; ++k)
        {
            for (int j = 0; j < pm_split.int_y; ++j)
            {
                for (int i = 0; i < pm_split.int_x; ++i)
                {
                    const int rank_id = i + j * pm_split.int_x + k * pm_xy;
                    if (rank_id < 0 || rank_id >= pm_size)
                    {
                        Set_PME_Validation_Error(
                            error,
                            "a generated PM corner rank is out of "
                            "range");
                        return false;
                    }
                    checked.minimum_corners[rank_id] = {
                        box_length.x / pm_split.int_x * i,
                        box_length.y / pm_split.int_y * j,
                        box_length.z / pm_split.int_z * k};
                    checked.maximum_corners[rank_id] = {
                        box_length.x / pm_split.int_x * (i + 1),
                        box_length.y / pm_split.int_y * (j + 1),
                        box_length.z / pm_split.int_z * (k + 1)};
                }
            }
        }

        const int pp_per_pm_x = pp_split.int_x / pm_split.int_x;
        const int pp_per_pm_y = pp_split.int_y / pm_split.int_y;
        const int pp_per_pm_z = pp_split.int_z / pm_split.int_z;
        for (int k = 0; k < pp_split.int_z; ++k)
        {
            for (int j = 0; j < pp_split.int_y; ++j)
            {
                for (int i = 0; i < pp_split.int_x; ++i)
                {
                    const int pp_rank =
                        i + j * pp_split.int_x + k * static_cast<int>(pp_xy);
                    const int pm_rank =
                        i / pp_per_pm_x + j / pp_per_pm_y * pm_split.int_x +
                        k / pp_per_pm_z * pm_split.int_x * pm_split.int_y;
                    if (pp_rank < 0 || pp_rank >= pp_size || pm_rank < 0 ||
                        pm_rank >= pm_size)
                    {
                        Set_PME_Validation_Error(
                            error,
                            "a generated PM/PP correspondence rank is "
                            "out of range");
                        return false;
                    }
                    checked.pm_to_pp[pm_rank].push_back(pp_rank);
                }
            }
        }
        for (const std::vector<int>& ranks : checked.pm_to_pp)
        {
            if (ranks.size() != static_cast<std::size_t>(pp_per_pm))
            {
                Set_PME_Validation_Error(
                    error,
                    "the PM/PP correspondence is not uniformly "
                    "partitioned");
                return false;
            }
        }
        *layout = std::move(checked);
    }
    catch (const std::bad_alloc&)
    {
        Set_PME_Validation_Error(
            error,
            "host allocation failed while constructing the dynamic "
            "PM/PP rank layout");
        return false;
    }
    return true;
}

// Domain Decomposition
void Particle_Mesh::Domain_Decomposition(CONTROLLER* controller,
                                         VECTOR box_length,
                                         INT_VECTOR pp_split_num)
{
    if (!PM_MPI_size)
    {
        return;
    }

    std::string validation_error;
    const long long expected_pp_size =
        controller->MPI_size == 1
            ? 1
            : static_cast<long long>(controller->MPI_size) -
                  controller->CC_MPI_size - PM_MPI_size;
    if (!Try_Validate_PME_Process_Count(controller->MPI_size,
                                        controller->CC_MPI_size, PM_MPI_size,
                                        &validation_error) ||
        controller->PM_MPI_size != PM_MPI_size ||
        expected_pp_size != controller->PP_MPI_size ||
        controller->MPI_rank < 0 ||
        controller->MPI_rank >= controller->MPI_size)
    {
        if (validation_error.empty())
            validation_error =
                "controller PM/PP rank counts or MPI_rank are inconsistent";
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "Particle_Mesh_Ewald::Domain_Decomposition",
            "Reason:\n\tinvalid PME rank layout: MPI_size=%d, MPI_rank=%d, "
            "PP_MPI_size=%d, PM_MPI_size=%d, CC_MPI_size=%d, local "
            "PM_MPI_size=%d: %s\n",
            controller->MPI_size, controller->MPI_rank, controller->PP_MPI_size,
            controller->PM_MPI_size, controller->CC_MPI_size, PM_MPI_size,
            validation_error.c_str());
        return;
    }
    // 区域分割只在主进程上做一次。
    if (controller->MPI_rank != 0) return;

    PME_Domain_Layout layout;
    if (!Try_Build_PME_Domain_Layout(PM_MPI_size, controller->PP_MPI_size,
                                     pp_split_num, box_length, &layout,
                                     &validation_error))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "Particle_Mesh_Ewald::Domain_Decomposition",
            "Reason:\n\tcannot construct the PME domain layout for %d PM "
            "ranks, %d PP ranks, and PP split (%d, %d, %d): %s\n",
            PM_MPI_size, controller->PP_MPI_size, pp_split_num.int_x,
            pp_split_num.int_y, pp_split_num.int_z, validation_error.c_str());
        return;
    }

    pm_dom_dec_split_num = layout.pm_split;
    pm_pp_corres = std::move(layout.pm_to_pp);
    min_corner_set = std::move(layout.minimum_corners);
    max_corner_set = std::move(layout.maximum_corners);
    pm_pp_num.resize(pm_pp_corres.size());
    for (std::size_t pm = 0; pm < pm_pp_corres.size(); ++pm)
        pm_pp_num[pm] = static_cast<int>(pm_pp_corres[pm].size());
    controller->printf("    PM domain split: (%d, %d, %d)\n",
                       pm_dom_dec_split_num.int_x, pm_dom_dec_split_num.int_y,
                       pm_dom_dec_split_num.int_z);
}

void Particle_Mesh::Send_Recv_Dom_Dec(CONTROLLER* controller)
{
    if (!PM_MPI_size) return;

    const bool root_storage_is_valid =
        controller->MPI_rank != 0 ||
        (pm_pp_corres.size() == static_cast<std::size_t>(PM_MPI_size) &&
         pm_pp_num.size() == static_cast<std::size_t>(PM_MPI_size) &&
         min_corner_set.size() == static_cast<std::size_t>(PM_MPI_size) &&
         max_corner_set.size() == static_cast<std::size_t>(PM_MPI_size));
    if (controller->PM_MPI_size != PM_MPI_size ||
        controller->PP_MPI_size <= 0 || !root_storage_is_valid)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Send_Recv_Dom_Dec",
            "Reason:\n\tPME rank counts or dynamic domain-layout storage are "
            "inconsistent (PP=%d, controller PM=%d, local PM=%d)\n",
            controller->PP_MPI_size, controller->PM_MPI_size, PM_MPI_size);
        return;
    }
    // 如果PM与PP共用一个进程，则不需要通信
    if (controller->MPI_size == 1 && PM_MPI_size == 1)
    {
        strcpy(this->FFT_MPI_TYPE, "DISABLE");
        pp_corres_pm_rank = 0;
        pm_corres_pp_num = 1;
        pm_corres_pp_rank_set.assign(1, 0);
        pm_corres_pp_atom_number.assign(1, 0);
        pm_corres_pp_atom_number_prefix.assign(1, 0);
        min_corner = min_corner_set[0];
        max_corner = max_corner_set[0];
        return;
    }
#ifdef USE_MPI
    // PP进程与PM进程分割
    if (controller->MPI_rank == 0)
    {
        // 发送PP进程对应的PM进程号
        for (int pm_id = 0; pm_id < controller->PM_MPI_size; ++pm_id)
        {
            const int expected_count =
                static_cast<int>(pm_pp_corres[pm_id].size());
            if (pm_pp_num[pm_id] != expected_count || expected_count <= 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "Particle_Mesh::Send_Recv_Dom_Dec",
                    "Reason:\n\tPM rank %d has an inconsistent PP rank list "
                    "size (%d versus %d)\n",
                    pm_id, pm_pp_num[pm_id], expected_count);
                return;
            }
            const int pm_rank_tot = pm_id + controller->PP_MPI_size;
            for (int pp_rank : pm_pp_corres[pm_id])
            {
                if (pp_rank < 0 || pp_rank >= controller->PP_MPI_size)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorValueErrorCommand,
                        "Particle_Mesh::Send_Recv_Dom_Dec",
                        "Reason:\n\tPM rank %d contains out-of-range PP rank "
                        "%d\n",
                        pm_id, pp_rank);
                    return;
                }
                if (pp_rank != 0)
                {
                    MPI_Send(&pm_rank_tot, sizeof(int), MPI_BYTE, pp_rank,
                             PME_MPI_TAG_DOMAIN_PM_ASSIGNMENT, MPI_COMM_WORLD);
                }
                else
                {
                    pp_corres_pm_rank = pm_rank_tot;
                }
            }
        }
        // 发送PM进程对应的pp进程数与进程索引集合；发送域分割信息
        for (int pm_id = 0; pm_id < controller->PM_MPI_size; ++pm_id)
        {
            const int pm_rank_tot = pm_id + controller->PP_MPI_size;
            const int rank_bytes =
                pm_pp_num[pm_id] * static_cast<int>(sizeof(int));
            MPI_Send(&min_corner_set[pm_id], sizeof(VECTOR), MPI_BYTE,
                     pm_rank_tot, PME_MPI_TAG_DOMAIN_MIN_CORNER,
                     MPI_COMM_WORLD);
            MPI_Send(&max_corner_set[pm_id], sizeof(VECTOR), MPI_BYTE,
                     pm_rank_tot, PME_MPI_TAG_DOMAIN_MAX_CORNER,
                     MPI_COMM_WORLD);
            MPI_Send(&pm_pp_num[pm_id], sizeof(int), MPI_BYTE, pm_rank_tot,
                     PME_MPI_TAG_DOMAIN_PP_COUNT, MPI_COMM_WORLD);
            MPI_Send(pm_pp_corres[pm_id].data(), rank_bytes, MPI_BYTE,
                     pm_rank_tot, PME_MPI_TAG_DOMAIN_PP_RANKS, MPI_COMM_WORLD);
            MPI_Send(&pm_dom_dec_split_num, sizeof(INT_VECTOR), MPI_BYTE,
                     pm_rank_tot, PME_MPI_TAG_DOMAIN_SPLIT, MPI_COMM_WORLD);
        }
    }
    else
    {
        if (controller->MPI_rank < controller->PP_MPI_size)
        {
            MPI_Recv(&pp_corres_pm_rank, sizeof(int), MPI_BYTE, 0,
                     PME_MPI_TAG_DOMAIN_PM_ASSIGNMENT, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            if (pp_corres_pm_rank < controller->PP_MPI_size ||
                pp_corres_pm_rank >=
                    controller->PP_MPI_size + controller->PM_MPI_size)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "Particle_Mesh::Send_Recv_Dom_Dec",
                    "Reason:\n\tPP rank %d received out-of-range PM rank %d\n",
                    controller->MPI_rank, pp_corres_pm_rank);
                return;
            }
        }
        else
        {
            const int local_pm_rank =
                controller->MPI_rank - controller->PP_MPI_size;
            if (local_pm_rank < 0 || local_pm_rank >= controller->PM_MPI_size ||
                pm_rank != local_pm_rank)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "Particle_Mesh::Send_Recv_Dom_Dec",
                    "Reason:\n\tglobal rank %d has inconsistent local PM rank "
                    "%d (expected %d)\n",
                    controller->MPI_rank, pm_rank, local_pm_rank);
                return;
            }
            MPI_Recv(&min_corner, sizeof(VECTOR), MPI_BYTE, 0,
                     PME_MPI_TAG_DOMAIN_MIN_CORNER, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Recv(&max_corner, sizeof(VECTOR), MPI_BYTE, 0,
                     PME_MPI_TAG_DOMAIN_MAX_CORNER, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Recv(&pm_corres_pp_num, sizeof(int), MPI_BYTE, 0,
                     PME_MPI_TAG_DOMAIN_PP_COUNT, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            if (pm_corres_pp_num <= 0 ||
                pm_corres_pp_num > controller->PP_MPI_size ||
                pm_corres_pp_num > INT_MAX / static_cast<int>(sizeof(int)))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "Particle_Mesh::Send_Recv_Dom_Dec",
                    "Reason:\n\tPM rank %d received invalid PP rank count %d\n",
                    pm_rank, pm_corres_pp_num);
                return;
            }
            try
            {
                pm_corres_pp_rank_set.resize(
                    static_cast<std::size_t>(pm_corres_pp_num));
                pm_corres_pp_atom_number.assign(
                    static_cast<std::size_t>(pm_corres_pp_num), 0);
                pm_corres_pp_atom_number_prefix.assign(
                    static_cast<std::size_t>(pm_corres_pp_num), 0);
            }
            catch (const std::bad_alloc&)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorMallocFailed, "Particle_Mesh::Send_Recv_Dom_Dec",
                    "Reason:\n\tcannot allocate the received PP rank list\n");
                return;
            }
            MPI_Recv(pm_corres_pp_rank_set.data(),
                     pm_corres_pp_num * static_cast<int>(sizeof(int)), MPI_BYTE,
                     0, PME_MPI_TAG_DOMAIN_PP_RANKS, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            try
            {
                std::vector<unsigned char> seen(
                    static_cast<std::size_t>(controller->PP_MPI_size), 0);
                for (int pp_rank : pm_corres_pp_rank_set)
                {
                    if (pp_rank < 0 || pp_rank >= controller->PP_MPI_size ||
                        seen[pp_rank])
                    {
                        controller->Throw_Formatted_SPONGE_Error(
                            spongeErrorValueErrorCommand,
                            "Particle_Mesh::Send_Recv_Dom_Dec",
                            "Reason:\n\tPM rank %d received invalid or "
                            "duplicate PP rank %d\n",
                            pm_rank, pp_rank);
                        return;
                    }
                    seen[pp_rank] = 1;
                }
            }
            catch (const std::bad_alloc&)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorMallocFailed, "Particle_Mesh::Send_Recv_Dom_Dec",
                    "Reason:\n\tcannot validate the received PP rank list\n");
                return;
            }
            MPI_Recv(&pm_dom_dec_split_num, sizeof(INT_VECTOR), MPI_BYTE, 0,
                     PME_MPI_TAG_DOMAIN_SPLIT, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            if (!PME_Split_Product_Equals(pm_dom_dec_split_num,
                                          controller->PM_MPI_size))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "Particle_Mesh::Send_Recv_Dom_Dec",
                    "Reason:\n\tPM rank %d received invalid PM split (%d, "
                    "%d, %d) for %d PM ranks\n",
                    pm_rank, pm_dom_dec_split_num.int_x,
                    pm_dom_dec_split_num.int_y, pm_dom_dec_split_num.int_z,
                    controller->PM_MPI_size);
                return;
            }
            if (pm_dom_dec_split_num.int_y == 1 &&
                pm_dom_dec_split_num.int_z == 1)
            {
                strcpy(this->FFT_MPI_TYPE, "SLAB");
            }
            else
            {
                strcpy(this->FFT_MPI_TYPE, "BRICK");
            }
        }
    }
#endif
}

void Particle_Mesh::Find_Neighbor_Domain(CONTROLLER* controller)
{
    for (int direction = 0; direction < 6; ++direction)
    {
        neighbor_num[direction] = 0;
        neighbor_dir[direction].clear();
    }
    if (!PM_MPI_size) return;
    if (controller->PM_MPI_size == 1 ||
        controller->MPI_rank < controller->PP_MPI_size)
    {
        return;
    }

    const int nx = pm_dom_dec_split_num.int_x;
    const int ny = pm_dom_dec_split_num.int_y;
    const int nz = pm_dom_dec_split_num.int_z;
    if (!PME_Split_Product_Equals(pm_dom_dec_split_num,
                                  controller->PM_MPI_size) ||
        pm_rank < 0 || pm_rank >= controller->PM_MPI_size)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Find_Neighbor_Domain",
            "Reason:\n\tinvalid PM split (%d, %d, %d), PM rank %d, or PM "
            "rank count %d\n",
            nx, ny, nz, pm_rank, controller->PM_MPI_size);
        return;
    }

    const int rank_id = pm_rank;
    const int i = rank_id % nx;
    const int j = (rank_id / nx) % ny;
    const int k = rank_id / (nx * ny);

    if (nx > 1)
    {
        neighbor_dir[0].push_back((i + 1) % nx + j * nx + k * nx * ny);
        neighbor_dir[1].push_back((i - 1 + nx) % nx + j * nx + k * nx * ny);
    }

    if (ny > 1)
    {
        neighbor_dir[2].push_back(i + ((j + 1) % ny) * nx + k * nx * ny);
        neighbor_dir[3].push_back(i + ((j - 1 + ny) % ny) * nx + k * nx * ny);
    }

    if (nz > 1)
    {
        neighbor_dir[4].push_back(i + j * nx + ((k + 1) % nz) * nx * ny);
        neighbor_dir[5].push_back(i + j * nx + ((k - 1 + nz) % nz) * nx * ny);
    }
    for (int direction = 0; direction < 6; ++direction)
        neighbor_num[direction] =
            static_cast<int>(neighbor_dir[direction].size());
}

enum PME_ATOM_ID_VALIDATION_ERROR
{
    PME_ATOM_IDS_VALID = 0,
    PME_ATOM_ID_OUT_OF_RANGE = 1,
    PME_ATOM_ID_DUPLICATE = 2,
    PME_ATOM_ID_MISSING = 3,
};

static __device__ __forceinline__ void Record_PME_Atom_ID_Error(int* error,
                                                                int code,
                                                                int index,
                                                                int value)
{
#ifdef GPU_ARCH_NAME
    if (atomicCAS(error, PME_ATOM_IDS_VALID, code) == PME_ATOM_IDS_VALID)
    {
        error[1] = index;
        error[2] = value;
    }
#elif defined(__GNUC__) || defined(__clang__)
    int expected = PME_ATOM_IDS_VALID;
    if (__atomic_compare_exchange_n(error, &expected, code, false,
                                    __ATOMIC_RELAXED, __ATOMIC_RELAXED))
    {
        error[1] = index;
        error[2] = value;
    }
#else
#pragma omp critical(sponge_pme_atom_id_error)
    {
        if (error[0] == PME_ATOM_IDS_VALID)
        {
            error[1] = index;
            error[2] = value;
            error[0] = code;
        }
    }
#endif
}

static __device__ __forceinline__ int Claim_PME_Global_Atom_ID(
    int* candidate_inverse, int global_id, int local_id)
{
    int* slot = candidate_inverse + global_id;
#ifdef GPU_ARCH_NAME
    return atomicCAS(slot, -1, local_id);
#elif defined(__GNUC__) || defined(__clang__)
    int expected = -1;
    if (__atomic_compare_exchange_n(slot, &expected, local_id, false,
                                    __ATOMIC_RELAXED, __ATOMIC_RELAXED))
        return -1;
    return expected;
#else
    int observed = -1;
#pragma omp critical(sponge_pme_atom_id_claim)
    {
        observed = slot[0];
        if (observed == -1) slot[0] = local_id;
    }
    return observed;
#endif
}

static __global__ void Build_PME_Inverse_Atom_ID_Candidate(
    const int* local_to_global, int* candidate_inverse, int* error,
    int atom_numbers)
{
    SIMPLE_DEVICE_FOR(local_id, atom_numbers)
    {
        const int global_id = local_to_global[local_id];
        if (global_id < 0 || global_id >= atom_numbers)
        {
            Record_PME_Atom_ID_Error(error, PME_ATOM_ID_OUT_OF_RANGE, local_id,
                                     global_id);
        }
        else if (Claim_PME_Global_Atom_ID(candidate_inverse, global_id,
                                          local_id) != -1)
        {
            Record_PME_Atom_ID_Error(error, PME_ATOM_ID_DUPLICATE, local_id,
                                     global_id);
        }
    }
}

static __global__ void Check_PME_Inverse_Atom_ID_Candidate(
    const int* candidate_inverse, int* error, int atom_numbers)
{
    SIMPLE_DEVICE_FOR(global_id, atom_numbers)
    {
        if (candidate_inverse[global_id] < 0)
            Record_PME_Atom_ID_Error(error, PME_ATOM_ID_MISSING, global_id, -1);
    }
}

static __global__ void crd_local_to_global(VECTOR* l_crd, VECTOR* g_crd,
                                           int* atom_id_l_g, int N)
{
#ifdef USE_GPU
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)
#else
#pragma omp parallel for
    for (int i = 0; i < N; i++)
#endif
    {
        int g_id = atom_id_l_g[i];
        g_crd[g_id] = l_crd[i];
    }
}

void Particle_Mesh::Get_Atoms(CONTROLLER* controller, VECTOR* pme_crd,
                              float* pme_charge, int pp_atom_numbers,
                              VECTOR* pp_crd, float* pp_charge, int* atom_local,
                              bool atom_number_label, bool charge_label,
                              bool crd_label, bool id_label)
{
    if (!PM_MPI_size)
    {
        return;
    }
    // 若单进程PM与PP共享同一进程下，什么也不做。共享进程下直接传入dd.crd等计算，不再重复拷贝内存
    if (controller->MPI_size == 1 && PM_MPI_size == 1)
    {
        return;
    }
#ifdef USE_MPI
    if (atom_number_label && !id_label)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
            "Reason:\n\ta changed PP atom-count layout must be accompanied "
            "by a complete atom-ID mapping\n");
        return;
    }
    const bool is_pp_rank = controller->MPI_rank < controller->PP_MPI_size;
    if (is_pp_rank)
    {
        if (pp_corres_pm_rank < controller->PP_MPI_size ||
            pp_corres_pm_rank >=
                controller->PP_MPI_size + controller->PM_MPI_size ||
            pp_atom_numbers < 0 || pp_atom_numbers > max_atom_numbers ||
            !PME_MPI_Byte_Count_Fits(pp_atom_numbers, sizeof(VECTOR)) ||
            !PME_MPI_Byte_Count_Fits(pp_atom_numbers, sizeof(float)) ||
            !PME_MPI_Byte_Count_Fits(pp_atom_numbers, sizeof(int)))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
                "Reason:\n\tPP rank %d has invalid PM rank %d or local atom "
                "count %d (global PME capacity %d)\n",
                controller->MPI_rank, pp_corres_pm_rank, pp_atom_numbers,
                max_atom_numbers);
            return;
        }
        if (!atom_number_label && reported_pp_atom_numbers != pp_atom_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
                "Reason:\n\tPP rank %d now owns %d atoms, but its last "
                "validated PME partition reported %d; atom counts and IDs "
                "must be refreshed together\n",
                controller->MPI_rank, pp_atom_numbers,
                reported_pp_atom_numbers);
            return;
        }
    }
    else if (pm_corres_pp_num <= 0 ||
             pm_corres_pp_rank_set.size() !=
                 static_cast<std::size_t>(pm_corres_pp_num) ||
             pm_corres_pp_atom_number.size() !=
                 static_cast<std::size_t>(pm_corres_pp_num) ||
             pm_corres_pp_atom_number_prefix.size() !=
                 static_cast<std::size_t>(pm_corres_pp_num))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
            "Reason:\n\tPM rank correspondence vectors are not initialized "
            "to the received PP rank count\n");
        return;
    }

    // 先阻塞通信原子数
    if (is_pp_rank)
    {
        if (atom_number_label)
        {
            MPI_Send(&pp_atom_numbers, sizeof(int), MPI_BYTE, pp_corres_pm_rank,
                     PME_MPI_TAG_ATOM_COUNT, MPI_COMM_WORLD);
            reported_pp_atom_numbers = pp_atom_numbers;
        }
    }
    else
    {
        if (atom_number_label)
        {
            std::vector<int> received_counts;
            std::vector<int> received_prefixes;
            try
            {
                received_counts.resize(
                    static_cast<std::size_t>(pm_corres_pp_num));
            }
            catch (const std::bad_alloc&)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorMallocFailed, "Particle_Mesh::Get_Atoms",
                    "Reason:\n\tcannot stage PP atom counts for validation\n");
                return;
            }
            // 接收PP进程对应的原子数
            for (int i = 0; i < pm_corres_pp_num; ++i)
            {
                const int pp_rank = pm_corres_pp_rank_set[i];
                MPI_Recv(&received_counts[i], sizeof(int), MPI_BYTE, pp_rank,
                         PME_MPI_TAG_ATOM_COUNT, MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
            }
            std::string count_error;
            if (!Try_Build_PME_Atom_Count_Prefix(
                    received_counts, max_atom_numbers, &received_prefixes,
                    &count_error))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
                    "Reason:\n\tPM rank %d received an invalid PP atom-count "
                    "partition for global capacity %d: %s\n",
                    pm_rank, max_atom_numbers, count_error.c_str());
                return;
            }
            pm_corres_pp_atom_number = std::move(received_counts);
            pm_corres_pp_atom_number_prefix = std::move(received_prefixes);
            this->atom_numbers = max_atom_numbers;
            atom_mapping_is_valid = false;
        }
    }
    if (!is_pp_rank)
    {
        int expected_prefix = 0;
        for (int i = 0; i < pm_corres_pp_num; ++i)
        {
            const int count = pm_corres_pp_atom_number[i];
            if (pm_corres_pp_atom_number_prefix[i] != expected_prefix ||
                count < 0 || count > max_atom_numbers - expected_prefix ||
                !PME_MPI_Byte_Count_Fits(count, sizeof(VECTOR)) ||
                !PME_MPI_Byte_Count_Fits(count, sizeof(float)) ||
                !PME_MPI_Byte_Count_Fits(count, sizeof(int)))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
                    "Reason:\n\tPM rank %d has an invalid committed atom-count "
                    "entry %d (count %d, prefix %d, expected prefix %d)\n",
                    pm_rank, i, count, pm_corres_pp_atom_number_prefix[i],
                    expected_prefix);
                return;
            }
            expected_prefix += count;
        }
        if (expected_prefix != max_atom_numbers ||
            this->atom_numbers != max_atom_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
                "Reason:\n\tPM rank %d atom-count partition covers %d atoms, "
                "but the global PME count is %d (active count %d)\n",
                pm_rank, expected_prefix, max_atom_numbers, this->atom_numbers);
            return;
        }
    }
    // 通信坐标与电荷

    if (is_pp_rank)
    {
        D_MPI_GroupStart();
        if (crd_label && pp_atom_numbers > 0)
        {
            D_MPI_Send(pp_crd,
                       PME_MPI_Byte_Count(pp_atom_numbers, sizeof(VECTOR)),
                       D_MPI_BYTE, pp_corres_pm_rank, PME_MPI_TAG_COORDINATES,
                       controller->D_MPI_COMM_WORLD, pm_stream);
        }
        if (charge_label && pp_atom_numbers > 0)
        {
            D_MPI_Send(pp_charge,
                       PME_MPI_Byte_Count(pp_atom_numbers, sizeof(float)),
                       D_MPI_BYTE, pp_corres_pm_rank, PME_MPI_TAG_CHARGES,
                       controller->D_MPI_COMM_WORLD, pm_stream);
        }
        if (id_label && pp_atom_numbers > 0)
        {
            D_MPI_Send(atom_local,
                       PME_MPI_Byte_Count(pp_atom_numbers, sizeof(int)),
                       D_MPI_BYTE, pp_corres_pm_rank, PME_MPI_TAG_ATOM_IDS,
                       controller->D_MPI_COMM_WORLD, pm_stream);
        }
        D_MPI_GroupEnd();
#ifdef USE_GPU
        deviceStreamSynchronize(pm_stream);
#endif
    }
    else
    {
        if (id_label) atom_mapping_is_valid = false;
        D_MPI_GroupStart();
        for (int i = 0; i < pm_corres_pp_num; ++i)
        {
            int pp_rank = pm_corres_pp_rank_set[i];
            const int received_count = pm_corres_pp_atom_number[i];
            if (received_count == 0) continue;
            if (crd_label)
            {
                D_MPI_Recv(pme_crd + pm_corres_pp_atom_number_prefix[i],
                           PME_MPI_Byte_Count(received_count, sizeof(VECTOR)),
                           D_MPI_BYTE, pp_rank, PME_MPI_TAG_COORDINATES,
                           controller->D_MPI_COMM_WORLD, pm_stream);
            }
            if (charge_label)
            {
                D_MPI_Recv(pme_charge + pm_corres_pp_atom_number_prefix[i],
                           PME_MPI_Byte_Count(received_count, sizeof(float)),
                           D_MPI_BYTE, pp_rank, PME_MPI_TAG_CHARGES,
                           controller->D_MPI_COMM_WORLD, pm_stream);
            }
            if (id_label)
            {
                D_MPI_Recv(atom_id_l_g + pm_corres_pp_atom_number_prefix[i],
                           PME_MPI_Byte_Count(received_count, sizeof(int)),
                           D_MPI_BYTE, pp_rank, PME_MPI_TAG_ATOM_IDS,
                           controller->D_MPI_COMM_WORLD, pm_stream);
            }
        }
        D_MPI_GroupEnd();
#ifdef USE_GPU
        deviceStreamSynchronize(pm_stream);
#endif
        // Validate the entire staged local-to-global mapping before either
        // committing its inverse or scattering coordinates through it.
        if (id_label)
        {
            deviceMemset(atom_id_g_l_candidate, -1,
                         sizeof(int) * this->atom_numbers);
            deviceMemset(atom_id_validation_error, 0, 3 * sizeof(int));
            Launch_Device_Kernel(
                Build_PME_Inverse_Atom_ID_Candidate,
                (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, atom_id_l_g,
                atom_id_g_l_candidate, atom_id_validation_error,
                this->atom_numbers);
            Launch_Device_Kernel(
                Check_PME_Inverse_Atom_ID_Candidate,
                (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, atom_id_g_l_candidate,
                atom_id_validation_error, this->atom_numbers);

            int validation_error[3] = {PME_ATOM_IDS_VALID, -1, -1};
            deviceMemcpy(validation_error, atom_id_validation_error,
                         sizeof(validation_error), deviceMemcpyDeviceToHost);
            if (validation_error[0] != PME_ATOM_IDS_VALID)
            {
                const char* reason = "unknown validation error";
                if (validation_error[0] == PME_ATOM_ID_OUT_OF_RANGE)
                    reason = "out-of-range global atom ID";
                else if (validation_error[0] == PME_ATOM_ID_DUPLICATE)
                    reason = "duplicate global atom ID";
                else if (validation_error[0] == PME_ATOM_ID_MISSING)
                    reason = "missing global atom ID";
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
                    "Reason:\n\tPM rank %d received an invalid global atom-ID "
                    "mapping: %s (index %d, value %d, code %d)\n",
                    pm_rank, reason, validation_error[1], validation_error[2],
                    validation_error[0]);
                return;
            }
            deviceMemcpy(atom_id_g_l, atom_id_g_l_candidate,
                         sizeof(int) * this->atom_numbers,
                         deviceMemcpyDeviceToDevice);
            atom_mapping_is_valid = true;
        }
        if (crd_label && !atom_mapping_is_valid)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Get_Atoms",
                "Reason:\n\tPME coordinates cannot be scattered before a "
                "complete atom-ID mapping has been validated\n");
            return;
        }
        if (crd_label && this->atom_numbers > 0)
        {
            Launch_Device_Kernel(
                crd_local_to_global,
                (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, pme_crd, g_crd,
                atom_id_l_g, this->atom_numbers);
        }
    }
#endif
}

// 目前只做单进程PME，暂时不考虑ghost，同样也不考虑get_local

void Particle_Mesh::Send_Recv_Force(CONTROLLER* controller, VECTOR* frc,
                                    VECTOR* pp_frc, int pp_atom_numbers)
{
    if (!PM_MPI_size) return;
#ifdef USE_MPI

    if (controller->MPI_rank < controller->PP_MPI_size)
    {
        if (pp_atom_numbers < 0 || pp_atom_numbers > max_atom_numbers ||
            !PME_MPI_Byte_Count_Fits(pp_atom_numbers, sizeof(VECTOR)) ||
            reported_pp_atom_numbers != pp_atom_numbers ||
            pp_corres_pm_rank < controller->PP_MPI_size ||
            pp_corres_pm_rank >=
                controller->PP_MPI_size + controller->PM_MPI_size)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Send_Recv_Force",
                "Reason:\n\tPP rank %d has invalid atom count %d or PM rank "
                "%d\n",
                controller->MPI_rank, pp_atom_numbers, pp_corres_pm_rank);
            return;
        }
        if (pp_atom_numbers > 0)
        {
            D_MPI_GroupStart();
            D_MPI_Recv(frc, PME_MPI_Byte_Count(pp_atom_numbers, sizeof(VECTOR)),
                       D_MPI_BYTE, pp_corres_pm_rank, PME_MPI_TAG_FORCES,
                       controller->D_MPI_COMM_WORLD, pm_stream);
            D_MPI_GroupEnd();
#ifdef USE_GPU
            deviceStreamSynchronize(pm_stream);
#endif
            Launch_Device_Kernel(
                device_add_force,
                (pp_atom_numbers + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, pp_atom_numbers, 1,
                pp_frc, frc);
        }
    }
    else
    {
        if (pm_corres_pp_num <= 0 ||
            pm_corres_pp_rank_set.size() !=
                static_cast<std::size_t>(pm_corres_pp_num) ||
            pm_corres_pp_atom_number.size() !=
                static_cast<std::size_t>(pm_corres_pp_num) ||
            pm_corres_pp_atom_number_prefix.size() !=
                static_cast<std::size_t>(pm_corres_pp_num))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Send_Recv_Force",
                "Reason:\n\tPM rank correspondence vectors are "
                "inconsistent\n");
            return;
        }
        int expected_prefix = 0;
        for (int i = 0; i < pm_corres_pp_num; ++i)
        {
            const int pp_rank = pm_corres_pp_rank_set[i];
            const int prefix = pm_corres_pp_atom_number_prefix[i];
            const int count = pm_corres_pp_atom_number[i];
            if (pp_rank < 0 || pp_rank >= controller->PP_MPI_size ||
                prefix != expected_prefix || count < 0 ||
                count > max_atom_numbers - expected_prefix ||
                !PME_MPI_Byte_Count_Fits(count, sizeof(VECTOR)))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "Particle_Mesh::Send_Recv_Force",
                    "Reason:\n\tPM rank %d has invalid PP rank/count/prefix "
                    "entry (%d, %d, %d)\n",
                    pm_rank, pp_rank, count, prefix);
                return;
            }
            expected_prefix += count;
        }
        if (expected_prefix != max_atom_numbers ||
            atom_numbers != max_atom_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "Particle_Mesh::Send_Recv_Force",
                "Reason:\n\tPM rank %d force partition covers %d atoms, but "
                "the global PME count is %d (active count %d)\n",
                pm_rank, expected_prefix, max_atom_numbers, atom_numbers);
            return;
        }
        for (int i = 0; i < pm_corres_pp_num; ++i)
        {
            const int pp_rank = pm_corres_pp_rank_set[i];
            const int prefix = pm_corres_pp_atom_number_prefix[i];
            const int count = pm_corres_pp_atom_number[i];
            if (count == 0) continue;
            D_MPI_GroupStart();
            D_MPI_Send(frc + prefix, PME_MPI_Byte_Count(count, sizeof(VECTOR)),
                       D_MPI_BYTE, pp_rank, PME_MPI_TAG_FORCES,
                       controller->D_MPI_COMM_WORLD, pm_stream);
            D_MPI_GroupEnd();
#ifdef USE_GPU
            deviceStreamSynchronize(pm_stream);
#endif
        }
    }
#endif
}

void Particle_Mesh::Create_Stream() { deviceStreamCreate(&pm_stream); }

void Particle_Mesh::Destroy_Stream() { deviceStreamDestroy(pm_stream); }

static __global__ void MPI_PME_Excluded_Force_With_Atom_Energy_Correction(
    const int atom_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const float* charge, const float pme_beta,
    const int* excluded_list_start, const int* excluded_list,
    const int* excluded_atom_numbers, VECTOR* frc, float* atom_ene,
    float* this_ene, LTMatrix3* atom_virial, int need_energy, int need_virial,
    const PME_EXCLUSION_DEPENDENCY_STATE* exclusion_dependencies,
    const float factor)
{
    SIMPLE_DEVICE_FOR(local_i, atom_numbers)
    {
        int excluded_numbers = excluded_atom_numbers[local_i];
        if (excluded_numbers > 0)
        {
            int list_start = excluded_list_start[local_i];
            int list_end = list_start + excluded_numbers;
            int partner_index;

            float charge_i = charge[local_i];
            float charge_j;
            VECTOR r1 = crd[local_i], r2;
            VECTOR dr;
            float dr2;

            float frc_abs = 0.;
            VECTOR frc_lin;
            VECTOR frc_record = {0., 0., 0.};
            LTMatrix3 virial_record = {0, 0, 0, 0, 0, 0};
            float ene_lin = 0.;

            for (int i = list_start; i < list_end; i = i + 1)
            {
                partner_index = excluded_list[i];
                if (partner_index >= 0)
                {
                    r2 = crd[partner_index];
                    charge_j = charge[partner_index];
                }
                else
                {
                    const int dependency = -partner_index - 1;
                    r2 = exclusion_dependencies[dependency].crd;
                    charge_j = exclusion_dependencies[dependency].charge;
                }
                if (!PairwiseInteraction::Coulomb_Is_Active(charge_i, charge_j))
                {
                    continue;
                }

                dr = Get_Periodic_Displacement(r2, r1, cell, rcell);
                dr2 = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
                const PME_Excluded_Radial_Kernels radial =
                    Get_PME_Excluded_Radial_Kernels(dr2, pme_beta);
                const float charge_product = charge_i * charge_j;
                frc_abs = charge_product * radial.force;
                frc_lin = frc_abs * dr;
                // A single PP rank consumes the original triangular CSR and
                // therefore owns both force updates.  Multiple PP ranks use
                // the symmetric CSR; each owner computes only its own force,
                // including when the partner came from the dedicated remote
                // dependency exchange.
                if (factor > 0.6f && partner_index >= 0)
                    atomicAdd(frc + partner_index, -frc_lin);
                frc_record = frc_record + frc_lin;
                if (need_energy)
                    ene_lin -= factor * charge_product * radial.energy;
                if (need_virial)
                    virial_record =
                        virial_record -
                        factor * Get_Virial_From_Force_Dis(frc_lin, dr);
            }  // atom_j cycle
            atomicAdd(frc + local_i, frc_record);
            if (need_energy)
            {
                atomicAdd(atom_ene + local_i, ene_lin);
                this_ene[local_i] = ene_lin;
            }
            if (need_virial) atomicAdd(atom_virial + local_i, virial_record);
        }  // if need excluded
    }
}

void Particle_Mesh::MPI_PME_Excluded_Force_With_Atom_Energy(
    const int local_atom_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const float* charge, const int* excluded_list_start,
    const int* excluded_list, const int* excluded_atom_numbers,
    const PME_EXCLUSION_DEPENDENCY_STATE* exclusion_dependencies, VECTOR* frc,
    int need_energy, float* atom_ene, int need_virial, LTMatrix3* atom_virial)
{
    if (is_initialized && calculate_excluded_part)
    {
        if (need_energy)
        {
            // Step_Print reduces the complete allocation, not just the
            // current owned prefix.  Domain migration can shrink that prefix,
            // so clearing only local_atom_numbers would retain stale energy
            // from atoms previously owned by this PP rank.
            deviceMemset(d_correction_atom_energy, 0,
                         sizeof(float) * atom_numbers);
        }

        // A valid domain decomposition may leave a PP rank with no owned
        // atoms.  Such ranks still participate in the dependency collectives,
        // and must clear their printable energy state, but must not issue a
        // zero-block device launch.
        if (local_atom_numbers == 0) return;

        Launch_Device_Kernel(
            MPI_PME_Excluded_Force_With_Atom_Energy_Correction,
            (local_atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, local_atom_numbers, crd,
            cell, rcell, charge, beta, excluded_list_start, excluded_list,
            excluded_atom_numbers, frc, atom_ene, d_correction_atom_energy,
            atom_virial, need_energy, need_virial, exclusion_dependencies,
            exclude_factor);
    }
}

void Particle_Mesh::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    // 单进程, PM与PP共享同一进程情况
    if (CONTROLLER::MPI_size == 1 && CONTROLLER::PM_MPI_size == 1)
    {
        Sum_Of_List(d_correction_atom_energy, d_correction_ene, atom_numbers);
        Sum_Of_List(d_direct_atom_energy, d_direct_ene, atom_numbers);
        deviceMemcpy(&direct_ene, d_direct_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(&correction_ene, d_correction_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(&self_ene, d_self_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(&reciprocal_ene, d_reciprocal_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
        ee_ene = direct_ene + reciprocal_ene + self_ene + correction_ene;
        controller->Step_Print("PM", ee_ene, true);
        if (print_detail)
        {
            controller->Step_Print("PM_direct", direct_ene);
            controller->Step_Print("PM_reciprocal", reciprocal_ene);
            controller->Step_Print("PM_self", self_ene);
            controller->Step_Print("PM_correction", correction_ene);
        }
        return;
    }
    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        Sum_Of_List(d_correction_atom_energy, d_correction_ene, atom_numbers);
        Sum_Of_List(d_direct_atom_energy, d_direct_ene, atom_numbers);
        self_ene = 0;
        reciprocal_ene = 0;
        deviceMemcpy(&direct_ene, d_direct_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(&correction_ene, d_correction_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
    }
    else
    {
        direct_ene = 0;
        correction_ene = 0;
        deviceMemcpy(&self_ene, d_self_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(&reciprocal_ene, d_reciprocal_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
    }
    ee_ene = direct_ene + reciprocal_ene + self_ene + correction_ene;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &ee_ene, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
    controller->Step_Print("PM", ee_ene, true);
    if (print_detail)
    {
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &self_ene, 1, MPI_FLOAT, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &reciprocal_ene, 1, MPI_FLOAT, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &direct_ene, 1, MPI_FLOAT, MPI_SUM,
                      MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &correction_ene, 1, MPI_FLOAT, MPI_SUM,
                      MPI_COMM_WORLD);
#endif
        controller->Step_Print("PM_direct", direct_ene);
        controller->Step_Print("PM_reciprocal", reciprocal_ene);
        controller->Step_Print("PM_self", self_ene);
        controller->Step_Print("PM_correction", correction_ene);
    }
}

void Particle_Mesh::reset_global_force(
    int no_direct_interaction_virtual_atom_numbers)
{
    deviceMemset(g_frc, 0,
                 sizeof(VECTOR) * (atom_numbers +
                                   no_direct_interaction_virtual_atom_numbers));
}

static __global__ void add_global_to_local_force(const VECTOR* g_frc,
                                                 VECTOR* l_frc,
                                                 const int* atom_id_g_l, int N)
{
#ifdef USE_GPU
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N)
#else
#pragma omp parallel for
    for (int i = 0; i < N; i++)
#endif
    {
        int l_id = atom_id_g_l[i];
        l_frc[l_id] = l_frc[l_id] + g_frc[i];
    }
}

void Particle_Mesh::add_force_g_to_l(VECTOR* l_frc)
{
    Launch_Device_Kernel(add_global_to_local_force,
                         (atom_numbers + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, g_frc, l_frc,
                         atom_id_g_l, atom_numbers);
}
