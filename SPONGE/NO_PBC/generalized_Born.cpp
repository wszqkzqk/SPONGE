#include "generalized_Born.h"

#include "../Lennard_Jones_force/pair_activity.h"
#include "../xponge/load/native/gb.hpp"
#include "../xponge/xponge.h"

#include <cerrno>

namespace
{
float GB_Parse_Finite_Command(CONTROLLER* controller, const char* module_name,
                              const char* command, bool positive)
{
    const char* token = controller->Command(module_name, command);
    errno = 0;
    char* end = NULL;
    const float value = std::strtof(token, &end);
    if (errno == ERANGE || end == token || end[0] != '\0' ||
        !Float_Memory_Is_Finite(&value) ||
        !Float_Memory_Is_Zero_Or_Normal(&value) ||
        (positive && !(value > 0.0f)))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\t%s.%s must be representable as a finite %sfloat; "
            "got \"%s\"\n",
            module_name, command,
            positive ? "positive normal " : "zero or normal ", token);
    }
    return value;
}
}  // namespace

static __global__ void Effective_Born_Radii_Factor_Device(
    const int atom_numbers, const VECTOR* crd, const float cutoff_square,
    const float* self_radius, const float* other_radius,
    float* effective_radius, int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    int atom_j = blockDim.y * blockIdx.y + threadIdx.y;
    if (atom_i < atom_numbers && atom_j < atom_numbers && atom_i != atom_j)
    {
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
        for (int atom_j = 0; atom_j < atom_numbers; atom_j++)
        {
            if (atom_j == atom_i) continue;
#endif
        VECTOR dr = crd[atom_j] - crd[atom_i];
        float dr2 = dr * dr;
        if (dr.x == 0.0f && dr.y == 0.0f && dr.z == 0.0f)
        {
            PairwiseInteraction::Fail_Exact_Overlap(
                atom_i, atom_j,
                PairwiseInteraction::PAIR_COMPONENT_GENERALIZED_BORN,
                pair_overlap_error);
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }
        if (dr2 < cutoff_square)
        {
            float R = sqrtf(dr2);
            float self_radii = self_radius[atom_i];
            float other_radii = other_radius[atom_j];
            float inner_distance = R - other_radii;
            float outer_distance = R + other_radii;
            float U = 1;
            float L = 1;
            if (self_radii <= inner_distance)
            {
                L = 1.0 / inner_distance;
                U = 1.0 / outer_distance;
            }
            else if (self_radii < outer_distance)
            {
                L = 1.0 / self_radii;
                U = 1.0 / outer_distance;
            }

            float temp =
                0.125 / R *
                (4 * R * (L - U) + dr2 * (U * U - L * L) + 2 * logf(U / L) +
                 other_radii * other_radii * (L * L - U * U));
            atomicAdd(effective_radius + atom_i, temp);
        }
    }
}

static __global__ void Effective_Born_Radii_Device(const int atom_numbers,
                                                   const float* self_radius,
                                                   float* effective_radius)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
#endif
    {
        effective_radius[atom_i] =
            1.0 / (1.0 / self_radius[atom_i] - effective_radius[atom_i]);
    }
}

static __host__ __device__ __forceinline__ bool GB_Is_Positive_Normal_Float(
    const float value)
{
    // The backend name is visible in both host and device passes.  Gate the
    // intrinsic on the compiler's device-pass marker instead.
#if defined(__CUDA_ARCH__) || \
    (defined(__HIP_DEVICE_COMPILE__) && __HIP_DEVICE_COMPILE__)
    const unsigned int bits = __float_as_uint(value);
#else
    unsigned int bits = 0;
    static_assert(sizeof(bits) == sizeof(value));
    memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
#endif
    return (bits & 0x80000000U) == 0U &&
           (bits & 0x7f800000U) != 0U &&
           (bits & 0x7f800000U) != 0x7f800000U;
}

static __global__ void Validate_Effective_Born_Radii_Device(
    const int atom_numbers, const float* effective_radius,
    int* effective_radius_error)
{
#ifdef USE_GPU
    const int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
#endif
    {
        const float radius = effective_radius[atom_i];
        if (!GB_Is_Positive_Normal_Float(radius))
        {
#ifdef GPU_ARCH_NAME
            printf(
                "Fatal SPONGE generalized-Born error: atom %d has invalid "
                "effective radius %.9g; the pairwise-descreening denominator "
                "must produce a positive normal float.\n",
                atom_i, radius);
#if defined(USE_CUDA)
            asm volatile("trap;");
#elif defined(USE_HIP)
            __builtin_trap();
#endif
#else
#if defined(__GNUC__) || defined(__clang__)
            int expected = -1;
            __atomic_compare_exchange_n(effective_radius_error, &expected,
                                        atom_i, false, __ATOMIC_RELAXED,
                                        __ATOMIC_RELAXED);
#else
#pragma omp critical(sponge_gb_effective_radius_error)
            {
                if (effective_radius_error[0] == -1)
                    effective_radius_error[0] = atom_i;
            }
#endif
#endif
        }
    }
}

static __global__ void GB_inej_Force_Energy_Device(
    const int atom_numbers, const VECTOR* crd, const float* charge,
    const float* effective_radius, const float epsilon_1_minus_1,
    VECTOR* frc, float* atom_ene, float* de_da, float* this_ene)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    int atom_j = atom_i + 1 + blockDim.y * blockIdx.y + threadIdx.y;
    if (atom_i < atom_numbers && atom_j < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
        for (int atom_j = atom_i + 1; atom_j < atom_numbers; atom_j++)
#endif
    {
        VECTOR dr = crd[atom_j] - crd[atom_i];
        float dr2 = dr * dr;
        {
            SADfloat<3> R(sqrtf(dr2), 0);
            SADfloat<3> a1(effective_radius[atom_i], 1);
            SADfloat<3> a2(effective_radius[atom_j], 2);

            SADfloat<3> born_radii_ij_square = a1 * a2;
            SADfloat<3> R2 = R * R;
            SADfloat<3> D = expf(-0.25f * R2 / born_radii_ij_square);
            SADfloat<3> temp_ene = charge[atom_i] * charge[atom_j] *
                                   epsilon_1_minus_1 /
                                   sqrtf(R2 + born_radii_ij_square * D);

            VECTOR temp_frc = -temp_ene.dval[0] / R.val * dr;

            atomicAdd(&frc[atom_j].x, temp_frc.x);
            atomicAdd(&frc[atom_j].y, temp_frc.y);
            atomicAdd(&frc[atom_j].z, temp_frc.z);
            atomicAdd(&frc[atom_i].x, -temp_frc.x);
            atomicAdd(&frc[atom_i].y, -temp_frc.y);
            atomicAdd(&frc[atom_i].z, -temp_frc.z);

            atomicAdd(&de_da[atom_i], temp_ene.dval[1]);
            atomicAdd(&de_da[atom_j], temp_ene.dval[2]);

            atomicAdd(&atom_ene[atom_i], temp_ene.val);
            atomicAdd(&this_ene[atom_i], temp_ene.val);
        }
    }
}

static __global__ void GB_ieqj_Force_Energy_Device(
    const int atom_numbers, const float* charge,
    const float* effective_radius, const float epsilon_1_minus_1_half,
    float* atom_ene, float* de_da, float* this_ene)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
#endif
    {
        float f_1 = 1 / effective_radius[atom_i];
        float temp_ene = charge[atom_i];
        temp_ene *= temp_ene;
        temp_ene *= f_1 * epsilon_1_minus_1_half;

        atomicAdd(&de_da[atom_i], -temp_ene * f_1);
        atomicAdd(&atom_ene[atom_i], temp_ene);
        atomicAdd(&this_ene[atom_i], temp_ene);
    }
}

static __global__ void GB_accumulate_Force_Energy_Device(
    const int atom_numbers, const VECTOR* crd, const float cutoff_square,
    const float* self_radius, const float* other_radius,
    const float* effective_raius, const float* de_da, VECTOR* frc)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    int atom_j = blockDim.y * blockIdx.y + threadIdx.y;
    if (atom_i < atom_numbers && atom_j < atom_numbers && atom_i != atom_j)
    {
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
        for (int atom_j = 0; atom_j < atom_numbers; atom_j++)
        {
            if (atom_j == atom_i) continue;
#endif
        VECTOR dr = crd[atom_j] - crd[atom_i];
        float dr2 = dr * dr;
        if (dr2 < cutoff_square)
        {
            SADfloat<1> R(sqrtf(dr2), 0);
            float self_radii = self_radius[atom_i];
            float other_radii = other_radius[atom_j];
            SADfloat<1> inner_distance = R - other_radii;
            SADfloat<1> outer_distance = R + other_radii;
            SADfloat<1> U = 1.0f;
            SADfloat<1> L = 1.0f;
            if (self_radii <= inner_distance.val)
            {
                L = 1.0f / inner_distance;
                U = 1.0f / outer_distance;
            }
            else if (self_radii < outer_distance.val)
            {
                L = 1.0f / self_radii;
                U = 1.0f / outer_distance;
            }

            SADfloat<1> temp = 0.125f / R *
                               (4.0f * R * (L - U) + R * R * (U * U - L * L) +
                                2.0f * logf(U / L) +
                                other_radii * other_radii * (L * L - U * U));

            float reff = effective_raius[atom_i];
            VECTOR temp_frc =
                reff * reff * temp.dval[0] * de_da[atom_i] / R.val * dr;
            atomicAdd(&frc[atom_i].x, temp_frc.x);
            atomicAdd(&frc[atom_i].y, temp_frc.y);
            atomicAdd(&frc[atom_i].z, temp_frc.z);
            atomicAdd(&frc[atom_j].x, -temp_frc.x);
            atomicAdd(&frc[atom_j].y, -temp_frc.y);
            atomicAdd(&frc[atom_j].z, -temp_frc.z);
        }
    }
}

void GENERALIZED_BORN_INFORMATION::Malloc()
{
    Malloc_Safely((void**)&h_GB_energy_atom, sizeof(float) * atom_numbers);
    Malloc_Safely((void**)&h_GB_self_radius, sizeof(float) * atom_numbers);
    Malloc_Safely((void**)&h_GB_other_radius, sizeof(float) * atom_numbers);

    memset(h_GB_energy_atom, 0, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_GB_energy_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_GB_energy_atom,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_GB_self_radius,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_GB_other_radius,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_GB_effective_radius,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_dE_da, sizeof(float) * atom_numbers);
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_pair_overlap_error, 3 * sizeof(int));
    Device_Malloc_Safely((void**)&d_effective_radius_error, sizeof(int));
#endif
    deviceMemset(d_GB_energy_sum, 0, sizeof(float));
    deviceMemset(d_GB_energy_atom, 0, sizeof(float) * atom_numbers);
}

void GENERALIZED_BORN_INFORMATION::Initial(CONTROLLER* controller,
                                           int expected_atom_numbers,
                                           float default_radii_cutoff,
                                           const char* module_name)
{
    if (is_initialized)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tthe generalized-Born module cannot be initialized "
            "twice without first releasing its state\n");
        return;
    }
    this->controller = controller;
    const char* selected_module_name = module_name == NULL ? "gb" : module_name;
    if (strlen(selected_module_name) >= sizeof(this->module_name))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tthe generalized-Born module name is too long\n");
        return;
    }
    strcpy(this->module_name, selected_module_name);
    if (expected_atom_numbers <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tGB requires a positive runtime atom count; got %d\n",
            expected_atom_numbers);
        return;
    }
    controller[0].printf(
        "START INITIALIZING STANDARD GENERALIZED BORN INFORMATION:\n");
    relative_dielectric_constant = 78.5;
    if (controller->Command_Exist(this->module_name, "epsilon"))
    {
        relative_dielectric_constant = GB_Parse_Finite_Command(
            controller, this->module_name, "epsilon", true);
    }

    radii_offset = 0.09;
    if (controller->Command_Exist(this->module_name, "radii_offset"))
    {
        radii_offset = GB_Parse_Finite_Command(
            controller, this->module_name, "radii_offset", false);
    }

    if (!Float_Memory_Is_Finite(&default_radii_cutoff) ||
        !Float_Memory_Is_Normal(&default_radii_cutoff) ||
        !(default_radii_cutoff > 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tthe default GB radii cutoff must be a finite "
            "positive normal float; got %.9g\n",
            default_radii_cutoff);
    }
    radii_cutoff = default_radii_cutoff;
    if (controller->Command_Exist(this->module_name, "radii_cutoff"))
    {
        radii_cutoff = GB_Parse_Finite_Command(
            controller, this->module_name, "radii_cutoff", true);
    }
    const double radii_cutoff_square_double =
        static_cast<double>(radii_cutoff) * radii_cutoff;
    if (!Double_Memory_Is_Finite(&radii_cutoff_square_double) ||
        radii_cutoff_square_double >
            static_cast<double>(std::numeric_limits<float>::max()))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\t%s.radii_cutoff %.9g has a square outside the "
            "finite float range used by the descreening kernel\n",
            this->module_name, radii_cutoff);
        return;
    }
    radii_cutoff_square = static_cast<float>(radii_cutoff_square_double);
    if (!Float_Memory_Is_Normal(&radii_cutoff_square))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\t%s.radii_cutoff %.9g has a square that is not a "
            "normal float on FTZ backends\n",
            this->module_name, radii_cutoff);
        return;
    }

    const auto& gb_system = Xponge::system.generalized_born;
    Xponge::GeneralizedBorn local_gb;
    const Xponge::GeneralizedBorn* gb_to_use = NULL;
    if (module_name == NULL)
    {
        gb_to_use = &gb_system;
    }
    else if (controller->Command_Exist(this->module_name, "in_file"))
    {
        Xponge::Native_Load_Generalized_Born(
            &local_gb, controller, expected_atom_numbers, this->module_name);
        gb_to_use = &local_gb;
    }

    if (gb_to_use == NULL || gb_to_use->radius.empty())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tGB requires a non-empty radius and scale-factor "
            "table\n");
        return;
    }
    if (gb_to_use->radius.size() != gb_to_use->scale_factor.size() ||
        gb_to_use->radius.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tGB radius/scale table sizes are inconsistent or "
            "unrepresentable: %zu/%zu\n",
            gb_to_use->radius.size(), gb_to_use->scale_factor.size());
        return;
    }

    const int validated_atom_numbers =
        static_cast<int>(gb_to_use->radius.size());
    if (validated_atom_numbers != expected_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tGB table atom count %d differs from the runtime "
            "atom count %d\n",
            validated_atom_numbers, expected_atom_numbers);
        return;
    }
    std::vector<float> validated_self;
    std::vector<float> validated_other;
    try
    {
        validated_self.resize(validated_atom_numbers);
        validated_other.resize(validated_atom_numbers);
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tthe GB atom table exceeds the maximum supported "
            "container size\n");
        return;
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat,
            "GENERALIZED_BORN_INFORMATION::Initial",
            "Reason:\n\tfailed to allocate validated GB atom parameters\n");
        return;
    }
    for (int i = 0; i < validated_atom_numbers; i++)
    {
        const float radius = gb_to_use->radius[i];
        const float scale = gb_to_use->scale_factor[i];
        const double self = static_cast<double>(radius) - radii_offset;
        const double other = static_cast<double>(scale) * self;
        if (!Float_Memory_Is_Finite(&radius) ||
            !Float_Memory_Is_Normal(&radius) || !(radius > 0.0f) ||
            !Float_Memory_Is_Finite(&scale) ||
            !Float_Memory_Is_Zero_Or_Normal(&scale) || scale < 0.0f ||
            !Double_Memory_Is_Finite(&self) || !(self > 0.0) ||
            self > std::numeric_limits<float>::max() ||
            !Double_Memory_Is_Finite(&other) || other < 0.0 ||
            other > std::numeric_limits<float>::max())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat,
                "GENERALIZED_BORN_INFORMATION::Initial",
                "Reason:\n\tGB atom %d has invalid radius/scale/offset "
                "values %.9g/%.9g/%.9g\n",
                i, radius, scale, radii_offset);
            return;
        }
        validated_self[i] = static_cast<float>(self);
        validated_other[i] = static_cast<float>(other);
        if (!Float_Memory_Is_Normal(&validated_self[i]) ||
            !Float_Memory_Is_Zero_Or_Normal(&validated_other[i]) ||
            (other != 0.0 && validated_other[i] == 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat,
                "GENERALIZED_BORN_INFORMATION::Initial",
                "Reason:\n\tGB atom %d effective radii are not representable "
                "as finite normal floats\n",
                i);
            return;
        }
    }

    atom_numbers = validated_atom_numbers;
    Malloc();
    for (int i = 0; i < atom_numbers; i++)
    {
        h_GB_self_radius[i] = validated_self[i];
        h_GB_other_radius[i] = validated_other[i];
    }
    deviceMemcpy(d_GB_self_radius, h_GB_self_radius,
                 sizeof(float) * atom_numbers, deviceMemcpyHostToDevice);
    deviceMemcpy(d_GB_other_radius, h_GB_other_radius,
                 sizeof(float) * atom_numbers, deviceMemcpyHostToDevice);

    this->is_initialized = 1;

    if (is_initialized && !is_controller_printf_initialized)
    {
        controller[0].Step_Print_Initial(this->module_name, "%.2f");
        is_controller_printf_initialized = 1;
        controller[0].printf("    structure last modify date is %d\n",
                             last_modify_date);
    }
    controller[0].printf(
        "END INITIALIZING STANDARD GENERALIZED BORN INFORMATION\n\n");
}

void GENERALIZED_BORN_INFORMATION::Reset_Pair_Overlap_Error()
{
#ifndef GPU_ARCH_NAME
    deviceMemset(d_pair_overlap_error, 0, 3 * sizeof(int));
#endif
}

bool GENERALIZED_BORN_INFORMATION::Check_Pair_Overlap_Error(
    const char* error_by)
{
#ifndef GPU_ARCH_NAME
    int overlap_error[3] = {0, -1, -1};
    deviceMemcpy(overlap_error, d_pair_overlap_error, sizeof(overlap_error),
                 deviceMemcpyDeviceToHost);
    if (overlap_error[0] != PairwiseInteraction::PAIR_COMPONENT_NONE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s global atoms %d %d overlap exactly; "
            "generalized-Born pairwise descreening is undefined at zero "
            "distance\n",
            module_name, overlap_error[1], overlap_error[2]);
        return true;
    }
#else
    (void)error_by;
#endif
    return false;
}

void GENERALIZED_BORN_INFORMATION::Reset_Effective_Radius_Error()
{
#ifndef GPU_ARCH_NAME
    const int no_error = -1;
    deviceMemcpy(d_effective_radius_error, &no_error, sizeof(int),
                 deviceMemcpyHostToDevice);
#endif
}

bool GENERALIZED_BORN_INFORMATION::Check_Effective_Radius_Error(
    const char* error_by)
{
#ifndef GPU_ARCH_NAME
    int invalid_atom = -1;
    deviceMemcpy(&invalid_atom, d_effective_radius_error, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (invalid_atom >= 0)
    {
        float invalid_radius = 0.0f;
        deviceMemcpy(&invalid_radius, d_GB_effective_radius + invalid_atom,
                     sizeof(float), deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s atom %d has invalid effective Born radius %.9g; "
            "the pairwise-descreening denominator must produce a finite "
            "positive normal float\n",
            module_name, invalid_atom, invalid_radius);
        return true;
    }
#else
    (void)error_by;
#endif
    return false;
}

bool GENERALIZED_BORN_INFORMATION::Get_Effective_Born_Radius(
    int atom_numbers, const VECTOR* crd)
{
    if (!is_initialized) return false;
    if (atom_numbers != this->atom_numbers || crd == NULL ||
        d_GB_self_radius == NULL || d_GB_other_radius == NULL ||
        d_GB_effective_radius == NULL)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "GENERALIZED_BORN_INFORMATION::Get_Effective_Born_Radius",
            "Reason:\n\tinvalid GB effective-radius call: runtime/module "
            "atom counts %d/%d or required storage is inconsistent\n",
            atom_numbers, this->atom_numbers);
        return false;
    }
    deviceMemset(d_GB_effective_radius, 0, sizeof(float) * atom_numbers);
    Reset_Pair_Overlap_Error();

    dim3 blockSize = {
        CONTROLLER::device_warp,
        CONTROLLER::device_max_thread / CONTROLLER::device_warp};
    dim3 gridSize = {(atom_numbers + blockSize.x - 1) / blockSize.x,
                     (atom_numbers + blockSize.y - 1) / blockSize.y};
    Launch_Device_Kernel(
        Effective_Born_Radii_Factor_Device, gridSize, blockSize, 0, NULL,
        atom_numbers, crd, radii_cutoff_square, d_GB_self_radius,
        d_GB_other_radius, d_GB_effective_radius, d_pair_overlap_error);
    if (Check_Pair_Overlap_Error(
            "GENERALIZED_BORN_INFORMATION::Get_Effective_Born_Radius"))
    {
        return false;
    }

    Launch_Device_Kernel(
        Effective_Born_Radii_Device,
        (atom_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
        d_GB_self_radius, d_GB_effective_radius);
    Reset_Effective_Radius_Error();
    Launch_Device_Kernel(
        Validate_Effective_Born_Radii_Device,
        Positive_Int_Ceil_Div(atom_numbers, CONTROLLER::device_max_thread),
        CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
        d_GB_effective_radius, d_effective_radius_error);
    if (Check_Effective_Radius_Error(
            "GENERALIZED_BORN_INFORMATION::Get_Effective_Born_Radius"))
    {
        return false;
    }
    return true;
}

void GENERALIZED_BORN_INFORMATION::GB_Force_With_Atom_Energy(
    const int atom_numbers, const VECTOR* crd, const float* charge, VECTOR* frc,
    float* atom_energy)
{
    if (is_initialized)
    {
        if (atom_numbers != this->atom_numbers || crd == NULL ||
            charge == NULL || frc == NULL || atom_energy == NULL ||
            d_dE_da == NULL || d_GB_energy_atom == NULL)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "GENERALIZED_BORN_INFORMATION::GB_Force_With_Atom_Energy",
                "Reason:\n\tinvalid GB force call: runtime/module atom "
                "counts %d/%d or required storage is inconsistent\n",
                atom_numbers, this->atom_numbers);
            return;
        }
        // Effective radii must be derived from the coordinates consumed by
        // this force call.  A single API boundary prevents stale radii after
        // coordinate updates and lets overlap validation finish before any
        // caller-owned force or energy storage is modified.
        if (!Get_Effective_Born_Radius(atom_numbers, crd)) return;
        deviceMemset(d_dE_da, 0, sizeof(float) * atom_numbers);
        deviceMemset(d_GB_energy_atom, 0, sizeof(float) * atom_numbers);

        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        dim3 gridSize = {(atom_numbers + blockSize.x - 1) / blockSize.x,
                         (atom_numbers + blockSize.y - 1) / blockSize.y};

        Launch_Device_Kernel(
            GB_inej_Force_Energy_Device, gridSize, blockSize, 0, NULL,
            atom_numbers, crd, charge, d_GB_effective_radius,
            1.0 / relative_dielectric_constant - 1.0, frc, atom_energy,
            d_dE_da, d_GB_energy_atom);

        Launch_Device_Kernel(
            GB_ieqj_Force_Energy_Device,
            (atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers, charge,
            d_GB_effective_radius, 0.5 / relative_dielectric_constant - 0.5,
            atom_energy, d_dE_da, d_GB_energy_atom);

        Launch_Device_Kernel(
            GB_accumulate_Force_Energy_Device, gridSize, blockSize, 0, NULL,
            atom_numbers, crd, radii_cutoff_square, d_GB_self_radius,
            d_GB_other_radius, d_GB_effective_radius, d_dE_da, frc);
    }
}

void GENERALIZED_BORN_INFORMATION::Clear()
{
    free(h_GB_energy_atom);
    free(h_GB_self_radius);
    free(h_GB_other_radius);
    h_GB_energy_atom = NULL;
    h_GB_self_radius = NULL;
    h_GB_other_radius = NULL;
    Free_Single_Device_Pointer((void**)&d_GB_energy_atom);
    Free_Single_Device_Pointer((void**)&d_GB_energy_sum);
    Free_Single_Device_Pointer((void**)&d_GB_self_radius);
    Free_Single_Device_Pointer((void**)&d_GB_other_radius);
    Free_Single_Device_Pointer((void**)&d_GB_effective_radius);
    Free_Single_Device_Pointer((void**)&d_dE_da);
    Free_Single_Device_Pointer((void**)&d_pair_overlap_error);
    Free_Single_Device_Pointer((void**)&d_effective_radius_error);
    atom_numbers = 0;
    h_GB_energy_sum = 0.0f;
    relative_dielectric_constant = 78.5f;
    radii_offset = 0.09f;
    radii_cutoff = 25.0f;
    radii_cutoff_square = 625.0f;
    is_initialized = 0;
    is_controller_printf_initialized = 0;
    controller = NULL;
}

void GENERALIZED_BORN_INFORMATION::Step_Print(CONTROLLER* controller)
{
    if (is_initialized)
    {
        Sum_Of_List(d_GB_energy_atom, d_GB_energy_sum, atom_numbers);
        deviceMemcpy(&h_GB_energy_sum, d_GB_energy_sum, sizeof(float),
                     deviceMemcpyDeviceToHost);
        controller->Step_Print(this->module_name, h_GB_energy_sum, true);
    }
}
