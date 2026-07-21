#include "LJ_soft_core.h"

#include <cctype>
#include <cerrno>
#include <cstdint>

#include "../xponge/load/native/lj_soft.hpp"
#include "../xponge/xponge.h"
#include "pair_activity.h"

__global__ void Copy_LJ_Type_And_Mask_To_New_Crd(const int atom_numbers,
                                                 VECTOR_LJ_SOFT_TYPE* new_crd,
                                                 const int* LJ_type_A,
                                                 const int* LJ_type_B,
                                                 const int* mask)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
#endif
    {
        new_crd[atom_i].LJ_type = LJ_type_A[atom_i];
        new_crd[atom_i].LJ_type_B = LJ_type_B[atom_i];
        new_crd[atom_i].mask = mask[atom_i];
        new_crd[atom_i].global_atom = atom_i;
    }
}

static __global__ void device_add(float* variable, const float adder)
{
    variable[0] += adder;
}

__global__ void Copy_Crd_And_Charge_To_New_Crd(const int atom_numbers,
                                               const VECTOR* crd,
                                               VECTOR_LJ_SOFT_TYPE* new_crd,
                                               const float* charge)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
#endif
    {
        new_crd[atom_i].crd = crd[atom_i];
        new_crd[atom_i].charge = charge[atom_i];
        new_crd[atom_i].charge_A = charge[atom_i];
        new_crd[atom_i].charge_B = charge[atom_i];
        new_crd[atom_i].charge_BA = 0.0f;
    }
}

static __device__ __forceinline__ void Record_Charge_Endpoint_Error(
    int* endpoint_error, const int code, const int atom_i)
{
    if (endpoint_error == NULL) return;
#ifdef GPU_ARCH_NAME
    if (atomicCAS(endpoint_error,
                  PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE,
                  code) == PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE)
    {
        endpoint_error[1] = atom_i;
    }
#elif defined(__GNUC__) || defined(__clang__)
    int expected = PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE;
    if (__atomic_compare_exchange_n(endpoint_error, &expected, code, false,
                                    __ATOMIC_RELAXED, __ATOMIC_RELAXED))
    {
        endpoint_error[1] = atom_i;
    }
#else
#pragma omp critical(sponge_charge_endpoint_error)
    {
        if (endpoint_error[0] ==
            PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE)
        {
            endpoint_error[1] = atom_i;
            endpoint_error[0] = code;
        }
    }
#endif
}

__global__ void Copy_Crd_And_Charge_To_New_Crd(
    const int atom_numbers, const VECTOR* crd, VECTOR_LJ_SOFT_TYPE* new_crd,
    const float* charge, const float lambda, int* endpoint_error)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
#endif
    {
        const float this_charge_A = new_crd[atom_i].charge_A;
        const float this_charge_B = new_crd[atom_i].charge_B;
        const PairwiseInteraction::Charge_Endpoint_Validation validation =
            PairwiseInteraction::Validate_Charge_Endpoints(
                charge[atom_i], this_charge_A, this_charge_B, NULL, lambda);
        if (validation.error != PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE)
        {
            Record_Charge_Endpoint_Error(endpoint_error, validation.error,
                                         atom_i);
        }
        new_crd[atom_i].crd = crd[atom_i];
        new_crd[atom_i].charge = validation.current;
        new_crd[atom_i].charge_A = this_charge_A;
        new_crd[atom_i].charge_B = this_charge_B;
        new_crd[atom_i].charge_BA = validation.derivative;
    }
}
__global__ void Copy_Crd_To_New_Crd(const int atom_numbers, const VECTOR* crd,
                                    VECTOR_LJ_SOFT_TYPE* new_crd)
{
#ifdef USE_GPU
    int atom_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
#endif
    {
        new_crd[atom_i].crd = crd[atom_i];
    }
}

static __global__ void Total_C6_Get(int atom_numbers, int* atom_lj_type_A,
                                    int* atom_lj_type_B, float* d_lj_Ab,
                                    float* d_lj_Bb, double* d_factor,
                                    const float lambda)
{
    double temp_sum = 0.0;
    const float lambda_ = 1.0f - lambda;
#ifdef USE_GPU
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < atom_numbers;
         i += gridDim.x * blockDim.x)
#else
#pragma omp parallel for reduction(+ : temp_sum)
    for (int i = 0; i < atom_numbers; i++)
#endif
    {
        const int itype_A = atom_lj_type_A[i];
        const int itype_B = atom_lj_type_B[i];
        double temp_small_sum = 0.0;
#ifdef USE_GPU
        for (int j = blockIdx.y * blockDim.y + threadIdx.y; j < atom_numbers;
             j += gridDim.y * blockDim.y)
#else
        for (int j = 0; j < atom_numbers; j++)
#endif
        {
            const int atom_pair_LJ_type_A =
                Get_LJ_Type(itype_A, atom_lj_type_A[j]);
            const int atom_pair_LJ_type_B =
                Get_LJ_Type(itype_B, atom_lj_type_B[j]);

            temp_small_sum += lambda_ * d_lj_Ab[atom_pair_LJ_type_A];
            temp_small_sum += lambda * d_lj_Bb[atom_pair_LJ_type_B];
        }
        temp_sum += temp_small_sum;
    }
    atomicAdd(d_factor, temp_sum);
}

static __global__ void Total_C6_B_A_Get(int atom_numbers, int* atom_lj_type_A,
                                        int* atom_lj_type_B, float* d_lj_Ab,
                                        float* d_lj_Bb, double* d_factor)
{
    double temp_sum = 0.0;
#ifdef USE_GPU
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < atom_numbers;
         i += gridDim.x * blockDim.x)
#else
#pragma omp parallel for reduction(+ : temp_sum)
    for (int i = 0; i < atom_numbers; i++)
#endif
    {
        const int itype_A = atom_lj_type_A[i];
        const int itype_B = atom_lj_type_B[i];
        double temp_small_sum = 0.0;
#ifdef USE_GPU
        for (int j = blockIdx.y * blockDim.y + threadIdx.y; j < atom_numbers;
             j += gridDim.y * blockDim.y)
#else
        for (int j = 0; j < atom_numbers; j++)
#endif
        {
            const int atom_pair_LJ_type_A =
                Get_LJ_Type(itype_A, atom_lj_type_A[j]);
            const int atom_pair_LJ_type_B =
                Get_LJ_Type(itype_B, atom_lj_type_B[j]);

            temp_small_sum +=
                d_lj_Bb[atom_pair_LJ_type_B] - d_lj_Ab[atom_pair_LJ_type_A];
        }
        temp_sum += temp_small_sum;
    }
    atomicAdd(d_factor, temp_sum);
}

template <bool need_force, bool need_energy, bool need_virial,
          bool need_coulomb, bool need_du_dlambda>
static __global__ void Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA(
    const int atom_numbers, const int solvent_numbers, const ATOM_GROUP* nl,
    const VECTOR_LJ_SOFT_TYPE* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    const float* LJ_type_AA, const float* LJ_type_AB, const float* LJ_type_BA,
    const float* LJ_type_BB, const float cutoff, VECTOR* frc,
    const float pme_beta, float* atom_energy, LTMatrix3* atom_virial,
    float* atom_direct_cf_energy, float* atom_du_dlambda_lj,
    float* atom_du_dlambda_direct, const float lambda, const float alpha,
    const float p, const float input_sigma_6, const float input_sigma_6_min,
    float* this_energy, int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = blockDim.y * blockIdx.x + threadIdx.y;
    if (atom_i < atom_numbers - solvent_numbers)
#else
#pragma omp parallel for firstprivate(lambda)
    for (int atom_i = 0; atom_i < atom_numbers - solvent_numbers; atom_i++)
#endif
    {
        ATOM_GROUP nl_i = nl[atom_i];
        VECTOR_LJ_SOFT_TYPE r1 = crd[atom_i];
        VECTOR frc_record = {0., 0., 0.};
        LTMatrix3 virial_record = {0, 0, 0, 0, 0, 0};
        float energy_lj = 0.;
        float energy_coulomb = 0.;
        float du_dlambda_lj = 0.;
        float du_dlambda_direct = 0.;
#ifdef USE_GPU
        for (int j = threadIdx.x; j < nl_i.atom_numbers; j += blockDim.x)
#else
        for (int j = 0; j < nl_i.atom_numbers; j++)
#endif
        {
            int atom_j = nl_i.atom_serial[j];
            float ij_factor = atom_j < atom_numbers ? 1.0f : 0.5f;
            VECTOR_LJ_SOFT_TYPE r2 = crd[atom_j];
            VECTOR dr = Get_Periodic_Displacement(r2, r1, cell, rcell);
            float dr_abs = norm3df(dr.x, dr.y, dr.z);
            if (dr_abs < cutoff)
            {
                int atom_pair_LJ_type_A = Get_LJ_Type(r1.LJ_type, r2.LJ_type);
                int atom_pair_LJ_type_B =
                    Get_LJ_Type(r1.LJ_type_B, r2.LJ_type_B);
                float AA = LJ_type_AA[atom_pair_LJ_type_A];
                float AB = LJ_type_AB[atom_pair_LJ_type_A];
                float BA = LJ_type_BA[atom_pair_LJ_type_B];
                float BB = LJ_type_BB[atom_pair_LJ_type_B];
                const LJ_SOFT_CORE_PAIR_RESULT pair =
                    Evaluate_LJ_Soft_Core_Pair(
                        r1, r2, dr_abs, AA, AB, BA, BB, pme_beta, lambda, alpha,
                        p, input_sigma_6, input_sigma_6_min, need_force,
                        need_energy, need_coulomb, need_du_dlambda);
                if (!pair.any_interaction)
                {
                    continue;
                }
                if (pair.singular_components !=
                    PairwiseInteraction::PAIR_COMPONENT_NONE)
                {
                    PairwiseInteraction::Fail_Exact_Overlap(
                        r1.global_atom, r2.global_atom,
                        pair.singular_components, pair_overlap_error);
                    continue;
                }
                if (need_force)
                {
                    const VECTOR frc_lin = pair.force * dr;
                    frc_record = frc_record + frc_lin;
                    if (atom_j < atom_numbers)
                    {
                        atomicAdd(frc + atom_j, -frc_lin);
                    }
                    if (need_virial)
                    {
                        virial_record =
                            virial_record -
                            ij_factor * Get_Virial_From_Force_Dis(frc_lin, dr);
                    }
                }
                if (need_energy)
                {
                    energy_lj += ij_factor * pair.lj_energy;
                    if (need_coulomb)
                    {
                        energy_coulomb += ij_factor * pair.coulomb_energy;
                    }
                }
                if (need_du_dlambda)
                {
                    // Apply ownership exactly once to the complete pair
                    // derivative, including both endpoint and radial terms.
                    du_dlambda_lj += ij_factor * pair.du_dlambda_lj;
                    if (need_coulomb)
                    {
                        du_dlambda_direct +=
                            ij_factor * pair.du_dlambda_coulomb;
                    }
                }
            }
        }
        if (need_force)
        {
            Warp_Sum_To(frc + atom_i, frc_record, warpSize);
        }
        if (need_energy)
        {
            float energy_total = energy_lj;
            if (need_coulomb)
            {
                energy_total += energy_coulomb;
            }
            Warp_Sum_To(atom_energy + atom_i, energy_total, warpSize);
            Warp_Sum_To(this_energy + atom_i, energy_lj, warpSize);
        }
        if (need_coulomb && need_energy)
        {
            Warp_Sum_To(atom_direct_cf_energy + atom_i, energy_coulomb,
                        warpSize);
        }
        if (need_virial)
        {
            Warp_Sum_To(atom_virial + atom_i, virial_record, warpSize);
        }
        if (need_du_dlambda)
        {
            Warp_Sum_To(atom_du_dlambda_lj, du_dlambda_lj, warpSize);
            if (need_coulomb)
            {
                Warp_Sum_To(atom_du_dlambda_direct, du_dlambda_direct,
                            warpSize);
            }
        }
    }
}

static bool Parse_Finite_Soft_Core_Command(CONTROLLER* controller,
                                           const char* command,
                                           const char* error_by, float* value)
{
    const char* text = controller->Command(command);
    char* end = NULL;
    errno = 0;
    const float parsed = strtof(text, &end);
    while (end != NULL && *end != '\0' &&
           isspace(static_cast<unsigned char>(*end)))
    {
        ++end;
    }
    if (!Xponge::Native_Core_Is_Strict_Decimal(text) || text == end ||
        end == NULL || *end != '\0' || errno == ERANGE ||
        !Float_Memory_Is_Finite(&parsed))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand, error_by,
            "Reason:\n\t'%s' must be one finite decimal value, received '%s'\n",
            command, text);
        return false;
    }
    *value = parsed;
    return true;
}

static bool Validate_Soft_Core_Parameters(CONTROLLER* controller,
                                          const char* error_by,
                                          const float lambda, const float alpha,
                                          const float p, const float sigma,
                                          const float sigma_min,
                                          const float cutoff)
{
    const double sigma_6 = pow(static_cast<double>(sigma), 6.0);
    const double sigma_min_6 = pow(static_cast<double>(sigma_min), 6.0);
    const double alpha_p = static_cast<double>(alpha) * p;
    const float stored_sigma_6 = static_cast<float>(sigma_6);
    const float stored_sigma_min_6 = static_cast<float>(sigma_min_6);
    const float stored_alpha_p = static_cast<float>(alpha_p);
    const double lambda_shift =
        static_cast<double>(alpha) *
        pow(static_cast<double>(lambda), static_cast<double>(p)) * sigma_6;
    const double lambda_complement = 1.0 - static_cast<double>(lambda);
    const double complement_shift =
        static_cast<double>(alpha) *
        pow(lambda_complement, static_cast<double>(p)) * sigma_6;
    const float stored_lambda_shift = static_cast<float>(lambda_shift);
    const float stored_complement_shift = static_cast<float>(complement_shift);
    if (!Float_Memory_Is_Finite(&lambda) || lambda < 0.0f || lambda > 1.0f ||
        !Float_Memory_Is_Zero_Or_Normal(&lambda) ||
        !Float_Memory_Is_Finite(&alpha) || alpha <= 0.0f ||
        !Float_Memory_Is_Normal(&alpha) || !Float_Memory_Is_Finite(&p) ||
        p < 1.0f || !Float_Memory_Is_Normal(&p) ||
        !Float_Memory_Is_Finite(&sigma) || sigma <= 0.0f ||
        !Float_Memory_Is_Normal(&sigma) ||
        !Float_Memory_Is_Finite(&sigma_min) || sigma_min < 0.0f ||
        !Float_Memory_Is_Zero_Or_Normal(&sigma_min) ||
        !Float_Memory_Is_Finite(&cutoff) || cutoff <= 0.0f ||
        !Float_Memory_Is_Normal(&cutoff) ||
        !Double_Memory_Is_Finite(&sigma_6) || sigma_6 > FLT_MAX ||
        !Float_Memory_Is_Normal(&stored_sigma_6) ||
        !Double_Memory_Is_Finite(&sigma_min_6) || sigma_min_6 > FLT_MAX ||
        !Float_Memory_Is_Zero_Or_Normal(&stored_sigma_min_6) ||
        (sigma_min != 0.0f && stored_sigma_min_6 == 0.0f) ||
        !Double_Memory_Is_Finite(&alpha_p) || alpha_p > FLT_MAX ||
        !Float_Memory_Is_Normal(&stored_alpha_p) ||
        !Double_Memory_Is_Finite(&lambda_shift) ||
        !Float_Memory_Is_Finite(&stored_lambda_shift) ||
        !Float_Memory_Is_Zero_Or_Normal(&stored_lambda_shift) ||
        (lambda != 0.0f && stored_lambda_shift == 0.0f) ||
        !Double_Memory_Is_Finite(&complement_shift) ||
        !Float_Memory_Is_Finite(&stored_complement_shift) ||
        !Float_Memory_Is_Zero_Or_Normal(&stored_complement_shift) ||
        (lambda != 1.0f && stored_complement_shift == 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand, error_by,
            "Reason:\n\tsoft-core parameters must satisfy finite lambda in "
            "[0,1], alpha > 0, p >= 1, sigma > 0, sigma_min >= 0, cutoff > "
            "0, with every nonzero value, sigma^6, alpha*p, and endpoint "
            "softening shift representable as a normal float on FTZ "
            "backends; received "
            "lambda=%g alpha=%g p=%g sigma=%g sigma_min=%g cutoff=%g\n",
            lambda, alpha, p, sigma, sigma_min, cutoff);
        return false;
    }
    return true;
}

static bool Store_Soft_Core_Long_Range_Factor(CONTROLLER* controller,
                                              const char* factor_name,
                                              const double c6_sum,
                                              const double scale,
                                              float* destination)
{
    const double scaled = c6_sum * scale;
    const float stored = static_cast<float>(scaled);
    if (!Double_Memory_Is_Finite(&c6_sum) || !Double_Memory_Is_Finite(&scale) ||
        !Double_Memory_Is_Finite(&scaled) || fabs(scaled) > FLT_MAX ||
        !Float_Memory_Is_Finite(&stored) ||
        !Float_Memory_Is_Zero_Or_Normal(&stored) ||
        (scaled != 0.0 && stored == 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
            "Reason:\n\tsoft-core %s long-range aggregate is not "
            "representable as a finite zero or normal float: C6 sum=%g, "
            "scale=%g, scaled=%g\n",
            factor_name, c6_sum, scale, scaled);
        return false;
    }
    *destination = stored;
    return true;
}

void LJ_SOFT_CORE::Initial(CONTROLLER* controller, float cutoff,
                           char* module_name)
{
    Clear();
    this->controller = controller;
    if (module_name == NULL)
    {
        strcpy(this->module_name, "LJ_soft_core");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }
    controller->printf(
        "START INITIALIZING FEP SOFT CORE FOR LJ AND COULOMB:\n");
    const auto& lj_soft = Xponge::system.classical_force_field.lj_soft_core;
    Xponge::LJSoftCore local_lj_soft;
    const Xponge::LJSoftCore* lj_soft_to_use = NULL;
    if (module_name == NULL)
    {
        lj_soft_to_use = &lj_soft;
    }
    else if (controller->Command_Exist(this->module_name, "in_file"))
    {
        Xponge::Native_Load_LJ_Soft_Core(&local_lj_soft, controller,
                                         this->module_name);
        lj_soft_to_use = &local_lj_soft;
    }

    if (lj_soft_to_use != NULL && lj_soft_to_use->atom_numbers > 0)
    {
        if (controller->Command_Exist("lambda_lj"))
        {
            if (!Parse_Finite_Soft_Core_Command(controller, "lambda_lj",
                                                "LJ_SOFT_CORE::Initial",
                                                &this->lambda))
            {
                return;
            }
            controller->printf("    FEP lj lambda: %f\n", this->lambda);
        }
        else
        {
            char error_reason[CHAR_LENGTH_MAX];
            sprintf(error_reason,
                    "Reason:\n\t'lambda_lj' is required for the calculation of "
                    "LJ_soft_core\n");
            controller->Throw_SPONGE_Error(spongeErrorMissingCommand,
                                           "LJ_SOFT_CORE::Initial",
                                           error_reason);
        }

        if (controller->Command_Exist("soft_core_alpha"))
        {
            if (!Parse_Finite_Soft_Core_Command(controller, "soft_core_alpha",
                                                "LJ_SOFT_CORE::Initial",
                                                &this->alpha))
            {
                return;
            }
            controller->printf("    FEP soft core alpha: %f\n", this->alpha);
        }
        else
        {
            controller->printf(
                "    FEP soft core alpha is set to default value 0.5\n");
            this->alpha = 0.5;
        }

        if (controller->Command_Exist("soft_core_powfer"))
        {
            if (!Parse_Finite_Soft_Core_Command(controller, "soft_core_powfer",
                                                "LJ_SOFT_CORE::Initial",
                                                &this->p))
            {
                return;
            }
            controller->printf("    FEP soft core powfer: %f\n", this->p);
        }
        else
        {
            controller->printf(
                "    FEP soft core powfer is set to default value 1.0.\n");
            this->p = 1.0;
        }

        if (controller->Command_Exist("soft_core_sigma"))
        {
            if (!Parse_Finite_Soft_Core_Command(controller, "soft_core_sigma",
                                                "LJ_SOFT_CORE::Initial",
                                                &this->sigma))
            {
                return;
            }
            controller->printf("    FEP soft core sigma: %f\n", this->sigma);
        }
        else
        {
            controller->printf(
                "    FEP soft core sigma is set to default value 3.0\n");
            this->sigma = 3.0;
        }
        if (controller->Command_Exist("soft_core_sigma_min"))
        {
            if (!Parse_Finite_Soft_Core_Command(
                    controller, "soft_core_sigma_min", "LJ_SOFT_CORE::Initial",
                    &this->sigma_min))
            {
                return;
            }
            controller->printf("    FEP soft core sigma min: %f\n",
                               this->sigma_min);
        }
        else
        {
            controller->printf(
                "    FEP soft core sigma min is set to default value 0.0\n");
            this->sigma_min = 0.0;
        }

        if (!Validate_Soft_Core_Parameters(controller, "LJ_SOFT_CORE::Initial",
                                           lambda, alpha, p, sigma, sigma_min,
                                           cutoff))
        {
            return;
        }

        atom_numbers = lj_soft_to_use->atom_numbers;
        atom_type_numbers_A = lj_soft_to_use->atom_type_numbers_A;
        atom_type_numbers_B = lj_soft_to_use->atom_type_numbers_B;
        constexpr int max_lj_type_count = 65535;
        if (atom_type_numbers_A <= 0 ||
            atom_type_numbers_A > max_lj_type_count ||
            atom_type_numbers_B <= 0 || atom_type_numbers_B > max_lj_type_count)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                "Reason:\n\tsoft-core LJ endpoint type counts %d/%d are "
                "outside the representable range [1, %d]\n",
                atom_type_numbers_A, atom_type_numbers_B, max_lj_type_count);
            return;
        }
        const std::uint64_t type_count_A =
            static_cast<std::uint64_t>(atom_type_numbers_A);
        const std::uint64_t type_count_B =
            static_cast<std::uint64_t>(atom_type_numbers_B);
        const std::uint64_t expected_pair_count_A =
            type_count_A * (type_count_A + 1ULL) / 2ULL;
        const std::uint64_t expected_pair_count_B =
            type_count_B * (type_count_B + 1ULL) / 2ULL;
        if (expected_pair_count_A >
                static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
            expected_pair_count_B >
                static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
            lj_soft_to_use->LJ_AA.size() != expected_pair_count_A ||
            lj_soft_to_use->LJ_AB.size() != expected_pair_count_A ||
            lj_soft_to_use->LJ_BA.size() != expected_pair_count_B ||
            lj_soft_to_use->LJ_BB.size() != expected_pair_count_B)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                "Reason:\n\tsoft-core LJ table shape is inconsistent: "
                "types A/B=%d/%d, expected pairs A/B=%llu/%llu, table "
                "sizes AA/AB/BA/BB=%llu/%llu/%llu/%llu\n",
                atom_type_numbers_A, atom_type_numbers_B,
                static_cast<unsigned long long>(expected_pair_count_A),
                static_cast<unsigned long long>(expected_pair_count_B),
                static_cast<unsigned long long>(lj_soft_to_use->LJ_AA.size()),
                static_cast<unsigned long long>(lj_soft_to_use->LJ_AB.size()),
                static_cast<unsigned long long>(lj_soft_to_use->LJ_BA.size()),
                static_cast<unsigned long long>(lj_soft_to_use->LJ_BB.size()));
            return;
        }
        if (lj_soft_to_use->atom_LJ_type_A.size() !=
                static_cast<std::size_t>(atom_numbers) ||
            lj_soft_to_use->atom_LJ_type_B.size() !=
                static_cast<std::size_t>(atom_numbers) ||
            (!lj_soft_to_use->subsystem_division.empty() &&
             lj_soft_to_use->subsystem_division.size() !=
                 static_cast<std::size_t>(atom_numbers)))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                "Reason:\n\tsoft-core atom metadata shape is inconsistent: "
                "atoms=%d, type maps A/B=%llu/%llu, subsystem mask=%llu\n",
                atom_numbers,
                static_cast<unsigned long long>(
                    lj_soft_to_use->atom_LJ_type_A.size()),
                static_cast<unsigned long long>(
                    lj_soft_to_use->atom_LJ_type_B.size()),
                static_cast<unsigned long long>(
                    lj_soft_to_use->subsystem_division.size()));
            return;
        }
        for (std::uint64_t pair = 0; pair < expected_pair_count_A; pair++)
        {
            if (!Float_Memory_Is_Finite(lj_soft_to_use->LJ_AA.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_soft_to_use->LJ_AA.data() +
                                                pair) ||
                !Float_Memory_Is_Finite(lj_soft_to_use->LJ_AB.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_soft_to_use->LJ_AB.data() +
                                                pair))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                    "Reason:\n\tsoft-core endpoint A LJ pair %llu has a "
                    "non-finite or subnormal coefficient\n",
                    static_cast<unsigned long long>(pair));
                return;
            }
        }
        for (std::uint64_t pair = 0; pair < expected_pair_count_B; pair++)
        {
            if (!Float_Memory_Is_Finite(lj_soft_to_use->LJ_BA.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_soft_to_use->LJ_BA.data() +
                                                pair) ||
                !Float_Memory_Is_Finite(lj_soft_to_use->LJ_BB.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_soft_to_use->LJ_BB.data() +
                                                pair))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                    "Reason:\n\tsoft-core endpoint B LJ pair %llu has a "
                    "non-finite or subnormal coefficient\n",
                    static_cast<unsigned long long>(pair));
                return;
            }
        }
        for (int atom = 0; atom < atom_numbers; atom++)
        {
            const int type_A = lj_soft_to_use->atom_LJ_type_A[atom];
            const int type_B = lj_soft_to_use->atom_LJ_type_B[atom];
            if (type_A < 0 || type_A >= atom_type_numbers_A || type_B < 0 ||
                type_B >= atom_type_numbers_B)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                    "Reason:\n\tsoft-core atom %d has endpoint LJ types "
                    "%d/%d outside [0, %d)/[0, %d)\n",
                    atom, type_A, type_B, atom_type_numbers_A,
                    atom_type_numbers_B);
                return;
            }
        }
        const bool has_endpoint_a = !lj_soft_to_use->charge_A.empty();
        const bool has_endpoint_b = !lj_soft_to_use->charge_B.empty();
        if (has_endpoint_a != has_endpoint_b ||
            (has_endpoint_a && (lj_soft_to_use->charge_A.size() !=
                                    static_cast<std::size_t>(atom_numbers) ||
                                lj_soft_to_use->charge_B.size() !=
                                    static_cast<std::size_t>(atom_numbers) ||
                                Xponge::system.atoms.charge.size() !=
                                    static_cast<std::size_t>(atom_numbers))))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                "Reason:\n\tsoft-core charge endpoints and the current "
                "charge array must either all be absent or contain exactly "
                "one value per soft-core atom\n");
            return;
        }
        has_charge_endpoints = has_endpoint_a ? 1 : 0;
        controller->printf("    atom_numbers is %d\n", atom_numbers);
        controller->printf(
            "    atom_LJ_type_number_A is %d, atom_LJ_type_number_B is %d\n",
            atom_type_numbers_A, atom_type_numbers_B);
        pair_type_numbers_A = static_cast<int>(expected_pair_count_A);
        pair_type_numbers_B = static_cast<int>(expected_pair_count_B);
        LJ_Soft_Core_Malloc();

        for (int i = 0; i < pair_type_numbers_A; i++)
        {
            h_LJ_AA[i] = lj_soft_to_use->LJ_AA[i];
        }
        for (int i = 0; i < pair_type_numbers_A; i++)
        {
            h_LJ_AB[i] = lj_soft_to_use->LJ_AB[i];
        }
        for (int i = 0; i < pair_type_numbers_B; ++i)
        {
            h_LJ_BA[i] = lj_soft_to_use->LJ_BA[i];
        }
        for (int i = 0; i < pair_type_numbers_B; ++i)
        {
            h_LJ_BB[i] = lj_soft_to_use->LJ_BB[i];
        }
        for (int i = 0; i < atom_numbers; i++)
        {
            h_atom_LJ_type_A[i] = lj_soft_to_use->atom_LJ_type_A[i];
            h_atom_LJ_type_B[i] = lj_soft_to_use->atom_LJ_type_B[i];
            if (has_charge_endpoints)
            {
                h_charge_A[i] = lj_soft_to_use->charge_A[i];
                h_charge_B[i] = lj_soft_to_use->charge_B[i];
                const float expected = PairwiseInteraction::Interpolate_Charge(
                    h_charge_A[i], h_charge_B[i], lambda);
                const float derivative = h_charge_B[i] - h_charge_A[i];
                if (!Float_Memory_Is_Finite(h_charge_A + i) ||
                    !Float_Memory_Is_Zero_Or_Normal(h_charge_A + i) ||
                    !Float_Memory_Is_Finite(h_charge_B + i) ||
                    !Float_Memory_Is_Zero_Or_Normal(h_charge_B + i) ||
                    !Float_Memory_Is_Finite(&expected) ||
                    !Float_Memory_Is_Zero_Or_Normal(&expected) ||
                    !Float_Memory_Is_Finite(&derivative) ||
                    !Float_Memory_Is_Zero_Or_Normal(&derivative) ||
                    Xponge::system.atoms.charge[i] != expected)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorConflictingCommand, "LJ_SOFT_CORE::Initial",
                        "Reason:\n\tsoft-core charge endpoint contract "
                        "failed for global atom %d: current=%g, qA=%g, "
                        "qB=%g, expected=%g, lambda=%g\n",
                        i, Xponge::system.atoms.charge[i], h_charge_A[i],
                        h_charge_B[i], expected, lambda);
                    return;
                }
            }
        }

        if (!lj_soft_to_use->subsystem_division.empty())
        {
            controller->printf(
                "    Start reading subsystem division information:\n");
            for (int i = 0; i < atom_numbers; i++)
            {
                h_subsys_division[i] = lj_soft_to_use->subsystem_division[i];
            }
            controller->printf(
                "    End reading subsystem division information\n\n");
        }
        else
        {
            controller->printf("    subsystem mask is set to 0 as default\n");
            for (int i = 0; i < atom_numbers; i++)
            {
                h_subsys_division[i] = 0;
            }
        }

        Parameter_Host_To_Device();
        is_initialized = 1;
        sigma_6 = static_cast<float>(pow(static_cast<double>(sigma), 6.0));
        sigma_6_min =
            static_cast<float>(pow(static_cast<double>(sigma_min), 6.0));
    }
    if (is_initialized)
    {
        this->cutoff = cutoff;
        Device_Malloc_Safely((void**)&crd_with_parameters,
                             sizeof(VECTOR_LJ_SOFT_TYPE) * atom_numbers);
        Launch_Device_Kernel(
            Copy_LJ_Type_And_Mask_To_New_Crd,
            (atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
            crd_with_parameters, d_atom_LJ_type_A, d_atom_LJ_type_B,
            d_subsys_division);
        controller->printf("    Start initializing long range LJ correction\n");
        long_range_factor = 0;
        double h_factor = 0.0;
        double* d_factor = NULL;
        Device_Malloc_Safely((void**)&d_factor, sizeof(double));
        deviceMemset(d_factor, 0, sizeof(double));

        dim3 gridSize = {4, 4};
        dim3 blockSize = {32, 32};
        Launch_Device_Kernel(Total_C6_Get, gridSize, blockSize, 0, NULL,
                             atom_numbers, d_atom_LJ_type_A, d_atom_LJ_type_B,
                             d_LJ_AB, d_LJ_BB, d_factor, this->lambda);

        deviceMemcpy(&h_factor, d_factor, sizeof(double),
                     deviceMemcpyDeviceToHost);
        const double c6_sum = h_factor;
        deviceMemset(d_factor, 0, sizeof(double));

        Launch_Device_Kernel(Total_C6_B_A_Get, gridSize, blockSize, 0, NULL,
                             atom_numbers, d_atom_LJ_type_A, d_atom_LJ_type_B,
                             d_LJ_AB, d_LJ_BB, d_factor);
        deviceMemcpy(&h_factor, d_factor, sizeof(double),
                     deviceMemcpyDeviceToHost);
        const double c6_derivative_sum = h_factor;
        Free_Single_Device_Pointer((void**)&d_factor);

        const double cutoff_double = static_cast<double>(cutoff);
        const double long_range_scale =
            -2.0 / 3.0 * static_cast<double>(CONSTANT_Pi) /
            (cutoff_double * cutoff_double * cutoff_double) / 6.0;
        if (!Store_Soft_Core_Long_Range_Factor(controller, "energy", c6_sum,
                                               long_range_scale,
                                               &long_range_factor) ||
            !Store_Soft_Core_Long_Range_Factor(
                controller, "lambda derivative", c6_derivative_sum,
                long_range_scale, &long_range_factor_TI))
        {
            is_initialized = 0;
            return;
        }
        controller->printf("        long range correction factor is: %e\n",
                           long_range_factor);
        controller->printf("    End initializing long range LJ correction\n");
    }
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial("LJ_soft", "%.2f");
        controller->Step_Print_Initial("LJ_soft_short", "%.2f");
        controller->Step_Print_Initial("LJ_soft_long", "%.2f");
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
    }
    controller->printf(
        "END INITIALIZING LENNADR JONES SOFT CORE INFORMATION\n\n");
}

void LJ_SOFT_CORE::Clear()
{
    is_initialized = 0;
    local_metadata_is_ready = false;

    const auto free_host_and_device = [](void** host, void** device)
    {
        if (*device != NULL)
        {
            Free_Host_And_Device_Pointer(host, device);
        }
        else if (*host != NULL)
        {
            free(*host);
            *host = NULL;
        }
    };
    const auto free_device_copy = [](void** device)
    {
        if (*device != NULL)
        {
            Free_Host_And_Device_Pointer(NULL, device);
        }
    };
    const auto free_device_only = [](void** device)
    {
        if (*device != NULL)
        {
            Free_Single_Device_Pointer(device);
        }
    };

    free_host_and_device((void**)&h_atom_LJ_type_A, (void**)&d_atom_LJ_type_A);
    free_host_and_device((void**)&h_atom_LJ_type_B, (void**)&d_atom_LJ_type_B);
    free_host_and_device((void**)&h_LJ_AA, (void**)&d_LJ_AA);
    free_host_and_device((void**)&h_LJ_AB, (void**)&d_LJ_AB);
    free_host_and_device((void**)&h_LJ_BA, (void**)&d_LJ_BA);
    free_host_and_device((void**)&h_LJ_BB, (void**)&d_LJ_BB);
    free_host_and_device((void**)&h_subsys_division,
                         (void**)&d_subsys_division);
    free_host_and_device((void**)&h_charge_A, (void**)&d_charge_A);
    free_host_and_device((void**)&h_charge_B, (void**)&d_charge_B);

    free_host_and_device((void**)&h_LJ_energy_atom, (void**)&d_LJ_energy_atom);
    free_device_copy((void**)&d_LJ_energy_sum);
    free_host_and_device((void**)&h_LJ_energy_atom_intersys,
                         (void**)&d_LJ_energy_atom_intersys);
    free_host_and_device((void**)&h_LJ_energy_atom_intrasys,
                         (void**)&d_LJ_energy_atom_intrasys);
    free_device_copy((void**)&d_LJ_energy_sum_intersys);
    free_device_copy((void**)&d_LJ_energy_sum_intrasys);
    free_device_copy((void**)&d_direct_ene_sum_intersys);
    free_device_copy((void**)&d_direct_ene_sum_intrasys);
    free_host_and_device((void**)&h_sigma_of_dH_dlambda_lj,
                         (void**)&d_sigma_of_dH_dlambda_lj);
    free_host_and_device((void**)&h_sigma_of_dH_dlambda_direct,
                         (void**)&d_sigma_of_dH_dlambda_direct);

    free_device_only((void**)&crd_with_parameters);
    free_device_only((void**)&crd_with_LJ_parameters_local);
    free_device_only((void**)&d_pair_overlap_error);
    free_device_only((void**)&d_local_metadata_error);
    free_device_only((void**)&d_charge_endpoint_error);

    controller = NULL;
    atom_numbers = 0;
    atom_type_numbers_A = 0;
    atom_type_numbers_B = 0;
    pair_type_numbers_A = 0;
    pair_type_numbers_B = 0;
    local_atom_numbers = 0;
    ghost_numbers = 0;
    has_charge_endpoints = 0;
    h_LJ_energy_sum = 0.0f;
    h_LJ_energy_sum_intersys = 0.0f;
    h_LJ_energy_sum_intrasys = 0.0f;
    h_direct_ene_sum = 0.0f;
    h_direct_ene_sum_intersys = 0.0f;
    h_direct_ene_sum_intrasys = 0.0f;
    lambda = 0.0f;
    alpha = 0.0f;
    p = 0.0f;
    sigma = 0.0f;
    sigma_6 = 0.0f;
    sigma_min = 0.0f;
    sigma_6_min = 0.0f;
    cutoff = 10.0f;
    h_LJ_long_energy = 0.0f;
    long_range_factor = 0.0f;
    long_range_factor_TI = 0.0f;
}

void LJ_SOFT_CORE::LJ_Soft_Core_Malloc()
{
    Malloc_Safely((void**)&h_LJ_energy_atom, sizeof(float) * atom_numbers);
    Malloc_Safely((void**)&h_atom_LJ_type_A, sizeof(int) * atom_numbers);
    Malloc_Safely((void**)&h_atom_LJ_type_B, sizeof(int) * atom_numbers);
    Malloc_Safely((void**)&h_LJ_AA, sizeof(float) * pair_type_numbers_A);
    Malloc_Safely((void**)&h_LJ_AB, sizeof(float) * pair_type_numbers_A);
    Malloc_Safely((void**)&h_LJ_BA, sizeof(float) * pair_type_numbers_B);
    Malloc_Safely((void**)&h_LJ_BB, sizeof(float) * pair_type_numbers_B);
    Malloc_Safely((void**)&h_subsys_division, sizeof(int) * atom_numbers);
    if (has_charge_endpoints)
    {
        Malloc_Safely((void**)&h_charge_A, sizeof(float) * atom_numbers);
        Malloc_Safely((void**)&h_charge_B, sizeof(float) * atom_numbers);
    }

    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_sum, &h_LJ_energy_sum,
                                  sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_atom, h_LJ_energy_atom,
                                  sizeof(float) * atom_numbers);

    Malloc_Safely((void**)&h_LJ_energy_atom_intersys,
                  sizeof(float) * atom_numbers);
    Malloc_Safely((void**)&h_LJ_energy_atom_intrasys,
                  sizeof(float) * atom_numbers);

    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_atom_intersys,
                                  h_LJ_energy_atom_intersys,
                                  sizeof(float) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_atom_intrasys,
                                  h_LJ_energy_atom_intrasys,
                                  sizeof(float) * atom_numbers);

    Device_Malloc_And_Copy_Safely((void**)&d_direct_ene_sum_intersys,
                                  &h_direct_ene_sum_intersys, sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_direct_ene_sum_intrasys,
                                  &h_direct_ene_sum_intrasys, sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_sum_intersys,
                                  &h_LJ_energy_sum_intersys, sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_sum_intrasys,
                                  &h_LJ_energy_sum_intrasys, sizeof(float));

    Malloc_Safely((void**)&h_sigma_of_dH_dlambda_lj, sizeof(float));
    Malloc_Safely((void**)&h_sigma_of_dH_dlambda_direct, sizeof(float));

    Device_Malloc_And_Copy_Safely((void**)&d_sigma_of_dH_dlambda_lj,
                                  h_sigma_of_dH_dlambda_lj, sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&d_sigma_of_dH_dlambda_direct,
                                  h_sigma_of_dH_dlambda_direct, sizeof(float));
}

void LJ_SOFT_CORE::Parameter_Host_To_Device()
{
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_AA, h_LJ_AA,
                                  sizeof(float) * pair_type_numbers_A);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_AB, h_LJ_AB,
                                  sizeof(float) * pair_type_numbers_A);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_BA, h_LJ_BA,
                                  sizeof(float) * pair_type_numbers_B);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_BB, h_LJ_BB,
                                  sizeof(float) * pair_type_numbers_B);

    Device_Malloc_And_Copy_Safely((void**)&d_atom_LJ_type_A, h_atom_LJ_type_A,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_atom_LJ_type_B, h_atom_LJ_type_B,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_subsys_division, h_subsys_division,
                                  sizeof(int) * atom_numbers);
    if (has_charge_endpoints)
    {
        Device_Malloc_And_Copy_Safely((void**)&d_charge_A, h_charge_A,
                                      sizeof(float) * atom_numbers);
        Device_Malloc_And_Copy_Safely((void**)&d_charge_B, h_charge_B,
                                      sizeof(float) * atom_numbers);
    }
    Device_Malloc_Safely((void**)&crd_with_LJ_parameters_local,
                         sizeof(VECTOR_LJ_SOFT_TYPE) * atom_numbers);
    Device_Malloc_Safely((void**)&d_local_metadata_error, sizeof(int));
    Device_Malloc_Safely((void**)&d_charge_endpoint_error, 2 * sizeof(int));
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_pair_overlap_error, 3 * sizeof(int));
#endif
}

void LJ_SOFT_CORE::LJ_Soft_Core_PME_Direct_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const int local_atom_numbers,
    const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
    const float* charge, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, const float pme_beta,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_lj_virial, float* atom_direct_pme_energy)
{
    if (is_initialized)
    {
        const char* error_by =
            "LJ_SOFT_CORE::"
            "LJ_Soft_Core_PME_Direct_Force_With_Atom_Energy_And_Virial";
        if (!Validate_Local_State(error_by, atom_numbers, local_atom_numbers,
                                  ghost_numbers, solvent_numbers))
        {
            return;
        }
        if (!Prepare_Local_Coordinates(
                error_by, this->local_atom_numbers + this->ghost_numbers, crd,
                charge))
        {
            return;
        }

        if (need_atom_energy)
        {
            deviceMemset(d_LJ_energy_atom, 0, sizeof(float) * atom_numbers);
            deviceMemset(atom_direct_pme_energy, 0,
                         sizeof(float) * atom_numbers);
        }

        if (atom_numbers == 0 || local_atom_numbers == 0) return;

        Reset_Pair_Overlap_Error();

        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        dim3 gridSize = (local_atom_numbers + blockSize.y - 1) / blockSize.y;

        auto f =
            Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA<true, false, false,
                                                            true, false>;

        if (!need_atom_energy && !need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA<
                true, false, false, true, false>;
        }
        else if (need_atom_energy && !need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA<
                true, true, false, true, false>;
        }
        else if (!need_atom_energy && need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA<
                true, false, true, true, false>;
        }
        else
        {
            f = Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA<
                true, true, true, true, false>;
        }
        Launch_Device_Kernel(
            f, gridSize, blockSize, 0, NULL, local_atom_numbers,
            solvent_numbers, nl, crd_with_LJ_parameters_local, cell, rcell,
            d_LJ_AA, d_LJ_AB, d_LJ_BA, d_LJ_BB, cutoff, frc, pme_beta,
            atom_energy, atom_lj_virial, atom_direct_pme_energy, NULL, NULL,
            lambda, alpha, p, sigma_6, sigma_6_min, d_LJ_energy_atom,
            d_pair_overlap_error);
        Check_Pair_Overlap_Error(error_by);
    }
}

float LJ_SOFT_CORE::Get_Partial_H_Partial_Lambda_With_Columb_Direct(
    const int solvent_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const float* charge, const ATOM_GROUP* nl,
    const float pme_beta)
{
    if (is_initialized)
    {
        if (!Validate_Local_State(
                "LJ_SOFT_CORE::"
                "Get_Partial_H_Partial_Lambda_With_Columb_Direct",
                atom_numbers, local_atom_numbers, ghost_numbers,
                solvent_numbers))
        {
            return NAN;
        }
        deviceMemset(d_sigma_of_dH_dlambda_lj, 0, sizeof(float));
        deviceMemset(d_sigma_of_dH_dlambda_direct, 0, sizeof(float));
        Reset_Pair_Overlap_Error();

        const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
        if (!Prepare_Local_Coordinates(
                "LJ_SOFT_CORE::"
                "Get_Partial_H_Partial_Lambda_With_Columb_Direct",
                local_coordinate_numbers, crd, charge))
        {
            return NAN;
        }

        if (local_atom_numbers > 0)
        {
            dim3 blockSize = {
                CONTROLLER::device_warp,
                CONTROLLER::device_max_thread / CONTROLLER::device_warp};
            dim3 gridSize =
                (local_atom_numbers + blockSize.y - 1) / blockSize.y;
            auto f = Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA<
                false, false, false, true, true>;

            if (!has_charge_endpoints)
            {
                f = Lennard_Jones_And_Direct_Coulomb_Soft_Core_CUDA<
                    false, false, false, false, true>;
            }
            Launch_Device_Kernel(
                f, gridSize, blockSize, 0, NULL, local_atom_numbers,
                solvent_numbers, nl, crd_with_LJ_parameters_local, cell, rcell,
                d_LJ_AA, d_LJ_AB, d_LJ_BA, d_LJ_BB, cutoff, NULL, pme_beta,
                NULL, NULL, NULL, d_sigma_of_dH_dlambda_lj,
                d_sigma_of_dH_dlambda_direct, lambda, alpha, p, sigma_6,
                sigma_6_min, NULL, d_pair_overlap_error);
            if (Check_Pair_Overlap_Error(
                    "LJ_SOFT_CORE::"
                    "Get_Partial_H_Partial_Lambda_With_Columb_Direct"))
            {
                return NAN;
            }
        }

        deviceMemcpy(h_sigma_of_dH_dlambda_lj, d_sigma_of_dH_dlambda_lj,
                     sizeof(float), deviceMemcpyDeviceToHost);
        deviceMemcpy(h_sigma_of_dH_dlambda_direct, d_sigma_of_dH_dlambda_direct,
                     sizeof(float), deviceMemcpyDeviceToHost);
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, h_sigma_of_dH_dlambda_lj, 1, MPI_FLOAT,
                      MPI_SUM, CONTROLLER::pp_comm);
        MPI_Allreduce(MPI_IN_PLACE, h_sigma_of_dH_dlambda_direct, 1, MPI_FLOAT,
                      MPI_SUM, CONTROLLER::pp_comm);
#endif
        return *h_sigma_of_dH_dlambda_lj + *h_sigma_of_dH_dlambda_direct +
               long_range_factor_TI / cell.a11 / cell.a22 / cell.a33;
    }
    else
    {
        return NAN;
    }
}

void LJ_SOFT_CORE::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized || CONTROLLER::MPI_rank >= CONTROLLER::PP_MPI_size)
        return;
    Sum_Of_List(d_LJ_energy_atom, d_LJ_energy_sum, atom_numbers);
    deviceMemcpy(&h_LJ_energy_sum, d_LJ_energy_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &h_LJ_energy_sum, 1, MPI_FLOAT, MPI_SUM,
                  CONTROLLER::pp_comm);
#endif
    controller->Step_Print("LJ_soft_short", h_LJ_energy_sum);
    controller->Step_Print("LJ_soft_long", h_LJ_long_energy);
    controller->Step_Print("LJ_soft", h_LJ_energy_sum + h_LJ_long_energy, true);
}

static __global__ void Long_Range_Virial_Correction(LTMatrix3* d_virial,
                                                    const float factor)
{
    d_virial->a11 += factor;
    d_virial->a22 += factor;
    d_virial->a33 += factor;
}

void LJ_SOFT_CORE::Long_Range_Correction(int need_pressure, LTMatrix3* d_virial,
                                         int need_potential, float* d_potential,
                                         const float volume)
{
    if (is_initialized && CONTROLLER::PP_MPI_rank == 0)
    {
        if (need_pressure > 0)
        {
            Launch_Device_Kernel(Long_Range_Virial_Correction, 1, 1, 0, NULL,
                                 d_virial, 2 * long_range_factor / volume);
        }
        if (need_potential > 0)
        {
            Launch_Device_Kernel(device_add, 1, 1, 0, NULL, d_potential,
                                 long_range_factor / volume);
            h_LJ_long_energy = long_range_factor / volume;
        }
    }
}

static __global__ void get_local_device(
    int* atom_local, int local_atom_numbers, int ghost_numbers,
    int global_atom_numbers, int* d_atom_LJ_type_A, int* d_atom_LJ_type_B,
    int* d_mask, const float* charge_A, const float* charge_B,
    const int has_charge_endpoints,
    VECTOR_LJ_SOFT_TYPE* crd_with_LJ_parameters_local, int* invalid_local_index)
{
    SIMPLE_DEVICE_FOR(i, local_atom_numbers + ghost_numbers)
    {
        int atom_i = atom_local[i];
        if (atom_i < 0 || atom_i >= global_atom_numbers)
        {
            atomicExch(invalid_local_index, i);
        }
        else
        {
            crd_with_LJ_parameters_local[i].LJ_type = d_atom_LJ_type_A[atom_i];
            crd_with_LJ_parameters_local[i].LJ_type_B =
                d_atom_LJ_type_B[atom_i];
            crd_with_LJ_parameters_local[i].mask = d_mask[atom_i];
            crd_with_LJ_parameters_local[i].global_atom = atom_i;
            if (has_charge_endpoints)
            {
                crd_with_LJ_parameters_local[i].charge_A = charge_A[atom_i];
                crd_with_LJ_parameters_local[i].charge_B = charge_B[atom_i];
                crd_with_LJ_parameters_local[i].charge_BA =
                    charge_B[atom_i] - charge_A[atom_i];
            }
        }
    }
}

void LJ_SOFT_CORE::Get_Local(int* atom_local, int local_atom_numbers,
                             int ghost_numbers)
{
    if (!is_initialized) return;
    local_metadata_is_ready = false;
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > atom_numbers ||
        ghost_numbers > atom_numbers - local_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "LJ_SOFT_CORE::Get_Local",
            "Reason:\n\t%s received invalid local/ghost atom counts %d/%d "
            "for %d global atoms\n",
            module_name, local_atom_numbers, ghost_numbers, atom_numbers);
        return;
    }
    this->local_atom_numbers = local_atom_numbers;
    this->ghost_numbers = ghost_numbers;
    const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
    if (local_coordinate_numbers == 0)
    {
        local_metadata_is_ready = true;
        return;
    }
    deviceMemset(d_local_metadata_error, -1, sizeof(int));
    Launch_Device_Kernel(
        get_local_device,
        (local_coordinate_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, atom_local, local_atom_numbers,
        ghost_numbers, atom_numbers, d_atom_LJ_type_A, d_atom_LJ_type_B,
        d_subsys_division, d_charge_A, d_charge_B, has_charge_endpoints,
        crd_with_LJ_parameters_local, d_local_metadata_error);
    int invalid_local_index = -1;
    deviceMemcpy(&invalid_local_index, d_local_metadata_error, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (invalid_local_index >= 0)
    {
        int invalid_global_atom = -1;
        deviceMemcpy(&invalid_global_atom, atom_local + invalid_local_index,
                     sizeof(int), deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "LJ_SOFT_CORE::Get_Local",
            "Reason:\n\t%s local coordinate %d maps to global atom %d "
            "outside [0, %d)\n",
            module_name, invalid_local_index, invalid_global_atom,
            atom_numbers);
        return;
    }
    local_metadata_is_ready = true;
}

bool LJ_SOFT_CORE::Validate_Local_State(const char* error_by,
                                        int global_atom_numbers,
                                        int local_atom_numbers,
                                        int ghost_numbers, int solvent_numbers)
{
    if (global_atom_numbers != atom_numbers ||
        local_atom_numbers != this->local_atom_numbers ||
        ghost_numbers != this->ghost_numbers || !local_metadata_is_ready ||
        solvent_numbers < 0 || solvent_numbers > local_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s local nonbond metadata mismatch: call has "
            "global/local/ghost/solvent counts %d/%d/%d/%d, initialized "
            "state has %d/%d/%d and ready=%d\n",
            module_name, global_atom_numbers, local_atom_numbers, ghost_numbers,
            solvent_numbers, atom_numbers, this->local_atom_numbers,
            this->ghost_numbers, static_cast<int>(local_metadata_is_ready));
        return false;
    }
    return true;
}

bool LJ_SOFT_CORE::Prepare_Local_Coordinates(const char* error_by,
                                             const int coordinate_numbers,
                                             const VECTOR* crd,
                                             const float* charge)
{
    if (coordinate_numbers < 0 || coordinate_numbers > atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s received invalid local coordinate count %d for "
            "%d global atoms\n",
            module_name, coordinate_numbers, atom_numbers);
        return false;
    }
    if (coordinate_numbers == 0)
    {
        return true;
    }
    if (crd == NULL || charge == NULL)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s received null coordinate/current-charge storage "
            "for %d local coordinates\n",
            module_name, coordinate_numbers);
        return false;
    }

    const int blocks =
        (coordinate_numbers + CONTROLLER::device_max_thread - 1) /
        CONTROLLER::device_max_thread;
    if (!has_charge_endpoints)
    {
        Launch_Device_Kernel(Copy_Crd_And_Charge_To_New_Crd, blocks,
                             CONTROLLER::device_max_thread, 0, NULL,
                             coordinate_numbers, crd,
                             crd_with_LJ_parameters_local, charge);
        return true;
    }

    deviceMemset(d_charge_endpoint_error, 0, 2 * sizeof(int));
    Launch_Device_Kernel(Copy_Crd_And_Charge_To_New_Crd, blocks,
                         CONTROLLER::device_max_thread, 0, NULL,
                         coordinate_numbers, crd, crd_with_LJ_parameters_local,
                         charge, lambda, d_charge_endpoint_error);
    int endpoint_error[2] = {PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE,
                             -1};
    deviceMemcpy(endpoint_error, d_charge_endpoint_error,
                 sizeof(endpoint_error), deviceMemcpyDeviceToHost);
    if (endpoint_error[0] == PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE)
    {
        return true;
    }

    float current = NAN;
    float endpoint_a = NAN;
    float endpoint_b = NAN;
    float derivative = NAN;
    int global_atom = -1;
    if (endpoint_error[1] < 0 || endpoint_error[1] >= coordinate_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s charge endpoint validation returned invalid "
            "local coordinate index %d for %d coordinates\n",
            module_name, endpoint_error[1], coordinate_numbers);
        return false;
    }
    else
    {
        VECTOR_LJ_SOFT_TYPE endpoint_record{};
        deviceMemcpy(&current, charge + endpoint_error[1], sizeof(float),
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(&endpoint_record,
                     crd_with_LJ_parameters_local + endpoint_error[1],
                     sizeof(endpoint_record), deviceMemcpyDeviceToHost);
        endpoint_a = endpoint_record.charge_A;
        endpoint_b = endpoint_record.charge_B;
        derivative = endpoint_record.charge_BA;
        global_atom = endpoint_record.global_atom;
    }
    if (global_atom < 0 || global_atom >= atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s local charge coordinate %d maps to invalid global "
            "atom %d outside [0, %d)\n",
            module_name, endpoint_error[1], global_atom, atom_numbers);
        return false;
    }
    const char* reason =
        endpoint_error[0] ==
                PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONFINITE
            ? "contains a non-finite charge or endpoint"
        : endpoint_error[0] ==
                PairwiseInteraction::CHARGE_ENDPOINT_ERROR_CURRENT_MISMATCH
            ? "does not equal fma(lambda, qB-qA, qA)"
            : "has qB-qA inconsistent with its explicit endpoints";
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown, error_by,
        "Reason:\n\t%s local charge coordinate %d (global atom %d) %s: "
        "current=%g, qA=%g, qB=%g, qB-qA=%g, lambda=%g\n",
        module_name, endpoint_error[1], global_atom, reason, current,
        endpoint_a, endpoint_b, derivative, lambda);
    return false;
}

void LJ_SOFT_CORE::Reset_Pair_Overlap_Error()
{
#ifndef GPU_ARCH_NAME
    deviceMemset(d_pair_overlap_error, 0, 3 * sizeof(int));
#endif
}

bool LJ_SOFT_CORE::Check_Pair_Overlap_Error(const char* error_by)
{
#ifndef GPU_ARCH_NAME
    int overlap_error[3] = {0, -1, -1};
    deviceMemcpy(overlap_error, d_pair_overlap_error, sizeof(overlap_error),
                 deviceMemcpyDeviceToHost);
    if (overlap_error[0] != PairwiseInteraction::PAIR_COMPONENT_NONE)
    {
        const char* component =
            overlap_error[0] ==
                    (PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES |
                     PairwiseInteraction::PAIR_COMPONENT_COULOMB)
                ? "LJ and Coulomb components"
            : overlap_error[0] &
                    PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES
                ? "LJ component"
                : "Coulomb component";
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s global atoms %d %d overlap exactly with active "
            "%s; an active hard nonbond pair has undefined inverse-distance "
            "energy and force\n",
            module_name, overlap_error[1], overlap_error[2], component);
        return true;
    }
#else
    (void)error_by;
#endif
    return false;
}
