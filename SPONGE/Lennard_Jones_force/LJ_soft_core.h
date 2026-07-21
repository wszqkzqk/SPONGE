#pragma once
#include "../common.h"
#include "../control.h"
#include "Lennard_Jones_force.h"
#include "pair_activity.h"

// 用于计算LJ_Force时使用的坐标和记录的原子LJ种类序号与原子电荷
#ifndef UINT_VECTOR_LJ_FEP_TYPE_DEFINE
#define UINT_VECTOR_LJ_FEP_TYPE_DEFINE

__device__ __host__ __forceinline__ float Get_Soft_Core_Sigma(
    const float A, const float B, const float input_sigma_6,
    const float input_sigma_6_min);
__device__ __host__ __forceinline__ float Get_Soft_Core_Distance(
    const float A, const float B, const float sigma, const float dr_abs,
    const float alpha, const float p, const float one_minus_lambda);
__device__ __host__ __forceinline__ float Get_Soft_Core_dU_dlambda(
    const float F, const float sigma, const float dr_soft_core,
    const float alpha, const float p, const float lambda);
struct VECTOR_LJ_SOFT_TYPE
{
    VECTOR crd;
    int LJ_type;
    int LJ_type_B;
    int mask;
    float charge;
    float charge_A;
    float charge_B;
    float charge_BA;
    int global_atom;
    friend __device__ __host__ __forceinline__ VECTOR Get_Periodic_Displacement(
        VECTOR_LJ_SOFT_TYPE vec_a, VECTOR_LJ_SOFT_TYPE vec_b, LTMatrix3 cell,
        LTMatrix3 rcell)
    {
        return Get_Periodic_Displacement(vec_a.crd, vec_b.crd, cell, rcell);
    }
    friend __device__ __host__ __forceinline__ float Get_LJ_Energy(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float A, const float B)
    {
        float dr_6 = powf(dr_abs, -6.0f);
        return (0.083333333f * A * dr_6 - 0.166666667f * B) * dr_6;
    }
    friend __device__ __host__ __forceinline__ float Get_LJ_Force(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float A, const float B)
    {
        return (B - A * powf(dr_abs, -6.0f)) * powf(dr_abs, -8.0f);
    }
    friend __device__ __host__ __forceinline__ float Get_LJ_Virial(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float A, const float B)
    {
        float dr_6 = powf(dr_abs, -6.0f);
        return -(B - A * dr_6) * dr_6;
    }
    friend __device__ __host__ __forceinline__ float Get_Direct_Coulomb_Energy(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float pme_beta)
    {
        const double coefficient = static_cast<double>(r1.charge) * r2.charge;
        return static_cast<float>(coefficient * erfc(pme_beta * dr_abs) /
                                  dr_abs);
    }
    friend __device__ __host__ __forceinline__ float Get_Direct_Coulomb_Force(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float pme_beta)
    {
        const double beta_dr = static_cast<double>(pme_beta) * dr_abs;
        const double coefficient = static_cast<double>(r1.charge) * r2.charge;
        return static_cast<float>(
            coefficient * pow(static_cast<double>(dr_abs), -3.0) *
            (beta_dr * TWO_DIVIDED_BY_SQRT_PI * exp(-beta_dr * beta_dr) +
             erfc(beta_dr)));
    }
    friend __device__ __host__ __forceinline__ float Get_Direct_Coulomb_Virial(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float pme_beta)
    {
        const double beta_dr = static_cast<double>(pme_beta) * dr_abs;
        const double coefficient = static_cast<double>(r1.charge) * r2.charge;
        return static_cast<float>(
            coefficient / dr_abs *
            (beta_dr * TWO_DIVIDED_BY_SQRT_PI * exp(-beta_dr * beta_dr) +
             erfc(beta_dr)));
    }
    friend __device__ __host__ __forceinline__ float
    Get_Direct_Coulomb_dU_dlambda(VECTOR_LJ_SOFT_TYPE r1,
                                  VECTOR_LJ_SOFT_TYPE r2, const float dr_abs,
                                  const float pme_beta)
    {
        const double derivative =
            static_cast<double>(r2.charge_BA) * r1.charge +
            static_cast<double>(r2.charge) * r1.charge_BA;
        return static_cast<float>(erfc(pme_beta * dr_abs) * derivative /
                                  dr_abs);
    }
    friend __device__ __host__ __forceinline__ float Get_Soft_Core_Sigma(
        const float A, const float B, const float input_sigma_6,
        const float input_sigma_6_min)
    {
        return (A == 0 || B == 0) ? input_sigma_6
                                  : fmaxf(0.5f * A / B, input_sigma_6_min);
    }
    friend __device__ __host__ __forceinline__ float Get_Soft_Core_Distance(
        const float A, const float B, const float sigma, const float dr_abs,
        const float alpha, const float p, const float one_minus_lambda)
    {
        float alpha_lambda_p = alpha * powf(one_minus_lambda, p);
        float dr6 = powf(dr_abs, 6.0f);
        return powf(dr6 + alpha_lambda_p * sigma, 0.16666667f);
    }
    friend __device__ __host__ __forceinline__ float Get_Soft_Core_LJ_Force(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float dr_soft_core, const float A, const float B)
    {
        return powf(dr_abs, 4.0f) * (B - A * powf(dr_soft_core, -6.0f)) *
               powf(dr_soft_core, -12.0f);
    }
    friend __device__ __host__ __forceinline__ float Get_Soft_Core_LJ_Virial(
        VECTOR_LJ_SOFT_TYPE r1, VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
        const float dr_soft_core, const float A, const float B)
    {
        return -powf(dr_abs, 6.0f) * (B - A * powf(dr_soft_core, -6.0f)) *
               powf(dr_soft_core, -12.0f);
    }

    friend __device__ __host__ __forceinline__ float
    Get_Soft_Core_Direct_Coulomb_Force(VECTOR_LJ_SOFT_TYPE r1,
                                       VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
                                       const float dr_soft_core,
                                       const float pme_beta)
    {
        const double beta_dr_soft_core =
            static_cast<double>(dr_soft_core) * pme_beta;
        const double coefficient = static_cast<double>(r1.charge) * r2.charge;
        return static_cast<float>(coefficient *
                                  pow(static_cast<double>(dr_abs), 4.0) *
                                  (exp(-beta_dr_soft_core * beta_dr_soft_core) *
                                       TWO_DIVIDED_BY_SQRT_PI * pme_beta +
                                   erfc(beta_dr_soft_core) / dr_soft_core) *
                                  pow(static_cast<double>(dr_soft_core), -6.0));
    }
    friend __device__ __host__ __forceinline__ float
    Get_Soft_Core_Direct_Coulomb_Virial(VECTOR_LJ_SOFT_TYPE r1,
                                        VECTOR_LJ_SOFT_TYPE r2, float dr_abs,
                                        const float dr_soft_core,
                                        const float pme_beta)
    {
        const double beta_dr_soft_core =
            static_cast<double>(dr_soft_core) * pme_beta;
        const double coefficient = static_cast<double>(r1.charge) * r2.charge;
        return static_cast<float>(coefficient *
                                  pow(static_cast<double>(dr_abs), 6.0) *
                                  (exp(-beta_dr_soft_core * beta_dr_soft_core) *
                                       TWO_DIVIDED_BY_SQRT_PI * pme_beta +
                                   erfc(beta_dr_soft_core) / dr_soft_core) *
                                  pow(static_cast<double>(dr_soft_core), -6.0));
    }
    friend __device__ __host__ __forceinline__ float Get_Soft_Core_dU_dlambda(
        const float F, const float sigma, const float dr_soft_core,
        const float alpha, const float p, const float lambda)
    {
        const float lambda_power_derivative =
            lambda == 0.0f ? (p == 1.0f ? 1.0f : 0.0f) : powf(lambda, p - 1.0f);
        return 0.16666667f * p * alpha * (1 - lambda) *
               lambda_power_derivative * F * powf(dr_soft_core, -4.0f) * sigma;
    }
};

struct LJ_SOFT_CORE_PAIR_RESULT
{
    float force = 0.0f;
    float lj_energy = 0.0f;
    float coulomb_energy = 0.0f;
    float du_dlambda_lj = 0.0f;
    float du_dlambda_coulomb = 0.0f;
    int singular_components = PairwiseInteraction::PAIR_COMPONENT_NONE;
    bool any_interaction = false;
};

// Evaluate one pair under the Hamiltonian used by both production dynamics and
// thermodynamic integration:
//
//   U_LJ = (1-lambda) U_A(r_A) + lambda U_B(r_B)
//   U_C  = q_i(lambda) q_j(lambda)
//          [(1-lambda) K(r_A) + lambda K(r_B)]
//
// The Coulomb bracket is replaced by K(r) when both endpoint pair
// coefficients are active.  Explicit endpoint charges decide that discrete
// choice; current interpolated charges decide the instantaneous coefficient.
// This keeps force, energy, virial, and dH/dlambda on one Hamiltonian and also
// handles a pair that is inactive at both endpoints but active in between.
__host__ __device__ __forceinline__ LJ_SOFT_CORE_PAIR_RESULT
Evaluate_LJ_Soft_Core_Pair(const VECTOR_LJ_SOFT_TYPE r1,
                           const VECTOR_LJ_SOFT_TYPE r2, const float dr_abs,
                           const float AA, const float AB, const float BA,
                           const float BB, const float pme_beta,
                           const float lambda, const float alpha, const float p,
                           const float input_sigma_6,
                           const float input_sigma_6_min, const bool need_force,
                           const bool need_energy, const bool need_coulomb,
                           const bool need_du_dlambda)
{
    LJ_SOFT_CORE_PAIR_RESULT result;
    const float lambda_ = 1.0f - lambda;
    const bool lj_A_active =
        PairwiseInteraction::Lennard_Jones_Is_Active(AA, AB);
    const bool lj_B_active =
        PairwiseInteraction::Lennard_Jones_Is_Active(BA, BB);
    const bool coulomb_active =
        need_coulomb &&
        PairwiseInteraction::Coulomb_Is_Active(r1.charge, r2.charge);
    const bool coulomb_derivative_active =
        need_coulomb && need_du_dlambda &&
        PairwiseInteraction::Coulomb_Derivative_Is_Active(
            r1.charge, r2.charge, r1.charge_BA, r2.charge_BA);
    const PairwiseInteraction::Coulomb_Endpoint_Activity coulomb_endpoints =
        need_coulomb ? PairwiseInteraction::Classify_Coulomb_Endpoints(
                           r1.charge_A, r2.charge_A, r1.charge_B, r2.charge_B)
                     : PairwiseInteraction::Coulomb_Endpoint_Activity{};

    const bool lj_soft_core_needed = lj_A_active != lj_B_active;
    const bool coulomb_soft_core_needed =
        need_coulomb && coulomb_endpoints.Changes();
    result.any_interaction = lj_A_active || lj_B_active || coulomb_active ||
                             coulomb_derivative_active;
    if (!result.any_interaction)
    {
        return result;
    }

    float sigma_A = input_sigma_6;
    float sigma_B = input_sigma_6;
    float dr_softcore_A = dr_abs;
    float dr_softcore_B = dr_abs;
    if (lj_soft_core_needed || coulomb_soft_core_needed)
    {
        sigma_A = Get_Soft_Core_Sigma(AA, AB, input_sigma_6, input_sigma_6_min);
        sigma_B = Get_Soft_Core_Sigma(BA, BB, input_sigma_6, input_sigma_6_min);
        dr_softcore_A =
            Get_Soft_Core_Distance(AA, AB, sigma_A, dr_abs, alpha, p, lambda);
        dr_softcore_B =
            Get_Soft_Core_Distance(BB, BA, sigma_B, dr_abs, alpha, p, lambda_);
    }
    const float dr_lj_A = lj_soft_core_needed ? dr_softcore_A : dr_abs;
    const float dr_lj_B = lj_soft_core_needed ? dr_softcore_B : dr_abs;
    const float dr_coulomb_A =
        coulomb_soft_core_needed ? dr_softcore_A : dr_abs;
    const float dr_coulomb_B =
        coulomb_soft_core_needed ? dr_softcore_B : dr_abs;

    if (dr_abs == 0.0f)
    {
        const bool produces_value = need_force || need_energy;
        const bool lj_uses_A =
            lj_A_active &&
            (need_du_dlambda || (produces_value && lambda_ != 0.0f));
        const bool lj_uses_B =
            lj_B_active &&
            (need_du_dlambda || (produces_value && lambda != 0.0f));
        const bool coulomb_uses_A =
            (coulomb_active && (need_du_dlambda ||
                                (produces_value && (!coulomb_soft_core_needed ||
                                                    lambda_ != 0.0f)))) ||
            (coulomb_derivative_active &&
             (!coulomb_soft_core_needed || lambda_ != 0.0f));
        const bool coulomb_uses_B =
            (coulomb_active && (need_du_dlambda ||
                                (produces_value && (!coulomb_soft_core_needed ||
                                                    lambda != 0.0f)))) ||
            (coulomb_derivative_active &&
             (!coulomb_soft_core_needed || lambda != 0.0f));
        const bool singular_lj =
            (dr_lj_A == 0.0f && lj_uses_A) || (dr_lj_B == 0.0f && lj_uses_B);
        const bool singular_coulomb =
            (dr_coulomb_A == 0.0f && coulomb_uses_A) ||
            (dr_coulomb_B == 0.0f && coulomb_uses_B);
        result.singular_components =
            PairwiseInteraction::Components(singular_lj, singular_coulomb);
        if (result.singular_components !=
            PairwiseInteraction::PAIR_COMPONENT_NONE)
        {
            return result;
        }
    }

    if (need_force)
    {
        if (lj_A_active && lambda_ != 0.0f)
        {
            result.force += lambda_ * Get_Soft_Core_LJ_Force(r1, r2, dr_abs,
                                                             dr_lj_A, AA, AB);
        }
        if (lj_B_active && lambda != 0.0f)
        {
            result.force += lambda * Get_Soft_Core_LJ_Force(r1, r2, dr_abs,
                                                            dr_lj_B, BA, BB);
        }
        if (coulomb_active)
        {
            float coulomb_force = 0.0f;
            if (coulomb_soft_core_needed)
            {
                if (lambda_ != 0.0f)
                {
                    coulomb_force +=
                        lambda_ * Get_Soft_Core_Direct_Coulomb_Force(
                                      r1, r2, dr_abs, dr_coulomb_A, pme_beta);
                }
                if (lambda != 0.0f)
                {
                    coulomb_force +=
                        lambda * Get_Soft_Core_Direct_Coulomb_Force(
                                     r1, r2, dr_abs, dr_coulomb_B, pme_beta);
                }
            }
            else
            {
                coulomb_force =
                    Get_Direct_Coulomb_Force(r1, r2, dr_abs, pme_beta);
            }
            result.force -= coulomb_force;
        }
    }

    if (need_energy)
    {
        if (lj_A_active && lambda_ != 0.0f)
        {
            result.lj_energy +=
                lambda_ * Get_LJ_Energy(r1, r2, dr_lj_A, AA, AB);
        }
        if (lj_B_active && lambda != 0.0f)
        {
            result.lj_energy += lambda * Get_LJ_Energy(r1, r2, dr_lj_B, BA, BB);
        }
        if (coulomb_active)
        {
            if (coulomb_soft_core_needed)
            {
                if (lambda_ != 0.0f)
                {
                    result.coulomb_energy +=
                        lambda_ * Get_Direct_Coulomb_Energy(
                                      r1, r2, dr_coulomb_A, pme_beta);
                }
                if (lambda != 0.0f)
                {
                    result.coulomb_energy +=
                        lambda * Get_Direct_Coulomb_Energy(r1, r2, dr_coulomb_B,
                                                           pme_beta);
                }
            }
            else
            {
                result.coulomb_energy =
                    Get_Direct_Coulomb_Energy(r1, r2, dr_abs, pme_beta);
            }
        }
    }

    if (need_du_dlambda)
    {
        if (lj_B_active)
        {
            result.du_dlambda_lj += Get_LJ_Energy(r1, r2, dr_lj_B, BA, BB);
            if (lj_soft_core_needed)
            {
                result.du_dlambda_lj -= Get_Soft_Core_dU_dlambda(
                    Get_LJ_Force(r1, r2, dr_lj_B, BA, BB), sigma_B, dr_lj_B,
                    alpha, p, lambda_);
            }
        }
        if (lj_A_active)
        {
            result.du_dlambda_lj -= Get_LJ_Energy(r1, r2, dr_lj_A, AA, AB);
            if (lj_soft_core_needed)
            {
                result.du_dlambda_lj += Get_Soft_Core_dU_dlambda(
                    Get_LJ_Force(r1, r2, dr_lj_A, AA, AB), sigma_A, dr_lj_A,
                    alpha, p, lambda);
            }
        }

        if (coulomb_active && coulomb_soft_core_needed)
        {
            result.du_dlambda_coulomb +=
                Get_Direct_Coulomb_Energy(r1, r2, dr_coulomb_B, pme_beta) -
                Get_Direct_Coulomb_Energy(r1, r2, dr_coulomb_A, pme_beta);
            result.du_dlambda_coulomb +=
                Get_Soft_Core_dU_dlambda(
                    Get_Direct_Coulomb_Force(r1, r2, dr_coulomb_B, pme_beta),
                    sigma_B, dr_coulomb_B, alpha, p, lambda_) -
                Get_Soft_Core_dU_dlambda(
                    Get_Direct_Coulomb_Force(r1, r2, dr_coulomb_A, pme_beta),
                    sigma_A, dr_coulomb_A, alpha, p, lambda);
        }
        if (coulomb_derivative_active)
        {
            if (coulomb_soft_core_needed)
            {
                if (lambda != 0.0f)
                {
                    result.du_dlambda_coulomb +=
                        lambda * Get_Direct_Coulomb_dU_dlambda(
                                     r1, r2, dr_coulomb_B, pme_beta);
                }
                if (lambda_ != 0.0f)
                {
                    result.du_dlambda_coulomb +=
                        lambda_ * Get_Direct_Coulomb_dU_dlambda(
                                      r1, r2, dr_coulomb_A, pme_beta);
                }
            }
            else
            {
                result.du_dlambda_coulomb +=
                    Get_Direct_Coulomb_dU_dlambda(r1, r2, dr_abs, pme_beta);
            }
        }
    }
    return result;
}

__global__ void Copy_LJ_Type_And_Mask_To_New_Crd(const int atom_numbers,
                                                 VECTOR_LJ_SOFT_TYPE* new_crd,
                                                 const int* LJ_type_A,
                                                 const int* LJ_type_B,
                                                 const int* mask);
__global__ void Copy_Crd_And_Charge_To_New_Crd(const int atom_numbers,
                                               const VECTOR* crd,
                                               VECTOR_LJ_SOFT_TYPE* new_crd,
                                               const float* charge);
__global__ void Copy_Crd_And_Charge_To_New_Crd(
    const int atom_numbers, const VECTOR* crd, VECTOR_LJ_SOFT_TYPE* new_crd,
    const float* charge, const float lambda, int* endpoint_error);
__global__ void Copy_Crd_To_New_Crd(const int atom_numbers, const VECTOR* crd,
                                    VECTOR_LJ_SOFT_TYPE* new_crd);
#endif

struct LJ_SOFT_CORE
{
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260216;
    CONTROLLER* controller = NULL;

    int atom_numbers = 0;
    int atom_type_numbers_A = 0;
    int atom_type_numbers_B = 0;
    int pair_type_numbers_A = 0;
    int pair_type_numbers_B = 0;

    int* h_atom_LJ_type_A = NULL;
    int* h_atom_LJ_type_B = NULL;
    int* d_atom_LJ_type_A = NULL;
    int* d_atom_LJ_type_B = NULL;

    float* h_LJ_AA = NULL;
    float* h_LJ_AB = NULL;
    float* h_LJ_BA = NULL;
    float* h_LJ_BB = NULL;
    float* d_LJ_AA = NULL;
    float* d_LJ_AB = NULL;
    float* d_LJ_BA = NULL;
    float* d_LJ_BB = NULL;

    float* h_LJ_energy_atom = NULL;
    float h_LJ_energy_sum = 0;
    float* d_LJ_energy_atom = NULL;
    float* d_LJ_energy_sum = NULL;

    int* d_subsys_division = NULL;
    int* h_subsys_division = NULL;

    float* h_LJ_energy_atom_intersys = NULL;
    float* h_LJ_energy_atom_intrasys = NULL;
    float h_LJ_energy_sum_intersys = 0;
    float h_LJ_energy_sum_intrasys = 0;
    float* d_LJ_energy_atom_intersys = NULL;
    float* d_LJ_energy_atom_intrasys = NULL;
    float* d_LJ_energy_sum_intersys = NULL;
    float* d_LJ_energy_sum_intrasys = NULL;

    float* d_direct_ene_sum_intersys = NULL;
    float* d_direct_ene_sum_intrasys = NULL;
    float h_direct_ene_sum = 0.0;
    float h_direct_ene_sum_intersys = 0.0;
    float h_direct_ene_sum_intrasys = 0.0;

    float lambda;
    float alpha;
    float p;
    float sigma_6;
    float sigma;
    float sigma_min;
    float sigma_6_min;

    int has_charge_endpoints = 0;
    float* h_charge_A = NULL;
    float* h_charge_B = NULL;
    float* d_charge_A = NULL;
    float* d_charge_B = NULL;

    float* h_sigma_of_dH_dlambda_lj = NULL;
    float* d_sigma_of_dH_dlambda_lj = NULL;

    float* h_sigma_of_dH_dlambda_direct = NULL;
    float* d_sigma_of_dH_dlambda_direct = NULL;

    float cutoff = 10.0;
    VECTOR_LJ_SOFT_TYPE* crd_with_parameters = NULL;
    float h_LJ_long_energy = 0.0;
    float long_range_factor = 0.0;
    float long_range_factor_TI = 0.0;

    void Initial(CONTROLLER* controller, float cutoff,
                 char* module_name = NULL);

    void LJ_Soft_Core_Malloc();

    void Clear();

    void Parameter_Host_To_Device();

    void LJ_Soft_Core_PME_Direct_Force_With_Atom_Energy_And_Virial(
        const int atom_numbers, const int local_atom_numbers,
        const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
        const float* charge, VECTOR* frc, const LTMatrix3 cell,
        const LTMatrix3 rcell, const ATOM_GROUP* nl, const float pme_beta,
        const int need_atom_energy, float* atom_energy, const int need_virial,
        LTMatrix3* atom_lj_virial, float* atom_direct_pme_energy);

    void Step_Print(CONTROLLER* controller);

    void Long_Range_Correction(int need_pressure, LTMatrix3* d_virial,
                               int need_potential, float* d_potential,
                               const float volume);

    float Get_Partial_H_Partial_Lambda_With_Columb_Direct(
        const int solvent_numbers, const VECTOR* crd, const LTMatrix3 cell,
        const LTMatrix3 rcell, const float* charge, const ATOM_GROUP* nl,
        const float pme_beta);

    /*
        以下用于区域分解
    */
    int local_atom_numbers = 0;
    int ghost_numbers = 0;
    VECTOR_LJ_SOFT_TYPE* crd_with_LJ_parameters_local =
        NULL;  // 局域原子的坐标，电荷LJ_type打包
    int* d_pair_overlap_error = NULL;
    int* d_local_metadata_error = NULL;
    int* d_charge_endpoint_error = NULL;
    bool local_metadata_is_ready = false;
    void Get_Local(int* atom_local, int local_atom_numbers,
                   int ghost_numbers);  // 获取局域粒子信息
    bool Validate_Local_State(const char* error_by, int global_atom_numbers,
                              int local_atom_numbers, int ghost_numbers,
                              int solvent_numbers);
    bool Prepare_Local_Coordinates(const char* error_by,
                                   const int coordinate_numbers,
                                   const VECTOR* crd, const float* charge);
    void Reset_Pair_Overlap_Error();
    bool Check_Pair_Overlap_Error(const char* error_by);
};
