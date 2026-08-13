#include "SITS.h"

#include <cctype>
#include <cerrno>
#include <cstdint>

#include "../Lennard_Jones_force/pair_activity.h"

template <bool need_force, bool need_energy, bool need_virial,
          bool need_coulomb>
static __global__ void Selective_Lennard_Jones_And_Direct_Coulomb_Device(
    const int local_atom_numbers, const int solvent_numbers,
    const ATOM_GROUP* nl, float* atom_ene_LJ, const float4* crd_q,
    const int2* type_g, const LTMatrix3 cell, const LTMatrix3 rcell,
    const float* LJ_type_A, const float* LJ_type_B, const int* atom_sys_mark,
    const float cutoff, VECTOR* frc, VECTOR* frc_enhancing,
    const float pme_beta, float* atom_energy, float* atom_energy_enhancing,
    LTMatrix3* atom_virial, LTMatrix3* atom_virial_enhancing,
    float* atom_direct_cf_energy, const float pwwp_factor,
    int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = blockDim.y * blockIdx.x + threadIdx.y;
    if (atom_i < local_atom_numbers - solvent_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < local_atom_numbers - solvent_numbers;
         atom_i++)
#endif
    {
        ATOM_GROUP nl_i = nl[atom_i];
        VECTOR_LJ r1 = Load_VECTOR_LJ(crd_q, type_g, atom_i);
        int atom_mark_i = atom_sys_mark[atom_i];
        VECTOR frc_record = {0.0f, 0.0f, 0.0f},
               frc_enhancing_record = {0.0f, 0.0f, 0.0f};
        LTMatrix3 virial_record = {0, 0, 0, 0, 0, 0},
                  virial_enhancing = {0, 0, 0, 0, 0, 0};
        float energy_lj = 0.0f, energy_enhancing = 0.0f, energy_coulomb = 0.0f;
#ifdef USE_GPU
        for (int j = threadIdx.x; j < nl_i.atom_numbers; j += blockDim.x)
#else
        for (int j = 0; j < nl_i.atom_numbers; j++)
#endif
        {
            int atom_j = nl_i.atom_serial[j];
            float ij_factor = atom_j < local_atom_numbers ? 1.0f : 0.5f;
            VECTOR_LJ r2 = Load_VECTOR_LJ(crd_q, type_g, atom_j);
            VECTOR dr = Get_Periodic_Displacement(r2, r1, cell, rcell);
            float dr_abs = norm3df(dr.x, dr.y, dr.z);
            if (dr_abs < cutoff)
            {
                int atom_mark_j = atom_sys_mark[atom_j] + atom_mark_i;
                int atom_pair_LJ_type = Get_LJ_Type(r1.LJ_type, r2.LJ_type);
                float A = LJ_type_A[atom_pair_LJ_type];
                float B = LJ_type_B[atom_pair_LJ_type];
                const PairwiseInteraction::Pair_Activity activity =
                    PairwiseInteraction::Classify(
                        A, B,
                        need_coulomb && PairwiseInteraction::Coulomb_Is_Active(
                                            r1.charge, r2.charge));
                if (!activity.Any())
                {
                    continue;
                }
                if (dr_abs == 0.0f)
                {
                    PairwiseInteraction::Fail_Exact_Overlap(
                        r1.global_atom, r2.global_atom,
                        PairwiseInteraction::Components(activity.lennard_jones,
                                                        activity.coulomb),
                        pair_overlap_error);
                    continue;
                }
                float factor = 0;
                if (atom_mark_j == 0)
                {
                    factor = 1;
                }
                else if (atom_mark_j == 1)
                {
                    factor = pwwp_factor;
                }
                if (need_force)
                {
                    float frc_abs = 0.0f;
                    if (activity.lennard_jones)
                    {
                        frc_abs = Get_LJ_Force(r1, r2, dr_abs, A, B);
                    }
                    if (activity.coulomb)
                    {
                        float frc_cf_abs =
                            Get_Direct_Coulomb_Force(r1, r2, dr_abs, pme_beta);
                        frc_abs = frc_abs - frc_cf_abs;
                    }
                    const VECTOR frc_lin = frc_abs * dr;
                    const VECTOR frc_enhancing_lin = factor * frc_lin;
                    frc_record = frc_record + frc_lin;
                    if (atom_j < local_atom_numbers)
                    {
                        atomicAdd(frc + atom_j, -frc_lin);
                        atomicAdd(frc_enhancing + atom_j, -frc_enhancing_lin);
                    }
                    frc_enhancing_record =
                        frc_enhancing_record + frc_enhancing_lin;
                    if (need_virial)
                    {
                        const LTMatrix3 pair_virial =
                            Get_Virial_From_Force_Dis(frc_lin, dr);
                        virial_record = virial_record - ij_factor * pair_virial;
                        virial_enhancing =
                            virial_enhancing - ij_factor * factor * pair_virial;
                    }
                }
                if (activity.coulomb && need_energy)
                {
                    float energy_lin =
                        Get_Direct_Coulomb_Energy(r1, r2, dr_abs, pme_beta);
                    energy_coulomb += ij_factor * energy_lin;
                    energy_enhancing += ij_factor * factor * energy_lin;
                }
                if (need_energy)
                {
                    if (activity.lennard_jones)
                    {
                        float energy_lin = Get_LJ_Energy(r1, r2, dr_abs, A, B);
                        energy_lj += ij_factor * energy_lin;
                        energy_enhancing += ij_factor * factor * energy_lin;
                    }
                }
            }
        }
        if (need_force)
        {
            Warp_Sum_To(frc + atom_i, frc_record, warpSize);
            Warp_Sum_To(frc_enhancing + atom_i, frc_enhancing_record, warpSize);
        }
        if (need_coulomb && need_energy)
        {
            Warp_Sum_To(atom_direct_cf_energy + atom_i, energy_coulomb,
                        warpSize);
        }
        if (need_energy)
        {
            float energy_total = energy_lj;
            if (need_coulomb)
            {
                energy_total += energy_coulomb;
            }
            Warp_Sum_To(atom_energy + atom_i, energy_total, warpSize);
#ifdef USE_GPU
            if (threadIdx.x == 0)
#endif
                atomicAdd(atom_ene_LJ + atom_i, energy_lj);
            Warp_Sum_To(atom_energy_enhancing + atom_i, energy_enhancing,
                        warpSize);
        }
        if (need_virial)
        {
            Warp_Sum_To(atom_virial + atom_i, virial_record, warpSize);
            Warp_Sum_To(atom_virial_enhancing + atom_i, virial_enhancing,
                        warpSize);
        }
    }
}

template <bool need_force, bool need_energy, bool need_virial,
          bool need_coulomb, bool need_du_dlambda>
static __global__ void
Selective_Lennard_Jones_And_Direct_Coulomb_Soft_Core_Device(
    const int local_atom_numbers, const int solvent_numbers,
    const ATOM_GROUP* nl, float* atom_ene_LJ, const VECTOR_LJ_SOFT_TYPE* crd,
    const LTMatrix3 cell, const LTMatrix3 rcell, const int* atom_sys_mark,
    const float* LJ_type_AA, const float* LJ_type_AB, const float* LJ_type_BA,
    const float* LJ_type_BB, const float cutoff, VECTOR* frc,
    VECTOR* frc_enhancing, const float pme_beta, float* atom_energy,
    float* atom_energy_enhancing, LTMatrix3* atom_virial,
    LTMatrix3* atom_virial_enhancing, float* atom_direct_cf_energy,
    float* atom_du_dlambda_lj, float* atom_du_dlambda_direct,
    float* atom_du_dlambda_enhancing, const float lambda, const float alpha,
    const float p, const float input_sigma_6, const float input_sigma_6_min,
    const float pwwp_factor, int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = blockDim.y * blockIdx.x + threadIdx.y;
    if (atom_i < local_atom_numbers - solvent_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < local_atom_numbers - solvent_numbers;
         atom_i++)
#endif
    {
        ATOM_GROUP nl_i = nl[atom_i];
        VECTOR_LJ_SOFT_TYPE r1 = crd[atom_i];
        VECTOR frc_record = {0., 0., 0.},
               frc_enhancing_record = {0.0f, 0.0f, 0.0f};
        LTMatrix3 virial_record = {0, 0, 0, 0, 0, 0},
                  virial_enhancing = {0, 0, 0, 0, 0, 0};
        float energy_total = 0., energy_enhancing = 0.0f;
        float energy_coulomb = 0.;
        float du_dlambda_lj = 0.;
        float du_dlambda_direct = 0.;
        // float du_dlambda_enhancing = 0.0f;
        int atom_mark_i = atom_sys_mark[atom_i];
#ifdef USE_GPU
        for (int j = threadIdx.x; j < nl_i.atom_numbers; j += blockDim.x)
#else
        for (int j = 0; j < nl_i.atom_numbers; j++)
#endif
        {
            int atom_j = nl_i.atom_serial[j];
            float ij_factor = atom_j < local_atom_numbers ? 1.0f : 0.5f;
            VECTOR_LJ_SOFT_TYPE r2 = crd[atom_j];
            VECTOR dr = Get_Periodic_Displacement(r2, r1, cell, rcell);
            float dr_abs = norm3df(dr.x, dr.y, dr.z);
            if (dr_abs < cutoff)
            {
                int atom_mark_j = atom_sys_mark[atom_j] + atom_mark_i;
                float factor = 0;
                if (atom_mark_j == 0)
                {
                    factor = 1;
                }
                else if (atom_mark_j == 1)
                {
                    factor = pwwp_factor;
                }
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
                    const VECTOR frc_enhancing_lin = factor * frc_lin;
                    frc_record = frc_record + frc_lin;
                    frc_enhancing_record =
                        frc_enhancing_record + frc_enhancing_lin;
                    if (atom_j < local_atom_numbers)
                    {
                        atomicAdd(frc + atom_j, -frc_lin);
                        atomicAdd(frc_enhancing + atom_j, -frc_enhancing_lin);
                    }
                    if (need_virial)
                    {
                        const LTMatrix3 pair_virial =
                            Get_Virial_From_Force_Dis(frc_lin, dr);
                        virial_record = virial_record - ij_factor * pair_virial;
                        virial_enhancing =
                            virial_enhancing - ij_factor * factor * pair_virial;
                    }
                }
                if (need_energy)
                {
                    const float pair_energy =
                        pair.lj_energy + pair.coulomb_energy;
                    energy_total += ij_factor * pair.lj_energy;
                    energy_coulomb += ij_factor * pair.coulomb_energy;
                    energy_enhancing += ij_factor * factor * pair_energy;
                }
                if (need_du_dlambda)
                {
                    du_dlambda_lj += ij_factor * pair.du_dlambda_lj;
                    du_dlambda_direct += ij_factor * pair.du_dlambda_coulomb;
                }
            }
        }
        if (need_force)
        {
            Warp_Sum_To(frc + atom_i, frc_record, warpSize);
            Warp_Sum_To(frc_enhancing + atom_i, frc_enhancing_record, warpSize);
        }
        if (need_coulomb && need_energy)
        {
            Warp_Sum_To(atom_direct_cf_energy + atom_i, energy_coulomb,
                        warpSize);
        }
        if (need_energy)
        {
            float full_energy = energy_total;
            if (need_coulomb)
            {
                full_energy += energy_coulomb;
            }
            Warp_Sum_To(atom_energy + atom_i, full_energy, warpSize);
#ifdef USE_GPU
            if (threadIdx.x == 0)
#endif
                atomicAdd(atom_ene_LJ + atom_i, energy_total);
            Warp_Sum_To(atom_energy_enhancing + atom_i, energy_enhancing,
                        warpSize);
        }
        if (need_virial)
        {
            Warp_Sum_To(atom_virial + atom_i, virial_record, warpSize);
            Warp_Sum_To(atom_virial_enhancing + atom_i, virial_enhancing,
                        warpSize);
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

static __device__ float log_add_log(float a, float b)
{
    return fmaxf(a, b) + logf(1.0 + expf(-fabsf(a - b)));
}

static __global__ void SITS_Record_Ene_Device(float* ene_record,
                                              const float* enhancing_energy,
                                              const float pe_a,
                                              const float pe_b)
{
    *ene_record = pe_a * *enhancing_energy + pe_b;
}

static __global__ void SITS_Update_gf_Device(const int kn, float* gf,
                                             const float* ene_record,
                                             const float* log_nk,
                                             const float* beta_k)
{
#ifdef USE_GPU
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < kn)
#else
#pragma omp parallel for
    for (int i = 0; i < kn; i++)
#endif
    {
        gf[i] = -beta_k[i] * ene_record[0] + log_nk[i];
    }
}

static __global__ void SITS_Update_gfsum_Device(const int kn, float* gfsum,
                                                const float* gf)
{
    float temp = -FLT_MAX;
    for (int i = 0; i < kn; i = i + 1)
    {
        temp = log_add_log(temp, gf[i]);
    }
    gfsum[0] = temp;
}

static __global__ void SITS_Update_log_pk_Device(const int kn, float* log_pk,
                                                 const float* gf,
                                                 const float* gfsum,
                                                 const int reset)
{
#ifdef USE_GPU
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < kn)
#else
#pragma omp parallel for
    for (int i = 0; i < kn; i++)
#endif
    {
        float gfi = gf[i];
        log_pk[i] =
            ((float)reset) * gfi +
            ((float)(1 - reset)) * log_add_log(log_pk[i], gfi - gfsum[0]);
    }
}

static __global__ void SITS_Update_log_mk_inverse_Device(
    const int kn, float* log_weight, float* log_mk_inverse, float* log_norm_old,
    float* log_norm, const float* log_pk, const float* log_nk)
{
#ifdef USE_GPU
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < kn - 1)
#else
#pragma omp parallel for
    for (int i = 0; i < kn - 1; i++)
#endif
    {
        log_weight[i] = (log_pk[i] + log_pk[i + 1]) * 0.5;
        log_mk_inverse[i] = log_nk[i] - log_nk[i + 1];
        log_norm_old[i] = log_norm[i];
        log_norm[i] = log_add_log(log_norm[i], log_weight[i]);
        log_mk_inverse[i] =
            log_add_log(log_mk_inverse[i] + log_norm_old[i] - log_norm[i],
                        log_pk[i + 1] - log_pk[i] + log_mk_inverse[i] +
                            log_weight[i] - log_norm[i]);
    }
}

static __global__ void SITS_Update_log_nk_inverse_Device(
    const int kn, float* log_nk_inverse, const float* log_mk_inverse)
{
    for (int i = 0; i < kn - 1; i++)
    {
        log_nk_inverse[i + 1] = log_nk_inverse[i] + log_mk_inverse[i];
    }
}

static __global__ void SITS_Update_nk_Device(const int kn, float* log_nk,
                                             float* nk,
                                             const float* log_nk_inverse)
{
#ifdef USE_GPU
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < kn)
#else
#pragma omp parallel for
    for (int i = 0; i < kn; i++)
#endif
    {
        log_nk[i] = -log_nk_inverse[i];
        nk[i] = exp(log_nk[i]);
    }
}

static __global__ void SITS_For_Enhanced_Force_Calculate_NkExpBetakU_Device(
    const int k_numbers, const float* beta_k, const float* log_nk,
    float* nkexpbetaku, const float* ene, const float beta0, const float pe_a,
    const float pe_b)
{
#ifdef USE_GPU
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < k_numbers)
#else
#pragma omp parallel for
    for (int i = 0; i < k_numbers; i++)
#endif
    {
        nkexpbetaku[i] =
            -(beta_k[i] - beta0) * (pe_a * ene[0] + pe_b) + log_nk[i];
    }
}

static __global__ void SITS_For_Enhanced_Force_Sum_Of_Above_And_Below_Device(
    const int k_numbers, const float* nkexpbetaku, const float* beta_k,
    float* d_bias, float pe_a, float pe_b, float* sum_of_above,
    float* sum_of_below, float* factor, float beta0, float fb_bias,
    const float* h_enhancing_energy)
{
    float above = -FLT_MAX;
    float below = -FLT_MAX;
    for (int i = 0; i < k_numbers; i++)
    {
        above = log_add_log(above, logf(beta_k[i]) + nkexpbetaku[i]);
        below = log_add_log(below, nkexpbetaku[i]);
    }
    sum_of_above[0] = above;
    sum_of_below[0] = below;
    factor[0] = expf(above - below - logf(beta0)) + fb_bias;
    d_bias[0] =
        -below / beta0 / pe_a + fb_bias * (h_enhancing_energy[0] + pe_b / pe_a);
}

static __global__ void SITS_For_Enhanced_Force_Protein_Water_Device(
    const int atom_numbers, VECTOR* md_frc, const VECTOR* enhancing_frc,
    float* md_ene, const float* bias, const int need_pressure,
    LTMatrix3* md_virial, const LTMatrix3* virial_enhancing,
    const float factor_minus_one)
{
#ifdef USE_GPU
    if (blockIdx.x == 0 && threadIdx.x == 0)
#endif
    {
        md_ene[0] = md_ene[0] + bias[0];
        if (need_pressure)
        {
            md_virial[0] =
                md_virial[0] + factor_minus_one * virial_enhancing[0];
        }
    }
#ifdef USE_GPU
    __syncthreads();
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < atom_numbers)
#else
#pragma omp parallel for
    for (int i = 0; i < atom_numbers; i++)
#endif
    {
        md_frc[i] = md_frc[i] + factor_minus_one * enhancing_frc[i];
    }
}

static __global__ void ESITS_Get_Current_Fb(const float* enhancing_energy,
                                            float* factor, const float pe_a,
                                            const float pe_b,
                                            const float low_temperature_ratio,
                                            const float high_temperature_ratio,
                                            float* d_bias)
{
    const float energy = enhancing_energy[0];
    if (energy > pe_b)
    {
        // Integrate dB/dU = factor - 1 with B(E) = 0.  Computing the
        // logarithm and the cancellation-prone expression in double
        // precision keeps the two branches continuous near U = E.
        const double energy_above_threshold =
            static_cast<double>(energy) - static_cast<double>(pe_b);
        const double smoothing_energy = static_cast<double>(pe_a);
        const double low_temperature_factor =
            static_cast<double>(low_temperature_ratio);
        const double high_temperature_factor =
            static_cast<double>(high_temperature_ratio);
        const double temperature_factor_difference =
            low_temperature_factor - high_temperature_factor;
        factor[0] = static_cast<float>(
            low_temperature_factor -
            temperature_factor_difference * smoothing_energy /
                (energy_above_threshold + smoothing_energy));
        d_bias[0] = static_cast<float>(
            (low_temperature_factor - 1.0) * energy_above_threshold -
            temperature_factor_difference * smoothing_energy *
                log1p(energy_above_threshold / smoothing_energy));
    }
    else
    {
        factor[0] = high_temperature_ratio;
        d_bias[0] = (high_temperature_ratio - 1.0f) * (energy - pe_b);
    }
}

static __global__ void AMD_Get_Current_Fb(const float* enhancing_energy,
                                          float* factor, const float pe_a,
                                          const float pe_b, float* d_bias)
{
    float ene = enhancing_energy[0];
    if (ene < pe_b)
    {
        ene = pe_b - ene;
        factor[0] = 1.0f - ene * (2 * pe_a + ene) / (pe_a + ene) / (pe_a + ene);
        d_bias[0] = ene * ene / (pe_a + ene);
    }
    else
    {
        factor[0] = 1.0f;
        d_bias[0] = 0;
    }
}

static __global__ void GAMD_Get_Current_Fb(const float* enhancing_energy,
                                           float* factor, const float pe_a,
                                           const float pe_b, float* d_bias)
{
    const float energy = enhancing_energy[0];
    if (energy < pe_b)
    {
        const float energy_below_threshold = pe_b - energy;
        factor[0] = 1.0f - pe_a * energy_below_threshold;
        d_bias[0] =
            0.5f * pe_a * energy_below_threshold * energy_below_threshold;
    }
    else
    {
        factor[0] = 1.0f;
        d_bias[0] = 0;
    }
}

static __global__ void SITS_Record_Fb_Reference(const float* enhancing_energy,
                                                const float* bias,
                                                float* reference_energy,
                                                float* reference_bias)
{
    reference_energy[0] = enhancing_energy[0];
    reference_bias[0] = bias[0];
}

static __global__ void SITS_Evaluate_Linearized_Fb(
    const float* enhancing_energy, const float* factor,
    const float* reference_energy, const float* reference_bias, float* bias)
{
    bias[0] = reference_bias[0] +
              (factor[0] - 1.0f) * (enhancing_energy[0] - reference_energy[0]);
}

static int SITS_Kernel_Block_Count(const int item_count)
{
    // item_count + 63 overflows for otherwise valid large int counts.  This
    // form computes ceil(item_count / 64) over the complete positive range.
    return item_count > 0 ? 1 + (item_count - 1) / 64 : 0;
}

static void SITS_Get_Current_Fb(const int atom_numbers,
                                const float* energy_enhancing, float* d_bias,
                                const int k_numbers, float* nkexpbetaku,
                                const float* beta_k, const float* log_nk,
                                const float beta0, float* sum_a, float* sum_b,
                                float* factor, const float fb_bias,
                                const float pe_a, const float pe_b,
                                const float pwwp_enhance_factor)
{
    Launch_Device_Kernel(SITS_For_Enhanced_Force_Calculate_NkExpBetakU_Device,
                         SITS_Kernel_Block_Count(k_numbers), 64, 0, NULL,
                         k_numbers, beta_k, log_nk, nkexpbetaku,
                         energy_enhancing, beta0, pe_a, pe_b);

    Launch_Device_Kernel(SITS_For_Enhanced_Force_Sum_Of_Above_And_Below_Device,
                         1, 1, 0, NULL, k_numbers, nkexpbetaku, beta_k, d_bias,
                         pe_a, pe_b, sum_a, sum_b, factor, beta0, fb_bias,
                         energy_enhancing);
}

static void Validate_SITS_Float(CONTROLLER* controller, const char* module_name,
                                const char* command_name, const float value,
                                const bool must_be_positive)
{
    if (!Float_Memory_Is_Finite(&value) ||
        (must_be_positive && !(value > 0.0f)))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "CLASSIC_SITS_INFORMATION::Initial",
            "Reason:\n\t%s_%s must be %sfinite (received %.9g)\n", module_name,
            command_name, must_be_positive ? "positive and " : "", value);
    }
}

static void Parse_SITS_Temperature_List(CONTROLLER* controller,
                                        const char* module_name,
                                        const char* input, const int count,
                                        float* beta_k)
{
    const char* cursor = input;
    float previous_temperature = 0.0f;
    for (int i = 0; i < count; i++)
    {
        while (std::isspace(static_cast<unsigned char>(*cursor))) cursor++;
        errno = 0;
        char* end = NULL;
        const float temperature = strtof(cursor, &end);
        if (end == cursor || errno == ERANGE ||
            !Float_Memory_Is_Finite(&temperature) || !(temperature > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\t%s_T entry %d must be a positive finite "
                "temperature\n",
                module_name, i);
        }
        if (i > 0 && !(temperature > previous_temperature))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\t%s_T temperatures must be strictly increasing "
                "(entry %d is %.9g after %.9g)\n",
                module_name, i, temperature, previous_temperature);
        }
        previous_temperature = temperature;
        beta_k[i] = 1.0f / (CONSTANT_kB * temperature);
        cursor = end;
        while (std::isspace(static_cast<unsigned char>(*cursor))) cursor++;
        if (i + 1 < count)
        {
            if (*cursor != '/')
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "CLASSIC_SITS_INFORMATION::Initial",
                    "Reason:\n\t%s_T must contain exactly %d "
                    "slash-separated temperatures\n",
                    module_name, count);
            }
            cursor++;
        }
        else if (*cursor != '\0')
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\t%s_T contains more than %d temperatures or an "
                "unexpected trailing token\n",
                module_name, count);
        }
    }
}

void CLASSIC_SITS_INFORMATION::Initial(CONTROLLER* controller,
                                       SITS_INFORMATION* sits)
{
    is_initialized = 1;
    sits_controller = sits;
    record_count = 0;
    fb_interval = 1;
    Device_Malloc_Safely((void**)&d_bias, sizeof(float));
    if (controller->Command_Exist(sits->module_name, "fb_interval"))
    {
        controller->Check_Int(sits->module_name, "fb_interval",
                              "CLASSIC_SITS_INFORMATION::Initial");
        fb_interval =
            atoi(controller->Command(sits->module_name, "fb_interval"));
    }
    if (fb_interval <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "CLASSIC_SITS_INFORMATION::Initial",
            "Reason:\n\t%s_fb_interval must be positive (received %d)\n",
            sits->module_name, fb_interval);
    }
    controller->printf("    SITS fb update interval set to %d\n", fb_interval);
    if (sits->sits_mode == SITS_MODE_AMD)
    {
        if (controller->Command_Exist(sits->module_name, "pe_a"))
        {
            controller->Check_Float(sits->module_name, "pe_a",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_a = atof(controller->Command(sits->module_name, "pe_a"));
        }
        else
        {
            controller->Throw_SPONGE_Error(
                spongeErrorMissingCommand, "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\tAlpha (pe_a) is required for the Accelerated MD");
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_a", pe_a, true);
        controller->printf("    AMD alpha (pe_a) set to %f\n", pe_a);

        if (controller->Command_Exist(sits->module_name, "pe_b"))
        {
            controller->Check_Float(sits->module_name, "pe_b",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_b = atof(controller->Command(sits->module_name, "pe_b"));
        }
        else
        {
            controller->Throw_SPONGE_Error(
                spongeErrorMissingCommand, "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\tE (pe_b) is required for the Accelerated MD");
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_b", pe_b, false);
        controller->printf("    AMD E (pe_b) set to %f\n", pe_b);

        k_numbers = 0;
        nk_fix = 1;
        record_interval = 1;
        update_interval = INT_MAX;
        Memory_Allocate();
    }
    else if (sits->sits_mode == SITS_MODE_GAMD)
    {
        if (controller->Command_Exist(sits->module_name, "pe_a"))
        {
            controller->Check_Float(sits->module_name, "pe_a",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_a = atof(controller->Command(sits->module_name, "pe_a"));
        }
        else
        {
            controller->Throw_SPONGE_Error(spongeErrorMissingCommand,
                                           "CLASSIC_SITS_INFORMATION::Initial",
                                           "Reason:\n\tk (pe_a) is required "
                                           "for the Gaussian Accelerated MD");
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_a", pe_a, true);
        controller->printf("    GAMD k (pe_a) set to %f\n", pe_a);

        if (controller->Command_Exist(sits->module_name, "pe_b"))
        {
            controller->Check_Float(sits->module_name, "pe_b",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_b = atof(controller->Command(sits->module_name, "pe_b"));
        }
        else
        {
            controller->Throw_SPONGE_Error(
                spongeErrorMissingCommand, "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\tE (pe_b) is required for the Accelerated MD");
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_b", pe_b, false);
        controller->printf("    GAMD E (pe_b) set to %f\n", pe_b);

        k_numbers = 0;
        nk_fix = 1;
        record_interval = 1;
        update_interval = INT_MAX;
        Memory_Allocate();
    }
    else if (sits->sits_mode == SITS_MODE_EMPIRICAL)
    {
        if (controller->Command_Exist(sits->module_name, "pe_a"))
        {
            controller->Check_Float(sits->module_name, "pe_a",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_a = atof(controller->Command(sits->module_name, "pe_a"));
        }
        else
        {
            pe_a = 1.0;
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_a", pe_a, true);
        controller->printf("    SITS_pe_a set to %f\n", pe_a);

        if (controller->Command_Exist(sits->module_name, "pe_b"))
        {
            controller->Check_Float(sits->module_name, "pe_b",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_b = atof(controller->Command(sits->module_name, "pe_b"));
        }
        else
        {
            pe_b = 0.0;
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_b", pe_b, false);
        controller->printf("    SITS_pe_b set to %f\n", pe_b);

        if (!controller->Command_Exist(sits->module_name, "T_low") ||
            !controller->Command_Exist(sits->module_name, "T_high"))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorMissingCommand, "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\tSITS_T_high and SITS_T_low are required for "
                "empirical SITS");
        }
        controller->Check_Float(sits->module_name, "T_low",
                                "CLASSIC_SITS_INFORMATION::Initial");
        controller->Check_Float(sits->module_name, "T_high",
                                "CLASSIC_SITS_INFORMATION::Initial");
        T_low = atof(controller->Command(sits->module_name, "T_low"));
        T_high = atof(controller->Command(sits->module_name, "T_high"));
        Validate_SITS_Float(controller, sits->module_name, "T_low", T_low,
                            true);
        Validate_SITS_Float(controller, sits->module_name, "T_high", T_high,
                            true);
        if (!(T_high > T_low))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\t%s_T_high (%.9g) must be greater than %s_T_low "
                "(%.9g)\n",
                sits->module_name, T_high, sits->module_name, T_low);
        }
        controller->printf("    SITS_T_high set to %f\n", T_high);
        controller->printf("    SITS_T_low set to %f\n", T_low);

        k_numbers = 0;
        nk_fix = 1;
        record_interval = 1;
        update_interval = INT_MAX;
        Memory_Allocate();
    }
    else if (sits->sits_mode != SITS_MODE_OBSERVATION)
    {
        if (controller->Command_Exist(sits->module_name, "k_numbers"))
        {
            controller->Check_Int(sits->module_name, "k_numbers",
                                  "CLASSIC_SITS_INFORMATION::Initial");
            k_numbers =
                atoi(controller->Command(sits->module_name, "k_numbers"));
        }
        else
        {
            k_numbers = 40;
        }
        if (k_numbers < 2)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\t%s_k_numbers must be at least 2 (received %d)\n",
                sits->module_name, k_numbers);
        }
        controller->printf("    k numbers is %d\n", k_numbers);
        Memory_Allocate();

        controller->printf("    Read %s temperature information.\n",
                           sits->module_name);
        float* beta_k_tmp;
        Malloc_Safely((void**)&beta_k_tmp, sizeof(float) * k_numbers);
        const bool has_T_low =
            controller->Command_Exist(sits->module_name, "T_low");
        const bool has_T_high =
            controller->Command_Exist(sits->module_name, "T_high");
        const bool has_T_list =
            controller->Command_Exist(sits->module_name, "T");
        if (has_T_low != has_T_high)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorMissingCommand, "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\tSITS_T_low and SITS_T_high must be provided "
                "together\n");
        }
        if (has_T_list && has_T_low)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorConflictingCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\tprovide either SITS_T or SITS_T_low/SITS_T_high, "
                "not both\n");
        }
        if (!has_T_list && !has_T_low)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorMissingCommand, "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\tSITS temperature ladder must be provided through "
                "SITS_T or SITS_T_low/SITS_T_high\n");
        }
        if (has_T_low)
        {
            controller->Check_Float(sits->module_name, "T_low",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            controller->Check_Float(sits->module_name, "T_high",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            T_low = atof(controller->Command(sits->module_name, "T_low"));
            T_high = atof(controller->Command(sits->module_name, "T_high"));
            Validate_SITS_Float(controller, sits->module_name, "T_low", T_low,
                                true);
            Validate_SITS_Float(controller, sits->module_name, "T_high", T_high,
                                true);
            if (!(T_high > T_low))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "CLASSIC_SITS_INFORMATION::Initial",
                    "Reason:\n\t%s_T_high (%.9g) must be greater than "
                    "%s_T_low (%.9g)\n",
                    sits->module_name, T_high, sits->module_name, T_low);
            }
            float T_space = (T_high - T_low) / (k_numbers - 1);
            for (int i = 0; i < k_numbers; ++i)
            {
                beta_k_tmp[i] = 1.0 / (CONSTANT_kB * (T_low + T_space * i));
            }
        }
        else
        {
            Parse_SITS_Temperature_List(
                controller, sits->module_name,
                controller->Command(sits->module_name, "T"), k_numbers,
                beta_k_tmp);
        }
        for (int i = 0; i < k_numbers; i++)
        {
            if (!Float_Memory_Is_Finite(beta_k_tmp + i) ||
                !(beta_k_tmp[i] > 0.0f))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "CLASSIC_SITS_INFORMATION::Initial",
                    "Reason:\n\t%s temperature entry %d produces an invalid "
                    "inverse temperature\n",
                    sits->module_name, i);
            }
        }
        deviceMemcpy(beta_k, beta_k_tmp, sizeof(float) * k_numbers,
                     deviceMemcpyHostToDevice);
        free(beta_k_tmp);
        if (controller->Command_Exist(sits->module_name, "record_interval"))
        {
            controller->Check_Int(sits->module_name, "record_interval",
                                  "CLASSIC_SITS_INFORMATION::Initial");
            record_interval =
                atoi(controller->Command(sits->module_name, "record_interval"));
        }
        else
        {
            record_interval = 1;
        }
        if (record_interval <= 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\t%s_record_interval must be positive (received "
                "%d)\n",
                sits->module_name, record_interval);
        }
        controller->printf("    SITS record interval set to %d\n",
                           record_interval);

        if (controller->Command_Exist(sits->module_name, "update_interval"))
        {
            controller->Check_Int(sits->module_name, "update_interval",
                                  "CLASSIC_SITS_INFORMATION::Initial");
            update_interval =
                atoi(controller->Command(sits->module_name, "update_interval"));
        }
        else
        {
            update_interval = 100;
        }
        if (update_interval <= 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "CLASSIC_SITS_INFORMATION::Initial",
                "Reason:\n\t%s_update_interval must be positive (received "
                "%d)\n",
                sits->module_name, update_interval);
        }
        controller->printf("    SITS update interval set to %d\n",
                           update_interval);

        if (controller->Command_Exist(sits->module_name, "pe_a"))
        {
            controller->Check_Float(sits->module_name, "pe_a",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_a = atof(controller->Command(sits->module_name, "pe_a"));
        }
        else
        {
            pe_a = 1.0;
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_a", pe_a, true);
        controller->printf("    SITS_pe_a set to %f\n", pe_a);

        if (controller->Command_Exist(sits->module_name, "pe_b"))
        {
            controller->Check_Float(sits->module_name, "pe_b",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            pe_b = atof(controller->Command(sits->module_name, "pe_b"));
        }
        else
        {
            pe_b = 0.0;
        }
        Validate_SITS_Float(controller, sits->module_name, "pe_b", pe_b, false);
        controller->printf("    SITS_pe_b set to %f\n", pe_b);

        if (controller->Command_Exist(sits->module_name, "fb_bias"))
        {
            controller->Check_Float(sits->module_name, "fb_bias",
                                    "CLASSIC_SITS_INFORMATION::Initial");
            fb_bias = atof(controller->Command(sits->module_name, "fb_bias"));
        }
        else
        {
            fb_bias = 0.0;
        }
        Validate_SITS_Float(controller, sits->module_name, "fb_bias", fb_bias,
                            false);
        controller->printf("    SITS_fb_bias set to %f\n", fb_bias);

        reset = 1;

        int nk_rest;
        if (sits->sits_mode == SITS_MODE_ITERATION)
        {
            nk_rest = 0;
        }
        else
        {
            nk_rest = 1;
        }
        if (controller->Command_Exist(sits->module_name, "nk_rest"))
        {
            nk_rest = controller->Get_Bool(sits->module_name, "nk_rest",
                                           "CLASSIC_SITS_INFORMATION::Initial");
        }
        float* beta_lin;
        Malloc_Safely((void**)&beta_lin, sizeof(float) * k_numbers);

        for (int i = 0; i < k_numbers; ++i) beta_lin[i] = -FLT_MAX;

        deviceMemcpy(log_norm_old, beta_lin, sizeof(float) * k_numbers,
                     deviceMemcpyHostToDevice);
        deviceMemcpy(log_norm, beta_lin, sizeof(float) * k_numbers,
                     deviceMemcpyHostToDevice);
        deviceMemset(log_nk_inverse, 0, sizeof(float) * k_numbers);

        if (nk_rest == 0)
        {
            for (int i = 0; i < k_numbers; ++i)
            {
                beta_lin[i] = 0.0;
            }
        }
        else
        {
            FILE* nk_read_file;
            if (controller->Command_Exist(sits->module_name, "nk_in_file"))
            {
                controller->printf(
                    "    Read Nk from %s\n",
                    controller->Original_Command(sits->module_name,
                                                 "nk_in_file"));
                Open_File_Safely(
                    &nk_read_file,
                    controller->Original_Command(sits->module_name,
                                                 "nk_in_file"),
                    "r");
                for (int i = 0; i < k_numbers; ++i)
                {
                    const int retval = fscanf(nk_read_file, "%f", beta_lin + i);
                    if (retval != 1 || !Float_Memory_Is_Finite(beta_lin + i) ||
                        !(beta_lin[i] > 0.0f))
                    {
                        fclose(nk_read_file);
                        controller->Throw_Formatted_SPONGE_Error(
                            spongeErrorBadFileFormat,
                            "CLASSIC_SITS_INFORMATION::Initial",
                            "Reason:\n\t%s_nk_in_file entry %d must be a "
                            "positive finite number\n",
                            sits->module_name, i);
                    }
                    beta_lin[i] = logf(beta_lin[i]);
                }
                char trailing_token[2] = {0, 0};
                if (fscanf(nk_read_file, " %1s", trailing_token) == 1)
                {
                    fclose(nk_read_file);
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat,
                        "CLASSIC_SITS_INFORMATION::Initial",
                        "Reason:\n\t%s_nk_in_file contains data after its %d "
                        "required entries\n",
                            sits->module_name, k_numbers);
                }
                if (ferror(nk_read_file))
                {
                    fclose(nk_read_file);
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat,
                        "CLASSIC_SITS_INFORMATION::Initial",
                        "Reason:\n\tI/O error while checking the end of "
                        "%s_nk_in_file\n",
                        sits->module_name);
                }
                fclose(nk_read_file);
            }
            else
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorMissingCommand,
                    "CLASSIC_SITS_INFORMATION::Initial",
                    "Reason:\n\tSITS_nk_in_file must be given when "
                    "SITS_nk_rest = 1 or SITS_mode = production\n");
            }
        }
        deviceMemcpy(log_nk, beta_lin, sizeof(float) * k_numbers,
                     deviceMemcpyHostToDevice);

        for (int i = 0; i < k_numbers; ++i)
        {
            beta_lin[i] = -beta_lin[i];
        }
        deviceMemcpy(log_nk_inverse, beta_lin, sizeof(float) * k_numbers,
                     deviceMemcpyHostToDevice);

        for (int i = 0; i < k_numbers; ++i)
        {
            beta_lin[i] = expf(-beta_lin[i]);
        }
        deviceMemcpy(Nk, beta_lin, sizeof(float) * k_numbers,
                     deviceMemcpyHostToDevice);

        free(beta_lin);
        Reset_List(factor, 1.0, 1);

        if (controller->Command_Exist(sits->module_name, "nk_fix"))
        {
            nk_fix = controller->Get_Bool(sits->module_name, "nk_fix",
                                          "CLASSIC_SITS_INFORMATION::Initial");
        }
        else if (sits->sits_mode == SITS_MODE_ITERATION)
        {
            nk_fix = 0;
        }
        else
        {
            nk_fix = 1;
        }
        controller->printf("    SITS nk fix is: %d\n", nk_fix);
        if (nk_fix == 0)
        {
            if (controller->Command_Exist(sits->module_name, "nk_rest_file"))
            {
                nk_rest_file_name = controller->Original_Command(
                    sits->module_name, "nk_rest_file");
            }
            else if (controller->Command_Exist("default_out_file_prefix"))
            {
                nk_rest_file_name =
                    std::string(controller->Original_Command(
                        "default_out_file_prefix")) +
                    "_" + sits->module_name + "_nk_rest.txt";
            }
            else
            {
                nk_rest_file_name =
                    std::string(sits->module_name) + "_nk_rest.txt";
            }
            controller->printf("    Restart Nk will be written in %s\n",
                               nk_rest_file_name.c_str());
            std::string default_name = sits->module_name;
            default_name += "_nk_traj.dat";
            if (CONTROLLER::MPI_rank == 0)
            {
                if (controller->Command_Exist(sits->module_name,
                                              "nk_traj_file"))
                {
                    nk_traj_file_name = controller->Original_Command(
                        sits->module_name, "nk_traj_file");
                }
                else if (controller->Command_Exist("default_out_file_prefix"))
                {
                    nk_traj_file_name =
                        std::string(controller->Original_Command(
                            "default_out_file_prefix")) +
                        "_nk_traj.dat";
                }
                else
                {
                    nk_traj_file_name = default_name;
                }
                nk_traj_file = controller->Get_Output_File(
                    true, sits->module_name, "nk_traj_file", "_nk_traj.dat",
                    default_name.c_str());
            }
        }
    }
}

void CLASSIC_SITS_INFORMATION::Memory_Allocate()
{
    Malloc_Safely((void**)&nk_record_cpu, sizeof(float) * k_numbers);
    Malloc_Safely((void**)&log_norm_record_cpu, sizeof(float) * k_numbers);

    Device_Malloc_Safely((void**)&ene_recorded, sizeof(float));
    Device_Malloc_Safely((void**)&gf, sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&gfsum, sizeof(float));
    Device_Malloc_Safely((void**)&log_weight, sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&log_mk_inverse, sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&log_norm_old, sizeof(float) * k_numbers);
    Device_Malloc_And_Copy_Safely((void**)&log_norm, log_norm_record_cpu,
                                  sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&log_pk, sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&log_nk_inverse, sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&log_nk, sizeof(float) * k_numbers);

    Device_Malloc_Safely((void**)&beta_k, sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&NkExpBetakU, sizeof(float) * k_numbers);
    Device_Malloc_And_Copy_Safely((void**)&Nk, nk_record_cpu,
                                  sizeof(float) * k_numbers);
    Device_Malloc_Safely((void**)&sum_a, sizeof(float));
    Device_Malloc_Safely((void**)&sum_b, sizeof(float));
    Device_Malloc_And_Copy_Safely((void**)&factor, &sits_controller->h_factor,
                                  sizeof(float));
    Device_Malloc_Safely((void**)&fb_reference_energy, sizeof(float));
    Device_Malloc_Safely((void**)&fb_reference_bias, sizeof(float));
    deviceMemset(fb_reference_energy, 0, sizeof(float));
    deviceMemset(fb_reference_bias, 0, sizeof(float));
}

void CLASSIC_SITS_INFORMATION::SITS_Record_Ene()
{
    Launch_Device_Kernel(SITS_Record_Ene_Device, 1, 1, 0, NULL, ene_recorded,
                         sits_controller->pw_select.select_energy[0], pe_a,
                         pe_b);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_gf()
{
    Launch_Device_Kernel(SITS_Update_gf_Device,
                         SITS_Kernel_Block_Count(k_numbers), 64, 0, NULL,
                         k_numbers, gf, ene_recorded, log_nk, beta_k);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_gfsum()
{
    Launch_Device_Kernel(SITS_Update_gfsum_Device, 1, 1, 0, NULL, k_numbers,
                         gfsum, gf);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_log_pk()
{
    Launch_Device_Kernel(SITS_Update_log_pk_Device,
                         SITS_Kernel_Block_Count(k_numbers), 64, 0, NULL,
                         k_numbers, log_pk, gf, gfsum, reset);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_log_mk_inverse()
{
    Launch_Device_Kernel(SITS_Update_log_mk_inverse_Device,
                         SITS_Kernel_Block_Count(k_numbers), 64, 0, NULL,
                         k_numbers,
                         log_weight, log_mk_inverse, log_norm_old, log_norm,
                         log_pk, log_nk);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_log_nk_inverse()
{
    Launch_Device_Kernel(SITS_Update_log_nk_inverse_Device, 1, 1, 0, NULL,
                         k_numbers, log_nk_inverse, log_mk_inverse);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_nk()
{
    Launch_Device_Kernel(SITS_Update_nk_Device,
                         SITS_Kernel_Block_Count(k_numbers), 64, 0, NULL,
                         k_numbers, log_nk, Nk, log_nk_inverse);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_Fb(float beta_0, int step)
{
    if (!is_initialized || sits_controller->sits_mode == SITS_MODE_OBSERVATION)
    {
        return;
    }
    if (step % fb_interval != 0)
    {
        // Holding only the last force factor while clearing the bias makes the
        // reported energy and applied force different Hamiltonians.  Between
        // exact feedback evaluations, use the tangent potential through the
        // last reference point.  Its derivative is precisely the cached
        // factor, so energy, force, and virial remain mutually consistent.
        Launch_Device_Kernel(SITS_Evaluate_Linearized_Fb, 1, 1, 0, NULL,
                             sits_controller->pw_select.select_energy[0],
                             factor, fb_reference_energy, fb_reference_bias,
                             d_bias);
        deviceMemcpy(&sits_controller->h_factor, factor, sizeof(float),
                     deviceMemcpyDeviceToHost);
        return;
    }
    if (sits_controller->sits_mode < SITS_MODE_EMPIRICAL)
    {
        SITS_Get_Current_Fb(sits_controller->atom_numbers,
                            sits_controller->pw_select.select_energy[0], d_bias,
                            k_numbers, NkExpBetakU, beta_k, log_nk, beta_0,
                            sum_a, sum_b, factor, fb_bias, pe_a, pe_b,
                            sits_controller->pwwp_enhance_factor);
    }
    else if (sits_controller->sits_mode == SITS_MODE_EMPIRICAL)
    {
        Launch_Device_Kernel(ESITS_Get_Current_Fb, 1, 1, 0, NULL,
                             sits_controller->pw_select.select_energy[0],
                             factor, pe_a, pe_b,
                             1.0f / (beta_0 * T_low * CONSTANT_kB),
                             1.0f / (beta_0 * T_high * CONSTANT_kB), d_bias);
    }
    else if (sits_controller->sits_mode == SITS_MODE_AMD)
    {
        Launch_Device_Kernel(AMD_Get_Current_Fb, 1, 1, 0, NULL,
                             sits_controller->pw_select.select_energy[0],
                             factor, pe_a, pe_b, d_bias);
    }
    else if (sits_controller->sits_mode == SITS_MODE_GAMD)
    {
        Launch_Device_Kernel(GAMD_Get_Current_Fb, 1, 1, 0, NULL,
                             sits_controller->pw_select.select_energy[0],
                             factor, pe_a, pe_b, d_bias);
    }
    Launch_Device_Kernel(SITS_Record_Fb_Reference, 1, 1, 0, NULL,
                         sits_controller->pw_select.select_energy[0], d_bias,
                         fb_reference_energy, fb_reference_bias);
    deviceMemcpy(&sits_controller->h_factor, factor, sizeof(float),
                 deviceMemcpyDeviceToHost);
}

void CLASSIC_SITS_INFORMATION::SITS_Update_Common(const float beta)
{
    if (sits_controller->sits_mode != SITS_MODE_EMPIRICAL)
    {
        SITS_Record_Ene();
        SITS_Update_gf();
        SITS_Update_gfsum();
        SITS_Update_log_pk();
        reset = 0;
        record_count++;
    }
}

void CLASSIC_SITS_INFORMATION::SITS_Update_Nk()
{
    if (sits_controller->sits_mode != SITS_MODE_EMPIRICAL)
    {
        SITS_Update_log_mk_inverse();
        SITS_Update_log_nk_inverse();
        SITS_Update_nk();

        record_count = 0;
        reset = 1;

        SITS_Write_Nk_Norm();
    }
}

static int SITS_Effective_IO_Error(int error_number)
{
    return error_number == 0 ? EIO : error_number;
}

static bool Open_SITS_Restart_Temporary(const std::string& restart_name,
                                        std::string* temporary_name,
                                        FILE** temporary_file,
                                        int* open_error)
{
    temporary_name->clear();
    *temporary_file = NULL;
    *open_error = 0;
#ifdef _WIN32
    const unsigned long process_id = GetCurrentProcessId();
#else
    const long process_id = static_cast<long>(getpid());
#endif
    // O_EXCL semantics are essential here: a fixed or pre-existing temporary
    // path could be a stale symlink to the live restart and would let fopen
    // truncate the very file this transaction is supposed to preserve.
    for (std::uint64_t attempt = 0;; ++attempt)
    {
        std::string candidate;
        try
        {
            candidate = restart_name + ".sponge-tmp." +
                        std::to_string(process_id) + "." +
                        std::to_string(attempt);
        }
        catch (const std::length_error&)
        {
#ifdef ENAMETOOLONG
            *open_error = ENAMETOOLONG;
#elif defined(EOVERFLOW)
            *open_error = EOVERFLOW;
#else
            *open_error = ERANGE;
#endif
            return false;
        }
        catch (const std::bad_alloc&)
        {
            *open_error = ENOMEM;
            return false;
        }
        errno = 0;
#ifdef _WIN32
        FILE* candidate_file = NULL;
        const errno_t status =
            fopen_s(&candidate_file, candidate.c_str(), "wx");
        const int candidate_error =
            status == 0 ? 0 : static_cast<int>(status);
#else
        FILE* candidate_file = fopen(candidate.c_str(), "wx");
        const int candidate_error = errno;
#endif
        if (candidate_file != NULL)
        {
            *temporary_name = std::move(candidate);
            *temporary_file = candidate_file;
            return true;
        }
        if (candidate_error != EEXIST)
        {
            *open_error = SITS_Effective_IO_Error(candidate_error);
            return false;
        }
        if (attempt == std::numeric_limits<std::uint64_t>::max())
        {
#ifdef EOVERFLOW
            *open_error = EOVERFLOW;
#else
            *open_error = ERANGE;
#endif
            return false;
        }
    }
}

static std::string Remove_SITS_Temporary(const std::string& temporary_name)
{
    errno = 0;
    if (std::remove(temporary_name.c_str()) == 0 || errno == ENOENT)
    {
        return std::string();
    }
    const int cleanup_error = SITS_Effective_IO_Error(errno);
    return std::string("; additionally failed to remove the temporary file: ") +
           strerror(cleanup_error);
}

#ifdef _WIN32
static std::string SITS_Windows_Error_Text(unsigned long error_number)
{
    char buffer[512] = {0};
    const unsigned long length = FormatMessageA(
        FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS, NULL,
        error_number, 0, buffer, static_cast<unsigned long>(sizeof(buffer)),
        NULL);
    if (length == 0)
    {
        return "unknown Windows error";
    }
    std::string result(buffer, length);
    while (!result.empty() &&
           (result.back() == '\r' || result.back() == '\n'))
    {
        result.pop_back();
    }
    return result;
}
#endif

void CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm()
{
#ifdef USE_MPI
    if (CONTROLLER::MPI_rank != 0) return;
#endif
    deviceMemcpy(nk_record_cpu, Nk, sizeof(float) * k_numbers,
                 deviceMemcpyDeviceToHost);
    for (int i = 0; i < k_numbers; ++i)
    {
        if (!Float_Memory_Is_Zero_Or_Normal(&nk_record_cpu[i]))
        {
            sits_controller->controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm",
                "Reason:\n\tSITS Nk entry %d is non-finite or subnormal and "
                "cannot be persisted safely\n",
                i);
            return;
        }
    }
    if (nk_traj_file != NULL)
    {
        const std::size_t expected = static_cast<std::size_t>(k_numbers);
        errno = 0;
        const std::size_t written =
            fwrite(nk_record_cpu, sizeof(float), expected, nk_traj_file);
        const int write_error = errno;
        errno = 0;
        const int flush_status = fflush(nk_traj_file);
        const int flush_error = errno;
        const bool stream_error = ferror(nk_traj_file) != 0;
        if (written != expected || flush_status != 0 || stream_error)
        {
            const char* failed_operation =
                written != expected ? "write" : "flush";
            const int io_error = SITS_Effective_IO_Error(
                written != expected ? write_error : flush_error);
            sits_controller->controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed,
                "CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm",
                "Reason:\n\tfailed to %s the complete SITS Nk trajectory "
                "record in '%s': %s (wrote %zu of %zu float values)\n",
                failed_operation,
                nk_traj_file_name.empty() ? "<unnamed>"
                                          : nk_traj_file_name.c_str(),
                strerror(io_error), written, expected);
            return;
        }
    }

    std::string temporary_name;
    int open_error = 0;
    if (!Open_SITS_Restart_Temporary(nk_rest_file_name, &temporary_name,
                                     &nk_rest_file, &open_error))
    {
        sits_controller->controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed,
            "CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm",
            "Reason:\n\tfailed to create an exclusive temporary SITS "
            "restart next to '%s': %s; the previous restart was preserved\n",
            nk_rest_file_name.c_str(), strerror(open_error));
        return;
    }
    const char* failed_operation = NULL;
    int io_error = 0;
    for (int i = 0; i < k_numbers; ++i)
    {
        errno = 0;
        if (fprintf(nk_rest_file, "%e%c", nk_record_cpu[i],
                    i + 1 == k_numbers ? '\n' : ' ') < 0)
        {
            failed_operation = "write";
            io_error = SITS_Effective_IO_Error(errno);
            break;
        }
    }
    if (failed_operation == NULL)
    {
        errno = 0;
        if (fflush(nk_rest_file) != 0)
        {
            failed_operation = "flush";
            io_error = SITS_Effective_IO_Error(errno);
        }
    }
    if (failed_operation == NULL)
    {
        errno = 0;
#ifdef _WIN32
        if (_commit(_fileno(nk_rest_file)) != 0)
#else
        if (fsync(fileno(nk_rest_file)) != 0)
#endif
        {
            failed_operation = "synchronize";
            io_error = SITS_Effective_IO_Error(errno);
        }
    }
    if (failed_operation == NULL && ferror(nk_rest_file) != 0)
    {
        failed_operation = "write or flush";
        io_error = SITS_Effective_IO_Error(errno);
    }
    errno = 0;
    const int close_status = fclose(nk_rest_file);
    const int close_error = errno;
    nk_rest_file = NULL;
    if (failed_operation == NULL && close_status != 0)
    {
        failed_operation = "close";
        io_error = SITS_Effective_IO_Error(close_error);
    }
    if (failed_operation != NULL)
    {
        const std::string cleanup_detail =
            Remove_SITS_Temporary(temporary_name);
        sits_controller->controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed,
            "CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm",
            "Reason:\n\tfailed to %s the complete temporary SITS restart "
            "'%s': %s; the previous restart '%s' was preserved%s\n",
            failed_operation, temporary_name.c_str(), strerror(io_error),
            nk_rest_file_name.c_str(), cleanup_detail.c_str());
        return;
    }

    errno = 0;
#ifdef _WIN32
    const bool replace_failed =
        MoveFileExA(temporary_name.c_str(), nk_rest_file_name.c_str(),
                    MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) == 0;
    const unsigned long replace_error =
        replace_failed ? GetLastError() : ERROR_SUCCESS;
#else
    const bool replace_failed =
        std::rename(temporary_name.c_str(), nk_rest_file_name.c_str()) != 0;
    const int replace_error = errno;
#endif
    if (replace_failed)
    {
        const std::string cleanup_detail =
            Remove_SITS_Temporary(temporary_name);
#ifdef _WIN32
        const std::string replace_error_text =
            SITS_Windows_Error_Text(replace_error);
        sits_controller->controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed,
            "CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm",
            "Reason:\n\tfailed to atomically replace SITS restart '%s' with "
            "the completed temporary file (Windows error %lu: %s); the "
            "previous restart was preserved%s\n",
            nk_rest_file_name.c_str(), replace_error,
            replace_error_text.c_str(),
            cleanup_detail.c_str());
#else
        sits_controller->controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed,
            "CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm",
            "Reason:\n\tfailed to atomically replace SITS restart '%s' with "
            "the completed temporary file: %s; the previous restart was "
            "preserved%s\n",
            nk_rest_file_name.c_str(), strerror(replace_error),
            cleanup_detail.c_str());
#endif
        return;
    }
}

void SITS_INFORMATION::Initial(CONTROLLER* controller, int atom_numbers_,
                               const char* given_module_name)
{
    this->controller = controller;
    if (given_module_name == NULL)
    {
        strcpy(module_name, "SITS");
        strcpy(print_aa_kab_name, "SITS");
        strcpy(print_bias_name, "SITS");
        strcpy(print_fb_name, "SITS");
    }
    else
    {
        strcpy(module_name, given_module_name);
        strcpy(print_aa_kab_name, given_module_name);
        strcpy(print_bias_name, given_module_name);
        strcpy(print_fb_name, given_module_name);
    }
    strcat(print_aa_kab_name, "_AA_kAB");
    strcat(print_bias_name, "_bias");
    strcat(print_fb_name, "_fb");
    if (controller->Command_Exist(module_name, "mode"))
    {
        if (controller->Command_Choice(module_name, "mode", "observation"))
        {
            controller->printf(
                "START INITIALIZING %s\n    %s mode = observation\n",
                module_name, module_name);
            is_initialized = 1;
            sits_mode = SITS_MODE_OBSERVATION;
        }
        else if (controller->Command_Choice(module_name, "mode", "iteration"))
        {
            controller->printf(
                "START INITIALIZING %s\n    %s mode = iteration\n", module_name,
                module_name);
            is_initialized = 1;
            sits_mode = SITS_MODE_ITERATION;
        }
        else if (controller->Command_Choice(module_name, "mode", "production"))
        {
            controller->printf(
                "START INITIALIZING %s\n    %s mode = production\n",
                module_name, module_name);
            is_initialized = 1;
            sits_mode = SITS_MODE_PRODUCTION;
        }
        else if (controller->Command_Choice(module_name, "mode", "empirical"))
        {
            controller->printf(
                "START INITIALIZING %s\n    %s mode = empirical\n", module_name,
                module_name);
            is_initialized = 1;
            sits_mode = SITS_MODE_EMPIRICAL;
        }
        else if (controller->Command_Choice(module_name, "mode", "amd"))
        {
            controller->printf(
                "START INITIALIZING %s\n    %s mode = AMD (Accelerated MD)\n",
                module_name, module_name);
            is_initialized = 1;
            sits_mode = SITS_MODE_AMD;
        }
        else if (controller->Command_Choice(module_name, "mode", "gamd"))
        {
            controller->printf(
                "START INITIALIZING %s\n    %s mode = GAMD (Gaussian "
                "Accelerated MD)\n",
                module_name, module_name);
            is_initialized = 1;
            sits_mode = SITS_MODE_GAMD;
        }
        else
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "SITS_INFORMATION::Initial",
                "Reason:\n\tSITS mode must be one of observation, iteration, "
                "production, empirical, amd, or gamd\n");
        }
        atom_numbers = atom_numbers_;
        if (atom_numbers <= 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "SITS_INFORMATION::Initial",
                "Reason:\n\tSITS requires at least one atom\n");
        }
        controller->printf("\tAtom numbers is %d\n", atom_numbers);
        Memory_Allocate();

        pw_select.Initial();
        pw_select.Add_One_Energy(atom_numbers);
        pw_select.Add_One_Force(atom_numbers);
        pw_select.Add_One_Virial(atom_numbers);

        if (controller->Command_Exist(module_name, "cross_enhance_factor"))
        {
            controller->Check_Float(module_name, "cross_enhance_factor",
                                    "SITS_INFORMATION::Initial");
            pwwp_enhance_factor =
                atof(controller->Command(module_name, "cross_enhance_factor"));
        }
        else
        {
            pwwp_enhance_factor = 0.5;
        }
        if (!Float_Memory_Is_Finite(&pwwp_enhance_factor))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "SITS_INFORMATION::Initial",
                "Reason:\n\tSITS cross_enhance_factor must be finite\n");
            return;
        }
        controller->printf("\tpwwp enhance factor set to %f\n",
                           pwwp_enhance_factor);

        const bool has_atom_file =
            controller->Command_Exist(module_name, "atom_in_file");
        const bool has_atom_numbers =
            controller->Command_Exist(module_name, "atom_numbers");
        if (has_atom_file == has_atom_numbers)
        {
            controller->Throw_SPONGE_Error(
                has_atom_file ? spongeErrorConflictingCommand
                              : spongeErrorMissingCommand,
                "SITS_INFORMATION::Initial",
                "Reason:\n\texactly one of SITS_atom_in_file and "
                "SITS_atom_numbers must be provided\n");
        }
        this->selectively_applied = true;
        if (has_atom_file || has_atom_numbers)
        {
            controller->printf("    Set atom atribution information\n");
            int* atom_sys_mark_cpu;
            Malloc_Safely((void**)&atom_sys_mark_cpu,
                          sizeof(int) * atom_numbers);
            if (has_atom_file)
            {
                for (int i = 0; i < atom_numbers; i++)
                {
                    atom_sys_mark_cpu[i] = 1;
                }
                controller->printf("    reading %s_atom_in_file\n",
                                   module_name);
                FILE* fr = NULL;
                int temp_atom;
                Open_File_Safely(
                    &fr,
                    controller->Original_Command(module_name, "atom_in_file"),
                    "r");
                int read_status = 0;
                while ((read_status = fscanf(fr, "%d", &temp_atom)) == 1)
                {
                    if (temp_atom < 0 || temp_atom >= atom_numbers)
                    {
                        fclose(fr);
                        free(atom_sys_mark_cpu);
                        controller->Throw_Formatted_SPONGE_Error(
                            spongeErrorBadFileFormat,
                            "SITS_INFORMATION::Initial",
                            "Reason:\n\t%s_atom_in_file contains atom index "
                            "%d outside [0, %d)\n",
                            module_name, temp_atom, atom_numbers);
                    }
                    atom_sys_mark_cpu[temp_atom] = 0;
                }
                if (ferror(fr))
                {
                    fclose(fr);
                    free(atom_sys_mark_cpu);
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "SITS_INFORMATION::Initial",
                        "Reason:\n\tI/O error while reading "
                        "%s_atom_in_file\n",
                        module_name);
                }
                if (read_status != EOF)
                {
                    fclose(fr);
                    free(atom_sys_mark_cpu);
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "SITS_INFORMATION::Initial",
                        "Reason:\n\t%s_atom_in_file contains a non-integer "
                        "token\n",
                        module_name);
                }
                fclose(fr);
            }
            else if (strcmp(controller->Command(module_name, "atom_numbers"),
                            "ITS") == 0 ||
                     strcmp(controller->Command(module_name, "atom_numbers"),
                            "ALL") == 0)
            {
                this->selectively_applied = false;
                for (int i = 0; i < atom_numbers; i++)
                {
                    atom_sys_mark_cpu[i] = 0;
                }
            }
            else
            {
                controller->Check_Int(module_name, "atom_numbers",
                                      "SITS_INFORMATION::Initial");
                int protein_numbers =
                    atoi(controller->Command(module_name, "atom_numbers"));
                if (protein_numbers < 0 || protein_numbers > atom_numbers)
                {
                    free(atom_sys_mark_cpu);
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorValueErrorCommand,
                        "SITS_INFORMATION::Initial",
                        "Reason:\n\t%s_atom_numbers is %d but must be within "
                        "[0, %d]\n",
                        module_name, protein_numbers, atom_numbers);
                }
                for (int i = 0; i < protein_numbers; i++)
                {
                    atom_sys_mark_cpu[i] = 0;
                }
                for (int i = protein_numbers; i < atom_numbers; i++)
                {
                    atom_sys_mark_cpu[i] = 1;
                }
            }
            deviceMemcpy(atom_sys_mark, atom_sys_mark_cpu,
                         sizeof(int) * atom_numbers, deviceMemcpyHostToDevice);
            free(atom_sys_mark_cpu);
        }
        h_factor = 1.0f;
        classic_sits.Initial(controller, this);

        controller->Step_Print_Initial(print_aa_kab_name, "%.2f");
        controller->Step_Print_Initial(print_bias_name, "%.4f");
        controller->Step_Print_Initial(print_fb_name, "%.4f");

        controller->printf("END INTIALIZING %s\n\n", module_name);
    }
    else
    {
        is_initialized = 0;
        return;
    }
}

void SITS_INFORMATION::Memory_Allocate()
{
    Device_Malloc_Safely((void**)&atom_sys_mark, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&atom_sys_mark_local,
                         sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_local_metadata_error, sizeof(int));
}

void SITS_INFORMATION::Reset_Force_Energy(int* md_need_potential)
{
    if (!is_initialized) return;
    md_need_potential[0] += 1;

    deviceMemset(pw_select.select_atom_energy[0], 0,
                 sizeof(float) * atom_numbers);
    deviceMemset(pw_select.select_energy[0], 0, sizeof(float));
    deviceMemset(classic_sits.d_bias, 0, sizeof(float));
    deviceMemset(pw_select.select_force[0], 0, sizeof(VECTOR) * atom_numbers);
    deviceMemset(pw_select.select_atom_virial_tensor[0], 0,
                 sizeof(LTMatrix3) * atom_numbers);
    deviceMemset(pw_select.select_virial_tensor[0], 0, sizeof(LTMatrix3));
}

void SITS_INFORMATION::Update_And_Enhance(
    const int step, float* d_total_potential, int need_pressure,
    LTMatrix3* d_total_virial, VECTOR* frc, float beta0, bool update_statistics)
{
    if (!is_initialized) return;
    if (!Float_Memory_Is_Finite(&beta0) || !(beta0 > 0.0f))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "SITS_INFORMATION::Update_And_Enhance",
            "Reason:\n\tthe inverse target temperature supplied to SITS "
            "must be finite and positive\n");
        return;
    }
    if (selectively_applied)
    {
        Sum_Of_List(pw_select.select_atom_energy[0], pw_select.select_energy[0],
                    atom_numbers);
#ifdef USE_MPI
        if (CONTROLLER::PP_MPI_size != 1)
            D_MPI_Allreduce_IN_PLACE(pw_select.select_energy[0], 1, D_MPI_FLOAT,
                                     D_MPI_SUM, CONTROLLER::d_pp_comm, NULL);
#endif
        if (need_pressure)
        {
            // Keep this rank's virial contribution local.  The pressure path
            // performs the one global stress reduction after every local
            // contribution has been assembled.  Reducing here as well would
            // inject the same global selective virial on every PP rank.
            Sum_Of_List(pw_select.select_atom_virial_tensor[0],
                        pw_select.select_virial_tensor[0], atom_numbers);
        }
    }
    else
    {
        deviceMemcpy(pw_select.select_energy[0], d_total_potential,
                     sizeof(float), deviceMemcpyDeviceToDevice);
        deviceMemcpy(pw_select.select_force[0], frc,
                     sizeof(VECTOR) * local_atom_numbers,
                     deviceMemcpyDeviceToDevice);
        if (need_pressure)
        {
            deviceMemcpy(pw_select.select_virial_tensor[0], d_total_virial,
                         sizeof(LTMatrix3), deviceMemcpyDeviceToDevice);
        }
    }
    float enhancing_energy = 0.0f;
    deviceMemcpy(&enhancing_energy, pw_select.select_energy[0], sizeof(float),
                 deviceMemcpyDeviceToHost);
    if (!Float_Memory_Is_Finite(&enhancing_energy))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "SITS_INFORMATION::Update_And_Enhance",
            "Reason:\n\tthe SITS enhancing energy is non-finite\n");
        return;
    }
    if (update_statistics && sits_mode != SITS_MODE_OBSERVATION &&
        !classic_sits.nk_fix && step % classic_sits.record_interval == 0)
    {
        classic_sits.SITS_Update_Common(beta0);
        if (classic_sits.record_count % classic_sits.update_interval == 0)
        {
            classic_sits.SITS_Update_Nk();
        }
    }
    if (sits_mode != SITS_MODE_OBSERVATION)
    {
        classic_sits.SITS_Update_Fb(beta0, step);
    }
    float bias = 0.0f;
    deviceMemcpy(&bias, classic_sits.d_bias, sizeof(float),
                 deviceMemcpyDeviceToHost);
    if (!Float_Memory_Is_Finite(&h_factor) || !Float_Memory_Is_Finite(&bias))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "SITS_INFORMATION::Update_And_Enhance",
            "Reason:\n\tthe SITS bias or force factor is non-finite\n");
        return;
    }
    if (sits_mode == SITS_MODE_GAMD && h_factor < 0.0f)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "SITS_INFORMATION::Update_And_Enhance",
            "Reason:\n\tthe GaMD force factor is negative "
            "(1 - k * (E - U) = %.9g for k = %.9g, E = %.9g, "
            "U = %.9g); require k * (E - U) <= 1 whenever U < E\n",
            h_factor, classic_sits.pe_a, classic_sits.pe_b, enhancing_energy);
        return;
    }
    // Force buffers use domain-local indexing.  The selection arrays are
    // globally sized only to provide enough storage for every possible local
    // layout; walking atom_numbers here overruns dd.frc on more than one PP
    // rank.  Keep one block for an empty rank so the scalar/tensor update in
    // thread zero still executes.
    const int update_grid_size =
        (local_atom_numbers + CONTROLLER::device_max_thread - 1) /
        CONTROLLER::device_max_thread;
    Launch_Device_Kernel(SITS_For_Enhanced_Force_Protein_Water_Device,
                         update_grid_size > 0 ? update_grid_size : 1,
                         CONTROLLER::device_max_thread, 0, NULL,
                         local_atom_numbers, frc, pw_select.select_force[0],
                         d_total_potential, classic_sits.d_bias, need_pressure,
                         d_total_virial, pw_select.select_virial_tensor[0],
                         h_factor - 1);
}

void SITS_INFORMATION::Save_State(SITS_STATE_SNAPSHOT* snapshot,
                                  const float* d_effective_potential)
{
    if (!is_initialized || snapshot == NULL) return;
    deviceMemcpy(&snapshot->enhancing_energy, pw_select.select_energy[0],
                 sizeof(float), deviceMemcpyDeviceToHost);
    deviceMemcpy(&snapshot->bias, classic_sits.d_bias, sizeof(float),
                 deviceMemcpyDeviceToHost);
    snapshot->factor = h_factor;
    if (sits_mode != SITS_MODE_OBSERVATION)
    {
        deviceMemcpy(&snapshot->fb_reference_energy,
                     classic_sits.fb_reference_energy, sizeof(float),
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(&snapshot->fb_reference_bias,
                     classic_sits.fb_reference_bias, sizeof(float),
                     deviceMemcpyDeviceToHost);
    }
    if (d_effective_potential != NULL)
    {
        deviceMemcpy(&snapshot->effective_potential, d_effective_potential,
                     sizeof(float), deviceMemcpyDeviceToHost);
    }
}

void SITS_INFORMATION::Restore_State(const SITS_STATE_SNAPSHOT& snapshot,
                                     float* d_effective_potential)
{
    if (!is_initialized) return;
    deviceMemcpy(pw_select.select_energy[0], &snapshot.enhancing_energy,
                 sizeof(float), deviceMemcpyHostToDevice);
    deviceMemcpy(classic_sits.d_bias, &snapshot.bias, sizeof(float),
                 deviceMemcpyHostToDevice);
    h_factor = snapshot.factor;
    if (sits_mode != SITS_MODE_OBSERVATION)
    {
        deviceMemcpy(classic_sits.factor, &snapshot.factor, sizeof(float),
                     deviceMemcpyHostToDevice);
        deviceMemcpy(classic_sits.fb_reference_energy,
                     &snapshot.fb_reference_energy, sizeof(float),
                     deviceMemcpyHostToDevice);
        deviceMemcpy(classic_sits.fb_reference_bias,
                     &snapshot.fb_reference_bias, sizeof(float),
                     deviceMemcpyHostToDevice);
    }
    if (d_effective_potential != NULL)
    {
        deviceMemcpy(d_effective_potential, &snapshot.effective_potential,
                     sizeof(float), deviceMemcpyHostToDevice);
    }
}

void SITS_INFORMATION::SITS_LJ_Direct_CF_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const int local_atom_numbers,
    const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
    const float* charge, LENNARD_JONES_INFORMATION* lj_info, VECTOR* md_frc,
    const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
    const float cutoff, const float pme_beta, const int need_potential,
    float* atom_energy, const int need_pressure, LTMatrix3* atom_virial,
    float* coulomb_atom_ene)
{
    if (is_initialized && lj_info->is_initialized)
    {
        if (!Validate_Local_State(
                "SITS_INFORMATION::"
                "SITS_LJ_Direct_CF_Force_With_Atom_Energy_And_Virial",
                atom_numbers, local_atom_numbers, ghost_numbers))
        {
            return;
        }
        if (!lj_info->Validate_Local_State(
                "SITS_INFORMATION::"
                "SITS_LJ_Direct_CF_Force_With_Atom_Energy_And_Virial",
                atom_numbers, local_atom_numbers, ghost_numbers,
                solvent_numbers))
        {
            return;
        }
        Launch_Device_Kernel(
            Repack_LJ_Crd,
            (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL,
            local_atom_numbers + ghost_numbers, crd, lj_info->d_lj_crd_q);

        if (need_potential)
        {
            deviceMemset(coulomb_atom_ene, 0,
                         sizeof(float) * (local_atom_numbers + ghost_numbers));
            deviceMemset(lj_info->d_LJ_energy_atom, 0,
                         sizeof(float) * (local_atom_numbers + ghost_numbers));
            deviceMemset(classic_sits.d_bias, 0, sizeof(float));
        }
        if (!local_atom_numbers) return;
        lj_info->Reset_Pair_Overlap_Error();

        auto f = Selective_Lennard_Jones_And_Direct_Coulomb_Device<true, false,
                                                                   false, true>;
        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        dim3 gridSize = (local_atom_numbers + blockSize.y - 1) / blockSize.y;

        if (need_potential && !need_pressure)
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Device<true, true,
                                                                  false, true>;
        }
        else if (need_potential && need_pressure)
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Device<true, true,
                                                                  true, true>;
        }
        else if (!need_potential && need_pressure)
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Device<true, false,
                                                                  true, true>;
        }
        else
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Device<true, false,
                                                                  false, true>;
        }

        Launch_Device_Kernel(
            f, gridSize, blockSize, 0, NULL, local_atom_numbers,
            solvent_numbers, nl, lj_info->d_LJ_energy_atom, lj_info->d_lj_crd_q,
            lj_info->d_lj_type_g, cell, rcell, lj_info->d_LJ_A, lj_info->d_LJ_B,
            atom_sys_mark_local, cutoff, md_frc, pw_select.select_force[0],
            pme_beta, atom_energy, pw_select.select_atom_energy[0], atom_virial,
            pw_select.select_atom_virial_tensor[0], coulomb_atom_ene,
            pwwp_enhance_factor, lj_info->d_pair_overlap_error);
        lj_info->Check_Pair_Overlap_Error(
            "SITS_INFORMATION::"
            "SITS_LJ_Direct_CF_Force_With_Atom_Energy_And_Virial");
    }
}

void SITS_INFORMATION::
    SITS_LJ_Soft_Core_Direct_CF_Force_With_Atom_Energy_And_Virial(
        const int atom_numbers, const int local_atom_numbers,
        const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
        const float* charge, LJ_SOFT_CORE* lj_info, VECTOR* md_frc,
        const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
        const float cutoff, const float pme_beta, const int need_potential,
        float* atom_energy, const int need_pressure, LTMatrix3* atom_virial,
        float* coulomb_atom_ene)
{
    if (is_initialized && lj_info->is_initialized)
    {
        const char* error_by =
            "SITS_INFORMATION::"
            "SITS_LJ_Soft_Core_Direct_CF_Force_With_Atom_Energy_And_Virial";
        if (!Validate_Local_State(error_by, atom_numbers, local_atom_numbers,
                                  ghost_numbers))
        {
            return;
        }
        if (!lj_info->Validate_Local_State(error_by, atom_numbers,
                                           local_atom_numbers, ghost_numbers,
                                           solvent_numbers))
        {
            return;
        }
        if (!lj_info->Prepare_Local_Coordinates(
                error_by, local_atom_numbers + ghost_numbers, crd, charge))
        {
            return;
        }

        if (need_potential)
        {
            deviceMemset(coulomb_atom_ene, 0,
                         sizeof(float) * (local_atom_numbers + ghost_numbers));
            deviceMemset(lj_info->d_LJ_energy_atom, 0,
                         sizeof(float) * (local_atom_numbers + ghost_numbers));
            deviceMemset(classic_sits.d_bias, 0, sizeof(float));
        }
        if (!local_atom_numbers) return;
        lj_info->Reset_Pair_Overlap_Error();

        auto f = Selective_Lennard_Jones_And_Direct_Coulomb_Soft_Core_Device<
            true, false, false, true, false>;
        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        dim3 gridSize = (local_atom_numbers + blockSize.y - 1) / blockSize.y;

        if (need_potential && !need_pressure)
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Soft_Core_Device<
                true, true, false, true, false>;
        }
        else if (need_potential && need_pressure)
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Soft_Core_Device<
                true, true, true, true, false>;
        }
        else if (!need_potential && need_pressure)
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Soft_Core_Device<
                true, false, true, true, false>;
        }
        else
        {
            f = Selective_Lennard_Jones_And_Direct_Coulomb_Soft_Core_Device<
                true, false, false, true, false>;
        }
        Launch_Device_Kernel(
            f, gridSize, blockSize, 0, NULL, local_atom_numbers,
            solvent_numbers, nl, lj_info->d_LJ_energy_atom,
            lj_info->crd_with_LJ_parameters_local, cell, rcell,
            atom_sys_mark_local, lj_info->d_LJ_AA, lj_info->d_LJ_AB,
            lj_info->d_LJ_BA, lj_info->d_LJ_BB, cutoff, md_frc,
            pw_select.select_force[0], pme_beta, atom_energy,
            pw_select.select_atom_energy[0], atom_virial,
            pw_select.select_atom_virial_tensor[0], coulomb_atom_ene, NULL,
            NULL, NULL, lj_info->lambda, lj_info->alpha, lj_info->p,
            lj_info->sigma_6, lj_info->sigma_6_min, pwwp_enhance_factor,
            lj_info->d_pair_overlap_error);
        lj_info->Check_Pair_Overlap_Error(error_by);
    }
}

void SITS_INFORMATION::Step_Print(CONTROLLER* controller, const float beta0)
{
    if (!is_initialized) return;
    float bias;
    deviceMemcpy(&bias, classic_sits.d_bias, sizeof(float),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&h_enhancing_energy, pw_select.select_energy[0], sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print(print_aa_kab_name, h_enhancing_energy);
    controller->Step_Print(print_bias_name, bias);
    controller->Step_Print(print_fb_name, h_factor);
}

static __global__ void Check_Solvent_Atom_Included(int atom_numbers,
                                                   int solvent_numbers,
                                                   int* atom_sys_mark,
                                                   int* errored)
{
#ifdef USE_GPU
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    if (i < solvent_numbers)
#else
#pragma omp parallel for
    for (int i = 0; i < solvent_numbers; i++)
#endif
    {
        if (!atom_sys_mark[atom_numbers - solvent_numbers + i]) errored[0] = 1;
    }
}

void SITS_INFORMATION::Check_Solvent(CONTROLLER* controller, int atom_numbers,
                                     int solvent_numbers)
{
    if (!is_initialized || !selectively_applied || solvent_numbers == 0) return;
    if (atom_numbers != this->atom_numbers || solvent_numbers < 0 ||
        solvent_numbers > atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "SITS_INFORMATION::Check_Solvent",
            "Reason:\n\tinvalid global/solvent atom counts %d/%d for a "
            "SITS system initialized with %d atoms\n",
            atom_numbers, solvent_numbers, this->atom_numbers);
    }
    int *errored, h_errored;
    Device_Malloc_Safely((void**)&errored, sizeof(int));
    deviceMemset(errored, 0, sizeof(int));
    Launch_Device_Kernel(Check_Solvent_Atom_Included,
                         (solvent_numbers + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
                         solvent_numbers, atom_sys_mark, errored);

    deviceMemcpy(&h_errored, errored, sizeof(int), deviceMemcpyDeviceToHost);
    if (h_errored == 1)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "SITS_INFORMATION::Check_Solvent",
            "Reason:\n\tYou are trying to apply SITS to the solvents. If YOU "
            "KNOW WHAT YOU ARE DOING, set the command 'solvent_LJ' to 0 to run "
            "the simulation.");
    }
    Free_Single_Device_Pointer((void**)&errored);
}

void SELECT::Initial()
{
    select_atom_energy.clear();
    select_energy.clear();
    select_force.clear();
    select_atom_virial_tensor.clear();
    select_virial_tensor.clear();
}

int SELECT::Add_One_Energy(int atom_numbers)
{
    float* tmp_atom_energy;
    float* tmp_energy;
    Device_Malloc_Safely((void**)&tmp_atom_energy,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&tmp_energy, sizeof(float));
    select_atom_energy.push_back(tmp_atom_energy);
    select_energy.push_back(tmp_energy);
    return select_energy.size() - 1;
}

int SELECT::Add_One_Force(int atom_numbers)
{
    VECTOR* tmp_force;
    Device_Malloc_Safely((void**)&tmp_force, sizeof(VECTOR) * atom_numbers);
    select_force.push_back(tmp_force);
    return (select_force.size() - 1);
}

int SELECT::Add_One_Virial(int atom_numbers)
{
    LTMatrix3* tmp_atom_virial_tensor;
    LTMatrix3* tmp_virial_tensor;
    Device_Malloc_Safely((void**)&tmp_atom_virial_tensor,
                         sizeof(LTMatrix3) * atom_numbers);
    Device_Malloc_Safely((void**)&tmp_virial_tensor, sizeof(LTMatrix3));
    select_atom_virial_tensor.push_back(tmp_atom_virial_tensor);
    select_virial_tensor.push_back(tmp_virial_tensor);
    return select_virial_tensor.size() - 1;
}

static __global__ void get_local_device(int* atom_local, int local_atom_numbers,
                                        int ghost_numbers,
                                        int global_atom_numbers,
                                        int* atom_sys_mark,
                                        int* atom_sys_mark_local,
                                        int* invalid_local_index)
{
    int total = local_atom_numbers + ghost_numbers;
    SIMPLE_DEVICE_FOR(i, total)
    {
        const int global_atom = atom_local[i];
        if (global_atom < 0 || global_atom >= global_atom_numbers)
        {
            atomicExch(invalid_local_index, i);
        }
        else
        {
            atom_sys_mark_local[i] = atom_sys_mark[global_atom];
        }
    }
}

void SITS_INFORMATION::Get_Local(int* atom_local, int local_atom_numbers_,
                                 int ghost_numbers_)
{
    if (is_initialized)
    {
        local_metadata_is_ready = false;
        if (local_atom_numbers_ < 0 || ghost_numbers_ < 0 ||
            local_atom_numbers_ > atom_numbers ||
            ghost_numbers_ > atom_numbers - local_atom_numbers_)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "SITS_INFORMATION::Get_Local",
                "Reason:\n\t%s received invalid local/ghost atom counts "
                "%d/%d for %d global atoms\n",
                module_name, local_atom_numbers_, ghost_numbers_, atom_numbers);
            return;
        }
        local_atom_numbers = local_atom_numbers_;
        ghost_numbers = ghost_numbers_;
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
            CONTROLLER::device_max_thread, 0, NULL, atom_local,
            local_atom_numbers, ghost_numbers, atom_numbers, atom_sys_mark,
            atom_sys_mark_local, d_local_metadata_error);
        int invalid_local_index = -1;
        deviceMemcpy(&invalid_local_index, d_local_metadata_error, sizeof(int),
                     deviceMemcpyDeviceToHost);
        if (invalid_local_index >= 0)
        {
            int invalid_global_atom = -1;
            deviceMemcpy(&invalid_global_atom, atom_local + invalid_local_index,
                         sizeof(int), deviceMemcpyDeviceToHost);
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "SITS_INFORMATION::Get_Local",
                "Reason:\n\t%s local coordinate %d maps to global atom %d "
                "outside [0, %d)\n",
                module_name, invalid_local_index, invalid_global_atom,
                atom_numbers);
            return;
        }
        local_metadata_is_ready = true;
    }
}

bool SITS_INFORMATION::Validate_Local_State(const char* error_by,
                                            int global_atom_numbers,
                                            int local_atom_numbers,
                                            int ghost_numbers)
{
    if (global_atom_numbers != atom_numbers ||
        local_atom_numbers != this->local_atom_numbers ||
        ghost_numbers != this->ghost_numbers || !local_metadata_is_ready)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s local selection metadata mismatch: call has "
            "global/local/ghost counts %d/%d/%d, initialized state has "
            "%d/%d/%d and ready=%d\n",
            module_name, global_atom_numbers, local_atom_numbers, ghost_numbers,
            atom_numbers, this->local_atom_numbers, this->ghost_numbers,
            static_cast<int>(local_metadata_is_ready));
        return false;
    }
    return true;
}
