#include "bond_order.h"

#include "reaxff_input.h"

// Evaluate 1 / (1 + exp(x)) without ever forming an overflowing exp(x).
// Besides the value, form the derivative from -f * (1 - f) directly.  The
// generic SAD quotient would otherwise turn the mathematically valid
// saturation limit into 0 * inf / inf = NaN for strongly coordinated atoms.
template <int N>
static __device__ __forceinline__ SADfloat<N>
Stable_Inverse_One_Plus_Exp(const SADfloat<N>& x)
{
    SADfloat<N> result;
    if (x.val >= 0.0f)
    {
        const float exp_negative_x = expf(-x.val);
        result.val = exp_negative_x / (1.0f + exp_negative_x);
    }
    else
    {
        const float exp_x = expf(x.val);
        result.val = 1.0f / (1.0f + exp_x);
    }
    const float slope = -result.val * (1.0f - result.val);
    for (int derivative_i = 0; derivative_i < N; derivative_i++)
    {
        result.dval[derivative_i] = slope * x.dval[derivative_i];
    }
    return result;
}

// f2 = exp(a) + exp(b) appears only through
//
//     (valency + f2) / (valency + f2 + f3).
//
// Work with log(f2) and scale by f2 whenever it is larger than one.  This is
// algebraically identical, while keeping both branches and their SAD
// derivatives representable when the individual exponentials overflow.
template <int N>
static __device__ __forceinline__ SADfloat<N> Stable_BOC_F1_Ratio(
    float valency, const SADfloat<N>& log_f2, const SADfloat<N>& f3)
{
    if (log_f2.val > 0.0f)
    {
        const SADfloat<N> inverse_f2 = expf(-1.0f * log_f2);
        return (1.0f + valency * inverse_f2) /
               (1.0f + (valency + f3) * inverse_f2);
    }
    const SADfloat<N> f2 = expf(log_f2);
    return (valency + f2) / (valency + f2 + f3);
}

// Use neighbor list instead of O(N^2) all-pairs scan
static __global__ void Calculate_Uncorrected_Bond_Orders_Kernel(
    int atom_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, float cutoff, const int* atom_type, const float* r_s,
    const float* r_p, const float* r_pp, const float* bo_1, const float* bo_2,
    const float* bo_3, const float* bo_4, const float* bo_5, const float* bo_6,
    const float* ro_pi, const float* ro_pi2, const int atom_type_numbers,
    float bo_cut, float* total_bond_order, const ATOM_GROUP* nl, int* pair_i,
    int* pair_j, float* distances, int max_pairs, unsigned long long* num_pairs,
    int* geometry_error)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type_i = atom_type[i];
        if (type_i < 0 || type_i >= atom_type_numbers)
        {
            Record_ReaxFF_Geometry_Error(
                geometry_error, REAXFF_INVALID_ATOM_TYPE, i, -1, type_i);
#ifdef GPU_ARCH_NAME
            return;
#else
            continue;
#endif
        }
        VECTOR ri = crd[i];
        ATOM_GROUP nl_i = nl[i];

        for (int nn = 0; nn < nl_i.atom_numbers; nn++)
        {
            int j = nl_i.atom_serial[nn];
            if (j <= i) continue;  // only process each pair once

            int type_j = atom_type[j];
            if (type_j < 0 || type_j >= atom_type_numbers)
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_INVALID_ATOM_TYPE, i, j,
                                             type_i, type_j);
                continue;
            }

            VECTOR rj = crd[j];
            VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
            float r2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;

            if (!ReaxFF_Vector_Is_Finite(drij) || !ReaxFF_Float_Is_Finite(r2))
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_BOND_NONFINITE, i, j);
                continue;
            }

            if (r2 < cutoff * cutoff)
            {
                int idx = type_i * atom_type_numbers + type_j;

                const bool sigma_active = r_s[idx] > 0.0f;
                const bool pi_active = ro_pi[type_i] > 0.0f &&
                                       ro_pi[type_j] > 0.0f && r_p[idx] > 0.0f;
                const bool pi2_active = ro_pi2[type_i] > 0.0f &&
                                        ro_pi2[type_j] > 0.0f &&
                                        r_pp[idx] > 0.0f;
                if (r2 == 0.0f)
                {
                    if (sigma_active || pi_active || pi2_active)
                    {
                        Record_ReaxFF_Geometry_Error(geometry_error,
                                                     REAXFF_BOND_OVERLAP, i, j);
                    }
                    continue;
                }

                float r = sqrtf(r2);
                if (!ReaxFF_Float_Is_Finite(r) || !(r > 0.0f))
                {
                    Record_ReaxFF_Geometry_Error(geometry_error,
                                                 REAXFF_BOND_NONFINITE, i, j);
                    continue;
                }

                float ros = r_s[idx];
                float bo_s = 0.0f;
                if (ros > 0.0f)
                {
                    float C12 = bo_1[idx] * powf(r / ros, bo_2[idx]);
                    bo_s = (1.0f + bo_cut) * expf(C12);
                }

                float bo_p = 0.0f;
                if (ro_pi[type_i] > 0.0f && ro_pi[type_j] > 0.0f)
                {
                    float rop = r_p[idx];
                    if (rop > 0.0f)
                    {
                        float C34 = bo_3[idx] * powf(r / rop, bo_4[idx]);
                        bo_p = expf(C34);
                    }
                }

                float bo_p2 = 0.0f;
                if (ro_pi2[type_i] > 0.0f && ro_pi2[type_j] > 0.0f)
                {
                    float rop2 = r_pp[idx];
                    if (rop2 > 0.0f)
                    {
                        float C56 = bo_5[idx] * powf(r / rop2, bo_6[idx]);
                        bo_p2 = expf(C56);
                    }
                }

                float total_bo = bo_s + bo_p + bo_p2;

                if (!ReaxFF_Float_Is_Finite(bo_s) ||
                    !ReaxFF_Float_Is_Finite(bo_p) ||
                    !ReaxFF_Float_Is_Finite(bo_p2) ||
                    !ReaxFF_Float_Is_Finite(total_bo))
                {
                    Record_ReaxFF_Geometry_Error(geometry_error,
                                                 REAXFF_BOND_NONFINITE, i, j);
                    continue;
                }

                if (total_bo >= bo_cut)
                {
                    bo_s -= bo_cut;
                    if (bo_s < 0.0f) bo_s = 0.0f;
                    total_bo -= bo_cut;

                    atomicAdd(&total_bond_order[i], total_bo);
                    atomicAdd(&total_bond_order[j], total_bo);

                    unsigned long long pos = atomicAdd(num_pairs, 1ULL);
                    if (pos < static_cast<unsigned long long>(max_pairs))
                    {
                        pair_i[static_cast<int>(pos)] = i;
                        pair_j[static_cast<int>(pos)] = j;
                        distances[static_cast<int>(pos)] = r;
                    }
                }
            }
        }
    }
}

// Writes corrected BO and derivatives directly to sparse per-bond arrays
static __global__ void Apply_Bond_Order_Corrections_Kernel(
    int num_pairs, int* pair_i, int* pair_j, float* distances,
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    const int* atom_type, const float* r_s, const float* r_p, const float* r_pp,
    const float* bo_1, const float* bo_2, const float* bo_3, const float* bo_4,
    const float* bo_5, const float* bo_6, const float* ro_pi,
    const float* ro_pi2, const float* valency, const float* valency_val,
    const float* ovc, const float* v13cor, const float* p_boc3,
    const float* p_boc4, const float* p_boc5, const int atom_type_numbers,
    const int atom_numbers, float gp_boc1, float gp_boc2, float bo_cut,
    const float* total_bond_order, float* corrected_bo_s,
    float* corrected_bo_pi, float* corrected_bo_pi2, float* dbo_s_dr,
    float* dbo_pi_dr, float* dbo_pi2_dr, float* dbo_s_dDelta_i,
    float* dbo_pi_dDelta_i, float* dbo_pi2_dDelta_i, float* dbo_s_dDelta_j,
    float* dbo_pi_dDelta_j, float* dbo_pi2_dDelta_j, float* dbo_raw_total_dr,
    int* geometry_error)
{
    SIMPLE_DEVICE_FOR(idx, num_pairs)
    {
        int i = pair_i[idx];
        int j = pair_j[idx];
        float r = distances[idx];

        int type_i = atom_type[i];
        int type_j = atom_type[j];

        if (type_i < 0 || type_i >= atom_type_numbers || type_j < 0 ||
            type_j >= atom_type_numbers)
        {
            Record_ReaxFF_Geometry_Error(
                geometry_error, REAXFF_INVALID_ATOM_TYPE, i, j, type_i, type_j);
#ifdef GPU_ARCH_NAME
            return;
#else
            continue;
#endif
        }

        if (type_i >= 0 && type_i < atom_type_numbers && type_j >= 0 &&
            type_j < atom_type_numbers)
        {
            int pair_idx = type_i * atom_type_numbers + type_j;

            float ros = r_s[pair_idx];
            float bo_s_raw_val = 0.0f, dbo_s_raw_dr = 0.0f;
            if (ros > 0.0f)
            {
                float ratio = r / ros;
                float pow_ratio = powf(ratio, bo_2[pair_idx]);
                bo_s_raw_val =
                    (1.0f + bo_cut) * expf(bo_1[pair_idx] * pow_ratio);
                dbo_s_raw_dr = bo_s_raw_val * bo_1[pair_idx] * bo_2[pair_idx] *
                               powf(ratio, bo_2[pair_idx] - 1.0f) *
                               (1.0f / ros);
            }

            float bo_p_val = 0.0f, dbo_p_raw_dr = 0.0f;
            if (ro_pi[type_i] > 0.0f && ro_pi[type_j] > 0.0f)
            {
                float rop = r_p[pair_idx];
                if (rop > 0.0f)
                {
                    float ratio = r / rop;
                    float pow_ratio = powf(ratio, bo_4[pair_idx]);
                    bo_p_val = expf(bo_3[pair_idx] * pow_ratio);
                    dbo_p_raw_dr = bo_p_val * bo_3[pair_idx] * bo_4[pair_idx] *
                                   powf(ratio, bo_4[pair_idx] - 1.0f) *
                                   (1.0f / rop);
                }
            }

            float bo_p2_val = 0.0f, dbo_p2_raw_dr = 0.0f;
            if (ro_pi2[type_i] > 0.0f && ro_pi2[type_j] > 0.0f)
            {
                float rop2 = r_pp[pair_idx];
                if (rop2 > 0.0f)
                {
                    float ratio = r / rop2;
                    float pow_ratio = powf(ratio, bo_6[pair_idx]);
                    bo_p2_val = expf(bo_5[pair_idx] * pow_ratio);
                    dbo_p2_raw_dr =
                        bo_p2_val * bo_5[pair_idx] * bo_6[pair_idx] *
                        powf(ratio, bo_6[pair_idx] - 1.0f) * (1.0f / rop2);
                }
            }

            float total_bo_raw = bo_s_raw_val + bo_p_val + bo_p2_val;
            float dbo_raw_total_dr_val =
                dbo_s_raw_dr + dbo_p_raw_dr + dbo_p2_raw_dr;

            if (!ReaxFF_Float_Is_Finite(total_bo_raw) ||
                !ReaxFF_Float_Is_Finite(dbo_raw_total_dr_val) ||
                !ReaxFF_Float_Is_Finite(total_bond_order[i]) ||
                !ReaxFF_Float_Is_Finite(total_bond_order[j]))
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_BOND_NONFINITE, i, j);
#ifdef GPU_ARCH_NAME
                return;
#else
                continue;
#endif
            }

            if (total_bo_raw >= bo_cut)
            {
                SADfloat<5> bo_s_raw(bo_s_raw_val, 0);
                SADfloat<5> bo_p(bo_p_val, 1);
                SADfloat<5> bo_p2(bo_p2_val, 2);
                SADfloat<5> Delta_i(total_bond_order[i], 3);
                SADfloat<5> Delta_j(total_bond_order[j], 4);

                SADfloat<5> total_bo_orig = (bo_s_raw + bo_p + bo_p2) - bo_cut;
                SADfloat<5> bo_s = bo_s_raw - bo_cut;
                if (bo_s.val < 0) bo_s = SADfloat<5>(0.0f);

                float ovc_val = ovc[pair_idx];
                float v13cor_val = v13cor[pair_idx];

                SADfloat<5> f1(1.0f);
                if (ovc_val >= 0.001f)
                {
                    SADfloat<5> Deltap_i = Delta_i - valency[type_i];
                    SADfloat<5> Deltap_j = Delta_j - valency[type_j];

                    SADfloat<5> log_f2 = Log_Sum_Exp(
                        -gp_boc1 * Deltap_i, -gp_boc1 * Deltap_j);
                    SADfloat<5> f3 =
                        -1.0f / gp_boc2 *
                        (Log_Sum_Exp(-gp_boc2 * Deltap_i, -gp_boc2 * Deltap_j) -
                         0.6931471805599453f);

                    float val_i = valency[type_i];
                    float val_j = valency[type_j];

                    f1 = 0.5f *
                         (Stable_BOC_F1_Ratio(val_i, log_f2, f3) +
                          Stable_BOC_F1_Ratio(val_j, log_f2, f3));
                }

                SADfloat<5> f4(1.0f), f5(1.0f);
                if (v13cor_val >= 0.001f)
                {
                    SADfloat<5> Deltap_boc_i = Delta_i - valency_val[type_i];
                    SADfloat<5> Deltap_boc_j = Delta_j - valency_val[type_j];

                    float p_boc3_val = p_boc3[pair_idx];
                    float p_boc4_val = p_boc4[pair_idx];
                    float p_boc5_val = p_boc5[pair_idx];

                    const SADfloat<5> exponent_f4 =
                        -(p_boc4_val * total_bo_orig * total_bo_orig -
                          Deltap_boc_i) *
                            p_boc3_val +
                        p_boc5_val;
                    const SADfloat<5> exponent_f5 =
                        -(p_boc4_val * total_bo_orig * total_bo_orig -
                          Deltap_boc_j) *
                            p_boc3_val +
                        p_boc5_val;

                    f4 = Stable_Inverse_One_Plus_Exp(exponent_f4);
                    f5 = Stable_Inverse_One_Plus_Exp(exponent_f5);
                }

                SADfloat<5> A0 = f1 * f4 * f5;

                SADfloat<5> s_corrected_bo_pi = bo_p * A0 * f1;
                SADfloat<5> s_corrected_bo_pi2 = bo_p2 * A0 * f1;
                SADfloat<5> s_corrected_bo_s =
                    total_bo_orig * A0 -
                    (s_corrected_bo_pi + s_corrected_bo_pi2);
                if (s_corrected_bo_s.val < 0)
                    s_corrected_bo_s = SADfloat<5>(0.0f);

                const float corrected_s = s_corrected_bo_s.val;
                const float corrected_pi = s_corrected_bo_pi.val;
                const float corrected_pi2 = s_corrected_bo_pi2.val;
                const float derivative_s =
                    s_corrected_bo_s.dval[0] * dbo_s_raw_dr +
                    s_corrected_bo_s.dval[1] * dbo_p_raw_dr +
                    s_corrected_bo_s.dval[2] * dbo_p2_raw_dr;
                const float derivative_pi =
                    s_corrected_bo_pi.dval[0] * dbo_s_raw_dr +
                    s_corrected_bo_pi.dval[1] * dbo_p_raw_dr +
                    s_corrected_bo_pi.dval[2] * dbo_p2_raw_dr;
                const float derivative_pi2 =
                    s_corrected_bo_pi2.dval[0] * dbo_s_raw_dr +
                    s_corrected_bo_pi2.dval[1] * dbo_p_raw_dr +
                    s_corrected_bo_pi2.dval[2] * dbo_p2_raw_dr;
                const float results[] = {corrected_s,
                                         corrected_pi,
                                         corrected_pi2,
                                         derivative_s,
                                         derivative_pi,
                                         derivative_pi2,
                                         s_corrected_bo_s.dval[3],
                                         s_corrected_bo_pi.dval[3],
                                         s_corrected_bo_pi2.dval[3],
                                         s_corrected_bo_s.dval[4],
                                         s_corrected_bo_pi.dval[4],
                                         s_corrected_bo_pi2.dval[4],
                                         dbo_raw_total_dr_val};
                bool valid = true;
                for (int result_i = 0;
                     result_i <
                     static_cast<int>(sizeof(results) / sizeof(results[0]));
                     result_i++)
                {
                    valid = valid && ReaxFF_Float_Is_Finite(results[result_i]);
                }
                if (!valid)
                {
                    Record_ReaxFF_Geometry_Error(geometry_error,
                                                 REAXFF_BOND_NONFINITE, i, j);
#ifdef GPU_ARCH_NAME
                    return;
#else
                    continue;
#endif
                }

                // Write only after the complete pair state is known valid.
                corrected_bo_s[idx] = corrected_s;
                corrected_bo_pi[idx] = corrected_pi;
                corrected_bo_pi2[idx] = corrected_pi2;
                dbo_s_dr[idx] = derivative_s;
                dbo_pi_dr[idx] = derivative_pi;
                dbo_pi2_dr[idx] = derivative_pi2;
                dbo_s_dDelta_i[idx] = s_corrected_bo_s.dval[3];
                dbo_pi_dDelta_i[idx] = s_corrected_bo_pi.dval[3];
                dbo_pi2_dDelta_i[idx] = s_corrected_bo_pi2.dval[3];
                dbo_s_dDelta_j[idx] = s_corrected_bo_s.dval[4];
                dbo_pi_dDelta_j[idx] = s_corrected_bo_pi.dval[4];
                dbo_pi2_dDelta_j[idx] = s_corrected_bo_pi2.dval[4];
                dbo_raw_total_dr[idx] = dbo_raw_total_dr_val;
            }
            else
            {
                corrected_bo_s[idx] = 0.0f;
                corrected_bo_pi[idx] = 0.0f;
                corrected_bo_pi2[idx] = 0.0f;

                dbo_s_dr[idx] = 0.0f;
                dbo_pi_dr[idx] = 0.0f;
                dbo_pi2_dr[idx] = 0.0f;

                dbo_s_dDelta_i[idx] = 0.0f;
                dbo_pi_dDelta_i[idx] = 0.0f;
                dbo_pi2_dDelta_i[idx] = 0.0f;

                dbo_s_dDelta_j[idx] = 0.0f;
                dbo_pi_dDelta_j[idx] = 0.0f;
                dbo_pi2_dDelta_j[idx] = 0.0f;

                dbo_raw_total_dr[idx] = 0.0f;
            }
        }
    }
}

// Reduce corrected bond orders per atom using CSR bond list
static __global__ void Reduce_Total_Corrected_Bond_Order_Kernel(
    int atom_numbers, const int* bond_count, const int* bond_offset,
    const int* bond_idx, const float* bo_s, const float* bo_pi,
    const float* bo_pi2, float* total_bo)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        float sum = 0.0f;
        int start = bond_offset[i];
        int end = start + bond_count[i];
        for (int k = start; k < end; k++)
        {
            int b = bond_idx[k];
            sum += bo_s[b] + bo_pi[b] + bo_pi2[b];
        }
        total_bo[i] = sum;
    }
}

// --- CSR build kernels ---
static __global__ void Count_Bonds_Per_Atom_Kernel(int num_bonds,
                                                   const int* bond_i,
                                                   const int* bond_j,
                                                   int* bond_count)
{
    SIMPLE_DEVICE_FOR(b, num_bonds)
    {
        atomicAdd(&bond_count[bond_i[b]], 1);
        atomicAdd(&bond_count[bond_j[b]], 1);
    }
}

static __global__ void Exclusive_Prefix_Sum_Kernel(int n, const int* input,
                                                   int* output)
{
    // Simple sequential prefix sum (launched with 1 thread)
    SIMPLE_DEVICE_FOR(dummy, 1)
    {
        output[0] = 0;
        for (int i = 0; i < n; i++)
        {
            output[i + 1] = output[i] + input[i];
        }
    }
}

static __global__ void Fill_Bond_CSR_Kernel(int num_bonds, const int* bond_i,
                                            const int* bond_j,
                                            const int* bond_offset,
                                            int* fill_count, int* bond_nbr,
                                            int* bond_idx)
{
    SIMPLE_DEVICE_FOR(b, num_bonds)
    {
        int i = bond_i[b];
        int j = bond_j[b];
        // Entry for atom i: neighbor is j
        int pos_i = bond_offset[i] + atomicAdd(&fill_count[i], 1);
        bond_nbr[pos_i] = j;
        bond_idx[pos_i] = b;
        // Entry for atom j: neighbor is i
        int pos_j = bond_offset[j] + atomicAdd(&fill_count[j], 1);
        bond_nbr[pos_j] = i;
        bond_idx[pos_j] = b;
    }
}

// --- Force projection kernels (sparse) ---

// In sparse mode, dE_dBO is accumulated to a single bond index by all
// consumers, so no need to sum [i*N+j] + [j*N+i].
static __global__ void Calculate_CdDelta_Prime_Kernel(
    int num_pairs, const int* pair_i, const int* pair_j, const float* dE_dBO_s,
    const float* dE_dBO_pi, const float* dE_dBO_pi2, const float* CdDelta,
    const float* dbo_s_dDelta_i, const float* dbo_pi_dDelta_i,
    const float* dbo_pi2_dDelta_i, const float* dbo_s_dDelta_j,
    const float* dbo_pi_dDelta_j, const float* dbo_pi2_dDelta_j,
    float* CdDelta_prime)
{
    SIMPLE_DEVICE_FOR(idx, num_pairs)
    {
        int i = pair_i[idx];
        int j = pair_j[idx];

        float de_dbo_s_total = dE_dBO_s[idx];
        float de_dbo_pi_total = dE_dBO_pi[idx];
        float de_dbo_pi2_total = dE_dBO_pi2[idx];

        float eff_cdd = CdDelta[i] + CdDelta[j];

        float term_i = (de_dbo_s_total + eff_cdd) * dbo_s_dDelta_i[idx] +
                       (de_dbo_pi_total + eff_cdd) * dbo_pi_dDelta_i[idx] +
                       (de_dbo_pi2_total + eff_cdd) * dbo_pi2_dDelta_i[idx];
        atomicAdd(&CdDelta_prime[i], term_i);

        float term_j = (de_dbo_s_total + eff_cdd) * dbo_s_dDelta_j[idx] +
                       (de_dbo_pi_total + eff_cdd) * dbo_pi_dDelta_j[idx] +
                       (de_dbo_pi2_total + eff_cdd) * dbo_pi2_dDelta_j[idx];
        atomicAdd(&CdDelta_prime[j], term_j);
    }
}

static __global__ void REAXFF_Force_Projection_Kernel(
    int num_pairs, const int* pair_i, const int* pair_j, const float* distances,
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    const float* dE_dBO_s, const float* dE_dBO_pi, const float* dE_dBO_pi2,
    const float* CdDelta, const float* dbo_s_dr, const float* dbo_pi_dr,
    const float* dbo_pi2_dr, const float* dbo_raw_total_dr,
    const float* CdDelta_prime, VECTOR* frc, LTMatrix3* atom_virial,
    int* geometry_error)
{
    SIMPLE_DEVICE_FOR(idx, num_pairs)
    {
        int i = pair_i[idx];
        int j = pair_j[idx];
        float r_val = distances[idx];
        if (!(r_val > 0.0f) || !ReaxFF_Float_Is_Finite(r_val))
        {
            Record_ReaxFF_Geometry_Error(
                geometry_error,
                r_val == 0.0f ? REAXFF_BOND_OVERLAP : REAXFF_BOND_NONFINITE, i,
                j);
        }
        else
        {
            float de_dbo_s_total = dE_dBO_s[idx];
            float de_dbo_pi_total = dE_dBO_pi[idx];
            float de_dbo_pi2_total = dE_dBO_pi2[idx];

            float eff_cdd = CdDelta[i] + CdDelta[j];

            float de_dr = (de_dbo_s_total + eff_cdd) * dbo_s_dr[idx] +
                          (de_dbo_pi_total + eff_cdd) * dbo_pi_dr[idx] +
                          (de_dbo_pi2_total + eff_cdd) * dbo_pi2_dr[idx];

            de_dr +=
                (CdDelta_prime[i] + CdDelta_prime[j]) * dbo_raw_total_dr[idx];

            float force_mag = -de_dr;

            VECTOR ri = crd[i];
            VECTOR rj = crd[j];
            VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);

            float fx = force_mag * drij.x / r_val;
            float fy = force_mag * drij.y / r_val;
            float fz = force_mag * drij.z / r_val;

            if (!ReaxFF_Float_Is_Finite(de_dr) ||
                !ReaxFF_Float_Is_Finite(force_mag) ||
                !ReaxFF_Vector_Is_Finite(drij) || !ReaxFF_Float_Is_Finite(fx) ||
                !ReaxFF_Float_Is_Finite(fy) || !ReaxFF_Float_Is_Finite(fz))
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_BOND_NONFINITE, i, j);
#ifdef GPU_ARCH_NAME
                return;
#else
                continue;
#endif
            }

            atomicAdd(&frc[i].x, fx);
            atomicAdd(&frc[i].y, fy);
            atomicAdd(&frc[i].z, fz);
            atomicAdd(&frc[j].x, -fx);
            atomicAdd(&frc[j].y, -fy);
            atomicAdd(&frc[j].z, -fz);

            if (atom_virial)
            {
                VECTOR fij = {fx, fy, fz};
                atomicAdd(atom_virial + i,
                          Get_Virial_From_Force_Dis(fij, drij));
            }
        }
    }
}

// ============================================================
// Implementation
// ============================================================

void REAXFF_BOND_ORDER::Initial(CONTROLLER* controller, int atom_numbers,
                                const char* parameter_in_file,
                                const char* type_in_file, const float cutoff,
                                float* cutoff_full)
{
    if (parameter_in_file == NULL || type_in_file == NULL) return;

    controller->printf("START INITIALIZING REAXFF_BOND_ORDER\n");
    REAXFF_INPUT_ERROR input_error;
    REAXFF_FORCE_FIELD_IR force_field;
    if (!ReaxFF_Parse_Force_Field_File(parameter_in_file, &force_field,
                                       &input_error))
    {
        const std::string reason = input_error.Describe();
        controller->Throw_Formatted_SPONGE_Error(
            input_error.kind == REAXFF_INPUT_OPEN_ERROR
                ? spongeErrorOpenFileFailed
                : spongeErrorBadFileFormat,
            "REAXFF_BOND_ORDER::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }
    std::vector<int> atom_type;
    if (!ReaxFF_Parse_Type_File_Path(type_in_file, atom_numbers, force_field,
                                     &atom_type, NULL, &input_error))
    {
        const std::string reason = input_error.Describe();
        controller->Throw_Formatted_SPONGE_Error(
            input_error.kind == REAXFF_INPUT_OPEN_ERROR
                ? spongeErrorOpenFileFailed
                : spongeErrorBadFileFormat,
            "REAXFF_BOND_ORDER::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }
    if (force_field.general_parameters.size() < 2)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_BOND_ORDER::Initial",
            "Reason:\n\tforce field %s does not contain the two required "
            "bond-order general parameters",
            parameter_in_file);
        return;
    }

    const float staged_gp_boc1 = force_field.general_parameters[0];
    const float staged_gp_boc2 = force_field.general_parameters[1];
    float staged_gp_bo_cut = 0.001f;
    if (force_field.general_parameters.size() > 29)
        staged_gp_bo_cut = 0.01f * force_field.general_parameters[29];
    if (!Float_Memory_Is_Finite(&staged_gp_bo_cut) ||
        staged_gp_bo_cut < 0.0f)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_BOND_ORDER::Initial",
            "Reason:\n\tbond-order cutoff must remain finite and "
            "non-negative in file %s",
            parameter_in_file);
        return;
    }

    const int n_atom_types = static_cast<int>(force_field.atom_types.size());
    int pair_parameter_count = 0;
    if (!ReaxFF_Checked_Dense_Table_Count(
            n_atom_types, 2, sizeof(float), &pair_parameter_count))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_BOND_ORDER::Initial",
            "Reason:\n\tatom type count %d exceeds the supported bond-order "
            "pair-table extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    std::vector<float> ro_sigma(n_atom_types);
    std::vector<float> ro_pi(n_atom_types);
    std::vector<float> ro_pi2(n_atom_types);
    std::vector<float> valency(n_atom_types);
    std::vector<float> valency_val(n_atom_types);
    std::vector<float> b_o_131(n_atom_types);
    std::vector<float> b_o_132(n_atom_types);
    std::vector<float> b_o_133(n_atom_types);
    for (int i = 0; i < n_atom_types; i++)
    {
        const REAXFF_ATOM_TYPE_IR& atom = force_field.atom_types[i];
        ro_sigma[i] = atom.values[0][0];
        valency[i] = atom.values[0][1];
        ro_pi[i] = atom.values[0][6];
        ro_pi2[i] = atom.values[2][0];
        b_o_131[i] = atom.values[2][3];
        b_o_132[i] = atom.values[2][4];
        b_o_133[i] = atom.values[2][5];
        valency_val[i] = atom.values[3][3];
    }

    std::vector<float> bo_1(pair_parameter_count, 0.0f);
    std::vector<float> bo_2(pair_parameter_count, 0.0f);
    std::vector<float> bo_3(pair_parameter_count, 0.0f);
    std::vector<float> bo_4(pair_parameter_count, 0.0f);
    std::vector<float> bo_5(pair_parameter_count, 0.0f);
    std::vector<float> bo_6(pair_parameter_count, 0.0f);
    std::vector<float> ovc(pair_parameter_count, 0.0f);
    std::vector<float> v13cor(pair_parameter_count, 0.0f);
    std::vector<float> p_boc3(pair_parameter_count, 0.0f);
    std::vector<float> p_boc4(pair_parameter_count, 0.0f);
    std::vector<float> p_boc5(pair_parameter_count, 0.0f);
    std::vector<float> r_s(pair_parameter_count, 0.0f);
    std::vector<float> r_p(pair_parameter_count, 0.0f);
    std::vector<float> r_pp(pair_parameter_count, 0.0f);

    auto geometric_mix = [](float lhs, float rhs)
    {
        const double product =
            static_cast<double>(lhs) * static_cast<double>(rhs);
        return static_cast<float>(sqrt(product));
    };
    auto arithmetic_mix = [](float lhs, float rhs)
    {
        return static_cast<float>(
            0.5 * (static_cast<double>(lhs) + static_cast<double>(rhs)));
    };
    for (int i = 0; i < n_atom_types; i++)
    {
        for (int j = 0; j < n_atom_types; j++)
        {
            const int idx = i * n_atom_types + j;
            p_boc3[idx] = geometric_mix(b_o_132[i], b_o_132[j]);
            p_boc4[idx] = geometric_mix(b_o_131[i], b_o_131[j]);
            p_boc5[idx] = geometric_mix(b_o_133[i], b_o_133[j]);
            r_s[idx] = arithmetic_mix(ro_sigma[i], ro_sigma[j]);
            r_p[idx] = arithmetic_mix(ro_pi[i], ro_pi[j]);
            r_pp[idx] = arithmetic_mix(ro_pi2[i], ro_pi2[j]);
            if (!Float_Memory_Is_Finite(&p_boc3[idx]) ||
                !Float_Memory_Is_Finite(&p_boc4[idx]) ||
                !Float_Memory_Is_Finite(&p_boc5[idx]) ||
                !Float_Memory_Is_Finite(&r_s[idx]) ||
                !Float_Memory_Is_Finite(&r_p[idx]) ||
                !Float_Memory_Is_Finite(&r_pp[idx]))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "REAXFF_BOND_ORDER::Initial",
                    "Reason:\n\tatom parameters produce a non-finite "
                    "bond-order pair parameter for types %d/%d in file %s",
                    i + 1, j + 1, parameter_in_file);
                return;
            }
        }
    }

    for (const REAXFF_BOND_IR& bond : force_field.bonds)
    {
        const int idx1 = bond.type1 - 1;
        const int idx2 = bond.type2 - 1;
        const int forward = idx1 * n_atom_types + idx2;
        const int reverse = idx2 * n_atom_types + idx1;
        bo_1[forward] = bo_1[reverse] = bond.line2[4];
        bo_2[forward] = bo_2[reverse] = bond.line2[5];
        bo_3[forward] = bo_3[reverse] = bond.line2[1];
        bo_4[forward] = bo_4[reverse] = bond.line2[2];
        bo_5[forward] = bo_5[reverse] = bond.line1[4];
        bo_6[forward] = bo_6[reverse] = bond.line1[6];
        ovc[forward] = ovc[reverse] = bond.line2[6];
        v13cor[forward] = v13cor[reverse] = bond.line1[5];
    }
    for (const REAXFF_OFF_DIAGONAL_IR& entry : force_field.off_diagonal)
    {
        const int idx1 = entry.type1 - 1;
        const int idx2 = entry.type2 - 1;
        const int forward = idx1 * n_atom_types + idx2;
        const int reverse = idx2 * n_atom_types + idx1;
        if (entry.values[3] > 0.0f)
            r_s[forward] = r_s[reverse] = entry.values[3];
        if (entry.values[4] > 0.0f)
            r_p[forward] = r_p[reverse] = entry.values[4];
        if (entry.values[5] > 0.0f)
            r_pp[forward] = r_pp[reverse] = entry.values[5];
    }
    for (int pair = 0; pair < pair_parameter_count; pair++)
    {
        if (!Float_Memory_Is_Finite(&r_s[pair]) ||
            !Float_Memory_Is_Finite(&r_p[pair]) ||
            !Float_Memory_Is_Finite(&r_pp[pair]))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "REAXFF_BOND_ORDER::Initial",
                "Reason:\n\tfinal bond-order radii for pair index %d are "
                "non-finite in file %s",
                pair + 1, parameter_in_file);
            return;
        }
    }

    constexpr int max_representable_bonds = INT_MAX / 2;
    long long initial_capacity =
        std::max(1LL, static_cast<long long>(atom_numbers) * 32LL);
    if (controller->Command_Exist("REAXFF", "initial_bond_capacity"))
    {
        controller->Check_Int("REAXFF", "initial_bond_capacity",
                              "REAXFF_BOND_ORDER::Initial");
        initial_capacity =
            atoll(controller->Command("REAXFF", "initial_bond_capacity"));
        if (initial_capacity <= 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "REAXFF_BOND_ORDER::Initial",
                "REAXFF.initial_bond_capacity must be positive, but got %lld.",
                initial_capacity);
            return;
        }
    }
    if (initial_capacity > max_representable_bonds)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "REAXFF_BOND_ORDER::Initial",
            "Initial bond capacity %lld exceeds the CSR representation limit "
            "%d for %d atoms.",
            initial_capacity, max_representable_bonds, atom_numbers);
        return;
    }
    const int staged_max_bonds = static_cast<int>(initial_capacity);

    auto allocate_float_copy = [](const std::vector<float>& values)
    {
        float* result = NULL;
        Malloc_Safely((void**)&result, sizeof(float) * values.size());
        memcpy(result, values.data(), sizeof(float) * values.size());
        return result;
    };
    float* staged_ro_sigma = allocate_float_copy(ro_sigma);
    float* staged_ro_pi = allocate_float_copy(ro_pi);
    float* staged_ro_pi2 = allocate_float_copy(ro_pi2);
    float* staged_valency = allocate_float_copy(valency);
    float* staged_valency_val = allocate_float_copy(valency_val);
    float* staged_b_o_131 = allocate_float_copy(b_o_131);
    float* staged_b_o_132 = allocate_float_copy(b_o_132);
    float* staged_b_o_133 = allocate_float_copy(b_o_133);
    float* staged_bo_1 = allocate_float_copy(bo_1);
    float* staged_bo_2 = allocate_float_copy(bo_2);
    float* staged_bo_3 = allocate_float_copy(bo_3);
    float* staged_bo_4 = allocate_float_copy(bo_4);
    float* staged_bo_5 = allocate_float_copy(bo_5);
    float* staged_bo_6 = allocate_float_copy(bo_6);
    float* staged_ovc = allocate_float_copy(ovc);
    float* staged_v13cor = allocate_float_copy(v13cor);
    float* staged_p_boc3 = allocate_float_copy(p_boc3);
    float* staged_p_boc4 = allocate_float_copy(p_boc4);
    float* staged_p_boc5 = allocate_float_copy(p_boc5);
    float* staged_r_s = allocate_float_copy(r_s);
    float* staged_r_p = allocate_float_copy(r_p);
    float* staged_r_pp = allocate_float_copy(r_pp);
    int* staged_atom_type = NULL;
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);

    this->controller = controller;
    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    gp_boc1 = staged_gp_boc1;
    gp_boc2 = staged_gp_boc2;
    gp_bo_cut = staged_gp_bo_cut;
    max_bonds = staged_max_bonds;
    h_ro_sigma = staged_ro_sigma;
    h_ro_pi = staged_ro_pi;
    h_ro_pi2 = staged_ro_pi2;
    h_valency = staged_valency;
    h_valency_val = staged_valency_val;
    h_b_o_131 = staged_b_o_131;
    h_b_o_132 = staged_b_o_132;
    h_b_o_133 = staged_b_o_133;
    h_bo_1 = staged_bo_1;
    h_bo_2 = staged_bo_2;
    h_bo_3 = staged_bo_3;
    h_bo_4 = staged_bo_4;
    h_bo_5 = staged_bo_5;
    h_bo_6 = staged_bo_6;
    h_ovc = staged_ovc;
    h_v13cor = staged_v13cor;
    h_p_boc3 = staged_p_boc3;
    h_p_boc4 = staged_p_boc4;
    h_p_boc5 = staged_p_boc5;
    h_r_s = staged_r_s;
    h_r_p = staged_r_p;
    h_r_pp = staged_r_pp;
    h_atom_type = staged_atom_type;

    Device_Malloc_And_Copy_Safely((void**)&d_ro_sigma, h_ro_sigma,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_ro_pi, h_ro_pi,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_ro_pi2, h_ro_pi2,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_bo_1, h_bo_1,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_bo_2, h_bo_2,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_bo_3, h_bo_3,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_bo_4, h_bo_4,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_bo_5, h_bo_5,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_bo_6, h_bo_6,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_r_s, h_r_s,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_r_p, h_r_p,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_r_pp, h_r_pp,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_valency, h_valency,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_valency_val, h_valency_val,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_ovc, h_ovc,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_v13cor, h_v13cor,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_p_boc3, h_p_boc3,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_p_boc4, h_p_boc4,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_p_boc5, h_p_boc5,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);

    Device_Malloc_Safely((void**)&d_total_bond_order,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_total_corrected_bond_order,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_CdDelta_prime,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_bond_count, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_bond_offset,
                         sizeof(int) * (atom_numbers + 1));
    Device_Malloc_Safely((void**)&d_fill_count, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_geometry_error,
                         sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);
    Device_Malloc_Safely((void**)&d_num_pairs_ptr, sizeof(unsigned long long));
    Allocate_Bond_Storage(max_bonds);

    is_initialized = 1;
    controller->printf("  Sparse bond storage: max_bonds = %d\n", max_bonds);
    controller->printf("END INITIALIZING REAXFF_BOND_ORDER\n\n");
}
void REAXFF_BOND_ORDER::Allocate_Bond_Storage(int capacity)
{
    const size_t unsigned_capacity =
        capacity > 0 ? static_cast<size_t>(capacity) : 0;
    if (unsigned_capacity == 0 ||
        unsigned_capacity >
            std::numeric_limits<size_t>::max() / (2 * sizeof(int)) ||
        unsigned_capacity > std::numeric_limits<size_t>::max() / sizeof(float))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "REAXFF_BOND_ORDER::Allocate_Bond_Storage",
            "Bond capacity %d cannot be represented as allocation byte "
            "counts on this platform.",
            capacity);
    }
    Device_Malloc_Safely((void**)&d_bond_nbr, sizeof(int) * 2 * capacity);
    Device_Malloc_Safely((void**)&d_bond_idx, sizeof(int) * 2 * capacity);

    Device_Malloc_Safely((void**)&d_pair_i, sizeof(int) * capacity);
    Device_Malloc_Safely((void**)&d_pair_j, sizeof(int) * capacity);
    Device_Malloc_Safely((void**)&d_pair_distances, sizeof(float) * capacity);

    Device_Malloc_Safely((void**)&d_corrected_bo_s, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_corrected_bo_pi, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_corrected_bo_pi2, sizeof(float) * capacity);

    Device_Malloc_Safely((void**)&d_dE_dBO_s, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dE_dBO_pi, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dE_dBO_pi2, sizeof(float) * capacity);

    Device_Malloc_Safely((void**)&d_dbo_s_dr, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_pi_dr, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_pi2_dr, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_s_dDelta_i, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_pi_dDelta_i, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_pi2_dDelta_i, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_s_dDelta_j, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_pi_dDelta_j, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_pi2_dDelta_j, sizeof(float) * capacity);
    Device_Malloc_Safely((void**)&d_dbo_raw_total_dr, sizeof(float) * capacity);
}

void REAXFF_BOND_ORDER::Release_Bond_Storage()
{
    auto release = [](void** pointer)
    {
        if (pointer != NULL && pointer[0] != NULL)
        {
            Free_Single_Device_Pointer(pointer);
        }
    };

    release((void**)&d_bond_nbr);
    release((void**)&d_bond_idx);
    release((void**)&d_pair_i);
    release((void**)&d_pair_j);
    release((void**)&d_pair_distances);
    release((void**)&d_corrected_bo_s);
    release((void**)&d_corrected_bo_pi);
    release((void**)&d_corrected_bo_pi2);
    release((void**)&d_dE_dBO_s);
    release((void**)&d_dE_dBO_pi);
    release((void**)&d_dE_dBO_pi2);
    release((void**)&d_dbo_s_dr);
    release((void**)&d_dbo_pi_dr);
    release((void**)&d_dbo_pi2_dr);
    release((void**)&d_dbo_s_dDelta_i);
    release((void**)&d_dbo_pi_dDelta_i);
    release((void**)&d_dbo_pi2_dDelta_i);
    release((void**)&d_dbo_s_dDelta_j);
    release((void**)&d_dbo_pi_dDelta_j);
    release((void**)&d_dbo_pi2_dDelta_j);
    release((void**)&d_dbo_raw_total_dr);
}

void REAXFF_BOND_ORDER::Ensure_Bond_Capacity(
    unsigned long long required_capacity)
{
    if (required_capacity <= static_cast<unsigned long long>(max_bonds)) return;

    constexpr int max_representable_bonds = INT_MAX / 2;
    if (required_capacity >
        static_cast<unsigned long long>(max_representable_bonds))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "REAXFF_BOND_ORDER::Ensure_Bond_Capacity",
            "Required bond count %llu exceeds the signed-int CSR limit %d.",
            required_capacity, max_representable_bonds);
    }

    const long long geometric_capacity =
        static_cast<long long>(max_bonds) +
        std::max(1LL, static_cast<long long>(max_bonds) / 2LL);
    const unsigned long long next_capacity =
        std::max(required_capacity,
                 static_cast<unsigned long long>(std::min(
                     geometric_capacity,
                     static_cast<long long>(max_representable_bonds))));

    controller->printf(
        "  Expanding sparse bond storage from %d to %llu for %llu bonds\n",
        max_bonds, next_capacity, required_capacity);
    Release_Bond_Storage();
    max_bonds = static_cast<int>(next_capacity);
    Allocate_Bond_Storage(max_bonds);
}

void REAXFF_BOND_ORDER::Calculate_Uncorrected_Bond_Orders_GPU(
    int atom_numbers, const VECTOR* d_crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, float cutoff, const ATOM_GROUP* d_nl, int* d_pair_i,
    int* d_pair_j, float* d_distances, unsigned long long* d_num_pairs_ptr)
{
    if (!is_initialized) return;

    dim3 blockSize = {CONTROLLER::device_max_thread};
    dim3 gridSize = {(atom_numbers + blockSize.x - 1) / blockSize.x};

    deviceMemset(d_total_bond_order, 0, sizeof(float) * atom_numbers);

    unsigned long long h_num_pairs = 0;
    deviceMemcpy(d_num_pairs_ptr, &h_num_pairs, sizeof(h_num_pairs),
                 deviceMemcpyHostToDevice);

    Launch_Device_Kernel(
        Calculate_Uncorrected_Bond_Orders_Kernel, gridSize, blockSize, 0, NULL,
        atom_numbers, d_crd, cell, rcell, cutoff, d_atom_type, d_r_s, d_r_p,
        d_r_pp, d_bo_1, d_bo_2, d_bo_3, d_bo_4, d_bo_5, d_bo_6, d_ro_pi,
        d_ro_pi2, atom_type_numbers, gp_bo_cut, d_total_bond_order, d_nl,
        d_pair_i, d_pair_j, d_distances, max_bonds, d_num_pairs_ptr,
        d_geometry_error);
}

void REAXFF_BOND_ORDER::Calculate_Corrected_Bond_Orders_GPU(
    int atom_numbers, const VECTOR* d_crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, float cutoff, int num_pairs, int* d_pair_i,
    int* d_pair_j, float* d_distances)
{
    if (!is_initialized) return;

    if (num_pairs <= 0) return;

    dim3 blockSize = {CONTROLLER::device_max_thread};
    dim3 gridSize = {(num_pairs + blockSize.x - 1) / blockSize.x};

    Launch_Device_Kernel(
        Apply_Bond_Order_Corrections_Kernel, gridSize, blockSize, 0, NULL,
        num_pairs, d_pair_i, d_pair_j, d_distances, d_crd, cell, rcell,
        d_atom_type, d_r_s, d_r_p, d_r_pp, d_bo_1, d_bo_2, d_bo_3, d_bo_4,
        d_bo_5, d_bo_6, d_ro_pi, d_ro_pi2, d_valency, d_valency_val, d_ovc,
        d_v13cor, d_p_boc3, d_p_boc4, d_p_boc5, atom_type_numbers, atom_numbers,
        gp_boc1, gp_boc2, gp_bo_cut, d_total_bond_order, d_corrected_bo_s,
        d_corrected_bo_pi, d_corrected_bo_pi2, d_dbo_s_dr, d_dbo_pi_dr,
        d_dbo_pi2_dr, d_dbo_s_dDelta_i, d_dbo_pi_dDelta_i, d_dbo_pi2_dDelta_i,
        d_dbo_s_dDelta_j, d_dbo_pi_dDelta_j, d_dbo_pi2_dDelta_j,
        d_dbo_raw_total_dr, d_geometry_error);
}

void REAXFF_BOND_ORDER::Build_Bond_CSR(int atom_numbers, int num_bonds)
{
    if (num_bonds <= 0)
    {
        deviceMemset(d_bond_count, 0, sizeof(int) * atom_numbers);
        deviceMemset(d_bond_offset, 0, sizeof(int) * (atom_numbers + 1));
        return;
    }

    dim3 blockSize = {CONTROLLER::device_max_thread};
    dim3 gridSize_bonds = {(num_bonds + blockSize.x - 1) / blockSize.x};

    // Phase 1: Count bonds per atom
    deviceMemset(d_bond_count, 0, sizeof(int) * atom_numbers);
    Launch_Device_Kernel(Count_Bonds_Per_Atom_Kernel, gridSize_bonds, blockSize,
                         0, NULL, num_bonds, d_pair_i, d_pair_j, d_bond_count);

    // Phase 2: Exclusive prefix sum
    Launch_Device_Kernel(Exclusive_Prefix_Sum_Kernel, dim3(1), dim3(1), 0, NULL,
                         atom_numbers, d_bond_count, d_bond_offset);

    // Phase 3: Fill CSR
    deviceMemset(d_fill_count, 0, sizeof(int) * atom_numbers);
    Launch_Device_Kernel(Fill_Bond_CSR_Kernel, gridSize_bonds, blockSize, 0,
                         NULL, num_bonds, d_pair_i, d_pair_j, d_bond_offset,
                         d_fill_count, d_bond_nbr, d_bond_idx);
}

void REAXFF_BOND_ORDER::Calculate_Corrected_Bond_Order(
    int atom_numbers, const VECTOR* d_crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* fnl_d_nl, float cutoff)
{
    if (!is_initialized) return;

    if (h_atom_type == NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF_BOND_ORDER::Calculate_Corrected_Bond_Order",
            "The host atom-type table is NULL; bond-order evaluation cannot "
            "produce a valid Hamiltonian.");
        return;
    }

    unsigned long long required_pairs = 0;
    while (true)
    {
        Reset_Geometry_Error();
        Calculate_Uncorrected_Bond_Orders_GPU(
            atom_numbers, d_crd, cell, rcell, cutoff, fnl_d_nl, d_pair_i,
            d_pair_j, d_pair_distances, d_num_pairs_ptr);
        Check_Geometry_Error(
            "REAXFF_BOND_ORDER::Calculate_Corrected_Bond_Order");
        deviceMemcpy(&required_pairs, d_num_pairs_ptr, sizeof(required_pairs),
                     deviceMemcpyDeviceToHost);

        if (required_pairs <= static_cast<unsigned long long>(max_bonds)) break;

        Ensure_Bond_Capacity(required_pairs);
        // The overflowing pass accumulated all raw bond orders but stored only
        // the old capacity.  Re-run from a cleared accumulator so pair storage
        // and per-atom totals describe the same complete interaction set.
    }

    h_num_pairs = static_cast<int>(required_pairs);
    const int num_pairs = h_num_pairs;

    Build_Bond_CSR(atom_numbers, num_pairs);
    deviceMemset(d_total_corrected_bond_order, 0, sizeof(float) * atom_numbers);
    if (num_pairs > 0)
    {
        Calculate_Corrected_Bond_Orders_GPU(atom_numbers, d_crd, cell, rcell,
                                            cutoff, num_pairs, d_pair_i,
                                            d_pair_j, d_pair_distances);
        Check_Geometry_Error(
            "REAXFF_BOND_ORDER::Calculate_Corrected_Bond_Order");

        // Reduce corrected BO per atom using CSR
        dim3 blockSize = {CONTROLLER::device_max_thread};
        dim3 gridSize = {(atom_numbers + blockSize.x - 1) / blockSize.x};
        Launch_Device_Kernel(Reduce_Total_Corrected_Bond_Order_Kernel, gridSize,
                             blockSize, 0, NULL, atom_numbers, d_bond_count,
                             d_bond_offset, d_bond_idx, d_corrected_bo_s,
                             d_corrected_bo_pi, d_corrected_bo_pi2,
                             d_total_corrected_bond_order);
    }
}

void REAXFF_BOND_ORDER::Calculate_Forces(int atom_numbers, const VECTOR* d_crd,
                                         VECTOR* d_frc, const LTMatrix3 cell,
                                         const LTMatrix3 rcell, float cutoff,
                                         float* d_CdDelta, int need_virial,
                                         LTMatrix3* atom_virial)
{
    if (!is_initialized || h_num_pairs <= 0) return;

    dim3 blockSize = {CONTROLLER::device_max_thread};
    dim3 gridSize = {(h_num_pairs + blockSize.x - 1) / blockSize.x};

    Reset_Geometry_Error();

    Launch_Device_Kernel(Calculate_CdDelta_Prime_Kernel, gridSize, blockSize, 0,
                         NULL, h_num_pairs, d_pair_i, d_pair_j, d_dE_dBO_s,
                         d_dE_dBO_pi, d_dE_dBO_pi2, d_CdDelta, d_dbo_s_dDelta_i,
                         d_dbo_pi_dDelta_i, d_dbo_pi2_dDelta_i,
                         d_dbo_s_dDelta_j, d_dbo_pi_dDelta_j,
                         d_dbo_pi2_dDelta_j, d_CdDelta_prime);

    Launch_Device_Kernel(
        REAXFF_Force_Projection_Kernel, gridSize, blockSize, 0, NULL,
        h_num_pairs, d_pair_i, d_pair_j, d_pair_distances, d_crd, cell, rcell,
        d_dE_dBO_s, d_dE_dBO_pi, d_dE_dBO_pi2, d_CdDelta, d_dbo_s_dr,
        d_dbo_pi_dr, d_dbo_pi2_dr, d_dbo_raw_total_dr, d_CdDelta_prime, d_frc,
        need_virial ? atom_virial : NULL, d_geometry_error);
    Check_Geometry_Error("REAXFF_BOND_ORDER::Calculate_Forces");
}

void REAXFF_BOND_ORDER::Reset_Geometry_Error()
{
    deviceMemset(d_geometry_error, 0, sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);
}

void REAXFF_BOND_ORDER::Check_Geometry_Error(const char* error_by)
{
    int error[REAXFF_GEOMETRY_ERROR_SIZE] = {0, -1, -1, -1, -1};
    deviceMemcpy(error, d_geometry_error, sizeof(error),
                 deviceMemcpyDeviceToHost);
    if (error[0] == REAXFF_GEOMETRY_OK) return;

    if (error[0] == REAXFF_INVALID_ATOM_TYPE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Invalid ReaxFF atom type while evaluating local atoms %d %d "
            "(types %d %d; valid range is [0, %d)).",
            error[1], error[2], error[3], error[4], atom_type_numbers);
        return;
    }

    const char* reason =
        error[0] == REAXFF_BOND_OVERLAP
            ? "overlap exactly while at least one ReaxFF bond-order component "
              "is active; the radial force direction is undefined"
            : "produce a non-finite or unrepresentable ReaxFF bond order, "
              "derivative, displacement, force, or virial";
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown, error_by,
        "Reason:\n\tlocal atoms %d %d %s\n", error[1], error[2], reason);
}

void REAXFF_BOND_ORDER::Clear_Derivatives(int atom_numbers, float* d_CdDelta)
{
    if (!is_initialized) return;
    if (h_num_pairs > 0)
    {
        deviceMemset(d_dE_dBO_s, 0, sizeof(float) * h_num_pairs);
        deviceMemset(d_dE_dBO_pi, 0, sizeof(float) * h_num_pairs);
        deviceMemset(d_dE_dBO_pi2, 0, sizeof(float) * h_num_pairs);
    }
    if (d_CdDelta)
    {
        deviceMemset(d_CdDelta, 0, sizeof(float) * atom_numbers);
    }
    deviceMemset(d_CdDelta_prime, 0, sizeof(float) * atom_numbers);
}

void REAXFF_BOND_ORDER::Calculate_Bond_Order(
    int atom_numbers, const VECTOR* d_crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* fnl_d_nl, float cutoff)
{
    Calculate_Corrected_Bond_Order(atom_numbers, d_crd, cell, rcell, fnl_d_nl,
                                   cutoff);
}
