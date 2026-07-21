#include "torsion.h"

#include "bond_order.h"  // for find_bond_index
#include "reaxff_geometry.h"
#include "reaxff_input.h"

static __global__ void Calculate_Torsion_Kernel(
    int atom_numbers, const VECTOR* crd, const int* atom_type,
    const float p_tor2, const float p_tor3, const float p_tor4,
    const float p_cot2, const float thb_cut, const float* Delta_boc,
    const REAXFF_TORSION_Info* torsion_info,
    const REAXFF_TORSION_Entry* torsion_entries, int atom_type_numbers,
    const float* bo_s, const float* bo_pi, const float* bo_pi2,
    float* d_dE_dBO_s, float* d_dE_dBO_pi, float* d_dE_dBO_pi2, float* CdDelta,
    const LTMatrix3 cell, const LTMatrix3 rcell, float* atom_energy,
    VECTOR* frc, LTMatrix3* atom_virial, float* d_energy_tor_sum,
    float* d_energy_cot_sum, const int* bond_count, const int* bond_offset,
    const int* bond_nbr, const int* bond_idx_arr, int* geometry_error)
{
    SIMPLE_DEVICE_FOR(j, atom_numbers)
    {
        int type_j = atom_type[j];
        if (type_j >= 0 && type_j < atom_type_numbers)
        {
            VECTOR rj = crd[j];
            float delta_j_val = Delta_boc[j];

            double en_tor = 0.0;
            double en_cot = 0.0;

            int bc_j = bond_count[j];
            int bo_j = bond_offset[j];

            for (int pk = 0; pk < bc_j; pk++)
            {
                int b_jk = bond_idx_arr[bo_j + pk];
                int k = bond_nbr[bo_j + pk];
                if (j >= k) continue;

                int type_k = atom_type[k];
                if (type_k < 0 || type_k >= atom_type_numbers)
                {
                    Record_ReaxFF_Geometry_Error(
                        geometry_error, REAXFF_INVALID_ATOM_TYPE, j, k,
                        type_j, type_k);
                    continue;
                }
                float bo_jk_val = bo_s[b_jk] + bo_pi[b_jk] + bo_pi2[b_jk];
                if (bo_jk_val <= thb_cut) continue;
                float bo_jk_pi_val = bo_pi[b_jk];

                VECTOR rk = crd[k];
                VECTOR djk = Get_Periodic_Displacement(rj, rk, cell, rcell);
                float r_jk = norm3df(djk.x, djk.y, djk.z);
                float delta_k_val = Delta_boc[k];

                int bc_k = bond_count[k];
                int bo_k = bond_offset[k];

                for (int pi = 0; pi < bc_j; pi++)
                {
                    int b_ij = bond_idx_arr[bo_j + pi];
                    int i = bond_nbr[bo_j + pi];
                    if (i == k) continue;
                    int type_i = atom_type[i];
                    if (type_i < 0 || type_i >= atom_type_numbers)
                    {
                        Record_ReaxFF_Geometry_Error(
                            geometry_error, REAXFF_INVALID_ATOM_TYPE, i, j,
                            type_i, type_j);
                        continue;
                    }
                    float bo_ij_val = bo_s[b_ij] + bo_pi[b_ij] + bo_pi2[b_ij];
                    if (bo_ij_val <= thb_cut) continue;

                    VECTOR ri = crd[i];
                    VECTOR dji = Get_Periodic_Displacement(rj, ri, cell, rcell);
                    float r_ij = norm3df(dji.x, dji.y, dji.z);

                    for (int pl = 0; pl < bc_k; pl++)
                    {
                        int b_kl = bond_idx_arr[bo_k + pl];
                        int l = bond_nbr[bo_k + pl];
                        if (l == j || l == i) continue;
                        int type_l = atom_type[l];
                        if (type_l < 0 || type_l >= atom_type_numbers)
                        {
                            Record_ReaxFF_Geometry_Error(
                                geometry_error, REAXFF_INVALID_ATOM_TYPE, i, j,
                                k, l);
                            continue;
                        }
                        float bo_kl_val =
                            bo_s[b_kl] + bo_pi[b_kl] + bo_pi2[b_kl];
                        if (bo_kl_val <= thb_cut) continue;
                        if (bo_ij_val * bo_jk_val * bo_kl_val <= thb_cut)
                            continue;

                        VECTOR rl = crd[l];
                        VECTOR dkl =
                            Get_Periodic_Displacement(rk, rl, cell, rcell);
                        float r_kl = norm3df(dkl.x, dkl.y, dkl.z);

                        int quartet_idx =
                            (((type_i * atom_type_numbers + type_j) *
                                  atom_type_numbers +
                              type_k) *
                                 atom_type_numbers +
                             type_l);
                        REAXFF_TORSION_Info info = torsion_info[quartet_idx];
                        if (info.entry_count <= 0) continue;

                        REAXFF_ANGLE_GEOMETRY angle_ijk;
                        REAXFF_ANGLE_GEOMETRY angle_jkl;
                        const VECTOR dkj = {-djk.x, -djk.y, -djk.z};
                        const bool geometry_valid =
                            Compute_ReaxFF_Angle_Geometry(dji, djk,
                                                         &angle_ijk) &&
                            Compute_ReaxFF_Angle_Geometry(dkj, dkl,
                                                         &angle_jkl);
                        bool direction_dependent = false;
                        for (int e = 0; e < info.entry_count; e++)
                        {
                            const REAXFF_TORSION_Entry* param =
                                &torsion_entries[info.start_idx + e];
                            direction_dependent =
                                direction_dependent || param->V1 != 0.0f ||
                                param->V2 != 0.0f || param->V3 != 0.0f ||
                                param->p_cot1 != 0.0f;
                        }
                        if (!geometry_valid || angle_ijk.is_collinear ||
                            angle_jkl.is_collinear)
                        {
                            // A parameter entry with no torsional or
                            // conjugation coefficient is the exact inactive
                            // branch and remains well defined at collinearity.
                            if (direction_dependent)
                            {
                                Record_ReaxFF_Geometry_Error(
                                    geometry_error, REAXFF_TORSION_UNDEFINED,
                                    i, j, k, l);
                            }
                            continue;
                        }
                        const float sin_ijk = angle_ijk.sine;
                        const float sin_jkl = angle_jkl.sine;

                        float unnorm_cos_phi =
                            -(dji.x * djk.x + dji.y * djk.y + dji.z * djk.z) *
                                (djk.x * dkl.x + djk.y * dkl.y +
                                 djk.z * dkl.z) +
                            (r_jk * r_jk) *
                                (dji.x * dkl.x + dji.y * dkl.y + dji.z * dkl.z);
                        VECTOR cross_jk_kl = {
                            (float)(djk.y * dkl.z - djk.z * dkl.y),
                            (float)(djk.z * dkl.x - djk.x * dkl.z),
                            (float)(djk.x * dkl.y - djk.y * dkl.x)};
                        float unnorm_sin_phi = -r_jk * (dji.x * cross_jk_kl.x +
                                                        dji.y * cross_jk_kl.y +
                                                        dji.z * cross_jk_kl.z);
                        float phi = atan2f(unnorm_sin_phi, unnorm_cos_phi);
                        float cos_phi = cosf(phi);
                        if (!ReaxFF_Float_Is_Finite(r_ij) ||
                            !ReaxFF_Float_Is_Finite(r_jk) ||
                            !ReaxFF_Float_Is_Finite(r_kl) ||
                            !ReaxFF_Float_Is_Finite(phi) ||
                            !ReaxFF_Float_Is_Finite(cos_phi) ||
                            (unnorm_sin_phi == 0.0f &&
                             unnorm_cos_phi == 0.0f))
                        {
                            Record_ReaxFF_Geometry_Error(
                                geometry_error, REAXFF_TORSION_NONFINITE, i, j,
                                k, l);
                            continue;
                        }

                        SADvector<9> s_dji;
                        s_dji.x = SADfloat<9>(dji.x, 0);
                        s_dji.y = SADfloat<9>(dji.y, 1);
                        s_dji.z = SADfloat<9>(dji.z, 2);
                        SADvector<9> s_djk;
                        s_djk.x = SADfloat<9>(djk.x, 3);
                        s_djk.y = SADfloat<9>(djk.y, 4);
                        s_djk.z = SADfloat<9>(djk.z, 5);
                        SADvector<9> s_dkl;
                        s_dkl.x = SADfloat<9>(dkl.x, 6);
                        s_dkl.y = SADfloat<9>(dkl.y, 7);
                        s_dkl.z = SADfloat<9>(dkl.z, 8);

                        SADfloat<9> s_r_ij = norm3df(s_dji.x, s_dji.y, s_dji.z);
                        SADfloat<9> s_r_jk = norm3df(s_djk.x, s_djk.y, s_djk.z);
                        SADfloat<9> s_r_kl = norm3df(s_dkl.x, s_dkl.y, s_dkl.z);

                        SADvector<9> s_cross_ijk = s_dji ^ s_djk;
                        SADfloat<9> s_sin_ijk =
                            norm3df(s_cross_ijk.x, s_cross_ijk.y,
                                    s_cross_ijk.z) /
                            (s_r_ij * s_r_jk);

                        SADvector<9> s_cross_jkl =
                            ((-1.0f) * s_djk) ^ s_dkl;
                        SADfloat<9> s_sin_jkl =
                            norm3df(s_cross_jkl.x, s_cross_jkl.y,
                                    s_cross_jkl.z) /
                            (s_r_jk * s_r_kl);

                        SADfloat<9> s_unnorm_cos_phi =
                            -(s_dji * s_djk) * (s_djk * s_dkl) +
                            (s_r_jk * s_r_jk) * (s_dji * s_dkl);
                        SADvector<9> s_cross_jk_kl = s_djk ^ s_dkl;
                        SADfloat<9> s_unnorm_sin_phi =
                            -s_r_jk * (s_dji * s_cross_jk_kl);
                        SADfloat<9> s_phi =
                            atan2f(s_unnorm_sin_phi, s_unnorm_cos_phi);
                        SADfloat<9> s_cos_phi = cosf(s_phi);

                        for (int e = 0; e < info.entry_count; e++)
                        {
                            const REAXFF_TORSION_Entry* param =
                                &torsion_entries[info.start_idx + e];

                            SADfloat<6> s_bo_ij(bo_ij_val, 0);
                            SADfloat<6> s_bo_jk(bo_jk_val, 1);
                            SADfloat<6> s_bo_kl(bo_kl_val, 2);
                            SADfloat<6> s_delta_j(delta_j_val, 3);
                            SADfloat<6> s_delta_k(delta_k_val, 4);
                            SADfloat<6> s_bo_jk_pi(bo_jk_pi_val, 5);

                            SADfloat<6> s_boa_ij = s_bo_ij - thb_cut;
                            SADfloat<6> s_boa_jk = s_bo_jk - thb_cut;
                            SADfloat<6> s_boa_kl = s_bo_kl - thb_cut;

                            SADfloat<6> s_fn10 =
                                (1.0f - expf(-p_tor2 * s_boa_ij)) *
                                (1.0f - expf(-p_tor2 * s_boa_jk)) *
                                (1.0f - expf(-p_tor2 * s_boa_kl));
                            SADfloat<6> s_exp_tor3 =
                                expf(-p_tor3 * (s_delta_j + s_delta_k));
                            SADfloat<6> s_exp_tor4 =
                                expf(p_tor4 * (s_delta_j + s_delta_k));
                            SADfloat<6> s_f11_DjDk =
                                (2.0f + s_exp_tor3) /
                                (1.0f + s_exp_tor3 + s_exp_tor4);

                            SADfloat<6> s_tor_diff =
                                2.0f - s_bo_jk_pi - s_f11_DjDk;
                            SADfloat<6> s_exp_tor1 =
                                expf(param->p_tor1 * s_tor_diff * s_tor_diff);
                            SADfloat<6> s_CV =
                                0.5f * (param->V1 * (1.0f + cos_phi) +
                                        param->V2 * s_exp_tor1 *
                                            (1.0f - cosf(2.0f * phi)) +
                                        param->V3 * (1.0f + cosf(3.0f * phi)));

                            SADfloat<6> s_en_tor =
                                s_fn10 * sin_ijk * sin_jkl * s_CV;

                            SADfloat<6> s_en_cot(0.0f);
                            if (param->p_cot1 != 0.0f)
                            {
                                SADfloat<6> s_cot_diff_ij = s_boa_ij - 1.5f;
                                SADfloat<6> s_cot_diff_jk = s_boa_jk - 1.5f;
                                SADfloat<6> s_cot_diff_kl = s_boa_kl - 1.5f;
                                SADfloat<6> s_fn12 =
                                    expf(-p_cot2 * s_cot_diff_ij *
                                         s_cot_diff_ij) *
                                    expf(-p_cot2 * s_cot_diff_jk *
                                         s_cot_diff_jk) *
                                    expf(-p_cot2 * s_cot_diff_kl *
                                         s_cot_diff_kl);
                                s_en_cot = (float)param->p_cot1 * s_fn12 *
                                           (1.0f + (cos_phi * cos_phi - 1.0f) *
                                                       sin_ijk * sin_jkl);
                            }

                            SADfloat<6> s_en_total = s_en_tor + s_en_cot;

                            bool bond_derivatives_finite =
                                ReaxFF_Float_Is_Finite(s_en_tor.val) &&
                                ReaxFF_Float_Is_Finite(s_en_cot.val) &&
                                ReaxFF_Float_Is_Finite(s_en_total.val);
                            for (int derivative_i = 0;
                                 bond_derivatives_finite && derivative_i < 6;
                                 derivative_i++)
                            {
                                bond_derivatives_finite =
                                    ReaxFF_Float_Is_Finite(
                                        s_en_total.dval[derivative_i]);
                            }
                            if (!bond_derivatives_finite)
                            {
                                Record_ReaxFF_Geometry_Error(
                                    geometry_error,
                                    REAXFF_TORSION_NONFINITE, i, j, k, l);
                                continue;
                            }

                            atomicAdd(&d_dE_dBO_s[b_ij], s_en_total.dval[0]);
                            atomicAdd(&d_dE_dBO_pi[b_ij], s_en_total.dval[0]);
                            atomicAdd(&d_dE_dBO_pi2[b_ij], s_en_total.dval[0]);

                            atomicAdd(&d_dE_dBO_s[b_jk], s_en_total.dval[1]);
                            atomicAdd(&d_dE_dBO_pi[b_jk],
                                      s_en_total.dval[1] + s_en_total.dval[5]);
                            atomicAdd(&d_dE_dBO_pi2[b_jk], s_en_total.dval[1]);

                            atomicAdd(&d_dE_dBO_s[b_kl], s_en_total.dval[2]);
                            atomicAdd(&d_dE_dBO_pi[b_kl], s_en_total.dval[2]);
                            atomicAdd(&d_dE_dBO_pi2[b_kl], s_en_total.dval[2]);

                            atomicAdd(&CdDelta[j], s_en_total.dval[3]);
                            atomicAdd(&CdDelta[k], s_en_total.dval[4]);

                            float boa_ij_val = bo_ij_val - thb_cut;
                            float boa_jk_val = bo_jk_val - thb_cut;
                            float boa_kl_val = bo_kl_val - thb_cut;
                            float fn10_val =
                                (1.0f - expf(-p_tor2 * boa_ij_val)) *
                                (1.0f - expf(-p_tor2 * boa_jk_val)) *
                                (1.0f - expf(-p_tor2 * boa_kl_val));
                            float exp_tor3_val =
                                expf(-p_tor3 * (delta_j_val + delta_k_val));
                            float exp_tor4_val =
                                expf(p_tor4 * (delta_j_val + delta_k_val));
                            float f11_DjDk_val =
                                (2.0f + exp_tor3_val) /
                                (1.0f + exp_tor3_val + exp_tor4_val);
                            float tor_diff_val =
                                2.0f - bo_jk_pi_val - f11_DjDk_val;
                            float exp_tor1_val = expf(
                                param->p_tor1 * tor_diff_val * tor_diff_val);

                            SADfloat<9> s_cv_dir =
                                0.5f *
                                (param->V1 * (1.0f + s_cos_phi) +
                                 param->V2 * exp_tor1_val *
                                     (1.0f - cosf(2.0f * s_phi)) +
                                 param->V3 * (1.0f + cosf(3.0f * s_phi)));
                            SADfloat<9> s_en_tor_dir =
                                fn10_val * s_sin_ijk * s_sin_jkl * s_cv_dir;

                            SADfloat<9> s_en_cot_dir(0.0f);
                            if (param->p_cot1 != 0.0f)
                            {
                                float cot_diff_ij_val = boa_ij_val - 1.5f;
                                float cot_diff_jk_val = boa_jk_val - 1.5f;
                                float cot_diff_kl_val = boa_kl_val - 1.5f;
                                float fn12_val =
                                    expf(-p_cot2 * cot_diff_ij_val *
                                         cot_diff_ij_val) *
                                    expf(-p_cot2 * cot_diff_jk_val *
                                         cot_diff_jk_val) *
                                    expf(-p_cot2 * cot_diff_kl_val *
                                         cot_diff_kl_val);
                                s_en_cot_dir =
                                    param->p_cot1 * fn12_val *
                                    (1.0f + (s_cos_phi * s_cos_phi - 1.0f) *
                                                s_sin_ijk * s_sin_jkl);
                            }
                            SADfloat<9> s_en_direct =
                                s_en_tor_dir + s_en_cot_dir;

                            bool direct_derivatives_finite =
                                ReaxFF_Float_Is_Finite(s_en_direct.val);
                            for (int derivative_i = 0;
                                 direct_derivatives_finite && derivative_i < 9;
                                 derivative_i++)
                            {
                                direct_derivatives_finite =
                                    ReaxFF_Float_Is_Finite(
                                        s_en_direct.dval[derivative_i]);
                            }
                            if (!direct_derivatives_finite)
                            {
                                Record_ReaxFF_Geometry_Error(
                                    geometry_error,
                                    REAXFF_TORSION_NONFINITE, i, j, k, l);
                                continue;
                            }

                            VECTOR dE_ddji = {s_en_direct.dval[0],
                                              s_en_direct.dval[1],
                                              s_en_direct.dval[2]};
                            VECTOR dE_ddjk = {s_en_direct.dval[3],
                                              s_en_direct.dval[4],
                                              s_en_direct.dval[5]};
                            VECTOR dE_ddkl = {s_en_direct.dval[6],
                                              s_en_direct.dval[7],
                                              s_en_direct.dval[8]};

                            VECTOR fi = dE_ddji;
                            VECTOR fj = {-(dE_ddji.x + dE_ddjk.x),
                                         -(dE_ddji.y + dE_ddjk.y),
                                         -(dE_ddji.z + dE_ddjk.z)};
                            VECTOR fk = {dE_ddjk.x - dE_ddkl.x,
                                         dE_ddjk.y - dE_ddkl.y,
                                         dE_ddjk.z - dE_ddkl.z};
                            VECTOR fl = dE_ddkl;

                            if (!ReaxFF_Vector_Is_Finite(fi) ||
                                !ReaxFF_Vector_Is_Finite(fj) ||
                                !ReaxFF_Vector_Is_Finite(fk) ||
                                !ReaxFF_Vector_Is_Finite(fl))
                            {
                                Record_ReaxFF_Geometry_Error(
                                    geometry_error,
                                    REAXFF_TORSION_NONFINITE, i, j, k, l);
                                continue;
                            }

                            LTMatrix3 interaction_virial = {0, 0, 0, 0, 0, 0};
                            if (atom_virial)
                            {
                                VECTOR dr_ji = {-dji.x, -dji.y, -dji.z};
                                VECTOR dr_jk = {-djk.x, -djk.y, -djk.z};
                                VECTOR dr_jl = {dr_jk.x - dkl.x,
                                                dr_jk.y - dkl.y,
                                                dr_jk.z - dkl.z};
                                interaction_virial =
                                    Get_Virial_From_Force_Dis(fi, dr_ji) +
                                    Get_Virial_From_Force_Dis(fk, dr_jk) +
                                    Get_Virial_From_Force_Dis(fl, dr_jl);
                                if (!ReaxFF_Matrix_Is_Finite(
                                        interaction_virial))
                                {
                                    Record_ReaxFF_Geometry_Error(
                                        geometry_error,
                                        REAXFF_TORSION_NONFINITE, i, j, k, l);
                                    continue;
                                }
                            }

                            atomicAdd(&frc[i].x, fi.x);
                            atomicAdd(&frc[i].y, fi.y);
                            atomicAdd(&frc[i].z, fi.z);
                            atomicAdd(&frc[j].x, fj.x);
                            atomicAdd(&frc[j].y, fj.y);
                            atomicAdd(&frc[j].z, fj.z);
                            atomicAdd(&frc[k].x, fk.x);
                            atomicAdd(&frc[k].y, fk.y);
                            atomicAdd(&frc[k].z, fk.z);
                            atomicAdd(&frc[l].x, fl.x);
                            atomicAdd(&frc[l].y, fl.y);
                            atomicAdd(&frc[l].z, fl.z);
                            if (atom_virial)
                            {
                                atomicAdd(atom_virial + j,
                                          interaction_virial);
                            }

                            en_tor += s_en_tor.val;
                            en_cot += s_en_cot.val;
                        }
                    }
                }
            }
            atomicAdd(d_energy_tor_sum, (float)en_tor);
            atomicAdd(d_energy_cot_sum, (float)en_cot);
            if (atom_energy)
                atomicAdd(&atom_energy[j], (float)(en_tor + en_cot));
        }
        else
        {
            Record_ReaxFF_Geometry_Error(
                geometry_error, REAXFF_INVALID_ATOM_TYPE, j, -1, type_j);
        }
    }
}

void REAXFF_TORSION::Initial(CONTROLLER* controller, int atom_numbers,
                             const char* module_name)
{
    if (module_name == NULL) module_name = "REAXFF";
    if (!controller->Command_Exist(module_name, "in_file")) return;

    const char* parameter_in_file =
        controller->Original_Command(module_name, "in_file");
    const char* type_in_file =
        controller->Original_Command(module_name, "type_in_file");
    if (parameter_in_file == NULL || type_in_file == NULL)
    {
        controller->printf(
            "REAXFF_TORSION IS NOT INITIALIZED (missing input files)\n\n");
        return;
    }

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
            "REAXFF_TORSION::Initial", "Reason:\n\t%s", reason.c_str());
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
            "REAXFF_TORSION::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }
    if (force_field.general_parameters.size() <= 29)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_TORSION::Initial",
            "Reason:\n\tgeneral parameter count is too small for torsion in "
            "file %s",
            parameter_in_file);
        return;
    }

    float staged_thb_cut = 0.001f;
    if (controller->Command_Exist(module_name, "thb_cutoff"))
    {
        controller->Check_Float(module_name, "thb_cutoff",
                                "REAXFF_TORSION::Initial");
        staged_thb_cut =
            atof(controller->Command(module_name, "thb_cutoff"));
    }
    if (!Float_Memory_Is_Finite(&staged_thb_cut) ||
        staged_thb_cut < 0.0f)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "REAXFF_TORSION::Initial",
            "REAXFF.thb_cutoff must be finite and non-negative, but got %g.",
            staged_thb_cut);
        return;
    }

    const int n_atom_types = static_cast<int>(force_field.atom_types.size());
    int quartet_parameter_count = 0;
    if (!ReaxFF_Checked_Dense_Table_Count(
            n_atom_types, 4, 3 * sizeof(REAXFF_TORSION_Info),
            &quartet_parameter_count))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_TORSION::Initial",
            "Reason:\n\tatom type count %d exceeds the supported four-body "
            "parameter table extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    auto make_entry = [](const REAXFF_TORSION_IR& source)
    {
        REAXFF_TORSION_Entry entry;
        entry.p_tor1 = source.values[3];
        entry.V1 = source.values[0];
        entry.V2 = source.values[1];
        entry.V3 = source.values[2];
        entry.p_tor2 = 0.0f;
        entry.p_cot1 = source.values[4];
        return entry;
    };
    auto quartet_index = [n_atom_types](int type1, int type2, int type3,
                                         int type4)
    {
        return (((type1)*n_atom_types + type2) * n_atom_types + type3) *
                   n_atom_types +
               type4;
    };

    std::map<int, REAXFF_TORSION_Entry> quartet_entries;
    std::set<int> explicit_quartets;
    for (const REAXFF_TORSION_IR& source : force_field.torsions)
    {
        if (source.types[0] == 0) continue;
        const int type1 = source.types[0] - 1;
        const int type2 = source.types[1] - 1;
        const int type3 = source.types[2] - 1;
        const int type4 = source.types[3] - 1;
        const int forward = quartet_index(type1, type2, type3, type4);
        const int reverse = quartet_index(type4, type3, type2, type1);
        const REAXFF_TORSION_Entry entry = make_entry(source);
        quartet_entries[forward] = entry;
        explicit_quartets.insert(forward);
        if (reverse != forward)
        {
            quartet_entries[reverse] = entry;
            explicit_quartets.insert(reverse);
        }
    }

    std::size_t wildcard_count = 0;
    for (const REAXFF_TORSION_IR& source : force_field.torsions)
        if (source.types[0] == 0) wildcard_count++;
    std::size_t wildcard_work = 0;
    std::size_t pair_count = 0;
    if (!ReaxFF_Checked_Size_Multiply(
            static_cast<std::size_t>(n_atom_types),
            static_cast<std::size_t>(n_atom_types), &pair_count) ||
        !ReaxFF_Checked_Size_Multiply(wildcard_count, pair_count,
                                      &wildcard_work) ||
        wildcard_work >
            static_cast<std::size_t>(quartet_parameter_count) * 4)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_TORSION::Initial",
            "Reason:\n\twildcard torsion expansion (%zu entries over %zu "
            "outer-type pairs) exceeds the supported preprocessing work in "
            "file %s",
            wildcard_count, pair_count, parameter_in_file);
        return;
    }

    for (const REAXFF_TORSION_IR& source : force_field.torsions)
    {
        if (source.types[0] != 0) continue;
        const int type2 = source.types[1] - 1;
        const int type3 = source.types[2] - 1;
        const REAXFF_TORSION_Entry entry = make_entry(source);
        for (int outer1 = 0; outer1 < n_atom_types; outer1++)
        {
            for (int outer2 = 0; outer2 < n_atom_types; outer2++)
            {
                const int forward =
                    quartet_index(outer1, type2, type3, outer2);
                if (explicit_quartets.count(forward) == 0)
                    quartet_entries[forward] = entry;
                const int reverse =
                    quartet_index(outer2, type3, type2, outer1);
                if (explicit_quartets.count(reverse) == 0)
                    quartet_entries[reverse] = entry;
            }
        }
    }

    std::vector<REAXFF_TORSION_Info> torsion_info(
        quartet_parameter_count, REAXFF_TORSION_Info{0, 0});
    std::vector<REAXFF_TORSION_Entry> sorted_entries;
    for (const auto& indexed_entry : quartet_entries)
    {
        if (indexed_entry.first < 0 ||
            indexed_entry.first >= quartet_parameter_count ||
            sorted_entries.size() >= static_cast<std::size_t>(INT_MAX))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "REAXFF_TORSION::Initial",
                "Reason:\n\ttorsion parameter index/count exceeds the "
                "supported representation in file %s",
                parameter_in_file);
            return;
        }
        REAXFF_TORSION_Info& info = torsion_info[indexed_entry.first];
        info.start_idx = static_cast<int>(sorted_entries.size());
        info.entry_count = 1;
        sorted_entries.push_back(indexed_entry.second);
    }

    REAXFF_TORSION_Info* staged_torsion_info = NULL;
    REAXFF_TORSION_Entry* staged_torsion_entries = NULL;
    int* staged_atom_type = NULL;
    Malloc_Safely((void**)&staged_torsion_info,
                  sizeof(REAXFF_TORSION_Info) * torsion_info.size());
    if (!sorted_entries.empty())
        Malloc_Safely((void**)&staged_torsion_entries,
                      sizeof(REAXFF_TORSION_Entry) * sorted_entries.size());
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    memcpy(staged_torsion_info, torsion_info.data(),
           sizeof(REAXFF_TORSION_Info) * torsion_info.size());
    if (!sorted_entries.empty())
        memcpy(staged_torsion_entries, sorted_entries.data(),
               sizeof(REAXFF_TORSION_Entry) * sorted_entries.size());
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);

    this->controller = controller;
    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    p_tor2 = force_field.general_parameters[23];
    p_tor3 = force_field.general_parameters[24];
    p_tor4 = force_field.general_parameters[25];
    p_cot2 = force_field.general_parameters[27];
    thb_cut = staged_thb_cut;
    h_torsion_info = staged_torsion_info;
    h_torsion_entries = staged_torsion_entries;
    h_atom_type = staged_atom_type;

    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely(
        (void**)&d_torsion_info, h_torsion_info,
        sizeof(REAXFF_TORSION_Info) * quartet_parameter_count);
    if (!sorted_entries.empty())
        Device_Malloc_And_Copy_Safely(
            (void**)&d_torsion_entries, h_torsion_entries,
            sizeof(REAXFF_TORSION_Entry) * sorted_entries.size());
    Device_Malloc_Safely((void**)&d_energy_tor_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_energy_cot_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_geometry_error,
                         sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);

    controller->Step_Print_Initial("REAXFF_TOR", "%14.7e");
    controller->Step_Print_Initial("REAXFF_CONJ", "%14.7e");
    is_initialized = 1;
}
void REAXFF_TORSION::Calculate_Torsion_Energy_And_Force(
    int atom_numbers, const VECTOR* crd, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, REAXFF_BOND_ORDER* bo_module,
    const float* Delta_boc, const int need_atom_energy, float* atom_energy,
    const int need_virial, LTMatrix3* atom_virial)
{
    if (!is_initialized) return;
    dim3 blockSize(32);
    dim3 gridSize((atom_numbers + blockSize.x - 1) / blockSize.x);
    deviceMemset(d_energy_tor_sum, 0, sizeof(float));
    deviceMemset(d_energy_cot_sum, 0, sizeof(float));
    deviceMemset(d_geometry_error, 0,
                 sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);

    Launch_Device_Kernel(
        Calculate_Torsion_Kernel, gridSize, blockSize, 0, NULL, atom_numbers,
        crd, d_atom_type, p_tor2, p_tor3, p_tor4, p_cot2, thb_cut, Delta_boc,
        d_torsion_info, d_torsion_entries, atom_type_numbers,
        bo_module->d_corrected_bo_s, bo_module->d_corrected_bo_pi,
        bo_module->d_corrected_bo_pi2, d_dE_dBO_s, d_dE_dBO_pi, d_dE_dBO_pi2,
        d_CdDelta, cell, rcell, atom_energy, frc,
        need_virial ? atom_virial : NULL, d_energy_tor_sum, d_energy_cot_sum,
        bo_module->d_bond_count, bo_module->d_bond_offset,
        bo_module->d_bond_nbr, bo_module->d_bond_idx, d_geometry_error);
    Check_Geometry_Error();
}

void REAXFF_TORSION::Check_Geometry_Error()
{
    int error[REAXFF_GEOMETRY_ERROR_SIZE] = {0, -1, -1, -1, -1};
    deviceMemcpy(error, d_geometry_error, sizeof(error),
                 deviceMemcpyDeviceToHost);
    if (error[0] == REAXFF_GEOMETRY_OK) return;
    if (error[0] == REAXFF_INVALID_ATOM_TYPE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF_TORSION::Calculate_Torsion_Energy_And_Force",
            "Invalid ReaxFF atom type in torsion evaluation record "
            "[%d, %d, %d, %d] (valid range is [0, %d)).",
            error[1], error[2], error[3], error[4], atom_type_numbers);
        return;
    }
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown,
        "REAXFF_TORSION::Calculate_Torsion_Energy_And_Force",
        "Reason:\n\tlocal atoms %d %d %d %d have an undefined zero-bond/"
        "collinear ReaxFF torsion geometry or produce a non-finite/"
        "unrepresentable torsion energy, derivative, force, or virial\n",
        error[1], error[2], error[3], error[4]);
}

void REAXFF_TORSION::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_tor, d_energy_tor_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&h_energy_cot, d_energy_cot_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print("REAXFF_TOR", h_energy_tor, true);
    controller->Step_Print("REAXFF_CONJ", h_energy_cot, true);
}
