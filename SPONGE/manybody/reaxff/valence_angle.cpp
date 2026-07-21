#include "valence_angle.h"

#include "bond_order.h"  // for find_bond_index
#include "reaxff_geometry.h"
#include "reaxff_input.h"

static __global__ void Calculate_Valence_Angle_Kernel(
    int atom_numbers, const VECTOR* crd, const int* atom_type,
    const float* Delta_boc, const float* Delta, const float* Delta_val,
    const float* p_val3, const float* p_val5,
    const REAXFF_VALENCE_ANGLE_PARAMS p, const REAXFF_THBP_Info* thbp_info,
    const REAXFF_THBP_Entry* thbp_entries, int atom_type_numbers,
    const float* bo_s, const float* bo_pi, const float* bo_pi2,
    const float* total_bo, const float* nlp, const float* vlpex,
    const float* dDelta_lp, const LTMatrix3 cell, const LTMatrix3 rcell,
    float* d_dE_dBO_s, float* d_dE_dBO_pi, float* d_dE_dBO_pi2, float* CdDelta,
    float* atom_energy, VECTOR* frc, LTMatrix3* atom_virial,
    float* d_energy_ang_sum, float* d_energy_pen_sum, float* d_energy_coa_sum,
    const int* bond_count, const int* bond_offset, const int* bond_nbr,
    const int* bond_idx_arr, int* geometry_error)
{
    SIMPLE_DEVICE_FOR(j, atom_numbers)
    {
        int type_j = atom_type[j];
        if (type_j >= 0 && type_j < atom_type_numbers)
        {
            VECTOR rj = crd[j];

            float p_val3_j = p_val3[type_j];
            float p_val5_j = p_val5[type_j];
            float delta_boc_j_val = Delta_boc[j];
            float delta_j_val = Delta[j];
            float delta_val_j_val = Delta_val[j];

            float SBOp = 0, prod_SBO = 1.0f;
            int bc_j = bond_count[j];
            int bo_j = bond_offset[j];
            for (int t = 0; t < bc_j; t++)
            {
                int b_jt = bond_idx_arr[bo_j + t];
                float bo_jt = bo_s[b_jt] + bo_pi[b_jt] + bo_pi2[b_jt];

                SBOp += (bo_pi[b_jt] + bo_pi2[b_jt]);
                float bo_jt_sq = bo_jt * bo_jt;
                float bo_jt_p4 = bo_jt_sq * bo_jt_sq;
                prod_SBO *= expf(-bo_jt_p4 * bo_jt_p4);
            }

            bool has_lp_corr = (vlpex[j] < 0.0f);
            float vlpadj = has_lp_corr ? nlp[j] : 0.0f;
            float SBO = SBOp + (1.0f - prod_SBO) *
                                   (-delta_boc_j_val - p.p_val8 * vlpadj);
            float dSBO1 =
                -8.0f * prod_SBO * (delta_boc_j_val + p.p_val8 * vlpadj);
            float dSBO2 = has_lp_corr ? (prod_SBO - 1.0f) *
                                            (1.0f - p.p_val8 * dDelta_lp[j])
                                      : (prod_SBO - 1.0f);
            float SBO2;
            float CSBO2;
            if (SBO <= 0.0f)
            {
                SBO2 = 0.0f;
                CSBO2 = 0.0f;
            }
            else if (SBO <= 1.0f)
            {
                SBO2 = powf(SBO, p.p_val9);
                CSBO2 = p.p_val9 * powf(SBO, p.p_val9 - 1.0f);
            }
            else if (SBO < 2.0f)
            {
                SBO2 = 2.0f - powf(2.0f - SBO, p.p_val9);
                CSBO2 = p.p_val9 * powf(2.0f - SBO, p.p_val9 - 1.0f);
            }
            else
            {
                SBO2 = 2.0f;
                CSBO2 = 0.0f;
            }

            for (int bi = 0; bi < bc_j; bi++)
            {
                int b_ij = bond_idx_arr[bo_j + bi];
                int i = bond_nbr[bo_j + bi];
                int type_i = atom_type[i];
                if (type_i < 0 || type_i >= atom_type_numbers)
                {
                    Record_ReaxFF_Geometry_Error(
                        geometry_error, REAXFF_INVALID_ATOM_TYPE, i, j,
                        type_i, type_j);
                    continue;
                }
                float bo_ij_val = bo_s[b_ij] + bo_pi[b_ij] + bo_pi2[b_ij];
                float boa_ij_val = bo_ij_val - p.thb_cut;

                if (boa_ij_val <= 0) continue;

                VECTOR ri = crd[i];
                VECTOR dji = Get_Periodic_Displacement(rj, ri, cell, rcell);

                for (int bk = bi + 1; bk < bc_j; bk++)
                {
                    int b_kj = bond_idx_arr[bo_j + bk];
                    int k = bond_nbr[bo_j + bk];
                    int type_k = atom_type[k];
                    if (type_k < 0 || type_k >= atom_type_numbers)
                    {
                        Record_ReaxFF_Geometry_Error(
                            geometry_error, REAXFF_INVALID_ATOM_TYPE, i, j, k,
                            type_k);
                        continue;
                    }
                    float bo_jk_val = bo_s[b_kj] + bo_pi[b_kj] + bo_pi2[b_kj];
                    float boa_jk_val = bo_jk_val - p.thb_cut;

                    if (boa_jk_val <= 0) continue;
                    if (bo_ij_val * bo_jk_val <= p.thb_cutsq) continue;

                    VECTOR rk = crd[k];
                    VECTOR djk = Get_Periodic_Displacement(rj, rk, cell, rcell);

                    REAXFF_ANGLE_GEOMETRY geometry;
                    if (!Compute_ReaxFF_Angle_Geometry(dji, djk, &geometry))
                    {
                        Record_ReaxFF_Geometry_Error(
                            geometry_error, REAXFF_ANGLE_UNDEFINED, i, j, k);
                        continue;
                    }
                    const float theta = geometry.theta;

                    int tri_info_idx = (type_i * atom_type_numbers + type_j) *
                                           atom_type_numbers +
                                       type_k;
                    REAXFF_THBP_Info info = thbp_info[tri_info_idx];

                    for (int e = 0; e < info.entry_count; e++)
                    {
                        const REAXFF_THBP_Entry* param =
                            &thbp_entries[info.start_idx + e];
                        // Exactly zero is the force-field sentinel for an
                        // inactive three-body entry.  A small nonzero
                        // coefficient remains a real interaction.
                        if (param->p_val1 == 0.0f) continue;

                        SADfloat<3> s_bo_ij(bo_ij_val, 0);
                        SADfloat<3> s_bo_jk(bo_jk_val, 1);
                        SADfloat<3> s_delta_j(delta_boc_j_val, 2);
                        SADfloat<3> s_delta_pen =
                            s_delta_j + (delta_j_val - delta_boc_j_val);

                        SADfloat<3> s_boa_ij = s_bo_ij - p.thb_cut;
                        SADfloat<3> s_boa_jk = s_bo_jk - p.thb_cut;

                        SADfloat<3> exp3ij = expf(
                            -p_val3_j * powf(s_boa_ij, (float)param->p_val4));
                        SADfloat<3> f7_ij = 1.0f - exp3ij;
                        SADfloat<3> exp3jk = expf(
                            -p_val3_j * powf(s_boa_jk, (float)param->p_val4));
                        SADfloat<3> f7_jk = 1.0f - exp3jk;

                        SADfloat<3> expval6 = expf(p.p_val6 * s_delta_j);
                        SADfloat<3> expval7 = expf(-param->p_val7 * s_delta_j);
                        SADfloat<3> trm8 = 1.0f + expval6 + expval7;
                        SADfloat<3> f8_Dj =
                            p_val5_j -
                            ((p_val5_j - 1.0f) * (2.0f + expval6) / trm8);

                        float theta_0 =
                            (180.0f -
                             param->theta_00 *
                                 (1.0f - expf(-p.p_val10 * (2.0f - SBO2)))) *
                            (CONSTANT_Pi / 180.0f);
                        float expval2theta =
                            expf(-param->p_val2 * (theta_0 - theta) *
                                 (theta_0 - theta));
                        float val12theta =
                            (param->p_val1 >= 0)
                                ? param->p_val1 * (1.0f - expval2theta)
                                : param->p_val1 * (-expval2theta);

                        SADfloat<3> s_en_ang =
                            f7_ij * f7_jk * f8_Dj * val12theta;

                        SADfloat<3> s_pen_diff_ij = s_boa_ij - 2.0f;
                        SADfloat<3> s_pen_diff_jk = s_boa_jk - 2.0f;
                        SADfloat<3> exp_pen2ij =
                            expf(-p.p_pen2 * s_pen_diff_ij * s_pen_diff_ij);
                        SADfloat<3> exp_pen2jk =
                            expf(-p.p_pen2 * s_pen_diff_jk * s_pen_diff_jk);

                        SADfloat<3> exp_pen3 = expf(-p.p_pen3 * s_delta_pen);
                        SADfloat<3> exp_pen4 = expf(p.p_pen4 * s_delta_pen);
                        SADfloat<3> trm_pen34 = 1.0f + exp_pen3 + exp_pen4;
                        SADfloat<3> f9_Dj = (2.0f + exp_pen3) / trm_pen34;
                        SADfloat<3> s_en_pen =
                            param->p_pen1 * f9_Dj * exp_pen2ij * exp_pen2jk;

                        SADfloat<3> exp_coa2 =
                            expf(p.p_coa2 * (s_delta_j + delta_val_j_val -
                                             delta_boc_j_val));
                        SADfloat<3> s_coa_diff_i = total_bo[i] - bo_ij_val;
                        SADfloat<3> s_coa_diff_k = total_bo[k] - bo_jk_val;
                        SADfloat<3> s_coa_diff_ij = s_boa_ij - 1.5f + p.thb_cut;
                        SADfloat<3> s_coa_diff_jk = s_boa_jk - 1.5f + p.thb_cut;

                        SADfloat<3> s_en_coa =
                            param->p_coa1 / (1.0f + exp_coa2) *
                            expf(-p.p_coa3 * s_coa_diff_i * s_coa_diff_i) *
                            expf(-p.p_coa3 * s_coa_diff_k * s_coa_diff_k) *
                            expf(-p.p_coa4 * s_coa_diff_ij * s_coa_diff_ij) *
                            expf(-p.p_coa4 * s_coa_diff_jk * s_coa_diff_jk);

                        SADfloat<3> s_en_total = s_en_ang + s_en_pen + s_en_coa;

                        atomicAdd(&d_dE_dBO_s[b_ij], s_en_total.dval[0]);
                        atomicAdd(&d_dE_dBO_pi[b_ij], s_en_total.dval[0]);
                        atomicAdd(&d_dE_dBO_pi2[b_ij], s_en_total.dval[0]);

                        atomicAdd(&d_dE_dBO_s[b_kj], s_en_total.dval[1]);
                        atomicAdd(&d_dE_dBO_pi[b_kj], s_en_total.dval[1]);
                        atomicAdd(&d_dE_dBO_pi2[b_kj], s_en_total.dval[1]);

                        atomicAdd(&CdDelta[j], s_en_total.dval[2]);

                        float dE_dtheta = -f7_ij.val * f7_jk.val * f8_Dj.val *
                                          param->p_val1 * expval2theta * 2.0f *
                                          param->p_val2 * (theta_0 - theta);
                        if (!ReaxFF_Float_Is_Finite(dE_dtheta) ||
                            !ReaxFF_Float_Is_Finite(s_en_total.val) ||
                            !ReaxFF_Float_Is_Finite(s_en_ang.val) ||
                            !ReaxFF_Float_Is_Finite(s_en_pen.val) ||
                            !ReaxFF_Float_Is_Finite(s_en_coa.val))
                        {
                            Record_ReaxFF_Geometry_Error(
                                geometry_error, REAXFF_ANGLE_NONFINITE, i, j,
                                k);
                            continue;
                        }
                        float Ctheta_0 = p.p_val10 * (CONSTANT_Pi / 180.0f) *
                                         param->theta_00 *
                                         expf(-p.p_val10 * (2.0f - SBO2));
                        float CEval5 = -dE_dtheta * Ctheta_0 * CSBO2;
                        float CEval6 = CEval5 * dSBO1;
                        float CEval7 = CEval5 * dSBO2;
                        atomicAdd(&CdDelta[j], CEval7);

                        for (int pt = 0; pt < bc_j; pt++)
                        {
                            int b_jt = bond_idx_arr[bo_j + pt];
                            float bo_jt =
                                bo_s[b_jt] + bo_pi[b_jt] + bo_pi2[b_jt];
                            if (bo_jt <= 0.0f) continue;
                            float bo_jt_2 = bo_jt * bo_jt;
                            float bo_jt_4 = bo_jt_2 * bo_jt_2;
                            float bo_jt_7 = bo_jt_4 * bo_jt_2 * bo_jt;
                            float dE_dbo_total = CEval6 * bo_jt_7;

                            atomicAdd(&d_dE_dBO_s[b_jt], dE_dbo_total);
                            atomicAdd(&d_dE_dBO_pi[b_jt],
                                      dE_dbo_total + CEval5);
                            atomicAdd(&d_dE_dBO_pi2[b_jt],
                                      dE_dbo_total + CEval5);
                        }

                        VECTOR fi = {0.0f, 0.0f, 0.0f};
                        VECTOR fk = {0.0f, 0.0f, 0.0f};
                        if (geometry.is_collinear)
                        {
                            // At an exact endpoint the Cartesian angle
                            // derivative is unique only when the angular
                            // derivative itself is exactly zero.
                            if (dE_dtheta != 0.0f)
                            {
                                Record_ReaxFF_Geometry_Error(
                                    geometry_error, REAXFF_ANGLE_UNDEFINED, i,
                                    j, k);
                                continue;
                            }
                        }
                        else
                        {
                            fi = {static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_du[0]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_du[1]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_du[2])};
                            fk = {static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_dv[0]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_dv[1]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_dv[2])};
                        }
                        VECTOR fj = {-(fi.x + fk.x), -(fi.y + fk.y),
                                     -(fi.z + fk.z)};
                        if (!ReaxFF_Vector_Is_Finite(fi) ||
                            !ReaxFF_Vector_Is_Finite(fj) ||
                            !ReaxFF_Vector_Is_Finite(fk))
                        {
                            Record_ReaxFF_Geometry_Error(
                                geometry_error, REAXFF_ANGLE_NONFINITE, i, j,
                                k);
                            continue;
                        }

                        LTMatrix3 interaction_virial = {0, 0, 0, 0, 0, 0};
                        if (atom_virial)
                        {
                            interaction_virial =
                                interaction_virial -
                                Get_Virial_From_Force_Dis(fi, dji) -
                                Get_Virial_From_Force_Dis(fk, djk);
                            if (!ReaxFF_Matrix_Is_Finite(interaction_virial))
                            {
                                Record_ReaxFF_Geometry_Error(
                                    geometry_error, REAXFF_ANGLE_NONFINITE, i,
                                    j, k);
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
                        if (atom_virial)
                        {
                            atomicAdd(atom_virial + j, interaction_virial);
                        }

                        if (atom_energy)
                        {
                            float en_total = s_en_total.val;
                            atomicAdd(&atom_energy[j], en_total);
                            atomicAdd(d_energy_ang_sum, s_en_ang.val);
                            atomicAdd(d_energy_pen_sum, s_en_pen.val);
                            atomicAdd(d_energy_coa_sum, s_en_coa.val);
                        }
                    }
                }
            }
        }
        else
        {
            Record_ReaxFF_Geometry_Error(
                geometry_error, REAXFF_INVALID_ATOM_TYPE, j, -1, type_j);
        }
    }
}

void REAXFF_VALENCE_ANGLE::Initial(CONTROLLER* controller, int atom_numbers,
                                   const char* module_name)
{
    if (module_name == NULL) module_name = "REAXFF";
    if (!controller->Command_Exist(module_name, "in_file")) return;

    controller->printf("START INITIALIZING REAXFF VALENCE ANGLE\n");
    const char* parameter_in_file =
        controller->Original_Command(module_name, "in_file");
    const char* type_in_file =
        controller->Original_Command(module_name, "type_in_file");
    if (parameter_in_file == NULL || type_in_file == NULL)
    {
        controller->printf(
            "REAXFF_VALENCE_ANGLE IS NOT INITIALIZED (missing input "
            "files)\n\n");
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
            "REAXFF_VALENCE_ANGLE::Initial", "Reason:\n\t%s",
            reason.c_str());
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
            "REAXFF_VALENCE_ANGLE::Initial", "Reason:\n\t%s",
            reason.c_str());
        return;
    }
    if (force_field.general_parameters.size() <= 38)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_VALENCE_ANGLE::Initial",
            "Reason:\n\tgeneral parameter count is too small for valence "
            "angle in file %s",
            parameter_in_file);
        return;
    }

    REAXFF_VALENCE_ANGLE_PARAMS staged_params = {};
    staged_params.p_coa2 = force_field.general_parameters[2];
    staged_params.p_val6 = force_field.general_parameters[14];
    staged_params.p_val9 = force_field.general_parameters[16];
    staged_params.p_val10 = force_field.general_parameters[17];
    staged_params.p_pen2 = force_field.general_parameters[19];
    staged_params.p_pen3 = force_field.general_parameters[20];
    staged_params.p_pen4 = force_field.general_parameters[21];
    staged_params.p_coa4 = force_field.general_parameters[30];
    staged_params.p_val8 = force_field.general_parameters[33];
    staged_params.p_coa3 = force_field.general_parameters[38];
    staged_params.thb_cut = 0.001f;
    if (controller->Command_Exist(module_name, "thb_cutoff"))
    {
        controller->Check_Float(
            module_name, "thb_cutoff", "REAXFF_VALENCE_ANGLE::Initial");
        staged_params.thb_cut =
            atof(controller->Command(module_name, "thb_cutoff"));
    }
    staged_params.thb_cutsq =
        staged_params.thb_cut * staged_params.thb_cut;
    if (!Float_Memory_Is_Finite(&staged_params.thb_cut) ||
        !Float_Memory_Is_Finite(&staged_params.thb_cutsq) ||
        staged_params.thb_cut < 0.0f)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "REAXFF_VALENCE_ANGLE::Initial",
            "REAXFF.thb_cutoff must be finite and non-negative, but got %g.",
            staged_params.thb_cut);
        return;
    }

    const int n_atom_types = static_cast<int>(force_field.atom_types.size());
    int triplet_parameter_count = 0;
    if (!ReaxFF_Checked_Dense_Table_Count(
            n_atom_types, 3, 3 * sizeof(REAXFF_THBP_Info),
            &triplet_parameter_count))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_VALENCE_ANGLE::Initial",
            "Reason:\n\tatom type count %d exceeds the supported three-body "
            "parameter table extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    std::vector<float> p_val3_values(n_atom_types);
    std::vector<float> p_val5_values(n_atom_types);
    std::vector<float> mass_values(n_atom_types);
    std::vector<float> valency_boc_values(n_atom_types);
    for (int i = 0; i < n_atom_types; i++)
    {
        const REAXFF_ATOM_TYPE_IR& atom = force_field.atom_types[i];
        mass_values[i] = atom.values[0][2];
        valency_boc_values[i] = atom.values[1][2];
        p_val3_values[i] = atom.values[3][1];
        p_val5_values[i] = atom.values[3][4];
    }

    std::map<int, std::vector<REAXFF_THBP_Entry>> triplet_entries;
    for (const REAXFF_ANGLE_IR& source : force_field.angles)
    {
        REAXFF_THBP_Entry entry;
        entry.theta_00 = source.values[0];
        entry.p_val1 = source.values[1];
        entry.p_val2 = source.values[2];
        entry.p_coa1 = source.values[3];
        entry.p_val7 = source.values[4];
        entry.p_pen1 = source.values[5];
        entry.p_val4 = source.values[6];

        const int idx1 = source.types[0] - 1;
        const int idx2 = source.types[1] - 1;
        const int idx3 = source.types[2] - 1;
        const int forward =
            (idx1 * n_atom_types + idx2) * n_atom_types + idx3;
        triplet_entries[forward].push_back(entry);
        if (idx1 != idx3)
        {
            const int reverse =
                (idx3 * n_atom_types + idx2) * n_atom_types + idx1;
            triplet_entries[reverse].push_back(entry);
        }
    }

    std::vector<REAXFF_THBP_Info> thbp_info(
        triplet_parameter_count, REAXFF_THBP_Info{0, 0});
    std::vector<REAXFF_THBP_Entry> sorted_entries;
    for (const auto& indexed_entries : triplet_entries)
    {
        if (indexed_entries.first < 0 ||
            indexed_entries.first >= triplet_parameter_count ||
            sorted_entries.size() >
                static_cast<std::size_t>(INT_MAX) -
                    indexed_entries.second.size())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "REAXFF_VALENCE_ANGLE::Initial",
                "Reason:\n\tthree-body parameter index/count exceeds the "
                "supported representation in file %s",
                parameter_in_file);
            return;
        }
        REAXFF_THBP_Info& info = thbp_info[indexed_entries.first];
        info.start_idx = static_cast<int>(sorted_entries.size());
        info.entry_count =
            static_cast<int>(indexed_entries.second.size());
        sorted_entries.insert(sorted_entries.end(),
                              indexed_entries.second.begin(),
                              indexed_entries.second.end());
    }

    auto allocate_float_copy = [](const std::vector<float>& values)
    {
        float* result = NULL;
        Malloc_Safely((void**)&result, sizeof(float) * values.size());
        memcpy(result, values.data(), sizeof(float) * values.size());
        return result;
    };
    float* staged_p_val3 = allocate_float_copy(p_val3_values);
    float* staged_p_val5 = allocate_float_copy(p_val5_values);
    float* staged_mass = allocate_float_copy(mass_values);
    float* staged_valency_boc = allocate_float_copy(valency_boc_values);
    REAXFF_THBP_Info* staged_thbp_info = NULL;
    REAXFF_THBP_Entry* staged_thbp_entries = NULL;
    int* staged_atom_type = NULL;
    Malloc_Safely((void**)&staged_thbp_info,
                  sizeof(REAXFF_THBP_Info) * thbp_info.size());
    if (!sorted_entries.empty())
        Malloc_Safely((void**)&staged_thbp_entries,
                      sizeof(REAXFF_THBP_Entry) * sorted_entries.size());
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    memcpy(staged_thbp_info, thbp_info.data(),
           sizeof(REAXFF_THBP_Info) * thbp_info.size());
    if (!sorted_entries.empty())
        memcpy(staged_thbp_entries, sorted_entries.data(),
               sizeof(REAXFF_THBP_Entry) * sorted_entries.size());
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);

    this->controller = controller;
    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    params = staged_params;
    h_p_val3 = staged_p_val3;
    h_p_val5 = staged_p_val5;
    h_mass = staged_mass;
    h_valency_boc = staged_valency_boc;
    h_thbp_info = staged_thbp_info;
    h_thbp_entries = staged_thbp_entries;
    h_atom_type = staged_atom_type;

    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_p_val3, h_p_val3,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_p_val5, h_p_val5,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_mass, h_mass,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_valency_boc, h_valency_boc,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely(
        (void**)&d_thbp_info, h_thbp_info,
        sizeof(REAXFF_THBP_Info) * triplet_parameter_count);
    if (!sorted_entries.empty())
        Device_Malloc_And_Copy_Safely(
            (void**)&d_thbp_entries, h_thbp_entries,
            sizeof(REAXFF_THBP_Entry) * sorted_entries.size());

    Device_Malloc_Safely((void**)&d_energy_ang_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_energy_pen_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_energy_coa_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_geometry_error,
                         sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);

    controller->Step_Print_Initial("REAXFF_ANG", "%14.7e");
    controller->Step_Print_Initial("REAXFF_PEN", "%14.7e");
    controller->Step_Print_Initial("REAXFF_COA", "%14.7e");
    is_initialized = 1;
    controller->printf("END INITIALIZING REAXFF VALENCE ANGLE\n\n");
}
void REAXFF_VALENCE_ANGLE::Calculate_Valence_Angle_Energy_And_Force(
    int atom_numbers, const VECTOR* crd, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, REAXFF_BOND_ORDER* bo_module,
    const float* Delta, const float* Delta_boc, const float* Delta_val,
    const float* nlp, const float* vlpex, const float* dDelta_lp,
    float* CdDelta, const int need_atom_energy, float* atom_energy,
    const int need_virial, LTMatrix3* atom_virial)
{
    if (!is_initialized) return;

    dim3 blockSize(32);
    dim3 gridSize((atom_numbers + blockSize.x - 1) / blockSize.x);

    deviceMemset(d_energy_ang_sum, 0, sizeof(float));
    deviceMemset(d_energy_pen_sum, 0, sizeof(float));
    deviceMemset(d_energy_coa_sum, 0, sizeof(float));
    deviceMemset(d_geometry_error, 0,
                 sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);

    Launch_Device_Kernel(
        Calculate_Valence_Angle_Kernel, gridSize, blockSize, 0, NULL,
        atom_numbers, crd, d_atom_type, Delta_boc, Delta, Delta_val, d_p_val3,
        d_p_val5, params, d_thbp_info, d_thbp_entries, atom_type_numbers,
        bo_module->d_corrected_bo_s, bo_module->d_corrected_bo_pi,
        bo_module->d_corrected_bo_pi2, bo_module->d_total_corrected_bond_order,
        nlp, vlpex, dDelta_lp, cell, rcell, d_dE_dBO_s, d_dE_dBO_pi,
        d_dE_dBO_pi2, CdDelta, need_atom_energy ? atom_energy : NULL, frc,
        need_virial ? atom_virial : NULL, d_energy_ang_sum, d_energy_pen_sum,
        d_energy_coa_sum, bo_module->d_bond_count, bo_module->d_bond_offset,
        bo_module->d_bond_nbr, bo_module->d_bond_idx, d_geometry_error);

    Check_Geometry_Error();

#ifdef USE_GPU
    deviceError_t err = deviceGetLastError();
    if (err != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF_VALENCE_ANGLE::Calculate_Valence_Angle_Energy_And_Force",
            "Device failure in ReaxFF valence-angle kernel: %s",
            deviceGetErrorString(err));
    }
#endif
}

void REAXFF_VALENCE_ANGLE::Check_Geometry_Error()
{
    int error[REAXFF_GEOMETRY_ERROR_SIZE] = {0, -1, -1, -1, -1};
    deviceMemcpy(error, d_geometry_error, sizeof(error),
                 deviceMemcpyDeviceToHost);
    if (error[0] == REAXFF_GEOMETRY_OK) return;
    if (error[0] == REAXFF_INVALID_ATOM_TYPE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF_VALENCE_ANGLE::Calculate_Valence_Angle_Energy_And_Force",
            "Invalid ReaxFF atom type in valence-angle evaluation record "
            "[%d, %d, %d, %d] (valid range is [0, %d)).",
            error[1], error[2], error[3], error[4], atom_type_numbers);
        return;
    }
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown,
        "REAXFF_VALENCE_ANGLE::Calculate_Valence_Angle_Energy_And_Force",
        "Reason:\n\tlocal atoms %d %d %d have an undefined zero-arm/"
        "collinear ReaxFF valence-angle geometry or produce a non-finite/"
        "unrepresentable angle energy, derivative, force, or virial\n",
        error[1], error[2], error[3]);
}

void REAXFF_VALENCE_ANGLE::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_ang, d_energy_ang_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&h_energy_pen, d_energy_pen_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&h_energy_coa, d_energy_coa_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print("REAXFF_ANG", h_energy_ang, true);
    controller->Step_Print("REAXFF_PEN", h_energy_pen, true);
    controller->Step_Print("REAXFF_COA", h_energy_coa, true);
}
