#include "over_under.h"

#include "reaxff_input.h"

static __global__ void Calculate_Delta_Kernel(
    int atom_numbers, const int* atom_type, const float* total_corrected_bo,
    const float* valency, const float* valency_e, const float* valency_boc,
    const float* valency_val, const float* mass, float p_lp1, float* d_Delta,
    float* d_Delta_boc, float* d_Delta_val, float* d_Delta_lp, float* d_nlp,
    float* d_vlpex, float* d_Delta_lp_temp, float* d_dDelta_lp)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type = atom_type[i];
        float total_bo = total_corrected_bo[i];

        float val = valency[type];
        float val_e = valency_e[type];
        float val_boc = valency_boc[type];
        float val_val = valency_val[type];
        float m = mass[type];

        d_Delta[i] = total_bo - val;
        d_Delta_boc[i] = total_bo - val_boc;
        d_Delta_val[i] = total_bo - val_val;

        float Delta_e = total_bo - val_e;

        float vlpex = Delta_e - 2.0f * (int)(Delta_e / 2.0f);
        d_vlpex[i] = vlpex;

        float explp1 = expf(-p_lp1 * (2.0f + vlpex) * (2.0f + vlpex));
        float nlp = explp1 - (int)(Delta_e / 2.0f);
        d_nlp[i] = nlp;

        float nlp_opt = 0.5f * (val_e - val);
        d_Delta_lp[i] = nlp_opt - nlp;

        float Clp = 2.0f * p_lp1 * explp1 * (2.0f + vlpex);
        d_dDelta_lp[i] = Clp;

        if (m > 21.0f)
        {
            d_Delta_lp_temp[i] = 0.0f;
        }
        else
        {
            d_Delta_lp_temp[i] = d_Delta_lp[i];
        }
    }
}

static __global__ void Calculate_Energy_Force_Prep_Kernel(
    int atom_numbers, const int* atom_type, const float* mass,
    const float* Delta, const float* Delta_lp, const float* Delta_lp_temp,
    const float* dDelta_lp, const float* bo_s, const float* bo_pi,
    const float* bo_pi2, const float* p_ovun1, const float* De_s,
    int atom_type_numbers, const float* p_lp2, const float* valency,
    float p_ovun3, float p_ovun4, float p_ovun6, float p_ovun7, float p_ovun8,
    const float* p_ovun2, const float* p_ovun5, float* d_dE_dBO_s,
    float* d_dE_dBO_pi, float* d_dE_dBO_pi2, float* CdDelta, float* atom_energy,
    float* d_energy_ovun_sum, float* d_energy_elp_sum, const int* bond_count,
    const int* bond_offset, const int* bond_nbr, const int* bond_idx_arr)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type_i = atom_type[i];
        float m = mass[type_i];
        float dfvl = (m > 21.0f) ? 0.0f : 1.0f;
        float val = valency[type_i];

        float s_ovun1 = 0.0f;
        float s_ovun2 = 0.0f;

        int start = bond_offset[i];
        int end = start + bond_count[i];
        for (int kk = start; kk < end; kk++)
        {
            int j = bond_nbr[kk];
            int b = bond_idx_arr[kk];
            float b_s = bo_s[b];
            float b_pi = bo_pi[b];
            float b_pi2 = bo_pi2[b];
            float bo_total = b_s + b_pi + b_pi2;

            if (bo_total == 0.0f) continue;

            int type_j = atom_type[j];
            int pair_idx = type_i * atom_type_numbers + type_j;

            s_ovun1 += p_ovun1[pair_idx] * De_s[pair_idx] * bo_total;

            float dfvl_j = (mass[type_j] > 21.0f) ? 0.0f : 1.0f;
            s_ovun2 += (Delta[j] - dfvl_j * Delta_lp_temp[j]) * (b_pi + b_pi2);
        }

        float en_ovun = 0.0f;
        float en_lp = 0.0f;
        float cdd = 0.0f;

        float delta_lp_i = Delta_lp[i];
        float p_lp2_val = p_lp2[type_i];
        float expvd2 = expf(-75.0f * delta_lp_i);
        float inv_expvd2 = 1.0f / (1.0f + expvd2);

        en_lp = p_lp2_val * delta_lp_i * inv_expvd2;
        float dElp = p_lp2_val * inv_expvd2 + 75.0f * p_lp2_val * delta_lp_i *
                                                  expvd2 * inv_expvd2 *
                                                  inv_expvd2;
        cdd += dElp * dDelta_lp[i];

        float exp_ovun1 = p_ovun3 * expf(p_ovun4 * s_ovun2);
        float inv_exp_ovun1 = 1.0f / (1.0f + exp_ovun1);
        float Delta_lpcorr =
            Delta[i] - (dfvl * Delta_lp_temp[i]) * inv_exp_ovun1;

        float p_ovun2_val = p_ovun2[type_i];
        float exp_ovun2 = expf(p_ovun2_val * Delta_lpcorr);
        float inv_exp_ovun2 = 1.0f / (1.0f + exp_ovun2);

        float DlpVi = 1.0f / (Delta_lpcorr + val + 1e-8f);
        float CEover1 = Delta_lpcorr * DlpVi * inv_exp_ovun2;

        float e_ov = s_ovun1 * CEover1;
        en_ovun += e_ov;

        float CEover2 =
            s_ovun1 * DlpVi * inv_exp_ovun2 *
            (1.0f -
             Delta_lpcorr * (DlpVi + p_ovun2_val * exp_ovun2 * inv_exp_ovun2));
        float CEover3 = CEover2 * (1.0f - dfvl * dDelta_lp[i] * inv_exp_ovun1);
        float CEover4 = CEover2 * (dfvl * Delta_lp_temp[i]) * p_ovun4 *
                        exp_ovun1 * inv_exp_ovun1 * inv_exp_ovun1;

        cdd += CEover3;

        float p_ovun5_val = p_ovun5[type_i];
        float exp_ovun2n = 1.0f / exp_ovun2;
        float exp_ovun6 = expf(p_ovun6 * Delta_lpcorr);
        float exp_ovun8 = p_ovun7 * expf(p_ovun8 * s_ovun2);
        float inv_exp_ovun2n = 1.0f / (1.0f + exp_ovun2n);
        float inv_exp_ovun8 = 1.0f / (1.0f + exp_ovun8);

        float e_un =
            -p_ovun5_val * (1.0f - exp_ovun6) * inv_exp_ovun2n * inv_exp_ovun8;
        en_ovun += e_un;

        float CEunder1 = inv_exp_ovun2n *
                         (p_ovun5_val * p_ovun6 * exp_ovun6 * inv_exp_ovun8 +
                          p_ovun2_val * e_un * exp_ovun2n);
        float CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
        float CEunder3 =
            CEunder1 * (1.0f - dfvl * dDelta_lp[i] * inv_exp_ovun1);

        cdd += CEunder3;

        atomicAdd(&CdDelta[i], cdd);

        float CEunder4 = CEunder1 * (dfvl * Delta_lp_temp[i]) * p_ovun4 *
                             exp_ovun1 * inv_exp_ovun1 * inv_exp_ovun1 +
                         CEunder2;
        float CE_sum_4 = CEover4 + CEunder4;

        for (int kk = start; kk < end; kk++)
        {
            int j = bond_nbr[kk];
            int b = bond_idx_arr[kk];
            float b_s = bo_s[b];
            float b_pi = bo_pi[b];
            float b_pi2 = bo_pi2[b];
            float bo_total = b_s + b_pi + b_pi2;

            if (bo_total == 0.0f) continue;

            int type_j = atom_type[j];
            int pair_idx = type_i * atom_type_numbers + type_j;

            float de_dbo_s = CEover1 * p_ovun1[pair_idx] * De_s[pair_idx];
            float de_dbo_pi =
                de_dbo_s + CE_sum_4 * (Delta[j] - dfvl * Delta_lp_temp[j]);
            float de_dbo_pi2 = de_dbo_pi;

            atomicAdd(&d_dE_dBO_s[b], de_dbo_s);
            atomicAdd(&d_dE_dBO_pi[b], de_dbo_pi);
            atomicAdd(&d_dE_dBO_pi2[b], de_dbo_pi2);

            float dfvl_j = (mass[type_j] > 21.0f) ? 0.0f : 1.0f;
            float term =
                CE_sum_4 * (1.0f - dfvl_j * dDelta_lp[j]) * (b_pi + b_pi2);
            atomicAdd(&CdDelta[j], term);
        }

        if (atom_energy) atomicAdd(&atom_energy[i], en_ovun + en_lp);
        atomicAdd(d_energy_ovun_sum, en_ovun);
        atomicAdd(d_energy_elp_sum, en_lp);
    }
}

void REAXFF_OVER_UNDER::Initial(CONTROLLER* controller, int atom_numbers,
                                const char* module_name)
{
    if (module_name == NULL) module_name = "REAXFF";
    if (!controller->Command_Exist(module_name, "in_file")) return;

    controller->printf("START INITIALIZING REAXFF OVER/UNDER COORD\n");

    const char* parameter_in_file =
        controller->Original_Command(module_name, "in_file");
    const char* type_in_file =
        controller->Original_Command(module_name, "type_in_file");
    if (parameter_in_file == NULL || type_in_file == NULL)
    {
        controller->printf(
            "REAXFF_OVER_UNDER IS NOT INITIALIZED (missing input files)\n\n");
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
            "REAXFF_OVER_UNDER::Initial", "Reason:\n\t%s", reason.c_str());
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
            "REAXFF_OVER_UNDER::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }
    if (force_field.general_parameters.size() <= 32)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_OVER_UNDER::Initial",
            "Reason:\n\tgeneral parameter count is too small for "
            "over/under coordination in file %s",
            parameter_in_file);
        return;
    }

    const int n_atom_types = static_cast<int>(force_field.atom_types.size());
    int pair_parameter_count = 0;
    if (!ReaxFF_Checked_Dense_Table_Count(
            n_atom_types, 2, sizeof(float), &pair_parameter_count))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_OVER_UNDER::Initial",
            "Reason:\n\tatom type count %d exceeds the supported pair-table "
            "extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    std::vector<float> valency(n_atom_types);
    std::vector<float> valency_e(n_atom_types);
    std::vector<float> valency_boc(n_atom_types);
    std::vector<float> valency_val(n_atom_types);
    std::vector<float> mass(n_atom_types);
    std::vector<float> p_lp2_values(n_atom_types);
    std::vector<float> p_ovun2_values(n_atom_types);
    std::vector<float> p_ovun5_values(n_atom_types);
    for (int i = 0; i < n_atom_types; i++)
    {
        const REAXFF_ATOM_TYPE_IR& atom = force_field.atom_types[i];
        valency[i] = atom.values[0][1];
        mass[i] = atom.values[0][2];
        valency_e[i] = atom.values[0][7];
        valency_boc[i] = atom.values[1][2];
        p_ovun5_values[i] = atom.values[1][3];
        p_lp2_values[i] = atom.values[2][1];
        p_ovun2_values[i] = atom.values[3][0];
        valency_val[i] = atom.values[3][3];
    }
    std::vector<float> p_ovun1_values(pair_parameter_count, 0.0f);
    std::vector<float> de_s_values(pair_parameter_count, 0.0f);
    for (const REAXFF_BOND_IR& bond : force_field.bonds)
    {
        const int idx1 = bond.type1 - 1;
        const int idx2 = bond.type2 - 1;
        de_s_values[idx1 * n_atom_types + idx2] =
            de_s_values[idx2 * n_atom_types + idx1] = bond.line1[0];
        p_ovun1_values[idx1 * n_atom_types + idx2] =
            p_ovun1_values[idx2 * n_atom_types + idx1] = bond.line1[7];
    }

    auto allocate_float_copy = [](const std::vector<float>& values)
    {
        float* result = NULL;
        Malloc_Safely((void**)&result, sizeof(float) * values.size());
        memcpy(result, values.data(), sizeof(float) * values.size());
        return result;
    };
    float* staged_valency = allocate_float_copy(valency);
    float* staged_valency_e = allocate_float_copy(valency_e);
    float* staged_valency_boc = allocate_float_copy(valency_boc);
    float* staged_valency_val = allocate_float_copy(valency_val);
    float* staged_mass = allocate_float_copy(mass);
    float* staged_p_lp2 = allocate_float_copy(p_lp2_values);
    float* staged_p_ovun2 = allocate_float_copy(p_ovun2_values);
    float* staged_p_ovun5 = allocate_float_copy(p_ovun5_values);
    float* staged_p_ovun1 = allocate_float_copy(p_ovun1_values);
    float* staged_de_s = allocate_float_copy(de_s_values);
    int* staged_atom_type = NULL;
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);

    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    p_lp1 = force_field.general_parameters[15];
    p_lp3 = force_field.general_parameters[5];
    p_ovun3 = force_field.general_parameters[32];
    p_ovun4 = force_field.general_parameters[31];
    p_ovun6 = force_field.general_parameters[6];
    p_ovun7 = force_field.general_parameters[8];
    p_ovun8 = force_field.general_parameters[9];
    h_valency = staged_valency;
    h_valency_e = staged_valency_e;
    h_valency_boc = staged_valency_boc;
    h_valency_val = staged_valency_val;
    h_mass = staged_mass;
    h_p_lp2 = staged_p_lp2;
    h_p_ovun2 = staged_p_ovun2;
    h_p_ovun5 = staged_p_ovun5;
    h_p_ovun1 = staged_p_ovun1;
    h_De_s = staged_de_s;
    h_atom_type = staged_atom_type;

    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_p_lp2, h_p_lp2,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_p_ovun2, h_p_ovun2,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_p_ovun5, h_p_ovun5,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_valency, h_valency,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_valency_e, h_valency_e,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_valency_boc, h_valency_boc,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_valency_val, h_valency_val,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_mass, h_mass,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_p_ovun1, h_p_ovun1,
                                  sizeof(float) * n_atom_types * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_De_s, h_De_s,
                                  sizeof(float) * n_atom_types * n_atom_types);

    Device_Malloc_Safely((void**)&d_Delta, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_Delta_boc, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_Delta_val, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_Delta_lp, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_nlp, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_vlpex, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_Delta_lp_temp,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_dDelta_lp, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_CdDelta, sizeof(float) * atom_numbers);

    Device_Malloc_Safely((void**)&d_energy_elp_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_energy_ovun_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_energy_atom, sizeof(float) * atom_numbers);

    controller->Step_Print_Initial("REAXFF_ELP", "%14.7e");
    controller->Step_Print_Initial("REAXFF_OVUN", "%14.7e");
    is_initialized = 1;
    controller->printf("END INITIALIZING REAXFF OVER/UNDER COORD\n\n");
}

void REAXFF_OVER_UNDER::Calculate_Over_Under_Energy_And_Force(
    int atom_numbers, const VECTOR* crd, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, REAXFF_BOND_ORDER* bo_module,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_virial)
{
    if (!is_initialized) return;

    dim3 blockSize(128);
    dim3 gridSize((atom_numbers + blockSize.x - 1) / blockSize.x);

    Launch_Device_Kernel(Calculate_Delta_Kernel, gridSize, blockSize, 0, NULL,
                         atom_numbers, d_atom_type,
                         bo_module->d_total_corrected_bond_order, d_valency,
                         d_valency_e, d_valency_boc, d_valency_val, d_mass,
                         p_lp1, d_Delta, d_Delta_boc, d_Delta_val, d_Delta_lp,
                         d_nlp, d_vlpex, d_Delta_lp_temp, d_dDelta_lp);

    deviceMemset(d_energy_elp_sum, 0, sizeof(float));
    deviceMemset(d_energy_ovun_sum, 0, sizeof(float));

    Launch_Device_Kernel(
        Calculate_Energy_Force_Prep_Kernel, gridSize, blockSize, 0, NULL,
        atom_numbers, d_atom_type, d_mass, d_Delta, d_Delta_lp, d_Delta_lp_temp,
        d_dDelta_lp, bo_module->d_corrected_bo_s, bo_module->d_corrected_bo_pi,
        bo_module->d_corrected_bo_pi2, d_p_ovun1, d_De_s, atom_type_numbers,
        d_p_lp2, d_valency, p_ovun3, p_ovun4, p_ovun6, p_ovun7, p_ovun8,
        d_p_ovun2, d_p_ovun5, d_dE_dBO_s, d_dE_dBO_pi, d_dE_dBO_pi2, d_CdDelta,
        need_atom_energy ? atom_energy : NULL, d_energy_ovun_sum,
        d_energy_elp_sum, bo_module->d_bond_count, bo_module->d_bond_offset,
        bo_module->d_bond_nbr, bo_module->d_bond_idx);
}

void REAXFF_OVER_UNDER::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_ovun, d_energy_ovun_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print("REAXFF_OVUN", h_energy_ovun, true);
}

void REAXFF_OVER_UNDER::Step_Print_ELP(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_lp, d_energy_elp_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print("REAXFF_ELP", h_energy_lp, true);
}
