#include "bond.h"

#include "reaxff_input.h"

#include "bond_order.h"  // for find_bond_index

static const int De_s = 0;
static const int De_p = 1;
static const int De_PP = 2;
static const int p_be1 = 3;
static const int p_be2 = 4;
static const int PARAM_STRIDE = 5;

template <int N>
__device__ __forceinline__ SADfloat<N> reax_bond_energy_sad(SADfloat<N> BO_s,
                                                            SADfloat<N> BO_pi,
                                                            SADfloat<N> BO_pi2,
                                                            const float* param)
{
    float De_s_val = param[De_s];
    float De_p_val = param[De_p];
    float De_PP_val = param[De_PP];
    float p_be1_val = param[p_be1];
    float p_be2_val = param[p_be2];

    SADfloat<N> pow_BOs_be2 = powf(BO_s, p_be2_val);
    SADfloat<N> exp_be12 = expf(p_be1_val * (1.0f - pow_BOs_be2));

    SADfloat<N> ebond =
        -De_s_val * BO_s * exp_be12 - De_p_val * BO_pi - De_PP_val * BO_pi2;
    return ebond;
}

static __global__ void REAXFF_Bond_Force_CUDA(
    const int atom_numbers, const VECTOR* crd, VECTOR* frc,
    const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
    int* atom_types, float* params, int ntypes, float* bo_s, float* bo_pi,
    float* bo_pi2, float* d_dE_dBO_s, float* d_dE_dBO_pi, float* d_dE_dBO_pi2,
    float* atom_energy, LTMatrix3* atom_virial, float* d_energy_sum,
    const int* bond_count, const int* bond_offset, const int* bond_nbr,
    const int* bond_idx_arr)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type_i = atom_types[i];
        VECTOR ri = crd[i];
        ATOM_GROUP nl_i = nl[i];
        float en_i = 0;

        for (int jj = 0; jj < nl_i.atom_numbers; jj++)
        {
            int j = nl_i.atom_serial[jj];

            if (j <= i) continue;

            int b = find_bond_index(i, j, bond_count, bond_offset, bond_nbr,
                                    bond_idx_arr);
            if (b < 0) continue;

            int param_idx = type_i * ntypes + atom_types[j];
            const float* param = params + param_idx * PARAM_STRIDE;

            float BO_s_ij = bo_s[b];
            float BO_pi_ij = bo_pi[b];
            float BO_pi2_ij = bo_pi2[b];

            if (BO_s_ij + BO_pi_ij + BO_pi2_ij == 0.0f) continue;

            SADfloat<3> BO_s_sad(BO_s_ij, 0);
            SADfloat<3> BO_pi_sad(BO_pi_ij, 1);
            SADfloat<3> BO_pi2_sad(BO_pi2_ij, 2);

            SADfloat<3> energy_sad =
                reax_bond_energy_sad(BO_s_sad, BO_pi_sad, BO_pi2_sad, param);

            atomicAdd(&d_dE_dBO_s[b], energy_sad.dval[0]);
            atomicAdd(&d_dE_dBO_pi[b], energy_sad.dval[1]);
            atomicAdd(&d_dE_dBO_pi2[b], energy_sad.dval[2]);

            if (atom_energy)
            {
                en_i += energy_sad.val;
                atomicAdd(d_energy_sum, energy_sad.val);
            }
        }

        if (atom_energy)
        {
            atom_energy[i] += en_i;
        }
    }
}

void REAXFF_BOND::Initial(CONTROLLER* controller, int atom_numbers,
                          const char* module_name, bool* need_full_nl_flag)
{
    if (module_name == NULL) module_name = "REAXFF";
    if (!controller->Command_Exist(module_name, "in_file")) return;

    controller->printf("START INITIALIZING REAXFF BOND FORCE\n");
    const char* parameter_in_file =
        controller->Original_Command(module_name, "in_file");
    const char* type_in_file =
        controller->Original_Command(module_name, "type_in_file");
    if (parameter_in_file == NULL || type_in_file == NULL)
    {
        controller->printf(
            "REAXFF_BOND IS NOT INITIALIZED (missing input files)\n\n");
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
            "REAXFF_BOND::Initial", "Reason:\n\t%s", reason.c_str());
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
            "REAXFF_BOND::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }

    const int n_atom_types = static_cast<int>(force_field.atom_types.size());
    int pair_parameter_count = 0;
    std::size_t parameter_count = 0;
    std::size_t parameter_bytes = 0;
    if (!ReaxFF_Checked_Dense_Table_Count(
            n_atom_types, 2, sizeof(float) * PARAM_STRIDE,
            &pair_parameter_count) ||
        !ReaxFF_Checked_Size_Multiply(
            static_cast<std::size_t>(pair_parameter_count),
            static_cast<std::size_t>(PARAM_STRIDE), &parameter_count) ||
        !ReaxFF_Checked_Size_Multiply(parameter_count, sizeof(float),
                                      &parameter_bytes))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_BOND::Initial",
            "Reason:\n\tatom type count %d exceeds the supported bond "
            "parameter table extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    std::vector<float> twobody_params(parameter_count, 0.0f);
    for (const REAXFF_BOND_IR& bond : force_field.bonds)
    {
        const int idx1 = bond.type1 - 1;
        const int idx2 = bond.type2 - 1;
        const int pair_idx1 = (idx1 * n_atom_types + idx2) * PARAM_STRIDE;
        const int pair_idx2 = (idx2 * n_atom_types + idx1) * PARAM_STRIDE;
        twobody_params[pair_idx1 + De_s] =
            twobody_params[pair_idx2 + De_s] = bond.line1[0];
        twobody_params[pair_idx1 + De_p] =
            twobody_params[pair_idx2 + De_p] = bond.line1[1];
        twobody_params[pair_idx1 + De_PP] =
            twobody_params[pair_idx2 + De_PP] = bond.line1[2];
        twobody_params[pair_idx1 + p_be1] =
            twobody_params[pair_idx2 + p_be1] = bond.line1[3];
        twobody_params[pair_idx1 + p_be2] =
            twobody_params[pair_idx2 + p_be2] = bond.line2[0];
    }

    float* staged_twobody_params = NULL;
    int* staged_atom_type = NULL;
    Malloc_Safely((void**)&staged_twobody_params, parameter_bytes);
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    memcpy(staged_twobody_params, twobody_params.data(), parameter_bytes);
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);

    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    h_twobody_params = staged_twobody_params;
    h_atom_type = staged_atom_type;

    Device_Malloc_And_Copy_Safely((void**)&d_twobody_params,
                                  h_twobody_params, parameter_bytes);

    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_energy_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_energy_atom, sizeof(float) * atom_numbers);
    deviceMemset(d_energy_sum, 0, sizeof(float));
    deviceMemset(d_energy_atom, 0, sizeof(float) * atom_numbers);

    if (need_full_nl_flag != NULL) *need_full_nl_flag = true;

    is_initialized = true;
    controller->Step_Print_Initial("REAXFF_BOND", "%14.7e");
    controller->printf("END INITIALIZING REAXFF BOND FORCE\n\n");
}

void REAXFF_BOND::REAXFF_Bond_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const VECTOR* crd, VECTOR* frc,
    const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_virial)
{
    if (!is_initialized) return;

    if (need_atom_energy)
    {
        deviceMemset(d_energy_sum, 0, sizeof(float));
        if (atom_energy)
            deviceMemset(d_energy_atom, 0, sizeof(float) * atom_numbers);
    }

    dim3 blockSize(128);
    dim3 gridSize((atom_numbers + blockSize.x - 1) / blockSize.x);

    Launch_Device_Kernel(REAXFF_Bond_Force_CUDA, gridSize, blockSize, 0, NULL,
                         atom_numbers, crd, frc, cell, rcell, nl, d_atom_type,
                         d_twobody_params, atom_type_numbers, d_bo_s, d_bo_pi,
                         d_bo_pi2, d_dE_dBO_s, d_dE_dBO_pi, d_dE_dBO_pi2,
                         atom_energy, atom_virial, d_energy_sum, d_bond_count,
                         d_bond_offset, d_bond_nbr, d_bond_idx);
}

void REAXFF_BOND::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_sum, d_energy_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print("REAXFF_BOND", h_energy_sum, true);
}
