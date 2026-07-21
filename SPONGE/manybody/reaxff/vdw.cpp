#include "vdw.h"

#include "reaxff_geometry.h"
#include "reaxff_input.h"

static const int p_rvdw = 0;
static const int p_epsilon = 1;
static const int p_alpha = 2;
static const int p_gamma_w = 3;
static const int PARAM_STRIDE = 8;

template <int N>
__device__ __forceinline__ SADfloat<N> reax_vdw_energy_sad(SADfloat<N> r,
                                                           const float* param,
                                                           float cutoff,
                                                           float p_vdw1)
{
    float rvdw = param[p_rvdw];
    float epsilon = param[p_epsilon];
    float alpha = param[p_alpha];
    float gamma_w = param[p_gamma_w];

    if (r.val > cutoff) return SADfloat<N>(0.0f);
    if (epsilon == 0.0f) return SADfloat<N>(0.0f);

    SADfloat<N> x = r / SADfloat<N>(cutoff);
    SADfloat<N> x4 = x * x * x * x;
    SADfloat<N> x5 = x4 * x;
    SADfloat<N> x6 = x5 * x;
    SADfloat<N> x7 = x6 * x;

    SADfloat<N> tap = SADfloat<N>(1.0f) - SADfloat<N>(35.0f) * x4 +
                      SADfloat<N>(84.0f) * x5 - SADfloat<N>(70.0f) * x6 +
                      SADfloat<N>(20.0f) * x7;

    float inv_gamma = 1.0f / gamma_w;
    SADfloat<N> inv_gamma_p = powf(SADfloat<N>(inv_gamma), SADfloat<N>(p_vdw1));
    SADfloat<N> r_p = powf(r, SADfloat<N>(p_vdw1));
    SADfloat<N> shielded_r =
        powf(r_p + inv_gamma_p, SADfloat<N>(1.0f / p_vdw1));

    SADfloat<N> exp_term = alpha * (1.0f - shielded_r / rvdw);
    SADfloat<N> term1 = expf(exp_term);
    SADfloat<N> term2 = -2.0f * expf(0.5f * exp_term);

    return tap * epsilon * (term1 + term2);
}

static __global__ void REAXFF_VDW_Force_CUDA(
    const int atom_numbers, const VECTOR* crd, VECTOR* frc,
    const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
    const int* atom_types, const float* params, int ntypes, float cutoff,
    float p_vdw1, float* atom_energy, LTMatrix3* atom_virial,
    float* d_energy_sum, int* geometry_error)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type_i = atom_types[i];
        if (type_i < 0 || type_i >= ntypes)
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
        if (!ReaxFF_Vector_Is_Finite(ri))
        {
            Record_ReaxFF_Geometry_Error(geometry_error, REAXFF_VDW_NONFINITE,
                                         i);
#ifdef GPU_ARCH_NAME
            return;
#else
            continue;
#endif
        }
        ATOM_GROUP nl_i = nl[i];
        VECTOR fi = {0, 0, 0};
        LTMatrix3 vi = {0, 0, 0, 0, 0, 0};
        float en_i = 0;

        for (int jj = 0; jj < nl_i.atom_numbers; jj++)
        {
            int j = nl_i.atom_serial[jj];
            if (j <= i) continue;
            int type_j = atom_types[j];
            if (type_j < 0 || type_j >= ntypes)
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_INVALID_ATOM_TYPE, i, j,
                                             type_i, type_j);
                continue;
            }

            VECTOR rj = crd[j];
            VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
            if (!ReaxFF_Vector_Is_Finite(drij))
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_VDW_NONFINITE, i, j);
                continue;
            }
            float rij = norm3df(drij.x, drij.y, drij.z);
            if (!ReaxFF_Float_Is_Finite(rij) || rij < 0.0f)
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_VDW_NONFINITE, i, j);
                continue;
            }

            if (rij >= cutoff) continue;

            int param_idx = (type_i * ntypes + type_j) * PARAM_STRIDE;
            const float* param = params + param_idx;

            SADfloat<1> rij_sad(rij, 0);
            SADfloat<1> energy_sad =
                reax_vdw_energy_sad(rij_sad, param, cutoff, p_vdw1);

            float force_mag = 0.0f;
            VECTOR fij = {0.0f, 0.0f, 0.0f};
            if (rij == 0.0f)
            {
                // With p_vdw1 > 1, dE/dr = O(r^(p_vdw1-1)) and the
                // Cartesian force has the unique zero limit.  For another
                // positive exponent, an exact overlap is valid only if the
                // evaluated radial derivative itself vanishes.
                if (!(p_vdw1 > 1.0f) && energy_sad.dval[0] != 0.0f)
                {
                    Record_ReaxFF_Geometry_Error(geometry_error,
                                                 REAXFF_VDW_OVERLAP, i, j);
                    continue;
                }
            }
            else
            {
                force_mag = -energy_sad.dval[0] / rij;
                fij = {force_mag * drij.x, force_mag * drij.y,
                       force_mag * drij.z};
            }

            LTMatrix3 pair_virial = {0, 0, 0, 0, 0, 0};
            if (atom_virial)
            {
                pair_virial = Get_Virial_From_Force_Dis(fij, drij);
            }
            const VECTOR next_fi = {fi.x + fij.x, fi.y + fij.y, fi.z + fij.z};
            const LTMatrix3 next_vi = vi + pair_virial;
            const float next_en_i = en_i + energy_sad.val;
            if (!ReaxFF_Float_Is_Finite(energy_sad.val) ||
                !ReaxFF_Float_Is_Finite(energy_sad.dval[0]) ||
                !ReaxFF_Float_Is_Finite(force_mag) ||
                !ReaxFF_Vector_Is_Finite(fij) ||
                !ReaxFF_Vector_Is_Finite(next_fi) ||
                (atom_virial && (!ReaxFF_Matrix_Is_Finite(pair_virial) ||
                                 !ReaxFF_Matrix_Is_Finite(next_vi))) ||
                (atom_energy && !ReaxFF_Float_Is_Finite(next_en_i)))
            {
                Record_ReaxFF_Geometry_Error(geometry_error,
                                             REAXFF_VDW_NONFINITE, i, j);
                continue;
            }

            fi = next_fi;

            atomicAdd(&frc[j].x, -fij.x);
            atomicAdd(&frc[j].y, -fij.y);
            atomicAdd(&frc[j].z, -fij.z);

            if (atom_virial)
            {
                vi = next_vi;
            }

            if (atom_energy)
            {
                en_i = next_en_i;
                atomicAdd(d_energy_sum, energy_sad.val);
            }
        }

        atomicAdd(&frc[i].x, fi.x);
        atomicAdd(&frc[i].y, fi.y);
        atomicAdd(&frc[i].z, fi.z);

        if (atom_energy)
        {
            atom_energy[i] += en_i;
        }
        if (atom_virial)
        {
            atomicAdd(atom_virial + i, vi);
        }
    }
}

void REAXFF_VDW::Initial(CONTROLLER* controller, int atom_numbers,
                         const char* module_name, bool* need_full_nl_flag)
{
    if (module_name == NULL) module_name = "REAXFF";
    if (!controller->Command_Exist(module_name, "in_file")) return;

    controller->printf("START INITIALIZING REAXFF VDW FORCE\n");
    const char* parameter_in_file =
        controller->Original_Command(module_name, "in_file");
    const char* type_in_file =
        controller->Original_Command(module_name, "type_in_file");
    if (parameter_in_file == NULL || type_in_file == NULL)
    {
        controller->printf(
            "REAXFF_VDW IS NOT INITIALIZED (missing input files)\n\n");
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
            "REAXFF_VDW::Initial", "Reason:\n\t%s", reason.c_str());
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
            "REAXFF_VDW::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }
    if (force_field.general_parameters.size() <= 28)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_VDW::Initial",
            "Reason:\n\tmissing general parameter p_vdw1 at index 29 in "
            "file %s",
            parameter_in_file);
        return;
    }
    const float staged_p_vdw1 = force_field.general_parameters[28];
    if (!(staged_p_vdw1 > 0.0f) ||
        !Float_Memory_Is_Finite(&staged_p_vdw1))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_VDW::Initial",
            "Reason:\n\tgeneral parameter p_vdw1 must be positive and "
            "finite in file %s",
            parameter_in_file);
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
            spongeErrorBadFileFormat, "REAXFF_VDW::Initial",
            "Reason:\n\tatom type count %d exceeds the supported van der "
            "Waals parameter table extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    std::vector<float> twobody_params(parameter_count, 0.0f);
    auto geometric_mix = [](float lhs, float rhs, double scale)
    {
        const double product =
            static_cast<double>(lhs) * static_cast<double>(rhs);
        return static_cast<float>(scale * sqrt(product));
    };
    for (int i = 0; i < n_atom_types; i++)
    {
        for (int j = 0; j < n_atom_types; j++)
        {
            const REAXFF_ATOM_TYPE_IR& atom_i = force_field.atom_types[i];
            const REAXFF_ATOM_TYPE_IR& atom_j = force_field.atom_types[j];
            const float rvdw_ij = geometric_mix(
                atom_i.values[0][3], atom_j.values[0][3], 2.0);
            const float epsilon_ij = geometric_mix(
                atom_i.values[0][4], atom_j.values[0][4], 1.0);
            const float alpha_ij = geometric_mix(
                atom_i.values[1][0], atom_j.values[1][0], 1.0);
            const float gamma_w_ij = geometric_mix(
                atom_i.values[1][1], atom_j.values[1][1], 1.0);
            if (!Float_Memory_Is_Finite(&rvdw_ij) ||
                !Float_Memory_Is_Finite(&epsilon_ij) ||
                !Float_Memory_Is_Finite(&alpha_ij) ||
                !Float_Memory_Is_Finite(&gamma_w_ij) ||
                (epsilon_ij != 0.0f &&
                 (!(rvdw_ij > 0.0f) || !(gamma_w_ij > 0.0f))))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "REAXFF_VDW::Initial",
                    "Reason:\n\tatom parameters produce an invalid active "
                    "van der Waals pair for types %d/%d in file %s",
                    i + 1, j + 1, parameter_in_file);
                return;
            }
            const int idx = (i * n_atom_types + j) * PARAM_STRIDE;
            twobody_params[idx + p_rvdw] = rvdw_ij;
            twobody_params[idx + p_epsilon] = epsilon_ij;
            twobody_params[idx + p_alpha] = alpha_ij;
            twobody_params[idx + p_gamma_w] = gamma_w_ij;
        }
    }
    for (const REAXFF_OFF_DIAGONAL_IR& entry : force_field.off_diagonal)
    {
        const int idx1 = entry.type1 - 1;
        const int idx2 = entry.type2 - 1;
        const int pair_idx1 = (idx1 * n_atom_types + idx2) * PARAM_STRIDE;
        const int pair_idx2 = (idx2 * n_atom_types + idx1) * PARAM_STRIDE;
        if (entry.values[0] > 0.0f)
            twobody_params[pair_idx1 + p_epsilon] =
                twobody_params[pair_idx2 + p_epsilon] = entry.values[0];
        if (entry.values[1] > 0.0f)
            twobody_params[pair_idx1 + p_rvdw] =
                twobody_params[pair_idx2 + p_rvdw] = 2.0f * entry.values[1];
        if (entry.values[2] > 0.0f)
            twobody_params[pair_idx1 + p_alpha] =
                twobody_params[pair_idx2 + p_alpha] = entry.values[2];
    }
    for (int pair = 0; pair < pair_parameter_count; pair++)
    {
        const float* params = twobody_params.data() + pair * PARAM_STRIDE;
        if (!Float_Memory_Is_Finite(params + p_rvdw) ||
            !Float_Memory_Is_Finite(params + p_epsilon) ||
            !Float_Memory_Is_Finite(params + p_alpha) ||
            !Float_Memory_Is_Finite(params + p_gamma_w) ||
            (params[p_epsilon] != 0.0f &&
             (!(params[p_rvdw] > 0.0f) || !(params[p_gamma_w] > 0.0f))))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "REAXFF_VDW::Initial",
                "Reason:\n\tfinal van der Waals parameters for pair index %d "
                "are invalid in file %s",
                pair + 1, parameter_in_file);
            return;
        }
    }

    float* staged_twobody_params = NULL;
    int* staged_atom_type = NULL;
    Malloc_Safely((void**)&staged_twobody_params, parameter_bytes);
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    memcpy(staged_twobody_params, twobody_params.data(), parameter_bytes);
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);

    this->controller = controller;
    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    this->p_vdw1 = staged_p_vdw1;
    h_twobody_params = staged_twobody_params;
    h_atom_type = staged_atom_type;

    Device_Malloc_And_Copy_Safely((void**)&d_twobody_params,
                                  h_twobody_params, parameter_bytes);

    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_energy_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_energy_atom, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_geometry_error,
                         sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);
    deviceMemset(d_energy_sum, 0, sizeof(float));
    deviceMemset(d_energy_atom, 0, sizeof(float) * atom_numbers);

    if (need_full_nl_flag != NULL) *need_full_nl_flag = true;
    is_initialized = true;
    controller->Step_Print_Initial("REAXFF_VDW", "%14.7e");
    controller->printf("END INITIALIZING REAXFF VDW FORCE\n\n");
}

void REAXFF_VDW::REAXFF_VDW_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const VECTOR* crd, VECTOR* frc,
    const LTMatrix3 cell, const LTMatrix3 rcell, const ATOM_GROUP* nl,
    const float cutoff, const int need_atom_energy, float* atom_energy,
    const int need_virial, LTMatrix3* atom_virial)
{
    if (!is_initialized) return;

    float checked_cutoff = cutoff;
    if (!(checked_cutoff > 0.0f) || !Float_Memory_Is_Finite(&checked_cutoff))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF_VDW::REAXFF_VDW_Force_With_Atom_Energy_And_Virial",
            "The van der Waals cutoff must be positive and finite, but got "
            "%g.",
            cutoff);
        return;
    }

    if (need_atom_energy)
    {
        deviceMemset(d_energy_sum, 0, sizeof(float));
        if (atom_energy)
            deviceMemset(d_energy_atom, 0, sizeof(float) * atom_numbers);
    }

    dim3 blockSize(128);
    dim3 gridSize((atom_numbers + blockSize.x - 1) / blockSize.x);
    deviceMemset(d_geometry_error, 0, sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);

    Launch_Device_Kernel(
        REAXFF_VDW_Force_CUDA, gridSize, blockSize, 0, NULL, atom_numbers, crd,
        frc, cell, rcell, nl, d_atom_type, d_twobody_params, atom_type_numbers,
        cutoff, this->p_vdw1, atom_energy, need_virial ? atom_virial : NULL,
        d_energy_sum, d_geometry_error);
    Check_Geometry_Error();
}

void REAXFF_VDW::Check_Geometry_Error()
{
    int error[REAXFF_GEOMETRY_ERROR_SIZE] = {0, -1, -1, -1, -1};
    deviceMemcpy(error, d_geometry_error, sizeof(error),
                 deviceMemcpyDeviceToHost);
    if (error[0] == REAXFF_GEOMETRY_OK) return;

    if (error[0] == REAXFF_INVALID_ATOM_TYPE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF_VDW::REAXFF_VDW_Force_With_Atom_Energy_And_Virial",
            "Invalid ReaxFF atom type while evaluating local atoms %d %d "
            "(types %d %d; valid range is [0, %d)).",
            error[1], error[2], error[3], error[4], atom_type_numbers);
        return;
    }

    const char* reason =
        error[0] == REAXFF_VDW_OVERLAP
            ? "overlap exactly while the shielded radial derivative is "
              "nonzero; the Cartesian force direction is undefined"
            : "produce a non-finite or unrepresentable ReaxFF van der Waals "
              "energy, derivative, force, or virial";
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown,
        "REAXFF_VDW::REAXFF_VDW_Force_With_Atom_Energy_And_Virial",
        "Reason:\n\tlocal atoms %d %d %s\n", error[1], error[2], reason);
}

void REAXFF_VDW::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_sum, d_energy_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print("REAXFF_VDW", h_energy_sum, true);
}
