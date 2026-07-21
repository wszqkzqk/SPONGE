#include "hydrogen_bond.h"

#include "bond_order.h"  // for find_bond_index
#include "reaxff_geometry.h"
#include "reaxff_input.h"

static __global__ void Calculate_HB_Kernel(
    int atom_numbers, const VECTOR* crd, const int* atom_type,
    const int* is_hydrogen, const REAXFF_HB_Info* hb_info,
    const REAXFF_HB_Entry* hb_entries, int atom_type_numbers, const float* bo_s,
    const float* bo_pi, const float* bo_pi2, float* d_dE_dBO_s,
    float* d_dE_dBO_pi, float* d_dE_dBO_pi2, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, float* atom_energy,
    VECTOR* frc, LTMatrix3* atom_virial, float* d_energy_hb_sum,
    const int* bond_count, const int* bond_offset, const int* bond_nbr,
    const int* bond_idx_arr, int* geometry_error)
{
    SIMPLE_DEVICE_FOR(h, atom_numbers)
    {
        if (is_hydrogen[h])
        {
            int type_h = atom_type[h];
            if (type_h < 0 || type_h >= atom_type_numbers)
            {
                Record_ReaxFF_Geometry_Error(
                    geometry_error, REAXFF_INVALID_ATOM_TYPE, h, -1, type_h);
#ifdef GPU_ARCH_NAME
                return;
#else
                continue;
#endif
            }
            VECTOR rh = crd[h];
            ATOM_GROUP nl_h = nl[h];

            double en_hb = 0.0;

            int bc_h = bond_count[h];
            int bo_h = bond_offset[h];
            for (int pd = 0; pd < bc_h; pd++)
            {
                int b_dh = bond_idx_arr[bo_h + pd];
                int d = bond_nbr[bo_h + pd];
                int type_d = atom_type[d];
                if (type_d < 0 || type_d >= atom_type_numbers)
                {
                    Record_ReaxFF_Geometry_Error(
                        geometry_error, REAXFF_INVALID_ATOM_TYPE, d, h, type_d);
                    continue;
                }
                float bo_dh_val = bo_s[b_dh] + bo_pi[b_dh] + bo_pi2[b_dh];
                if (!ReaxFF_Float_Is_Finite(bo_dh_val))
                {
                    Record_ReaxFF_Geometry_Error(geometry_error,
                                                 REAXFF_HB_NONFINITE, d, h);
                    continue;
                }
                if (bo_dh_val < 0.01f) continue;

                VECTOR rd = crd[d];
                VECTOR ddh = Get_Periodic_Displacement(rd, rh, cell, rcell);

                for (int pa = 0; pa < nl_h.atom_numbers; pa++)
                {
                    int a = nl_h.atom_serial[pa];
                    if (a == d) continue;
                    int type_a = atom_type[a];
                    if (type_a < 0 || type_a >= atom_type_numbers)
                    {
                        Record_ReaxFF_Geometry_Error(geometry_error,
                                                     REAXFF_INVALID_ATOM_TYPE,
                                                     d, h, a, type_a);
                        continue;
                    }

                    int hb_idx = ((type_d * atom_type_numbers + type_h) *
                                      atom_type_numbers +
                                  type_a);
                    REAXFF_HB_Info info = hb_info[hb_idx];
                    if (info.entry_count == 0) continue;

                    bool interaction_active = false;
                    for (int e = 0; e < info.entry_count; e++)
                    {
                        interaction_active =
                            interaction_active ||
                            hb_entries[info.start_idx + e].p_hb1 != 0.0f;
                    }

                    VECTOR ra = crd[a];
                    VECTOR dah = Get_Periodic_Displacement(ra, rh, cell, rcell);
                    REAXFF_ANGLE_GEOMETRY geometry;
                    if (!Compute_ReaxFF_Angle_Geometry(ddh, dah, &geometry))
                    {
                        if (interaction_active)
                        {
                            Record_ReaxFF_Geometry_Error(
                                geometry_error, REAXFF_HB_UNDEFINED, d, h, a);
                        }
                        continue;
                    }
                    const float r_dh = norm3df(ddh.x, ddh.y, ddh.z);
                    float r_ah = norm3df(dah.x, dah.y, dah.z);
                    if (!(r_dh > 0.0f) || !(r_ah > 0.0f) ||
                        !ReaxFF_Float_Is_Finite(r_dh) ||
                        !ReaxFF_Float_Is_Finite(r_ah))
                    {
                        if (interaction_active)
                        {
                            Record_ReaxFF_Geometry_Error(
                                geometry_error, REAXFF_HB_UNDEFINED, d, h, a);
                        }
                        continue;
                    }
                    if (r_ah > 7.5f) continue;

                    const float half_sine_squared =
                        geometry.half_sine_squared;
                    // The geometry helper forms sin(theta/2)^2 from the cross
                    // product near zero, where 1-cos(theta) would round away.
                    const float sin_p4 =
                        half_sine_squared * half_sine_squared;

                    for (int e = 0; e < info.entry_count; e++)
                    {
                        const REAXFF_HB_Entry* param =
                            &hb_entries[info.start_idx + e];
                        if (param->p_hb1 == 0.0f) continue;

                        SADfloat<1> s_bo_dh(bo_dh_val, 0);
                        SADfloat<1> s_f_hb =
                            1.0f - expf(-(float)param->p_hb2 * s_bo_dh);

                        const double r_ah_double =
                            static_cast<double>(r_ah);
                        const double r0_hb =
                            static_cast<double>(param->r0_hb);
                        const double p_hb3 =
                            static_cast<double>(param->p_hb3);
                        const double radial_shape =
                            r0_hb / r_ah_double + r_ah_double / r0_hb - 2.0;
                        const double exp_hb3_double =
                            exp(-p_hb3 * radial_shape);
                        const float exp_hb3 =
                            static_cast<float>(exp_hb3_double);

                        SADfloat<1> s_en_total =
                            (float)param->p_hb1 * s_f_hb * exp_hb3 * sin_p4;

                        const double radial_log_derivative =
                            -p_hb3 *
                            (-r0_hb / (r_ah_double * r_ah_double) +
                             1.0 / r0_hb);
                        const double dE_dr_ah_double =
                            static_cast<double>(param->p_hb1) *
                            static_cast<double>(s_f_hb.val) *
                            static_cast<double>(sin_p4) * exp_hb3_double *
                            radial_log_derivative;
                        const float dE_dr_ah =
                            static_cast<float>(dE_dr_ah_double);

                        float f_ah = -dE_dr_ah;
                        VECTOR f_a_rad = {f_ah * dah.x / r_ah,
                                          f_ah * dah.y / r_ah,
                                          f_ah * dah.z / r_ah};

                        float dE_dsinp4 =
                            (float)param->p_hb1 * s_f_hb.val * exp_hb3;
                        // d(sin(theta/2)^4)/dtheta =
                        // sin(theta/2)^2 * sin(theta).  Project through the
                        // stable atan2 angle Jacobian instead of a rounded
                        // float cosine, preserving tiny nonzero forces.
                        const float dE_dtheta =
                            dE_dsinp4 * half_sine_squared * geometry.sine;

                        VECTOR fd = {0.0f, 0.0f, 0.0f};
                        VECTOR fa = {0.0f, 0.0f, 0.0f};
                        if (!geometry.is_collinear)
                        {
                            fd = {static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_du[0]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_du[1]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_du[2])};
                            fa = {static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_dv[0]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_dv[1]),
                                  static_cast<float>(
                                      static_cast<double>(dE_dtheta) *
                                      geometry.dtheta_dv[2])};
                        }

                        VECTOR fh;
                        fh.x = -(fd.x + fa.x);
                        fh.y = -(fd.y + fa.y);
                        fh.z = -(fd.z + fa.z);

                        const VECTOR force_d = {-fd.x, -fd.y, -fd.z};
                        const VECTOR force_a = {f_a_rad.x - fa.x,
                                                f_a_rad.y - fa.y,
                                                f_a_rad.z - fa.z};
                        const VECTOR force_h = {-f_a_rad.x - fh.x,
                                                -f_a_rad.y - fh.y,
                                                -f_a_rad.z - fh.z};
                        LTMatrix3 pair_virial = {0, 0, 0, 0, 0, 0};
                        if (atom_virial)
                        {
                            pair_virial =
                                Get_Virial_From_Force_Dis(force_d, ddh) +
                                Get_Virial_From_Force_Dis(force_a, dah);
                        }
                        const double next_en_hb =
                            en_hb + static_cast<double>(s_en_total.val);
                        const float next_en_hb_float =
                            static_cast<float>(next_en_hb);

                        if (!ReaxFF_Float_Is_Finite(s_en_total.val) ||
                            !ReaxFF_Float_Is_Finite(s_en_total.dval[0]) ||
                            !ReaxFF_Double_Is_Finite(radial_shape) ||
                            !ReaxFF_Double_Is_Finite(exp_hb3_double) ||
                            !ReaxFF_Double_Is_Finite(radial_log_derivative) ||
                            !ReaxFF_Double_Is_Finite(dE_dr_ah_double) ||
                            !ReaxFF_Float_Is_Finite(dE_dr_ah) ||
                            !ReaxFF_Float_Is_Finite(dE_dtheta) ||
                            !ReaxFF_Vector_Is_Finite(f_a_rad) ||
                            !ReaxFF_Vector_Is_Finite(fd) ||
                            !ReaxFF_Vector_Is_Finite(fa) ||
                            !ReaxFF_Vector_Is_Finite(fh) ||
                            !ReaxFF_Vector_Is_Finite(force_d) ||
                            !ReaxFF_Vector_Is_Finite(force_a) ||
                            !ReaxFF_Vector_Is_Finite(force_h) ||
                            (atom_virial &&
                             !ReaxFF_Matrix_Is_Finite(pair_virial)) ||
                            !ReaxFF_Double_Is_Finite(next_en_hb) ||
                            !ReaxFF_Float_Is_Finite(next_en_hb_float))
                        {
                            Record_ReaxFF_Geometry_Error(
                                geometry_error, REAXFF_HB_NONFINITE, d, h, a);
                            continue;
                        }

                        atomicAdd(&d_dE_dBO_s[b_dh], s_en_total.dval[0]);
                        atomicAdd(&d_dE_dBO_pi[b_dh], s_en_total.dval[0]);
                        atomicAdd(&d_dE_dBO_pi2[b_dh], s_en_total.dval[0]);
                        atomicAdd(&frc[d].x, force_d.x);
                        atomicAdd(&frc[d].y, force_d.y);
                        atomicAdd(&frc[d].z, force_d.z);
                        atomicAdd(&frc[a].x, force_a.x);
                        atomicAdd(&frc[a].y, force_a.y);
                        atomicAdd(&frc[a].z, force_a.z);
                        atomicAdd(&frc[h].x, force_h.x);
                        atomicAdd(&frc[h].y, force_h.y);
                        atomicAdd(&frc[h].z, force_h.z);
                        if (atom_virial)
                        {
                            atomicAdd(atom_virial + h, pair_virial);
                        }

                        en_hb = next_en_hb;
                    }
                }
            }
            atomicAdd(d_energy_hb_sum, (float)en_hb);
            if (atom_energy) atomicAdd(&atom_energy[h], (float)en_hb);
        }
    }
}

void REAXFF_HYDROGEN_BOND::Initial(CONTROLLER* controller, int atom_numbers,
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
            "REAXFF_HYDROGEN_BOND IS NOT INITIALIZED (missing input "
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
            "REAXFF_HYDROGEN_BOND::Initial", "Reason:\n\t%s",
            reason.c_str());
        return;
    }
    std::vector<int> atom_type;
    std::vector<int> is_hydrogen;
    if (!ReaxFF_Parse_Type_File_Path(type_in_file, atom_numbers, force_field,
                                     &atom_type, &is_hydrogen, &input_error))
    {
        const std::string reason = input_error.Describe();
        controller->Throw_Formatted_SPONGE_Error(
            input_error.kind == REAXFF_INPUT_OPEN_ERROR
                ? spongeErrorOpenFileFailed
                : spongeErrorBadFileFormat,
            "REAXFF_HYDROGEN_BOND::Initial", "Reason:\n\t%s",
            reason.c_str());
        return;
    }

    const int n_atom_types = static_cast<int>(force_field.atom_types.size());
    int triplet_parameter_count = 0;
    if (!ReaxFF_Checked_Dense_Table_Count(
            n_atom_types, 3, 3 * sizeof(REAXFF_HB_Info),
            &triplet_parameter_count))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_HYDROGEN_BOND::Initial",
            "Reason:\n\tatom type count %d exceeds the supported hydrogen-bond "
            "parameter table extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    std::map<int, std::vector<REAXFF_HB_Entry>> triplet_entries;
    for (std::size_t i = 0; i < force_field.hydrogen_bonds.size(); i++)
    {
        const REAXFF_HYDROGEN_BOND_IR& source =
            force_field.hydrogen_bonds[i];
        if (source.values[1] != 0.0f && !(source.values[0] > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "REAXFF_HYDROGEN_BOND::Initial",
                "Reason:\n\tactive hydrogen bond entry %zu must have positive "
                "r0 in file %s",
                i + 1, parameter_in_file);
            return;
        }
        REAXFF_HB_Entry entry = {
            source.values[0], source.values[1], source.values[2],
            source.values[3]};
        const int idx1 = source.types[0] - 1;
        const int idx2 = source.types[1] - 1;
        const int idx3 = source.types[2] - 1;
        const int triplet =
            (idx1 * n_atom_types + idx2) * n_atom_types + idx3;
        triplet_entries[triplet].push_back(entry);
    }

    std::vector<REAXFF_HB_Info> hb_info(
        triplet_parameter_count, REAXFF_HB_Info{0, 0});
    std::vector<REAXFF_HB_Entry> sorted_entries;
    for (const auto& indexed_entries : triplet_entries)
    {
        if (indexed_entries.first < 0 ||
            indexed_entries.first >= triplet_parameter_count ||
            sorted_entries.size() >
                static_cast<std::size_t>(INT_MAX) -
                    indexed_entries.second.size())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "REAXFF_HYDROGEN_BOND::Initial",
                "Reason:\n\thydrogen-bond parameter index/count exceeds the "
                "supported representation in file %s",
                parameter_in_file);
            return;
        }
        REAXFF_HB_Info& info = hb_info[indexed_entries.first];
        info.start_idx = static_cast<int>(sorted_entries.size());
        info.entry_count =
            static_cast<int>(indexed_entries.second.size());
        sorted_entries.insert(sorted_entries.end(),
                              indexed_entries.second.begin(),
                              indexed_entries.second.end());
    }

    REAXFF_HB_Info* staged_hb_info = NULL;
    REAXFF_HB_Entry* staged_hb_entries = NULL;
    int* staged_atom_type = NULL;
    int* staged_is_hydrogen = NULL;
    Malloc_Safely((void**)&staged_hb_info,
                  sizeof(REAXFF_HB_Info) * hb_info.size());
    if (!sorted_entries.empty())
        Malloc_Safely((void**)&staged_hb_entries,
                      sizeof(REAXFF_HB_Entry) * sorted_entries.size());
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    Malloc_Safely((void**)&staged_is_hydrogen, sizeof(int) * atom_numbers);
    memcpy(staged_hb_info, hb_info.data(),
           sizeof(REAXFF_HB_Info) * hb_info.size());
    if (!sorted_entries.empty())
        memcpy(staged_hb_entries, sorted_entries.data(),
               sizeof(REAXFF_HB_Entry) * sorted_entries.size());
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);
    memcpy(staged_is_hydrogen, is_hydrogen.data(),
           sizeof(int) * atom_numbers);

    this->controller = controller;
    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    h_hb_info = staged_hb_info;
    h_hb_entries = staged_hb_entries;
    h_atom_type = staged_atom_type;
    h_is_hydrogen = staged_is_hydrogen;

    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_is_hydrogen_global, h_is_hydrogen,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_is_hydrogen, sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely(
        (void**)&d_hb_info, h_hb_info,
        sizeof(REAXFF_HB_Info) * triplet_parameter_count);
    if (!sorted_entries.empty())
        Device_Malloc_And_Copy_Safely(
            (void**)&d_hb_entries, h_hb_entries,
            sizeof(REAXFF_HB_Entry) * sorted_entries.size());
    Device_Malloc_Safely((void**)&d_energy_hb_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_geometry_error,
                         sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);

    controller->Step_Print_Initial("REAXFF_HB", "%14.7e");
    is_initialized = 1;
}
void REAXFF_HYDROGEN_BOND::Calculate_HB_Energy_And_Force(
    int atom_numbers, const VECTOR* crd, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, REAXFF_BOND_ORDER* bo_module,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_virial)
{
    if (!is_initialized) return;
    dim3 blockSize(32);
    dim3 gridSize((atom_numbers + blockSize.x - 1) / blockSize.x);
    deviceMemset(d_energy_hb_sum, 0, sizeof(float));
    deviceMemset(d_geometry_error, 0, sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);

    Launch_Device_Kernel(
        Calculate_HB_Kernel, gridSize, blockSize, 0, NULL, atom_numbers, crd,
        d_atom_type, d_is_hydrogen, d_hb_info, d_hb_entries, atom_type_numbers,
        bo_module->d_corrected_bo_s, bo_module->d_corrected_bo_pi,
        bo_module->d_corrected_bo_pi2, d_dE_dBO_s, d_dE_dBO_pi, d_dE_dBO_pi2,
        cell, rcell, nl, atom_energy, frc, need_virial ? atom_virial : NULL,
        d_energy_hb_sum, bo_module->d_bond_count, bo_module->d_bond_offset,
        bo_module->d_bond_nbr, bo_module->d_bond_idx, d_geometry_error);
    Check_Geometry_Error();
}

void REAXFF_HYDROGEN_BOND::Check_Geometry_Error()
{
    int error[REAXFF_GEOMETRY_ERROR_SIZE] = {0, -1, -1, -1, -1};
    deviceMemcpy(error, d_geometry_error, sizeof(error),
                 deviceMemcpyDeviceToHost);
    if (error[0] == REAXFF_GEOMETRY_OK) return;
    if (error[0] == REAXFF_INVALID_ATOM_TYPE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF_HYDROGEN_BOND::Calculate_HB_Energy_And_Force",
            "Invalid ReaxFF atom type in hydrogen-bond evaluation record "
            "[%d, %d, %d, %d] (valid range is [0, %d)).",
            error[1], error[2], error[3], error[4], atom_type_numbers);
        return;
    }
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown,
        "REAXFF_HYDROGEN_BOND::Calculate_HB_Energy_And_Force",
        "Reason:\n\tlocal donor/hydrogen/acceptor atoms %d %d %d have an "
        "undefined zero-arm or non-finite/unrepresentable ReaxFF hydrogen-"
        "bond geometry, energy, derivative, force, or virial\n",
        error[1], error[2], error[3]);
}

void REAXFF_HYDROGEN_BOND::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    deviceMemcpy(&h_energy_hb, d_energy_hb_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print("REAXFF_HB", h_energy_hb, true);
}
