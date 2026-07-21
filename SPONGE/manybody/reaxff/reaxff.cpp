#include "reaxff.h"

#include "reaxff_evaluation.h"

void REAXFF::Initial(CONTROLLER* controller, int atom_numbers, float cutoff,
                     float* cutoff_full, bool* need_full_nl_flag)
{
    is_initialized = 0;
    if (!controller->Command_Exist("REAXFF", "in_file"))
    {
        return;
    }

    this->controller = controller;
    this->atom_numbers = atom_numbers;

    const char* parameter_in_file =
        controller->Original_Command("REAXFF", "in_file");
    const char* type_in_file =
        controller->Original_Command("REAXFF", "type_in_file");

    eeq.Initial(controller, atom_numbers, parameter_in_file, type_in_file);
    bond_order.Initial(controller, atom_numbers, parameter_in_file,
                       type_in_file, cutoff, cutoff_full);
    bond.Initial(controller, atom_numbers, "REAXFF", need_full_nl_flag);
    vdw.Initial(controller, atom_numbers, "REAXFF", need_full_nl_flag);
    ovun.Initial(controller, atom_numbers, "REAXFF");
    angle.Initial(controller, atom_numbers, "REAXFF");
    torsion.Initial(controller, atom_numbers, "REAXFF");
    hb.Initial(controller, atom_numbers, "REAXFF");

    Validate_Module_Consistency();
    Wire_Shared_State();

    if (atom_numbers <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Initial",
            "ReaxFF requires a positive atom count, but received %d.",
            atom_numbers);
    }
    Ensure_Staging_Capacity(atom_numbers);
    Device_Malloc_Safely((void**)&d_seen_global, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_evaluation_error,
                         sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);
    is_initialized = 1;
}

void REAXFF::Ensure_Staging_Capacity(int coordinate_count)
{
    if (coordinate_count <= staging_capacity) return;
    if (coordinate_count <= 0 ||
        static_cast<std::size_t>(coordinate_count) >
            std::numeric_limits<std::size_t>::max() / sizeof(LTMatrix3))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "REAXFF::Ensure_Staging_Capacity",
            "ReaxFF coordinate capacity %d cannot be represented safely.",
            coordinate_count);
    }

    auto release = [](void** pointer)
    {
        if (pointer != NULL && pointer[0] != NULL)
            Free_Single_Device_Pointer(pointer);
    };
    release((void**)&d_staged_frc);
    release((void**)&d_staged_energy);
    release((void**)&d_staged_virial);
    release((void**)&d_staged_charge);
    Device_Malloc_Safely((void**)&d_staged_frc,
                         sizeof(VECTOR) * coordinate_count);
    Device_Malloc_Safely((void**)&d_staged_energy,
                         sizeof(float) * coordinate_count);
    Device_Malloc_Safely((void**)&d_staged_virial,
                         sizeof(LTMatrix3) * coordinate_count);
    Device_Malloc_Safely((void**)&d_staged_charge,
                         sizeof(float) * coordinate_count);
    staging_capacity = coordinate_count;
}

void REAXFF::Validate_Module_Consistency()
{
    struct MODULE_VIEW
    {
        const char* name;
        int initialized;
        int atoms;
        int types;
        const int* atom_type;
    };
    const MODULE_VIEW modules[] = {
        {"EEQ", eeq.is_initialized, eeq.atom_numbers, eeq.atom_type_numbers,
         eeq.h_atom_type},
        {"bond order", bond_order.is_initialized, bond_order.atom_numbers,
         bond_order.atom_type_numbers, bond_order.h_atom_type},
        {"bond", bond.is_initialized, bond.atom_numbers,
         bond.atom_type_numbers, bond.h_atom_type},
        {"van der Waals", vdw.is_initialized, vdw.atom_numbers,
         vdw.atom_type_numbers, vdw.h_atom_type},
        {"over/under coordination", ovun.is_initialized, ovun.atom_numbers,
         ovun.atom_type_numbers, ovun.h_atom_type},
        {"valence angle", angle.is_initialized, angle.atom_numbers,
         angle.atom_type_numbers, angle.h_atom_type},
        {"torsion", torsion.is_initialized, torsion.atom_numbers,
         torsion.atom_type_numbers, torsion.h_atom_type},
        {"hydrogen bond", hb.is_initialized, hb.atom_numbers,
         hb.atom_type_numbers, hb.h_atom_type}};

    for (const MODULE_VIEW& module : modules)
    {
        if (!module.initialized || module.atoms != atom_numbers ||
            module.types != eeq.atom_type_numbers ||
            (atom_numbers > 0 && module.atom_type == NULL))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "REAXFF::Validate_Module_Consistency",
                "ReaxFF %s state is inconsistent: initialized=%d, atom "
                "count=%d/%d, type count=%d/%d, type table=%p.",
                module.name, module.initialized, module.atoms, atom_numbers,
                module.types, eeq.atom_type_numbers,
                static_cast<const void*>(module.atom_type));
        }
        for (int atom = 0; atom < atom_numbers; atom++)
        {
            if (module.atom_type[atom] != eeq.h_atom_type[atom])
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "REAXFF::Validate_Module_Consistency",
                    "ReaxFF %s maps global atom %d to type %d, while the "
                    "authoritative EEQ table maps it to type %d.",
                    module.name, atom, module.atom_type[atom],
                    eeq.h_atom_type[atom]);
            }
        }
    }
    if (hb.h_is_hydrogen == NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "REAXFF::Validate_Module_Consistency",
            "The ReaxFF hydrogen identity table is NULL.");
    }
    for (int atom = 0; atom < atom_numbers; atom++)
    {
        if (hb.h_is_hydrogen[atom] != 0 && hb.h_is_hydrogen[atom] != 1)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "REAXFF::Validate_Module_Consistency",
                "Global atom %d has invalid ReaxFF hydrogen flag %d.", atom,
                hb.h_is_hydrogen[atom]);
        }
    }
}

void REAXFF::Wire_Shared_State()
{
    bond.d_bo_s = bond_order.d_corrected_bo_s;
    bond.d_bo_pi = bond_order.d_corrected_bo_pi;
    bond.d_bo_pi2 = bond_order.d_corrected_bo_pi2;
    bond.d_dE_dBO_s = bond_order.d_dE_dBO_s;
    bond.d_dE_dBO_pi = bond_order.d_dE_dBO_pi;
    bond.d_dE_dBO_pi2 = bond_order.d_dE_dBO_pi2;
    bond.d_dbo_s_dr = bond_order.d_dbo_s_dr;
    bond.d_dbo_pi_dr = bond_order.d_dbo_pi_dr;
    bond.d_dbo_pi2_dr = bond_order.d_dbo_pi2_dr;
    bond.d_dbo_s_dDelta_i = bond_order.d_dbo_s_dDelta_i;
    bond.d_dbo_pi_dDelta_i = bond_order.d_dbo_pi_dDelta_i;
    bond.d_dbo_pi2_dDelta_i = bond_order.d_dbo_pi2_dDelta_i;
    bond.d_dbo_s_dDelta_j = bond_order.d_dbo_s_dDelta_j;
    bond.d_dbo_pi_dDelta_j = bond_order.d_dbo_pi_dDelta_j;
    bond.d_dbo_pi2_dDelta_j = bond_order.d_dbo_pi2_dDelta_j;
    bond.d_dbo_raw_total_dr = bond_order.d_dbo_raw_total_dr;
    bond.d_bond_count = bond_order.d_bond_count;
    bond.d_bond_offset = bond_order.d_bond_offset;
    bond.d_bond_nbr = bond_order.d_bond_nbr;
    bond.d_bond_idx = bond_order.d_bond_idx;

    ovun.d_dE_dBO_s = bond_order.d_dE_dBO_s;
    ovun.d_dE_dBO_pi = bond_order.d_dE_dBO_pi;
    ovun.d_dE_dBO_pi2 = bond_order.d_dE_dBO_pi2;
    ovun.d_dbo_s_dr = bond_order.d_dbo_s_dr;
    ovun.d_dbo_pi_dr = bond_order.d_dbo_pi_dr;
    ovun.d_dbo_pi2_dr = bond_order.d_dbo_pi2_dr;
    ovun.d_dbo_s_dDelta_i = bond_order.d_dbo_s_dDelta_i;
    ovun.d_dbo_pi_dDelta_i = bond_order.d_dbo_pi_dDelta_i;
    ovun.d_dbo_pi2_dDelta_i = bond_order.d_dbo_pi2_dDelta_i;
    ovun.d_dbo_s_dDelta_j = bond_order.d_dbo_s_dDelta_j;
    ovun.d_dbo_pi_dDelta_j = bond_order.d_dbo_pi_dDelta_j;
    ovun.d_dbo_pi2_dDelta_j = bond_order.d_dbo_pi2_dDelta_j;
    ovun.d_dbo_raw_total_dr = bond_order.d_dbo_raw_total_dr;

    bond.d_CdDelta = ovun.d_CdDelta;

    angle.d_dE_dBO_s = bond_order.d_dE_dBO_s;
    angle.d_dE_dBO_pi = bond_order.d_dE_dBO_pi;
    angle.d_dE_dBO_pi2 = bond_order.d_dE_dBO_pi2;
    angle.d_CdDelta = ovun.d_CdDelta;
    angle.d_dbo_s_dr = bond_order.d_dbo_s_dr;
    angle.d_dbo_pi_dr = bond_order.d_dbo_pi_dr;
    angle.d_dbo_pi2_dr = bond_order.d_dbo_pi2_dr;
    angle.d_dbo_s_dDelta_i = bond_order.d_dbo_s_dDelta_i;
    angle.d_dbo_pi_dDelta_i = bond_order.d_dbo_pi_dDelta_i;
    angle.d_dbo_pi2_dDelta_i = bond_order.d_dbo_pi2_dDelta_i;
    angle.d_dbo_s_dDelta_j = bond_order.d_dbo_s_dDelta_j;
    angle.d_dbo_pi_dDelta_j = bond_order.d_dbo_pi_dDelta_j;
    angle.d_dbo_pi2_dDelta_j = bond_order.d_dbo_pi2_dDelta_j;
    angle.d_dbo_raw_total_dr = bond_order.d_dbo_raw_total_dr;

    torsion.d_dE_dBO_s = bond_order.d_dE_dBO_s;
    torsion.d_dE_dBO_pi = bond_order.d_dE_dBO_pi;
    torsion.d_dE_dBO_pi2 = bond_order.d_dE_dBO_pi2;
    torsion.d_CdDelta = ovun.d_CdDelta;

    hb.d_dE_dBO_s = bond_order.d_dE_dBO_s;
    hb.d_dE_dBO_pi = bond_order.d_dE_dBO_pi;
    hb.d_dE_dBO_pi2 = bond_order.d_dE_dBO_pi2;
}

void REAXFF::Step_Print(CONTROLLER* controller, const float* d_charge)
{
    bond.Step_Print(controller);
    vdw.Step_Print(controller);
    eeq.Step_Print(controller);
    if (eeq.is_initialized)
    {
        eeq.Print_Charges(d_charge);
    }
    ovun.Step_Print_ELP(controller);
    ovun.Step_Print(controller);
    angle.Step_Print(controller);
    torsion.Step_Print(controller);
    hb.Step_Print(controller);

    if (bond.is_initialized && vdw.is_initialized && eeq.is_initialized)
    {
        const float total_reaxff =
            bond.h_energy_sum + vdw.h_energy_sum + eeq.h_energy +
            ovun.h_energy_lp + ovun.h_energy_ovun + angle.h_energy_ang +
            angle.h_energy_pen + angle.h_energy_coa + torsion.h_energy_tor +
            torsion.h_energy_cot + hb.h_energy_hb;
        controller->Step_Print("REAXFF", total_reaxff);
    }
}

static bool ReaxFF_Host_Matrix_Is_Finite(const LTMatrix3& matrix)
{
    return Float_Memory_Is_Finite(&matrix.a11) &&
           Float_Memory_Is_Finite(&matrix.a21) &&
           Float_Memory_Is_Finite(&matrix.a22) &&
           Float_Memory_Is_Finite(&matrix.a31) &&
           Float_Memory_Is_Finite(&matrix.a32) &&
           Float_Memory_Is_Finite(&matrix.a33);
}

void REAXFF::Preflight_Evaluation(DOMAIN_INFORMATION* dd,
                                  MD_INFORMATION* md_info,
                                  NEIGHBOR_LIST* neighbor_list)
{
    if (dd == NULL || md_info == NULL || neighbor_list == NULL)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Preflight_Evaluation",
            "ReaxFF received null evaluation state (DD=%p, MD=%p, NL=%p).",
            static_cast<void*>(dd), static_cast<void*>(md_info),
            static_cast<void*>(neighbor_list));
    }
    int coordinate_count = 0;
    const ReaxFFEvaluation::COORDINATE_COUNT_STATUS count_status =
        ReaxFFEvaluation::Checked_Coordinate_Count(
            dd->atom_numbers, dd->ghost_numbers, sizeof(LTMatrix3),
            &coordinate_count);
    if (count_status != ReaxFFEvaluation::COORDINATE_COUNT_OK)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "REAXFF::Preflight_Evaluation",
            "Invalid ReaxFF owned/ghost coordinate counts %d/%d: checked "
            "owned+ghost staging extent failed with status %d (negative count, "
            "int addition overflow, or size_t byte-count overflow).",
            dd->atom_numbers, dd->ghost_numbers,
            static_cast<int>(count_status));
    }
    if (CONTROLLER::PP_MPI_size != 1 || dd->atom_numbers != atom_numbers ||
        dd->ghost_numbers != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Preflight_Evaluation",
            "ReaxFF requires one PP rank, exactly %d owned atoms, and no "
            "ghosts; received PP=%d, owned=%d, ghosts=%d (checked coordinate "
            "extent=%d). Distributed/ghost evaluation is rejected before any "
            "layout-driven staging allocation and before the first Hamiltonian "
            "kernel because EEQ "
            "charges and bond-order coordination do not yet implement the "
            "required cross-rank global interaction semantics; evaluating only "
            "the local view would change the Hamiltonian.",
            atom_numbers, CONTROLLER::PP_MPI_size, dd->atom_numbers,
            dd->ghost_numbers, coordinate_count);
    }
    Ensure_Staging_Capacity(coordinate_count);
    if (dd->crd == NULL || dd->frc == NULL || dd->d_energy == NULL ||
        dd->d_virial == NULL || dd->d_charge == NULL ||
        dd->atom_local == NULL || md_info->d_charge == NULL ||
        d_staged_frc == NULL || d_staged_energy == NULL ||
        d_staged_virial == NULL || d_staged_charge == NULL ||
        d_seen_global == NULL || d_evaluation_error == NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Preflight_Evaluation",
            "ReaxFF evaluation or staging storage is NULL.");
    }

    const float cutoff = md_info->nb.cutoff;
    const LTMatrix3 cell = md_info->pbc.cell;
    const LTMatrix3 rcell = md_info->pbc.rcell;
    const double determinant = static_cast<double>(cell.a11) *
                               static_cast<double>(cell.a22) *
                               static_cast<double>(cell.a33);
    if (!(cutoff > 0.0f) || !Float_Memory_Is_Finite(&cutoff) ||
        !ReaxFF_Host_Matrix_Is_Finite(cell) ||
        !ReaxFF_Host_Matrix_Is_Finite(rcell) || !(determinant != 0.0) ||
        !Double_Memory_Is_Finite(&determinant))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Preflight_Evaluation",
            "ReaxFF requires a positive finite cutoff and finite, nonsingular "
            "periodic cell/reciprocal-cell matrices (cutoff=%g, det=%g).",
            cutoff, determinant);
    }

    FULL_NEIGHBOR_LIST& full = neighbor_list->full_neighbor_list;
    if (!neighbor_list->is_initialized ||
        neighbor_list->active_local_atom_numbers != atom_numbers ||
        neighbor_list->max_neighbor_numbers < 0 ||
        neighbor_list->d_nl == NULL || !full.is_initialized ||
        full.active_owned_atom_numbers != atom_numbers ||
        full.max_neighbor_numbers < 0 || full.d_nl == NULL ||
        full.last_build_error != FULL_NEIGHBOR_LIST::BUILD_OK)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Preflight_Evaluation",
            "ReaxFF requires complete validated half/full neighbor lists for "
            "%d atoms (half active/capacity=%d/%d, full active/capacity=%d/%d, "
            "full error=%d).",
            atom_numbers, neighbor_list->active_local_atom_numbers,
            neighbor_list->max_neighbor_numbers, full.active_owned_atom_numbers,
            full.max_neighbor_numbers, full.last_build_error);
    }
    if (eeq.nprev < 0 || eeq.nprev > REAXFF_EEQ::HIST_SIZE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Preflight_Evaluation",
            "EEQ predictor history count %d is outside [0, %d].", eeq.nprev,
            REAXFF_EEQ::HIST_SIZE);
    }

    deviceMemset(d_seen_global, 0, sizeof(int) * atom_numbers);
    deviceMemset(d_evaluation_error, 0,
                 sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);
    const dim3 block_size(128);
    const dim3 grid_size((atom_numbers + block_size.x - 1) / block_size.x);
    Launch_Device_Kernel(
        ReaxFFEvaluation::Preflight_Kernel, grid_size, block_size, 0, NULL,
        atom_numbers, dd->crd, dd->atom_local, eeq.d_atom_type,
        eeq.d_atom_type_global, eeq.atom_type_numbers, bond_order.d_atom_type,
        bond.d_atom_type, vdw.d_atom_type, ovun.d_atom_type, angle.d_atom_type,
        torsion.d_atom_type, hb.d_atom_type, hb.d_is_hydrogen,
        hb.d_is_hydrogen_global, neighbor_list->d_nl,
        neighbor_list->max_neighbor_numbers, full.d_nl,
        full.max_neighbor_numbers, d_seen_global, d_evaluation_error);
    Check_Evaluation_Error("REAXFF::Preflight_Evaluation");
}

void REAXFF::Check_Evaluation_Error(const char* error_by)
{
    int error[REAXFF_GEOMETRY_ERROR_SIZE] = {0, -1, -1, -1, -1};
    deviceMemcpy(error, d_evaluation_error, sizeof(error),
                 deviceMemcpyDeviceToHost);
    if (error[0] == ReaxFFEvaluation::EVALUATION_OK) return;

    const char* reason = "an unknown validation failure";
    switch (error[0])
    {
        case ReaxFFEvaluation::INVALID_LOCAL_TO_GLOBAL:
            reason = "an out-of-range local-to-global atom ID";
            break;
        case ReaxFFEvaluation::DUPLICATE_GLOBAL_ATOM:
            reason = "a duplicate local-to-global atom ID";
            break;
        case ReaxFFEvaluation::NONFINITE_COORDINATE:
            reason = "a non-finite coordinate";
            break;
        case ReaxFFEvaluation::INVALID_OR_INCONSISTENT_TYPE:
            reason = "an invalid or module-inconsistent atom type";
            break;
        case ReaxFFEvaluation::INVALID_HYDROGEN_FLAG:
            reason = "an invalid or stale hydrogen identity flag";
            break;
        case ReaxFFEvaluation::INVALID_NEIGHBOR_ROW:
            reason = "an invalid neighbor-list row";
            break;
        case ReaxFFEvaluation::INVALID_NEIGHBOR_ATOM:
            reason = "an out-of-range/self neighbor-list entry";
            break;
        case ReaxFFEvaluation::NONFINITE_STAGED_RESULT:
            reason = "a non-finite staged force, energy, virial, or charge";
            break;
        case ReaxFFEvaluation::NONFINITE_COMMITTED_RESULT:
            reason = "a non-finite value after combining staged and existing "
                     "force/energy/virial";
            break;
        case ReaxFFEvaluation::NONFINITE_HISTORY:
            reason = "a non-finite EEQ predictor value/history frame";
            break;
    }
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown, error_by,
        "Transactional ReaxFF validation found %s at record [%d, %d, %d, %d]; "
        "no staged ReaxFF force, energy, virial, charge, or history was "
        "published.",
        reason, error[1], error[2], error[3], error[4]);
}

void REAXFF::Validate_Staged_Result(DOMAIN_INFORMATION* dd,
                                    MD_INFORMATION* md_info,
                                    bool commit_sampling_state)
{
    if (!Float_Memory_Is_Finite(&eeq.pending_energy))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Validate_Staged_Result",
            "EEQ produced a non-finite staged total energy; no ReaxFF state "
            "was published.");
    }
    struct COMPONENT_ENERGY_VIEW
    {
        const char* name;
        const float* device_value;
    };
    const COMPONENT_ENERGY_VIEW component_energies[] = {
        {"bond", bond.d_energy_sum},
        {"van der Waals", vdw.d_energy_sum},
        {"lone-pair", ovun.d_energy_elp_sum},
        {"over/under-coordination", ovun.d_energy_ovun_sum},
        {"valence-angle", angle.d_energy_ang_sum},
        {"penalty", angle.d_energy_pen_sum},
        {"three-body conjugation", angle.d_energy_coa_sum},
        {"torsion", torsion.d_energy_tor_sum},
        {"four-body conjugation", torsion.d_energy_cot_sum},
        {"hydrogen-bond", hb.d_energy_hb_sum}};
    for (const COMPONENT_ENERGY_VIEW& component : component_energies)
    {
        if (component.device_value == NULL)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "REAXFF::Validate_Staged_Result",
                "ReaxFF %s staged component-energy storage is NULL; no "
                "ReaxFF state was published.",
                component.name);
        }
        float value = 0.0f;
        deviceMemcpy(&value, component.device_value, sizeof(value),
                     deviceMemcpyDeviceToHost);
        if (!Float_Memory_Is_Finite(&value))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "REAXFF::Validate_Staged_Result",
                "ReaxFF %s staged component energy is non-finite; no force, "
                "energy, virial, charge, or EEQ history was published.",
                component.name);
        }
    }
    deviceMemset(d_evaluation_error, 0,
                 sizeof(int) * REAXFF_GEOMETRY_ERROR_SIZE);
    const dim3 block_size(128);
    const dim3 grid_size((atom_numbers + block_size.x - 1) / block_size.x);
    Launch_Device_Kernel(
        ReaxFFEvaluation::Validate_Staging_Kernel, grid_size, block_size, 0,
        NULL, atom_numbers, d_staged_frc, d_staged_energy, d_staged_virial,
        d_staged_charge, dd->frc, dd->d_energy, dd->d_virial,
        md_info->need_potential != 0, md_info->need_pressure != 0, eeq.d_s,
        eeq.d_t, eeq.d_s_hist, eeq.d_t_hist, eeq.nprev,
        commit_sampling_state, d_evaluation_error);
    Check_Evaluation_Error("REAXFF::Validate_Staged_Result");
}

void REAXFF::Commit_Staged_Result(DOMAIN_INFORMATION* dd,
                                  MD_INFORMATION* md_info,
                                  bool commit_sampling_state)
{
    const int history_count = eeq.nprev;
    const dim3 block_size(128);
    const dim3 grid_size((atom_numbers + block_size.x - 1) / block_size.x);
    Launch_Device_Kernel(
        ReaxFFEvaluation::Commit_Kernel, grid_size, block_size, 0, NULL,
        atom_numbers, dd->atom_local, d_staged_frc, d_staged_energy,
        d_staged_virial, d_staged_charge, dd->frc, dd->d_energy, dd->d_virial,
        dd->d_charge, md_info->d_charge, md_info->need_potential != 0,
        md_info->need_pressure != 0, eeq.d_s, eeq.d_t, eeq.d_s_hist,
        eeq.d_t_hist, history_count, REAXFF_EEQ::HIST_SIZE,
        commit_sampling_state);

    // The copy is intentionally after the publication kernel: it synchronizes
    // the default stream before host-side history/energy metadata is advanced.
    int synchronized_error = 0;
    deviceMemcpy(&synchronized_error, d_evaluation_error, sizeof(int),
                 deviceMemcpyDeviceToHost);
#ifdef USE_GPU
    const deviceError_t device_error = deviceGetLastError();
    if (device_error != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Commit_Staged_Result",
            "Device failure in the proven ReaxFF publication kernel: %s.",
            deviceGetErrorString(device_error));
    }
#endif
    if (commit_sampling_state && eeq.nprev < REAXFF_EEQ::HIST_SIZE)
        eeq.nprev++;
    eeq.h_energy = eeq.pending_energy;
}

void REAXFF::Calculate_Force(DOMAIN_INFORMATION* dd, MD_INFORMATION* md_info,
                             NEIGHBOR_LIST* neighbor_list,
                             bool commit_sampling_state)
{
    if (!is_initialized) return;
    Preflight_Evaluation(dd, md_info, neighbor_list);

    deviceMemset(d_staged_frc, 0, sizeof(VECTOR) * dd->atom_numbers);
    deviceMemset(d_staged_energy, 0, sizeof(float) * dd->atom_numbers);
    deviceMemset(d_staged_virial, 0, sizeof(LTMatrix3) * dd->atom_numbers);
    deviceMemset(d_staged_charge, 0, sizeof(float) * dd->atom_numbers);

    bond_order.Calculate_Bond_Order(
        dd->atom_numbers, dd->crd, md_info->pbc.cell, md_info->pbc.rcell,
        neighbor_list->full_neighbor_list.d_nl, md_info->nb.cutoff);
    // Bond-order storage grows on demand.  Refresh every shared alias before
    // any consumer can observe pointers released by an expansion.
    Wire_Shared_State();

    if (bond_order.is_initialized)
    {
        bond_order.Clear_Derivatives(dd->atom_numbers, ovun.d_CdDelta);
    }

    eeq.Calculate_Charges(
        dd->atom_numbers, dd->atom_local, d_staged_charge, dd->crd,
        md_info->pbc.cell, md_info->pbc.rcell,
        neighbor_list->full_neighbor_list.d_nl, md_info->nb.cutoff,
        d_staged_energy, d_staged_frc, md_info->need_pressure,
        md_info->need_pressure ? d_staged_virial : NULL);

    bond.REAXFF_Bond_Force_With_Atom_Energy_And_Virial(
        dd->atom_numbers, dd->crd, d_staged_frc, md_info->pbc.cell,
        md_info->pbc.rcell, neighbor_list->d_nl, 1, d_staged_energy,
        md_info->need_pressure,
        md_info->need_pressure ? d_staged_virial : NULL);
    vdw.REAXFF_VDW_Force_With_Atom_Energy_And_Virial(
        dd->atom_numbers, dd->crd, d_staged_frc, md_info->pbc.cell,
        md_info->pbc.rcell, neighbor_list->d_nl, md_info->nb.cutoff,
        1, d_staged_energy, md_info->need_pressure,
        md_info->need_pressure ? d_staged_virial : NULL);
    ovun.Calculate_Over_Under_Energy_And_Force(
        dd->atom_numbers, dd->crd, d_staged_frc, md_info->pbc.cell,
        md_info->pbc.rcell, &bond_order, 1, d_staged_energy,
        md_info->need_pressure,
        md_info->need_pressure ? d_staged_virial : NULL);
    angle.Calculate_Valence_Angle_Energy_And_Force(
        dd->atom_numbers, dd->crd, d_staged_frc, md_info->pbc.cell,
        md_info->pbc.rcell, neighbor_list->full_neighbor_list.d_nl, &bond_order,
        ovun.d_Delta, ovun.d_Delta_boc, ovun.d_Delta_val, ovun.d_nlp,
        ovun.d_vlpex, ovun.d_dDelta_lp, ovun.d_CdDelta, 1, d_staged_energy,
        md_info->need_pressure,
        md_info->need_pressure ? d_staged_virial : NULL);
    torsion.Calculate_Torsion_Energy_And_Force(
        dd->atom_numbers, dd->crd, d_staged_frc, md_info->pbc.cell,
        md_info->pbc.rcell, neighbor_list->full_neighbor_list.d_nl, &bond_order,
        ovun.d_Delta_boc, 1, d_staged_energy, md_info->need_pressure,
        md_info->need_pressure ? d_staged_virial : NULL);
    hb.Calculate_HB_Energy_And_Force(
        dd->atom_numbers, dd->crd, d_staged_frc, md_info->pbc.cell,
        md_info->pbc.rcell, neighbor_list->full_neighbor_list.d_nl, &bond_order,
        1, d_staged_energy, md_info->need_pressure,
        md_info->need_pressure ? d_staged_virial : NULL);

    if (bond_order.is_initialized)
    {
        bond_order.Calculate_Forces(dd->atom_numbers, dd->crd, d_staged_frc,
                                    md_info->pbc.cell, md_info->pbc.rcell,
                                    md_info->nb.cutoff, ovun.d_CdDelta,
                                    md_info->need_pressure,
                                    md_info->need_pressure ? d_staged_virial
                                                           : NULL);
    }

    Validate_Staged_Result(dd, md_info, commit_sampling_state);
    Commit_Staged_Result(dd, md_info, commit_sampling_state);
}
