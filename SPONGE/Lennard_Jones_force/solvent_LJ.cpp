#include "solvent_LJ.h"

#include "pair_activity.h"

#ifdef USE_GPU
enum SOLVENT_DISPATCH_ERROR
{
    SOLVENT_DISPATCH_VALID = 0,
    SOLVENT_DISPATCH_BAD_RESIDUE_START = 1,
    SOLVENT_DISPATCH_BAD_NEIGHBOR_LIST = 2,
    SOLVENT_DISPATCH_BAD_NEIGHBOR_INDEX = 3
};

static __device__ __forceinline__ void Record_Solvent_Dispatch_Error(
    int* error, const int code, const int residue, const int entry,
    const int value)
{
    if (atomicCAS(error, SOLVENT_DISPATCH_VALID, code) ==
        SOLVENT_DISPATCH_VALID)
    {
        error[1] = residue;
        error[2] = entry;
        error[3] = value;
    }
}

template <int WAT_POINTS, bool need_force, bool need_energy, bool need_virial,
          bool need_coulomb>
static __global__ void Lennard_Jones_And_Direct_Coulomb_Device(
    const int atom_numbers, const int coordinate_numbers, const ATOM_GROUP* nl,
    const int solvent_start_residue, const int res_numbers,
    const int* d_res_start, const VECTOR_LJ* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const float* LJ_type_A, const float* LJ_type_B,
    const float cutoff, VECTOR* frc, const float pme_beta, float* atom_energy,
    LTMatrix3* atom_lj_virial, float* atom_direct_cf_energy, float* this_energy,
    int* pair_overlap_error, int* dispatch_error)
{
    __shared__ float r1s_x[128];
    __shared__ float r1s_y[128];
    __shared__ float r1s_z[128];
    __shared__ int r1s_lj_type[128];
    __shared__ float r1s_charge[128];
    __shared__ int r1s_global_atom[128];
    int residue_i =
        blockDim.y * blockIdx.x + threadIdx.y + solvent_start_residue;
    if (residue_i < res_numbers)
    {
        VECTOR frc_record[4] = {{0.0f, 0.0f, 0.0f},
                                {0.0f, 0.0f, 0.0f},
                                {0.0f, 0.0f, 0.0f},
                                {0.0f, 0.0f, 0.0f}};
        VECTOR frc_record_j;
        int atom_i = d_res_start[residue_i];
        const int expected_atom_i =
            atom_numbers - (res_numbers - residue_i) * WAT_POINTS;
        if (atom_i < 0 || atom_i > atom_numbers - WAT_POINTS ||
            atom_i != expected_atom_i)
        {
            if (threadIdx.x == 0)
            {
                Record_Solvent_Dispatch_Error(
                    dispatch_error, SOLVENT_DISPATCH_BAD_RESIDUE_START,
                    residue_i, expected_atom_i, atom_i);
            }
            return;
        }
        if (threadIdx.x < WAT_POINTS)
        {
            const int shared_idx = threadIdx.y * WAT_POINTS + threadIdx.x;
            VECTOR_LJ r1 = crd[atom_i + threadIdx.x];
            r1s_x[shared_idx] = r1.crd.x;
            r1s_y[shared_idx] = r1.crd.y;
            r1s_z[shared_idx] = r1.crd.z;
            r1s_lj_type[shared_idx] = r1.LJ_type;
            r1s_charge[shared_idx] = r1.charge;
            r1s_global_atom[shared_idx] = r1.global_atom;
        }
        deviceSyncWarp();
        ATOM_GROUP nl_i = nl[atom_i];
        if (nl_i.atom_numbers < 0 ||
            (nl_i.atom_numbers > 0 && nl_i.atom_serial == NULL))
        {
            if (threadIdx.x == 0)
            {
                Record_Solvent_Dispatch_Error(
                    dispatch_error, SOLVENT_DISPATCH_BAD_NEIGHBOR_LIST,
                    residue_i, nl_i.atom_numbers, 0);
            }
            return;
        }
        LTMatrix3 virial_record = {0, 0, 0, 0, 0, 0};
        float energy_lj = 0.;
        float energy_coulomb = 0.;
        float energy_total = 0.0f;
        for (int j = threadIdx.x; j < nl_i.atom_numbers; j += blockDim.x)
        {
            int atom_j = nl_i.atom_serial[j];
            if (atom_j < 0 || atom_j >= coordinate_numbers)
            {
                Record_Solvent_Dispatch_Error(
                    dispatch_error, SOLVENT_DISPATCH_BAD_NEIGHBOR_INDEX,
                    residue_i, j, atom_j);
                continue;
            }
            float ij_factor = atom_j < atom_numbers ? 1.0f : 0.5f;
            VECTOR_LJ r2 = crd[atom_j];
            frc_record_j = {0.0f, 0.0f, 0.0f};
            for (int i = 0; i < WAT_POINTS; i++)
            {
                const int shared_idx = threadIdx.y * WAT_POINTS + i;
                VECTOR_LJ r1 = {
                    {r1s_x[shared_idx], r1s_y[shared_idx], r1s_z[shared_idx]},
                    r1s_lj_type[shared_idx],
                    r1s_charge[shared_idx],
                    r1s_global_atom[shared_idx]};
                VECTOR dr = Get_Periodic_Displacement(r2, r1, cell, rcell);
                float dr_abs = norm3df(dr.x, dr.y, dr.z);
                if (dr_abs < cutoff)
                {
                    int atom_pair_LJ_type = Get_LJ_Type(r1.LJ_type, r2.LJ_type);
                    float A = LJ_type_A[atom_pair_LJ_type];
                    float B = LJ_type_B[atom_pair_LJ_type];
                    const PairwiseInteraction::Pair_Activity activity =
                        PairwiseInteraction::Classify(
                            A, B,
                            need_coulomb &&
                                PairwiseInteraction::Coulomb_Is_Active(
                                    r1.charge, r2.charge));
                    if (!activity.Any())
                    {
                        continue;
                    }
                    if (dr_abs == 0.0f)
                    {
                        PairwiseInteraction::Fail_Exact_Overlap(
                            r1.global_atom, r2.global_atom,
                            PairwiseInteraction::Components(
                                activity.lennard_jones, activity.coulomb),
                            pair_overlap_error);
                        continue;
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
                            float frc_cf_abs = Get_Direct_Coulomb_Force(
                                r1, r2, dr_abs, pme_beta);
                            frc_abs = frc_abs - frc_cf_abs;
                        }
                        VECTOR frc_lin = frc_abs * dr;
                        frc_record[i] = frc_record[i] + frc_lin;
                        frc_record_j = frc_record_j - frc_lin;
                        if (need_virial)
                        {
                            virial_record =
                                virial_record -
                                ij_factor *
                                    Get_Virial_From_Force_Dis(frc_lin, dr);
                        }
                    }
                    if (need_energy)
                    {
                        if (activity.lennard_jones)
                        {
                            energy_lj +=
                                ij_factor * Get_LJ_Energy(r1, r2, dr_abs, A, B);
                        }
                        if (activity.coulomb)
                        {
                            energy_coulomb +=
                                ij_factor * Get_Direct_Coulomb_Energy(
                                                r1, r2, dr_abs, pme_beta);
                        }
                    }
                }
            }
            if (need_force && atom_j < atom_numbers)
            {
                atomicAdd(frc + atom_j, frc_record_j);
            }
        }
        energy_total = energy_lj + energy_coulomb;
        if (need_force)
        {
            for (int i = 0; i < WAT_POINTS; i++)
            {
                Warp_Sum_To(frc + atom_i + i, frc_record[i], warpSize);
            }
        }
        if (need_energy)
        {
            Warp_Sum_To(atom_energy + atom_i, energy_total, warpSize);
            Warp_Sum_To(this_energy + atom_i, energy_lj, warpSize);
            if (need_coulomb)
                Warp_Sum_To(atom_direct_cf_energy + atom_i, energy_coulomb,
                            warpSize);
        }
        if (need_virial)
        {
            Warp_Sum_To(atom_lj_virial + atom_i, virial_record, warpSize);
        }
    }
}
#endif

static bool Validate_Hard_Solvent_State(
    CONTROLLER* controller, const LENNARD_JONES_INFORMATION* lj_info,
    const MD_INFORMATION* md_info, const char* error_by)
{
    if (lj_info == NULL || md_info == NULL || !lj_info->is_initialized ||
        lj_info->controller == NULL || lj_info->atom_numbers <= 0 ||
        lj_info->atom_numbers != md_info->atom_numbers ||
        lj_info->atom_type_numbers <= 0 || lj_info->h_atom_LJ_type == NULL ||
        lj_info->h_LJ_A == NULL || lj_info->h_LJ_B == NULL ||
        lj_info->d_LJ_A == NULL || lj_info->d_LJ_B == NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\tthe optimized solvent kernel requires one complete "
            "initialized hard-LJ state whose atom count matches MD\n");
        return false;
    }

    const long long expected_pair_types =
        static_cast<long long>(lj_info->atom_type_numbers) *
        (static_cast<long long>(lj_info->atom_type_numbers) + 1LL) / 2LL;
    if (expected_pair_types <= 0 ||
        expected_pair_types > std::numeric_limits<int>::max() ||
        lj_info->pair_type_numbers != static_cast<int>(expected_pair_types))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\thard-LJ type metadata is inconsistent: %d atom "
            "types imply %lld pair types, but %d are stored\n",
            lj_info->atom_type_numbers, expected_pair_types,
            lj_info->pair_type_numbers);
        return false;
    }
    for (int atom = 0; atom < lj_info->atom_numbers; atom++)
    {
        const int type = lj_info->h_atom_LJ_type[atom];
        if (type < 0 || type >= lj_info->atom_type_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, error_by,
                "Reason:\n\thard-LJ atom %d has type %d outside [0, %d)\n",
                atom, type, lj_info->atom_type_numbers);
            return false;
        }
    }
    for (int pair = 0; pair < lj_info->pair_type_numbers; pair++)
    {
        if (!Float_Memory_Is_Finite(lj_info->h_LJ_A + pair) ||
            !Float_Memory_Is_Finite(lj_info->h_LJ_B + pair))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, error_by,
                "Reason:\n\thard-LJ pair type %d has a non-finite "
                "coefficient\n",
                pair);
            return false;
        }
    }
    if (!Float_Memory_Is_Finite(&lj_info->cutoff) || !(lj_info->cutoff > 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\tthe optimized solvent cutoff must be finite and "
            "positive; got %.9g\n",
            lj_info->cutoff);
        return false;
    }
    return true;
}

void SOLVENT_LENNARD_JONES::Initial(CONTROLLER* controller,
                                    LENNARD_JONES_INFORMATION* lj,
                                    LJ_SOFT_CORE* lj_soft,
                                    MD_INFORMATION* md_info,
                                    bool default_enable,
                                    const char* module_name)
{
    const char* error_by = "SOLVENT_LENNARD_JONES::Initial";
    this->lj_info = lj;
    this->lj_soft_info = lj_soft;
    is_initialized = 0;
    solvent_numbers = 0;
    local_solvent_numbers = 0;
    solvent_start = 0;
    solvent_start_local = 0;
    if (module_name == NULL)
    {
        strcpy(this->module_name, "solvent_LJ");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }
    if (controller == NULL || md_info == NULL || lj == NULL || lj_soft == NULL)
    {
        if (controller != NULL)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, error_by,
                "Reason:\n\tthe solvent dispatcher received a null MD or "
                "nonbond state\n");
        }
        return;
    }

    bool enable = false;
#ifdef USE_GPU
    if (!controller->Command_Exist(this->module_name))
    {
        enable = default_enable && (md_info->ug.ug_numbers >= 10);
    }
    else
    {
        enable = controller->Get_Bool(this->module_name,
                                      "SOLVENT_LENNARD_JONES::Initial");
    }

    // This optimized kernel accepts exactly one hard-LJ state.  A soft-core
    // state carries two type maps and two coefficient tables, so converting
    // it to VECTOR_LJ necessarily discards the B Hamiltonian and lambda.
    // Keep the feature fully supported by routing all atoms through the
    // general soft-core kernel instead (local_solvent_numbers stays zero).
    if (enable && lj_soft->is_initialized)
    {
        controller->printf(
            "    soft-core LJ is active; solvent interactions use the "
            "general soft-core dispatch\n");
        enable = false;
    }
    if (enable && !lj->is_initialized)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, error_by,
            "Reason:\n\tthe optimized solvent dispatcher was enabled "
            "without an initialized hard-LJ state\n");
        return;
    }
    if (enable && (md_info->ug.ug_numbers <= 0 || md_info->ug.ug == NULL ||
                   md_info->h_mass == NULL))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\tthe optimized solvent dispatcher requires a "
            "non-empty residue topology and atom masses\n");
        return;
    }
    if (enable &&
        !Validate_Hard_Solvent_State(controller, lj, md_info, error_by))
    {
        return;
    }
    if (enable)
    {
        water_points = md_info->ug.ug[md_info->ug.ug_numbers - 1].atom_numbers;
        if (water_points != 3 && water_points != 4)
        {
            enable = false;
        }
    }
#endif
    if (enable)
    {
        controller->printf("START INITIALIZING SOLVENT LJ:\n");
        solvent_start = md_info->ug.ug_numbers;
        for (int i = md_info->ug.ug_numbers - 1; i >= 0; i -= 1)
        {
            int res_atom_numbers = md_info->ug.ug[i].atom_numbers;
            if (res_atom_numbers != water_points || res_atom_numbers < 3)
            {
                break;
            }
            float mass_O = md_info->h_mass[md_info->ug.ug[i].atom_serial[0]];
            float mass_H1 = md_info->h_mass[md_info->ug.ug[i].atom_serial[1]];
            float mass_H2 = md_info->h_mass[md_info->ug.ug[i].atom_serial[2]];
            if ((mass_O > 15.9f && mass_O < 16.1f) &&
                (mass_H1 > 1.007f && mass_H1 < 1.009f) &&
                (mass_H2 > 1.007f && mass_H2 < 1.009f))
            {
                solvent_numbers += res_atom_numbers;
                solvent_start -= 1;
            }
            else
            {
                break;
            }
        }
        controller->printf("    the solvent is %d-point\n", water_points);
        controller->printf(
            "    the number of solvent atoms is %d (started from Residue "
            "#%d)\n",
            solvent_numbers, solvent_start);
        if (solvent_numbers > 0)
        {
            is_initialized = 1;
            Device_Malloc_Safely((void**)&d_solvent_start_local, sizeof(int));
            Device_Malloc_Safely((void**)&d_local_solvent_numbers, sizeof(int));
            Device_Malloc_Safely((void**)&d_dispatch_error, 4 * sizeof(int));
        }
    }
    if (!is_initialized)
    {
        solvent_numbers = 0;
        local_solvent_numbers = 0;
        controller->printf(
            "    optimized hard-LJ solvent dispatch is inactive; "
            "interactions remain in the general nonbond dispatch\n\n");
    }
    else if (!is_controller_printf_initialized)
    {
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
        controller->printf("END INITIALIZING SOLVENT LJ\n\n");
    }
    else
    {
        controller->printf("END INITIALIZING SOLVENT LJ\n\n");
    }
}

/*
    从输入读入local信息
    atom_numbers: local_原子数
    residue_numbers: local_残基数
    d_res_start: local_残基起始位置
*/
void SOLVENT_LENNARD_JONES::LJ_PME_Direct_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const int residue_numbers, const int* d_res_start,
    const VECTOR* crd, const float* charge, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, const float pme_beta,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_lj_virial, float* atom_direct_pme_energy)
{
    if (!is_initialized) return;

    const char* error_by =
        "SOLVENT_LENNARD_JONES::"
        "LJ_PME_Direct_Force_With_Atom_Energy_And_Virial";
    CONTROLLER* controller = lj_info == NULL ? NULL : lj_info->controller;
    if (controller == NULL || lj_info == NULL || !lj_info->is_initialized ||
        lj_soft_info == NULL || lj_soft_info->is_initialized)
    {
        if (controller != NULL)
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, error_by,
                "Reason:\n\tthe hard-only optimized solvent state changed "
                "after initialization; soft-core states must use the general "
                "dispatch\n");
        return;
    }
    if (!lj_info->Validate_Local_State(error_by, lj_info->atom_numbers,
                                       atom_numbers, lj_info->ghost_numbers, 0))
    {
        return;
    }
    if (atom_numbers < 0 || residue_numbers < 0 || solvent_start_local < 0 ||
        solvent_start_local > residue_numbers || water_points < 3 ||
        water_points > 4 || lj_info->ghost_numbers < 0 ||
        atom_numbers > std::numeric_limits<int>::max() - lj_info->ghost_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s received invalid owned/ghost/residue/solvent "
            "metadata %d/%d/%d/%d for %d-point solvent\n",
            module_name, atom_numbers, lj_info->ghost_numbers, residue_numbers,
            solvent_start_local, water_points);
        return;
    }
    const int selected_residues = residue_numbers - solvent_start_local;
    if (selected_residues == 0) return;
    if (selected_residues > atom_numbers / water_points ||
        local_solvent_numbers != selected_residues * water_points ||
        d_res_start == NULL || crd == NULL || charge == NULL || nl == NULL ||
        frc == NULL || lj_info->crd_with_LJ_parameters_local == NULL ||
        lj_info->d_LJ_A == NULL || lj_info->d_LJ_B == NULL ||
        d_dispatch_error == NULL ||
        (need_atom_energy &&
         (atom_energy == NULL || atom_direct_pme_energy == NULL ||
          lj_info->d_LJ_energy_atom == NULL)) ||
        (need_virial && atom_lj_virial == NULL))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s optimized dispatch has inconsistent selected "
            "residues/atoms (%d/%d) or a required buffer is null\n",
            module_name, selected_residues, local_solvent_numbers);
        return;
    }
#ifdef USE_GPU
    if (CONTROLLER::device_warp <= 0 ||
        CONTROLLER::device_max_thread < CONTROLLER::device_warp)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\tthe device thread geometry cannot form a "
            "solvent kernel warp\n");
        return;
    }
    const int coordinate_numbers = atom_numbers + lj_info->ghost_numbers;
    int* pair_overlap_error = lj_info->d_pair_overlap_error;
    lj_info->Reset_Pair_Overlap_Error();
    deviceMemset(d_dispatch_error, 0, 4 * sizeof(int));
    dim3 blockSize = {
        static_cast<unsigned int>(CONTROLLER::device_warp),
        static_cast<unsigned int>(
            min(8u, static_cast<unsigned int>(CONTROLLER::device_max_thread /
                                              CONTROLLER::device_warp)))};
    dim3 gridSize(static_cast<unsigned int>(
        Positive_Int_Ceil_Div(selected_residues,
                              static_cast<int>(blockSize.y))));

    switch (water_points)
    {
        case 3:
            if (lj_info->is_initialized)
            {
                auto f = Lennard_Jones_And_Direct_Coulomb_Device<3, true, false,
                                                                 false, true>;
                if (!need_atom_energy && !need_virial)
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<3, true, false,
                                                                false, true>;
                }
                else if (need_atom_energy && !need_virial)
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<3, true, true,
                                                                false, true>;
                }
                else if (!need_atom_energy && need_virial)
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<3, true, false,
                                                                true, true>;
                }
                else
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<3, true, true,
                                                                true, true>;
                }
                Launch_Device_Kernel(
                    f, gridSize, blockSize, 0, NULL, atom_numbers,
                    coordinate_numbers, nl, solvent_start_local,
                    residue_numbers, d_res_start,
                    lj_info->crd_with_LJ_parameters_local, cell, rcell,
                    lj_info->d_LJ_A, lj_info->d_LJ_B, lj_info->cutoff, frc,
                    pme_beta, atom_energy, atom_lj_virial,
                    atom_direct_pme_energy, lj_info->d_LJ_energy_atom,
                    pair_overlap_error, d_dispatch_error);
            }
            break;
        case 4:
            if (lj_info->is_initialized)
            {
                auto f = Lennard_Jones_And_Direct_Coulomb_Device<4, true, false,
                                                                 false, true>;
                if (!need_atom_energy && !need_virial)
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<4, true, false,
                                                                false, true>;
                }
                else if (need_atom_energy && !need_virial)
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<4, true, true,
                                                                false, true>;
                }
                else if (!need_atom_energy && need_virial)
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<4, true, false,
                                                                true, true>;
                }
                else
                {
                    f = Lennard_Jones_And_Direct_Coulomb_Device<4, true, true,
                                                                true, true>;
                }
                Launch_Device_Kernel(
                    f, gridSize, blockSize, 0, NULL, atom_numbers,
                    coordinate_numbers, nl, solvent_start_local,
                    residue_numbers, d_res_start,
                    lj_info->crd_with_LJ_parameters_local, cell, rcell,
                    lj_info->d_LJ_A, lj_info->d_LJ_B, lj_info->cutoff, frc,
                    pme_beta, atom_energy, atom_lj_virial,
                    atom_direct_pme_energy, lj_info->d_LJ_energy_atom,
                    pair_overlap_error, d_dispatch_error);
            }
            break;
        default:
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, error_by,
                "Reason:\n\toptimized solvent dispatch received "
                "unsupported water size %d\n",
                water_points);
            return;
    }

    int dispatch_error[4] = {SOLVENT_DISPATCH_VALID, -1, -1, -1};
    deviceMemcpy(dispatch_error, d_dispatch_error, sizeof(dispatch_error),
                 deviceMemcpyDeviceToHost);
    if (dispatch_error[0] != SOLVENT_DISPATCH_VALID)
    {
        const char* reason =
            dispatch_error[0] == SOLVENT_DISPATCH_BAD_RESIDUE_START
                ? "invalid or non-contiguous solvent residue start"
            : dispatch_error[0] == SOLVENT_DISPATCH_BAD_NEIGHBOR_LIST
                ? "invalid solvent neighbor-list record"
                : "neighbor index outside the owned-plus-ghost "
                  "coordinates";
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s at residue %d, entry/expected %d, value %d\n",
            reason, dispatch_error[1], dispatch_error[2], dispatch_error[3]);
        return;
    }
    lj_info->Check_Pair_Overlap_Error(error_by);
#endif
}

static __global__ void get_local_device(const int local_res_numbers,
                                        const int* d_res_len,
                                        const int water_points,
                                        int* d_solvent_start_local,
                                        int* d_local_solvent_numbers,
                                        const int atom_numbers, float* d_mass)
{
#ifdef USE_GPU
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx != 0) return;
#endif
    d_solvent_start_local[0] = local_res_numbers;
    d_local_solvent_numbers[0] = 0;
    int cnt = 0;
    for (int i = local_res_numbers - 1; i >= 0; i -= 1)
    {
        if (d_res_len[i] != water_points || cnt >= atom_numbers / water_points)
        {
            break;
        }
        int res_start = atom_numbers - water_points * (cnt + 1);
        float atom_O = d_mass[res_start];
        float atom_H1 = d_mass[res_start + 1];
        float atom_H2 = d_mass[res_start + 2];
        if ((atom_O > 15.9f && atom_O < 16.1f) &&
            (atom_H1 > 1.007f && atom_H1 < 1.009f) &&
            (atom_H2 > 1.007f && atom_H2 < 1.009f))
        {
            d_solvent_start_local[0] -= 1;
            d_local_solvent_numbers[0] += water_points;
            cnt++;
        }
        else
        {
            break;
        }
    }
}

void SOLVENT_LENNARD_JONES::Get_Local(const int local_res_numbers,
                                      const int* d_res_len,
                                      const int atom_numbers, float* local_mass)
{
    if (!is_initialized) return;

    CONTROLLER* controller = lj_info == NULL ? NULL : lj_info->controller;
    if (controller == NULL || lj_info == NULL || !lj_info->is_initialized ||
        lj_soft_info == NULL || lj_soft_info->is_initialized ||
        local_res_numbers < 0 || atom_numbers < 0 ||
        local_res_numbers > atom_numbers || water_points < 3 ||
        water_points > 4 || d_solvent_start_local == NULL ||
        d_local_solvent_numbers == NULL ||
        (local_res_numbers > 0 && d_res_len == NULL) ||
        (atom_numbers > 0 && local_mass == NULL))
    {
        if (controller != NULL)
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "SOLVENT_LENNARD_JONES::Get_Local",
                "Reason:\n\tthe hard-only optimized solvent dispatcher "
                "received invalid local topology or buffers\n");
        return;
    }
    if (local_res_numbers == 0 || atom_numbers == 0)
    {
        solvent_start_local = local_res_numbers;
        local_solvent_numbers = 0;
        deviceMemset(d_solvent_start_local, 0, sizeof(int));
        deviceMemset(d_local_solvent_numbers, 0, sizeof(int));
        if (local_res_numbers != 0)
            deviceMemcpy(d_solvent_start_local, &local_res_numbers, sizeof(int),
                         deviceMemcpyHostToDevice);
        return;
    }

    Launch_Device_Kernel(get_local_device, 1, 1, 0, NULL, local_res_numbers,
                         d_res_len, water_points, d_solvent_start_local,
                         d_local_solvent_numbers, atom_numbers, local_mass);
    deviceMemcpy(&solvent_start_local, d_solvent_start_local, sizeof(int),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&local_solvent_numbers, d_local_solvent_numbers, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (solvent_start_local < 0 || solvent_start_local > local_res_numbers ||
        local_solvent_numbers < 0 || local_solvent_numbers > atom_numbers ||
        local_solvent_numbers % water_points != 0 ||
        local_solvent_numbers / water_points !=
            local_res_numbers - solvent_start_local)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "SOLVENT_LENNARD_JONES::Get_Local",
            "Reason:\n\tinvalid local solvent result start/count %d/%d for "
            "%d residues and %d atoms\n",
            solvent_start_local, local_solvent_numbers, local_res_numbers,
            atom_numbers);
    }
}
