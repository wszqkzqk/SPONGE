#include "Coulomb_Force_No_PBC.h"

#include "../Lennard_Jones_force/pair_activity.h"

static __global__ void Coulomb_Force_Device(
    const int atom_numbers, const VECTOR* crd, const float* charge,
    const int* excluded_list_start, const int* excluded_list,
    const int* excluded_atom_numbers, VECTOR* frc, int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = blockDim.y * blockIdx.y + threadIdx.y;
    int atom_j = atom_i + 1 + blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers && atom_j < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
        for (int atom_j = atom_i + 1; atom_j < atom_numbers; atom_j++)
#endif
    {
        int tocal = 1;
        const int* start = excluded_list + excluded_list_start[atom_i];
        for (int k = 0; tocal == 1 && k < excluded_atom_numbers[atom_i]; k += 1)
        {
            if (start[k] == atom_j) tocal = 0;
        }
        if (tocal == 1)
        {
            VECTOR dr = crd[atom_j] - crd[atom_i];
            const float dr_abs = norm3df(dr.x, dr.y, dr.z);
            // The direct NOPBC Hamiltonian contains every non-excluded pair.
            {
                do
                {
                    const float charge_i = charge[atom_i];
                    const float charge_j = charge[atom_j];
                    if (!PairwiseInteraction::Coulomb_Is_Active(charge_i,
                                                                charge_j))
                    {
                        break;
                    }
                    if (dr.x == 0.0f && dr.y == 0.0f && dr.z == 0.0f)
                    {
                        PairwiseInteraction::Fail_Exact_Overlap(
                            atom_i, atom_j,
                            PairwiseInteraction::PAIR_COMPONENT_COULOMB,
                            pair_overlap_error);
                        break;
                    }
                    float dr_1 = 1.0f / dr_abs;
                    float dr_2 = dr_1 * dr_1;
                    float dr_3 = dr_1 * dr_2;
                    float chargeij = charge_i * charge_j;
                    float frc_abs = -chargeij * dr_3;
                    VECTOR temp_frc = frc_abs * dr;

                    atomicAdd(&frc[atom_j].x, -temp_frc.x);
                    atomicAdd(&frc[atom_j].y, -temp_frc.y);
                    atomicAdd(&frc[atom_j].z, -temp_frc.z);
                    atomicAdd(&frc[atom_i].x, temp_frc.x);
                    atomicAdd(&frc[atom_i].y, temp_frc.y);
                    atomicAdd(&frc[atom_i].z, temp_frc.z);
                } while (false);
            }
        }
    }
}

static __global__ void Coulomb_Force_Energy_Device(
    const int atom_numbers, const VECTOR* crd, const float* charge,
    const int* excluded_list_start, const int* excluded_list,
    const int* excluded_atom_numbers, float* atom_ene, VECTOR* frc,
    float* this_ene, int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = blockDim.y * blockIdx.y + threadIdx.y;
    int atom_j = atom_i + 1 + blockDim.x * blockIdx.x + threadIdx.x;
    if (atom_i < atom_numbers && atom_j < atom_numbers)
#else
#pragma omp parallel for
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
        for (int atom_j = atom_i + 1; atom_j < atom_numbers; atom_j++)
#endif
    {
        int tocal = 1;
        const int* start = excluded_list + excluded_list_start[atom_i];
        for (int k = 0; tocal == 1 && k < excluded_atom_numbers[atom_i]; k += 1)
        {
            if (start[k] == atom_j) tocal = 0;
        }
        if (tocal == 1)
        {
            VECTOR dr = crd[atom_j] - crd[atom_i];
            const float dr_abs = norm3df(dr.x, dr.y, dr.z);
            // See Coulomb_Force_Device: NOPBC has no distance cutoff.
            {
                do
                {
                    const float charge_i = charge[atom_i];
                    const float charge_j = charge[atom_j];
                    if (!PairwiseInteraction::Coulomb_Is_Active(charge_i,
                                                                charge_j))
                    {
                        break;
                    }
                    if (dr.x == 0.0f && dr.y == 0.0f && dr.z == 0.0f)
                    {
                        PairwiseInteraction::Fail_Exact_Overlap(
                            atom_i, atom_j,
                            PairwiseInteraction::PAIR_COMPONENT_COULOMB,
                            pair_overlap_error);
                        break;
                    }
                    float dr_1 = 1.0f / dr_abs;
                    float dr_2 = dr_1 * dr_1;
                    float dr_3 = dr_1 * dr_2;
                    float chargeij = charge_i * charge_j;
                    float temp_ene = chargeij * dr_1;
                    float frc_abs = -chargeij * dr_3;
                    VECTOR temp_frc = frc_abs * dr;

                    atomicAdd(&frc[atom_j].x, -temp_frc.x);
                    atomicAdd(&frc[atom_j].y, -temp_frc.y);
                    atomicAdd(&frc[atom_j].z, -temp_frc.z);
                    atomicAdd(&frc[atom_i].x, temp_frc.x);
                    atomicAdd(&frc[atom_i].y, temp_frc.y);
                    atomicAdd(&frc[atom_i].z, temp_frc.z);

                    atomicAdd(&atom_ene[atom_i], temp_ene);
                    atomicAdd(&this_ene[atom_i], temp_ene);
                } while (false);
            }
        }
    }
}

void COULOMB_FORCE_NO_PBC_INFORMATION::Malloc()
{
    Malloc_Safely((void**)&h_Coulomb_energy_atom, sizeof(float) * atom_numbers);
    memset(h_Coulomb_energy_atom, 0, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_Coulomb_energy_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_Coulomb_energy_atom,
                         sizeof(float) * atom_numbers);
    deviceMemset(d_Coulomb_energy_sum, 0, sizeof(float));
    deviceMemset(d_Coulomb_energy_atom, 0, sizeof(float) * atom_numbers);
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_pair_overlap_error, 3 * sizeof(int));
#endif
}

void COULOMB_FORCE_NO_PBC_INFORMATION::Initial(CONTROLLER* controller,
                                               int atom_numbers,
                                               const char* module_name)
{
    if (is_initialized)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "COULOMB_FORCE_NO_PBC_INFORMATION::Initial",
            "Reason:\n\tthe NOPBC Coulomb module cannot be initialized "
            "twice without first releasing its state\n");
        return;
    }
    this->controller = controller;
    const char* selected_module_name =
        module_name == NULL ? "Coulomb" : module_name;
    if (strlen(selected_module_name) >= sizeof(this->module_name))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "COULOMB_FORCE_NO_PBC_INFORMATION::Initial",
            "Reason:\n\tthe NOPBC Coulomb module name is too long\n");
        return;
    }
    strcpy(this->module_name, selected_module_name);
    controller[0].printf("START INITIALIZING COULOMB INFORMATION:\n");
    if (atom_numbers <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand,
            "COULOMB_FORCE_NO_PBC_INFORMATION::Initial",
            "Reason:\n\tinvalid NOPBC Coulomb atom count: %d\n",
            atom_numbers);
        return;
    }
    this->atom_numbers = atom_numbers;
    this->is_initialized = 1;
    Malloc();
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller[0].Step_Print_Initial(this->module_name, "%.2f");
        is_controller_printf_initialized = 1;
        controller[0].printf("    structure last modify date is %d\n",
                             last_modify_date);
    }
    controller[0].printf("END INITIALIZING COULOMB INFORMATION\n\n");
}

void COULOMB_FORCE_NO_PBC_INFORMATION::Reset_Pair_Overlap_Error()
{
#ifndef GPU_ARCH_NAME
    deviceMemset(d_pair_overlap_error, 0, 3 * sizeof(int));
#endif
}

bool COULOMB_FORCE_NO_PBC_INFORMATION::Check_Pair_Overlap_Error(
    const char* error_by)
{
#ifndef GPU_ARCH_NAME
    int overlap_error[3] = {0, -1, -1};
    deviceMemcpy(overlap_error, d_pair_overlap_error, sizeof(overlap_error),
                 deviceMemcpyDeviceToHost);
    if (overlap_error[0] != PairwiseInteraction::PAIR_COMPONENT_NONE)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s global atoms %d %d overlap exactly with an "
            "active Coulomb component; an active hard nonbond pair has "
            "undefined inverse-distance energy and force\n",
            module_name, overlap_error[1], overlap_error[2]);
        return true;
    }
#else
    (void)error_by;
#endif
    return false;
}

void COULOMB_FORCE_NO_PBC_INFORMATION::Coulomb_Force_With_Atom_Energy(
    const int atom_numbers, const VECTOR* crd, const float* charge, VECTOR* frc,
    const int need_atom_energy, float* atom_energy,
    const int* excluded_list_start, const int* excluded_list,
    const int* excluded_atom_numbers)
{
    if (is_initialized)
    {
        if (atom_numbers != this->atom_numbers || crd == NULL ||
            charge == NULL || frc == NULL || excluded_list_start == NULL ||
            excluded_list == NULL || excluded_atom_numbers == NULL ||
            (need_atom_energy && atom_energy == NULL))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "COULOMB_FORCE_NO_PBC_INFORMATION::"
                "Coulomb_Force_With_Atom_Energy",
                "Reason:\n\tinvalid NOPBC Coulomb call: runtime/module atom "
                "counts %d/%d or required storage is inconsistent\n",
                atom_numbers, this->atom_numbers);
            return;
        }
        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        dim3 gridSize = {(atom_numbers + blockSize.x - 1) / blockSize.x,
                         (atom_numbers + blockSize.y - 1) / blockSize.y};
        Reset_Pair_Overlap_Error();
        if (need_atom_energy == 0)
        {
            Launch_Device_Kernel(
                Coulomb_Force_Device, gridSize, blockSize, 0, NULL,
                atom_numbers, crd, charge, excluded_list_start, excluded_list,
                excluded_atom_numbers, frc, d_pair_overlap_error);
        }
        else
        {
            deviceMemset(d_Coulomb_energy_atom, 0,
                         sizeof(float) * atom_numbers);
            Launch_Device_Kernel(
                Coulomb_Force_Energy_Device, gridSize, blockSize, 0, NULL,
                atom_numbers, crd, charge, excluded_list_start, excluded_list,
                excluded_atom_numbers, atom_energy, frc,
                d_Coulomb_energy_atom, d_pair_overlap_error);
        }
        Check_Pair_Overlap_Error(
            "COULOMB_FORCE_NO_PBC_INFORMATION::"
            "Coulomb_Force_With_Atom_Energy");
    }
}

void COULOMB_FORCE_NO_PBC_INFORMATION::Clear()
{
    free(h_Coulomb_energy_atom);
    h_Coulomb_energy_atom = NULL;
    Free_Single_Device_Pointer((void**)&d_Coulomb_energy_atom);
    Free_Single_Device_Pointer((void**)&d_Coulomb_energy_sum);
    Free_Single_Device_Pointer((void**)&d_pair_overlap_error);
    atom_numbers = 0;
    h_Coulomb_energy_sum = 0.0f;
    is_initialized = 0;
    is_controller_printf_initialized = 0;
    controller = NULL;
}

void COULOMB_FORCE_NO_PBC_INFORMATION::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized) return;
    Sum_Of_List(d_Coulomb_energy_atom, d_Coulomb_energy_sum, atom_numbers);
    deviceMemcpy(&h_Coulomb_energy_sum, d_Coulomb_energy_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
    controller->Step_Print(module_name, h_Coulomb_energy_sum, true);
}
