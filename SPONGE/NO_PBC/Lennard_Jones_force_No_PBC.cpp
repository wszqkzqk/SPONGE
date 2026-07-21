#include "Lennard_Jones_force_No_PBC.h"

#include <cstdint>

#include "../Lennard_Jones_force/Lennard_Jones_force.h"
#include "../Lennard_Jones_force/pair_activity.h"
#include "../xponge/load/native/lj.hpp"
#include "../xponge/xponge.h"

static __global__ void LJ_Force_Device(
    const int atom_numbers, const VECTOR* crd, const int* LJ_types,
    const float* LJ_A, const float* LJ_B, const int* excluded_list_start,
    const int* excluded_list, const int* excluded_atom_numbers, VECTOR* frc,
    int* pair_overlap_error)
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
            // The isolated-system path evaluates every non-excluded pair.
            // Periodic neighbor cutoffs do not define the NOPBC Hamiltonian.
            {
                do
                {
                    const int type_ij =
                        Get_LJ_Type(LJ_types[atom_i], LJ_types[atom_j]);
                    const float coefficient_a = LJ_A[type_ij];
                    const float coefficient_b = LJ_B[type_ij];
                    if (!PairwiseInteraction::Lennard_Jones_Is_Active(
                            coefficient_a, coefficient_b))
                    {
                        break;
                    }
                    if (dr.x == 0.0f && dr.y == 0.0f && dr.z == 0.0f)
                    {
                        PairwiseInteraction::Fail_Exact_Overlap(
                            atom_i, atom_j,
                            PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES,
                            pair_overlap_error);
                        break;
                    }
                    const float dr_1 = 1.0f / dr_abs;
                    float dr_2 = dr_1 * dr_1;
                    float dr_4 = dr_2 * dr_2;
                    float dr_6 = dr_4 * dr_2;
                    float dr_8 = dr_4 * dr_4;

                    float frc_abs =
                        (-coefficient_a * dr_6 + coefficient_b) * dr_8;
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

static __global__ void LJ_Force_Energy_Device(
    const int atom_numbers, const VECTOR* crd, const int* LJ_types,
    const float* LJ_A, const float* LJ_B, const int* excluded_list_start,
    const int* excluded_list, const int* excluded_atom_numbers, float* atom_ene,
    VECTOR* frc, float* this_ene,
    int* pair_overlap_error)
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
            // See LJ_Force_Device: this direct path has no distance cutoff.
            {
                do
                {
                    const int type_ij =
                        Get_LJ_Type(LJ_types[atom_i], LJ_types[atom_j]);
                    const float coefficient_a = LJ_A[type_ij];
                    const float coefficient_b = LJ_B[type_ij];
                    if (!PairwiseInteraction::Lennard_Jones_Is_Active(
                            coefficient_a, coefficient_b))
                    {
                        break;
                    }
                    if (dr.x == 0.0f && dr.y == 0.0f && dr.z == 0.0f)
                    {
                        PairwiseInteraction::Fail_Exact_Overlap(
                            atom_i, atom_j,
                            PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES,
                            pair_overlap_error);
                        break;
                    }
                    const float dr_1 = 1.0f / dr_abs;
                    float dr_2 = dr_1 * dr_1;
                    float dr_4 = dr_2 * dr_2;
                    float dr_6 = dr_4 * dr_2;
                    float dr_8 = dr_4 * dr_4;

                    float temp_ene = (0.083333333 * coefficient_a * dr_6 -
                                      0.166666666 * coefficient_b) *
                                     dr_6;
                    float frc_abs =
                        (-coefficient_a * dr_6 + coefficient_b) * dr_8;
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

void LENNARD_JONES_NO_PBC_INFORMATION::LJ_Malloc()
{
    Malloc_Safely((void**)&h_LJ_energy_atom, sizeof(float) * atom_numbers);
    Malloc_Safely((void**)&h_atom_LJ_type, sizeof(int) * atom_numbers);
    Malloc_Safely((void**)&h_LJ_A, sizeof(float) * pair_type_numbers);
    Malloc_Safely((void**)&h_LJ_B, sizeof(float) * pair_type_numbers);

    memset(h_LJ_energy_atom, 0, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_LJ_energy_sum, sizeof(float));
    Device_Malloc_Safely((void**)&d_LJ_energy_atom,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_LJ_type, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_LJ_A, sizeof(float) * pair_type_numbers);
    Device_Malloc_Safely((void**)&d_LJ_B, sizeof(float) * pair_type_numbers);
    deviceMemset(d_LJ_energy_sum, 0, sizeof(float));
    deviceMemset(d_LJ_energy_atom, 0, sizeof(float) * atom_numbers);
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_pair_overlap_error, 3 * sizeof(int));
#endif
}

void LENNARD_JONES_NO_PBC_INFORMATION::Initial(CONTROLLER* controller,
                                               const char* module_name)
{
    if (is_initialized)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
            "Reason:\n\tthe NOPBC Lennard-Jones module cannot be "
            "initialized twice without first releasing its state\n");
        return;
    }
    this->controller = controller;
    const char* selected_module_name = module_name == NULL ? "LJ" : module_name;
    if (strlen(selected_module_name) >= sizeof(this->module_name))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
            "Reason:\n\tthe NOPBC Lennard-Jones module name is too long\n");
        return;
    }
    strcpy(this->module_name, selected_module_name);
    controller[0].printf("START INITIALIZING LENNADR JONES INFORMATION:\n");
    const auto& system_lj = Xponge::system.classical_force_field.lj;
    Xponge::LennardJones local_lj;
    const Xponge::LennardJones* lj_to_use = NULL;
    if (module_name == NULL)
    {
        lj_to_use = &system_lj;
    }
    else if (controller->Command_Exist(this->module_name, "in_file"))
    {
        Xponge::Native_Load_LJ(&local_lj, controller,
                               Xponge::Load_Get_Atom_Numbers(&Xponge::system),
                               this->module_name);
        lj_to_use = &local_lj;
    }
    if (lj_to_use != NULL && !lj_to_use->atom_type.empty())
    {
        if (lj_to_use->atom_type.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
                "Reason:\n\tLJ atom count cannot be represented by the "
                "runtime index type\n");
            return;
        }
        atom_numbers = static_cast<int>(lj_to_use->atom_type.size());
        const int expected_atom_numbers =
            Xponge::Load_Get_Atom_Numbers(&Xponge::system);
        if (expected_atom_numbers > 0 && atom_numbers != expected_atom_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
                "Reason:\n\tNOPBC LJ table atom count %d differs from the "
                "loaded system atom count %d\n",
                atom_numbers, expected_atom_numbers);
            return;
        }
        atom_type_numbers = lj_to_use->atom_type_numbers;
        if (atom_type_numbers <= 0 || atom_type_numbers > 65535)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
                "Reason:\n\tLJ atom type count %d is outside the "
                "representable range [1, 65535]\n",
                atom_type_numbers);
            return;
        }
        const std::uint64_t type_count =
            static_cast<std::uint64_t>(atom_type_numbers);
        const std::uint64_t expected_pair_count =
            type_count * (type_count + 1ULL) / 2ULL;
        if (expected_pair_count >
                static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
            lj_to_use->pair_A.size() != expected_pair_count ||
            lj_to_use->pair_B.size() != expected_pair_count)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
                "Reason:\n\tLJ table shape is inconsistent: atoms=%d, "
                "types=%d, expected pairs=%llu, A/B pairs=%llu/%llu\n",
                atom_numbers, atom_type_numbers,
                static_cast<unsigned long long>(expected_pair_count),
                static_cast<unsigned long long>(lj_to_use->pair_A.size()),
                static_cast<unsigned long long>(lj_to_use->pair_B.size()));
            return;
        }
        pair_type_numbers = static_cast<int>(expected_pair_count);
        for (int i = 0; i < pair_type_numbers; i++)
        {
            if (!Float_Memory_Is_Finite(lj_to_use->pair_A.data() + i) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_to_use->pair_A.data() + i) ||
                !Float_Memory_Is_Finite(lj_to_use->pair_B.data() + i) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_to_use->pair_B.data() + i))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand,
                    "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
                    "Reason:\n\tLJ pair %d contains a non-finite or "
                    "subnormal coefficient\n",
                    i);
                return;
            }
        }
        for (int i = 0; i < atom_numbers; i++)
        {
            const int type = lj_to_use->atom_type[i];
            if (type < 0 || type >= atom_type_numbers)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand,
                    "LENNARD_JONES_NO_PBC_INFORMATION::Initial",
                    "Reason:\n\tLJ atom %d has type %d outside [0, %d)\n", i,
                    type, atom_type_numbers);
                return;
            }
        }
        controller[0].printf("    atom_numbers is %d\n", atom_numbers);
        controller[0].printf("    atom_LJ_type_number is %d\n",
                             atom_type_numbers);
        LJ_Malloc();
        for (int i = 0; i < pair_type_numbers; i++)
        {
            h_LJ_A[i] = lj_to_use->pair_A[i];
            h_LJ_B[i] = lj_to_use->pair_B[i];
        }
        for (int i = 0; i < atom_numbers; i++)
        {
            h_atom_LJ_type[i] = lj_to_use->atom_type[i];
        }
        Parameter_Host_To_Device();
        is_initialized = 1;
    }
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller[0].Step_Print_Initial(this->module_name, "%.2f");
        is_controller_printf_initialized = 1;
        controller[0].printf("    structure last modify date is %d\n",
                             last_modify_date);
    }
    controller[0].printf("END INITIALIZING LENNADR JONES INFORMATION\n\n");
}

void LENNARD_JONES_NO_PBC_INFORMATION::Parameter_Host_To_Device()
{
    deviceMemcpy(d_LJ_B, h_LJ_B, sizeof(float) * pair_type_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_LJ_A, h_LJ_A, sizeof(float) * pair_type_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_atom_LJ_type, h_atom_LJ_type, sizeof(int) * atom_numbers,
                 deviceMemcpyHostToDevice);
}

void LENNARD_JONES_NO_PBC_INFORMATION::Reset_Pair_Overlap_Error()
{
#ifndef GPU_ARCH_NAME
    deviceMemset(d_pair_overlap_error, 0, 3 * sizeof(int));
#endif
}

bool LENNARD_JONES_NO_PBC_INFORMATION::Check_Pair_Overlap_Error(
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
            "active LJ component; an active hard nonbond pair has undefined "
            "inverse-distance energy and force\n",
            module_name, overlap_error[1], overlap_error[2]);
        return true;
    }
#else
    (void)error_by;
#endif
    return false;
}

void LENNARD_JONES_NO_PBC_INFORMATION::LJ_Force_With_Atom_Energy(
    const int atom_numbers, const VECTOR* crd, VECTOR* frc,
    const int need_atom_energy, float* atom_energy,
    const int* excluded_list_start, const int* excluded_list,
    const int* excluded_atom_numbers)
{
    if (is_initialized)
    {
        if (atom_numbers != this->atom_numbers || crd == NULL || frc == NULL ||
            excluded_list_start == NULL || excluded_list == NULL ||
            excluded_atom_numbers == NULL ||
            (need_atom_energy && atom_energy == NULL))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "LENNARD_JONES_NO_PBC_INFORMATION::"
                "LJ_Force_With_Atom_Energy",
                "Reason:\n\tinvalid NOPBC LJ call: runtime/module atom "
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
        if (!need_atom_energy)
        {
            Launch_Device_Kernel(LJ_Force_Device, gridSize, blockSize, 0, NULL,
                                 atom_numbers, crd, d_atom_LJ_type, d_LJ_A,
                                 d_LJ_B, excluded_list_start, excluded_list,
                                 excluded_atom_numbers, frc,
                                 d_pair_overlap_error);
        }
        else
        {
            deviceMemset(d_LJ_energy_atom, 0, sizeof(float) * atom_numbers);
            Launch_Device_Kernel(LJ_Force_Energy_Device, gridSize, blockSize, 0,
                                 NULL, atom_numbers, crd, d_atom_LJ_type,
                                 d_LJ_A, d_LJ_B, excluded_list_start,
                                 excluded_list, excluded_atom_numbers,
                                 atom_energy, frc, d_LJ_energy_atom,
                                 d_pair_overlap_error);
        }
        Check_Pair_Overlap_Error(
            "LENNARD_JONES_NO_PBC_INFORMATION::LJ_Force_With_Atom_Energy");
    }
}

void LENNARD_JONES_NO_PBC_INFORMATION::Clear()
{
    free(h_atom_LJ_type);
    free(h_LJ_A);
    free(h_LJ_B);
    free(h_LJ_energy_atom);
    h_atom_LJ_type = NULL;
    h_LJ_A = NULL;
    h_LJ_B = NULL;
    h_LJ_energy_atom = NULL;
    Free_Single_Device_Pointer((void**)&d_atom_LJ_type);
    Free_Single_Device_Pointer((void**)&d_LJ_A);
    Free_Single_Device_Pointer((void**)&d_LJ_B);
    Free_Single_Device_Pointer((void**)&d_LJ_energy_atom);
    Free_Single_Device_Pointer((void**)&d_LJ_energy_sum);
    Free_Single_Device_Pointer((void**)&d_pair_overlap_error);
    atom_numbers = 0;
    atom_type_numbers = 0;
    pair_type_numbers = 0;
    h_LJ_energy_sum = 0.0f;
    is_initialized = 0;
    is_controller_printf_initialized = 0;
    controller = NULL;
}

void LENNARD_JONES_NO_PBC_INFORMATION::Step_Print(CONTROLLER* controller)
{
    if (is_initialized)
    {
        Sum_Of_List(d_LJ_energy_atom, d_LJ_energy_sum, atom_numbers);
        deviceMemcpy(&h_LJ_energy_sum, d_LJ_energy_sum, sizeof(float),
                     deviceMemcpyDeviceToHost);
        controller->Step_Print(module_name, h_LJ_energy_sum, true);
    }
}
