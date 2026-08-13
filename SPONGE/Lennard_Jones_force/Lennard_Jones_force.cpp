#include "Lennard_Jones_force.h"

#include <cstdint>

#include "../xponge/load/native/lj.hpp"
#include "../xponge/xponge.h"
#include "pair_activity.h"
// #include "assert.h"

// 由LJ坐标和转化系数求距离
__global__ void Copy_LJ_Type_To_New_Crd(const int atom_numbers,
                                        VECTOR_LJ* new_crd, const int* LJ_type)
{
    SIMPLE_DEVICE_FOR(atom_i, atom_numbers)
    {
        new_crd[atom_i].LJ_type = LJ_type[atom_i];
        new_crd[atom_i].global_atom = atom_i;
    }
}

__global__ void Repack_LJ_Crd(const int atom_numbers, const VECTOR* crd,
                              float4* lj_crd_q)
{
    SIMPLE_DEVICE_FOR(atom_i, atom_numbers)
    {
        const VECTOR r = crd[atom_i];
        lj_crd_q[atom_i].x = r.x;
        lj_crd_q[atom_i].y = r.y;
        lj_crd_q[atom_i].z = r.z;
    }
}

static __global__ void device_add(float* variable, const float adder)
{
    variable[0] += adder;
}

template <bool need_force, bool need_energy, bool need_virial,
          bool need_coulomb>
static __global__ void Lennard_Jones_And_Direct_Coulomb_Device(
    const int local_atom_numbers, const int solvent_numbers,
    const ATOM_GROUP* nl, const float4* crd_q, const int2* type_g,
    const LTMatrix3 cell, const LTMatrix3 rcell, const float* LJ_type_A,
    const float* LJ_type_B, const float cutoff, VECTOR* frc,
    const float pme_beta, float* atom_energy, LTMatrix3* atom_virial,
    float* atom_direct_cf_energy, float* atom_LJ_ene, int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = 0 + blockDim.y * blockIdx.x + threadIdx.y;
    if (atom_i < local_atom_numbers - solvent_numbers)
#else
#pragma omp parallel for schedule(dynamic)
    for (int atom_i = 0; atom_i < local_atom_numbers - solvent_numbers;
         atom_i++)
#endif
    {
        VECTOR frc_record = {0.0f, 0.0f, 0.0f};
        LTMatrix3 virial = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        float energy_lj = 0.0f;
        float energy_coulomb = 0.0f;
        float energy_total = 0.0f;
        ATOM_GROUP nl_i = nl[atom_i];
        VECTOR_LJ r1 = Load_VECTOR_LJ(crd_q, type_g, atom_i);
#ifdef USE_GPU
        for (int j = threadIdx.x; j < nl_i.atom_numbers; j += blockDim.x)
#else
        for (int j = 0; j < nl_i.atom_numbers; j += 1)
#endif
        {
            int atom_j = nl_i.atom_serial[j];
            float ij_factor = atom_j < local_atom_numbers ? 1.0f : 0.5f;
            VECTOR_LJ r2 = Load_VECTOR_LJ(crd_q, type_g, atom_j);
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
                        need_coulomb && PairwiseInteraction::Coulomb_Is_Active(
                                            r1.charge, r2.charge));
                if (!activity.Any())
                {
                    continue;
                }
                if (dr_abs == 0.0f)
                {
                    PairwiseInteraction::Fail_Exact_Overlap(
                        r1.global_atom, r2.global_atom,
                        PairwiseInteraction::Components(activity.lennard_jones,
                                                        activity.coulomb),
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
                        float frc_cf_abs =
                            Get_Direct_Coulomb_Force(r1, r2, dr_abs, pme_beta);
                        frc_abs = frc_abs - frc_cf_abs;
                    }
                    VECTOR frc_lin = frc_abs * dr;
                    frc_record = frc_record + frc_lin;
                    if (atom_j < local_atom_numbers)
                    {
                        atomicAdd(frc + atom_j, -frc_lin);
                    }
                    if (need_virial)
                    {
                        virial = virial - ij_factor * Get_Virial_From_Force_Dis(
                                                          frc_lin, dr);
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
                            ij_factor *
                            Get_Direct_Coulomb_Energy(r1, r2, dr_abs, pme_beta);
                    }
                }
            }
        }
        energy_total = energy_lj + energy_coulomb;
        if (need_force)
        {
            Warp_Sum_To(frc + atom_i, frc_record, warpSize);
        }
        if (need_energy)
        {
            Warp_Sum_To(atom_energy + atom_i, energy_total, warpSize);
            Warp_Sum_To(atom_LJ_ene + atom_i, energy_lj, warpSize);
            if (need_coulomb)
                Warp_Sum_To(atom_direct_cf_energy + atom_i, energy_coulomb,
                            warpSize);
        }
        if (need_virial)
        {
            Warp_Sum_To(atom_virial + atom_i, virial, warpSize);
        }
    }
}

void LENNARD_JONES_INFORMATION::LJ_Malloc()
{
    Malloc_Safely((void**)&h_atom_LJ_type, sizeof(int) * atom_numbers);
    Malloc_Safely((void**)&h_LJ_A, sizeof(float) * pair_type_numbers);
    Malloc_Safely((void**)&h_LJ_B, sizeof(float) * pair_type_numbers);
    Malloc_Safely((void**)&h_LJ_energy_atom, sizeof(float) * atom_numbers);
}

void LENNARD_JONES_INFORMATION::Initial(CONTROLLER* controller, float cutoff,
                                        const char* module_name)
{
    this->controller = controller;
    if (module_name == NULL)
    {
        strcpy(this->module_name, "LJ");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }
    controller->printf("START INITIALIZING LENNADR JONES INFORMATION:\n");
    const auto& lj = Xponge::system.classical_force_field.lj;
    Xponge::LennardJones local_lj;
    const Xponge::LennardJones* lj_to_use = NULL;
    if (module_name == NULL)
    {
        lj_to_use = &lj;
    }
    else if (controller->Command_Exist(this->module_name, "in_file"))
    {
        Xponge::Native_Load_LJ(&local_lj, controller,
                               Xponge::Load_Get_Atom_Numbers(&Xponge::system),
                               this->module_name);
        lj_to_use = &local_lj;
    }
    if (lj_to_use != NULL)
    {
        if (lj_to_use->atom_type.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_INFORMATION::Initial",
                "Reason:\n\tLJ atom count cannot be represented by the "
                "runtime index type\n");
            return;
        }
        atom_numbers = static_cast<int>(lj_to_use->atom_type.size());
        atom_type_numbers = lj_to_use->atom_type_numbers;
        if (atom_type_numbers < 0 || atom_type_numbers > 65535)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_INFORMATION::Initial",
                "Reason:\n\tLJ atom type count %d is outside the "
                "representable range [0, 65535]\n",
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
            lj_to_use->pair_B.size() != expected_pair_count ||
            (atom_numbers > 0 && atom_type_numbers == 0))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_INFORMATION::Initial",
                "Reason:\n\tLJ table shape is inconsistent: atoms=%d, "
                "types=%d, expected pairs=%llu, A/B pairs=%llu/%llu\n",
                atom_numbers, atom_type_numbers,
                static_cast<unsigned long long>(expected_pair_count),
                static_cast<unsigned long long>(lj_to_use->pair_A.size()),
                static_cast<unsigned long long>(lj_to_use->pair_B.size()));
            return;
        }
        pair_type_numbers = static_cast<int>(expected_pair_count);
        for (int atom = 0; atom < atom_numbers; atom++)
        {
            const int type = lj_to_use->atom_type[atom];
            if (type < 0 || type >= atom_type_numbers)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand,
                    "LENNARD_JONES_INFORMATION::Initial",
                    "Reason:\n\tLJ atom %d has type %d outside [0, %d)\n", atom,
                    type, atom_type_numbers);
                return;
            }
        }
        for (int pair = 0; pair < pair_type_numbers; pair++)
        {
            if (!Float_Memory_Is_Finite(lj_to_use->pair_A.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_to_use->pair_A.data() +
                                                pair) ||
                !Float_Memory_Is_Finite(lj_to_use->pair_B.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_to_use->pair_B.data() +
                                                pair))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand,
                    "LENNARD_JONES_INFORMATION::Initial",
                    "Reason:\n\tLJ pair %d contains a non-finite or "
                    "subnormal coefficient\n",
                    pair);
                return;
            }
        }
    }
    if (atom_numbers > 0)
    {
        controller->printf("    atom_numbers is %d\n", atom_numbers);
        controller->printf("    atom_LJ_type_number is %d\n",
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
    if (is_initialized)
    {
        this->cutoff = cutoff;
        Device_Malloc_Safely((void**)&crd_with_LJ_parameters,
                             sizeof(VECTOR_LJ) * atom_numbers);
        Launch_Device_Kernel(
            Copy_LJ_Type_To_New_Crd,
            (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
            crd_with_LJ_parameters, d_atom_LJ_type);
        controller->printf("    Start initializing long range LJ correction\n");
        // 全对求和 Σ_i Σ_j B[type_i, type_j] 等于按类型直方图的
        // Σ_a count_a · Σ_b count_b · B[pair(a,b)]，后者按固定顺序双精度
        // 累加，结果确定；此前的全对 kernel 复杂度为 O(N²) 且依赖 float
        // 原子加顺序，本身就有运行间波动。
        std::vector<int64_t> type_count(atom_type_numbers, 0);
        for (int i = 0; i < atom_numbers; i++)
        {
            type_count[h_atom_LJ_type[i]] += 1;
        }
        double c6_sum = 0.0;
        for (int itype = 0; itype < atom_type_numbers; itype++)
        {
            if (type_count[itype] == 0) continue;
            double inner_sum = 0.0;
            for (int jtype = 0; jtype < atom_type_numbers; jtype++)
            {
                if (type_count[jtype] == 0) continue;
                inner_sum += static_cast<double>(type_count[jtype]) *
                             static_cast<double>(
                                 h_LJ_B[Get_LJ_Type(itype, jtype)]);
            }
            c6_sum += static_cast<double>(type_count[itype]) * inner_sum;
        }
        long_range_factor = static_cast<float>(c6_sum);
        printf("        Total C6 factor is %e\n", long_range_factor);

        long_range_factor *=
            -2.0f / 3.0f * CONSTANT_Pi / cutoff / cutoff / cutoff / 6.0f;
        controller->printf("        long range correction factor is: %e\n",
                           long_range_factor);
        controller->printf("    End initializing long range LJ correction\n");
    }
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial("LJ_short", "%.2f");
        controller->Step_Print_Initial("LJ_long", "%.2f");
        controller->Step_Print_Initial("LJ", "%.2f");
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
    }
    controller->printf("END INITIALIZING LENNADR JONES INFORMATION\n\n");
}

static __global__ void get_local_device(int* atom_local, int local_atom_numbers,
                                        int ghost_numbers,
                                        int global_atom_numbers,
                                        int* d_atom_LJ_type,
                                        const float* charge, float4* d_lj_crd_q,
                                        int2* d_lj_type_g,
                                        int* invalid_local_index)
{
    SIMPLE_DEVICE_FOR(i, local_atom_numbers + ghost_numbers)
    {
        int atom_i = atom_local[i];
        if (atom_i < 0 || atom_i >= global_atom_numbers)
        {
            atomicExch(invalid_local_index, i);
        }
        else
        {
            d_lj_type_g[i].x = d_atom_LJ_type[atom_i];
            d_lj_type_g[i].y = atom_i;
            d_lj_crd_q[i].w = charge[i];
        }
    }
}

void LENNARD_JONES_INFORMATION::Get_Local(int* atom_local,
                                          int local_atom_numbers,
                                          int ghost_numbers,
                                          const float* charge)
{
    if (!is_initialized) return;
    local_metadata_is_ready = false;
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > atom_numbers ||
        ghost_numbers > atom_numbers - local_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "LENNARD_JONES_INFORMATION::Get_Local",
            "Reason:\n\t%s received invalid local/ghost atom counts %d/%d "
            "for %d global atoms\n",
            module_name, local_atom_numbers, ghost_numbers, atom_numbers);
        return;
    }
    this->local_atom_numbers = local_atom_numbers;
    this->ghost_numbers = ghost_numbers;
    const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
    if (local_coordinate_numbers == 0)
    {
        local_metadata_is_ready = true;
        return;
    }
    deviceMemset(d_local_metadata_error, -1, sizeof(int));
    Launch_Device_Kernel(
        get_local_device,
        (local_coordinate_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, atom_local, local_atom_numbers,
        ghost_numbers, atom_numbers, d_atom_LJ_type, charge, d_lj_crd_q,
        d_lj_type_g, d_local_metadata_error);
    int invalid_local_index = -1;
    deviceMemcpy(&invalid_local_index, d_local_metadata_error, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (invalid_local_index >= 0)
    {
        int invalid_global_atom = -1;
        deviceMemcpy(&invalid_global_atom, atom_local + invalid_local_index,
                     sizeof(int), deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "LENNARD_JONES_INFORMATION::Get_Local",
            "Reason:\n\t%s local coordinate %d maps to global atom %d "
            "outside [0, %d)\n",
            module_name, invalid_local_index, invalid_global_atom,
            atom_numbers);
        return;
    }
    local_metadata_is_ready = true;
}

bool LENNARD_JONES_INFORMATION::Validate_Local_State(const char* error_by,
                                                     int global_atom_numbers,
                                                     int local_atom_numbers,
                                                     int ghost_numbers,
                                                     int solvent_numbers)
{
    if (global_atom_numbers != atom_numbers ||
        local_atom_numbers != this->local_atom_numbers ||
        ghost_numbers != this->ghost_numbers || !local_metadata_is_ready ||
        solvent_numbers < 0 || solvent_numbers > local_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s local nonbond metadata mismatch: call has "
            "global/local/ghost/solvent counts %d/%d/%d/%d, initialized "
            "state has %d/%d/%d and ready=%d\n",
            module_name, global_atom_numbers, local_atom_numbers, ghost_numbers,
            solvent_numbers, atom_numbers, this->local_atom_numbers,
            this->ghost_numbers, static_cast<int>(local_metadata_is_ready));
        return false;
    }
    return true;
}

void LENNARD_JONES_INFORMATION::Reset_Pair_Overlap_Error()
{
#ifndef GPU_ARCH_NAME
    deviceMemset(d_pair_overlap_error, 0, 3 * sizeof(int));
#endif
}

bool LENNARD_JONES_INFORMATION::Check_Pair_Overlap_Error(const char* error_by)
{
#ifndef GPU_ARCH_NAME
    int overlap_error[3] = {0, -1, -1};
    deviceMemcpy(overlap_error, d_pair_overlap_error, sizeof(overlap_error),
                 deviceMemcpyDeviceToHost);
    if (overlap_error[0] != PairwiseInteraction::PAIR_COMPONENT_NONE)
    {
        const char* component =
            overlap_error[0] ==
                    (PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES |
                     PairwiseInteraction::PAIR_COMPONENT_COULOMB)
                ? "LJ and Coulomb components"
            : overlap_error[0] &
                    PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES
                ? "LJ component"
                : "Coulomb component";
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s global atoms %d %d overlap exactly with active "
            "%s; an active hard nonbond pair has undefined inverse-distance "
            "energy and force\n",
            module_name, overlap_error[1], overlap_error[2], component);
        return true;
    }
#else
    (void)error_by;
#endif
    return false;
}

static __global__ void Long_Range_Virial_Correction(LTMatrix3* d_virial,
                                                    const float factor)
{
    d_virial[0].a11 += factor;
    d_virial[0].a22 += factor;
    d_virial[0].a33 += factor;
}

void LENNARD_JONES_INFORMATION::Long_Range_Correction(int need_pressure,
                                                      LTMatrix3* d_virial,
                                                      int need_potential,
                                                      float* d_potential,
                                                      const float volume)
{
    if (is_initialized && CONTROLLER::PP_MPI_rank == 0)
    {
        if (need_pressure)
        {
            Launch_Device_Kernel(Long_Range_Virial_Correction, 1, 1, 0, 0,
                                 d_virial, 2 * long_range_factor / volume);
        }
        if (need_potential)
        {
            Launch_Device_Kernel(device_add, 1, 1, 0, 0, d_potential,
                                 long_range_factor / volume);

            h_LJ_long_energy = long_range_factor / volume;
        }
    }
}

void LENNARD_JONES_INFORMATION::Parameter_Host_To_Device()
{
    Device_Malloc_And_Copy_Safely((void**)&d_atom_LJ_type, h_atom_LJ_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_A, h_LJ_A,
                                  sizeof(float) * pair_type_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_B, h_LJ_B,
                                  sizeof(float) * pair_type_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_sum, h_LJ_energy_atom,
                                  sizeof(float));
    Device_Malloc_Safely((void**)&d_LJ_energy_atom,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_lj_crd_q, sizeof(float4) * atom_numbers);
    Device_Malloc_Safely((void**)&d_lj_type_g, sizeof(int2) * atom_numbers);
    Device_Malloc_Safely((void**)&d_local_metadata_error, sizeof(int));
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_pair_overlap_error, 3 * sizeof(int));
#endif
}

void LENNARD_JONES_INFORMATION::LJ_PME_Direct_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const int local_atom_numbers,
    const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
    const float* charge, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, const float pme_beta,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_virial, float* atom_direct_pme_energy)
{
    if (is_initialized)
    {
        if (!Validate_Local_State(
                "LENNARD_JONES_INFORMATION::"
                "LJ_PME_Direct_Force_With_Atom_Energy_And_Virial",
                atom_numbers, local_atom_numbers, ghost_numbers,
                solvent_numbers))
        {
            return;
        }
        Launch_Device_Kernel(
            Repack_LJ_Crd,
            (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL,
            this->local_atom_numbers + this->ghost_numbers, crd, d_lj_crd_q);
        if (need_atom_energy)
        {
            deviceMemset(atom_direct_pme_energy, 0,
                         sizeof(float) * this->atom_numbers);
            deviceMemset(d_LJ_energy_atom, 0,
                         sizeof(float) * this->atom_numbers);
        }

        if (atom_numbers == 0 || local_atom_numbers == 0) return;

        Reset_Pair_Overlap_Error();

        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        dim3 gridSize = (local_atom_numbers + blockSize.y - 1) / blockSize.y;
        auto f =
            Lennard_Jones_And_Direct_Coulomb_Device<true, false, false, true>;
        if (!need_atom_energy && !need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, false, false,
                                                        true>;
        }
        else if (need_atom_energy && !need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, true, false,
                                                        true>;
        }
        else if (!need_atom_energy && need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, false, true,
                                                        true>;
        }
        else
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, true, true, true>;
        }
        Launch_Device_Kernel(
            f, gridSize, blockSize, 0, NULL, local_atom_numbers,
            solvent_numbers, nl, d_lj_crd_q, d_lj_type_g, cell, rcell,
            d_LJ_A, d_LJ_B, cutoff, frc, pme_beta, atom_energy, atom_virial,
            atom_direct_pme_energy, d_LJ_energy_atom, d_pair_overlap_error);
        Check_Pair_Overlap_Error(
            "LENNARD_JONES_INFORMATION::"
            "LJ_PME_Direct_Force_With_Atom_Energy_And_Virial");
    }
}

void LENNARD_JONES_INFORMATION::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized || CONTROLLER::MPI_rank >= CONTROLLER::PP_MPI_size)
        return;
    Sum_Of_List(d_LJ_energy_atom, d_LJ_energy_sum, atom_numbers);
    deviceMemcpy(&h_LJ_energy_sum, d_LJ_energy_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &h_LJ_energy_sum, 1, MPI_FLOAT, MPI_SUM,
                  CONTROLLER::pp_comm);
#endif
    controller->Step_Print("LJ_short", h_LJ_energy_sum);
    controller->Step_Print("LJ_long", h_LJ_long_energy);
    controller->Step_Print("LJ", h_LJ_energy_sum + h_LJ_long_energy, true);
}
