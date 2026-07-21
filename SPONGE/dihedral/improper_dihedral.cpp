#include "improper_dihedral.h"

#include "../xponge/load/native/improper_dihedral.hpp"
#include "../xponge/xponge.h"
#include "torsion_geometry.h"
#include "torsion_validation.h"

static __device__ __forceinline__ void Improper_Dihedral_Fail_Invalid_Geometry(
    int global_term, int atom_i, int atom_j, int atom_k, int atom_l,
    int* invalid_geometry_term, bool accumulator_failure = false)
{
#ifdef GPU_ARCH_NAME
    if (accumulator_failure)
    {
        printf(
            "Fatal SPONGE improper dihedral error: global term %d (local "
            "atoms %d %d %d %d) would produce a non-finite "
            "force/energy/virial accumulator.\n",
            global_term, atom_i, atom_j, atom_k, atom_l);
    }
    else
    {
        printf(
            "Fatal SPONGE improper dihedral error: global term %d (local "
            "atoms %d %d %d %d) has undefined, non-finite, or "
            "unrepresentable torsion geometry/energy/force/virial.\n",
            global_term, atom_i, atom_j, atom_k, atom_l);
    }
#if defined(USE_CUDA)
    asm volatile("trap;");
#elif defined(USE_HIP)
    __builtin_trap();
#endif
#else
    atomicExch(invalid_geometry_term,
               accumulator_failure ? -global_term - 2 : global_term);
#endif
}

static __global__ void Dihedral_Force_With_Atom_Energy_And_Virial_Device(
    const int dihedral_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const int local_atom_numbers, const int* atom_a,
    const int* atom_b, const int* atom_c, const int* atom_d, const float* pk,
    const float* phi0, VECTOR* frc, int need_atom_energy, float* ene,
    float* di_ene, int need_virial, LTMatrix3* virial, const int* global_index,
    int* invalid_geometry_term)
{
#ifdef USE_GPU
    int dihedral_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (dihedral_i < dihedral_numbers)
#else
#pragma omp parallel for
    for (int dihedral_i = 0; dihedral_i < dihedral_numbers; dihedral_i++)
#endif
    {
        int atom_i = atom_a[dihedral_i];
        int atom_j = atom_b[dihedral_i];
        int atom_k = atom_c[dihedral_i];
        int atom_l = atom_d[dihedral_i];
        int global_term = global_index[dihedral_i];

        const double temp_pk = static_cast<double>(pk[dihedral_i]);
        const double temp_phi0 = static_cast<double>(phi0[dihedral_i]);

        TORSION_GEOMETRY geometry;
        if (!Compute_Torsion_Geometry(crd[atom_i], crd[atom_j], crd[atom_k],
                                      crd[atom_l], cell, rcell, &geometry))
        {
            Improper_Dihedral_Fail_Invalid_Geometry(global_term, atom_i, atom_j,
                                                    atom_k, atom_l,
                                                    invalid_geometry_term);
            di_ene[dihedral_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        const double delta_phi =
            remainder(geometry.phi - temp_phi0,
                      6.283185307179586476925286766559005768);
        const double energy_double = temp_pk * delta_phi * delta_phi;
        const double derivative = 2.0 * temp_pk * delta_phi;
        float energy = 0.0f;
        VECTOR fi;
        VECTOR fj;
        VECTOR fk;
        VECTOR fl;
        if (!Torsion_Double_Is_Finite(delta_phi) ||
            !Torsion_Double_Is_Finite(derivative) ||
            !Torsion_Checked_Narrow(energy_double, &energy) ||
            !Torsion_Checked_Scale_And_Narrow(-derivative,
                                               geometry.dphi_dri, &fi) ||
            !Torsion_Checked_Scale_And_Narrow(-derivative,
                                               geometry.dphi_drj, &fj) ||
            !Torsion_Checked_Scale_And_Narrow(-derivative,
                                               geometry.dphi_drk, &fk) ||
            !Torsion_Checked_Scale_And_Narrow(-derivative,
                                               geometry.dphi_drl, &fl))
        {
            Improper_Dihedral_Fail_Invalid_Geometry(global_term, atom_i, atom_j,
                                                    atom_k, atom_l,
                                                    invalid_geometry_term);
            di_ene[dihedral_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }
        LTMatrix3 term_virial = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        if (need_virial && atom_i < local_atom_numbers)
        {
            const TORSION_DOUBLE_VECTOR displacement_i =
                Torsion_To_Double_Vector(geometry.drij);
            const TORSION_DOUBLE_VECTOR displacement_k =
                Torsion_To_Double_Vector(geometry.drkj);
            const TORSION_DOUBLE_VECTOR displacement_l = {
                static_cast<double>(geometry.drkj.x) -
                    static_cast<double>(geometry.drkl.x),
                static_cast<double>(geometry.drkj.y) -
                    static_cast<double>(geometry.drkl.y),
                static_cast<double>(geometry.drkj.z) -
                    static_cast<double>(geometry.drkl.z)};
            if (!Torsion_Checked_Term_Virial(
                    fi, displacement_i, fk, displacement_k, fl,
                    displacement_l, &term_virial))
            {
                Improper_Dihedral_Fail_Invalid_Geometry(global_term, atom_i,
                                                        atom_j, atom_k, atom_l,
                                                        invalid_geometry_term);
                di_ene[dihedral_i] = 0.0f;
#ifdef USE_GPU
                return;
#else
                continue;
#endif
            }
        }
        bool accumulation_is_valid = true;
        if (atom_i < local_atom_numbers)
        {
            accumulation_is_valid = Torsion_Finite_Atomic_Add(frc + atom_i, fi);
        }
        if (accumulation_is_valid && atom_j < local_atom_numbers)
        {
            accumulation_is_valid = Torsion_Finite_Atomic_Add(frc + atom_j, fj);
        }
        if (accumulation_is_valid && atom_k < local_atom_numbers)
        {
            accumulation_is_valid = Torsion_Finite_Atomic_Add(frc + atom_k, fk);
        }
        if (accumulation_is_valid && atom_l < local_atom_numbers)
        {
            accumulation_is_valid = Torsion_Finite_Atomic_Add(frc + atom_l, fl);
        }
        if (accumulation_is_valid && need_atom_energy &&
            atom_i < local_atom_numbers)
        {
            accumulation_is_valid = Finite_Atomic_Add(ene + atom_i, energy);
        }
        if (accumulation_is_valid && need_virial && atom_i < local_atom_numbers)
        {
            accumulation_is_valid =
                Torsion_Finite_Atomic_Add(virial + atom_i, term_virial);
        }
        if (!accumulation_is_valid)
        {
            Improper_Dihedral_Fail_Invalid_Geometry(
                global_term, atom_i, atom_j, atom_k, atom_l,
                invalid_geometry_term, true);
            di_ene[dihedral_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }
        if (need_atom_energy)
        {
            di_ene[dihedral_i] = atom_i < local_atom_numbers ? energy : 0.0f;
        }
    }
}

void IMPROPER_DIHEDRAL::Initial(CONTROLLER* controller, const char* module_name)
{
    this->controller = controller;
    if (module_name == NULL)
    {
        strcpy(this->module_name, "improper_dihedral");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }

    char file_name_suffix[CHAR_LENGTH_MAX];
    sprintf(file_name_suffix, "in_file");
    const auto& impropers = Xponge::system.classical_force_field.impropers;
    Xponge::Torsions local_impropers;
    const Xponge::Torsions* impropers_to_use = NULL;
    const char* init_source = NULL;

    if (module_name == NULL)
    {
        impropers_to_use = &impropers;
        init_source = "Xponge::system";
    }
    else if (controller[0].Command_Exist(this->module_name, file_name_suffix))
    {
        Xponge::Native_Load_Impropers(&local_impropers, controller,
                                      this->module_name);
        impropers_to_use = &local_impropers;
    }

    dihedral_numbers = 0;
    if (impropers_to_use != NULL)
    {
        const std::string validation_error = Validate_Torsions(
            *impropers_to_use, Xponge::system.atoms.mass.size(), "improper",
            true);
        if (!validation_error.empty())
        {
            const std::string reason =
                "Reason:\n\tinvalid improper "
                "dihedral data: " +
                validation_error + "\n";
            controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                           "IMPROPER_DIHEDRAL::Initial",
                                           reason.c_str());
        }
        dihedral_numbers = static_cast<int>(impropers_to_use->atom_a.size());
    }
    if (dihedral_numbers > 0)
    {
        if (module_name == NULL)
        {
            controller[0].printf("START INITIALIZING IMPROPER DIHEDRAL (%s):\n",
                                 init_source);
        }
        else
        {
            controller[0].printf(
                "START INITIALIZING IMPROPER DIHEDRAL (%s_%s):\n",
                this->module_name, file_name_suffix);
        }
        controller[0].printf("    dihedral_numbers is %d\n", dihedral_numbers);
        Memory_Allocate();
        for (int i = 0; i < dihedral_numbers; i++)
        {
            h_atom_a[i] = impropers_to_use->atom_a[i];
            h_atom_b[i] = impropers_to_use->atom_b[i];
            h_atom_c[i] = impropers_to_use->atom_c[i];
            h_atom_d[i] = impropers_to_use->atom_d[i];
            h_pk[i] = impropers_to_use->pk[i];
            h_phi0[i] = impropers_to_use->gamc[i];
        }
        Parameter_Host_To_Device();
        is_initialized = 1;
    }
    else
    {
        controller[0].printf("IMPROPER DIHEDRAL IS NOT INITIALIZED\n\n");
    }

    if (is_initialized && !is_controller_printf_initialized)
    {
        controller[0].Step_Print_Initial(this->module_name, "%.2f");
        is_controller_printf_initialized = 1;
        controller[0].printf("    structure last modify date is %d\n",
                             last_modify_date);
    }
    if (is_initialized)
    {
        controller[0].printf("END INITIALIZING IMPROPER DIHEDRAL\n\n");
    }
}

void IMPROPER_DIHEDRAL::Memory_Allocate()
{
    Malloc_Safely((void**)&h_atom_a, sizeof(int) * dihedral_numbers);
    Malloc_Safely((void**)&h_atom_b, sizeof(int) * dihedral_numbers);
    Malloc_Safely((void**)&h_atom_c, sizeof(int) * dihedral_numbers);
    Malloc_Safely((void**)&h_atom_d, sizeof(int) * dihedral_numbers);
    Malloc_Safely((void**)&h_pk, sizeof(float) * dihedral_numbers);
    Malloc_Safely((void**)&h_phi0, sizeof(float) * dihedral_numbers);
    Malloc_Safely((void**)&h_dihedral_ene, sizeof(float) * dihedral_numbers);
    memset(h_dihedral_ene, 0, sizeof(float) * dihedral_numbers);
    Malloc_Safely((void**)&h_sigma_of_dihedral_ene, sizeof(float));
    memset(h_sigma_of_dihedral_ene, 0, sizeof(float));
}

void IMPROPER_DIHEDRAL::Parameter_Host_To_Device()
{
    Device_Malloc_Safely((void**)&d_atom_a, sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_atom_b, sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_atom_c, sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_atom_d, sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_pk, sizeof(float) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_phi0, sizeof(float) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_dihedral_ene,
                         sizeof(float) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_sigma_of_dihedral_ene, sizeof(float));

    Device_Malloc_Safely((void**)&d_atom_a_local,
                         sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_atom_b_local,
                         sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_atom_c_local,
                         sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_atom_d_local,
                         sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_pk_local, sizeof(float) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_phi0_local,
                         sizeof(float) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_global_index_local,
                         sizeof(int) * dihedral_numbers);
    Device_Malloc_Safely((void**)&d_num_dihe_local, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_term, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_atom, sizeof(int));
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_invalid_geometry_term, sizeof(int));
#endif
    deviceMemset(d_num_dihe_local, 0, sizeof(int));

    deviceMemcpy(d_atom_a, h_atom_a, sizeof(int) * dihedral_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_atom_b, h_atom_b, sizeof(int) * dihedral_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_atom_c, h_atom_c, sizeof(int) * dihedral_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_atom_d, h_atom_d, sizeof(int) * dihedral_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_pk, h_pk, sizeof(float) * dihedral_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_phi0, h_phi0, sizeof(float) * dihedral_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemset(d_dihedral_ene, 0, sizeof(float) * dihedral_numbers);
    deviceMemset(d_sigma_of_dihedral_ene, 0, sizeof(float));
}

static __global__ void get_local_device(
    int dihedral_numbers, int local_coordinate_numbers, const int* d_atom_a,
    const int* d_atom_b, const int* d_atom_c, const int* d_atom_d,
    const char* atom_local_label, const int* atom_local_id, int* d_atom_a_local,
    int* d_atom_b_local, int* d_atom_c_local, int* d_atom_d_local,
    const float* d_pk, const float* d_phi0, float* d_pk_local,
    float* d_phi0_local, int* d_global_index_local, int* d_num_dihe_local,
    int* d_invalid_local_term, int* d_invalid_local_atom)
{
#ifdef USE_GPU
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx != 0) return;
#endif
    d_num_dihe_local[0] = 0;
    d_invalid_local_term[0] = -1;
    d_invalid_local_atom[0] = -1;
    for (int i = 0; i < dihedral_numbers; i++)
    {
        const int global_atoms[4] = {d_atom_a[i], d_atom_b[i], d_atom_c[i],
                                     d_atom_d[i]};
        if (atom_local_label[global_atoms[0]] == 1 ||
            atom_local_label[global_atoms[1]] == 1 ||
            atom_local_label[global_atoms[2]] == 1 ||
            atom_local_label[global_atoms[3]] == 1)
        {
            int local_atoms[4];
            for (int atom = 0; atom < 4; atom++)
            {
                local_atoms[atom] = atom_local_id[global_atoms[atom]];
                if (local_atoms[atom] < 0 ||
                    local_atoms[atom] >= local_coordinate_numbers)
                {
                    d_invalid_local_term[0] = i;
                    d_invalid_local_atom[0] = global_atoms[atom];
                    return;
                }
            }
            const int local_term = d_num_dihe_local[0];
            d_atom_a_local[local_term] = local_atoms[0];
            d_atom_b_local[local_term] = local_atoms[1];
            d_atom_c_local[local_term] = local_atoms[2];
            d_atom_d_local[local_term] = local_atoms[3];
            d_pk_local[local_term] = d_pk[i];
            d_phi0_local[local_term] = d_phi0[i];
            d_global_index_local[local_term] = i;
            d_num_dihe_local[0]++;
        }
    }
}

void IMPROPER_DIHEDRAL::Get_Local(int* atom_local, int local_atom_numbers,
                                  int ghost_numbers, char* atom_local_label,
                                  int* atom_local_id)
{
    if (!is_initialized) return;
    (void)atom_local;
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > std::numeric_limits<int>::max() - ghost_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "IMPROPER_DIHEDRAL::Get_Local",
            "Reason:\n\t%s received invalid local/ghost atom counts %d/%d\n",
            module_name, local_atom_numbers, ghost_numbers);
        return;
    }
    const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
    num_dihe_local = 0;
    this->local_atom_numbers = local_atom_numbers;
    Launch_Device_Kernel(get_local_device, 1, 1, 0, NULL, dihedral_numbers,
                         local_coordinate_numbers, d_atom_a, d_atom_b, d_atom_c,
                         d_atom_d, atom_local_label, atom_local_id,
                         d_atom_a_local, d_atom_b_local, d_atom_c_local,
                         d_atom_d_local, d_pk, d_phi0, d_pk_local, d_phi0_local,
                         d_global_index_local, d_num_dihe_local,
                         d_invalid_local_term, d_invalid_local_atom);
    int invalid_term = -1;
    int invalid_atom = -1;
    deviceMemcpy(&invalid_term, d_invalid_local_term, sizeof(int),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&invalid_atom, d_invalid_local_atom, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (invalid_term >= 0)
    {
        int invalid_local_id = -1;
        deviceMemcpy(&invalid_local_id, atom_local_id + invalid_atom,
                     sizeof(int), deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "IMPROPER_DIHEDRAL::Get_Local",
            "Reason:\n\t%s improper dihedral term %d (global atoms %d %d %d "
            "%d) maps global atom %d to local index %d outside the valid "
            "owned/ghost range [0, %d) on this domain\n",
            module_name, invalid_term, h_atom_a[invalid_term],
            h_atom_b[invalid_term], h_atom_c[invalid_term],
            h_atom_d[invalid_term], invalid_atom, invalid_local_id,
            local_coordinate_numbers);
        return;
    }
    deviceMemcpy(&num_dihe_local, d_num_dihe_local, sizeof(int),
                 deviceMemcpyDeviceToHost);
}

void IMPROPER_DIHEDRAL::Dihedral_Force_With_Atom_Energy_And_Virial(
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* frc,
    int need_atom_energy, float* atom_energy, int need_virial,
    LTMatrix3* atom_virial)
{
    if (is_initialized && num_dihe_local > 0)
    {
#ifndef GPU_ARCH_NAME
        deviceMemset(d_invalid_geometry_term, -1, sizeof(int));
#endif
        Launch_Device_Kernel(
            Dihedral_Force_With_Atom_Energy_And_Virial_Device,
            (num_dihe_local - 1) / CONTROLLER::device_max_thread + 1,
            CONTROLLER::device_max_thread, 0, NULL, this->num_dihe_local, crd,
            cell, rcell, this->local_atom_numbers, this->d_atom_a_local,
            this->d_atom_b_local, this->d_atom_c_local, this->d_atom_d_local,
            this->d_pk_local, this->d_phi0_local, frc, need_atom_energy,
            atom_energy, d_dihedral_ene, need_virial, atom_virial,
            this->d_global_index_local, d_invalid_geometry_term);
#ifndef GPU_ARCH_NAME
        int invalid_code = -1;
        deviceMemcpy(&invalid_code, d_invalid_geometry_term, sizeof(int),
                     deviceMemcpyDeviceToHost);
        if (invalid_code != -1)
        {
            const bool accumulator_failure = invalid_code < -1;
            const int invalid_term =
                accumulator_failure ? -(invalid_code + 2) : invalid_code;
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "IMPROPER_DIHEDRAL::"
                "Dihedral_Force_With_Atom_Energy_And_Virial",
                "Reason:\n\t%s improper dihedral term %d (global atoms %d "
                "%d %d %d) %s\n",
                module_name, invalid_term, h_atom_a[invalid_term],
                h_atom_b[invalid_term], h_atom_c[invalid_term],
                h_atom_d[invalid_term],
                accumulator_failure
                    ? "would produce a non-finite force/energy/virial "
                      "accumulator"
                    : "has undefined, non-finite, or unrepresentable torsion "
                      "geometry/energy/force/virial");
            return;
        }
#endif
    }
}

void IMPROPER_DIHEDRAL::Step_Print(CONTROLLER* controller, bool print_sum)
{
    if (is_initialized && CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        Sum_Of_List(d_dihedral_ene, d_sigma_of_dihedral_ene,
                    num_dihe_local);  // 修改为局部求和
        deviceMemcpy(h_sigma_of_dihedral_ene, d_sigma_of_dihedral_ene,
                     sizeof(float), deviceMemcpyDeviceToHost);
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, h_sigma_of_dihedral_ene, 1, MPI_FLOAT,
                      MPI_SUM, CONTROLLER::pp_comm);
#endif
        if (!Float_Memory_Is_Finite(h_sigma_of_dihedral_ene))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "IMPROPER_DIHEDRAL::Step_Print",
                "Reason:\n\t%s improper dihedral total energy is non-finite "
                "after local/MPI reduction\n",
                module_name);
            return;
        }
        if (CONTROLLER::MPI_rank == 0)
        {
            controller->Step_Print(this->module_name, h_sigma_of_dihedral_ene,
                                   print_sum);
        }
    }
}
