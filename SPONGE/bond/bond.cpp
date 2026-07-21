#include "bond.h"

#include "../xponge/load/native/bond.hpp"
#include "../xponge/xponge.h"

// 由于，大部分情况下bond的energy和virial计算耗时不显著，为简化bond模块的逻辑复杂度，
// 将bond
// 对原子上的力（frc）、原子上的能量（atom_energy）、原子上的维力值（atom_virial）
// 一并计算。
// 对于简易和轻度修改，可以不用考虑能量与维力值的计算。
//   在不使用涉及维力系数的模拟中，可以不用计算正确的维力值
//   在能量数值不影响模拟的过程中，可以不用计算正确的能量值
//   只有力是最基本的计算要求

static __device__ __forceinline__ bool Bond_Float_Is_Finite(float value)
{
#ifdef GPU_ARCH_NAME
    // CUDA/HIP fast-math must not be able to assume that an arithmetic result
    // is finite and fold this validation away.
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
#elif defined(__GNUC__) || defined(__clang__)
    unsigned int bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires 32-bit IEEE-754 floats");
    memcpy(&bits, &value, sizeof(value));
    __asm__ __volatile__("" : "+r"(bits));
    return (bits & 0x7f800000U) != 0x7f800000U;
#else
    return Float_Memory_Is_Finite(&value);
#endif
}

static __device__ __forceinline__ bool Bond_Vector_Is_Finite(
    const VECTOR& value)
{
    return Bond_Float_Is_Finite(value.x) && Bond_Float_Is_Finite(value.y) &&
           Bond_Float_Is_Finite(value.z);
}

static __device__ __forceinline__ bool Bond_Virial_Is_Finite(
    const LTMatrix3& value)
{
    return Bond_Float_Is_Finite(value.a11) && Bond_Float_Is_Finite(value.a21) &&
           Bond_Float_Is_Finite(value.a22) && Bond_Float_Is_Finite(value.a31) &&
           Bond_Float_Is_Finite(value.a32) && Bond_Float_Is_Finite(value.a33);
}

static __device__ __forceinline__ bool Bond_Finite_Atomic_Add(
    VECTOR* accumulator, const VECTOR& value)
{
    return Finite_Atomic_Add(&accumulator->x, value.x) &&
           Finite_Atomic_Add(&accumulator->y, value.y) &&
           Finite_Atomic_Add(&accumulator->z, value.z);
}

static __device__ __forceinline__ bool Bond_Finite_Atomic_Add(
    LTMatrix3* accumulator, const LTMatrix3& value)
{
    return Finite_Atomic_Add(&accumulator->a11, value.a11) &&
           Finite_Atomic_Add(&accumulator->a21, value.a21) &&
           Finite_Atomic_Add(&accumulator->a22, value.a22) &&
           Finite_Atomic_Add(&accumulator->a31, value.a31) &&
           Finite_Atomic_Add(&accumulator->a32, value.a32) &&
           Finite_Atomic_Add(&accumulator->a33, value.a33);
}

static __device__ __forceinline__ void Bond_Fail_Invalid_Term(
    int global_term, int global_atom_i, int global_atom_j,
    int* invalid_geometry_term, bool accumulator_failure = false)
{
#ifdef GPU_ARCH_NAME
    if (accumulator_failure)
    {
        printf(
            "Fatal SPONGE bond error: global term %d (global atoms %d %d) "
            "would produce a non-finite force/energy/virial accumulator.\n",
            global_term, global_atom_i, global_atom_j);
    }
    else
    {
        printf(
            "Fatal SPONGE bond error: global term %d (global atoms %d %d) "
            "has undefined zero-distance geometry or a non-finite "
            "geometry/energy/force/virial.\n",
            global_term, global_atom_i, global_atom_j);
    }
#if defined(USE_CUDA)
    asm volatile("trap;");
#elif defined(USE_HIP)
    __builtin_trap();
#endif
#else
    // -1 means no failure. A validated term count is at most INT_MAX, so both
    // the geometry and accumulator encodings remain representable.
    atomicExch(invalid_geometry_term,
               accumulator_failure ? -global_term - 2 : global_term);
#endif
}

static __global__ void Bond_Force_With_Atom_Energy_And_Virial_Device(
    const int bond_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const int local_atom_numbers, const int* atom_a,
    const int* atom_b, const float* bond_k, const float* bond_r0, VECTOR* frc,
    int need_atom_energy, float* atom_energy, int need_virial,
    LTMatrix3* atom_virial, float* bond_ene, const int* global_index,
    const int* global_atom_a, const int* global_atom_b,
    int* invalid_geometry_term)
{
#ifdef USE_GPU
    int bond_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (bond_i < bond_numbers)
#else
#pragma omp parallel for
    for (int bond_i = 0; bond_i < bond_numbers; bond_i++)
#endif
    {
        // 获取第bond_i根键的两个连接的原子编号
        // 和键强度、平衡长度
        int atom_i = atom_a[bond_i];
        int atom_j = atom_b[bond_i];
        int global_term = global_index[bond_i];
        int global_i = global_atom_a[global_term];
        int global_j = global_atom_b[global_term];
        float k = bond_k[bond_i];
        float r0 = bond_r0[bond_i];
        bond_ene[bond_i] = 0.0f;

        // An exact-zero coefficient makes the complete term identically zero,
        // including at coincident or otherwise unused coordinates. Avoiding
        // the geometry entirely is also what prevents 0 * undefined from
        // contaminating the force path.
        if (k == 0.0f)
        {
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        // 获取该对原子的考虑周期性边界的最短位置矢量（dr），和最短距离abs_r
        VECTOR dr =
            Get_Periodic_Displacement(crd[atom_i], crd[atom_j], cell, rcell);
        const double r_squared =
            static_cast<double>(dr.x) * static_cast<double>(dr.x) +
            static_cast<double>(dr.y) * static_cast<double>(dr.y) +
            static_cast<double>(dr.z) * static_cast<double>(dr.z);
        if (!Bond_Vector_Is_Finite(dr) || !(r_squared >= 0.0))
        {
            Bond_Fail_Invalid_Term(global_term, global_i, global_j,
                                   invalid_geometry_term);
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        if (r_squared == 0.0)
        {
            // k*r^2 (r0 == 0) is differentiable at the origin with exactly
            // zero energy and force. For r0 != 0, |r| has no unique Cartesian
            // gradient at the origin, so choosing an arbitrary zero direction
            // would silently change the model.
            if (r0 == 0.0f)
            {
#ifdef USE_GPU
                return;
#else
                continue;
#endif
            }
            Bond_Fail_Invalid_Term(global_term, global_i, global_j,
                                   invalid_geometry_term);
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        const double distance = sqrt(r_squared);
        const double displacement = distance - static_cast<double>(r0);
        const double force_magnitude =
            -2.0 * static_cast<double>(k) * displacement;
        const float energy = static_cast<float>(static_cast<double>(k) *
                                                displacement * displacement);
        const VECTOR f = {
            static_cast<float>(force_magnitude *
                               (static_cast<double>(dr.x) / distance)),
            static_cast<float>(force_magnitude *
                               (static_cast<double>(dr.y) / distance)),
            static_cast<float>(force_magnitude *
                               (static_cast<double>(dr.z) / distance))};
        if (!Bond_Float_Is_Finite(energy) || !Bond_Vector_Is_Finite(f))
        {
            Bond_Fail_Invalid_Term(global_term, global_i, global_j,
                                   invalid_geometry_term);
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        LTMatrix3 term_virial = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        if (need_virial && atom_i < local_atom_numbers)
        {
            term_virial = Get_Virial_From_Force_Dis(f, dr);
            if (!Bond_Virial_Is_Finite(term_virial))
            {
                Bond_Fail_Invalid_Term(global_term, global_i, global_j,
                                       invalid_geometry_term);
#ifdef USE_GPU
                return;
#else
                continue;
#endif
            }
        }

        bool accumulation_is_valid = true;
        if (atom_j < local_atom_numbers)
        {
            accumulation_is_valid = Bond_Finite_Atomic_Add(frc + atom_j, -f);
        }
        if (accumulation_is_valid && atom_i < local_atom_numbers)
        {
            accumulation_is_valid = Bond_Finite_Atomic_Add(frc + atom_i, f);
        }
        // 将计算得到的能量和维力值加到该bond中的其中一个原子身上
        // 原理上，该bond能量是不可分的。但是，一般情况，bond相连的两个原子
        // 总是被看作属于一个分子来讨论，因此可以直接将能量和维力值加到其中一个原子上
        if (accumulation_is_valid && need_atom_energy &&
            atom_i < local_atom_numbers)
        {
            accumulation_is_valid =
                Finite_Atomic_Add(atom_energy + atom_i, energy);
        }
        if (accumulation_is_valid && need_virial && atom_i < local_atom_numbers)
        {
            accumulation_is_valid =
                Bond_Finite_Atomic_Add(atom_virial + atom_i, term_virial);
        }
        if (!accumulation_is_valid)
        {
            Bond_Fail_Invalid_Term(global_term, global_i, global_j,
                                   invalid_geometry_term, true);
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }
        if (need_atom_energy && atom_i < local_atom_numbers)
        {
            bond_ene[bond_i] = energy;
        }
    }
}

static std::string Validate_Bonds(const Xponge::Bonds& bonds,
                                  size_t atom_numbers)
{
    const size_t bond_numbers = bonds.atom_a.size();
    if (atom_numbers > static_cast<size_t>(std::numeric_limits<int>::max()))
    {
        return "system atom count exceeds the bond kernel index range";
    }
    if (bonds.atom_b.size() != bond_numbers || bonds.k.size() != bond_numbers ||
        bonds.r0.size() != bond_numbers)
    {
        return "bond arrays have inconsistent lengths";
    }
    if (bond_numbers > static_cast<size_t>(std::numeric_limits<int>::max()))
    {
        return "bond count exceeds the supported int range";
    }
    for (size_t i = 0; i < bond_numbers; i++)
    {
        if (bonds.atom_a[i] < 0 || bonds.atom_b[i] < 0 ||
            static_cast<size_t>(bonds.atom_a[i]) >= atom_numbers ||
            static_cast<size_t>(bonds.atom_b[i]) >= atom_numbers)
        {
            std::ostringstream message;
            message << "bond term " << i << " has atom indices "
                    << bonds.atom_a[i] << ' ' << bonds.atom_b[i]
                    << " outside [0, " << atom_numbers << ')';
            return message.str();
        }
        if (bonds.atom_a[i] == bonds.atom_b[i])
        {
            std::ostringstream message;
            message << "bond term " << i << " repeats global atom "
                    << bonds.atom_a[i] << " as both endpoints";
            return message.str();
        }
        if (!Float_Memory_Is_Finite(&bonds.k[i]) ||
            !Float_Memory_Is_Finite(&bonds.r0[i]))
        {
            std::ostringstream message;
            message << "bond term " << i << " (global atoms " << bonds.atom_a[i]
                    << ' ' << bonds.atom_b[i]
                    << ") has a non-finite force constant or equilibrium "
                       "distance";
            return message.str();
        }
        if (bonds.r0[i] < 0.0f)
        {
            std::ostringstream message;
            message << "bond term " << i << " (global atoms " << bonds.atom_a[i]
                    << ' ' << bonds.atom_b[i]
                    << ") has a negative equilibrium distance";
            return message.str();
        }
    }
    return "";
}

void BOND::Initial(CONTROLLER* controller, CONECT* connectivity,
                   PAIR_DISTANCE* con_dis, const char* module_name)
{
    this->controller = controller;
    // 给予bond模块一个默认名字：bond
    if (module_name == NULL)
    {
        strcpy(this->module_name, "bond");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }

    char file_name_suffix[CHAR_LENGTH_MAX];
    sprintf(file_name_suffix, "in_file");
    const auto& bonds = Xponge::system.classical_force_field.bonds;
    Xponge::Bonds local_bonds;
    const Xponge::Bonds* bonds_to_use = NULL;
    const char* init_source = NULL;
    if (module_name == NULL)
    {
        bonds_to_use = &bonds;
        init_source = "Xponge::system";
    }
    else if (controller->Command_Exist(this->module_name, file_name_suffix))
    {
        Xponge::Native_Load_Bonds(&local_bonds, controller, this->module_name);
        bonds_to_use = &local_bonds;
    }

    if (bonds_to_use != NULL)
    {
        const std::string validation_error =
            Validate_Bonds(*bonds_to_use, Xponge::system.atoms.mass.size());
        if (!validation_error.empty())
        {
            const std::string reason =
                "Reason:\n\tinvalid bond data: " + validation_error + "\n";
            controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                           "BOND::Initial", reason.c_str());
            return;
        }
        bond_numbers = static_cast<int>(bonds_to_use->atom_a.size());
    }
    if (bond_numbers > 0)
    {
        if (module_name == NULL)
        {
            controller->printf("START INITIALIZING BOND (%s):\n", init_source);
        }
        else
        {
            controller->printf("START INITIALIZING BOND (%s_%s):\n",
                               this->module_name, file_name_suffix);
        }
        controller->printf("    bond_numbers is %d\n", bond_numbers);
        Memory_Allocate();
        for (int i = 0; i < bond_numbers; i++)
        {
            h_atom_a[i] = bonds_to_use->atom_a[i];
            h_atom_b[i] = bonds_to_use->atom_b[i];
            h_k[i] = bonds_to_use->k[i];
            h_r0[i] = bonds_to_use->r0[i];
        }
        Parameter_Host_To_Device();
        is_initialized = 1;
    }
    else
    {
        controller->printf("BOND IS NOT INITIALIZED\n\n");
    }

    // 初始化了，且第一次加载用于间隔输出的信息
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial(this->module_name, "%.2f");
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
    }

    // 初始化完成
    if (is_initialized && connectivity)
    {
        for (int i = 0; i < bond_numbers; i += 1)
        {
            connectivity[0][h_atom_a[i]].insert(h_atom_b[i]);
            connectivity[0][h_atom_b[i]].insert(h_atom_a[i]);
            if (h_atom_a[i] < h_atom_b[i])
            {
                con_dis[0][std::pair<int, int>(h_atom_a[i], h_atom_b[i])] =
                    h_r0[i];
            }
            else
            {
                con_dis[0][std::pair<int, int>(h_atom_b[i], h_atom_a[i])] =
                    h_r0[i];
            }
        }
        controller->printf("END INITIALIZING BOND\n\n");
    }
}

void BOND::Memory_Allocate()
{
    Malloc_Safely((void**)&h_atom_a, sizeof(int) * this->bond_numbers);
    Malloc_Safely((void**)&h_atom_b, sizeof(int) * this->bond_numbers);
    Malloc_Safely((void**)&h_k, sizeof(float) * this->bond_numbers);
    Malloc_Safely((void**)&h_r0, sizeof(float) * this->bond_numbers);
    Malloc_Safely((void**)&h_bond_ene, sizeof(float) * this->bond_numbers);
    memset(h_bond_ene, 0, sizeof(float) * this->bond_numbers);
    Malloc_Safely((void**)&h_sigma_of_bond_ene, sizeof(float));
    memset(h_sigma_of_bond_ene, 0, sizeof(float));
}

void BOND::Parameter_Host_To_Device()
{
    Device_Malloc_Safely((void**)&d_atom_a, sizeof(int) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_atom_b, sizeof(int) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_k, sizeof(float) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_r0, sizeof(float) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_bond_ene,
                         sizeof(float) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_sigma_of_bond_ene, sizeof(float));

    Device_Malloc_Safely((void**)&d_atom_a_local,
                         sizeof(int) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_atom_b_local,
                         sizeof(int) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_k_local,
                         sizeof(float) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_r0_local,
                         sizeof(float) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_global_index_local,
                         sizeof(int) * this->bond_numbers);
    Device_Malloc_Safely((void**)&d_num_bond_local, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_term, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_atom, sizeof(int));
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_invalid_geometry_term, sizeof(int));
#endif
    deviceMemset(d_num_bond_local, 0, sizeof(int));

    deviceMemcpy(d_atom_a, h_atom_a, sizeof(int) * this->bond_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_atom_b, h_atom_b, sizeof(int) * this->bond_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_k, h_k, sizeof(float) * this->bond_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_r0, h_r0, sizeof(float) * this->bond_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemset(d_bond_ene, 0, sizeof(float) * this->bond_numbers);
    deviceMemset(d_sigma_of_bond_ene, 0, sizeof(float));
}

static __global__ void get_local_device(
    int bond_numbers, int local_coordinate_numbers, const int* d_atom_a,
    const int* d_atom_b, const char* atom_local_label, const int* atom_local_id,
    int* d_atom_a_local, int* d_atom_b_local, const float* d_k,
    const float* d_r0, float* d_k_local, float* d_r0_local,
    int* d_global_index_local, int* d_num_bond_local, int* d_invalid_local_term,
    int* d_invalid_local_atom)
{
#ifdef USE_GPU
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx != 0) return;
#endif
    d_num_bond_local[0] = 0;
    d_invalid_local_term[0] = -1;
    d_invalid_local_atom[0] = -1;
    for (int i = 0; i < bond_numbers; i++)
    {
        if (atom_local_label[d_atom_a[i]] == 1 ||
            atom_local_label[d_atom_b[i]] == 1)
        {
            const int local_atom_a = atom_local_id[d_atom_a[i]];
            const int local_atom_b = atom_local_id[d_atom_b[i]];
            if (local_atom_a < 0 || local_atom_a >= local_coordinate_numbers)
            {
                d_invalid_local_term[0] = i;
                d_invalid_local_atom[0] = d_atom_a[i];
                return;
            }
            if (local_atom_b < 0 || local_atom_b >= local_coordinate_numbers)
            {
                d_invalid_local_term[0] = i;
                d_invalid_local_atom[0] = d_atom_b[i];
                return;
            }
            const int local_term = d_num_bond_local[0];
            d_atom_a_local[local_term] = local_atom_a;
            d_atom_b_local[local_term] = local_atom_b;
            d_k_local[local_term] = d_k[i];
            d_r0_local[local_term] = d_r0[i];
            d_global_index_local[local_term] = i;
            d_num_bond_local[0]++;
        }
    }
}

void BOND::Get_Local(int* atom_local, int local_atom_numbers, int ghost_numbers,
                     char* atom_local_label, int* atom_local_id)
{
    if (!is_initialized) return;
    (void)atom_local;
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > std::numeric_limits<int>::max() - ghost_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "BOND::Get_Local",
            "Reason:\n\t%s received invalid local/ghost atom counts %d/%d\n",
            module_name, local_atom_numbers, ghost_numbers);
        return;
    }
    const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
    num_bond_local = 0;
    this->local_atom_numbers = local_atom_numbers;
    Launch_Device_Kernel(
        get_local_device, 1, 1, 0, NULL, this->bond_numbers,
        local_coordinate_numbers, this->d_atom_a, this->d_atom_b,
        atom_local_label, atom_local_id, this->d_atom_a_local,
        this->d_atom_b_local, this->d_k, this->d_r0, this->d_k_local,
        this->d_r0_local, this->d_global_index_local, this->d_num_bond_local,
        this->d_invalid_local_term, this->d_invalid_local_atom);
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
            spongeErrorSimulationBreakDown, "BOND::Get_Local",
            "Reason:\n\t%s bond term %d (global atoms %d %d) maps global "
            "atom %d to local index %d outside the valid owned/ghost range "
            "[0, %d) on this domain\n",
            module_name, invalid_term, h_atom_a[invalid_term],
            h_atom_b[invalid_term], invalid_atom, invalid_local_id,
            local_coordinate_numbers);
        return;
    }
    deviceMemcpy(&this->num_bond_local, this->d_num_bond_local, sizeof(int),
                 deviceMemcpyDeviceToHost);
}

void BOND::Bond_Force_With_Atom_Energy_And_Virial(
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* frc,
    int need_atom_energy, float* atom_energy, int need_virial,
    LTMatrix3* atom_virial)
{
    if (is_initialized && num_bond_local > 0)
    {
#ifndef GPU_ARCH_NAME
        deviceMemset(d_invalid_geometry_term, -1, sizeof(int));
#endif
        Launch_Device_Kernel(
            Bond_Force_With_Atom_Energy_And_Virial_Device,
            (num_bond_local - 1) / CONTROLLER::device_max_thread + 1,
            CONTROLLER::device_max_thread, 0, NULL, this->num_bond_local, crd,
            cell, rcell, this->local_atom_numbers, this->d_atom_a_local,
            this->d_atom_b_local, this->d_k_local, this->d_r0_local, frc,
            need_atom_energy, atom_energy, need_virial, atom_virial,
            this->d_bond_ene, this->d_global_index_local, this->d_atom_a,
            this->d_atom_b, this->d_invalid_geometry_term);
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
                "BOND::Bond_Force_With_Atom_Energy_And_Virial",
                "Reason:\n\t%s bond term %d (global atoms %d %d) %s\n",
                module_name, invalid_term, h_atom_a[invalid_term],
                h_atom_b[invalid_term],
                accumulator_failure
                    ? "would produce a non-finite force/energy/virial "
                      "accumulator"
                    : "has undefined zero-distance geometry or a non-finite "
                      "geometry/energy/force/virial");
            return;
        }
#endif
    }
}

void BOND::Step_Print(CONTROLLER* controller)
{
    if (is_initialized && CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        Sum_Of_List(d_bond_ene, d_sigma_of_bond_ene,
                    num_bond_local);  // 局部求和
        deviceMemcpy(h_sigma_of_bond_ene, d_sigma_of_bond_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, h_sigma_of_bond_ene, 1, MPI_FLOAT, MPI_SUM,
                      CONTROLLER::pp_comm);
#endif
        if (!Float_Memory_Is_Finite(h_sigma_of_bond_ene))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "BOND::Step_Print",
                "Reason:\n\t%s bond total energy is non-finite after "
                "local/MPI reduction\n",
                module_name);
            return;
        }
        controller->Step_Print(this->module_name, h_sigma_of_bond_ene, true);
    }
}
