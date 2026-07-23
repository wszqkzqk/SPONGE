#include "angle.h"

#include "../xponge/load/native/angle.hpp"
#include "../xponge/xponge.h"
#include "angle_validation.h"

static __device__ __forceinline__ bool Angle_Float_Is_Finite(float value)
{
#ifdef GPU_ARCH_NAME
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

static __device__ __forceinline__ bool Angle_Double_Is_Finite(double value)
{
#ifdef GPU_ARCH_NAME
    const unsigned long long bits =
        static_cast<unsigned long long>(__double_as_longlong(value));
    return (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
#elif defined(__GNUC__) || defined(__clang__)
    unsigned long long bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires 64-bit IEEE-754 doubles");
    memcpy(&bits, &value, sizeof(value));
    __asm__ __volatile__("" : "+r"(bits));
    return (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
#else
    return Double_Memory_Is_Finite(&value);
#endif
}

static __device__ __forceinline__ bool Angle_Vector_Is_Finite(
    const VECTOR& value)
{
    return Angle_Float_Is_Finite(value.x) && Angle_Float_Is_Finite(value.y) &&
           Angle_Float_Is_Finite(value.z);
}

static __device__ __forceinline__ bool Angle_Virial_Is_Finite(
    const LTMatrix3& value)
{
    return Angle_Float_Is_Finite(value.a11) &&
           Angle_Float_Is_Finite(value.a21) &&
           Angle_Float_Is_Finite(value.a22) &&
           Angle_Float_Is_Finite(value.a31) &&
           Angle_Float_Is_Finite(value.a32) && Angle_Float_Is_Finite(value.a33);
}

static __device__ __forceinline__ bool Angle_Finite_Atomic_Add(
    VECTOR* accumulator, const VECTOR& value)
{
    return Finite_Atomic_Add(&accumulator->x, value.x) &&
           Finite_Atomic_Add(&accumulator->y, value.y) &&
           Finite_Atomic_Add(&accumulator->z, value.z);
}

static __device__ __forceinline__ bool Angle_Finite_Atomic_Add(
    LTMatrix3* accumulator, const LTMatrix3& value)
{
    return Finite_Atomic_Add(&accumulator->a11, value.a11) &&
           Finite_Atomic_Add(&accumulator->a21, value.a21) &&
           Finite_Atomic_Add(&accumulator->a22, value.a22) &&
           Finite_Atomic_Add(&accumulator->a31, value.a31) &&
           Finite_Atomic_Add(&accumulator->a32, value.a32) &&
           Finite_Atomic_Add(&accumulator->a33, value.a33);
}

static __device__ __forceinline__ void Angle_Fail_Invalid_Term(
    int global_term, int global_atom_i, int global_atom_j, int global_atom_k,
    int* invalid_geometry_term, bool accumulator_failure = false)
{
#ifdef GPU_ARCH_NAME
    if (accumulator_failure)
    {
        printf(
            "Fatal SPONGE angle error: global term %d (global atoms %d %d "
            "%d) would produce a non-finite force/energy/virial "
            "accumulator.\n",
            global_term, global_atom_i, global_atom_j, global_atom_k);
    }
    else
    {
        printf(
            "Fatal SPONGE angle error: global term %d (global atoms %d %d "
            "%d) has an undefined zero-arm/collinear geometry or a "
            "non-finite/unrepresentable geometry/energy/force/virial.\n",
            global_term, global_atom_i, global_atom_j, global_atom_k);
    }
#if defined(USE_CUDA)
    asm volatile("trap;");
#elif defined(USE_HIP)
    __builtin_trap();
#endif
#else
    // -1 means no failure. A validated interaction count is at most INT_MAX,
    // so both encodings remain representable for every valid global term.
    atomicExch(invalid_geometry_term,
               accumulator_failure ? -global_term - 2 : global_term);
#endif
}

struct ANGLE_GEOMETRY
{
    float theta;
    bool is_collinear;
    double dtheta_du[3];
    double dtheta_dv[3];
    VECTOR u;
    VECTOR v;
};

// atan2 keeps small nonzero angles distinct from exact endpoints. Double
// intermediates classify collinearity from the float coordinates without an
// empirical cosine clamp, and the cross-product Jacobian avoids 1/sin(theta).
static __device__ __forceinline__ bool Compute_Angle_Geometry(
    const VECTOR& ri, const VECTOR& rj, const VECTOR& rk, const LTMatrix3& cell,
    const LTMatrix3& rcell, ANGLE_GEOMETRY* geometry)
{
    geometry->u = Get_Periodic_Displacement(ri, rj, cell, rcell);
    geometry->v = Get_Periodic_Displacement(rk, rj, cell, rcell);
    if (!Angle_Vector_Is_Finite(geometry->u) ||
        !Angle_Vector_Is_Finite(geometry->v))
    {
        return false;
    }

    const double ux = static_cast<double>(geometry->u.x);
    const double uy = static_cast<double>(geometry->u.y);
    const double uz = static_cast<double>(geometry->u.z);
    const double vx = static_cast<double>(geometry->v.x);
    const double vy = static_cast<double>(geometry->v.y);
    const double vz = static_cast<double>(geometry->v.z);
    const double u_squared = ux * ux + uy * uy + uz * uz;
    const double v_squared = vx * vx + vy * vy + vz * vz;
    if (!(u_squared > 0.0) || !(v_squared > 0.0) ||
        !Angle_Double_Is_Finite(u_squared) ||
        !Angle_Double_Is_Finite(v_squared))
    {
        return false;
    }

    const double nx = uy * vz - uz * vy;
    const double ny = uz * vx - ux * vz;
    const double nz = ux * vy - uy * vx;
    const double normal_squared = nx * nx + ny * ny + nz * nz;
    const double dot = ux * vx + uy * vy + uz * vz;
    if (!(normal_squared >= 0.0) || !Angle_Double_Is_Finite(normal_squared) ||
        !Angle_Double_Is_Finite(dot))
    {
        return false;
    }

    const double normal_length = sqrt(normal_squared);
    const double theta_double = atan2(normal_length, dot);
    geometry->theta = static_cast<float>(theta_double);
    if (!Angle_Double_Is_Finite(normal_length) ||
        !Angle_Double_Is_Finite(theta_double) ||
        !Angle_Float_Is_Finite(geometry->theta) ||
        (theta_double != 0.0 && geometry->theta == 0.0f))
    {
        return false;
    }

    geometry->is_collinear = normal_squared == 0.0;
    if (geometry->is_collinear)
    {
        return dot != 0.0;
    }

    const double inverse_u = 1.0 / sqrt(u_squared);
    const double inverse_v = 1.0 / sqrt(v_squared);
    const double inverse_normal = 1.0 / normal_length;
    // FLT_MAX: std::numeric_limits<float>::max() is a host-only constexpr
    // under nvcc and cannot be called from this device function.
    const double maximum_float = static_cast<double>(FLT_MAX);
    if (!Angle_Double_Is_Finite(inverse_u) ||
        !Angle_Double_Is_Finite(inverse_v) ||
        !Angle_Double_Is_Finite(inverse_normal) || inverse_u > maximum_float ||
        inverse_v > maximum_float)
    {
        return false;
    }

    const double nhat_x = nx * inverse_normal;
    const double nhat_y = ny * inverse_normal;
    const double nhat_z = nz * inverse_normal;
    const double inverse_u_squared = inverse_u * inverse_u;
    const double inverse_v_squared = inverse_v * inverse_v;
    geometry->dtheta_du[0] = (uy * nhat_z - uz * nhat_y) * inverse_u_squared;
    geometry->dtheta_du[1] = (uz * nhat_x - ux * nhat_z) * inverse_u_squared;
    geometry->dtheta_du[2] = (ux * nhat_y - uy * nhat_x) * inverse_u_squared;
    geometry->dtheta_dv[0] = (nhat_y * vz - nhat_z * vy) * inverse_v_squared;
    geometry->dtheta_dv[1] = (nhat_z * vx - nhat_x * vz) * inverse_v_squared;
    geometry->dtheta_dv[2] = (nhat_x * vy - nhat_y * vx) * inverse_v_squared;
    for (int axis = 0; axis < 3; axis++)
    {
        if (!Angle_Double_Is_Finite(geometry->dtheta_du[axis]) ||
            !Angle_Double_Is_Finite(geometry->dtheta_dv[axis]))
        {
            return false;
        }
    }
    return true;
}

static __global__ void Angle_Force_With_Atom_Energy_And_Virial_Device(
    const int angle_numbers, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const int local_atom_numbers, const int* atom_a,
    const int* atom_b, const int* atom_c, const float* angle_k,
    const float* angle_theta0, VECTOR* frc, int need_atom_energy,
    float* atom_energy, float* angle_energy, int need_virial, LTMatrix3* virial,
    const int* global_index, const int* global_atom_a, const int* global_atom_b,
    const int* global_atom_c, int* invalid_geometry_term)
{
#ifdef USE_GPU
    int angle_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (angle_i < angle_numbers)
#else
#pragma omp parallel for
    for (int angle_i = 0; angle_i < angle_numbers; angle_i++)
#endif
    {
        const int atom_i = atom_a[angle_i];
        const int atom_j = atom_b[angle_i];
        const int atom_k = atom_c[angle_i];
        const int global_term = global_index[angle_i];
        const int global_i = global_atom_a[global_term];
        const int global_j = global_atom_b[global_term];
        const int global_k = global_atom_c[global_term];
        const float force_constant = angle_k[angle_i];
        const float theta0 = angle_theta0[angle_i];
        angle_energy[angle_i] = 0.0f;

        if (force_constant == 0.0f)
        {
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        ANGLE_GEOMETRY geometry;
        if (!Compute_Angle_Geometry(crd[atom_i], crd[atom_j], crd[atom_k], cell,
                                    rcell, &geometry))
        {
            Angle_Fail_Invalid_Term(global_term, global_i, global_j, global_k,
                                    invalid_geometry_term);
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        const float dtheta = geometry.theta - theta0;
        if (geometry.is_collinear)
        {
            // k*(theta-theta0)^2 has a unique zero Cartesian derivative at an
            // exact endpoint only when the same float theta gives dtheta==0.
            if (dtheta == 0.0f)
            {
#ifdef USE_GPU
                return;
#else
                continue;
#endif
            }
            Angle_Fail_Invalid_Term(global_term, global_i, global_j, global_k,
                                    invalid_geometry_term);
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        const double derivative = 2.0 * static_cast<double>(force_constant) *
                                  static_cast<double>(dtheta);
        const double energy_double = static_cast<double>(force_constant) *
                                     static_cast<double>(dtheta) *
                                     static_cast<double>(dtheta);
        const double fi_double[3] = {-derivative * geometry.dtheta_du[0],
                                     -derivative * geometry.dtheta_du[1],
                                     -derivative * geometry.dtheta_du[2]};
        const double fk_double[3] = {-derivative * geometry.dtheta_dv[0],
                                     -derivative * geometry.dtheta_dv[1],
                                     -derivative * geometry.dtheta_dv[2]};
        const VECTOR fi = {static_cast<float>(fi_double[0]),
                           static_cast<float>(fi_double[1]),
                           static_cast<float>(fi_double[2])};
        const VECTOR fk = {static_cast<float>(fk_double[0]),
                           static_cast<float>(fk_double[1]),
                           static_cast<float>(fk_double[2])};
        const VECTOR fj = -fi - fk;
        const float energy = static_cast<float>(energy_double);
        if (!Angle_Double_Is_Finite(derivative) ||
            !Angle_Double_Is_Finite(energy_double) ||
            !Angle_Float_Is_Finite(energy) || !Angle_Vector_Is_Finite(fi) ||
            !Angle_Vector_Is_Finite(fj) || !Angle_Vector_Is_Finite(fk))
        {
            Angle_Fail_Invalid_Term(global_term, global_i, global_j, global_k,
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
            term_virial = Get_Virial_From_Force_Dis(fi, geometry.u) +
                          Get_Virial_From_Force_Dis(fk, geometry.v);
            if (!Angle_Virial_Is_Finite(term_virial))
            {
                Angle_Fail_Invalid_Term(global_term, global_i, global_j,
                                        global_k, invalid_geometry_term);
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
            accumulation_is_valid = Angle_Finite_Atomic_Add(frc + atom_i, fi);
        }
        if (accumulation_is_valid && atom_j < local_atom_numbers)
        {
            accumulation_is_valid = Angle_Finite_Atomic_Add(frc + atom_j, fj);
        }
        if (accumulation_is_valid && atom_k < local_atom_numbers)
        {
            accumulation_is_valid = Angle_Finite_Atomic_Add(frc + atom_k, fk);
        }
        if (accumulation_is_valid && need_atom_energy &&
            atom_i < local_atom_numbers)
        {
            accumulation_is_valid =
                Finite_Atomic_Add(atom_energy + atom_i, energy);
        }
        if (accumulation_is_valid && need_virial && atom_i < local_atom_numbers)
        {
            accumulation_is_valid =
                Angle_Finite_Atomic_Add(virial + atom_i, term_virial);
        }
        if (!accumulation_is_valid)
        {
            Angle_Fail_Invalid_Term(global_term, global_i, global_j, global_k,
                                    invalid_geometry_term, true);
            angle_energy[angle_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }
        if (need_atom_energy)
        {
            angle_energy[angle_i] = atom_i < local_atom_numbers ? energy : 0.0f;
        }
    }
}

void ANGLE::Initial(CONTROLLER* controller, const char* module_name)
{
    this->controller = controller;
    if (module_name == NULL)
    {
        strcpy(this->module_name, "angle");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }

    char file_name_suffix[CHAR_LENGTH_MAX];
    sprintf(file_name_suffix, "in_file");
    const auto& angles = Xponge::system.classical_force_field.angles;
    Xponge::Angles local_angles;
    const Xponge::Angles* angles_to_use = NULL;
    const char* init_source = NULL;

    if (module_name == NULL)
    {
        angles_to_use = &angles;
        init_source = "Xponge::system";
    }
    else if (controller->Command_Exist(this->module_name, file_name_suffix))
    {
        Xponge::Native_Load_Angles(&local_angles, controller,
                                   this->module_name);
        angles_to_use = &local_angles;
        init_source =
            controller->Original_Command(this->module_name, file_name_suffix);
    }

    angle_numbers = 0;
    if (angles_to_use != NULL)
    {
        const std::string validation_error = Validate_Angles(
            *angles_to_use, Xponge::system.atoms.mass.size(), "angle");
        if (!validation_error.empty())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "ANGLE::Initial",
                "Reason:\n\tinvalid angle data: %s\n\tSource: %s\n",
                validation_error.c_str(), init_source);
            return;
        }
        angle_numbers = static_cast<int>(angles_to_use->atom_a.size());
    }
    if (angle_numbers > 0)
    {
        if (module_name == NULL)
        {
            controller->printf("START INITIALIZING ANGLE (%s):\n", init_source);
        }
        else
        {
            controller->printf("START INITIALIZING ANGLE (%s_%s):\n",
                               this->module_name, file_name_suffix);
        }
        controller->printf("    angle_numbers is %d\n", angle_numbers);
        Memory_Allocate();
        for (int i = 0; i < angle_numbers; i++)
        {
            h_atom_a[i] = angles_to_use->atom_a[i];
            h_atom_b[i] = angles_to_use->atom_b[i];
            h_atom_c[i] = angles_to_use->atom_c[i];
            h_angle_k[i] = angles_to_use->k[i];
            h_angle_theta0[i] = angles_to_use->theta0[i];
        }
        Parameter_Host_To_Device();
        is_initialized = 1;
    }
    else
    {
        controller->printf("ANGLE IS NOT INITIALIZED\n\n");
    }
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial(this->module_name, "%.2f");
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
    }
    if (is_initialized)
    {
        controller->printf("END INITIALIZING ANGLE\n\n");
    }
}

void ANGLE::Memory_Allocate()
{
    Malloc_Safely((void**)&h_atom_a, sizeof(int) * angle_numbers);
    Malloc_Safely((void**)&h_atom_b, sizeof(int) * angle_numbers);
    Malloc_Safely((void**)&h_atom_c, sizeof(int) * angle_numbers);
    Malloc_Safely((void**)&h_angle_k, sizeof(float) * angle_numbers);
    Malloc_Safely((void**)&h_angle_theta0, sizeof(float) * angle_numbers);
    Malloc_Safely((void**)&h_angle_ene, sizeof(float) * angle_numbers);
    memset(h_angle_ene, 0, sizeof(float) * angle_numbers);
    Malloc_Safely((void**)&h_sigma_of_angle_ene, sizeof(float));
    memset(h_sigma_of_angle_ene, 0, sizeof(float));
}

void ANGLE::Parameter_Host_To_Device()
{
    Device_Malloc_Safely((void**)&d_atom_a, sizeof(int) * angle_numbers);
    Device_Malloc_Safely((void**)&d_atom_b, sizeof(int) * angle_numbers);
    Device_Malloc_Safely((void**)&d_atom_c, sizeof(int) * angle_numbers);
    Device_Malloc_Safely((void**)&d_angle_k, sizeof(float) * angle_numbers);
    Device_Malloc_Safely((void**)&d_angle_theta0,
                         sizeof(float) * angle_numbers);
    Device_Malloc_Safely((void**)&d_angle_ene, sizeof(float) * angle_numbers);
    Device_Malloc_Safely((void**)&d_sigma_of_angle_ene, sizeof(float));

    Device_Malloc_Safely((void**)&d_atom_a_local, sizeof(int) * angle_numbers);
    Device_Malloc_Safely((void**)&d_atom_b_local, sizeof(int) * angle_numbers);
    Device_Malloc_Safely((void**)&d_atom_c_local, sizeof(int) * angle_numbers);
    Device_Malloc_Safely((void**)&d_angle_k_local,
                         sizeof(float) * angle_numbers);
    Device_Malloc_Safely((void**)&d_angle_theta0_local,
                         sizeof(float) * angle_numbers);
    Device_Malloc_Safely((void**)&d_global_index_local,
                         sizeof(int) * angle_numbers);
    Device_Malloc_Safely((void**)&d_num_angle_local, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_term, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_atom, sizeof(int));
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_invalid_geometry_term, sizeof(int));
#endif
    deviceMemset(d_num_angle_local, 0, sizeof(int));

    deviceMemcpy(d_atom_a, h_atom_a, sizeof(int) * angle_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_atom_b, h_atom_b, sizeof(int) * angle_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_atom_c, h_atom_c, sizeof(int) * angle_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_angle_k, h_angle_k, sizeof(float) * angle_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_angle_theta0, h_angle_theta0, sizeof(float) * angle_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemcpy(d_angle_ene, h_angle_ene, sizeof(float) * angle_numbers,
                 deviceMemcpyHostToDevice);
    deviceMemset(d_sigma_of_angle_ene, 0, sizeof(float));
}

static __global__ void get_local_device(
    int angle_numbers, int local_coordinate_numbers, const int* d_atom_a,
    const int* d_atom_b, const int* d_atom_c, const char* atom_local_label,
    const int* atom_local_id, int* d_atom_a_local, int* d_atom_b_local,
    int* d_atom_c_local, const float* d_angle_k, const float* d_angle_theta0,
    float* d_angle_k_local, float* d_angle_theta0_local,
    int* d_global_index_local, int* d_num_angle_local,
    int* d_invalid_local_term, int* d_invalid_local_atom)
{
#ifdef USE_GPU
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx != 0) return;
#endif
    d_num_angle_local[0] = 0;
    d_invalid_local_term[0] = -1;
    d_invalid_local_atom[0] = -1;
    for (int i = 0; i < angle_numbers; i++)
    {
        const int global_atoms[3] = {d_atom_a[i], d_atom_b[i], d_atom_c[i]};
        if (atom_local_label[global_atoms[0]] == 1 ||
            atom_local_label[global_atoms[1]] == 1 ||
            atom_local_label[global_atoms[2]] == 1)
        {
            int local_atoms[3];
            for (int atom = 0; atom < 3; atom++)
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
            const int local_term = d_num_angle_local[0];
            d_atom_a_local[local_term] = local_atoms[0];
            d_atom_b_local[local_term] = local_atoms[1];
            d_atom_c_local[local_term] = local_atoms[2];
            d_angle_k_local[local_term] = d_angle_k[i];
            d_angle_theta0_local[local_term] = d_angle_theta0[i];
            d_global_index_local[local_term] = i;
            d_num_angle_local[0]++;
        }
    }
}

void ANGLE::Get_Local(int* atom_local, int local_atom_numbers,
                      int ghost_numbers, char* atom_local_label,
                      int* atom_local_id)
{
    if (!is_initialized) return;
    (void)atom_local;
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > std::numeric_limits<int>::max() - ghost_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "ANGLE::Get_Local",
            "Reason:\n\t%s received invalid local/ghost atom counts %d/%d\n",
            module_name, local_atom_numbers, ghost_numbers);
        return;
    }
    const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
    num_angle_local = 0;
    this->local_atom_numbers = local_atom_numbers;
    Launch_Device_Kernel(get_local_device, 1, 1, 0, NULL, angle_numbers,
                         local_coordinate_numbers, d_atom_a, d_atom_b, d_atom_c,
                         atom_local_label, atom_local_id, d_atom_a_local,
                         d_atom_b_local, d_atom_c_local, d_angle_k,
                         d_angle_theta0, d_angle_k_local, d_angle_theta0_local,
                         d_global_index_local, d_num_angle_local,
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
            spongeErrorSimulationBreakDown, "ANGLE::Get_Local",
            "Reason:\n\t%s angle term %d (global atoms %d %d %d) maps "
            "global atom %d to local index %d outside the valid owned/ghost "
            "range [0, %d) on this domain\n",
            module_name, invalid_term, h_atom_a[invalid_term],
            h_atom_b[invalid_term], h_atom_c[invalid_term], invalid_atom,
            invalid_local_id, local_coordinate_numbers);
        return;
    }
    deviceMemcpy(&num_angle_local, d_num_angle_local, sizeof(int),
                 deviceMemcpyDeviceToHost);
}

void ANGLE::Angle_Force_With_Atom_Energy_And_Virial(
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* frc,
    int need_atom_energy, float* atom_energy, int need_virial,
    LTMatrix3* atom_virial_tensor)
{
    if (is_initialized && num_angle_local > 0)
    {
#ifndef GPU_ARCH_NAME
        deviceMemset(d_invalid_geometry_term, -1, sizeof(int));
#endif
        Launch_Device_Kernel(
            Angle_Force_With_Atom_Energy_And_Virial_Device,
            (num_angle_local + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, this->num_angle_local, crd,
            cell, rcell, this->local_atom_numbers, this->d_atom_a_local,
            this->d_atom_b_local, this->d_atom_c_local, this->d_angle_k_local,
            this->d_angle_theta0_local, frc, need_atom_energy, atom_energy,
            this->d_angle_ene, need_virial, atom_virial_tensor,
            this->d_global_index_local, this->d_atom_a, this->d_atom_b,
            this->d_atom_c, this->d_invalid_geometry_term);
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
                "ANGLE::Angle_Force_With_Atom_Energy_And_Virial",
                "Reason:\n\t%s angle term %d (global atoms %d %d %d) %s\n",
                module_name, invalid_term, h_atom_a[invalid_term],
                h_atom_b[invalid_term], h_atom_c[invalid_term],
                accumulator_failure
                    ? "would produce a non-finite force/energy/virial "
                      "accumulator"
                    : "has an undefined zero-arm/collinear geometry or a "
                      "non-finite/unrepresentable "
                      "geometry/energy/force/virial");
            return;
        }
#endif
    }
}

void ANGLE::Step_Print(CONTROLLER* controller)
{
    if (is_initialized && CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        Sum_Of_List(d_angle_ene, d_sigma_of_angle_ene,
                    num_angle_local);  // 修改为local求和
        deviceMemcpy(h_sigma_of_angle_ene, d_sigma_of_angle_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, h_sigma_of_angle_ene, 1, MPI_FLOAT, MPI_SUM,
                      CONTROLLER::pp_comm);
#endif
        if (!Float_Memory_Is_Finite(h_sigma_of_angle_ene))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "ANGLE::Step_Print",
                "Reason:\n\t%s angle total energy is non-finite after "
                "local/MPI reduction\n",
                module_name);
            return;
        }
        if (CONTROLLER::MPI_rank == 0)
        {
            controller->Step_Print(this->module_name, h_sigma_of_angle_ene,
                                   true);
        }
    }
}
