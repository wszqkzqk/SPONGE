#include "virtual_atoms.h"

#include "../xponge/ir/virtual_atoms.hpp"
#include "../xponge/load/native/virtual_atoms.hpp"
#include "../xponge/xponge.h"

struct CV_VIRTUAL_ATOM_DEFINITION
{
    std::string name;
    std::string type;
    int target = -1;
    int level = 0;
    std::vector<int> from;
    std::vector<float> weights;
};

static int Virtual_Atom_Block_Count(const int item_count)
{
    const int block_size = CONTROLLER::device_max_thread;
    return item_count / block_size + (item_count % block_size != 0);
}

static __device__ __forceinline__ bool Virtual_Atom_Float_Is_Finite(float value)
{
#ifdef GPU_ARCH_NAME
    // Do not use isfinite here.  CUDA/HIP fast-math compilation is allowed to
    // assume finite operands and can optimize that predicate away.
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
#else
    return Float_Memory_Is_Finite(&value);
#endif
}

static __device__ __forceinline__ bool Virtual_Atom_Vector_Is_Finite(
    const VECTOR& value)
{
    return Virtual_Atom_Float_Is_Finite(value.x) &&
           Virtual_Atom_Float_Is_Finite(value.y) &&
           Virtual_Atom_Float_Is_Finite(value.z);
}

static __global__ void v0_Coordinate_Refresh(const int virtual_numbers,
                                             const VIRTUAL_TYPE_0* v_info,
                                             VECTOR* crd, const LTMatrix3 cell,
                                             const LTMatrix3 rcell)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_0 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        float h_double = v_temp.h_double;
        VECTOR temp = crd[atom_1];
        temp.z = h_double - temp.z;
        crd[atom_v] = temp;
    }
}

static __global__ void v1_Coordinate_Refresh(const int virtual_numbers,
                                             const VIRTUAL_TYPE_1* v_info,
                                             VECTOR* crd, const LTMatrix3 cell,
                                             const LTMatrix3 rcell)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_1 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        int atom_2 = v_temp.from_2;
        float a = v_temp.a;
        VECTOR rv1 = a * Get_Periodic_Displacement(crd[atom_2], crd[atom_1],
                                                   cell, rcell);
        crd[atom_v] = crd[atom_1] + rv1;
    }
}

static __global__ void v2_Coordinate_Refresh(const int virtual_numbers,
                                             const VIRTUAL_TYPE_2* v_info,
                                             VECTOR* crd, const LTMatrix3 cell,
                                             const LTMatrix3 rcell)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_2 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        int atom_2 = v_temp.from_2;
        int atom_3 = v_temp.from_3;
        float a = v_temp.a;
        float b = v_temp.b;

        const VECTOR r1 = crd[atom_1];
        const VECTOR r2 = crd[atom_2];
        const VECTOR r3 = crd[atom_3];

        VECTOR rv1 = a * Get_Periodic_Displacement(r2, r1, cell, rcell) +
                     b * Get_Periodic_Displacement(r3, r1, cell, rcell);

        crd[atom_v] = crd[atom_1] + rv1;
    }
}

static __global__ void v3_Coordinate_Refresh(const int virtual_numbers,
                                             const VIRTUAL_TYPE_3* v_info,
                                             VECTOR* crd, const LTMatrix3 cell,
                                             const LTMatrix3 rcell,
                                             int* singularity)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_3 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        int atom_2 = v_temp.from_2;
        int atom_3 = v_temp.from_3;
        float d = v_temp.d;
        float k = v_temp.k;
        const VECTOR r1 = crd[atom_1];

        if (d == 0.0f)
        {
            if (!Virtual_Atom_Vector_Is_Finite(r1))
            {
                atomicExch(singularity, v_temp.global_virtual_atom);
            }
            crd[atom_v] = r1;
        }
        else
        {
            const VECTOR r2 = crd[atom_2];
            const VECTOR r3 = crd[atom_3];
            VECTOR r21 = Get_Periodic_Displacement(r2, r1, cell, rcell);
            VECTOR r32 = Get_Periodic_Displacement(r3, r2, cell, rcell);
            VECTOR direction = r21 + k * r32;
            float direction_squared = direction * direction;
            if (!Virtual_Atom_Vector_Is_Finite(r1) ||
                !Virtual_Atom_Vector_Is_Finite(r2) ||
                !Virtual_Atom_Vector_Is_Finite(r3) ||
                !Virtual_Atom_Vector_Is_Finite(r21) ||
                !Virtual_Atom_Vector_Is_Finite(r32) ||
                !Virtual_Atom_Vector_Is_Finite(direction) ||
                !Virtual_Atom_Float_Is_Finite(direction_squared) ||
                !(direction_squared > 0.0f))
            {
                atomicExch(singularity, v_temp.global_virtual_atom);
                crd[atom_v] = r1;
            }
            else
            {
                float inverse_direction =
                    rnorm3df(direction.x, direction.y, direction.z);
                VECTOR unit_direction = inverse_direction * direction;
                VECTOR displacement = d * unit_direction;
                VECTOR virtual_position = r1 + displacement;
                if (!Virtual_Atom_Float_Is_Finite(inverse_direction) ||
                    !Virtual_Atom_Vector_Is_Finite(unit_direction) ||
                    !Virtual_Atom_Vector_Is_Finite(displacement) ||
                    !Virtual_Atom_Vector_Is_Finite(virtual_position))
                {
                    atomicExch(singularity, v_temp.global_virtual_atom);
                    crd[atom_v] = r1;
                }
                else
                {
                    crd[atom_v] = virtual_position;
                }
            }
        }
    }
}

static __global__ void v5_Coordinate_Refresh(const int virtual_numbers,
                                             const VIRTUAL_TYPE_5* v_info,
                                             VECTOR* crd, const LTMatrix3 cell,
                                             const LTMatrix3 rcell,
                                             int* singularity)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_5 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_o = v_temp.from_1;
        int atom_h1 = v_temp.from_2;
        int atom_h2 = v_temp.from_3;
        float d = v_temp.d;
        const VECTOR r_o = crd[atom_o];

        if (d == 0.0f)
        {
            if (!Virtual_Atom_Vector_Is_Finite(r_o))
            {
                atomicExch(singularity, v_temp.global_virtual_atom);
            }
            crd[atom_v] = r_o;
        }
        else
        {
            VECTOR oh1 =
                Get_Periodic_Displacement(crd[atom_h1], r_o, cell, rcell);
            VECTOR oh2 =
                Get_Periodic_Displacement(crd[atom_h2], r_o, cell, rcell);
            float oh1_squared = oh1 * oh1;
            float oh2_squared = oh2 * oh2;
            if (!Virtual_Atom_Vector_Is_Finite(r_o) ||
                !Virtual_Atom_Vector_Is_Finite(oh1) ||
                !Virtual_Atom_Vector_Is_Finite(oh2) ||
                !Virtual_Atom_Float_Is_Finite(oh1_squared) ||
                !Virtual_Atom_Float_Is_Finite(oh2_squared) ||
                !(oh1_squared > 0.0f) || !(oh2_squared > 0.0f))
            {
                atomicExch(singularity, v_temp.global_virtual_atom);
                crd[atom_v] = r_o;
            }
            else
            {
                float inverse_oh1 = rnorm3df(oh1.x, oh1.y, oh1.z);
                float inverse_oh2 = rnorm3df(oh2.x, oh2.y, oh2.z);
                VECTOR unit_oh1 = inverse_oh1 * oh1;
                VECTOR unit_oh2 = inverse_oh2 * oh2;
                VECTOR bisector = unit_oh1 + unit_oh2;
                float bisector_squared = bisector * bisector;
                if (!Virtual_Atom_Float_Is_Finite(inverse_oh1) ||
                    !Virtual_Atom_Float_Is_Finite(inverse_oh2) ||
                    !Virtual_Atom_Vector_Is_Finite(unit_oh1) ||
                    !Virtual_Atom_Vector_Is_Finite(unit_oh2) ||
                    !Virtual_Atom_Vector_Is_Finite(bisector) ||
                    !Virtual_Atom_Float_Is_Finite(bisector_squared) ||
                    !(bisector_squared > 0.0f))
                {
                    atomicExch(singularity, v_temp.global_virtual_atom);
                    crd[atom_v] = r_o;
                }
                else
                {
                    float inverse_bisector =
                        rnorm3df(bisector.x, bisector.y, bisector.z);
                    VECTOR unit_bisector = inverse_bisector * bisector;
                    VECTOR displacement = d * unit_bisector;
                    VECTOR virtual_position = r_o + displacement;
                    if (!Virtual_Atom_Float_Is_Finite(inverse_bisector) ||
                        !Virtual_Atom_Vector_Is_Finite(unit_bisector) ||
                        !Virtual_Atom_Vector_Is_Finite(displacement) ||
                        !Virtual_Atom_Vector_Is_Finite(virtual_position))
                    {
                        atomicExch(singularity, v_temp.global_virtual_atom);
                        crd[atom_v] = r_o;
                    }
                    else
                    {
                        crd[atom_v] = virtual_position;
                    }
                }
            }
        }
    }
}

static __global__ void v4_Coordinate_Refresh(const int atom_numbers,
                                             const int virtual_atom,
                                             const int* from_atoms,
                                             const float* weight,
                                             VECTOR* coordinate)
{
    VECTOR new_position = {0, 0, 0};
#ifdef USE_GPU
    // One lane group handles one CV site.  Striding is required because a CV
    // center may contain more source atoms than the CUDA/HIP warp size.
    for (std::size_t i = static_cast<std::size_t>(threadIdx.x);
         i < static_cast<std::size_t>(atom_numbers);
         i += static_cast<std::size_t>(blockDim.x))
    {
        new_position = new_position + weight[i] * coordinate[from_atoms[i]];
    }
    for (int delta = warpSize >> 1; delta > 0; delta >>= 1)
    {
        new_position.x +=
            deviceShflDown(FULL_MASK, new_position.x, delta, warpSize);
        new_position.y +=
            deviceShflDown(FULL_MASK, new_position.y, delta, warpSize);
        new_position.z +=
            deviceShflDown(FULL_MASK, new_position.z, delta, warpSize);
    }
    if (threadIdx.x == 0)
    {
        coordinate[virtual_atom] = new_position;
    }
#else
    float px = 0.0f, py = 0.0f, pz = 0.0f;
#pragma omp parallel for reduction(+ : px, py, pz)
    for (int i = 0; i < atom_numbers; i++)
    {
        VECTOR p = weight[i] * coordinate[from_atoms[i]];
        px += p.x;
        py += p.y;
        pz += p.z;
    }
    new_position = {px, py, pz};
    coordinate[virtual_atom] = new_position;
#endif
}

static __global__ void v0_Force_Redistribute(
    const int virtual_numbers, const VIRTUAL_TYPE_0* v_info, const VECTOR* crd,
    const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* force)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_0 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        VECTOR force_v = force[atom_v];
        atomicAdd(&force[atom_1].x, force_v.x);
        atomicAdd(&force[atom_1].y, force_v.y);
        atomicAdd(&force[atom_1].z, -force_v.z);
        force_v.x = 0.0f;
        force_v.y = 0.0f;
        force_v.z = 0.0f;
        force[atom_v] = force_v;
    }
}

static __global__ void v1_Force_Redistribute(
    const int virtual_numbers, const VIRTUAL_TYPE_1* v_info, const VECTOR* crd,
    const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* force)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_1 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        int atom_2 = v_temp.from_2;
        float a = v_temp.a;
        VECTOR force_v = force[atom_v];
        atomicAdd(&force[atom_1].x, (1 - a) * force_v.x);
        atomicAdd(&force[atom_1].y, (1 - a) * force_v.y);
        atomicAdd(&force[atom_1].z, (1 - a) * force_v.z);

        atomicAdd(&force[atom_2].x, a * force_v.x);
        atomicAdd(&force[atom_2].y, a * force_v.y);
        atomicAdd(&force[atom_2].z, a * force_v.z);

        force_v.x = 0.0f;
        force_v.y = 0.0f;
        force_v.z = 0.0f;
        force[atom_v] = force_v;
    }
}

static __global__ void v2_Force_Redistribute(
    const int virtual_numbers, const VIRTUAL_TYPE_2* v_info, const VECTOR* crd,
    const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* force)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_2 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        int atom_2 = v_temp.from_2;
        int atom_3 = v_temp.from_3;
        float a = v_temp.a;
        float b = v_temp.b;
        VECTOR force_v = force[atom_v];
        atomicAdd(&force[atom_1].x, (1 - a - b) * force_v.x);
        atomicAdd(&force[atom_1].y, (1 - a - b) * force_v.y);
        atomicAdd(&force[atom_1].z, (1 - a - b) * force_v.z);

        atomicAdd(&force[atom_2].x, a * force_v.x);
        atomicAdd(&force[atom_2].y, a * force_v.y);
        atomicAdd(&force[atom_2].z, a * force_v.z);

        atomicAdd(&force[atom_3].x, b * force_v.x);
        atomicAdd(&force[atom_3].y, b * force_v.y);
        atomicAdd(&force[atom_3].z, b * force_v.z);

        force_v.x = 0.0f;
        force_v.y = 0.0f;
        force_v.z = 0.0f;
        force[atom_v] = force_v;
    }
}

static __global__ void v2_Force_Redistribute_No_Atomic(
    const int virtual_numbers, const VIRTUAL_TYPE_2* v_info, const VECTOR* crd,
    const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* force)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_2 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        int atom_2 = v_temp.from_2;
        int atom_3 = v_temp.from_3;
        float a = v_temp.a;
        float b = v_temp.b;
        VECTOR force_v = force[atom_v];

        force[atom_1].x += (1 - a - b) * force_v.x;
        force[atom_1].y += (1 - a - b) * force_v.y;
        force[atom_1].z += (1 - a - b) * force_v.z;

        force[atom_2].x += a * force_v.x;
        force[atom_2].y += a * force_v.y;
        force[atom_2].z += a * force_v.z;

        force[atom_3].x += b * force_v.x;
        force[atom_3].y += b * force_v.y;
        force[atom_3].z += b * force_v.z;

        force[atom_v] = {0, 0, 0};
    }
}

static __global__ void v3_Force_Redistribute(const int virtual_numbers,
                                             const VIRTUAL_TYPE_3* v_info,
                                             const VECTOR* crd,
                                             const LTMatrix3 cell,
                                             const LTMatrix3 rcell,
                                             VECTOR* force, int* singularity)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_3 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_1 = v_temp.from_1;
        int atom_2 = v_temp.from_2;
        int atom_3 = v_temp.from_3;
        float d = v_temp.d;
        float k = v_temp.k;
        VECTOR force_v = force[atom_v];

        bool valid = Virtual_Atom_Vector_Is_Finite(force_v);
        VECTOR force_1 = force_v;
        VECTOR force_2 = {0.0f, 0.0f, 0.0f};
        VECTOR force_3 = {0.0f, 0.0f, 0.0f};
        if (valid && d != 0.0f)
        {
            const VECTOR r1 = crd[atom_1];
            const VECTOR r2 = crd[atom_2];
            const VECTOR r3 = crd[atom_3];
            VECTOR r21 = Get_Periodic_Displacement(r2, r1, cell, rcell);
            VECTOR r32 = Get_Periodic_Displacement(r3, r2, cell, rcell);
            VECTOR direction = r21 + k * r32;
            float direction_squared = direction * direction;
            valid = Virtual_Atom_Vector_Is_Finite(r1) &&
                    Virtual_Atom_Vector_Is_Finite(r2) &&
                    Virtual_Atom_Vector_Is_Finite(r3) &&
                    Virtual_Atom_Vector_Is_Finite(r21) &&
                    Virtual_Atom_Vector_Is_Finite(r32) &&
                    Virtual_Atom_Vector_Is_Finite(direction) &&
                    Virtual_Atom_Float_Is_Finite(direction_squared) &&
                    direction_squared > 0.0f;
            if (valid)
            {
                float inverse_direction =
                    rnorm3df(direction.x, direction.y, direction.z);
                VECTOR unit_direction = inverse_direction * direction;
                float parallel_force = unit_direction * force_v;
                VECTOR perpendicular_force =
                    force_v - parallel_force * unit_direction;
                float scale = d * inverse_direction;
                VECTOR redistributed_force = scale * perpendicular_force;
                force_1 = force_v - redistributed_force;
                force_2 = (1.0f - k) * redistributed_force;
                force_3 = k * redistributed_force;
                valid = Virtual_Atom_Float_Is_Finite(inverse_direction) &&
                        Virtual_Atom_Vector_Is_Finite(unit_direction) &&
                        Virtual_Atom_Float_Is_Finite(parallel_force) &&
                        Virtual_Atom_Vector_Is_Finite(perpendicular_force) &&
                        Virtual_Atom_Float_Is_Finite(scale) &&
                        Virtual_Atom_Vector_Is_Finite(redistributed_force) &&
                        Virtual_Atom_Vector_Is_Finite(force_1) &&
                        Virtual_Atom_Vector_Is_Finite(force_2) &&
                        Virtual_Atom_Vector_Is_Finite(force_3);
            }
        }
        if (valid)
        {
            valid = Finite_Atomic_Add(&force[atom_1].x, force_1.x) &&
                    Finite_Atomic_Add(&force[atom_1].y, force_1.y) &&
                    Finite_Atomic_Add(&force[atom_1].z, force_1.z);
            if (valid && d != 0.0f)
            {
                valid = Finite_Atomic_Add(&force[atom_2].x, force_2.x) &&
                        Finite_Atomic_Add(&force[atom_2].y, force_2.y) &&
                        Finite_Atomic_Add(&force[atom_2].z, force_2.z) &&
                        Finite_Atomic_Add(&force[atom_3].x, force_3.x) &&
                        Finite_Atomic_Add(&force[atom_3].y, force_3.y) &&
                        Finite_Atomic_Add(&force[atom_3].z, force_3.z);
            }
        }
        if (!valid)
        {
            atomicExch(singularity, v_temp.global_virtual_atom);
        }
        force[atom_v] = {0.0f, 0.0f, 0.0f};
    }
}

static __global__ void v5_Force_Redistribute(const int virtual_numbers,
                                             const VIRTUAL_TYPE_5* v_info,
                                             const VECTOR* crd,
                                             const LTMatrix3 cell,
                                             const LTMatrix3 rcell,
                                             VECTOR* force, int* singularity)
{
#ifdef USE_GPU
    std::size_t i = static_cast<std::size_t>(threadIdx.x) +
                    static_cast<std::size_t>(blockIdx.x) * blockDim.x;
    if (i < static_cast<std::size_t>(virtual_numbers))
#else
#pragma omp parallel for
    for (int i = 0; i < virtual_numbers; i++)
#endif
    {
        VIRTUAL_TYPE_5 v_temp = v_info[i];
        int atom_v = v_temp.virtual_atom;
        int atom_o = v_temp.from_1;
        int atom_h1 = v_temp.from_2;
        int atom_h2 = v_temp.from_3;
        float d = v_temp.d;
        VECTOR force_v = force[atom_v];

        if (d == 0.0f)
        {
            if (!Virtual_Atom_Vector_Is_Finite(force_v))
            {
                atomicExch(singularity, v_temp.global_virtual_atom);
            }
            else
            {
                if (!Finite_Atomic_Add(&force[atom_o].x, force_v.x) ||
                    !Finite_Atomic_Add(&force[atom_o].y, force_v.y) ||
                    !Finite_Atomic_Add(&force[atom_o].z, force_v.z))
                {
                    atomicExch(singularity, v_temp.global_virtual_atom);
                }
            }
        }
        else
        {
            const VECTOR r_o = crd[atom_o];
            VECTOR oh1 =
                Get_Periodic_Displacement(crd[atom_h1], r_o, cell, rcell);
            VECTOR oh2 =
                Get_Periodic_Displacement(crd[atom_h2], r_o, cell, rcell);
            float oh1_squared = oh1 * oh1;
            float oh2_squared = oh2 * oh2;
            if (!Virtual_Atom_Vector_Is_Finite(force_v) ||
                !Virtual_Atom_Vector_Is_Finite(r_o) ||
                !Virtual_Atom_Vector_Is_Finite(oh1) ||
                !Virtual_Atom_Vector_Is_Finite(oh2) ||
                !Virtual_Atom_Float_Is_Finite(oh1_squared) ||
                !Virtual_Atom_Float_Is_Finite(oh2_squared) ||
                !(oh1_squared > 0.0f) || !(oh2_squared > 0.0f))
            {
                atomicExch(singularity, v_temp.global_virtual_atom);
            }
            else
            {
                float inverse_oh1 = rnorm3df(oh1.x, oh1.y, oh1.z);
                float inverse_oh2 = rnorm3df(oh2.x, oh2.y, oh2.z);
                VECTOR unit_oh1 = inverse_oh1 * oh1;
                VECTOR unit_oh2 = inverse_oh2 * oh2;
                VECTOR bisector = unit_oh1 + unit_oh2;
                float bisector_squared = bisector * bisector;
                if (!Virtual_Atom_Float_Is_Finite(inverse_oh1) ||
                    !Virtual_Atom_Float_Is_Finite(inverse_oh2) ||
                    !Virtual_Atom_Vector_Is_Finite(unit_oh1) ||
                    !Virtual_Atom_Vector_Is_Finite(unit_oh2) ||
                    !Virtual_Atom_Vector_Is_Finite(bisector) ||
                    !Virtual_Atom_Float_Is_Finite(bisector_squared) ||
                    !(bisector_squared > 0.0f))
                {
                    atomicExch(singularity, v_temp.global_virtual_atom);
                }
                else
                {
                    float inverse_bisector =
                        rnorm3df(bisector.x, bisector.y, bisector.z);
                    VECTOR unit_bisector = inverse_bisector * bisector;
                    float scale = d * inverse_bisector;
                    VECTOR q = scale * (force_v - (unit_bisector * force_v) *
                                                      unit_bisector);
                    VECTOR force_h1 =
                        inverse_oh1 * (q - (unit_oh1 * q) * unit_oh1);
                    VECTOR force_h2 =
                        inverse_oh2 * (q - (unit_oh2 * q) * unit_oh2);
                    VECTOR force_o = force_v - force_h1 - force_h2;
                    if (!Virtual_Atom_Float_Is_Finite(inverse_bisector) ||
                        !Virtual_Atom_Float_Is_Finite(scale) ||
                        !Virtual_Atom_Vector_Is_Finite(q) ||
                        !Virtual_Atom_Vector_Is_Finite(force_o) ||
                        !Virtual_Atom_Vector_Is_Finite(force_h1) ||
                        !Virtual_Atom_Vector_Is_Finite(force_h2))
                    {
                        atomicExch(singularity, v_temp.global_virtual_atom);
                    }
                    else
                    {
                        if (!Finite_Atomic_Add(&force[atom_o].x, force_o.x) ||
                            !Finite_Atomic_Add(&force[atom_o].y, force_o.y) ||
                            !Finite_Atomic_Add(&force[atom_o].z, force_o.z) ||
                            !Finite_Atomic_Add(&force[atom_h1].x, force_h1.x) ||
                            !Finite_Atomic_Add(&force[atom_h1].y, force_h1.y) ||
                            !Finite_Atomic_Add(&force[atom_h1].z, force_h1.z) ||
                            !Finite_Atomic_Add(&force[atom_h2].x, force_h2.x) ||
                            !Finite_Atomic_Add(&force[atom_h2].y, force_h2.y) ||
                            !Finite_Atomic_Add(&force[atom_h2].z, force_h2.z))
                        {
                            atomicExch(singularity, v_temp.global_virtual_atom);
                        }
                    }
                }
            }
        }

        force[atom_v] = {0.0f, 0.0f, 0.0f};
    }
}

static __global__ void v4_Force_Redistribute(const int atom_numbers,
                                             const int virtual_atom,
                                             const int* from_atoms,
                                             const float* weight, VECTOR* frc)
{
    const VECTOR new_force = frc[virtual_atom];
#ifdef USE_GPU
    for (std::size_t i = static_cast<std::size_t>(threadIdx.x);
         i < static_cast<std::size_t>(atom_numbers);
         i += static_cast<std::size_t>(blockDim.x))
#else
#pragma omp parallel for firstprivate(new_force)
    for (int i = 0; i < atom_numbers; i++)
#endif
    {
        float this_weight = weight[i];
        float* this_frc = &frc[from_atoms[i]].x;
        atomicAdd(this_frc, this_weight * new_force.x);
        atomicAdd(this_frc + 1, this_weight * new_force.y);
        atomicAdd(this_frc + 2, this_weight * new_force.z);
    }
#ifdef USE_GPU
    __syncthreads();
    if (threadIdx.x == 0)
#endif
    {
        frc[virtual_atom] = {0.0f, 0.0f, 0.0f};
    }
}

void VIRTUAL_INFORMATION::Initial(CONTROLLER* controller,
                                  COLLECTIVE_VARIABLE_CONTROLLER* cv_controller,
                                  int atom_numbers, int no_direct_vatom_numbers,
                                  CheckMap cv_vatom_name, float* h_mass,
                                  int* system_freedom, CONECT* connectivity,
                                  const char* module_name)
{
    this->controller = controller;
    has_local_layout = false;
    global_atom_numbers = atom_numbers;
    auto fail_initialization = [&](const std::string& reason)
    {
        std::string message =
            "Reason:\n\tinvalid virtual-atom initialization: " + reason + "\n";
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                       "VIRTUAL_INFORMATION::Initial",
                                       message.c_str());
    };
    auto checked_allocation_size = [&](std::size_t count,
                                       std::size_t element_size,
                                       const char* description)
    {
        if (element_size != 0 &&
            count > std::numeric_limits<std::size_t>::max() / element_size)
        {
            fail_initialization(std::string(description) +
                                " allocation size overflows size_t");
        }
        return count * element_size;
    };
    if (atom_numbers < 0)
    {
        fail_initialization("atom count is negative");
    }
    if (CONTROLLER::device_max_thread <= 0 || CONTROLLER::device_warp <= 0)
    {
        fail_initialization("device block or lane-group size is non-positive");
    }
    if (no_direct_vatom_numbers < 0)
    {
        fail_initialization("CV virtual-atom count is negative");
    }
    if (atom_numbers >
        std::numeric_limits<int>::max() - no_direct_vatom_numbers)
    {
        fail_initialization(
            "atom count plus CV virtual-atom count overflows int");
    }
    const int atom_numbers_with_cv = atom_numbers + no_direct_vatom_numbers;
    checked_allocation_size(static_cast<std::size_t>(atom_numbers_with_cv),
                            sizeof(int), "virtual-level array");
    const char* selected_module_name =
        module_name == NULL ? "virtual_atom" : module_name;
    if (strlen(selected_module_name) >= sizeof(this->module_name))
    {
        fail_initialization("module name exceeds the supported length");
    }
    strcpy(this->module_name, selected_module_name);
    const auto& system_virtual_atoms = Xponge::system.virtual_atoms.records;
    Xponge::VirtualAtoms local_virtual_atoms;
    const std::vector<Xponge::VirtualAtomRecord>* records_to_use = NULL;
    if (module_name == NULL)
    {
        records_to_use = &system_virtual_atoms;
    }
    else if (controller->Command_Exist(this->module_name, "in_file"))
    {
        Xponge::Native_Load_Virtual_Atoms(&local_virtual_atoms, controller,
                                          this->module_name);
        records_to_use = &local_virtual_atoms.records;
    }
    bool has_in_file = records_to_use != NULL && !records_to_use->empty();
    auto fail_cv_virtual_atom = [&](const std::string& reason)
    {
        std::string message =
            "Reason:\n\tinvalid CV virtual atom: " + reason + "\n";
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                       "VIRTUAL_INFORMATION::Initial",
                                       message.c_str());
    };
    if (cv_vatom_name.size() !=
        static_cast<std::size_t>(no_direct_vatom_numbers))
    {
        fail_cv_virtual_atom("definition count is inconsistent");
    }
    if (records_to_use != NULL &&
        (records_to_use->size() >
             static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
         records_to_use->size() > static_cast<std::size_t>(atom_numbers)))
    {
        fail_initialization(
            "force-field virtual-atom record count exceeds the supported atom "
            "or int range");
    }
    checked_allocation_size(records_to_use == NULL ? 0 : records_to_use->size(),
                            sizeof(Xponge::VirtualAtomRecord),
                            "virtual-atom records");
    checked_allocation_size(static_cast<std::size_t>(no_direct_vatom_numbers),
                            sizeof(CV_VIRTUAL_ATOM_DEFINITION),
                            "CV virtual-atom definitions");
    std::vector<CV_VIRTUAL_ATOM_DEFINITION> cv_definitions(
        static_cast<std::size_t>(no_direct_vatom_numbers));
    std::vector<int> cv_topological_order;
    cv_topological_order.reserve(
        static_cast<std::size_t>(no_direct_vatom_numbers));
    std::vector<bool> cv_index_seen(no_direct_vatom_numbers, false);
    for (const auto& item : cv_vatom_name)
    {
        if (item.second < 0 || item.second >= no_direct_vatom_numbers ||
            cv_index_seen[item.second])
        {
            fail_cv_virtual_atom(
                "definition indices are invalid or duplicated");
        }
        cv_index_seen[item.second] = true;
        CV_VIRTUAL_ATOM_DEFINITION& definition = cv_definitions[item.second];
        definition.name = item.first;
        definition.target = atom_numbers + item.second;
        definition.type =
            cv_controller->Command(item.first.c_str(), "vatom_type");
        if (definition.type != "center" && definition.type != "center_of_mass")
        {
            fail_cv_virtual_atom("'" + definition.name +
                                 "' has unsupported type '" + definition.type +
                                 "'");
        }
        definition.from =
            cv_controller->Ask_For_Indefinite_Length_Int_Parameter(
                item.first.c_str(), "atom");
        if (definition.from.empty())
        {
            fail_cv_virtual_atom("'" + definition.name +
                                 "' has no source atoms");
        }
        if (definition.from.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            fail_cv_virtual_atom("'" + definition.name +
                                 "' has too many source atoms");
        }
        checked_allocation_size(definition.from.size(), sizeof(int),
                                "CV virtual-atom source array");
        for (int source : definition.from)
        {
            if (source < 0 || source >= atom_numbers_with_cv)
            {
                fail_cv_virtual_atom(
                    "'" + definition.name + "' has source atom index " +
                    std::to_string(source) + " outside [0, " +
                    std::to_string(atom_numbers_with_cv) + ")");
            }
            if (source == definition.target)
            {
                fail_cv_virtual_atom("'" + definition.name +
                                     "' uses itself as a source");
            }
            if (definition.type == "center_of_mass" && source >= atom_numbers)
            {
                fail_cv_virtual_atom(
                    "'" + definition.name +
                    "' uses a massless CV virtual atom as a center-of-mass "
                    "source");
            }
        }
        if (definition.type == "center")
        {
            definition.weights =
                cv_controller->Ask_For_Indefinite_Length_Float_Parameter(
                    item.first.c_str(), "weight");
            if (definition.weights.size() != definition.from.size())
            {
                fail_cv_virtual_atom("'" + definition.name +
                                     "' has a weight/source count mismatch");
            }
            checked_allocation_size(definition.weights.size(), sizeof(float),
                                    "CV virtual-atom weight array");
            for (float weight : definition.weights)
            {
                if (!Xponge::Virtual_Atom_Parameter_Is_Finite(weight))
                {
                    fail_cv_virtual_atom("'" + definition.name +
                                         "' has a non-finite weight");
                }
            }
        }
        else
        {
            double total_mass = 0.0;
            for (int source : definition.from)
            {
                float mass = h_mass[source];
                if (!Xponge::Virtual_Atom_Parameter_Is_Finite(mass) ||
                    mass < 0.0f)
                {
                    fail_cv_virtual_atom("'" + definition.name +
                                         "' has a source with invalid mass");
                }
                total_mass += static_cast<double>(mass);
                if (!Double_Memory_Is_Finite(&total_mass))
                {
                    fail_cv_virtual_atom("'" + definition.name +
                                         "' has an overflowing total mass");
                }
            }
            if (total_mass <= 0.0f)
            {
                fail_cv_virtual_atom("'" + definition.name +
                                     "' has non-positive total mass");
            }
            for (int source : definition.from)
            {
                const double exact_weight =
                    static_cast<double>(h_mass[source]) / total_mass;
                const float stored_weight = static_cast<float>(exact_weight);
                if (!Xponge::Virtual_Atom_Parameter_Is_Finite(stored_weight) ||
                    (exact_weight != 0.0 && stored_weight == 0.0f))
                {
                    fail_cv_virtual_atom(
                        "'" + definition.name +
                        "' has a center-of-mass weight outside the supported "
                        "finite float range");
                }
                definition.weights.push_back(stored_weight);
            }
        }
    }

    checked_allocation_size(static_cast<std::size_t>(no_direct_vatom_numbers),
                            sizeof(int), "CV dependency indegree array");
    checked_allocation_size(static_cast<std::size_t>(no_direct_vatom_numbers),
                            sizeof(std::vector<int>),
                            "CV dependency consumer array");
    std::vector<int> cv_indegree(
        static_cast<std::size_t>(no_direct_vatom_numbers), 0);
    std::vector<std::vector<int>> cv_consumers(
        static_cast<std::size_t>(no_direct_vatom_numbers));
    for (int definition_index = 0; definition_index < no_direct_vatom_numbers;
         definition_index++)
    {
        std::vector<int> dependencies;
        for (int source : cv_definitions[definition_index].from)
        {
            if (source < atom_numbers)
            {
                continue;
            }
            int dependency = source - atom_numbers;
            if (std::find(dependencies.begin(), dependencies.end(),
                          dependency) != dependencies.end())
            {
                continue;
            }
            dependencies.push_back(dependency);
            if (cv_indegree[definition_index] ==
                std::numeric_limits<int>::max())
            {
                fail_cv_virtual_atom("dependency count overflows int");
            }
            cv_indegree[definition_index]++;
            cv_consumers[dependency].push_back(definition_index);
        }
    }
    for (int definition_index = 0; definition_index < no_direct_vatom_numbers;
         definition_index++)
    {
        if (cv_indegree[definition_index] == 0)
        {
            cv_topological_order.push_back(definition_index);
        }
    }
    for (std::size_t head = 0; head < cv_topological_order.size(); head++)
    {
        int definition_index = cv_topological_order[head];
        for (int consumer : cv_consumers[definition_index])
        {
            cv_indegree[consumer]--;
            if (cv_indegree[consumer] == 0)
            {
                cv_topological_order.push_back(consumer);
            }
        }
    }
    if (cv_topological_order.size() != cv_definitions.size())
    {
        for (int definition_index = 0;
             definition_index < no_direct_vatom_numbers; definition_index++)
        {
            if (cv_indegree[definition_index] > 0)
            {
                fail_cv_virtual_atom("dependency graph containing '" +
                                     cv_definitions[definition_index].name +
                                     "' has a cycle");
            }
        }
    }
    if (has_in_file || no_direct_vatom_numbers > 0)
    {
        controller->printf("START INITIALIZING VIRTUAL ATOM\n");
        const std::size_t virtual_level_bytes = checked_allocation_size(
            static_cast<std::size_t>(atom_numbers_with_cv), sizeof(int),
            "virtual-level array");
        Malloc_Safely((void**)&virtual_level, virtual_level_bytes);
        for (int i = 0; i < atom_numbers_with_cv; i++)
        {
            virtual_level[i] = 0;
        }

        int virtual_type;
        int virtual_atom;

        // 文件会从头到尾读三遍，分别确定每个原子的虚拟等级（因为可能存在坐标依赖于虚原子的虚原子，所以不得不如此做）
        // 第一遍确定虚拟原子的层级
        controller->printf("    Start reading virtual levels\n");
        if (has_in_file)
        {
            Xponge::VirtualAtomLayout layout;
            std::string validation_error;
            if (!Xponge::Validate_And_Build_Virtual_Atom_Layout(
                    *records_to_use, atom_numbers, &layout, &validation_error))
            {
                std::string reason =
                    "Reason:\n\tinvalid virtual-atom definition: " +
                    validation_error + "\n";
                controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                               "VIRTUAL_INFORMATION::Initial",
                                               reason.c_str());
            }
            for (int atom = 0; atom < atom_numbers; atom++)
            {
                virtual_level[atom] = layout.atom_levels[atom];
            }
            for (const auto& record : *records_to_use)
            {
                virtual_atom = record.virtual_atom;
                for (int source : record.from)
                {
                    connectivity[0][virtual_atom].insert(source);
                    connectivity[0][source].insert(virtual_atom);
                }
            }
        }
        // Add validated CV virtual atoms in dependency order.  A CV virtual
        // atom may depend on a force-field virtual atom or on an earlier CV
        // virtual atom, so levels must be propagated after both graphs are
        // known rather than inferred in map iteration order.
        for (int definition_index : cv_topological_order)
        {
            CV_VIRTUAL_ATOM_DEFINITION& definition =
                cv_definitions[definition_index];
            int source_level = 0;
            for (int source : definition.from)
            {
                source_level = std::max(source_level, virtual_level[source]);
            }
            if (source_level == std::numeric_limits<int>::max())
            {
                fail_cv_virtual_atom("dependency depth overflows int");
            }
            definition.level = source_level + 1;
            virtual_level[definition.target] = definition.level;
        }
        // 层级初始化
        max_level = 0;
        int total_virtual_atoms = 0;
        for (int i = 0; i < atom_numbers_with_cv; i++)
        {
            int vli = virtual_level[i];
            if (vli > 0)
            {
                if (total_virtual_atoms == std::numeric_limits<int>::max())
                {
                    fail_initialization("virtual-atom count overflows int");
                }
                total_virtual_atoms++;
            }
            if (vli > max_level)
            {
                for (int j = 0; j < vli - max_level; j++)
                {
                    VIRTUAL_LAYER_INFORMATION virtual_layer;
                    virtual_layer_info.push_back(virtual_layer);
                }
                max_level = vli;
            }
        }
        const int ff_virtual_atoms =
            total_virtual_atoms - no_direct_vatom_numbers;
        if (ff_virtual_atoms < 0 ||
            ff_virtual_atoms > std::numeric_limits<int>::max() / 3)
        {
            fail_initialization(
                "three times the force-field virtual-atom count overflows int");
        }
        const int freedom_reduction = 3 * ff_virtual_atoms;
        if (system_freedom[0] <
            std::numeric_limits<int>::min() + freedom_reduction)
        {
            fail_initialization(
                "virtual-atom freedom-of-motion adjustment overflows int");
        }
        system_freedom[0] -= freedom_reduction;
        controller->printf("        Virtual Atoms Max Level is %d\n",
                           max_level);
        controller->printf("        Virtual Atoms Number is %d\n",
                           total_virtual_atoms);
        controller->printf("            FF Virtual Atoms Number is %d\n",
                           ff_virtual_atoms);
        controller->printf("            CV Virtual Atoms Number is %d\n",
                           no_direct_vatom_numbers);
        controller->printf("    End reading virtual levels\n");
        // 第二遍确定虚拟原子每一层的个数
        controller->printf(
            "    Start reading virtual type numbers in different levels\n");
        if (has_in_file)
        {
            for (const auto& record : *records_to_use)
            {
                virtual_type = record.type;
                virtual_atom = record.virtual_atom;
                VIRTUAL_LAYER_INFORMATION* temp_vl =
                    &virtual_layer_info[virtual_level[virtual_atom] - 1];
                switch (virtual_type)
                {
                    case 0:
                        temp_vl->v0_info.virtual_numbers += 1;
                        break;
                    case 1:
                        temp_vl->v1_info.virtual_numbers += 1;
                        break;
                    case 2:
                        temp_vl->v2_info.virtual_numbers += 1;
                        break;
                    case 3:
                        temp_vl->v3_info.virtual_numbers += 1;
                        break;
                    case 5:
                        temp_vl->v5_info.virtual_numbers += 1;
                        break;
                    default:
                        break;
                }
            }
        }

        for (const CV_VIRTUAL_ATOM_DEFINITION& definition : cv_definitions)
        {
            VIRTUAL_LAYER_INFORMATION* temp_vl =
                &virtual_layer_info[definition.level - 1];
            temp_vl->v4_info.virtual_numbers += 1;
        }

        // 每层的每种虚拟原子初始化
        for (int layer = 0; layer < max_level; layer++)
        {
            controller->printf("        Virutual level %d:\n", layer);
            VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];
            if (temp_vl->v0_info.virtual_numbers > 0)
            {
                controller->printf(
                    "            Virtual type 0 atom numbers is %d\n",
                    temp_vl->v0_info.virtual_numbers);
                Malloc_Safely(
                    (void**)&temp_vl->v0_info.h_virtual_type_0,
                    checked_allocation_size(
                        static_cast<std::size_t>(
                            temp_vl->v0_info.virtual_numbers),
                        sizeof(VIRTUAL_TYPE_0), "type-0 virtual-atom array"));
            }
            if (temp_vl->v1_info.virtual_numbers > 0)
            {
                controller->printf(
                    "            Virtual type 1 atom numbers is %d\n",
                    temp_vl->v1_info.virtual_numbers);
                Malloc_Safely(
                    (void**)&temp_vl->v1_info.h_virtual_type_1,
                    checked_allocation_size(
                        static_cast<std::size_t>(
                            temp_vl->v1_info.virtual_numbers),
                        sizeof(VIRTUAL_TYPE_1), "type-1 virtual-atom array"));
            }
            if (temp_vl->v2_info.virtual_numbers > 0)
            {
                controller->printf(
                    "            Virtual type 2 atom numbers is %d\n",
                    temp_vl->v2_info.virtual_numbers);
                Malloc_Safely(
                    (void**)&temp_vl->v2_info.h_virtual_type_2,
                    checked_allocation_size(
                        static_cast<std::size_t>(
                            temp_vl->v2_info.virtual_numbers),
                        sizeof(VIRTUAL_TYPE_2), "type-2 virtual-atom array"));
            }
            if (temp_vl->v3_info.virtual_numbers > 0)
            {
                controller->printf(
                    "            Virtual type 3 atom numbers is %d\n",
                    temp_vl->v3_info.virtual_numbers);
                Malloc_Safely(
                    (void**)&temp_vl->v3_info.h_virtual_type_3,
                    checked_allocation_size(
                        static_cast<std::size_t>(
                            temp_vl->v3_info.virtual_numbers),
                        sizeof(VIRTUAL_TYPE_3), "type-3 virtual-atom array"));
            }
            if (temp_vl->v4_info.virtual_numbers > 0)
            {
                controller->printf(
                    "            Virtual type 4 atom numbers is %d\n",
                    temp_vl->v4_info.virtual_numbers);
                Malloc_Safely(
                    (void**)&temp_vl->v4_info.h_virtual_type_4,
                    checked_allocation_size(
                        static_cast<std::size_t>(
                            temp_vl->v4_info.virtual_numbers),
                        sizeof(VIRTUAL_TYPE_4), "type-4 virtual-atom array"));
            }
            if (temp_vl->v5_info.virtual_numbers > 0)
            {
                controller->printf(
                    "            Virtual type 5 atom numbers is %d\n",
                    temp_vl->v5_info.virtual_numbers);
                Malloc_Safely(
                    (void**)&temp_vl->v5_info.h_virtual_type_5,
                    checked_allocation_size(
                        static_cast<std::size_t>(
                            temp_vl->v5_info.virtual_numbers),
                        sizeof(VIRTUAL_TYPE_5), "type-5 virtual-atom array"));
            }
        }
        controller->printf(
            "    End reading virtual type numbers in different levels\n");
        // 第三遍将所有信息填入
        controller->printf(
            "    Start reading information for every virtual atom\n");
        if (has_in_file)
        {
            std::map<int, int> count0, count1, count2, count3, count5;
            for (int i = 0; i < virtual_layer_info.size(); i++)
            {
                count0[i] = 0;
                count1[i] = 0;
                count2[i] = 0;
                count3[i] = 0;
                count5[i] = 0;
            }
            for (const auto& record : *records_to_use)
            {
                virtual_type = record.type;
                virtual_atom = record.virtual_atom;
                int this_level = virtual_level[virtual_atom] - 1;
                VIRTUAL_LAYER_INFORMATION* temp_vl =
                    &virtual_layer_info[this_level];
                switch (virtual_type)
                {
                    case 0:
                        temp_vl->v0_info.h_virtual_type_0[count0[this_level]]
                            .virtual_atom = record.virtual_atom;
                        temp_vl->v0_info.h_virtual_type_0[count0[this_level]]
                            .from_1 = record.from[0];
                        temp_vl->v0_info.h_virtual_type_0[count0[this_level]]
                            .h_double = 2 * record.parameter[0];
                        count0[this_level]++;
                        break;

                    case 1:
                        temp_vl->v1_info.h_virtual_type_1[count1[this_level]]
                            .virtual_atom = record.virtual_atom;
                        temp_vl->v1_info.h_virtual_type_1[count1[this_level]]
                            .from_1 = record.from[0];
                        temp_vl->v1_info.h_virtual_type_1[count1[this_level]]
                            .from_2 = record.from[1];
                        temp_vl->v1_info.h_virtual_type_1[count1[this_level]]
                            .a = record.parameter[0];
                        count1[this_level]++;
                        break;

                    case 2:
                        temp_vl->v2_info.h_virtual_type_2[count2[this_level]]
                            .virtual_atom = record.virtual_atom;
                        temp_vl->v2_info.h_virtual_type_2[count2[this_level]]
                            .from_1 = record.from[0];
                        temp_vl->v2_info.h_virtual_type_2[count2[this_level]]
                            .from_2 = record.from[1];
                        temp_vl->v2_info.h_virtual_type_2[count2[this_level]]
                            .from_3 = record.from[2];
                        temp_vl->v2_info.h_virtual_type_2[count2[this_level]]
                            .a = record.parameter[0];
                        temp_vl->v2_info.h_virtual_type_2[count2[this_level]]
                            .b = record.parameter[1];
                        count2[this_level]++;
                        break;

                    case 3:
                        temp_vl->v3_info.h_virtual_type_3[count3[this_level]]
                            .virtual_atom = record.virtual_atom;
                        temp_vl->v3_info.h_virtual_type_3[count3[this_level]]
                            .global_virtual_atom = record.virtual_atom;
                        temp_vl->v3_info.h_virtual_type_3[count3[this_level]]
                            .from_1 = record.from[0];
                        temp_vl->v3_info.h_virtual_type_3[count3[this_level]]
                            .from_2 = record.from[1];
                        temp_vl->v3_info.h_virtual_type_3[count3[this_level]]
                            .from_3 = record.from[2];
                        temp_vl->v3_info.h_virtual_type_3[count3[this_level]]
                            .d = record.parameter[0];
                        temp_vl->v3_info.h_virtual_type_3[count3[this_level]]
                            .k = record.parameter[1];
                        count3[this_level]++;
                        break;

                    case 5:
                        temp_vl->v5_info.h_virtual_type_5[count5[this_level]]
                            .virtual_atom = record.virtual_atom;
                        temp_vl->v5_info.h_virtual_type_5[count5[this_level]]
                            .global_virtual_atom = record.virtual_atom;
                        temp_vl->v5_info.h_virtual_type_5[count5[this_level]]
                            .from_1 = record.from[0];
                        temp_vl->v5_info.h_virtual_type_5[count5[this_level]]
                            .from_2 = record.from[1];
                        temp_vl->v5_info.h_virtual_type_5[count5[this_level]]
                            .from_3 = record.from[2];
                        temp_vl->v5_info.h_virtual_type_5[count5[this_level]]
                            .d = record.parameter[0];
                        count5[this_level]++;
                        break;

                    default:
                        break;
                }
            }
        }
        std::map<int, int> count4;
        for (int i = 0; i < virtual_layer_info.size(); i++)
        {
            count4[i] = 0;
        }
        for (const CV_VIRTUAL_ATOM_DEFINITION& definition : cv_definitions)
        {
            int this_level = definition.level - 1;
            VIRTUAL_LAYER_INFORMATION* temp_vl =
                &virtual_layer_info[this_level];
            VIRTUAL_TYPE_4& record =
                temp_vl->v4_info.h_virtual_type_4[count4[this_level]];
            record.virtual_atom = definition.target;
            record.atom_numbers = static_cast<int>(definition.from.size());
            const std::size_t source_bytes =
                checked_allocation_size(definition.from.size(), sizeof(int),
                                        "CV virtual-atom source array");
            const std::size_t weight_bytes = checked_allocation_size(
                definition.weights.size(), sizeof(float),
                "CV virtual-atom weight array");
            Malloc_Safely((void**)&record.h_from, source_bytes);
            memcpy(record.h_from, definition.from.data(), source_bytes);
            Device_Malloc_And_Copy_Safely((void**)&record.d_from, record.h_from,
                                          source_bytes);
            Malloc_Safely((void**)&record.h_weight, weight_bytes);
            memcpy(record.h_weight, definition.weights.data(), weight_bytes);
            Device_Malloc_And_Copy_Safely((void**)&record.d_weight,
                                          record.h_weight, weight_bytes);
            count4[this_level]++;
        }
        // 每层的数据信息传到cuda上去
        Device_Malloc_Safely((void**)&d_invalid_local_layout, 6 * sizeof(int));
        for (int layer = 0; layer < max_level; layer++)
        {
            VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];
            if (temp_vl->v0_info.virtual_numbers > 0)
            {
                const std::size_t bytes = checked_allocation_size(
                    static_cast<std::size_t>(temp_vl->v0_info.virtual_numbers),
                    sizeof(VIRTUAL_TYPE_0), "type-0 device/local array");
                Device_Malloc_And_Copy_Safely(
                    (void**)&temp_vl->v0_info.d_virtual_type_0,
                    temp_vl->v0_info.h_virtual_type_0, bytes);
                Device_Malloc_Safely((void**)&temp_vl->v0_info.l_virtual_type_0,
                                     bytes);
                Device_Malloc_Safely((void**)&temp_vl->v0_info.d_local_numbers,
                                     sizeof(int));
            }
            if (temp_vl->v1_info.virtual_numbers > 0)
            {
                const std::size_t bytes = checked_allocation_size(
                    static_cast<std::size_t>(temp_vl->v1_info.virtual_numbers),
                    sizeof(VIRTUAL_TYPE_1), "type-1 device/local array");
                Device_Malloc_And_Copy_Safely(
                    (void**)&temp_vl->v1_info.d_virtual_type_1,
                    temp_vl->v1_info.h_virtual_type_1, bytes);
                Device_Malloc_Safely((void**)&temp_vl->v1_info.l_virtual_type_1,
                                     bytes);
                Device_Malloc_Safely((void**)&temp_vl->v1_info.d_local_numbers,
                                     sizeof(int));
            }
            if (temp_vl->v2_info.virtual_numbers > 0)
            {
                const std::size_t bytes = checked_allocation_size(
                    static_cast<std::size_t>(temp_vl->v2_info.virtual_numbers),
                    sizeof(VIRTUAL_TYPE_2), "type-2 device/local array");
                Device_Malloc_And_Copy_Safely(
                    (void**)&temp_vl->v2_info.d_virtual_type_2,
                    temp_vl->v2_info.h_virtual_type_2, bytes);
                Device_Malloc_Safely((void**)&temp_vl->v2_info.l_virtual_type_2,
                                     bytes);
                Device_Malloc_Safely((void**)&temp_vl->v2_info.d_local_numbers,
                                     sizeof(int));
            }
            if (temp_vl->v3_info.virtual_numbers > 0)
            {
                const std::size_t bytes = checked_allocation_size(
                    static_cast<std::size_t>(temp_vl->v3_info.virtual_numbers),
                    sizeof(VIRTUAL_TYPE_3), "type-3 device/local array");
                if (d_type3_singularity == NULL)
                {
                    Device_Malloc_Safely((void**)&d_type3_singularity,
                                         sizeof(int));
                }
                Device_Malloc_And_Copy_Safely(
                    (void**)&temp_vl->v3_info.d_virtual_type_3,
                    temp_vl->v3_info.h_virtual_type_3, bytes);
                Device_Malloc_Safely((void**)&temp_vl->v3_info.l_virtual_type_3,
                                     bytes);
                Device_Malloc_Safely((void**)&temp_vl->v3_info.d_local_numbers,
                                     sizeof(int));
            }
            if (temp_vl->v5_info.virtual_numbers > 0)
            {
                const std::size_t bytes = checked_allocation_size(
                    static_cast<std::size_t>(temp_vl->v5_info.virtual_numbers),
                    sizeof(VIRTUAL_TYPE_5), "type-5 device/local array");
                if (d_type5_singularity == NULL)
                {
                    Device_Malloc_Safely((void**)&d_type5_singularity,
                                         sizeof(int));
                }
                Device_Malloc_And_Copy_Safely(
                    (void**)&temp_vl->v5_info.d_virtual_type_5,
                    temp_vl->v5_info.h_virtual_type_5, bytes);
                Device_Malloc_Safely((void**)&temp_vl->v5_info.l_virtual_type_5,
                                     bytes);
                Device_Malloc_Safely((void**)&temp_vl->v5_info.d_local_numbers,
                                     sizeof(int));
            }
        }
        controller->printf(
            "    End reading information for every virtual atom\n");

        is_initialized = 1;
        if (is_initialized && !is_controller_printf_initialized)
        {
            is_controller_printf_initialized = 1;
            controller->printf("    structure last modify date is %d\n",
                               last_modify_date);
        }

        for (int layer = 0; layer < max_level; layer++)
        {
            std::vector<int> mark(atom_numbers, -1);
            VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];
            VIRTUAL_TYPE_2* v_info = temp_vl->v2_info.h_virtual_type_2;
            int virtual_numbers = temp_vl->v2_info.virtual_numbers;
            for (int i = 0; i < virtual_numbers; ++i)
            {
                for (auto x :
                     {v_info[i].from_1, v_info[i].from_2, v_info[i].from_3})
                {
                    if (mark[x] < 0)
                    {
                        mark[x] = i;
                    }
                    else if (mark[x] != i)
                    {
                        need_atomic = true;
                    }
                }
            }
        }
        controller->printf("END INITIALIZING VIRTUAL ATOM\n\n");
    }
    else
    {
        controller->printf("VIRTUAL ATOM IS NOT INITIALIZED\n\n");
    }
}

void VIRTUAL_INFORMATION::Coordinate_Refresh(VECTOR* crd, const LTMatrix3 cell,
                                             const LTMatrix3 rcell)
{
    if (is_initialized)
    {
        if (d_type3_singularity != NULL)
        {
            int no_singularity = -1;
            deviceMemcpy(d_type3_singularity, &no_singularity, sizeof(int),
                         deviceMemcpyHostToDevice);
        }
        if (d_type5_singularity != NULL)
        {
            int no_singularity = -1;
            deviceMemcpy(d_type5_singularity, &no_singularity, sizeof(int),
                         deviceMemcpyHostToDevice);
        }
        // 每层之间需要串行计算，层内并行计算
        for (int layer = 0; layer < max_level; layer++)
        {
            VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];
            int v0_numbers = has_local_layout
                                 ? temp_vl->v0_info.local_numbers
                                 : temp_vl->v0_info.virtual_numbers;
            VIRTUAL_TYPE_0* v0_info = has_local_layout
                                          ? temp_vl->v0_info.l_virtual_type_0
                                          : temp_vl->v0_info.d_virtual_type_0;
            if (v0_numbers > 0)
                Launch_Device_Kernel(v0_Coordinate_Refresh,
                                     Virtual_Atom_Block_Count(v0_numbers),
                                     CONTROLLER::device_max_thread, 0, NULL,
                                     v0_numbers, v0_info, crd, cell, rcell);

            int v1_numbers = has_local_layout
                                 ? temp_vl->v1_info.local_numbers
                                 : temp_vl->v1_info.virtual_numbers;
            VIRTUAL_TYPE_1* v1_info = has_local_layout
                                          ? temp_vl->v1_info.l_virtual_type_1
                                          : temp_vl->v1_info.d_virtual_type_1;
            if (v1_numbers > 0)
                Launch_Device_Kernel(v1_Coordinate_Refresh,
                                     Virtual_Atom_Block_Count(v1_numbers),
                                     CONTROLLER::device_max_thread, 0, NULL,
                                     v1_numbers, v1_info, crd, cell, rcell);

            int v2_numbers = has_local_layout
                                 ? temp_vl->v2_info.local_numbers
                                 : temp_vl->v2_info.virtual_numbers;
            VIRTUAL_TYPE_2* v2_info = has_local_layout
                                          ? temp_vl->v2_info.l_virtual_type_2
                                          : temp_vl->v2_info.d_virtual_type_2;
            if (v2_numbers > 0)
                Launch_Device_Kernel(v2_Coordinate_Refresh,
                                     Virtual_Atom_Block_Count(v2_numbers),
                                     CONTROLLER::device_max_thread, 0, NULL,
                                     v2_numbers, v2_info, crd, cell, rcell);

            int v3_numbers = has_local_layout
                                 ? temp_vl->v3_info.local_numbers
                                 : temp_vl->v3_info.virtual_numbers;
            VIRTUAL_TYPE_3* v3_info = has_local_layout
                                          ? temp_vl->v3_info.l_virtual_type_3
                                          : temp_vl->v3_info.d_virtual_type_3;
            if (v3_numbers > 0)
                Launch_Device_Kernel(
                    v3_Coordinate_Refresh, Virtual_Atom_Block_Count(v3_numbers),
                    CONTROLLER::device_max_thread, 0, NULL, v3_numbers, v3_info,
                    crd, cell, rcell, d_type3_singularity);

            int v5_numbers = has_local_layout
                                 ? temp_vl->v5_info.local_numbers
                                 : temp_vl->v5_info.virtual_numbers;
            VIRTUAL_TYPE_5* v5_info = has_local_layout
                                          ? temp_vl->v5_info.l_virtual_type_5
                                          : temp_vl->v5_info.d_virtual_type_5;
            if (v5_numbers > 0)
                Launch_Device_Kernel(
                    v5_Coordinate_Refresh, Virtual_Atom_Block_Count(v5_numbers),
                    CONTROLLER::device_max_thread, 0, NULL, v5_numbers, v5_info,
                    crd, cell, rcell, d_type5_singularity);
        }
        if (d_type3_singularity != NULL)
        {
            int singular_atom = -1;
            deviceMemcpy(&singular_atom, d_type3_singularity, sizeof(int),
                         deviceMemcpyDeviceToHost);
            if (singular_atom >= 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "VIRTUAL_INFORMATION::Coordinate_Refresh",
                    "Reason:\n\ttype-3 virtual atom %d has non-finite or "
                    "unrepresentable data, or zero construction direction "
                    "while its distance is nonzero\n",
                    singular_atom);
            }
        }
        if (d_type5_singularity != NULL)
        {
            int singular_atom = -1;
            deviceMemcpy(&singular_atom, d_type5_singularity, sizeof(int),
                         deviceMemcpyDeviceToHost);
            if (singular_atom >= 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "VIRTUAL_INFORMATION::Coordinate_Refresh",
                    "Reason:\n\ttype-5 virtual atom %d has non-finite or "
                    "unrepresentable data, or singular O-H or bisector "
                    "geometry\n",
                    singular_atom);
            }
        }
    }
}

void VIRTUAL_INFORMATION::Force_Redistribute(const VECTOR* crd,
                                             const LTMatrix3 cell,
                                             const LTMatrix3 rcell, VECTOR* frc)
{
    if (is_initialized)
    {
        if (d_type3_singularity != NULL)
        {
            int no_singularity = -1;
            deviceMemcpy(d_type3_singularity, &no_singularity, sizeof(int),
                         deviceMemcpyHostToDevice);
        }
        if (d_type5_singularity != NULL)
        {
            int no_singularity = -1;
            deviceMemcpy(d_type5_singularity, &no_singularity, sizeof(int),
                         deviceMemcpyHostToDevice);
        }
        // 每层之间需要串行逆向计算，层内并行计算
        for (int layer = max_level - 1; layer >= 0; layer--)
        {
            VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];
            if (temp_vl->v0_info.local_numbers > 0)
            {
                Launch_Device_Kernel(
                    v0_Force_Redistribute,
                    Virtual_Atom_Block_Count(temp_vl->v0_info.local_numbers),
                    CONTROLLER::device_max_thread, 0, NULL,
                    temp_vl->v0_info.local_numbers,
                    temp_vl->v0_info.l_virtual_type_0, crd, cell, rcell, frc);
            }
            if (temp_vl->v1_info.local_numbers > 0)
            {
                Launch_Device_Kernel(
                    v1_Force_Redistribute,
                    Virtual_Atom_Block_Count(temp_vl->v1_info.local_numbers),
                    CONTROLLER::device_max_thread, 0, NULL,
                    temp_vl->v1_info.local_numbers,
                    temp_vl->v1_info.l_virtual_type_1, crd, cell, rcell, frc);
            }
            if (temp_vl->v3_info.local_numbers > 0)
            {
                Launch_Device_Kernel(
                    v3_Force_Redistribute,
                    Virtual_Atom_Block_Count(temp_vl->v3_info.local_numbers),
                    CONTROLLER::device_max_thread, 0, NULL,
                    temp_vl->v3_info.local_numbers,
                    temp_vl->v3_info.l_virtual_type_3, crd, cell, rcell, frc,
                    d_type3_singularity);
            }
            if (temp_vl->v5_info.local_numbers > 0)
            {
                Launch_Device_Kernel(
                    v5_Force_Redistribute,
                    Virtual_Atom_Block_Count(temp_vl->v5_info.local_numbers),
                    CONTROLLER::device_max_thread, 0, NULL,
                    temp_vl->v5_info.local_numbers,
                    temp_vl->v5_info.l_virtual_type_5, crd, cell, rcell, frc,
                    d_type5_singularity);
            }

            if (temp_vl->v2_info.local_numbers > 0)
            {
                if (need_atomic)
                {
                    Launch_Device_Kernel(v2_Force_Redistribute,
                                         Virtual_Atom_Block_Count(
                                             temp_vl->v2_info.local_numbers),
                                         CONTROLLER::device_max_thread, 0, NULL,
                                         temp_vl->v2_info.local_numbers,
                                         temp_vl->v2_info.l_virtual_type_2, crd,
                                         cell, rcell, frc);
                }
                else
                {
                    Launch_Device_Kernel(v2_Force_Redistribute_No_Atomic,
                                         Virtual_Atom_Block_Count(
                                             temp_vl->v2_info.local_numbers),
                                         CONTROLLER::device_max_thread, 0, NULL,
                                         temp_vl->v2_info.local_numbers,
                                         temp_vl->v2_info.l_virtual_type_2, crd,
                                         cell, rcell, frc);
                }
            }
        }
        if (d_type3_singularity != NULL)
        {
            int singular_atom = -1;
            deviceMemcpy(&singular_atom, d_type3_singularity, sizeof(int),
                         deviceMemcpyDeviceToHost);
            if (singular_atom >= 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "VIRTUAL_INFORMATION::Force_Redistribute",
                    "Reason:\n\ttype-3 virtual atom %d has non-finite or "
                    "unrepresentable geometry/force data, a singular "
                    "construction direction, or an overflowing force "
                    "accumulator\n",
                    singular_atom);
            }
        }
        if (d_type5_singularity != NULL)
        {
            int singular_atom = -1;
            deviceMemcpy(&singular_atom, d_type5_singularity, sizeof(int),
                         deviceMemcpyDeviceToHost);
            if (singular_atom >= 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "VIRTUAL_INFORMATION::Force_Redistribute",
                    "Reason:\n\ttype-5 virtual atom %d has non-finite or "
                    "unrepresentable data, or singular O-H or bisector "
                    "geometry\n",
                    singular_atom);
            }
        }
    }
}

void VIRTUAL_INFORMATION::Coordinate_Refresh_CV(VECTOR* crd,
                                                const LTMatrix3 cell,
                                                const LTMatrix3 rcell)
{
    if (is_initialized)
    {
        // 每层之间需要串行计算，层内并行计算
        for (int layer = 0; layer < max_level; layer++)
        {
            VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];
            // 预留v4质心接口
            VIRTUAL_TYPE_4* temp_vl4;
            for (int iv4 = 0; iv4 < temp_vl->v4_info.virtual_numbers; iv4++)
            {
                temp_vl4 = temp_vl->v4_info.h_virtual_type_4 + iv4;
                Launch_Device_Kernel(
                    v4_Coordinate_Refresh, 1, CONTROLLER::device_warp, 0, NULL,
                    temp_vl4->atom_numbers, temp_vl4->virtual_atom,
                    temp_vl4->d_from, temp_vl4->d_weight, crd);
            }
        }
    }
}

void VIRTUAL_INFORMATION::Force_Redistribute_CV(const VECTOR* crd,
                                                const LTMatrix3 cell,
                                                const LTMatrix3 rcell,
                                                VECTOR* frc)
{
    if (is_initialized)
    {
        // 每层之间需要串行逆向计算，层内并行计算
        for (int layer = max_level - 1; layer >= 0; layer--)
        {
            VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];
            // 预留v4质心接口
            VIRTUAL_TYPE_4* temp_vl4;
            for (int iv4 = 0; iv4 < temp_vl->v4_info.virtual_numbers; iv4++)
            {
                temp_vl4 = temp_vl->v4_info.h_virtual_type_4 + iv4;
                Launch_Device_Kernel(
                    v4_Force_Redistribute, 1, CONTROLLER::device_warp, 0, NULL,
                    temp_vl4->atom_numbers, temp_vl4->virtual_atom,
                    temp_vl4->d_from, temp_vl4->d_weight, frc);
            }
        }
    }
}

static __device__ __forceinline__ bool Virtual_Atom_Get_Owned_Local_Id(
    const int global_atom, const int local_atom_numbers,
    const int* atom_local_id, const char* atom_local_label,
    const int virtual_type, const int virtual_atom, const int source_slot,
    int* local_atom, int* invalid_local_layout)
{
    const int candidate = atom_local_id[global_atom];
    const int label = static_cast<int>(atom_local_label[global_atom]);
    if (label != 1 || candidate < 0 || candidate >= local_atom_numbers)
    {
        if (invalid_local_layout[0] < 0)
        {
            invalid_local_layout[0] = virtual_type;
            invalid_local_layout[1] = virtual_atom;
            invalid_local_layout[2] = global_atom;
            invalid_local_layout[3] = candidate;
            invalid_local_layout[4] = label;
            invalid_local_layout[5] = source_slot;
        }
        return false;
    }
    *local_atom = candidate;
    return true;
}

static __global__ void get_local_device_V0(
    int virtual_numbers, int* local_numbers,
    const VIRTUAL_TYPE_0* d_virtual_type_0, VIRTUAL_TYPE_0* l_virtual_type_0,
    const int* atom_local_id, const char* atom_local_label,
    int local_atom_numbers, int* invalid_local_layout)
{
    local_numbers[0] = 0;
    for (int cluster = 0; cluster < virtual_numbers; cluster++)
    {
        int vatom = d_virtual_type_0[cluster].virtual_atom;
        int from1 = d_virtual_type_0[cluster].from_1;
        if (atom_local_label[vatom] == 1)
        {
            int local_vatom = -1, local_from1 = -1;
            bool valid = Virtual_Atom_Get_Owned_Local_Id(
                vatom, local_atom_numbers, atom_local_id, atom_local_label, 0,
                vatom, -1, &local_vatom, invalid_local_layout);
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from1, local_atom_numbers, atom_local_id, atom_local_label,
                    0, vatom, 0, &local_from1, invalid_local_layout) &&
                valid;
            if (valid)
            {
                l_virtual_type_0[local_numbers[0]] = d_virtual_type_0[cluster];
                l_virtual_type_0[local_numbers[0]].virtual_atom = local_vatom;
                l_virtual_type_0[local_numbers[0]].from_1 = local_from1;
                local_numbers[0] += 1;
            }
        }
    }
}

static __global__ void get_local_device_V1(
    int virtual_numbers, int* local_numbers,
    const VIRTUAL_TYPE_1* d_virtual_type_1, VIRTUAL_TYPE_1* l_virtual_type_1,
    const int* atom_local_id, const char* atom_local_label,
    int local_atom_numbers, int* invalid_local_layout)
{
    local_numbers[0] = 0;
    for (int cluster = 0; cluster < virtual_numbers; cluster++)
    {
        int vatom = d_virtual_type_1[cluster].virtual_atom;
        int from1 = d_virtual_type_1[cluster].from_1;
        int from2 = d_virtual_type_1[cluster].from_2;
        if (atom_local_label[vatom] == 1)
        {
            int local_vatom = -1, local_from1 = -1, local_from2 = -1;
            bool valid = Virtual_Atom_Get_Owned_Local_Id(
                vatom, local_atom_numbers, atom_local_id, atom_local_label, 1,
                vatom, -1, &local_vatom, invalid_local_layout);
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from1, local_atom_numbers, atom_local_id, atom_local_label,
                    1, vatom, 0, &local_from1, invalid_local_layout) &&
                valid;
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from2, local_atom_numbers, atom_local_id, atom_local_label,
                    1, vatom, 1, &local_from2, invalid_local_layout) &&
                valid;
            if (valid)
            {
                l_virtual_type_1[local_numbers[0]] = d_virtual_type_1[cluster];
                l_virtual_type_1[local_numbers[0]].virtual_atom = local_vatom;
                l_virtual_type_1[local_numbers[0]].from_1 = local_from1;
                l_virtual_type_1[local_numbers[0]].from_2 = local_from2;
                local_numbers[0] += 1;
            }
        }
    }
}

static __global__ void get_local_device_V2(
    int virtual_numbers, int* local_numbers,
    const VIRTUAL_TYPE_2* d_virtual_type_2, VIRTUAL_TYPE_2* l_virtual_type_2,
    const int* atom_local_id, const char* atom_local_label,
    int local_atom_numbers, int* invalid_local_layout)
{
    local_numbers[0] = 0;
    for (int cluster = 0; cluster < virtual_numbers; cluster++)
    {
        int vatom = d_virtual_type_2[cluster].virtual_atom;
        int from1 = d_virtual_type_2[cluster].from_1;
        int from2 = d_virtual_type_2[cluster].from_2;
        int from3 = d_virtual_type_2[cluster].from_3;
        if (atom_local_label[vatom] == 1)
        {
            int local_vatom = -1, local_from1 = -1, local_from2 = -1,
                local_from3 = -1;
            bool valid = Virtual_Atom_Get_Owned_Local_Id(
                vatom, local_atom_numbers, atom_local_id, atom_local_label, 2,
                vatom, -1, &local_vatom, invalid_local_layout);
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from1, local_atom_numbers, atom_local_id, atom_local_label,
                    2, vatom, 0, &local_from1, invalid_local_layout) &&
                valid;
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from2, local_atom_numbers, atom_local_id, atom_local_label,
                    2, vatom, 1, &local_from2, invalid_local_layout) &&
                valid;
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from3, local_atom_numbers, atom_local_id, atom_local_label,
                    2, vatom, 2, &local_from3, invalid_local_layout) &&
                valid;
            if (valid)
            {
                l_virtual_type_2[local_numbers[0]] = d_virtual_type_2[cluster];
                l_virtual_type_2[local_numbers[0]].virtual_atom = local_vatom;
                l_virtual_type_2[local_numbers[0]].from_1 = local_from1;
                l_virtual_type_2[local_numbers[0]].from_2 = local_from2;
                l_virtual_type_2[local_numbers[0]].from_3 = local_from3;
                local_numbers[0] += 1;
            }
        }
    }
}

static __global__ void get_local_device_V3(
    int virtual_numbers, int* local_numbers,
    const VIRTUAL_TYPE_3* d_virtual_type_3, VIRTUAL_TYPE_3* l_virtual_type_3,
    const int* atom_local_id, const char* atom_local_label,
    int local_atom_numbers, int* invalid_local_layout)
{
    local_numbers[0] = 0;
    for (int cluster = 0; cluster < virtual_numbers; cluster++)
    {
        int vatom = d_virtual_type_3[cluster].virtual_atom;
        int from1 = d_virtual_type_3[cluster].from_1;
        int from2 = d_virtual_type_3[cluster].from_2;
        int from3 = d_virtual_type_3[cluster].from_3;
        if (atom_local_label[vatom] == 1)
        {
            int local_vatom = -1, local_from1 = -1, local_from2 = -1,
                local_from3 = -1;
            bool valid = Virtual_Atom_Get_Owned_Local_Id(
                vatom, local_atom_numbers, atom_local_id, atom_local_label, 3,
                vatom, -1, &local_vatom, invalid_local_layout);
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from1, local_atom_numbers, atom_local_id, atom_local_label,
                    3, vatom, 0, &local_from1, invalid_local_layout) &&
                valid;
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from2, local_atom_numbers, atom_local_id, atom_local_label,
                    3, vatom, 1, &local_from2, invalid_local_layout) &&
                valid;
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from3, local_atom_numbers, atom_local_id, atom_local_label,
                    3, vatom, 2, &local_from3, invalid_local_layout) &&
                valid;
            if (valid)
            {
                l_virtual_type_3[local_numbers[0]] = d_virtual_type_3[cluster];
                l_virtual_type_3[local_numbers[0]].virtual_atom = local_vatom;
                l_virtual_type_3[local_numbers[0]].from_1 = local_from1;
                l_virtual_type_3[local_numbers[0]].from_2 = local_from2;
                l_virtual_type_3[local_numbers[0]].from_3 = local_from3;
                local_numbers[0] += 1;
            }
        }
    }
}

static __global__ void get_local_device_V5(
    int virtual_numbers, int* local_numbers,
    const VIRTUAL_TYPE_5* d_virtual_type_5, VIRTUAL_TYPE_5* l_virtual_type_5,
    const int* atom_local_id, const char* atom_local_label,
    int local_atom_numbers, int* invalid_local_layout)
{
    local_numbers[0] = 0;
    for (int cluster = 0; cluster < virtual_numbers; cluster++)
    {
        int vatom = d_virtual_type_5[cluster].virtual_atom;
        int from1 = d_virtual_type_5[cluster].from_1;
        int from2 = d_virtual_type_5[cluster].from_2;
        int from3 = d_virtual_type_5[cluster].from_3;
        if (atom_local_label[vatom] == 1)
        {
            int local_vatom = -1, local_from1 = -1, local_from2 = -1,
                local_from3 = -1;
            bool valid = Virtual_Atom_Get_Owned_Local_Id(
                vatom, local_atom_numbers, atom_local_id, atom_local_label, 5,
                vatom, -1, &local_vatom, invalid_local_layout);
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from1, local_atom_numbers, atom_local_id, atom_local_label,
                    5, vatom, 0, &local_from1, invalid_local_layout) &&
                valid;
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from2, local_atom_numbers, atom_local_id, atom_local_label,
                    5, vatom, 1, &local_from2, invalid_local_layout) &&
                valid;
            valid =
                Virtual_Atom_Get_Owned_Local_Id(
                    from3, local_atom_numbers, atom_local_id, atom_local_label,
                    5, vatom, 2, &local_from3, invalid_local_layout) &&
                valid;
            if (valid)
            {
                l_virtual_type_5[local_numbers[0]] = d_virtual_type_5[cluster];
                l_virtual_type_5[local_numbers[0]].virtual_atom = local_vatom;
                l_virtual_type_5[local_numbers[0]].from_1 = local_from1;
                l_virtual_type_5[local_numbers[0]].from_2 = local_from2;
                l_virtual_type_5[local_numbers[0]].from_3 = local_from3;
                local_numbers[0] += 1;
            }
        }
    }
}

// 预留get_local_device_V4接口

void VIRTUAL_INFORMATION::Get_Local(const int* atom_local_id,
                                    const char* atom_local_label,
                                    const int local_atom_numbers)
{
    if (!is_initialized) return;
    if (local_atom_numbers < 0 || local_atom_numbers > global_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "VIRTUAL_INFORMATION::Get_Local",
            "Reason:\n\towned local atom count %d is outside [0, %d]\n",
            local_atom_numbers, global_atom_numbers);
    }
    if (d_invalid_local_layout == NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "VIRTUAL_INFORMATION::Get_Local",
            "Reason:\n\tvirtual-atom local-layout validation storage is not "
            "initialized\n");
    }
    const int valid_local_layout[6] = {-1, -1, -1, -1, -1, -1};
    deviceMemcpy(d_invalid_local_layout, valid_local_layout,
                 sizeof(valid_local_layout), deviceMemcpyHostToDevice);
    auto validate_local_count =
        [&](int virtual_type, int local_count, int global_count)
    {
        if (local_count < 0 || local_count > global_count)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "VIRTUAL_INFORMATION::Get_Local",
                "Reason:\n\ttype-%d local virtual-atom count %d is outside "
                "[0, %d]\n",
                virtual_type, local_count, global_count);
        }
    };
    // 每层之间需要串行计算，层内并行计算
    for (int layer = 0; layer < max_level; layer++)
    {
        VIRTUAL_LAYER_INFORMATION* temp_vl = &virtual_layer_info[layer];

        if (temp_vl->v0_info.virtual_numbers > 0)
        {
            Launch_Device_Kernel(get_local_device_V0, 1, 1, 0, NULL,
                                 temp_vl->v0_info.virtual_numbers,
                                 temp_vl->v0_info.d_local_numbers,
                                 temp_vl->v0_info.d_virtual_type_0,
                                 temp_vl->v0_info.l_virtual_type_0,
                                 atom_local_id, atom_local_label,
                                 local_atom_numbers, d_invalid_local_layout);
            deviceMemcpy(&temp_vl->v0_info.local_numbers,
                         temp_vl->v0_info.d_local_numbers, sizeof(int),
                         deviceMemcpyDeviceToHost);
            validate_local_count(0, temp_vl->v0_info.local_numbers,
                                 temp_vl->v0_info.virtual_numbers);
        }

        if (temp_vl->v1_info.virtual_numbers > 0)
        {
            Launch_Device_Kernel(get_local_device_V1, 1, 1, 0, NULL,
                                 temp_vl->v1_info.virtual_numbers,
                                 temp_vl->v1_info.d_local_numbers,
                                 temp_vl->v1_info.d_virtual_type_1,
                                 temp_vl->v1_info.l_virtual_type_1,
                                 atom_local_id, atom_local_label,
                                 local_atom_numbers, d_invalid_local_layout);
            deviceMemcpy(&temp_vl->v1_info.local_numbers,
                         temp_vl->v1_info.d_local_numbers, sizeof(int),
                         deviceMemcpyDeviceToHost);
            validate_local_count(1, temp_vl->v1_info.local_numbers,
                                 temp_vl->v1_info.virtual_numbers);
        }

        if (temp_vl->v2_info.virtual_numbers > 0)
        {
            Launch_Device_Kernel(get_local_device_V2, 1, 1, 0, NULL,
                                 temp_vl->v2_info.virtual_numbers,
                                 temp_vl->v2_info.d_local_numbers,
                                 temp_vl->v2_info.d_virtual_type_2,
                                 temp_vl->v2_info.l_virtual_type_2,
                                 atom_local_id, atom_local_label,
                                 local_atom_numbers, d_invalid_local_layout);
            deviceMemcpy(&temp_vl->v2_info.local_numbers,
                         temp_vl->v2_info.d_local_numbers, sizeof(int),
                         deviceMemcpyDeviceToHost);
            validate_local_count(2, temp_vl->v2_info.local_numbers,
                                 temp_vl->v2_info.virtual_numbers);
        }

        if (temp_vl->v3_info.virtual_numbers > 0)
        {
            Launch_Device_Kernel(get_local_device_V3, 1, 1, 0, NULL,
                                 temp_vl->v3_info.virtual_numbers,
                                 temp_vl->v3_info.d_local_numbers,
                                 temp_vl->v3_info.d_virtual_type_3,
                                 temp_vl->v3_info.l_virtual_type_3,
                                 atom_local_id, atom_local_label,
                                 local_atom_numbers, d_invalid_local_layout);
            deviceMemcpy(&temp_vl->v3_info.local_numbers,
                         temp_vl->v3_info.d_local_numbers, sizeof(int),
                         deviceMemcpyDeviceToHost);
            validate_local_count(3, temp_vl->v3_info.local_numbers,
                                 temp_vl->v3_info.virtual_numbers);
        }

        if (temp_vl->v5_info.virtual_numbers > 0)
        {
            Launch_Device_Kernel(get_local_device_V5, 1, 1, 0, NULL,
                                 temp_vl->v5_info.virtual_numbers,
                                 temp_vl->v5_info.d_local_numbers,
                                 temp_vl->v5_info.d_virtual_type_5,
                                 temp_vl->v5_info.l_virtual_type_5,
                                 atom_local_id, atom_local_label,
                                 local_atom_numbers, d_invalid_local_layout);
            deviceMemcpy(&temp_vl->v5_info.local_numbers,
                         temp_vl->v5_info.d_local_numbers, sizeof(int),
                         deviceMemcpyDeviceToHost);
            validate_local_count(5, temp_vl->v5_info.local_numbers,
                                 temp_vl->v5_info.virtual_numbers);
        }

        // 预留v4质心接口
    }
    int invalid_layout[6] = {-1, -1, -1, -1, -1, -1};
    deviceMemcpy(invalid_layout, d_invalid_local_layout, sizeof(invalid_layout),
                 deviceMemcpyDeviceToHost);
    if (invalid_layout[0] >= 0)
    {
        const char* role = invalid_layout[5] < 0 ? "target" : "source";
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "VIRTUAL_INFORMATION::Get_Local",
            "Reason:\n\ttype-%d virtual atom %d has %s atom %d with "
            "ownership label %d and local id %d for %d owned atoms; every "
            "target and parent must be owned by the same domain because "
            "reverse virtual-site ghost-force communication is unavailable\n",
            invalid_layout[0], invalid_layout[1], role, invalid_layout[2],
            invalid_layout[4], invalid_layout[3], local_atom_numbers);
    }
    has_local_layout = true;
}

void VIRTUAL_INFORMATION::update_ug_connectivity(CONECT* connectivity)
{
    if (!is_initialized) return;
    for (int layer = 0; layer < max_level; layer++)
    {
        VIRTUAL_LAYER_INFORMATION& layer_info = virtual_layer_info[layer];
        for (int i = 0; i < layer_info.v0_info.virtual_numbers; i++)
        {
            int atomv = layer_info.v0_info.h_virtual_type_0[i].virtual_atom;
            int atom1 = layer_info.v0_info.h_virtual_type_0[i].from_1;
            (*connectivity)[atomv].insert(atom1);
            (*connectivity)[atom1].insert(atomv);
        }
        for (int i = 0; i < layer_info.v1_info.virtual_numbers; i++)
        {
            int atomv = layer_info.v1_info.h_virtual_type_1[i].virtual_atom;
            int atom1 = layer_info.v1_info.h_virtual_type_1[i].from_1;
            int atom2 = layer_info.v1_info.h_virtual_type_1[i].from_2;
            (*connectivity)[atomv].insert(atom1);
            (*connectivity)[atomv].insert(atom2);
            (*connectivity)[atom1].insert(atomv);
            (*connectivity)[atom1].insert(atom2);
            (*connectivity)[atom2].insert(atom1);
            (*connectivity)[atom2].insert(atomv);
        }
        for (int i = 0; i < layer_info.v2_info.virtual_numbers; i++)
        {
            int atomv = layer_info.v2_info.h_virtual_type_2[i].virtual_atom;
            int atom1 = layer_info.v2_info.h_virtual_type_2[i].from_1;
            int atom2 = layer_info.v2_info.h_virtual_type_2[i].from_2;
            int atom3 = layer_info.v2_info.h_virtual_type_2[i].from_3;
            (*connectivity)[atomv].insert(atom1);
            (*connectivity)[atomv].insert(atom2);
            (*connectivity)[atomv].insert(atom3);
            (*connectivity)[atom1].insert(atomv);
            (*connectivity)[atom1].insert(atom2);
            (*connectivity)[atom1].insert(atom3);
            (*connectivity)[atom2].insert(atom1);
            (*connectivity)[atom2].insert(atomv);
            (*connectivity)[atom2].insert(atom3);
            (*connectivity)[atom3].insert(atom1);
            (*connectivity)[atom3].insert(atom2);
            (*connectivity)[atom3].insert(atomv);
        }
        for (int i = 0; i < layer_info.v3_info.virtual_numbers; i++)
        {
            int atomv = layer_info.v3_info.h_virtual_type_3[i].virtual_atom;
            int atom1 = layer_info.v3_info.h_virtual_type_3[i].from_1;
            int atom2 = layer_info.v3_info.h_virtual_type_3[i].from_2;
            int atom3 = layer_info.v3_info.h_virtual_type_3[i].from_3;
            (*connectivity)[atomv].insert(atom1);
            (*connectivity)[atomv].insert(atom2);
            (*connectivity)[atomv].insert(atom3);
            (*connectivity)[atom1].insert(atomv);
            (*connectivity)[atom1].insert(atom2);
            (*connectivity)[atom1].insert(atom3);
            (*connectivity)[atom2].insert(atom1);
            (*connectivity)[atom2].insert(atomv);
            (*connectivity)[atom2].insert(atom3);
            (*connectivity)[atom3].insert(atom1);
            (*connectivity)[atom3].insert(atom2);
            (*connectivity)[atom3].insert(atomv);
        }
        for (int i = 0; i < layer_info.v5_info.virtual_numbers; i++)
        {
            int atomv = layer_info.v5_info.h_virtual_type_5[i].virtual_atom;
            int atom1 = layer_info.v5_info.h_virtual_type_5[i].from_1;
            int atom2 = layer_info.v5_info.h_virtual_type_5[i].from_2;
            int atom3 = layer_info.v5_info.h_virtual_type_5[i].from_3;
            (*connectivity)[atomv].insert(atom1);
            (*connectivity)[atomv].insert(atom2);
            (*connectivity)[atomv].insert(atom3);
            (*connectivity)[atom1].insert(atomv);
            (*connectivity)[atom1].insert(atom2);
            (*connectivity)[atom1].insert(atom3);
            (*connectivity)[atom2].insert(atom1);
            (*connectivity)[atom2].insert(atomv);
            (*connectivity)[atom2].insert(atom3);
            (*connectivity)[atom3].insert(atom1);
            (*connectivity)[atom3].insert(atom2);
            (*connectivity)[atom3].insert(atomv);
        }
    }
}
