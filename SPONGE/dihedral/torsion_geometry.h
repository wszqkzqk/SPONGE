#pragma once

#include "../common.h"

struct TORSION_DOUBLE_VECTOR
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct TORSION_GEOMETRY
{
    double phi = 0.0;
    TORSION_DOUBLE_VECTOR dphi_dri;
    TORSION_DOUBLE_VECTOR dphi_drj;
    TORSION_DOUBLE_VECTOR dphi_drk;
    TORSION_DOUBLE_VECTOR dphi_drl;
    VECTOR drij = {0.0f, 0.0f, 0.0f};
    VECTOR drkj = {0.0f, 0.0f, 0.0f};
    VECTOR drkl = {0.0f, 0.0f, 0.0f};
};

static __device__ __forceinline__ bool Torsion_Float_Is_Finite(float value)
{
#ifdef GPU_ARCH_NAME
    // Use the representation directly: CUDA/HIP fast-math must not be able to
    // infer that an arithmetic value is finite and fold the validation away.
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
#elif defined(__GNUC__) || defined(__clang__)
    unsigned int bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires 32-bit IEEE-754 floats");
    memcpy(&bits, &value, sizeof(value));
    // -ffast-math otherwise lets the compiler assume that a float argument is
    // finite and fold this test away.  The empty tied asm is a zero-instruction
    // optimizer barrier; the actual test remains an inlined integer mask.
    __asm__ __volatile__("" : "+r"(bits));
    return (bits & 0x7f800000U) != 0x7f800000U;
#else
    return Float_Memory_Is_Finite(&value);
#endif
}

static __device__ __forceinline__ bool Torsion_Vector_Is_Finite(
    const VECTOR& value)
{
    return Torsion_Float_Is_Finite(value.x) &&
           Torsion_Float_Is_Finite(value.y) && Torsion_Float_Is_Finite(value.z);
}

static __device__ __forceinline__ bool Torsion_Double_Is_Finite(double value)
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

static __device__ __forceinline__ bool Torsion_Double_Vector_Is_Finite(
    const TORSION_DOUBLE_VECTOR& value)
{
    return Torsion_Double_Is_Finite(value.x) &&
           Torsion_Double_Is_Finite(value.y) &&
           Torsion_Double_Is_Finite(value.z);
}

static __device__ __forceinline__ TORSION_DOUBLE_VECTOR
Torsion_To_Double_Vector(const VECTOR& value)
{
    return {static_cast<double>(value.x), static_cast<double>(value.y),
            static_cast<double>(value.z)};
}

static __device__ __forceinline__ double Torsion_Double_Dot(
    const TORSION_DOUBLE_VECTOR& first, const TORSION_DOUBLE_VECTOR& second)
{
    return first.x * second.x + first.y * second.y + first.z * second.z;
}

static __device__ __forceinline__ TORSION_DOUBLE_VECTOR Torsion_Double_Cross(
    const TORSION_DOUBLE_VECTOR& first, const TORSION_DOUBLE_VECTOR& second)
{
    return {first.y * second.z - first.z * second.y,
            first.z * second.x - first.x * second.z,
            first.x * second.y - first.y * second.x};
}

// Narrow only when the rounded float is finite and does not erase a nonzero
// double result.
static __device__ __forceinline__ bool Torsion_Checked_Narrow(
    double value, float* result)
{
    // FLT_MAX: std::numeric_limits<float>::max() is a host-only constexpr
    // under nvcc and cannot be called from this device function.
    const double maximum_float = static_cast<double>(FLT_MAX);
    if (!Torsion_Double_Is_Finite(value) || value > maximum_float ||
        value < -maximum_float)
    {
        return false;
    }
    const float narrowed = static_cast<float>(value);
    if (!Torsion_Float_Is_Finite(narrowed) ||
        (value != 0.0 && narrowed == 0.0f))
    {
        return false;
    }
    *result = narrowed;
    return true;
}

static __device__ __forceinline__ bool Torsion_Checked_Scale_And_Narrow(
    double scale, const TORSION_DOUBLE_VECTOR& value, VECTOR* result)
{
    VECTOR narrowed;
    if (!Torsion_Checked_Narrow(scale * value.x, &narrowed.x) ||
        !Torsion_Checked_Narrow(scale * value.y, &narrowed.y) ||
        !Torsion_Checked_Narrow(scale * value.z, &narrowed.z))
    {
        return false;
    }
    *result = narrowed;
    return true;
}

// Form the complete term virial in double, then narrow each final component.
// This avoids rejecting a representable sum merely because one float product
// overflowed before cancellation.
static __device__ __forceinline__ bool Torsion_Checked_Term_Virial(
    const VECTOR& force_i, const TORSION_DOUBLE_VECTOR& displacement_i,
    const VECTOR& force_k, const TORSION_DOUBLE_VECTOR& displacement_k,
    const VECTOR& force_l, const TORSION_DOUBLE_VECTOR& displacement_l,
    LTMatrix3* result)
{
    const TORSION_DOUBLE_VECTOR forces[3] = {
        Torsion_To_Double_Vector(force_i), Torsion_To_Double_Vector(force_k),
        Torsion_To_Double_Vector(force_l)};
    const TORSION_DOUBLE_VECTOR displacements[3] = {
        displacement_i, displacement_k, displacement_l};
    double values[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int pair = 0; pair < 3; pair++)
    {
        values[0] += forces[pair].x * displacements[pair].x;
        values[1] += forces[pair].x * displacements[pair].y +
                     forces[pair].y * displacements[pair].x;
        values[2] += forces[pair].y * displacements[pair].y;
        values[3] += forces[pair].x * displacements[pair].z +
                     forces[pair].z * displacements[pair].x;
        values[4] += forces[pair].y * displacements[pair].z +
                     forces[pair].z * displacements[pair].y;
        values[5] += forces[pair].z * displacements[pair].z;
    }

    LTMatrix3 narrowed;
    if (!Torsion_Checked_Narrow(values[0], &narrowed.a11) ||
        !Torsion_Checked_Narrow(values[1], &narrowed.a21) ||
        !Torsion_Checked_Narrow(values[2], &narrowed.a22) ||
        !Torsion_Checked_Narrow(values[3], &narrowed.a31) ||
        !Torsion_Checked_Narrow(values[4], &narrowed.a32) ||
        !Torsion_Checked_Narrow(values[5], &narrowed.a33))
    {
        return false;
    }
    *result = narrowed;
    return true;
}

static __device__ __forceinline__ bool Torsion_Finite_Atomic_Add(
    VECTOR* accumulator, const VECTOR& value)
{
    return Finite_Atomic_Add(&accumulator->x, value.x) &&
           Finite_Atomic_Add(&accumulator->y, value.y) &&
           Finite_Atomic_Add(&accumulator->z, value.z);
}

static __device__ __forceinline__ bool Torsion_Finite_Atomic_Add(
    LTMatrix3* accumulator, const LTMatrix3& value)
{
    return Finite_Atomic_Add(&accumulator->a11, value.a11) &&
           Finite_Atomic_Add(&accumulator->a21, value.a21) &&
           Finite_Atomic_Add(&accumulator->a22, value.a22) &&
           Finite_Atomic_Add(&accumulator->a31, value.a31) &&
           Finite_Atomic_Add(&accumulator->a32, value.a32) &&
           Finite_Atomic_Add(&accumulator->a33, value.a33);
}

// Compute the signed SPONGE/GROMACS torsion and its Cartesian Jacobian.
// atan2 keeps the angle and force well-conditioned near 0 and pi.  All
// geometry after the float minimum-image displacement is evaluated in double:
// a nonzero normal can have a square below the float range even though the
// final Cartesian force is readily representable.
static __device__ __forceinline__ bool Compute_Torsion_Geometry(
    const VECTOR& ri, const VECTOR& rj, const VECTOR& rk, const VECTOR& rl,
    const LTMatrix3& cell, const LTMatrix3& rcell, TORSION_GEOMETRY* geometry)
{
    geometry->drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
    geometry->drkj = Get_Periodic_Displacement(rk, rj, cell, rcell);
    geometry->drkl = Get_Periodic_Displacement(rk, rl, cell, rcell);
    if (!Torsion_Vector_Is_Finite(geometry->drij) ||
        !Torsion_Vector_Is_Finite(geometry->drkj) ||
        !Torsion_Vector_Is_Finite(geometry->drkl))
    {
        return false;
    }

    const TORSION_DOUBLE_VECTOR arm_i =
        Torsion_To_Double_Vector(geometry->drij);
    const TORSION_DOUBLE_VECTOR central =
        Torsion_To_Double_Vector(geometry->drkj);
    const TORSION_DOUBLE_VECTOR arm_l =
        Torsion_To_Double_Vector(geometry->drkl);
    const TORSION_DOUBLE_VECTOR normal_1 =
        Torsion_Double_Cross(arm_i, central);
    const TORSION_DOUBLE_VECTOR normal_2 =
        Torsion_Double_Cross(arm_l, central);
    const double normal_1_squared = Torsion_Double_Dot(normal_1, normal_1);
    const double normal_2_squared = Torsion_Double_Dot(normal_2, normal_2);
    const double central_bond_squared = Torsion_Double_Dot(central, central);
    if (!(normal_1_squared > 0.0) || !(normal_2_squared > 0.0) ||
        !(central_bond_squared > 0.0) ||
        !Torsion_Double_Is_Finite(normal_1_squared) ||
        !Torsion_Double_Is_Finite(normal_2_squared) ||
        !Torsion_Double_Is_Finite(central_bond_squared))
    {
        return false;
    }

    const double normal_1_length = sqrt(normal_1_squared);
    const double normal_2_length = sqrt(normal_2_squared);
    const double central_bond_length = sqrt(central_bond_squared);
    const double inverse_normal_1 = 1.0 / normal_1_length;
    const double inverse_normal_2 = 1.0 / normal_2_length;
    const double inverse_central_bond = 1.0 / central_bond_length;
    if (!Torsion_Double_Is_Finite(normal_1_length) ||
        !Torsion_Double_Is_Finite(normal_2_length) ||
        !Torsion_Double_Is_Finite(central_bond_length) ||
        !Torsion_Double_Is_Finite(inverse_normal_1) ||
        !Torsion_Double_Is_Finite(inverse_normal_2) ||
        !Torsion_Double_Is_Finite(inverse_central_bond))
    {
        return false;
    }

    const TORSION_DOUBLE_VECTOR unit_normal_1 = {
        normal_1.x * inverse_normal_1, normal_1.y * inverse_normal_1,
        normal_1.z * inverse_normal_1};
    const TORSION_DOUBLE_VECTOR unit_normal_2 = {
        normal_2.x * inverse_normal_2, normal_2.y * inverse_normal_2,
        normal_2.z * inverse_normal_2};
    const double cosine = -Torsion_Double_Dot(unit_normal_1, unit_normal_2);
    const double sine =
        Torsion_Double_Dot(Torsion_Double_Cross(unit_normal_2, unit_normal_1),
                           central) *
        inverse_central_bond;
    if (!Torsion_Double_Is_Finite(cosine) ||
        !Torsion_Double_Is_Finite(sine) || (cosine == 0.0 && sine == 0.0))
    {
        return false;
    }
    geometry->phi = atan2(sine, cosine);

    const double scale_i = central_bond_length * inverse_normal_1;
    const double scale_l = central_bond_length * inverse_normal_2;
    geometry->dphi_dri = {scale_i * unit_normal_1.x,
                          scale_i * unit_normal_1.y,
                          scale_i * unit_normal_1.z};
    geometry->dphi_drl = {scale_l * unit_normal_2.x,
                          scale_l * unit_normal_2.y,
                          scale_l * unit_normal_2.z};

    const double projection_i =
        Torsion_Double_Dot(arm_i, central) / central_bond_squared;
    const double projection_l =
        Torsion_Double_Dot(arm_l, central) / central_bond_squared;
    geometry->dphi_drj = {
        (projection_i - 1.0) * geometry->dphi_dri.x -
            projection_l * geometry->dphi_drl.x,
        (projection_i - 1.0) * geometry->dphi_dri.y -
            projection_l * geometry->dphi_drl.y,
        (projection_i - 1.0) * geometry->dphi_dri.z -
            projection_l * geometry->dphi_drl.z};
    geometry->dphi_drk = {
        -projection_i * geometry->dphi_dri.x +
            (projection_l - 1.0) * geometry->dphi_drl.x,
        -projection_i * geometry->dphi_dri.y +
            (projection_l - 1.0) * geometry->dphi_drl.y,
        -projection_i * geometry->dphi_dri.z +
            (projection_l - 1.0) * geometry->dphi_drl.z};

    return Torsion_Double_Is_Finite(geometry->phi) &&
           Torsion_Double_Vector_Is_Finite(geometry->dphi_dri) &&
           Torsion_Double_Vector_Is_Finite(geometry->dphi_drj) &&
           Torsion_Double_Vector_Is_Finite(geometry->dphi_drk) &&
           Torsion_Double_Vector_Is_Finite(geometry->dphi_drl);
}
