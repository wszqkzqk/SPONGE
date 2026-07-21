#pragma once

#include "../../common.h"

enum REAXFF_GEOMETRY_ERROR_CODE
{
    REAXFF_GEOMETRY_OK = 0,
    REAXFF_BOND_OVERLAP = 1,
    REAXFF_BOND_NONFINITE = 2,
    REAXFF_ANGLE_UNDEFINED = 3,
    REAXFF_ANGLE_NONFINITE = 4,
    REAXFF_TORSION_UNDEFINED = 5,
    REAXFF_TORSION_NONFINITE = 6,
    REAXFF_HB_UNDEFINED = 7,
    REAXFF_HB_NONFINITE = 8,
    REAXFF_VDW_OVERLAP = 9,
    REAXFF_VDW_NONFINITE = 10,
    REAXFF_INVALID_ATOM_TYPE = 11
};

constexpr int REAXFF_GEOMETRY_ERROR_SIZE = 5;

static __device__ __forceinline__ bool ReaxFF_Float_Is_Finite(float value)
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

static __device__ __forceinline__ bool ReaxFF_Double_Is_Finite(double value)
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

static __device__ __forceinline__ bool ReaxFF_Vector_Is_Finite(
    const VECTOR& value)
{
    return ReaxFF_Float_Is_Finite(value.x) &&
           ReaxFF_Float_Is_Finite(value.y) &&
           ReaxFF_Float_Is_Finite(value.z);
}

static __device__ __forceinline__ bool ReaxFF_Matrix_Is_Finite(
    const LTMatrix3& value)
{
    return ReaxFF_Float_Is_Finite(value.a11) &&
           ReaxFF_Float_Is_Finite(value.a21) &&
           ReaxFF_Float_Is_Finite(value.a22) &&
           ReaxFF_Float_Is_Finite(value.a31) &&
           ReaxFF_Float_Is_Finite(value.a32) &&
           ReaxFF_Float_Is_Finite(value.a33);
}

// Store the first error.  The record is [code, i, j, k, l].  Device kernels
// return from the failing interaction and the host checks immediately after
// launch, before a partially evaluated force field can be consumed.
static __device__ __forceinline__ void Record_ReaxFF_Geometry_Error(
    int* error, int code, int i, int j = -1, int k = -1, int l = -1)
{
    if (error == NULL) return;
#ifdef GPU_ARCH_NAME
    if (atomicCAS(error, REAXFF_GEOMETRY_OK, -1) == REAXFF_GEOMETRY_OK)
    {
        error[1] = i;
        error[2] = j;
        error[3] = k;
        error[4] = l;
        __threadfence();
        atomicExch(error, code);
    }
#else
#pragma omp critical(sponge_reaxff_geometry_error)
    {
        if (error[0] == REAXFF_GEOMETRY_OK)
        {
            error[1] = i;
            error[2] = j;
            error[3] = k;
            error[4] = l;
            error[0] = code;
        }
    }
#endif
}

struct REAXFF_ANGLE_GEOMETRY
{
    float theta = 0.0f;
    float sine = 0.0f;
    float cosine = 1.0f;
    // sin(theta / 2)^2, evaluated without subtracting a cosine close to one.
    // Hydrogen-bond terms use this directly so small but representable angles
    // do not collapse to an exactly zero interaction after float rounding.
    float half_sine_squared = 0.0f;
    bool is_collinear = false;
    double dtheta_du[3] = {0.0, 0.0, 0.0};
    double dtheta_dv[3] = {0.0, 0.0, 0.0};
};

// u and v point from the central atom to the two outer atoms.  atan2 and the
// cross-product Jacobian preserve every representable non-collinear angle;
// exact endpoints are reported separately because their Cartesian derivative
// is not generally unique.
static __device__ __forceinline__ bool Compute_ReaxFF_Angle_Geometry(
    const VECTOR& u, const VECTOR& v, REAXFF_ANGLE_GEOMETRY* geometry)
{
    if (geometry == NULL || !ReaxFF_Vector_Is_Finite(u) ||
        !ReaxFF_Vector_Is_Finite(v))
    {
        return false;
    }
    const double ux = static_cast<double>(u.x);
    const double uy = static_cast<double>(u.y);
    const double uz = static_cast<double>(u.z);
    const double vx = static_cast<double>(v.x);
    const double vy = static_cast<double>(v.y);
    const double vz = static_cast<double>(v.z);
    const double u2 = ux * ux + uy * uy + uz * uz;
    const double v2 = vx * vx + vy * vy + vz * vz;
    if (!(u2 > 0.0) || !(v2 > 0.0) || !ReaxFF_Double_Is_Finite(u2) ||
        !ReaxFF_Double_Is_Finite(v2))
    {
        return false;
    }

    const double nx = uy * vz - uz * vy;
    const double ny = uz * vx - ux * vz;
    const double nz = ux * vy - uy * vx;
    const double n2 = nx * nx + ny * ny + nz * nz;
    const double dot = ux * vx + uy * vy + uz * vz;
    if (!(n2 >= 0.0) || !ReaxFF_Double_Is_Finite(n2) ||
        !ReaxFF_Double_Is_Finite(dot))
    {
        return false;
    }

    const double n = sqrt(n2);
    const double norm_product = sqrt(u2) * sqrt(v2);
    const double theta = atan2(n, dot);
    const double sine = n / norm_product;
    const double cosine = dot / norm_product;
    double half_sine_squared = 0.0;
    if (dot >= 0.0)
    {
        // (1-cos(theta))/2, rewritten with the cross product to avoid
        // catastrophic cancellation for theta close to zero.
        const double denominator =
            2.0 * norm_product * (norm_product + dot);
        if (!(denominator > 0.0) ||
            !ReaxFF_Double_Is_Finite(denominator))
        {
            return false;
        }
        half_sine_squared = n2 / denominator;
    }
    else
    {
        // Near pi the direct expression is well conditioned, including the
        // exact antiparallel endpoint.
        half_sine_squared = 0.5 * (1.0 - cosine);
    }
    half_sine_squared = half_sine_squared > 1.0
                            ? 1.0
                            : (half_sine_squared < 0.0 ? 0.0
                                                       : half_sine_squared);
    geometry->theta = static_cast<float>(theta);
    geometry->sine = static_cast<float>(sine);
    geometry->cosine = static_cast<float>(
        cosine > 1.0 ? 1.0 : (cosine < -1.0 ? -1.0 : cosine));
    geometry->half_sine_squared =
        static_cast<float>(half_sine_squared);
    geometry->is_collinear = n2 == 0.0;
    if (!ReaxFF_Double_Is_Finite(n) || !ReaxFF_Double_Is_Finite(theta) ||
        !ReaxFF_Double_Is_Finite(sine) ||
        !ReaxFF_Double_Is_Finite(cosine) ||
        !ReaxFF_Double_Is_Finite(half_sine_squared) ||
        !ReaxFF_Float_Is_Finite(geometry->theta) ||
        !ReaxFF_Float_Is_Finite(geometry->sine) ||
        !ReaxFF_Float_Is_Finite(geometry->cosine) ||
        !ReaxFF_Float_Is_Finite(geometry->half_sine_squared) ||
        (theta != 0.0 && geometry->theta == 0.0f))
    {
        return false;
    }
    if (geometry->is_collinear) return dot != 0.0;
    if (!(geometry->sine > 0.0f)) return false;

    const double inverse_n = 1.0 / n;
    const double inverse_u2 = 1.0 / u2;
    const double inverse_v2 = 1.0 / v2;
    const double nhat_x = nx * inverse_n;
    const double nhat_y = ny * inverse_n;
    const double nhat_z = nz * inverse_n;
    geometry->dtheta_du[0] = (uy * nhat_z - uz * nhat_y) * inverse_u2;
    geometry->dtheta_du[1] = (uz * nhat_x - ux * nhat_z) * inverse_u2;
    geometry->dtheta_du[2] = (ux * nhat_y - uy * nhat_x) * inverse_u2;
    geometry->dtheta_dv[0] = (nhat_y * vz - nhat_z * vy) * inverse_v2;
    geometry->dtheta_dv[1] = (nhat_z * vx - nhat_x * vz) * inverse_v2;
    geometry->dtheta_dv[2] = (nhat_x * vy - nhat_y * vx) * inverse_v2;
    for (int axis = 0; axis < 3; axis++)
    {
        if (!ReaxFF_Double_Is_Finite(geometry->dtheta_du[axis]) ||
            !ReaxFF_Double_Is_Finite(geometry->dtheta_dv[axis]))
        {
            return false;
        }
    }
    return true;
}
