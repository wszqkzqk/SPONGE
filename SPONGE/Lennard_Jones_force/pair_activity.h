#pragma once

namespace PairwiseInteraction
{

enum Pair_Component
{
    PAIR_COMPONENT_NONE = 0,
    PAIR_COMPONENT_LENNARD_JONES = 1,
    PAIR_COMPONENT_COULOMB = 2,
    PAIR_COMPONENT_GENERALIZED_BORN = 4
};

// These tests intentionally use exact floating-point zero. Pair parameters
// are coefficients, not geometric tolerances: a nonzero coefficient must
// continue through the normal interaction path regardless of its magnitude.
__host__ __device__ __forceinline__ bool Lennard_Jones_Is_Active(
    const float coefficient_a, const float coefficient_b)
{
    return coefficient_a != 0.0f || coefficient_b != 0.0f;
}

__host__ __device__ __forceinline__ bool Coulomb_Is_Active(const float charge_i,
                                                           const float charge_j)
{
    // Activity is a property of the two input coefficients, not of a rounded
    // float product.  In particular, two representable nonzero charges remain
    // active even when their product underflows in single precision.
    return charge_i != 0.0f && charge_j != 0.0f;
}

__host__ __device__ __forceinline__ bool Coulomb_Derivative_Is_Active(
    const float charge_i, const float charge_j, const float charge_derivative_i,
    const float charge_derivative_j)
{
    const double derivative =
        static_cast<double>(charge_derivative_j) * charge_i +
        static_cast<double>(charge_j) * charge_derivative_i;
    return derivative != 0.0;
}

struct Coulomb_Endpoint_Activity
{
    bool state_a = false;
    bool state_b = false;
    bool ever_active = false;

    __host__ __device__ __forceinline__ bool Changes() const
    {
        // A pair can be inactive at both endpoints but active in between, for
        // example (q_i: 1 -> 0, q_j: 0 -> 1).  Such a pair still needs the
        // endpoint soft-core distances.
        return ever_active && !(state_a && state_b);
    }
};

__host__ __device__ __forceinline__ float Interpolate_Charge(
    const float charge_a, const float charge_b, const float lambda)
{
    // This is the canonical interpolation used by both the production
    // Hamiltonian and its lambda derivative.
    return fmaf(lambda, charge_b - charge_a, charge_a);
}

enum Charge_Endpoint_Error
{
    CHARGE_ENDPOINT_ERROR_NONE = 0,
    CHARGE_ENDPOINT_ERROR_NONFINITE = 1,
    CHARGE_ENDPOINT_ERROR_CURRENT_MISMATCH = 2,
    CHARGE_ENDPOINT_ERROR_DERIVATIVE_MISMATCH = 3
};

struct Charge_Endpoint_Validation
{
    float current = 0.0f;
    float derivative = 0.0f;
    int error = CHARGE_ENDPOINT_ERROR_NONE;
};

// This check must remain valid in Release builds compiled with -ffast-math.
// Inspecting IEEE-754 exponent bits avoids an isfinite() expression that the
// compiler is otherwise allowed to fold to true under finite-math-only.
#ifdef GPU_ARCH_NAME
static __device__ __forceinline__ bool Charge_Float_Is_Finite(
    const float value)
{
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
}
#else
static __host__ __device__ __forceinline__ bool Charge_Float_Is_Finite(
    const float value)
{
    unsigned int bits = 0;
    static_assert(sizeof(bits) == sizeof(value) &&
                      std::numeric_limits<float>::is_iec559,
                  "SPONGE requires 32-bit IEEE-754 floats");
    memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return (bits & 0x7f800000U) != 0x7f800000U;
}
#endif

#ifdef GPU_ARCH_NAME
static __device__ __forceinline__ Charge_Endpoint_Validation
#else
static __host__ __device__ __forceinline__ Charge_Endpoint_Validation
#endif
Validate_Charge_Endpoints(const float supplied_current, const float charge_a,
                          const float charge_b,
                          const float* supplied_derivative,
                          const float lambda)
{
    Charge_Endpoint_Validation result;
    result.current = Interpolate_Charge(charge_a, charge_b, lambda);
    result.derivative = charge_b - charge_a;
    if (!Charge_Float_Is_Finite(supplied_current) ||
        !Charge_Float_Is_Finite(charge_a) ||
        !Charge_Float_Is_Finite(charge_b) ||
        !Charge_Float_Is_Finite(result.current) ||
        !Charge_Float_Is_Finite(result.derivative) ||
        (supplied_derivative != NULL &&
         !Charge_Float_Is_Finite(*supplied_derivative)))
    {
        result.error = CHARGE_ENDPOINT_ERROR_NONFINITE;
    }
    else if (supplied_current != result.current)
    {
        result.error = CHARGE_ENDPOINT_ERROR_CURRENT_MISMATCH;
    }
    else if (supplied_derivative != NULL &&
             *supplied_derivative != result.derivative)
    {
        result.error = CHARGE_ENDPOINT_ERROR_DERIVATIVE_MISMATCH;
    }
    return result;
}

__host__ __device__ __forceinline__ Coulomb_Endpoint_Activity
Classify_Coulomb_Endpoints(const float charge_a_i, const float charge_a_j,
                           const float charge_b_i, const float charge_b_j)
{
    const bool atom_i_ever_active =
        charge_a_i != 0.0f || charge_b_i != 0.0f;
    const bool atom_j_ever_active =
        charge_a_j != 0.0f || charge_b_j != 0.0f;
    return {Coulomb_Is_Active(charge_a_i, charge_a_j),
            Coulomb_Is_Active(charge_b_i, charge_b_j),
            atom_i_ever_active && atom_j_ever_active};
}

struct Pair_Activity
{
    bool lennard_jones;
    bool coulomb;

    __host__ __device__ __forceinline__ bool Any() const
    {
        return lennard_jones || coulomb;
    }
};

__host__ __device__ __forceinline__ Pair_Activity
Classify(const float coefficient_a, const float coefficient_b,
         const bool coulomb_is_active)
{
    return {Lennard_Jones_Is_Active(coefficient_a, coefficient_b),
            coulomb_is_active};
}

__host__ __device__ __forceinline__ int Components(const bool lennard_jones,
                                                   const bool coulomb)
{
    return (lennard_jones ? PAIR_COMPONENT_LENNARD_JONES
                          : PAIR_COMPONENT_NONE) |
           (coulomb ? PAIR_COMPONENT_COULOMB : PAIR_COMPONENT_NONE);
}

// Active hard interactions are undefined at exact overlap.  GPU kernels print
// the global pair and stop at the point of failure.  CPU kernels record one
// complete pair in a three-int buffer (component mask, atom i, atom j); the
// owning host wrapper turns it into a normal SPONGE fatal error after the
// OpenMP loop has joined.
static __device__ __forceinline__ void Fail_Exact_Overlap(
    const int global_atom_i, const int global_atom_j, const int components,
    int* overlap_error)
{
#ifdef GPU_ARCH_NAME
    if (components & PAIR_COMPONENT_GENERALIZED_BORN)
    {
        printf(
            "Fatal SPONGE generalized-Born error: global atoms %d %d "
            "overlap exactly; pairwise descreening is undefined.\n",
            global_atom_i, global_atom_j);
    }
    else if (components ==
             (PAIR_COMPONENT_LENNARD_JONES | PAIR_COMPONENT_COULOMB))
    {
        printf(
            "Fatal SPONGE nonbond error: global atoms %d %d overlap exactly "
            "with active LJ and Coulomb components.\n",
            global_atom_i, global_atom_j);
    }
    else if (components & PAIR_COMPONENT_LENNARD_JONES)
    {
        printf(
            "Fatal SPONGE nonbond error: global atoms %d %d overlap exactly "
            "with an active LJ component.\n",
            global_atom_i, global_atom_j);
    }
    else
    {
        printf(
            "Fatal SPONGE nonbond error: global atoms %d %d overlap exactly "
            "with an active Coulomb component.\n",
            global_atom_i, global_atom_j);
    }
#if defined(USE_CUDA)
    asm volatile("trap;");
#elif defined(USE_HIP)
    __builtin_trap();
#endif
#else
    if (overlap_error == NULL || components == PAIR_COMPONENT_NONE)
    {
        return;
    }
#if defined(__GNUC__) || defined(__clang__)
    int expected = PAIR_COMPONENT_NONE;
    if (__atomic_compare_exchange_n(overlap_error, &expected, components, false,
                                    __ATOMIC_RELAXED, __ATOMIC_RELAXED))
    {
        overlap_error[1] = global_atom_i;
        overlap_error[2] = global_atom_j;
    }
#else
#pragma omp critical(sponge_nonbond_overlap_error)
    {
        if (overlap_error[0] == PAIR_COMPONENT_NONE)
        {
            overlap_error[1] = global_atom_i;
            overlap_error[2] = global_atom_j;
            overlap_error[0] = components;
        }
    }
#endif
#endif
}

}  // namespace PairwiseInteraction
