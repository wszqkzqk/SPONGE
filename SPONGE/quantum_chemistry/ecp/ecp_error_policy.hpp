#pragma once

#include <cfloat>
#include <cstdint>
#include <cstring>
#include <limits>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define QC_ECP_ERROR_HOST_DEVICE __host__ __device__
#else
#define QC_ECP_ERROR_HOST_DEVICE
#endif

namespace qc_ecp_error_policy
{
inline constexpr double CONTRACTED_ERROR_TOLERANCE = 2.0e-8;

QC_ECP_ERROR_HOST_DEVICE inline double Abs(double value)
{
    return value < 0.0 ? -value : value;
}

// SPONGE uses finite-math optimizations. Inspect IEEE-754 bits so fault checks
// remain effective for NaN/Inf in both optimized host code and device kernels.
QC_ECP_ERROR_HOST_DEVICE inline bool Double_Is_Finite(double value)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    const unsigned long long bits =
        static_cast<unsigned long long>(__double_as_longlong(value));
    return (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
#else
    std::uint64_t bits = 0;
    static_assert(
        sizeof(bits) == sizeof(value) && std::numeric_limits<double>::is_iec559,
        "SPONGE requires 64-bit IEEE-754 doubles");
    std::memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return (bits & UINT64_C(0x7ff0000000000000)) !=
           UINT64_C(0x7ff0000000000000);
#endif
}

QC_ECP_ERROR_HOST_DEVICE inline bool Contracted_Error_Is_Acceptable(
    double magnitude, double absolute_error_estimate)
{
    return Double_Is_Finite(magnitude) &&
           Double_Is_Finite(absolute_error_estimate) && magnitude >= 0.0 &&
           absolute_error_estimate >= 0.0 &&
           absolute_error_estimate <=
               CONTRACTED_ERROR_TOLERANCE * (1.0 + magnitude);
}

QC_ECP_ERROR_HOST_DEVICE inline bool Series_Estimate_Is_Small_Enough(
    bool primitive_is_usable, double contraction_scale,
    double primitive_error_estimate)
{
    return primitive_is_usable && Double_Is_Finite(contraction_scale) &&
           Double_Is_Finite(primitive_error_estimate) &&
           primitive_error_estimate >= 0.0 &&
           Abs(contraction_scale) * primitive_error_estimate <=
               0.25 * CONTRACTED_ERROR_TOLERANCE;
}

struct Matrix_Storage_Assessment
{
    bool accepted = false;
    double stored_value = 0.0;
    double storage_error = DBL_MAX;
    double storage_half_ulp = 0.0;
    double total_error_estimate = DBL_MAX;
    double allowed_error = 0.0;
};

QC_ECP_ERROR_HOST_DEVICE inline std::uint32_t Float_Bits(float value)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return static_cast<std::uint32_t>(__float_as_uint(value));
#else
    std::uint32_t bits = 0;
    static_assert(
        sizeof(bits) == sizeof(value) && std::numeric_limits<float>::is_iec559,
        "SPONGE requires 32-bit IEEE-754 floats");
    std::memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return bits;
#endif
}

QC_ECP_ERROR_HOST_DEVICE inline double Double_Power_Of_Two_From_Exponent_Field(
    unsigned int exponent_field)
{
    const std::uint64_t bits =
        static_cast<std::uint64_t>(exponent_field) << 52U;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return __longlong_as_double(static_cast<long long>(bits));
#else
    std::uint64_t protected_bits = bits;
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(protected_bits));
#endif
    double value = 0.0;
    std::memcpy(&value, &protected_bits, sizeof(value));
    return value;
#endif
}

// Return half the spacing of the IEEE-754 binary32 binade containing stored.
// Zero and every subnormal share the minimum spacing 2^-149, hence half-ulp is
// 2^-150.  For a normal float with biased exponent e, half-ulp is 2^(e-151).
// Constructing the result through the binary64 exponent keeps this exact under
// finite-math optimization and also represents the subnormal-float allowance
// as a normal double.
QC_ECP_ERROR_HOST_DEVICE inline double Float_Storage_Half_ULP(float stored)
{
    const std::uint32_t exponent =
        (Float_Bits(stored) >> 23U) & UINT32_C(0xff);
    if (exponent == UINT32_C(0xff)) return DBL_MAX;
    const std::uint32_t effective_exponent = exponent == 0 ? 1 : exponent;
    return Double_Power_Of_Two_From_Exponent_Field(
        static_cast<unsigned int>(effective_exponent + UINT32_C(872)));
}

// The QC matrix workspace stores float values. The integration estimate must
// independently meet its original budget; the final stored observable then
// receives only the unavoidable half-ulp float representation allowance.
QC_ECP_ERROR_HOST_DEVICE inline Matrix_Storage_Assessment
Assess_Matrix_Storage(double pre_storage_value,
                      double integration_error_estimate)
{
    Matrix_Storage_Assessment result;
    result.stored_value = pre_storage_value;
    if (!Double_Is_Finite(pre_storage_value) ||
        !Double_Is_Finite(integration_error_estimate) ||
        integration_error_estimate < 0.0)
        return result;

    const float stored = static_cast<float>(pre_storage_value);
    result.stored_value = static_cast<double>(stored);
    if (!Double_Is_Finite(result.stored_value)) return result;

    result.storage_error = Abs(result.stored_value - pre_storage_value);
    result.storage_half_ulp = Float_Storage_Half_ULP(stored);
    result.total_error_estimate =
        integration_error_estimate + result.storage_error;
    const double pre_storage_magnitude = Abs(pre_storage_value);
    const double integration_allowed_error =
        CONTRACTED_ERROR_TOLERANCE * (1.0 + pre_storage_magnitude);
    result.allowed_error =
        integration_allowed_error + result.storage_half_ulp;
    result.accepted =
        Contracted_Error_Is_Acceptable(pre_storage_magnitude,
                                       integration_error_estimate) &&
        Double_Is_Finite(result.storage_error) &&
        Double_Is_Finite(result.storage_half_ulp) &&
        result.storage_error <= result.storage_half_ulp &&
        Double_Is_Finite(result.total_error_estimate) &&
        Double_Is_Finite(result.allowed_error) &&
        result.total_error_estimate <= result.allowed_error;
    return result;
}

struct Gradient_Assessment
{
    bool accepted = false;
    double signed_observable = 0.0;
    double absolute_error_estimate = DBL_MAX;
};

// Relative scaling is based on the final signed ECP observable. An L1 sum of
// contribution magnitudes would incorrectly relax the absolute tolerance when
// large positive and negative terms cancel.
QC_ECP_ERROR_HOST_DEVICE inline Gradient_Assessment Assess_Gradient(
    double signed_observable, double absolute_error_estimate)
{
    Gradient_Assessment result;
    result.signed_observable = signed_observable;
    result.absolute_error_estimate = absolute_error_estimate;
    result.accepted =
        Double_Is_Finite(signed_observable) &&
        Contracted_Error_Is_Acceptable(Abs(signed_observable),
                                       absolute_error_estimate);
    return result;
}

// Triangle-inequality propagation for
//   c * (2 alpha V(l+1) - l V(l-1)).
QC_ECP_ERROR_HOST_DEVICE inline double Raising_Lowering_Error_Estimate(
    double absolute_contraction, double gaussian_exponent, int angular_momentum,
    double raised_error_estimate, double lowered_error_estimate)
{
    return Abs(absolute_contraction) *
           (2.0 * Abs(gaussian_exponent) * raised_error_estimate +
            static_cast<double>(angular_momentum) * lowered_error_estimate);
}
}  // namespace qc_ecp_error_policy

#undef QC_ECP_ERROR_HOST_DEVICE
