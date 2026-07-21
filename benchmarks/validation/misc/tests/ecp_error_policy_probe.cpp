#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "quantum_chemistry/ecp/ecp_error_policy.hpp"

namespace
{
double Double_From_Bits(std::uint64_t bits)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

bool Nearly_Equal(double left, double right, double tolerance = 1.0e-15)
{
    return std::fabs(left - right) <=
           tolerance * (1.0 + std::fabs(right));
}
}  // namespace

int main()
{
    using namespace qc_ecp_error_policy;

    // Halfway between adjacent floats at 1.0: the integration estimate is
    // zero, but the final observable still carries a 2^-24 storage error.
    const double halfway = 1.0 + std::ldexp(1.0, -24);
    const Matrix_Storage_Assessment storage =
        Assess_Matrix_Storage(halfway, 0.0);
    if (!storage.accepted || storage.stored_value != 1.0 ||
        !Nearly_Equal(storage.storage_error, std::ldexp(1.0, -24), 0.0) ||
        !Nearly_Equal(storage.storage_half_ulp, std::ldexp(1.0, -24), 0.0) ||
        storage.storage_error > storage.storage_half_ulp ||
        !Nearly_Equal(storage.total_error_estimate, storage.storage_error,
                      0.0))
    {
        std::fprintf(stderr, "float storage rounding was not accounted for\n");
        return EXIT_FAILURE;
    }

    // A relative epsilon expression is exact only at a binade's lower edge.
    // Near 2.0 it approaches one full ulp, whereas the representation budget
    // remains exactly 2^-24 throughout [1,2).  Normal/subnormal transition and
    // zero all share binary32 spacing 2^-149.
    const float below_two = std::nextafter(2.0f, 0.0f);
    const Matrix_Storage_Assessment one = Assess_Matrix_Storage(1.0, 0.0);
    const Matrix_Storage_Assessment upper_binade =
        Assess_Matrix_Storage(static_cast<double>(below_two), 0.0);
    const Matrix_Storage_Assessment minimum_normal = Assess_Matrix_Storage(
        static_cast<double>(std::numeric_limits<float>::min()), 0.0);
    const Matrix_Storage_Assessment minimum_subnormal = Assess_Matrix_Storage(
        static_cast<double>(std::numeric_limits<float>::denorm_min()), 0.0);
    const Matrix_Storage_Assessment zero = Assess_Matrix_Storage(0.0, 0.0);
    if (!one.accepted || !upper_binade.accepted || !minimum_normal.accepted ||
        !minimum_subnormal.accepted || !zero.accepted ||
        one.storage_half_ulp != std::ldexp(1.0, -24) ||
        upper_binade.storage_half_ulp != std::ldexp(1.0, -24) ||
        minimum_normal.storage_half_ulp != std::ldexp(1.0, -150) ||
        minimum_subnormal.storage_half_ulp != std::ldexp(1.0, -150) ||
        zero.storage_half_ulp != std::ldexp(1.0, -150))
    {
        std::fprintf(stderr, "float half-ulp binade policy is incorrect\n");
        return EXIT_FAILURE;
    }

    const Matrix_Storage_Assessment minimum_normal_halfway =
        Assess_Matrix_Storage(std::ldexp(1.0, -126) -
                                  std::ldexp(1.0, -150),
                              0.0);
    const Matrix_Storage_Assessment even_subnormal_halfway =
        Assess_Matrix_Storage(3.0 * std::ldexp(1.0, -150), 0.0);
    const Matrix_Storage_Assessment zero_halfway =
        Assess_Matrix_Storage(std::ldexp(1.0, -150), 0.0);
    if (!minimum_normal_halfway.accepted ||
        !even_subnormal_halfway.accepted || !zero_halfway.accepted ||
        minimum_normal_halfway.storage_error !=
            minimum_normal_halfway.storage_half_ulp ||
        even_subnormal_halfway.storage_error !=
            even_subnormal_halfway.storage_half_ulp ||
        zero_halfway.storage_error != zero_halfway.storage_half_ulp)
    {
        std::fprintf(stderr,
                     "float halfway cast exceeded its exact half-ulp\n");
        return EXIT_FAILURE;
    }

    // The half-ulp is an additional representation allowance, not a way to
    // relax the integration certificate itself.  This estimate is only 1%
    // above the original integration budget and remains below budget+half-ulp.
    const double integration_budget =
        CONTRACTED_ERROR_TOLERANCE * (1.0 + 1.0);
    const Matrix_Storage_Assessment excessive_integration =
        Assess_Matrix_Storage(1.0, 1.01 * integration_budget);
    if (excessive_integration.accepted ||
        !(excessive_integration.total_error_estimate <
          excessive_integration.allowed_error))
    {
        std::fprintf(stderr,
                     "float allowance relaxed the integration budget\n");
        return EXIT_FAILURE;
    }

    // +M and -M cancel in the signed observable. The old L1 scale (2M) would
    // accept this 0.4 absolute estimate, while the final-observable policy must
    // reject it against the near-zero absolute budget.
    const double large = 1.0e8;
    const double cancelled = large - large;
    const double cancellation_error = 0.4;
    const Gradient_Assessment cancellation =
        Assess_Gradient(cancelled, cancellation_error);
    const bool old_l1_would_accept = Contracted_Error_Is_Acceptable(
        std::fabs(large) + std::fabs(-large), cancellation_error);
    if (cancellation.accepted || !old_l1_would_accept)
    {
        std::fprintf(stderr, "signed cancellation used an L1 relative scale\n");
        return EXIT_FAILURE;
    }

    // A rejected scratch result must not be published. This mirrors the
    // production transaction: the kernel writes only ECP scratch and the
    // caller commits it after Assess_Gradient succeeds.
    double published_gradient = 7.25;
    const double scratch_gradient = 3.5;
    const Gradient_Assessment rejected_scratch =
        Assess_Gradient(scratch_gradient, 1.0);
    if (rejected_scratch.accepted) published_gradient += scratch_gradient;
    if (published_gradient != 7.25)
    {
        std::fprintf(stderr, "rejected ECP scratch polluted the gradient\n");
        return EXIT_FAILURE;
    }

    const double propagated = Raising_Lowering_Error_Estimate(
        1.0e6, 0.8, 3, 2.0e-13, 5.0e-13);
    const double expected =
        1.0e6 * (2.0 * 0.8 * 2.0e-13 + 3.0 * 5.0e-13);
    if (!Nearly_Equal(propagated, expected))
    {
        std::fprintf(stderr, "raising/lowering coefficient was not propagated\n");
        return EXIT_FAILURE;
    }

    // Deterministic accepted-series sample 7682 from the independent
    // series-vs-quadrature sweep. Its primitive estimate straddles the
    // production fallback threshold at roughly 1.00065e6 contraction.
    const double primitive_estimate = 4.9967497894065997e-15;
    if (!Series_Estimate_Is_Small_Enough(true, 1.0e6,
                                        primitive_estimate) ||
        Series_Estimate_Is_Small_Enough(true, 1.01e6,
                                        primitive_estimate))
    {
        std::fprintf(stderr,
                     "coefficient amplification did not select fallback\n");
        return EXIT_FAILURE;
    }

    const double positive_infinity =
        Double_From_Bits(UINT64_C(0x7ff0000000000000));
    const double quiet_nan =
        Double_From_Bits(UINT64_C(0x7ff8000000000042));
    if (Assess_Gradient(positive_infinity, 0.0).accepted ||
        Assess_Gradient(0.0, quiet_nan).accepted ||
        Assess_Matrix_Storage(positive_infinity, 0.0).accepted)
    {
        std::fprintf(stderr, "non-finite certificate input was accepted\n");
        return EXIT_FAILURE;
    }

    std::printf("storage_error=%.17g\n", storage.storage_error);
    std::printf("cancelled_observable=%.17g error=%.17g\n", cancelled,
                cancellation_error);
    std::printf("amplified_error=%.17g\n", propagated);
    std::printf("fallback_threshold_bracket=[1e6,1.01e6]\n");
    return EXIT_SUCCESS;
}
