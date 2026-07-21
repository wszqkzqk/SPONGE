#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "quantum_chemistry/gradient/ri_metric_response.hpp"

namespace
{

bool Close(double actual, double expected, double tolerance = 1.0e-14)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &actual, sizeof(bits));
    const bool finite = ((bits >> 52U) & 0x7ffU) != 0x7ffU;
    return finite && std::abs(actual - expected) <= tolerance;
}

}  // namespace

int main()
{
    // Two discarded modes (including a negative roundoff eigenvalue), followed
    // by an exactly degenerate retained pair.
    const double eigval[] = {-1.0e-12, 1.0e-12, 4.0, 4.0};
    const int naux = 4;
    const int naux_eff = 2;

    const double degenerate = QC_RI_Truncated_InvSqrt_Loewner(
        eigval, naux, naux_eff, 2, 3);
    if (!Close(degenerate, -0.0625))
    {
        std::fprintf(stderr, "degenerate retained coefficient: %.17g\n",
                     degenerate);
        return EXIT_FAILURE;
    }

    const double discarded = QC_RI_Truncated_InvSqrt_Loewner(
        eigval, naux, naux_eff, 0, 1);
    if (!Close(discarded, 0.0))
    {
        std::fprintf(stderr, "discarded coefficient: %.17g\n", discarded);
        return EXIT_FAILURE;
    }

    const double mixed = QC_RI_Truncated_InvSqrt_Loewner(
        eigval, naux, naux_eff, 1, 2);
    const double mixed_expected = (0.0 - 0.5) / (eigval[1] - eigval[2]);
    if (!Close(mixed, mixed_expected))
    {
        std::fprintf(stderr, "mixed coefficient: %.17g expected %.17g\n", mixed,
                     mixed_expected);
        return EXIT_FAILURE;
    }

    const double mixed_transpose = QC_RI_Truncated_InvSqrt_Loewner(
        eigval, naux, naux_eff, 2, 1);
    if (!Close(mixed_transpose, mixed_expected))
    {
        std::fprintf(stderr, "mixed transpose coefficient: %.17g\n",
                     mixed_transpose);
        return EXIT_FAILURE;
    }

    const double distinct_eigval[] = {1.0, 4.0};
    const double distinct = QC_RI_Truncated_InvSqrt_Loewner(
        distinct_eigval, 2, 2, 0, 1);
    if (!Close(distinct, -1.0 / 6.0))
    {
        std::fprintf(stderr, "distinct retained coefficient: %.17g\n", distinct);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
