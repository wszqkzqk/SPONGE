#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "manybody/reaxff/eeq.h"

namespace
{

bool Close(double actual, double expected, double relative, double absolute)
{
    return std::fabs(actual - expected) <=
           absolute +
               relative * std::fmax(std::fabs(actual), std::fabs(expected));
}

double Reference_Force_X(float qi, float qj, float distance, float shield,
                         float cutoff)
{
    const double r = distance;
    const double x = r / cutoff;
    const double x2 = x * x;
    const double x3 = x2 * x;
    const double x4 = x2 * x2;
    const double x5 = x4 * x;
    const double x6 = x5 * x;
    const double x7 = x6 * x;
    const double taper = 20.0 * x7 - 70.0 * x6 + 84.0 * x5 - 35.0 * x4 + 1.0;
    const double dtaper_dr =
        (140.0 * x6 - 420.0 * x5 + 420.0 * x4 - 140.0 * x3) / cutoff;
    const double u = r * r * r + shield;
    const double inv_u_cbrt = 1.0 / std::cbrt(u);
    const double dH_dr =
        static_cast<double>(ReaxFFEEQ::COULOMB_CONSTANT) *
        (dtaper_dr * inv_u_cbrt - taper * r * r * inv_u_cbrt / u);
    return -static_cast<double>(qi) * static_cast<double>(qj) * dH_dr;
}

}  // namespace

int main()
{
    constexpr float qi = 5.0e-11f;
    constexpr float qj = -4.0e-11f;
    constexpr float shield = 1.0f;
    constexpr float cutoff = 4.0f;
    const VECTOR displacement = {1.0f, 0.0f, 0.0f};

    const auto tiny =
        ReaxFFEEQ::Evaluate_Pair_Force(displacement, qi, qj, shield, cutoff);
    const double expected = Reference_Force_X(qi, qj, 1.0f, shield, cutoff);
    if (!tiny.active || !std::isfinite(tiny.force.x) || tiny.force.x == 0.0f ||
        !Close(tiny.force.x, expected, 2.0e-5, 1.0e-30))
    {
        std::fprintf(stderr,
                     "tiny nonzero EEQ charges lost their pair force: "
                     "actual %.9g expected %.17g\n",
                     tiny.force.x, expected);
        return EXIT_FAILURE;
    }

    const auto neutral =
        ReaxFFEEQ::Evaluate_Pair_Force(displacement, 0.0f, qj, shield, cutoff);
    if (neutral.active || neutral.force.x != 0.0f || neutral.force.y != 0.0f ||
        neutral.force.z != 0.0f)
    {
        std::fprintf(stderr, "an exactly neutral EEQ pair remained active\n");
        return EXIT_FAILURE;
    }

    const VECTOR overlap = {0.0f, 0.0f, 0.0f};
    const auto coincident =
        ReaxFFEEQ::Evaluate_Pair_Force(overlap, qi, qj, shield, cutoff);
    if (!coincident.active || coincident.force.x != 0.0f ||
        coincident.force.y != 0.0f || coincident.force.z != 0.0f)
    {
        std::fprintf(stderr,
                     "shielded EEQ exact-overlap zero-force limit failed\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
