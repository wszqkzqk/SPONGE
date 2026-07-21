#define NO_GLOBAL_CONTROLLER

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

// eri_common.hpp relies on the integral types and device compatibility macros
// pulled in by ri_workspace.hpp, so this order is intentional.
// clang-format off
#include "quantum_chemistry/integrals/ri/ri_workspace.hpp"
#include "quantum_chemistry/integrals/eri/common/eri_common.hpp"
// clang-format on

static void Add_Neighborhood(std::vector<float>& points, float center)
{
    points.push_back(std::nextafter(center, 0.0f));
    points.push_back(center);
    points.push_back(
        std::nextafter(center, std::numeric_limits<float>::infinity()));
}

int main()
{
    std::vector<float> points = {
        0.0f,     std::nextafter(0.0f, 1.0f),
        1.0e-20f, 1.0e-12f,
        1.0e-7f,  1.0e-2f,
        1.0e-1f,  1.0f,
        2.0f,     4.0f,
    };

    // The first four neighborhoods exercise error-controlled transitions for
    // different requested orders.  30 is the old fixed Miller/upward branch
    // boundary, where the previous implementations were discontinuous.
    Add_Neighborhood(points, 7.44f);
    Add_Neighborhood(points, 11.57f);
    Add_Neighborhood(points, 13.65f);
    Add_Neighborhood(points, 17.92f);
    points.push_back(20.0f);
    points.push_back(28.003f);  // worst observed point for max_m + 25 Miller
    Add_Neighborhood(points, 30.0f);
    points.push_back(64.0f);
    points.push_back(100.0f);
    points.push_back(1000.0f);

    for (float t : points)
    {
        std::uint32_t bits = 0;
        std::memcpy(&bits, &t, sizeof(bits));
        // The requested regression range is m=0..15.  m=16 is included as
        // well because the highest angular-momentum ERI kernel requests it.
        for (int max_m = 0; max_m <= 16; max_m++)
        {
            double one_e[17] = {};
            double ri[17] = {};
            double eri[17] = {};
            compute_boys_double(one_e, t, max_m);
            QC_RI_Compute_Boys_Double(ri, t, max_m);
            eri_boys(eri, t, max_m);
            for (int order = 0; order <= max_m; ++order)
                std::printf("%08x %d %d %.17e %.17e %.17e\n", bits,
                            max_m, order, one_e[order], ri[order],
                            eri[order]);
        }
    }
    return 0;
}
