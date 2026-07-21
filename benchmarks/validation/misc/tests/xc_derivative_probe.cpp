#define NO_GLOBAL_CONTROLLER

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>

#include "quantum_chemistry/dft/xc.hpp"

static double RKS_Energy(QC_METHOD method, double rho, double sigma)
{
    return QC_Local_Exc_Density(method, rho, sigma);
}

static double UKS_Energy(QC_METHOD method, const std::array<double, 5>& x)
{
    return QC_Local_Exc_Density_UKS(method, x[0], x[1], x[2], x[3], x[4]);
}

template <typename Function>
static double Five_Point_Derivative(Function function, double x, double step)
{
    return (-function(x + 2.0 * step) + 8.0 * function(x + step) -
            8.0 * function(x - step) + function(x - 2.0 * step)) /
           (12.0 * step);
}

static std::uint64_t Double_Bits(double value)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static std::uint64_t Next_Random(std::uint64_t& state)
{
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
}

static std::uint64_t Random_Positive_Finite_Bits(std::uint64_t& state)
{
    const std::uint64_t fraction =
        Next_Random(state) & UINT64_C(0x000fffffffffffff);
    const std::uint64_t exponent = Next_Random(state) % 2047u;
    const std::uint64_t bits = (exponent << 52) | fraction;
    return bits == 0u ? 1u : bits;
}

static double Random_Bounded_Finite(std::uint64_t& state)
{
    const std::uint64_t fraction =
        Next_Random(state) & UINT64_C(0x000fffffffffffff);
    // Unbiased exponents in [-400, 400] keep every exact three-vector norm
    // representable while still exercising a very wide dynamic range.
    const std::uint64_t exponent = 623u + Next_Random(state) % 801u;
    const std::uint64_t sign = Next_Random(state) & (UINT64_C(1) << 63);
    return QC_XC_Double_From_Bits(sign | (exponent << 52) | fraction);
}

int main()
{
    const QC_METHOD methods[5] = {QC_METHOD::LDA, QC_METHOD::PBE,
                                  QC_METHOD::BLYP, QC_METHOD::PBE0,
                                  QC_METHOD::B3LYP};

    const double rho = 0.7;
    const double sigma = 0.04;
    for (QC_METHOD method : methods)
    {
        double energy, vrho, vsigma;
        QC_VXC_Analytical_RKS(method, rho, sigma, energy, vrho, vsigma);
        const double vrho_reference = Five_Point_Derivative(
            [&](double value) { return RKS_Energy(method, value, sigma); },
            rho, 2.0e-4);
        const double vsigma_reference = Five_Point_Derivative(
            [&](double value) { return RKS_Energy(method, rho, value); },
            sigma, 1.0e-4);
        std::printf("rks %d %.17e %.17e %.17e %.17e %.17e\n",
                    static_cast<int>(method), energy, vrho, vrho_reference,
                    vsigma, vsigma_reference);

        const std::array<double, 5> point = {0.42, 0.28, 0.05, 0.01, 0.03};
        double uenergy, vra, vrb, vsaa, vsab, vsbb;
        QC_VXC_Analytical_UKS(method, point[0], point[1], point[2], point[3],
                              point[4], uenergy, vra, vrb, vsaa, vsab, vsbb);
        const double analytical[5] = {vra, vrb, vsaa, vsab, vsbb};
        std::printf("uks %d %.17e", static_cast<int>(method), uenergy);
        for (int component = 0; component < 5; ++component)
        {
            const double step = component < 2 ? 1.0e-4 : 5.0e-5;
            const double reference = Five_Point_Derivative(
                [&](double value)
                {
                    std::array<double, 5> varied = point;
                    varied[component] = value;
                    return UKS_Energy(method, varied);
                },
                point[component], step);
            std::printf(" %.17e %.17e", analytical[component], reference);
        }
        std::printf("\n");

        double closed_energy, closed_vra, closed_vrb, closed_vsaa,
            closed_vsab, closed_vsbb;
        QC_VXC_Analytical_UKS(method, 0.5 * rho, 0.5 * rho, 0.25 * sigma,
                              0.25 * sigma, 0.25 * sigma, closed_energy,
                              closed_vra, closed_vrb, closed_vsaa,
                              closed_vsab, closed_vsbb);
        std::printf("closed %d %.17e %.17e %.17e %.17e %.17e %.17e\n",
                    static_cast<int>(method), energy, closed_energy, vrho,
                    0.5 * (closed_vra + closed_vrb), vsigma,
                    0.25 * (closed_vsaa + closed_vsab + closed_vsbb));

        double swapped_energy, swapped_vra, swapped_vrb, swapped_vsaa,
            swapped_vsab, swapped_vsbb;
        QC_VXC_Analytical_UKS(method, point[1], point[0], point[4], point[3],
                              point[2], swapped_energy, swapped_vra,
                              swapped_vrb, swapped_vsaa, swapped_vsab,
                              swapped_vsbb);
        std::printf("swap %d %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                    static_cast<int>(method), uenergy, swapped_energy, vra,
                    swapped_vrb, vrb, swapped_vra,
                    vsaa - swapped_vsbb + vsbb - swapped_vsaa +
                        vsab - swapped_vsab);
    }

    const double tail_density[] = {1.0e-30, 1.0e-60, 1.0e-120,
                                   1.0e-180, 1.0e-240, 1.0e-300};
    for (double tail_rho : tail_density)
    {
        for (QC_METHOD method : methods)
        {
            double energy, vrho, vsigma;
            QC_VXC_Analytical_RKS(method, tail_rho, tail_rho, energy, vrho,
                                  vsigma);
            std::printf("tail %d %.1e %.17e %.17e %.17e %d%d%d\n",
                        static_cast<int>(method), tail_rho, energy, vrho,
                        vsigma, QC_XC_Double_Is_Finite(energy) ? 1 : 0,
                        QC_XC_Double_Is_Finite(vrho) ? 1 : 0,
                        QC_XC_Double_Is_Finite(vsigma) ? 1 : 0);
        }

        const QC_XC_Dual<1> tail_dual =
            QC_XC_Dual<1>::Variable(tail_rho, 0);
        const QC_XC_Dual<1> vwn = QC_XC_VWN5_Impl(tail_dual);
        std::printf("vwn-tail %.1e %.17e %.17e\n", tail_rho, vwn.value,
                    vwn.deriv[0]);
    }

    for (QC_METHOD method : {QC_METHOD::PBE, QC_METHOD::PBE0})
    {
        double energy, vra, vrb, vsaa, vsab, vsbb;
        QC_VXC_Analytical_UKS(method, 0.7, 0.0, 0.04, 0.0, 0.0, energy,
                              vra, vrb, vsaa, vsab, vsbb);
        const unsigned mask = QC_XC_UKS_Expected_Infinite_Output_Mask(
            method, 0.7, 0.0, 0.04, 0.0, 0.0);
        std::printf("endpoint %d %u %.17e %.17e %016llx %.17e %.17e %.17e\n",
                    static_cast<int>(method), mask, energy, vra,
                    static_cast<unsigned long long>(Double_Bits(vrb)), vsaa,
                    vsab, vsbb);

        QC_VXC_Analytical_UKS(method, 0.7, 0.0, 0.0, 0.0, 0.0, energy, vra,
                              vrb, vsaa, vsab, vsbb);
        std::printf("endpoint-zero-gradient %d %u %.17e %.17e %.17e %.17e "
                    "%.17e %.17e\n",
                    static_cast<int>(method),
                    QC_XC_UKS_Expected_Infinite_Output_Mask(
                        method, 0.7, 0.0, 0.0, 0.0, 0.0),
                    energy, vra, vrb, vsaa, vsab, vsbb);

        const double smallest_sigma =
            std::numeric_limits<double>::denorm_min();
        double zero_energy, zero_vrho, zero_vsigma;
        double tiny_energy, tiny_vrho, tiny_vsigma;
        QC_VXC_Analytical_RKS(method, 0.7, 0.0, zero_energy, zero_vrho,
                              zero_vsigma);
        QC_VXC_Analytical_RKS(method, 0.7, smallest_sigma, tiny_energy,
                              tiny_vrho, tiny_vsigma);
        std::printf("tiny-sigma %d %.17e %.17e %.17e %.17e %.17e %.17e\n",
                    static_cast<int>(method), zero_energy, tiny_energy,
                    zero_vrho, tiny_vrho, zero_vsigma, tiny_vsigma);

        double u_zero_energy, u_zero_vra, u_zero_vrb, u_zero_vsaa,
            u_zero_vsab, u_zero_vsbb;
        double u_tiny_energy, u_tiny_vra, u_tiny_vrb, u_tiny_vsaa,
            u_tiny_vsab, u_tiny_vsbb;
        QC_VXC_Analytical_UKS(method, 0.42, 0.28, 0.0, 0.0, 0.0,
                              u_zero_energy, u_zero_vra, u_zero_vrb,
                              u_zero_vsaa, u_zero_vsab, u_zero_vsbb);
        QC_VXC_Analytical_UKS(method, 0.42, 0.28, smallest_sigma, 0.0, 0.0,
                              u_tiny_energy, u_tiny_vra, u_tiny_vrb,
                              u_tiny_vsaa, u_tiny_vsab, u_tiny_vsbb);
        std::printf("tiny-sigma-uks %d %.17e %.17e %.17e %.17e %.17e "
                    "%.17e %.17e %.17e\n",
                    static_cast<int>(method), u_zero_energy, u_tiny_energy,
                    u_zero_vra, u_tiny_vra, u_zero_vrb, u_tiny_vrb,
                    u_zero_vsaa, u_tiny_vsaa);

        for (double zero_tail_rho : {1.0e-235, 1.0e-300})
        {
            double zero_tail_energy, zero_tail_vrho, zero_tail_vsigma;
            QC_VXC_Analytical_RKS(method, zero_tail_rho, 0.0,
                                  zero_tail_energy, zero_tail_vrho,
                                  zero_tail_vsigma);
            std::printf("zero-tail %d %.1e %016llx %016llx %016llx\n",
                        static_cast<int>(method), zero_tail_rho,
                        static_cast<unsigned long long>(
                            Double_Bits(zero_tail_energy)),
                        static_cast<unsigned long long>(
                            Double_Bits(zero_tail_vrho)),
                        static_cast<unsigned long long>(
                            Double_Bits(zero_tail_vsigma)));
        }
    }

    const double nan = QC_XC_Quiet_NaN();
    const double infinity = QC_XC_Positive_Infinity();
    std::printf("valid %d%d%d%d%d%d%d%d\n",
                QC_XC_Double_Is_Finite(1.0) ? 1 : 0,
                QC_XC_Double_Is_Finite(nan) ? 1 : 0,
                QC_XC_Double_Is_Finite(infinity) ? 1 : 0,
                QC_XC_Inputs_Valid_RKS(1.0e-300, 0.0) ? 1 : 0,
                QC_XC_Inputs_Valid_RKS(-1.0, 0.0) ? 1 : 0,
                QC_XC_Inputs_Valid_RKS(0.0, 1.0) ? 1 : 0,
                QC_XC_Inputs_Valid_UKS(0.7, 0.0, 0.04, 0.0, 0.0) ? 1 : 0,
                QC_XC_Inputs_Valid_UKS(0.7, 0.0, 0.04, 1.0e-9, 0.0) ? 1
                                                                          : 0);

    const double gram_diagonal = std::numeric_limits<double>::max() / 8.0;
    const double outside = std::nextafter(2.0, 3.0);
    std::printf("gram %d%d%d%d%d %016llx\n",
                QC_XC_Inputs_Valid_UKS(0.5, 0.5, 4.0, 2.0, 1.0) ? 1 : 0,
                QC_XC_Inputs_Valid_UKS(0.5, 0.5, 4.0, -2.0, 1.0) ? 1 : 0,
                QC_XC_Inputs_Valid_UKS(0.5, 0.5, 4.0, outside, 1.0) ? 1 : 0,
                QC_XC_Inputs_Valid_UKS(0.5, 0.5, 4.0, -outside, 1.0) ? 1 : 0,
                QC_XC_Inputs_Valid_UKS(0.5, 0.5, gram_diagonal,
                                       gram_diagonal, gram_diagonal)
                    ? 1
                    : 0,
                static_cast<unsigned long long>(Double_Bits(outside)));

    const double gram_saa = 1.3748930833837857e-158;
    const double gram_sbb = 2.27790958946488e-246;
    const double gram_invalid =
        QC_XC_Double_From_Bits(UINT64_C(0x160bbe25fe73ee69));
    const double gram_adjacent_valid =
        QC_XC_Double_From_Bits(UINT64_C(0x160bbe25fe73ee68));
    std::printf(
        "gram-counterexample %d%d%d%d\n",
        QC_XC_Inputs_Valid_UKS(0.5, 0.5, gram_saa,
                               gram_adjacent_valid, gram_sbb)
            ? 1
            : 0,
        QC_XC_Inputs_Valid_UKS(0.5, 0.5, gram_saa,
                               -gram_adjacent_valid, gram_sbb)
            ? 1
            : 0,
        QC_XC_Inputs_Valid_UKS(0.5, 0.5, gram_saa, gram_invalid, gram_sbb)
            ? 1
            : 0,
        QC_XC_Inputs_Valid_UKS(0.5, 0.5, gram_saa, -gram_invalid, gram_sbb)
            ? 1
            : 0);

    std::uint64_t random_state = UINT64_C(0x83a5d731c9472e6b);
    for (int sample = 0; sample < 256; ++sample)
    {
        const std::uint64_t aa_bits =
            Random_Positive_Finite_Bits(random_state);
        const std::uint64_t ab_bits =
            Random_Positive_Finite_Bits(random_state);
        const std::uint64_t bb_bits =
            Random_Positive_Finite_Bits(random_state);
        const double aa = QC_XC_Double_From_Bits(aa_bits);
        const double ab = QC_XC_Double_From_Bits(ab_bits);
        const double bb = QC_XC_Double_From_Bits(bb_bits);
        std::printf("gram-exact %016llx %016llx %016llx %d\n",
                    static_cast<unsigned long long>(aa_bits),
                    static_cast<unsigned long long>(ab_bits),
                    static_cast<unsigned long long>(bb_bits),
                    QC_XC_Gradient_Gram_Is_PSD(aa, ab, bb) ? 1 : 0);
    }

    for (int sample = 0; sample < 256; ++sample)
    {
        double ga[3] = {Random_Bounded_Finite(random_state),
                        Random_Bounded_Finite(random_state),
                        Random_Bounded_Finite(random_state)};
        double gb[3] = {Random_Bounded_Finite(random_state),
                        Random_Bounded_Finite(random_state),
                        Random_Bounded_Finite(random_state)};
        if (sample % 16 == 0)
            for (int dir = 0; dir < 3; ++dir) gb[dir] = ga[dir];
        else if (sample % 16 == 1)
            for (int dir = 0; dir < 3; ++dir) gb[dir] = -ga[dir];
        else if (sample % 16 == 2)
            ga[0] = ga[1] = ga[2] = 0.0;
        else if (sample % 16 == 3)
            gb[0] = gb[1] = gb[2] = 0.0;

        double saa = 0.0, sab = 0.0, sbb = 0.0;
        const bool built = QC_XC_Build_Gradient_Gram(
            ga[0], ga[1], ga[2], gb[0], gb[1], gb[2], saa, sab, sbb);
        const bool accepted =
            built && QC_XC_Inputs_Valid_UKS(0.5, 0.5, saa, sab, sbb);
        std::printf(
            "gram-builder %016llx %016llx %016llx %016llx %016llx "
            "%016llx %d%d %016llx %016llx %016llx\n",
            static_cast<unsigned long long>(Double_Bits(ga[0])),
            static_cast<unsigned long long>(Double_Bits(ga[1])),
            static_cast<unsigned long long>(Double_Bits(ga[2])),
            static_cast<unsigned long long>(Double_Bits(gb[0])),
            static_cast<unsigned long long>(Double_Bits(gb[1])),
            static_cast<unsigned long long>(Double_Bits(gb[2])), built ? 1 : 0,
            accepted ? 1 : 0,
            static_cast<unsigned long long>(Double_Bits(saa)),
            static_cast<unsigned long long>(Double_Bits(sab)),
            static_cast<unsigned long long>(Double_Bits(sbb)));
    }

    for (QC_METHOD method : methods)
    {
        double energy, vrho, vsigma;
        QC_VXC_Analytical_RKS(method, 0.0, 0.0, energy, vrho, vsigma);
        std::printf("vacuum %d %.17e %.17e %.17e\n",
                    static_cast<int>(method), energy, vrho, vsigma);
    }
    return 0;
}
