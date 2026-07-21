#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "cart2sph.hpp"
#include "cartesian_components.hpp"

namespace qc_analytical_norms
{
inline bool Float_Is_Positive_Normal(float value)
{
    std::uint32_t bits = 0U;
    std::memcpy(&bits, &value, sizeof(bits));
    const std::uint32_t exponent = bits & 0x7f800000U;
    return (bits & 0x80000000U) == 0U && exponent != 0U &&
           exponent != 0x7f800000U;
}

inline bool Double_Is_Positive_Finite(double value)
{
    std::uint64_t bits = 0ULL;
    static_assert(sizeof(bits) == sizeof(value), "unexpected double size");
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x8000000000000000ULL) == 0ULL &&
           (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL &&
           (bits & 0x7fffffffffffffffULL) != 0ULL;
}

inline double Odd_Double_Factorial(int power)
{
    if (power < 0)
        throw std::invalid_argument(
            "double-factorial Cartesian power must be non-negative");
    double result = 1.0;
    for (int k = 1; k <= power; ++k)
    {
        result *= static_cast<double>(2 * k - 1);
        if (!Double_Is_Positive_Finite(result))
            throw std::overflow_error("Cartesian double factorial overflow");
    }
    return result;
}

inline std::vector<float> Build(
    const std::vector<int>& l_list, const std::vector<int>& shell_sizes,
    const std::vector<int>& shell_offsets, const std::vector<float>& exponents,
    const std::vector<float>& contraction_coefficients,
    const std::vector<int>& ao_offsets_cart,
    const std::vector<int>& ao_offsets_effective, bool is_spherical,
    const std::vector<float>& cart2sph, int nao_cart, int nao_effective)
{
    const std::size_t nbas = l_list.size();
    if (shell_sizes.size() != nbas || shell_offsets.size() != nbas ||
        ao_offsets_cart.size() != nbas ||
        ao_offsets_effective.size() != nbas)
        throw std::invalid_argument(
            "analytical-norm shell metadata lengths are inconsistent");
    if (exponents.size() != contraction_coefficients.size())
        throw std::invalid_argument(
            "analytical-norm exponent/coefficient lengths differ");
    if (nao_cart < 0 || nao_effective < 0)
        throw std::invalid_argument(
            "analytical-norm AO dimensions must be non-negative");

    std::size_t expected_cart = 0U;
    std::size_t expected_effective = 0U;
    for (std::size_t shell = 0; shell < nbas; ++shell)
    {
        const std::pair<int, int> dimensions =
            qc_cart2sph::Int_Dimensions(l_list[shell]);
        const std::size_t dim_cart =
            static_cast<std::size_t>(dimensions.first);
        const std::size_t dim_effective = static_cast<std::size_t>(
            is_spherical ? dimensions.second : dimensions.first);
        if (expected_cart >
                static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                    dim_cart ||
            expected_effective >
                static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                    dim_effective)
            throw std::overflow_error("analytical-norm AO offset overflow");
        if (ao_offsets_cart[shell] != static_cast<int>(expected_cart) ||
            ao_offsets_effective[shell] !=
                static_cast<int>(expected_effective))
            throw std::invalid_argument(
                "analytical-norm AO offsets do not match shell order");
        expected_cart += dim_cart;
        expected_effective += dim_effective;
    }
    if (expected_cart != static_cast<std::size_t>(nao_cart) ||
        expected_effective != static_cast<std::size_t>(nao_effective))
        throw std::invalid_argument(
            "analytical-norm AO dimensions do not match shell metadata");
    if (is_spherical)
    {
        if (expected_effective != 0U &&
            expected_cart > std::numeric_limits<std::size_t>::max() /
                                expected_effective)
            throw std::overflow_error(
                "analytical-norm Cartesian transform size overflow");
        if (cart2sph.size() != expected_cart * expected_effective)
            throw std::invalid_argument(
                "analytical-norm Cartesian transform has wrong size");
    }

    std::vector<float> norms(expected_effective, 0.0f);
    const double pi = std::acos(-1.0);
    for (std::size_t shell = 0; shell < nbas; ++shell)
    {
        const int l = l_list[shell];
        const std::pair<int, int> dimensions =
            qc_cart2sph::Int_Dimensions(l);
        const int dim_cart = dimensions.first;
        const int dim_effective =
            is_spherical ? dimensions.second : dimensions.first;
        const int primitive_count = shell_sizes[shell];
        const int primitive_offset = shell_offsets[shell];
        if (primitive_count <= 0 || primitive_offset < 0 ||
            static_cast<std::size_t>(primitive_offset) > exponents.size() ||
            static_cast<std::size_t>(primitive_count) >
                exponents.size() -
                    static_cast<std::size_t>(primitive_offset))
            throw std::invalid_argument(
                "analytical-norm primitive range is invalid at shell " +
                std::to_string(shell));

        double radial_overlap = 0.0;
        for (int p = 0; p < primitive_count; ++p)
        {
            const float alpha_f = exponents[primitive_offset + p];
            const float coeff_p_f =
                contraction_coefficients[primitive_offset + p];
            if (!Float_Is_Positive_Normal(alpha_f) ||
                !qc_cart2sph::Float_Is_Zero_Or_Normal(coeff_p_f))
                throw std::runtime_error(
                    "invalid exponent/coefficient float at shell " +
                    std::to_string(shell));
            const double alpha = static_cast<double>(alpha_f);
            const double coeff_p = static_cast<double>(coeff_p_f);
            for (int q = 0; q < primitive_count; ++q)
            {
                const float beta_f = exponents[primitive_offset + q];
                const float coeff_q_f =
                    contraction_coefficients[primitive_offset + q];
                if (!Float_Is_Positive_Normal(beta_f) ||
                    !qc_cart2sph::Float_Is_Zero_Or_Normal(coeff_q_f))
                    throw std::runtime_error(
                        "invalid exponent/coefficient float at shell " +
                        std::to_string(shell));
                const double gamma = alpha + static_cast<double>(beta_f);
                double term =
                    coeff_p * static_cast<double>(coeff_q_f) *
                    std::pow(pi / gamma, 1.5);
                for (int degree = 0; degree < l; ++degree)
                    term *= 0.5 / gamma;
                radial_overlap += term;
            }
        }
        if (!Double_Is_Positive_Finite(radial_overlap))
            throw std::runtime_error(
                "contracted radial self-overlap is not finite and positive "
                "at shell " +
                std::to_string(shell));

        const int cart_offset = ao_offsets_cart[shell];
        const int effective_offset = ao_offsets_effective[shell];
        for (int component = 0; component < dim_effective; ++component)
        {
            double angular_overlap = 0.0;
            if (!is_spherical)
            {
                int lx = 0, ly = 0, lz = 0;
                if (!qc_cartesian::Component(l, component, lx, ly, lz))
                    throw std::runtime_error(
                        "invalid Cartesian component in analytical norm");
                angular_overlap = Odd_Double_Factorial(lx) *
                                  Odd_Double_Factorial(ly) *
                                  Odd_Double_Factorial(lz);
            }
            else
            {
                for (int cart_a = 0; cart_a < dim_cart; ++cart_a)
                {
                    const float transform_a =
                        cart2sph[(static_cast<std::size_t>(cart_offset) +
                                  static_cast<std::size_t>(cart_a)) *
                                     expected_effective +
                                 static_cast<std::size_t>(effective_offset +
                                                          component)];
                    if (!qc_cart2sph::Float_Is_Zero_Or_Normal(transform_a))
                        throw std::runtime_error(
                            "Cartesian-to-spherical coefficient is neither "
                            "zero nor normal");
                    if (transform_a == 0.0f) continue;
                    int lx_a = 0, ly_a = 0, lz_a = 0;
                    if (!qc_cartesian::Component(l, cart_a, lx_a, ly_a,
                                                 lz_a))
                        throw std::runtime_error(
                            "invalid Cartesian component in analytical norm");
                    for (int cart_b = 0; cart_b < dim_cart; ++cart_b)
                    {
                        const float transform_b =
                            cart2sph
                                [(static_cast<std::size_t>(cart_offset) +
                                  static_cast<std::size_t>(cart_b)) *
                                     expected_effective +
                                 static_cast<std::size_t>(effective_offset +
                                                          component)];
                        if (!qc_cart2sph::Float_Is_Zero_Or_Normal(transform_b))
                            throw std::runtime_error(
                                "Cartesian-to-spherical coefficient is "
                                "neither zero nor normal");
                        if (transform_b == 0.0f) continue;
                        int lx_b = 0, ly_b = 0, lz_b = 0;
                        if (!qc_cartesian::Component(l, cart_b, lx_b, ly_b,
                                                     lz_b))
                            throw std::runtime_error(
                                "invalid Cartesian component in analytical "
                                "norm");
                        const int x_degree = lx_a + lx_b;
                        const int y_degree = ly_a + ly_b;
                        const int z_degree = lz_a + lz_b;
                        if ((x_degree & 1) != 0 || (y_degree & 1) != 0 ||
                            (z_degree & 1) != 0)
                            continue;
                        angular_overlap +=
                            static_cast<double>(transform_a) *
                            static_cast<double>(transform_b) *
                            Odd_Double_Factorial(x_degree / 2) *
                            Odd_Double_Factorial(y_degree / 2) *
                            Odd_Double_Factorial(z_degree / 2);
                    }
                }
            }
            if (!Double_Is_Positive_Finite(angular_overlap))
                throw std::runtime_error(
                    "angular self-overlap is not finite and positive at shell " +
                    std::to_string(shell) + ", component " +
                    std::to_string(component));
            const double diagonal = radial_overlap * angular_overlap;
            if (!Double_Is_Positive_Finite(diagonal))
                throw std::runtime_error(
                    "AO self-overlap is not finite and positive at shell " +
                    std::to_string(shell) + ", component " +
                    std::to_string(component));
            const float norm = static_cast<float>(1.0 / std::sqrt(diagonal));
            if (!Float_Is_Positive_Normal(norm))
                throw std::runtime_error(
                    "AO normalization is not a positive normal float at shell " +
                    std::to_string(shell) + ", component " +
                    std::to_string(component));
            norms[static_cast<std::size_t>(effective_offset + component)] =
                norm;
        }
    }
    for (float norm : norms)
        if (!Float_Is_Positive_Normal(norm))
            throw std::runtime_error(
                "analytical-norm output contains an uninitialized value");
    return norms;
}
}  // namespace qc_analytical_norms
