#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Host-side construction of normalized real solid harmonics in the Cartesian
// monomial basis used by the integral kernels.  For l >= 2, spherical columns
// are ordered m=-l,...,+l.  The p shell retains the program's historical
// (px, py, pz) ordering.
namespace qc_cart2sph
{
inline bool Float_Is_Zero_Or_Normal(float value)
{
    std::uint32_t bits = 0U;
    static_assert(sizeof(bits) == sizeof(value), "unexpected float size");
    std::memcpy(&bits, &value, sizeof(bits));
    const std::uint32_t magnitude = bits & 0x7fffffffU;
    const std::uint32_t exponent = bits & 0x7f800000U;
    return magnitude == 0U ||
           (exponent != 0U && exponent != 0x7f800000U);
}

inline float Checked_Float_Coefficient(long double coefficient)
{
    if (!std::isfinite(coefficient) ||
        std::fabs(coefficient) >
            static_cast<long double>(std::numeric_limits<float>::max()))
        throw std::overflow_error(
            "Cartesian-to-spherical coefficient is not representable as "
            "float");
    const float value = static_cast<float>(coefficient);
    if (coefficient != 0.0L && value == 0.0f)
        throw std::overflow_error(
            "nonzero Cartesian-to-spherical coefficient underflowed to "
            "float zero");
    if (!Float_Is_Zero_Or_Normal(value))
        throw std::overflow_error(
            "Cartesian-to-spherical coefficient is neither zero nor a "
            "normal float");
    return value;
}

inline std::size_t Cartesian_Dimension(int l)
{
    if (l < 0)
        throw std::invalid_argument(
            "Cartesian-to-spherical angular momentum must be non-negative");

    const unsigned long long lp1 =
        static_cast<unsigned long long>(l) + 1ULL;
    const unsigned long long lp2 = lp1 + 1ULL;
    const unsigned long long dimension = lp1 * lp2 / 2ULL;
    if (dimension > std::numeric_limits<std::size_t>::max())
        throw std::overflow_error(
            "Cartesian shell dimension exceeds host address space");
    return static_cast<std::size_t>(dimension);
}

inline std::size_t Spherical_Dimension(int l)
{
    if (l < 0)
        throw std::invalid_argument(
            "Cartesian-to-spherical angular momentum must be non-negative");

    const unsigned long long dimension =
        2ULL * static_cast<unsigned long long>(l) + 1ULL;
    if (dimension > std::numeric_limits<std::size_t>::max())
        throw std::overflow_error(
            "spherical shell dimension exceeds host address space");
    return static_cast<std::size_t>(dimension);
}

inline std::pair<int, int> Int_Dimensions(int l)
{
    const std::size_t dim_cart = Cartesian_Dimension(l);
    const std::size_t dim_sph = Spherical_Dimension(l);
    if (dim_cart > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        dim_sph > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error(
            "shell dimension is not representable as int");
    return {static_cast<int>(dim_cart), static_cast<int>(dim_sph)};
}

inline std::size_t Cartesian_Component_Index(int l, int lx, int ly, int lz)
{
    if (l < 0 || lx < 0 || ly < 0 || lz < 0 || lx > l || ly > l ||
        lz > l || lx != l - ly - lz)
        throw std::invalid_argument("invalid Cartesian angular component");

    // Components are ordered by descending lx and then descending ly:
    // xx, xy, xz, yy, yz, zz for a d shell.
    const std::size_t remaining = static_cast<std::size_t>(l - lx);
    return remaining * (remaining + 1U) / 2U +
           static_cast<std::size_t>(lz);
}

inline long double Log_Factorial(int n)
{
    return std::lgamma(static_cast<long double>(n) + 1.0L);
}

inline long double Exp_Checked(long double value, const char* quantity)
{
    const long double result = std::exp(value);
    if (!std::isfinite(result))
        throw std::overflow_error(std::string(quantity) +
                                  " overflowed while building solid harmonic");
    return result;
}

inline std::vector<float> Build_Block(int l)
{
    const std::size_t dim_cart = Cartesian_Dimension(l);
    const std::size_t dim_sph = Spherical_Dimension(l);
    if (dim_sph != 0U &&
        dim_cart > std::numeric_limits<std::size_t>::max() / dim_sph)
        throw std::overflow_error(
            "Cartesian-to-spherical shell matrix size overflow");

    // Preserve the historical one-ULP rounding of the s coefficient.  Higher
    // legacy shells already round identically to the general construction.
    if (l == 0) return {0.28209479f};

    std::vector<long double> coefficients(dim_cart * dim_sph, 0.0L);
    const long double pi = std::acos(-1.0L);
    const long double log_two = std::log(2.0L);

    for (std::size_t column = 0; column < dim_sph; ++column)
    {
        // p functions have always been exposed as px,py,pz.  Other real
        // spherical shells follow m=-l,...,+l.
        int m = static_cast<int>(column) - l;
        if (l == 1)
        {
            static const int p_order[3] = {1, -1, 0};
            m = p_order[column];
        }
        const int abs_m = (m < 0) ? -m : m;

        long double log_normalization =
            0.5L * (std::log((2.0L * static_cast<long double>(l) + 1.0L) /
                             (4.0L * pi)) +
                    Log_Factorial(l - abs_m) -
                    Log_Factorial(l + abs_m));
        if (abs_m > 0) log_normalization += 0.5L * log_two;
        const long double normalization =
            Exp_Checked(log_normalization, "solid-harmonic normalization");

        const int max_k = (l - abs_m) / 2;
        for (int k = 0; k <= max_k; ++k)
        {
            const long double log_legendre =
                Log_Factorial(2 * l - 2 * k) -
                static_cast<long double>(l) * log_two - Log_Factorial(k) -
                Log_Factorial(l - k) - Log_Factorial(l - abs_m - 2 * k);
            long double legendre =
                Exp_Checked(log_legendre, "associated-Legendre coefficient");
            if ((k & 1) != 0) legendre = -legendre;

            for (int radial_x = 0; radial_x <= k; ++radial_x)
            {
                for (int radial_y = 0; radial_y <= k - radial_x;
                     ++radial_y)
                {
                    const int radial_z = k - radial_x - radial_y;
                    const long double radial_multinomial = Exp_Checked(
                        Log_Factorial(k) - Log_Factorial(radial_x) -
                            Log_Factorial(radial_y) -
                            Log_Factorial(radial_z),
                        "radial multinomial coefficient");

                    for (int y_power = 0; y_power <= abs_m; ++y_power)
                    {
                        int phase = 0;
                        if (m > 0 && (y_power & 1) == 0)
                            phase = ((y_power / 2) & 1) ? -1 : 1;
                        else if (m < 0 && (y_power & 1) != 0)
                            phase = (((y_power - 1) / 2) & 1) ? -1 : 1;
                        else if (m == 0)
                            phase = 1;
                        if (phase == 0) continue;

                        const long double azimuthal_binomial = Exp_Checked(
                            Log_Factorial(abs_m) - Log_Factorial(y_power) -
                                Log_Factorial(abs_m - y_power),
                            "azimuthal binomial coefficient");
                        const int lx =
                            2 * radial_x + abs_m - y_power;
                        const int ly = 2 * radial_y + y_power;
                        const int lz = l - abs_m - 2 * k + 2 * radial_z;
                        const std::size_t row =
                            Cartesian_Component_Index(l, lx, ly, lz);
                        coefficients[row * dim_sph + column] +=
                            static_cast<long double>(phase) * normalization *
                            legendre * radial_multinomial *
                            azimuthal_binomial;
                    }
                }
            }
        }
    }

    std::vector<float> block(coefficients.size());
    std::vector<unsigned char> nonzero_columns(dim_sph, 0U);
    for (std::size_t row = 0; row < dim_cart; ++row)
    {
        for (std::size_t column = 0; column < dim_sph; ++column)
        {
            const long double coefficient =
                coefficients[row * dim_sph + column];
            const float value = Checked_Float_Coefficient(coefficient);
            block[row * dim_sph + column] = value;
            if (value != 0.0f) nonzero_columns[column] = 1U;
        }
    }
    for (std::size_t column = 0; column < dim_sph; ++column)
    {
        if (nonzero_columns[column] == 0U)
            throw std::runtime_error(
                "Cartesian-to-spherical construction produced an empty "
                "spherical component");
    }
    return block;
}

inline std::vector<float> Build_Matrix(const std::vector<int>& l_list,
                                       int nao_cart, int nao_sph)
{
    if (nao_cart < 0 || nao_sph < 0)
        throw std::invalid_argument(
            "Cartesian-to-spherical matrix dimensions must be non-negative");

    std::size_t expected_cart = 0U;
    std::size_t expected_sph = 0U;
    for (int l : l_list)
    {
        const std::size_t dim_cart = Cartesian_Dimension(l);
        const std::size_t dim_sph = Spherical_Dimension(l);
        if (expected_cart >
                std::numeric_limits<std::size_t>::max() - dim_cart ||
            expected_sph >
                std::numeric_limits<std::size_t>::max() - dim_sph)
            throw std::overflow_error(
                "Cartesian-to-spherical basis dimension overflow");
        expected_cart += dim_cart;
        expected_sph += dim_sph;
    }
    if (expected_cart != static_cast<std::size_t>(nao_cart) ||
        expected_sph != static_cast<std::size_t>(nao_sph))
        throw std::invalid_argument(
            "Cartesian-to-spherical dimensions do not match shell list");
    if (expected_sph != 0U &&
        expected_cart >
            std::numeric_limits<std::size_t>::max() / expected_sph)
        throw std::overflow_error(
            "Cartesian-to-spherical matrix allocation size overflow");

    std::vector<float> matrix(expected_cart * expected_sph, 0.0f);
    std::size_t cart_offset = 0U;
    std::size_t sph_offset = 0U;
    for (int l : l_list)
    {
        const std::size_t dim_cart = Cartesian_Dimension(l);
        const std::size_t dim_sph = Spherical_Dimension(l);
        const std::vector<float> block = Build_Block(l);
        for (std::size_t row = 0; row < dim_cart; ++row)
            for (std::size_t column = 0; column < dim_sph; ++column)
                matrix[(cart_offset + row) * expected_sph + sph_offset +
                       column] = block[row * dim_sph + column];
        cart_offset += dim_cart;
        sph_offset += dim_sph;
    }
    if (cart_offset != expected_cart || sph_offset != expected_sph)
        throw std::runtime_error(
            "Cartesian-to-spherical block offsets are inconsistent");
    return matrix;
}

// Contract a matrix expressed in the normalized spherical AO basis with the
// Cartesian-to-spherical transform used by the integral code.  The transform
// is stored row-major as U[nao_cart][nao_sph], so
//
//   M_cart = U (N M_sph N) U^T.
//
// Keeping this operation next to Build_Matrix makes the storage convention
// explicit and prevents callers from accidentally reinterpreting U as its
// transpose when nao_cart != nao_sph.
inline std::vector<float> Transform_Weighted_Matrix_To_Cartesian(
    int nao_cart, int nao_sph, const std::vector<float>& transform_cart_sph,
    const std::vector<float>& spherical_norms,
    const std::vector<float>& matrix_sph)
{
    if (nao_cart < 0 || nao_sph < 0)
        throw std::invalid_argument(
            "spherical-to-Cartesian matrix dimensions must be non-negative");

    const std::size_t nc = static_cast<std::size_t>(nao_cart);
    const std::size_t ns = static_cast<std::size_t>(nao_sph);
    if ((ns != 0U && nc > std::numeric_limits<std::size_t>::max() / ns) ||
        (ns != 0U && ns > std::numeric_limits<std::size_t>::max() / ns) ||
        (nc != 0U && nc > std::numeric_limits<std::size_t>::max() / nc))
        throw std::overflow_error(
            "spherical-to-Cartesian matrix size overflow");

    if (transform_cart_sph.size() != nc * ns ||
        spherical_norms.size() != ns || matrix_sph.size() != ns * ns)
        throw std::invalid_argument(
            "spherical-to-Cartesian matrix input dimensions are inconsistent");

    std::vector<double> weighted_sph(ns * ns, 0.0);
    for (std::size_t p = 0; p < ns; ++p)
        for (std::size_t q = 0; q < ns; ++q)
            weighted_sph[p * ns + q] =
                static_cast<double>(spherical_norms[p]) *
                static_cast<double>(matrix_sph[p * ns + q]) *
                static_cast<double>(spherical_norms[q]);

    std::vector<double> left(nc * ns, 0.0);
    for (std::size_t a = 0; a < nc; ++a)
        for (std::size_t q = 0; q < ns; ++q)
            for (std::size_t p = 0; p < ns; ++p)
                left[a * ns + q] +=
                    static_cast<double>(transform_cart_sph[a * ns + p]) *
                    weighted_sph[p * ns + q];

    std::vector<float> matrix_cart(nc * nc, 0.0f);
    for (std::size_t a = 0; a < nc; ++a)
        for (std::size_t b = 0; b < nc; ++b)
        {
            double value = 0.0;
            for (std::size_t q = 0; q < ns; ++q)
                value += left[a * ns + q] *
                         static_cast<double>(
                             transform_cart_sph[b * ns + q]);
            matrix_cart[a * nc + b] = static_cast<float>(value);
        }
    return matrix_cart;
}
}  // namespace qc_cart2sph
