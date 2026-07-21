#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "quantum_chemistry/structure/analytical_norms.hpp"
#include "quantum_chemistry/structure/cart2sph.hpp"

int main()
{
    try
    {
        for (int l = 0; l <= 6; ++l)
        {
            const std::size_t dim_cart =
                qc_cart2sph::Cartesian_Dimension(l);
            const std::size_t dim_sph =
                qc_cart2sph::Spherical_Dimension(l);
            const std::vector<float> block = qc_cart2sph::Build_Block(l);
            std::cout << "BLOCK " << l << ' ' << dim_cart << ' ' << dim_sph;
            std::cout << std::setprecision(9) << std::scientific;
            for (float value : block) std::cout << ' ' << value;
            std::cout << '\n';
        }

        const std::vector<int> l_list = {5, 6};
        const int nao_cart = 21 + 28;
        const int nao_sph = 11 + 13;
        const std::vector<float> transform =
            qc_cart2sph::Build_Matrix(l_list, nao_cart, nao_sph);
        const std::vector<float> block5 = qc_cart2sph::Build_Block(5);
        const std::vector<float> block6 = qc_cart2sph::Build_Block(6);
        for (int row = 0; row < nao_cart; ++row)
        {
            for (int column = 0; column < nao_sph; ++column)
            {
                float expected = 0.0f;
                if (row < 21 && column < 11)
                    expected = block5[static_cast<std::size_t>(row) * 11 +
                                      column];
                else if (row >= 21 && column >= 11)
                    expected = block6[static_cast<std::size_t>(row - 21) * 13 +
                                      column - 11];
                if (transform[static_cast<std::size_t>(row) * nao_sph +
                              column] != expected)
                    throw std::runtime_error(
                        "combined Cartesian-to-spherical block is misplaced");
            }
        }
        const std::vector<float> norms = qc_analytical_norms::Build(
            l_list, {1, 1}, {0, 1}, {0.8f, 1.2f}, {1.0f, 1.0f}, {0, 21},
            {0, 11}, true, transform, nao_cart, nao_sph);
        if (norms.size() != static_cast<std::size_t>(nao_sph) ||
            !std::all_of(norms.begin(), norms.end(), [](float value) {
                return qc_analytical_norms::Float_Is_Positive_Normal(value);
            }))
            throw std::runtime_error(
                "high-angular-momentum analytical norms are invalid");

        // A deliberately rectangular transform catches row/column storage
        // confusion: U is [Cartesian=3][spherical=2], not its transpose.
        const std::vector<float> rectangular_u = {1.0f, 2.0f, 3.0f,
                                                   4.0f, 5.0f, 6.0f};
        const std::vector<float> rectangular_norms = {2.0f, 3.0f};
        const std::vector<float> spherical_matrix = {7.0f, 11.0f, 13.0f,
                                                      17.0f};
        const std::vector<float> cartesian_matrix =
            qc_cart2sph::Transform_Weighted_Matrix_To_Cartesian(
                3, 2, rectangular_u, rectangular_norms, spherical_matrix);
        const std::vector<float> expected_cartesian = {
            928.0f, 2040.0f, 3152.0f, 2016.0f, 4428.0f,
            6840.0f, 3104.0f, 6816.0f, 10528.0f};
        if (cartesian_matrix != expected_cartesian)
            throw std::runtime_error(
                "weighted spherical matrix used the wrong transform layout");

        bool rejected_bad_dimensions = false;
        try
        {
            (void)qc_cart2sph::Build_Matrix(l_list, nao_cart - 1, nao_sph);
        }
        catch (const std::invalid_argument&)
        {
            rejected_bad_dimensions = true;
        }
        bool rejected_negative_l = false;
        try
        {
            (void)qc_cart2sph::Build_Block(-1);
        }
        catch (const std::invalid_argument&)
        {
            rejected_negative_l = true;
        }
        bool rejected_float_underflow = false;
        try
        {
            const long double nonzero_below_float = std::ldexp(
                static_cast<long double>(std::numeric_limits<float>::min()),
                -100);
            (void)qc_cart2sph::Checked_Float_Coefficient(
                nonzero_below_float);
        }
        catch (const std::overflow_error&)
        {
            rejected_float_underflow = true;
        }

        std::vector<float> zero_column_transform = transform;
        for (int row = 21; row < nao_cart; ++row)
            zero_column_transform[static_cast<std::size_t>(row) * nao_sph +
                                  11] = 0.0f;
        bool rejected_zero_column = false;
        try
        {
            (void)qc_analytical_norms::Build(
                l_list, {1, 1}, {0, 1}, {0.8f, 1.2f}, {1.0f, 1.0f},
                {0, 21}, {0, 11}, true, zero_column_transform, nao_cart,
                nao_sph);
        }
        catch (const std::runtime_error&)
        {
            rejected_zero_column = true;
        }

        std::vector<float> subnormal_transform = transform;
        subnormal_transform[0] = std::numeric_limits<float>::denorm_min();
        bool rejected_subnormal_transform = false;
        try
        {
            (void)qc_analytical_norms::Build(
                l_list, {1, 1}, {0, 1}, {0.8f, 1.2f}, {1.0f, 1.0f},
                {0, 21}, {0, 11}, true, subnormal_transform, nao_cart,
                nao_sph);
        }
        catch (const std::runtime_error&)
        {
            rejected_subnormal_transform = true;
        }

        bool rejected_bad_density_dimensions = false;
        try
        {
            (void)qc_cart2sph::Transform_Weighted_Matrix_To_Cartesian(
                3, 2, rectangular_u, {1.0f}, spherical_matrix);
        }
        catch (const std::invalid_argument&)
        {
            rejected_bad_density_dimensions = true;
        }

        if (!rejected_bad_dimensions || !rejected_negative_l ||
            !rejected_float_underflow || !rejected_zero_column ||
            !rejected_subnormal_transform ||
            !rejected_bad_density_dimensions)
            throw std::runtime_error(
                "invalid Cartesian-to-spherical input was not rejected");
        std::cout << "VALIDATION 1\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
