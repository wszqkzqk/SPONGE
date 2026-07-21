#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>

#include "quantum_chemistry/guess/minao.h"
#include "quantum_chemistry/guess/minao_occupancy.hpp"
#include "quantum_chemistry/structure/ecp.h"
#include "quantum_chemistry/structure/electron_configuration.hpp"
#include "quantum_chemistry/structure/input_contract.hpp"
#include "quantum_chemistry/structure/molecule.h"

void deviceMemcpy(void* to, const void* from, size_t size,
                  deviceMemcpyKind /*kind*/)
{
    std::memcpy(to, from, size);
}

void deviceMemset(void* to, int value, size_t size)
{
    std::memset(to, value, size);
}

void deviceFree(void* pointer) { std::free(pointer); }

namespace
{

int Electron_Count(const std::array<int, 4>& electrons_by_l)
{
    int count = 0;
    for (int occupancy : electrons_by_l) count += occupancy;
    return count;
}

double Sum(const std::vector<double>& values)
{
    double result = 0.0;
    for (double value : values) result += value;
    return result;
}

double Double_From_Bits(std::uint64_t bits)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

bool Global_Occupancies_Are_Valid(
    const sponge_qc_minao::Global_AO_Occupancies& occupations,
    double target_alpha, double target_beta, bool unrestricted)
{
    if (std::fabs(Sum(occupations.alpha) - target_alpha) > 1.0e-12)
        return false;
    if (unrestricted)
    {
        if (occupations.alpha.size() != occupations.beta.size() ||
            std::fabs(Sum(occupations.beta) - target_beta) > 1.0e-12)
            return false;
    }
    else if (!occupations.beta.empty())
    {
        return false;
    }
    for (std::size_t i = 0; i < occupations.alpha.size(); ++i)
    {
        const double alpha = occupations.alpha[i];
        const double beta = unrestricted ? occupations.beta[i] : 0.0;
        if (alpha < 0.0 || beta < 0.0 || alpha > (unrestricted ? 1.0 : 2.0) ||
            beta > 1.0 || alpha + beta > 2.0)
            return false;
    }
    return true;
}

template <typename Function>
bool Rejects(Function&& function)
{
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        return true;
    }
    return false;
}

}  // namespace

int main()
{
    if (sponge_qc_input::Parse_Exact_Int("-2147483648") !=
            std::numeric_limits<int>::min() ||
        sponge_qc_input::Parse_Exact_Int("2147483647") !=
            std::numeric_limits<int>::max() ||
        sponge_qc_input::Parse_Finite_Nonnegative_Float("1e-12") <= 0.0f ||
        !Rejects([] { sponge_qc_input::Parse_Exact_Int("2147483648"); }) ||
        !Rejects([] { sponge_qc_input::Parse_Exact_Int("1x"); }) ||
        !Rejects(
            [] { sponge_qc_input::Parse_Finite_Nonnegative_Float("1e999"); }) ||
        !Rejects(
            []
            { sponge_qc_input::Parse_Finite_Nonnegative_Float("1e-999"); }) ||
        !Rejects([]
                 { sponge_qc_input::Parse_Finite_Nonnegative_Float("nan"); }) ||
        !Rejects([] { sponge_qc_input::Parse_Finite_Nonnegative_Float("-1"); }))
        return EXIT_FAILURE;

    const auto restricted_configuration =
        sponge_qc_electrons::Resolve_Electron_Configuration(2, 1, false);
    const auto unrestricted_configuration =
        sponge_qc_electrons::Resolve_Electron_Configuration(9, 2, true);
    if (restricted_configuration.total != 2 ||
        restricted_configuration.alpha != 1 ||
        restricted_configuration.beta != 0 ||
        restricted_configuration.occupation_factor != 2.0f ||
        unrestricted_configuration.total != 9 ||
        unrestricted_configuration.alpha != 5 ||
        unrestricted_configuration.beta != 4 ||
        unrestricted_configuration.occupation_factor != 1.0f ||
        sponge_qc_electrons::Electron_Count_From_Charge(
            std::numeric_limits<int>::min()) != 2147483648LL ||
        !Rejects(
            []
            {
                sponge_qc_electrons::Resolve_Electron_Configuration(
                    sponge_qc_electrons::Electron_Count_From_Charge(
                        std::numeric_limits<int>::min()),
                    1, false);
            }) ||
        !Rejects(
            []
            {
                const long long electrons =
                    sponge_qc_electrons::Add_Effective_Nuclear_Charge(-2, 1);
                sponge_qc_electrons::Resolve_Electron_Configuration(electrons,
                                                                    1, false);
            }) ||
        !Rejects(
            []
            {
                sponge_qc_electrons::Resolve_Electron_Configuration(0, 1,
                                                                    false);
            }) ||
        !Rejects(
            []
            {
                sponge_qc_electrons::Resolve_Electron_Configuration(2, 4, true);
            }) ||
        !Rejects(
            []
            {
                sponge_qc_electrons::Resolve_Electron_Configuration(2, 2, true);
            }) ||
        !Rejects(
            []
            {
                const auto configuration =
                    sponge_qc_electrons::Resolve_Electron_Configuration(4, 1,
                                                                        false);
                sponge_qc_electrons::Validate_AO_Capacity(configuration, 1);
            }))
        return EXIT_FAILURE;

    if (!QC_ECP_L_Max_Is_Supported(0) ||
        !QC_ECP_L_Max_Is_Supported(QC_ECP_MAX_SEMILOCAL_L + 1) ||
        QC_ECP_L_Max_Is_Supported(-1) ||
        QC_ECP_L_Max_Is_Supported(QC_ECP_MAX_SEMILOCAL_L + 2))
        return EXIT_FAILURE;

    for (int atomic_number = 1;
         atomic_number <= sponge_qc_elements::MAX_ATOMIC_NUMBER;
         ++atomic_number)
    {
        const char* symbol =
            sponge_qc_elements::Symbol_From_Atomic_Number(atomic_number);
        if (symbol == nullptr || sponge_qc_elements::Atomic_Number_From_Symbol(
                                     symbol) != atomic_number)
            return EXIT_FAILURE;
    }
    if (sponge_qc_elements::Atomic_Number_From_Symbol("NotAnElement") != 0 ||
        sponge_qc_elements::Symbol_From_Atomic_Number(0) != nullptr ||
        sponge_qc_elements::Symbol_From_Atomic_Number(87) != nullptr)
        return EXIT_FAILURE;

    QC_MOLECULE molecule{};
    molecule.natm = 1;
    molecule.h_atomic_numbers = {19};
    molecule.h_Z = {9};
    molecule.h_ecp_n_core = {10};
    if (molecule.Atomic_Number(0) != 19 ||
        molecule.Effective_Nuclear_Charge(0) != 9)
        return EXIT_FAILURE;

    const double fluorine_radius =
        sponge_qc_elements::Cordero_Covalent_Radius_Angstrom(9);
    const double potassium_radius =
        sponge_qc_elements::Cordero_Covalent_Radius_Angstrom(
            molecule.Atomic_Number(0));
    const double rubidium_radius =
        sponge_qc_elements::Cordero_Covalent_Radius_Angstrom(37);
    const double radon_radius =
        sponge_qc_elements::Cordero_Covalent_Radius_Angstrom(86);
    if (fluorine_radius != 0.57 || potassium_radius != 2.03 ||
        rubidium_radius != 2.20 || radon_radius != 1.50 ||
        potassium_radius == fluorine_radius)
        return EXIT_FAILURE;

    bool rejected_low = false;
    bool rejected_high = false;
    try
    {
        (void)sponge_qc_elements::Cordero_Covalent_Radius_Angstrom(0);
    }
    catch (const std::domain_error&)
    {
        rejected_low = true;
    }
    try
    {
        (void)sponge_qc_elements::Cordero_Covalent_Radius_Angstrom(87);
    }
    catch (const std::domain_error&)
    {
        rejected_high = true;
    }
    if (!rejected_low || !rejected_high) return EXIT_FAILURE;

    const std::array<int, 4> potassium_valence =
        sponge_qc_minao::Explicit_Electrons_By_Angular_Momentum(19, 10);
    const std::array<int, 4> iron_valence =
        sponge_qc_minao::Explicit_Electrons_By_Angular_Momentum(26, 10);
    if (potassium_valence != std::array<int, 4>{3, 6, 0, 0} ||
        Electron_Count(potassium_valence) != 9 ||
        iron_valence != std::array<int, 4>{4, 6, 6, 0} ||
        Electron_Count(iron_valence) != 16)
        return EXIT_FAILURE;

    // A charged restricted singlet must publish the total density target
    // (occupation_factor*n_alpha), not mistake n_beta=0 for zero electrons.
    const sponge_qc_minao::Global_AO_Occupancies charged_singlet =
        sponge_qc_minao::Allocate_Global_AO_Occupancies({2.0, 1.0, 1.0, 0.0},
                                                        false, 1, 0, 2.0);
    if (!Global_Occupancies_Are_Valid(charged_singlet, 2.0, 0.0, false) ||
        charged_singlet.alpha[0] == 1.0)
    {
        std::fprintf(stderr,
                     "charged restricted MINAO target was not projected\n");
        return EXIT_FAILURE;
    }

    // K with a 10-electron ECP contributes nine explicit electrons.  Adding
    // H gives a neutral preference of ten; KH+ is a nine-electron doublet and
    // therefore requires exactly five alpha and four beta electrons without
    // any spin occupation exceeding one.
    const std::vector<double> kh_ecp_preference = {1.5, 1.5, 2.0,
                                                   2.0, 2.0, 1.0};
    if (Electron_Count(potassium_valence) + 1 != 10) return EXIT_FAILURE;
    const sponge_qc_minao::Global_AO_Occupancies charged_doublet =
        sponge_qc_minao::Allocate_Global_AO_Occupancies(kh_ecp_preference, true,
                                                        5, 4, 1.0);
    if (!Global_Occupancies_Are_Valid(charged_doublet, 5.0, 4.0, true))
    {
        std::fprintf(stderr,
                     "charged ECP UKS MINAO spin targets are incorrect\n");
        return EXIT_FAILURE;
    }

    const std::vector<double> projection_stress = {0.0, 0.2, 0.7,
                                                   1.1, 1.8, 2.0};
    for (int alpha = 0; alpha <= 6; ++alpha)
    {
        const auto restricted_projection =
            sponge_qc_minao::Allocate_Global_AO_Occupancies(
                projection_stress, false, alpha, 0, 2.0);
        if (!Global_Occupancies_Are_Valid(restricted_projection, 2.0 * alpha,
                                          0.0, false))
            return EXIT_FAILURE;
        for (int beta = 0; beta <= 6; ++beta)
        {
            const auto unrestricted_projection =
                sponge_qc_minao::Allocate_Global_AO_Occupancies(
                    projection_stress, true, alpha, beta, 1.0);
            if (!Global_Occupancies_Are_Valid(unrestricted_projection, alpha,
                                              beta, true))
                return EXIT_FAILURE;
        }
    }

    bool rejected_nan = false;
    bool rejected_infinity = false;
    bool rejected_restricted_factor = false;
    bool rejected_unrestricted_factor = false;
    try
    {
        (void)sponge_qc_minao::Allocate_Global_AO_Occupancies(
            {Double_From_Bits(UINT64_C(0x7ff8000000000000))}, false, 0, 0, 2.0);
    }
    catch (const std::domain_error&)
    {
        rejected_nan = true;
    }
    try
    {
        (void)sponge_qc_minao::Allocate_Global_AO_Occupancies(
            {0.0}, false, 0, 0, Double_From_Bits(UINT64_C(0x7ff0000000000000)));
    }
    catch (const std::domain_error&)
    {
        rejected_infinity = true;
    }
    try
    {
        (void)sponge_qc_minao::Allocate_Global_AO_Occupancies({1.0}, false, 0,
                                                              0, 1.0);
    }
    catch (const std::domain_error&)
    {
        rejected_restricted_factor = true;
    }
    try
    {
        (void)sponge_qc_minao::Allocate_Global_AO_Occupancies({1.0}, true, 0, 0,
                                                              2.0);
    }
    catch (const std::domain_error&)
    {
        rejected_unrestricted_factor = true;
    }
    if (!rejected_nan || !rejected_infinity || !rejected_restricted_factor ||
        !rejected_unrestricted_factor)
    {
        std::fprintf(stderr,
                     "MINAO accepted a non-finite or inconsistent runtime "
                     "contract\n");
        return EXIT_FAILURE;
    }

    // Exercise the production MINAO publication path, including its full
    // density memset.  Two s shells on H share the one alpha electron.
    QC_MOLECULE minao_molecule{};
    minao_molecule.natm = 1;
    minao_molecule.nelectron = 1;
    minao_molecule.nbas = 2;
    minao_molecule.nao_cart = 2;
    minao_molecule.nao_sph = 2;
    minao_molecule.nao = 2;
    minao_molecule.nao2 = 4;
    minao_molecule.is_spherical = 0;
    minao_molecule.h_atomic_numbers = {1};
    minao_molecule.h_Z = {1};
    minao_molecule.h_ecp_n_core = {0};
    minao_molecule.h_bas.assign(16, 0);
    minao_molecule.h_l_list = {0, 0};
    minao_molecule.h_ao_offsets = {0, 1};
    minao_molecule.h_ao_offsets_sph = {0, 1};

    QC_SCF_Runtime_State minao_runtime{};
    minao_runtime.unrestricted = true;
    minao_runtime.n_alpha = 1;
    minao_runtime.n_beta = 0;
    minao_runtime.occ_factor = 1.0f;
    float alpha_density[4] = {91.0f, 92.0f, 93.0f, 94.0f};
    float beta_density[4] = {-91.0f, -92.0f, -93.0f, -94.0f};
    QC_Build_Minao_Guess(minao_molecule, minao_runtime, alpha_density,
                         beta_density);
    if (std::fabs(alpha_density[0] - 0.5f) > 1.0e-6f ||
        alpha_density[1] != 0.0f || alpha_density[2] != 0.0f ||
        std::fabs(alpha_density[3] - 0.5f) > 1.0e-6f)
        return EXIT_FAILURE;
    for (float value : beta_density)
        if (value != 0.0f) return EXIT_FAILURE;

    minao_molecule.h_l_list[0] = -1;
    if (!Rejects(
            [&]
            {
                QC_Build_Minao_Guess(minao_molecule, minao_runtime,
                                     alpha_density, beta_density);
            }))
        return EXIT_FAILURE;

    std::printf("K %d %d %.2f %.2f\n", molecule.Atomic_Number(0),
                molecule.Effective_Nuclear_Charge(0), potassium_radius,
                fluorine_radius);
    return EXIT_SUCCESS;
}
