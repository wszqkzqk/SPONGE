#pragma once

#include <array>
#include <stdexcept>
#include <string>

namespace sponge_qc_elements
{

constexpr int MAX_ATOMIC_NUMBER = 86;
constexpr double ANGSTROM_PER_BOHR = 0.52917721092;

inline const std::array<const char*, MAX_ATOMIC_NUMBER + 1>& Symbols()
{
    static const std::array<const char*, MAX_ATOMIC_NUMBER + 1> symbols = {
        "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
        "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",
        "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
        "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
        "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
        "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
        "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn"};
    return symbols;
}

inline const char* Symbol_From_Atomic_Number(int atomic_number)
{
    if (atomic_number < 1 || atomic_number > MAX_ATOMIC_NUMBER) return nullptr;
    return Symbols()[static_cast<std::size_t>(atomic_number)];
}

inline int Atomic_Number_From_Symbol(const std::string& symbol)
{
    const auto& symbols = Symbols();
    for (int atomic_number = 1; atomic_number <= MAX_ATOMIC_NUMBER;
         ++atomic_number)
    {
        if (symbol == symbols[static_cast<std::size_t>(atomic_number)])
            return atomic_number;
    }
    return 0;
}

inline double Cordero_Covalent_Radius_Angstrom(int atomic_number)
{
    // Single-bond covalent radii from Cordero et al., Dalton Trans. 2008,
    // 2832-2838.  The table covers the complete H-Rn identity range used by
    // the element lookup and universal auxiliary metadata; individual
    // orbital bases may support only a subset.  There is intentionally no
    // guessed fallback because it would silently change the atom-centred DFT
    // quadrature.
    static const std::array<double, MAX_ATOMIC_NUMBER + 1> radii = {
        0.00, 0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57,
        0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03,
        1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32,
        1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90,
        1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42,
        1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03,
        2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90,
        1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36,
        1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50};
    if (atomic_number < 1 || atomic_number > MAX_ATOMIC_NUMBER)
    {
        throw std::domain_error(
            "DFT radial scale is unavailable for atomic number " +
            std::to_string(atomic_number));
    }
    return radii[static_cast<std::size_t>(atomic_number)];
}

inline double Cordero_Covalent_Radius_Bohr(int atomic_number)
{
    return Cordero_Covalent_Radius_Angstrom(atomic_number) /
           ANGSTROM_PER_BOHR;
}

}  // namespace sponge_qc_elements
