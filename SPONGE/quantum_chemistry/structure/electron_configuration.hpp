#pragma once

#include <limits>
#include <stdexcept>

namespace sponge_qc_electrons
{

struct Electron_Configuration
{
    int total = 0;
    int alpha = 0;
    int beta = 0;
    float occupation_factor = 0.0f;
};

inline long long Electron_Count_From_Charge(int molecular_charge) noexcept
{
    // Cast before negation: -INT_MIN is undefined for int, but every int is
    // exactly representable as long long.
    return -static_cast<long long>(molecular_charge);
}

inline long long Add_Effective_Nuclear_Charge(long long electron_count,
                                              int effective_charge)
{
    if (effective_charge <= 0)
        throw std::domain_error(
            "effective nuclear charge must be a positive integer");
    if (electron_count >
        std::numeric_limits<long long>::max() - effective_charge)
        throw std::overflow_error(
            "molecular electron count overflows long long");
    return electron_count + effective_charge;
}

inline Electron_Configuration Resolve_Electron_Configuration(
    long long electron_count, int multiplicity, bool unrestricted)
{
    // Zero-electron systems are intentionally outside the current SCF
    // contract.  Supporting them would also require auditing the RI,
    // eigensolver, energy-weighted-density, and gradient zero-occupancy paths.
    if (electron_count <= 0)
        throw std::domain_error(
            "the QC region must contain at least one explicit electron");
    if (electron_count > std::numeric_limits<int>::max())
        throw std::overflow_error(
            "molecular electron count exceeds the int storage contract");
    if (multiplicity < 1)
        throw std::domain_error("spin multiplicity must be at least one");

    const long long spin_e = static_cast<long long>(multiplicity) - 1LL;
    if (spin_e > electron_count)
        throw std::domain_error(
            "spin multiplicity requires more unpaired electrons than exist");
    if (((electron_count + spin_e) & 1LL) != 0)
        throw std::domain_error(
            "electron count and spin multiplicity have inconsistent parity");

    const long long alpha = (electron_count + spin_e) / 2LL;
    const long long beta = (electron_count - spin_e) / 2LL;
    if (!unrestricted && spin_e != 0)
        throw std::domain_error(
            "restricted SCF requires a closed-shell singlet");

    Electron_Configuration result;
    result.total = static_cast<int>(electron_count);
    result.alpha = static_cast<int>(alpha);
    result.beta = unrestricted ? static_cast<int>(beta) : 0;
    result.occupation_factor = unrestricted ? 1.0f : 2.0f;
    return result;
}

inline void Validate_AO_Capacity(const Electron_Configuration& configuration,
                                 int nao)
{
    if (nao <= 0)
        throw std::domain_error(
            "orbital basis must provide at least one atomic orbital");
    if (configuration.alpha > nao || configuration.beta > nao)
        throw std::domain_error(
            "occupied spin orbitals exceed the orbital-basis capacity");
}

}  // namespace sponge_qc_electrons
