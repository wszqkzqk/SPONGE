#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace sponge_qc_minao
{

inline bool Double_Is_Finite(double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires a 64-bit double representation");
    std::memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return (bits & UINT64_C(0x7ff0000000000000)) !=
           UINT64_C(0x7ff0000000000000);
}

constexpr int MAX_SUBSHELLS = 8;
constexpr int MAX_SUPPORTED_ATOMIC_NUMBER = 36;

struct Atom_Config
{
    int n_subshells;
    std::array<int, MAX_SUBSHELLS> l;
    std::array<int, MAX_SUBSHELLS> occupancy;
};

// Ground-state configurations ordered from the innermost subshell outward.
// This order lets an ECP remove exactly its n_core innermost electrons rather
// than reinterpreting the remaining charge as a different element.
inline const std::array<Atom_Config, MAX_SUPPORTED_ATOMIC_NUMBER + 1>&
Atom_Configurations()
{
    static const std::array<Atom_Config, MAX_SUPPORTED_ATOMIC_NUMBER + 1>
        configurations = {{
            {0, {}, {}},
            {1, {0}, {1}},
            {1, {0}, {2}},
            {2, {0, 0}, {2, 1}},
            {2, {0, 0}, {2, 2}},
            {3, {0, 0, 1}, {2, 2, 1}},
            {3, {0, 0, 1}, {2, 2, 2}},
            {3, {0, 0, 1}, {2, 2, 3}},
            {3, {0, 0, 1}, {2, 2, 4}},
            {3, {0, 0, 1}, {2, 2, 5}},
            {3, {0, 0, 1}, {2, 2, 6}},
            {4, {0, 0, 1, 0}, {2, 2, 6, 1}},
            {4, {0, 0, 1, 0}, {2, 2, 6, 2}},
            {5, {0, 0, 1, 0, 1}, {2, 2, 6, 2, 1}},
            {5, {0, 0, 1, 0, 1}, {2, 2, 6, 2, 2}},
            {5, {0, 0, 1, 0, 1}, {2, 2, 6, 2, 3}},
            {5, {0, 0, 1, 0, 1}, {2, 2, 6, 2, 4}},
            {5, {0, 0, 1, 0, 1}, {2, 2, 6, 2, 5}},
            {5, {0, 0, 1, 0, 1}, {2, 2, 6, 2, 6}},
            {6, {0, 0, 1, 0, 1, 0}, {2, 2, 6, 2, 6, 1}},
            {6, {0, 0, 1, 0, 1, 0}, {2, 2, 6, 2, 6, 2}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 1}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 2}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 3}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 1, 5}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 5}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 6}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 7}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 8}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 1, 10}},
            {7, {0, 0, 1, 0, 1, 0, 2}, {2, 2, 6, 2, 6, 2, 10}},
            {8, {0, 0, 1, 0, 1, 0, 2, 1}, {2, 2, 6, 2, 6, 2, 10, 1}},
            {8, {0, 0, 1, 0, 1, 0, 2, 1}, {2, 2, 6, 2, 6, 2, 10, 2}},
            {8, {0, 0, 1, 0, 1, 0, 2, 1}, {2, 2, 6, 2, 6, 2, 10, 3}},
            {8, {0, 0, 1, 0, 1, 0, 2, 1}, {2, 2, 6, 2, 6, 2, 10, 4}},
            {8, {0, 0, 1, 0, 1, 0, 2, 1}, {2, 2, 6, 2, 6, 2, 10, 5}},
            {8, {0, 0, 1, 0, 1, 0, 2, 1}, {2, 2, 6, 2, 6, 2, 10, 6}},
        }};
    return configurations;
}

inline std::array<int, 4> Explicit_Electrons_By_Angular_Momentum(
    int atomic_number, int core_electron_count)
{
    if (atomic_number < 1 || atomic_number > MAX_SUPPORTED_ATOMIC_NUMBER)
        throw std::domain_error(
            "MINAO configuration is unavailable for atomic number " +
            std::to_string(atomic_number));
    if (core_electron_count < 0 || core_electron_count >= atomic_number)
        throw std::domain_error("invalid MINAO ECP core-electron count");

    const Atom_Config& config =
        Atom_Configurations()[static_cast<std::size_t>(atomic_number)];
    std::array<int, 4> electrons_by_l{};
    int remaining_core = core_electron_count;
    int explicit_electron_count = 0;
    for (int subshell = 0; subshell < config.n_subshells; ++subshell)
    {
        const int occupancy =
            config.occupancy[static_cast<std::size_t>(subshell)];
        const int removed = std::min(remaining_core, occupancy);
        remaining_core -= removed;
        const int explicit_occupancy = occupancy - removed;
        explicit_electron_count += explicit_occupancy;
        const int l = config.l[static_cast<std::size_t>(subshell)];
        electrons_by_l[static_cast<std::size_t>(l)] += explicit_occupancy;
    }
    if (remaining_core != 0 ||
        explicit_electron_count != atomic_number - core_electron_count)
        throw std::logic_error(
            "MINAO core removal does not conserve explicit electrons");
    return electrons_by_l;
}

struct Global_AO_Occupancies
{
    // Restricted calculations store the total (doubly occupiable) density in
    // alpha and leave beta empty.  Unrestricted calculations store one
    // occupation-bounded density per spin channel.
    std::vector<double> alpha;
    std::vector<double> beta;
};

namespace detail
{

struct Projection_Event
{
    double shift;
    double preferred;
    bool enters_free_set;
};

// Euclidean projection onto {x : sum(x)=target, 0<=x_i<=upper_bound}.
// The Lagrange-multiplier solution is x_i=clip(preferred_i+shift,0,upper).
// Sweep its exact breakpoints instead of iteratively rescaling a density;
// this preserves the neutral-atom preference while enforcing the molecular
// electron constraint as part of the construction itself.
inline std::vector<double> Project_Occupancies(
    const std::vector<double>& preferred, double target, double upper_bound)
{
    if (!(upper_bound > 0.0) || !Double_Is_Finite(upper_bound) ||
        !Double_Is_Finite(target))
        throw std::domain_error("invalid MINAO occupation projection target");
    const std::size_t count = preferred.size();
    const double capacity = upper_bound * static_cast<double>(count);
    if (target < 0.0 || target > capacity)
        throw std::domain_error(
            "MINAO target electron count exceeds the AO occupation capacity");
    for (double value : preferred)
        if (!Double_Is_Finite(value) || value < 0.0)
            throw std::domain_error(
                "MINAO preferred AO occupation is not finite and non-negative");

    std::vector<double> projected(count, 0.0);
    if (target == 0.0 || count == 0) return projected;
    if (target == capacity)
    {
        std::fill(projected.begin(), projected.end(), upper_bound);
        return projected;
    }

    std::vector<Projection_Event> events;
    events.reserve(2 * count);
    for (double value : preferred)
    {
        events.push_back({-value, value, true});
        events.push_back({upper_bound - value, value, false});
    }
    std::sort(events.begin(), events.end(),
              [](const Projection_Event& lhs, const Projection_Event& rhs)
              {
                  if (lhs.shift != rhs.shift) return lhs.shift < rhs.shift;
                  return lhs.enters_free_set && !rhs.enters_free_set;
              });

    std::size_t event = 0;
    int free_count = 0;
    int capped_count = 0;
    double free_preferred_sum = 0.0;
    double solution_shift = 0.0;
    bool found = false;
    while (event < events.size())
    {
        const double current_shift = events[event].shift;
        while (event < events.size() && events[event].shift == current_shift)
        {
            if (events[event].enters_free_set)
            {
                ++free_count;
                free_preferred_sum += events[event].preferred;
            }
            else
            {
                --free_count;
                ++capped_count;
                free_preferred_sum -= events[event].preferred;
            }
            ++event;
        }

        const double fixed_sum = upper_bound * capped_count;
        const double current_sum =
            fixed_sum + free_preferred_sum + free_count * current_shift;
        if (target == current_sum)
        {
            solution_shift = current_shift;
            found = true;
            break;
        }
        if (event == events.size()) break;

        const double next_shift = events[event].shift;
        const double next_sum =
            fixed_sum + free_preferred_sum + free_count * next_shift;
        if (target >= current_sum && target <= next_sum)
        {
            if (free_count <= 0)
            {
                if (target != current_sum)
                    throw std::logic_error(
                        "MINAO occupation projection crossed an empty free "
                        "set");
                solution_shift = current_shift;
            }
            else
            {
                solution_shift =
                    (target - fixed_sum - free_preferred_sum) / free_count;
            }
            found = true;
            break;
        }
    }
    if (!found)
        throw std::logic_error(
            "MINAO occupation projection could not locate its multiplier");

    double projected_sum = 0.0;
    for (std::size_t i = 0; i < count; ++i)
    {
        projected[i] =
            std::max(0.0, std::min(upper_bound, preferred[i] + solution_shift));
        projected_sum += projected[i];
    }

    // Close only the last few ulps introduced while evaluating the analytic
    // projection.  A material mismatch is an algorithm error, not something
    // to hide by scaling all occupations after the fact.
    double residual = target - projected_sum;
    const double rounding_bound =
        128.0 * std::numeric_limits<double>::epsilon() *
        std::max({1.0, capacity, static_cast<double>(count)});
    if (std::fabs(residual) > rounding_bound)
        throw std::logic_error(
            "MINAO occupation projection has a non-roundoff residual");
    for (std::size_t i = 0; i < count && residual != 0.0; ++i)
    {
        const double lower_delta = -projected[i];
        const double upper_delta = upper_bound - projected[i];
        const double delta =
            std::max(lower_delta, std::min(upper_delta, residual));
        projected[i] += delta;
        residual -= delta;
    }
    if (std::fabs(residual) > rounding_bound)
        throw std::logic_error(
            "MINAO occupation projection could not close roundoff");
    return projected;
}

}  // namespace detail

inline Global_AO_Occupancies Allocate_Global_AO_Occupancies(
    const std::vector<double>& preferred_total, bool unrestricted, int n_alpha,
    int n_beta, double occupation_factor)
{
    if (n_alpha < 0 || n_beta < 0 || !Double_Is_Finite(occupation_factor) ||
        !(occupation_factor > 0.0))
        throw std::domain_error("invalid MINAO molecular occupation target");

    Global_AO_Occupancies result;
    if (!unrestricted)
    {
        if (n_beta != 0 || occupation_factor != 2.0)
            throw std::domain_error(
                "restricted MINAO requires occupation factor two and no beta "
                "orbital target");
        result.alpha = detail::Project_Occupancies(
            preferred_total, occupation_factor * n_alpha, 2.0);
        return result;
    }

    if (occupation_factor != 1.0)
        throw std::domain_error(
            "unrestricted MINAO requires occupation factor one");

    std::vector<double> preferred_spin(preferred_total.size(), 0.0);
    for (std::size_t i = 0; i < preferred_total.size(); ++i)
        preferred_spin[i] = 0.5 * preferred_total[i];
    result.alpha = detail::Project_Occupancies(preferred_spin, n_alpha, 1.0);
    result.beta = detail::Project_Occupancies(preferred_spin, n_beta, 1.0);
    return result;
}

}  // namespace sponge_qc_minao
