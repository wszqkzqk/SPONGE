#pragma once

#include <limits>
#include <string>

#include "../common.h"
#include "../xponge/ir/forcefield.h"

static inline std::string Validate_Angles(const Xponge::Angles& angles,
                                          std::size_t atom_numbers,
                                          const char* interaction_name)
{
    const std::size_t term_numbers = angles.atom_a.size();
    const std::string prefix = interaction_name;
    if (term_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return prefix + " interaction count exceeds the kernel index range";
    }
    if (atom_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return "system atom count exceeds the angle kernel index range";
    }
    if (angles.atom_b.size() != term_numbers ||
        angles.atom_c.size() != term_numbers ||
        angles.k.size() != term_numbers || angles.theta0.size() != term_numbers)
    {
        return prefix + " interaction arrays have inconsistent lengths";
    }

    for (std::size_t term = 0; term < term_numbers; term++)
    {
        const int atoms[3] = {angles.atom_a[term], angles.atom_b[term],
                              angles.atom_c[term]};
        for (int atom : atoms)
        {
            if (atom < 0 || static_cast<std::size_t>(atom) >= atom_numbers)
            {
                return prefix + " interaction " + std::to_string(term) +
                       " has atom index " + std::to_string(atom) +
                       " outside [0, " + std::to_string(atom_numbers) + ")";
            }
        }
        if (atoms[0] == atoms[1] || atoms[0] == atoms[2] ||
            atoms[1] == atoms[2])
        {
            return prefix + " interaction " + std::to_string(term) +
                   " repeats atom index in (" + std::to_string(atoms[0]) +
                   ", " + std::to_string(atoms[1]) + ", " +
                   std::to_string(atoms[2]) +
                   "); a harmonic angle requires three distinct atoms";
        }
        if (!Float_Memory_Is_Finite(&angles.k[term]))
        {
            return prefix + " interaction " + std::to_string(term) +
                   " has a non-finite force constant";
        }
        if (!Float_Memory_Is_Finite(&angles.theta0[term]))
        {
            return prefix + " interaction " + std::to_string(term) +
                   " has a non-finite equilibrium angle";
        }
        if (angles.theta0[term] < 0.0f || angles.theta0[term] > CONSTANT_Pi)
        {
            return prefix + " interaction " + std::to_string(term) +
                   " has equilibrium angle " +
                   std::to_string(angles.theta0[term]) + " outside [0, pi]";
        }
    }
    return "";
}

static inline std::string Validate_Urey_Bradley(
    const Xponge::UreyBradley& urey_bradley, std::size_t atom_numbers)
{
    const std::size_t term_numbers = urey_bradley.atom_a.size();
    if (term_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return "Urey-Bradley interaction count exceeds the kernel index "
               "range";
    }
    if (atom_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return "system atom count exceeds the Urey-Bradley kernel index "
               "range";
    }
    if (urey_bradley.atom_b.size() != term_numbers ||
        urey_bradley.atom_c.size() != term_numbers ||
        urey_bradley.angle_k.size() != term_numbers ||
        urey_bradley.angle_theta0.size() != term_numbers ||
        urey_bradley.bond_k.size() != term_numbers ||
        urey_bradley.bond_r0.size() != term_numbers)
    {
        return "Urey-Bradley interaction arrays have inconsistent lengths";
    }

    for (std::size_t term = 0; term < term_numbers; term++)
    {
        const int atoms[3] = {urey_bradley.atom_a[term],
                              urey_bradley.atom_b[term],
                              urey_bradley.atom_c[term]};
        for (int atom : atoms)
        {
            if (atom < 0 || static_cast<std::size_t>(atom) >= atom_numbers)
            {
                return "Urey-Bradley interaction " + std::to_string(term) +
                       " has atom index " + std::to_string(atom) +
                       " outside [0, " + std::to_string(atom_numbers) + ")";
            }
        }
        if (atoms[0] == atoms[1] || atoms[0] == atoms[2] ||
            atoms[1] == atoms[2])
        {
            return "Urey-Bradley interaction " + std::to_string(term) +
                   " repeats atom index in (" + std::to_string(atoms[0]) +
                   ", " + std::to_string(atoms[1]) + ", " +
                   std::to_string(atoms[2]) +
                   "); a harmonic angle requires three distinct atoms";
        }

        const float parameters[4] = {
            urey_bradley.angle_k[term], urey_bradley.angle_theta0[term],
            urey_bradley.bond_k[term], urey_bradley.bond_r0[term]};
        const char* parameter_names[4] = {
            "angle force constant", "equilibrium angle", "bond force constant",
            "equilibrium distance"};
        for (int parameter = 0; parameter < 4; parameter++)
        {
            if (!Float_Memory_Is_Finite(&parameters[parameter]))
            {
                return "Urey-Bradley interaction " + std::to_string(term) +
                       " has a non-finite " + parameter_names[parameter];
            }
        }
        if (urey_bradley.angle_theta0[term] < 0.0f ||
            urey_bradley.angle_theta0[term] > CONSTANT_Pi)
        {
            return "Urey-Bradley interaction " + std::to_string(term) +
                   " has equilibrium angle " +
                   std::to_string(urey_bradley.angle_theta0[term]) +
                   " outside [0, pi]";
        }
        if (urey_bradley.bond_r0[term] < 0.0f)
        {
            return "Urey-Bradley interaction " + std::to_string(term) +
                   " has a negative equilibrium distance";
        }
    }
    return "";
}
