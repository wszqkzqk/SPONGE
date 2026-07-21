#pragma once

#include <limits>
#include <string>

#include "../common.h"
#include "../xponge/ir/forcefield.h"

static inline bool Torsion_Host_Float_Is_Finite(float value)
{
    return Float_Memory_Is_Finite(&value);
}

static inline std::string Validate_Torsions(const Xponge::Torsions& torsions,
                                            std::size_t atom_numbers,
                                            const char* interaction_name,
                                            bool harmonic_improper)
{
    const std::size_t term_numbers = torsions.atom_a.size();
    const std::string prefix = std::string(interaction_name) + " dihedral";
    if (term_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return prefix + " interaction count exceeds the kernel index range";
    }
    if (atom_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return "system atom count exceeds the dihedral kernel index range";
    }
    if (torsions.atom_b.size() != term_numbers ||
        torsions.atom_c.size() != term_numbers ||
        torsions.atom_d.size() != term_numbers ||
        torsions.pk.size() != term_numbers ||
        torsions.pn.size() != term_numbers ||
        torsions.ipn.size() != term_numbers ||
        torsions.gamc.size() != term_numbers ||
        torsions.gams.size() != term_numbers)
    {
        return prefix + " interaction arrays have inconsistent lengths";
    }

    for (std::size_t term = 0; term < term_numbers; term++)
    {
        const int atoms[4] = {torsions.atom_a[term], torsions.atom_b[term],
                              torsions.atom_c[term], torsions.atom_d[term]};
        for (int atom : atoms)
        {
            if (atom < 0 || static_cast<std::size_t>(atom) >= atom_numbers)
            {
                return prefix + " interaction " + std::to_string(term) +
                       " has atom index " + std::to_string(atom) +
                       " outside the system";
            }
        }
        for (int first = 0; first < 4; first++)
        {
            for (int second = first + 1; second < 4; second++)
            {
                if (atoms[first] == atoms[second])
                {
                    return prefix + " interaction " + std::to_string(term) +
                           " repeats atom " + std::to_string(atoms[first]);
                }
            }
        }

        const float parameters[4] = {torsions.pk[term], torsions.pn[term],
                                     torsions.gamc[term], torsions.gams[term]};
        const char* parameter_names[4] = {"pk", "pn", "gamc", "gams"};
        for (int parameter = 0; parameter < 4; parameter++)
        {
            if (!Torsion_Host_Float_Is_Finite(parameters[parameter]))
            {
                return prefix + " interaction " + std::to_string(term) +
                       " has non-finite " + parameter_names[parameter];
            }
        }
        if (harmonic_improper)
        {
            if (torsions.pk[term] < 0.0f)
            {
                return prefix + " interaction " + std::to_string(term) +
                       " has a negative harmonic force constant";
            }
            if (torsions.pn[term] != 0.0f || torsions.ipn[term] != 0 ||
                torsions.gams[term] != 0.0f)
            {
                return prefix + " interaction " + std::to_string(term) +
                       " has nonzero unused periodic parameters";
            }
        }
    }
    return "";
}
