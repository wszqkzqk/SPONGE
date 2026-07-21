#pragma once

#include "torsion_parse.hpp"

namespace Xponge
{

static void Native_Load_Dihedrals(Torsions* dihedrals, CONTROLLER* controller,
                                  const char* module_name = "dihedral")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    const Native_Torsion_Input_Context context = {
        input_path, "dihedral_in_file", "Xponge::Native_Load_Dihedrals",
        controller};
    std::ifstream input(input_path);
    if (!input.is_open())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Dihedrals",
            "Reason:\n\tfailed to open dihedral_in_file\n\tInput file: %s\n",
            input_path);
    }
    int dihedral_numbers = Native_Torsion_Parse_Int(
        Native_Torsion_Read_Token(&input, context, "interaction count", -1),
        context, "interaction count", -1);
    if (dihedral_numbers < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Dihedrals",
            "Reason:\n\tdihedral_in_file has a negative interaction count\n"
            "\tInput file: %s\n",
            input_path);
    }
    Torsions parsed;
    try
    {
        for (int i = 0; i < dihedral_numbers; i++)
        {
            int integer_values[5];
            const char* integer_fields[5] = {
                "atom A", "atom B", "atom C", "atom D", "multiplicity"};
            for (int field = 0; field < 5; field++)
            {
                integer_values[field] = Native_Torsion_Parse_Int(
                    Native_Torsion_Read_Token(&input, context,
                                              integer_fields[field], i),
                    context, integer_fields[field], i);
            }
            const double pk_double = Native_Torsion_Parse_Double(
                Native_Torsion_Read_Token(&input, context, "force constant",
                                          i),
                context, "force constant", i);
            const double phase_double = Native_Torsion_Parse_Double(
                Native_Torsion_Read_Token(&input, context, "phase", i),
                context, "phase", i);
            const float pk = Native_Torsion_Checked_Float(
                pk_double, context, "force constant", "parameter", i);
            const float phase = Native_Torsion_Checked_Float(
                phase_double, context, "phase", "parameter", i);
            const int multiplicity = integer_values[4];
            const float pn = static_cast<float>(multiplicity);
            if (static_cast<double>(pn) !=
                static_cast<double>(multiplicity))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat,
                    "Xponge::Native_Load_Dihedrals",
                    "Reason:\n\tdihedral_in_file interaction %d "
                    "multiplicity %d cannot be represented exactly by the "
                    "force kernel\n\tInput file: %s\n",
                    i, multiplicity, input_path);
            }
            const double gamc_double = static_cast<double>(cosf(phase)) *
                                       static_cast<double>(pk);
            const double gams_double = static_cast<double>(sinf(phase)) *
                                       static_cast<double>(pk);
            const float gamc = Native_Torsion_Checked_Float(
                gamc_double, context, "cosine", "coefficient", i);
            const float gams = Native_Torsion_Checked_Float(
                gams_double, context, "sine", "coefficient", i);
            parsed.atom_a.push_back(integer_values[0]);
            parsed.atom_b.push_back(integer_values[1]);
            parsed.atom_c.push_back(integer_values[2]);
            parsed.atom_d.push_back(integer_values[3]);
            parsed.pk.push_back(pk);
            parsed.pn.push_back(pn);
            parsed.ipn.push_back(multiplicity);
            parsed.gamc.push_back(gamc);
            parsed.gams.push_back(gams);
        }
    }
    catch (const std::bad_alloc&)
    {
        Native_Torsion_Throw_Storage_Error(context, dihedral_numbers);
        return;
    }
    catch (const std::length_error&)
    {
        Native_Torsion_Throw_Storage_Error(context, dihedral_numbers);
        return;
    }
    Native_Torsion_Ensure_End(&input, context);
    Native_Torsion_Close(&input, context);
    *dihedrals = std::move(parsed);
}

static void Native_Load_Dihedrals(System* system, CONTROLLER* controller)
{
    Native_Load_Dihedrals(&system->classical_force_field.dihedrals, controller);
}

}  // namespace Xponge
