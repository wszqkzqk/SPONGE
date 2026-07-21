#pragma once

#include "torsion_parse.hpp"

namespace Xponge
{

static void Native_Load_Impropers(Torsions* impropers, CONTROLLER* controller,
                                  const char* module_name = "improper_dihedral")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    const Native_Torsion_Input_Context context = {
        input_path, "improper_dihedral_in_file",
        "Xponge::Native_Load_Impropers", controller};
    std::ifstream input(input_path);
    if (!input.is_open())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Impropers",
            "Reason:\n\tfailed to open improper_dihedral_in_file\n"
            "\tInput file: %s\n",
            input_path);
    }
    int improper_numbers = Native_Torsion_Parse_Int(
        Native_Torsion_Read_Token(&input, context, "interaction count", -1),
        context, "interaction count", -1);
    if (improper_numbers < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Impropers",
            "Reason:\n\timproper_dihedral_in_file has a negative "
            "interaction count\n\tInput file: %s\n",
            input_path);
    }
    Torsions parsed;
    try
    {
        for (int i = 0; i < improper_numbers; i++)
        {
            int atoms[4];
            const char* integer_fields[4] = {"atom A", "atom B", "atom C",
                                             "atom D"};
            for (int field = 0; field < 4; field++)
            {
                atoms[field] = Native_Torsion_Parse_Int(
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
            parsed.atom_a.push_back(atoms[0]);
            parsed.atom_b.push_back(atoms[1]);
            parsed.atom_c.push_back(atoms[2]);
            parsed.atom_d.push_back(atoms[3]);
            parsed.pk.push_back(pk);
            parsed.pn.push_back(0.0f);
            parsed.ipn.push_back(0);
            parsed.gamc.push_back(phase);
            parsed.gams.push_back(0.0f);
        }
    }
    catch (const std::bad_alloc&)
    {
        Native_Torsion_Throw_Storage_Error(context, improper_numbers);
        return;
    }
    catch (const std::length_error&)
    {
        Native_Torsion_Throw_Storage_Error(context, improper_numbers);
        return;
    }
    Native_Torsion_Ensure_End(&input, context);
    Native_Torsion_Close(&input, context);
    *impropers = std::move(parsed);
}

static void Native_Load_Impropers(System* system, CONTROLLER* controller)
{
    Native_Load_Impropers(&system->classical_force_field.impropers, controller);
}

}  // namespace Xponge
