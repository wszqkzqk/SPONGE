#pragma once

#include "torsion_parse.hpp"

namespace Xponge
{

static void Native_Load_Bonds(Bonds* bonds, CONTROLLER* controller,
                              const char* module_name = "bond")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    const Native_Torsion_Input_Context context = {
        input_path, "bond_in_file", "Xponge::Native_Load_Bonds", controller};
    std::ifstream input(input_path);
    if (!input.is_open())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Bonds",
            "Reason:\n\tfailed to open bond_in_file\n\tInput file: %s\n",
            input_path);
    }
    const int bond_numbers = Native_Torsion_Parse_Int(
        Native_Torsion_Read_Token(&input, context, "interaction count", -1),
        context, "interaction count", -1);
    if (bond_numbers < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Bonds",
            "Reason:\n\tbond_in_file has a negative interaction count\n"
            "\tInput file: %s\n",
            input_path);
    }
    // Build a temporary value one complete record at a time.  Resizing every
    // destination vector from an untrusted count can attempt many gigabytes
    // before discovering that the file is truncated, and it leaves a
    // partially replaced force field on any parse failure.
    Bonds parsed;
    try
    {
        for (int i = 0; i < bond_numbers; i++)
        {
            int atoms[2];
            const char* integer_fields[2] = {"atom A", "atom B"};
            for (int field = 0; field < 2; field++)
            {
                atoms[field] = Native_Torsion_Parse_Int(
                    Native_Torsion_Read_Token(&input, context,
                                              integer_fields[field], i),
                    context, integer_fields[field], i);
            }
            const double force_constant = Native_Torsion_Parse_Double(
                Native_Torsion_Read_Token(&input, context, "force constant",
                                          i),
                context, "force constant", i);
            const double equilibrium_distance = Native_Torsion_Parse_Double(
                Native_Torsion_Read_Token(&input, context,
                                          "equilibrium distance", i),
                context, "equilibrium distance", i);
            const float stored_force_constant = Native_Torsion_Checked_Float(
                force_constant, context, "force constant", "parameter", i);
            const float stored_equilibrium_distance =
                Native_Torsion_Checked_Float(
                    equilibrium_distance, context, "equilibrium distance",
                    "parameter", i);
            parsed.atom_a.push_back(atoms[0]);
            parsed.atom_b.push_back(atoms[1]);
            parsed.k.push_back(stored_force_constant);
            parsed.r0.push_back(stored_equilibrium_distance);
        }
    }
    catch (const std::bad_alloc&)
    {
        Native_Torsion_Throw_Storage_Error(context, bond_numbers);
        return;
    }
    catch (const std::length_error&)
    {
        Native_Torsion_Throw_Storage_Error(context, bond_numbers);
        return;
    }
    Native_Torsion_Ensure_End(&input, context);
    Native_Torsion_Close(&input, context);
    *bonds = std::move(parsed);
}

static void Native_Load_Bonds(System* system, CONTROLLER* controller)
{
    Native_Load_Bonds(&system->classical_force_field.bonds, controller);
}

}  // namespace Xponge
