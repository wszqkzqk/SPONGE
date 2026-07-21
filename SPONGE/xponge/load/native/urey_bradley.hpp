#pragma once

#include "torsion_parse.hpp"

namespace Xponge
{

static void Native_Load_Urey_Bradley(UreyBradley* urey_bradley,
                                     CONTROLLER* controller,
                                     const char* module_name = "urey_bradley")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    const Native_Torsion_Input_Context context = {
        input_path, "urey_bradley_in_file", "Xponge::Native_Load_Urey_Bradley",
        controller};
    std::ifstream input(input_path);
    if (!input.is_open())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Urey_Bradley",
            "Reason:\n\tfailed to open urey_bradley_in_file\n"
            "\tInput file: %s\n",
            input_path);
    }
    const int term_numbers = Native_Torsion_Parse_Int(
        Native_Torsion_Read_Token(&input, context, "interaction count", -1),
        context, "interaction count", -1);
    if (term_numbers < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Urey_Bradley",
            "Reason:\n\turey_bradley_in_file has a negative interaction "
            "count\n\tInput file: %s\n",
            input_path);
    }

    UreyBradley parsed;
    try
    {
        for (int i = 0; i < term_numbers; i++)
        {
            int atoms[3];
            const char* integer_fields[3] = {"atom A", "atom B", "atom C"};
            for (int field = 0; field < 3; field++)
            {
                atoms[field] = Native_Torsion_Parse_Int(
                    Native_Torsion_Read_Token(&input, context,
                                              integer_fields[field], i),
                    context, integer_fields[field], i);
            }
            const char* real_fields[4] = {
                "angle force constant", "equilibrium angle",
                "bond force constant", "equilibrium distance"};
            float real_values[4];
            for (int field = 0; field < 4; field++)
            {
                const double parsed_value = Native_Torsion_Parse_Double(
                    Native_Torsion_Read_Token(&input, context,
                                              real_fields[field], i),
                    context, real_fields[field], i);
                real_values[field] = Native_Torsion_Checked_Float(
                    parsed_value, context, real_fields[field], "parameter", i);
            }
            parsed.atom_a.push_back(atoms[0]);
            parsed.atom_b.push_back(atoms[1]);
            parsed.atom_c.push_back(atoms[2]);
            parsed.angle_k.push_back(real_values[0]);
            parsed.angle_theta0.push_back(real_values[1]);
            parsed.bond_k.push_back(real_values[2]);
            parsed.bond_r0.push_back(real_values[3]);
        }
    }
    catch (const std::bad_alloc&)
    {
        Native_Torsion_Throw_Storage_Error(context, term_numbers);
        return;
    }
    catch (const std::length_error&)
    {
        Native_Torsion_Throw_Storage_Error(context, term_numbers);
        return;
    }
    Native_Torsion_Ensure_End(&input, context);
    Native_Torsion_Close(&input, context);
    *urey_bradley = std::move(parsed);
}

static void Native_Load_Urey_Bradley(System* system, CONTROLLER* controller)
{
    Native_Load_Urey_Bradley(&system->classical_force_field.urey_bradley,
                             controller);
}

}  // namespace Xponge
