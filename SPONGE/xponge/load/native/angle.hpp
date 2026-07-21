#pragma once

#include "torsion_parse.hpp"

namespace Xponge
{

static void Native_Load_Angles(Angles* angles, CONTROLLER* controller,
                               const char* module_name = "angle")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    const Native_Torsion_Input_Context context = {
        input_path, "angle_in_file", "Xponge::Native_Load_Angles", controller};
    std::ifstream input(input_path);
    if (!input.is_open())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Angles",
            "Reason:\n\tfailed to open angle_in_file\n\tInput file: %s\n",
            input_path);
    }
    const int angle_numbers = Native_Torsion_Parse_Int(
        Native_Torsion_Read_Token(&input, context, "interaction count", -1),
        context, "interaction count", -1);
    if (angle_numbers < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Native_Load_Angles",
            "Reason:\n\tangle_in_file has a negative interaction count\n"
            "\tInput file: %s\n",
            input_path);
    }
    Angles parsed;
    try
    {
        for (int i = 0; i < angle_numbers; i++)
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
            const double force_constant = Native_Torsion_Parse_Double(
                Native_Torsion_Read_Token(&input, context, "force constant",
                                          i),
                context, "force constant", i);
            const double equilibrium_angle = Native_Torsion_Parse_Double(
                Native_Torsion_Read_Token(&input, context,
                                          "equilibrium angle", i),
                context, "equilibrium angle", i);
            const float stored_force_constant = Native_Torsion_Checked_Float(
                force_constant, context, "force constant", "parameter", i);
            const float stored_equilibrium_angle = Native_Torsion_Checked_Float(
                equilibrium_angle, context, "equilibrium angle", "parameter",
                i);
            parsed.atom_a.push_back(atoms[0]);
            parsed.atom_b.push_back(atoms[1]);
            parsed.atom_c.push_back(atoms[2]);
            parsed.k.push_back(stored_force_constant);
            parsed.theta0.push_back(stored_equilibrium_angle);
        }
    }
    catch (const std::bad_alloc&)
    {
        Native_Torsion_Throw_Storage_Error(context, angle_numbers);
        return;
    }
    catch (const std::length_error&)
    {
        Native_Torsion_Throw_Storage_Error(context, angle_numbers);
        return;
    }
    Native_Torsion_Ensure_End(&input, context);
    Native_Torsion_Close(&input, context);
    *angles = std::move(parsed);
}

static void Native_Load_Angles(System* system, CONTROLLER* controller)
{
    Native_Load_Angles(&system->classical_force_field.angles, controller);
}

}  // namespace Xponge
