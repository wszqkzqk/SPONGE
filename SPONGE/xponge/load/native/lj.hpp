#pragma once

#include "md_core_parse.hpp"

namespace Xponge
{

static void Native_Load_LJ(LennardJones* lj, CONTROLLER* controller,
                           int atom_numbers_hint = 0,
                           const char* module_name = "LJ")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    Native_Core_Parser parser(input_path, "LJ_in_file",
                              "Xponge::Native_Load_LJ", controller);
    LennardJones parsed;
    const int atom_numbers =
        parser.Validate_Atom_Count(parser.Read_Int("atom count"), "atom count");
    const int atom_type_numbers = parser.Read_Int("atom type count");
    const int pair_type_numbers = parser.Validate_Triangular_Type_Count(
        atom_type_numbers, "atom type count");
    if (atom_numbers > 0 && atom_type_numbers == 0)
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "LJ_in_file has atoms but no LJ atom types");
    }
    if (atom_numbers_hint > 0 && atom_numbers_hint != atom_numbers)
    {
        parser.Fail(spongeErrorConflictingCommand,
                    "LJ_in_file atom count " + std::to_string(atom_numbers) +
                        " differs from the previously loaded atom count " +
                        std::to_string(atom_numbers_hint));
    }
    parsed.atom_type_numbers = atom_type_numbers;
    for (int i = 0; i < pair_type_numbers; i++)
    {
        const std::string field = Native_Core_Entry_Field("LJ A", i);
        parsed.pair_A.push_back(parser.Checked_Float(
            parser.Read_Double(field) * 12.0, field + " scaled by 12"));
    }
    for (int i = 0; i < pair_type_numbers; i++)
    {
        const std::string field = Native_Core_Entry_Field("LJ B", i);
        parsed.pair_B.push_back(parser.Checked_Float(
            parser.Read_Double(field) * 6.0, field + " scaled by 6"));
    }
    for (int i = 0; i < atom_numbers; i++)
    {
        const std::string field = Native_Core_Entry_Field("atom LJ type", i);
        const int atom_type = parser.Read_Int(field);
        if (atom_type < 0 || atom_type >= atom_type_numbers)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "LJ_in_file " + field + " " +
                            std::to_string(atom_type) + " is outside [0, " +
                            std::to_string(atom_type_numbers) + ")");
        }
        parsed.atom_type.push_back(atom_type);
    }
    parser.Ensure_End();
    parser.Close();
    *lj = std::move(parsed);
}

static void Native_Load_LJ(System* system, CONTROLLER* controller)
{
    Native_Load_LJ(&system->classical_force_field.lj, controller,
                   Load_Get_Atom_Numbers(system));
}

}  // namespace Xponge
