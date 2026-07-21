#pragma once

#include "md_core_parse.hpp"

namespace Xponge
{

static void Native_Load_Generalized_Born(GeneralizedBorn* gb,
                                         CONTROLLER* controller,
                                         int atom_numbers_hint = 0,
                                         const char* module_name = "gb")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    Native_Core_Parser parser(input_path, "gb_in_file",
                              "Xponge::Native_Load_Generalized_Born",
                              controller);
    GeneralizedBorn parsed;
    const int atom_numbers =
        parser.Validate_Atom_Count(parser.Read_Int("atom count"), "atom count");
    if (atom_numbers_hint > 0 && atom_numbers_hint != atom_numbers)
    {
        parser.Fail(spongeErrorConflictingCommand,
                    "gb_in_file atom count " +
                        std::to_string(atom_numbers) +
                        " differs from the previously loaded atom count " +
                        std::to_string(atom_numbers_hint));
    }
    for (int i = 0; i < atom_numbers; i++)
    {
        parser.Append(&parsed.radius,
                      parser.Read_Float(
                          Native_Core_Entry_Field("GB radius", i)),
                      Native_Core_Entry_Field("GB radius", i));
        parser.Append(&parsed.scale_factor,
                      parser.Read_Float(
                          Native_Core_Entry_Field("GB scale factor", i)),
                      Native_Core_Entry_Field("GB scale factor", i));
    }
    parser.Ensure_End();
    parser.Close();
    *gb = std::move(parsed);
}

static void Native_Load_Generalized_Born(System* system, CONTROLLER* controller)
{
    Native_Load_Generalized_Born(&system->generalized_born, controller,
                                 Load_Get_Atom_Numbers(system));
}

}  // namespace Xponge
