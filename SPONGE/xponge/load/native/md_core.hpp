#pragma once

#include "md_core_parse.hpp"

namespace Xponge
{

static std::vector<std::string> Native_Core_Read_Header(
    Native_Core_Parser* parser, const char* input_name, int maximum_fields)
{
    std::vector<std::string> fields = parser->Read_Line_Tokens("header");
    if (fields.empty())
    {
        parser->Fail(
            spongeErrorBadFileFormat,
            std::string(input_name) + " header is missing the atom count");
    }
    if (fields.size() > static_cast<std::size_t>(maximum_fields))
    {
        parser->Fail(spongeErrorBadFileFormat,
                     std::string(input_name) + " header has unexpected field " +
                         std::to_string(maximum_fields + 1) + " ('" +
                         fields[maximum_fields] + "')");
    }
    return fields;
}

static int Native_Core_Parse_Atom_Header(Native_Core_Parser* parser,
                                         const std::vector<std::string>& fields,
                                         System* system, double* start_time)
{
    const int atom_count = parser->Validate_Atom_Count(
        parser->Parse_Int(fields[0], "header atom count"), "atom count");
    parser->Ensure_Atom_Count_Matches(system, atom_count, "header atom count");

    if (fields.size() >= 2)
    {
        const double parsed_time =
            parser->Parse_Double(fields[1], "header start time");
        if (start_time != nullptr) *start_time = parsed_time;
    }
    if (fields.size() >= 3)
    {
        parser->Parse_Int(fields[2], "header step");
    }
    if (start_time != nullptr && fields.size() < 2) *start_time = 0.0;
    return atom_count;
}

static void Native_Load_Mass(System* system, CONTROLLER* controller)
{
    if (controller->Command_Exist("mass_in_file"))
    {
        Native_Core_Parser parser(controller->Original_Command("mass_in_file"),
                                  "mass_in_file", "Xponge::Native_Load_Mass",
                                  controller);
        const int atom_count = parser.Validate_Atom_Count(
            parser.Read_Int("atom count"), "atom count");
        parser.Ensure_Atom_Count_Matches(system, atom_count, "atom count");

        std::vector<float> mass;
        for (int atom = 0; atom < atom_count; atom++)
        {
            const std::string field = Native_Core_Entry_Field("mass", atom);
            parser.Append(&mass, parser.Read_Float(field), field);
        }
        parser.Ensure_End();
        parser.Close();
        system->atoms.mass = std::move(mass);
        return;
    }

    const int atom_count = Load_Get_Atom_Numbers(system);
    if (atom_count > 0)
    {
        system->atoms.mass.assign(atom_count, 20.0f);
        return;
    }

    controller->Throw_SPONGE_Error(
        spongeErrorMissingCommand, "Xponge::Native_Load_Mass",
        "Reason:\n\tno atom_numbers found. No mass_in_file is provided\n");
}

static void Native_Load_Charge(System* system, CONTROLLER* controller)
{
    if (controller->Command_Exist("charge_in_file"))
    {
        Native_Core_Parser parser(
            controller->Original_Command("charge_in_file"),
                                  "charge_in_file",
                                  "Xponge::Native_Load_Charge", controller);
        const int atom_count = parser.Validate_Atom_Count(
            parser.Read_Int("atom count"), "atom count");
        parser.Ensure_Atom_Count_Matches(system, atom_count, "atom count");

        std::vector<float> charge;
        for (int atom = 0; atom < atom_count; atom++)
        {
            const std::string field = Native_Core_Entry_Field("charge", atom);
            parser.Append(&charge, parser.Read_Float(field), field);
        }
        parser.Ensure_End();
        parser.Close();
        system->atoms.charge = std::move(charge);
        return;
    }

    const int atom_count = Load_Get_Atom_Numbers(system);
    if (atom_count > 0)
    {
        system->atoms.charge.assign(atom_count, 0.0f);
        return;
    }

    controller->Throw_SPONGE_Error(
        spongeErrorMissingCommand, "Xponge::Native_Load_Charge",
        "Reason:\n\tno atom_numbers found. No charge_in_file is provided\n");
}

static void Native_Load_Coordinate_And_Velocity(System* system,
                                                CONTROLLER* controller)
{
    const char* error_by = "Xponge::Native_Load_Coordinate_And_Velocity";
    if (!controller->Command_Exist("coordinate_in_file"))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand, error_by,
            "Reason:\n\tno coordinate information found");
    }

    Native_Core_Parser coordinate_parser(
        controller->Original_Command("coordinate_in_file"),
        "coordinate_in_file",
        error_by, controller);
    const std::vector<std::string> coordinate_header =
        Native_Core_Read_Header(&coordinate_parser, "coordinate_in_file", 3);
    double start_time = 0.0;
    const int atom_count = Native_Core_Parse_Atom_Header(
        &coordinate_parser, coordinate_header, system, &start_time);
    const std::size_t component_count =
        static_cast<std::size_t>(atom_count) * 3;

    std::vector<float> coordinate;
    const char* axes[3] = {"x", "y", "z"};
    for (int atom = 0; atom < atom_count; atom++)
    {
        for (int axis = 0; axis < 3; axis++)
        {
            const std::string field = "coordinate " +
                                      std::string(axes[axis]) + " entry " +
                                      std::to_string(atom);
            coordinate_parser.Append(
                &coordinate, coordinate_parser.Read_Float(field), field);
        }
    }

    std::vector<float> box_length(3);
    std::vector<float> box_angle(3);
    for (int axis = 0; axis < 3; axis++)
    {
        box_length[axis] = coordinate_parser.Read_Float(
            "box length " + std::string(axes[axis]));
    }
    const char* angles[3] = {"alpha", "beta", "gamma"};
    for (int axis = 0; axis < 3; axis++)
    {
        box_angle[axis] = coordinate_parser.Read_Float(
            "box angle " + std::string(angles[axis]));
    }
    coordinate_parser.Ensure_End();
    coordinate_parser.Close();

    std::vector<float> velocity;
    if (controller->Command_Exist("velocity_in_file"))
    {
        Native_Core_Parser velocity_parser(
            controller->Original_Command("velocity_in_file"),
            "velocity_in_file",
            error_by, controller);
        const std::vector<std::string> velocity_header =
            Native_Core_Read_Header(&velocity_parser, "velocity_in_file", 3);
        const int velocity_atom_count = Native_Core_Parse_Atom_Header(
            &velocity_parser, velocity_header, system, nullptr);

        for (int atom = 0; atom < velocity_atom_count; atom++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                const std::string field =
                    "velocity " + std::string(axes[axis]) + " entry " +
                    std::to_string(atom);
                velocity_parser.Append(
                    &velocity, velocity_parser.Read_Float(field), field);
            }
        }
        velocity_parser.Ensure_End();
        velocity_parser.Close();
    }
    else
    {
        coordinate_parser.Assign(&velocity, component_count, 0.0f,
                                 "default velocity components");
    }

    // Coordinate, box, time, and velocity describe one restart state.  Do not
    // publish any part of that state until both input files have been parsed
    // and their trailing-data checks have succeeded.
    system->start_time = start_time;
    system->atoms.coordinate = std::move(coordinate);
    system->atoms.velocity = std::move(velocity);
    system->box.box_length = std::move(box_length);
    system->box.box_angle = std::move(box_angle);
}

static void Native_Load_Residues(System* system, CONTROLLER* controller)
{
    const int atom_count = Load_Get_Atom_Numbers(system);
    if (!controller->Command_Exist("residue_in_file"))
    {
        system->residues.atom_numbers.assign(atom_count, 1);
        return;
    }

    Native_Core_Parser parser(controller->Original_Command("residue_in_file"),
                              "residue_in_file", "Xponge::Native_Load_Residues",
                              controller);
    const std::vector<std::string> header = parser.Read_Line_Tokens("header");
    if (header.size() != 2)
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "residue_in_file header must contain exactly the atom "
                    "count and residue count");
    }
    const int residue_atom_count = parser.Validate_Atom_Count(
        parser.Parse_Int(header[0], "header atom count"), "atom count");
    parser.Ensure_Atom_Count_Matches(system, residue_atom_count,
                                     "header atom count");
    const int residue_count =
        parser.Parse_Int(header[1], "header residue count");
    if (residue_count < 0)
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "residue_in_file has a negative residue count");
    }
    if (residue_count > residue_atom_count)
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "residue_in_file residue count cannot fit strictly "
                    "positive residue lengths within the atom count");
    }

    std::vector<int> residue_lengths;
    std::size_t length_sum = 0;
    for (int residue = 0; residue < residue_count; residue++)
    {
        const std::string field =
            Native_Core_Entry_Field("residue length", residue);
        const int length = parser.Read_Int(field);
        if (length <= 0)
        {
            parser.Fail(
                spongeErrorBadFileFormat,
                "residue_in_file " + field + " must be strictly positive");
        }
        if (static_cast<std::size_t>(length) >
            std::numeric_limits<std::size_t>::max() - length_sum)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "residue_in_file residue-length sum overflows the "
                        "supported size range at entry " +
                            std::to_string(residue));
        }
        length_sum += static_cast<std::size_t>(length);
        parser.Append(&residue_lengths, length, field);
    }
    parser.Ensure_End();
    if (length_sum != static_cast<std::size_t>(residue_atom_count))
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "residue_in_file residue-length sum " +
                        std::to_string(length_sum) +
                        " does not equal atom count " +
                        std::to_string(residue_atom_count));
    }
    parser.Close();
    system->residues.atom_numbers = std::move(residue_lengths);
}

static void Native_Load_Exclusions(System* system, CONTROLLER* controller)
{
    const int atom_count = Load_Get_Atom_Numbers(system);
    if (!controller->Command_Exist("exclude_in_file"))
    {
        system->exclusions.excluded_atoms.assign(atom_count, {});
        return;
    }

    Native_Core_Parser parser(controller->Original_Command("exclude_in_file"),
                              "exclude_in_file",
                              "Xponge::Native_Load_Exclusions", controller);
    const std::vector<std::string> header = parser.Read_Line_Tokens("header");
    if (header.size() != 2)
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "exclude_in_file header must contain exactly the atom "
                    "count and declared exclusion total");
    }
    const int exclusion_atom_count = parser.Validate_Atom_Count(
        parser.Parse_Int(header[0], "header atom count"), "atom count");
    parser.Ensure_Atom_Count_Matches(system, exclusion_atom_count,
                                     "header atom count");
    const int declared_total =
        parser.Parse_Int(header[1], "header declared exclusion total");
    if (declared_total < 0)
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "exclude_in_file has a negative declared exclusion "
                    "total");
    }
    const std::uint64_t maximum_total =
        static_cast<std::uint64_t>(exclusion_atom_count) *
        static_cast<std::uint64_t>(
            exclusion_atom_count > 0 ? exclusion_atom_count - 1 : 0) /
        2;
    if (static_cast<std::uint64_t>(declared_total) > maximum_total)
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "exclude_in_file declared exclusion total exceeds the "
                    "triangular-format maximum for its atom count");
    }

    std::vector<std::vector<int>> excluded_atoms;
    std::size_t parsed_total = 0;
    for (int atom = 0; atom < atom_count; atom++)
    {
        const std::string row_field = "exclusion row " + std::to_string(atom);
        const std::vector<std::string> row = parser.Read_Line_Tokens(row_field);
        if (row.empty())
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "exclude_in_file " + row_field +
                            " is missing its exclusion count");
        }
        const int row_count = parser.Parse_Int(row[0], row_field + " count");
        if (row_count < 0)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "exclude_in_file " + row_field +
                            " has a negative exclusion count");
        }
        const int maximum_row_count = atom_count - atom - 1;
        if (row_count > maximum_row_count)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "exclude_in_file " + row_field + " count " +
                            std::to_string(row_count) +
                            " exceeds the remaining upper-triangle atoms");
        }
        if (row.size() != static_cast<std::size_t>(row_count) + 1)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "exclude_in_file " + row_field + " declares " +
                            std::to_string(row_count) + " IDs but contains " +
                            std::to_string(row.size() - 1));
        }
        if (static_cast<std::size_t>(row_count) >
            std::numeric_limits<std::size_t>::max() - parsed_total)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "exclude_in_file exclusion-count sum overflows the "
                        "supported size range at row " +
                            std::to_string(atom));
        }
        parsed_total += static_cast<std::size_t>(row_count);
        if (parsed_total > static_cast<std::size_t>(declared_total))
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "exclude_in_file exclusion-count sum exceeds the "
                        "declared total at row " +
                            std::to_string(atom));
        }

        std::vector<int> row_neighbors;
        int previous_neighbor = atom;
        for (int neighbor_entry = 0; neighbor_entry < row_count;
             neighbor_entry++)
        {
            const std::string field =
                row_field + " neighbor entry " + std::to_string(neighbor_entry);
            const int neighbor =
                parser.Parse_Int(row[neighbor_entry + 1], field);
            if (neighbor < 0 || neighbor >= atom_count)
            {
                parser.Fail(spongeErrorBadFileFormat,
                            "exclude_in_file " + field + " ID " +
                                std::to_string(neighbor) +
                                " is outside atom range [0, " +
                                std::to_string(atom_count) + ")");
            }
            if (neighbor == atom)
            {
                parser.Fail(
                    spongeErrorBadFileFormat,
                    "exclude_in_file " + field + " is a self exclusion");
            }
            if (neighbor < atom)
            {
                parser.Fail(spongeErrorBadFileFormat,
                            "exclude_in_file " + field +
                                " is a reverse/lower-triangle exclusion; "
                                "neighbor IDs must be greater than the row "
                                "atom");
            }
            if (neighbor <= previous_neighbor)
            {
                const char* reason = neighbor == previous_neighbor
                                         ? "duplicates the previous ID"
                                         : "is not strictly increasing";
                parser.Fail(spongeErrorBadFileFormat,
                            "exclude_in_file " + field + " " + reason);
            }
            parser.Append(&row_neighbors, neighbor, field);
            previous_neighbor = neighbor;
        }
        parser.Append(&excluded_atoms, std::move(row_neighbors), row_field);
    }
    parser.Ensure_End();
    if (parsed_total != static_cast<std::size_t>(declared_total))
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "exclude_in_file exclusion-count sum " +
                        std::to_string(parsed_total) +
                        " does not equal declared total " +
                        std::to_string(declared_total));
    }
    parser.Close();
    system->exclusions.excluded_atoms = std::move(excluded_atoms);
}

}  // namespace Xponge
