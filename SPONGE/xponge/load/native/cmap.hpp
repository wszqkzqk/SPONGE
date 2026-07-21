#pragma once

#include <climits>
#include <fstream>

#include "../common.hpp"

namespace Xponge
{

static void Native_Load_CMap(CMap* cmap, CONTROLLER* controller,
                             const char* module_name = "cmap")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    const char* error_by = "Xponge::Native_Load_CMap";
    std::ifstream input(input_path);
    if (!input.is_open())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed to open cmap_in_file\n\tInput file: %s\n",
            input_path);
    }
    auto read_token = [&](const std::string& field) -> std::string
    {
        std::string token;
        if (!(input >> token))
        {
            if (input.bad())
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tI/O error while reading %s from "
                    "cmap_in_file\n\tInput file: %s\n",
                    field.c_str(), input_path);
            }
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tcmap_in_file is truncated while reading %s\n"
                "\tInput file: %s\n",
                field.c_str(), input_path);
        }
        return token;
    };
    auto parse_int = [&](const std::string& field) -> int
    {
        const std::string token = read_token(field);
        try
        {
            std::size_t consumed = 0;
            const long long value = std::stoll(token, &consumed, 10);
            if (consumed != token.size() ||
                value <
                    static_cast<long long>(std::numeric_limits<int>::min()) ||
                value > static_cast<long long>(std::numeric_limits<int>::max()))
            {
                throw std::out_of_range("not a strict signed int");
            }
            return static_cast<int>(value);
        }
        catch (const std::exception&)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tcmap_in_file %s token '%s' is not a strict "
                "signed integer in range\n\tInput file: %s\n",
                field.c_str(), token.c_str(), input_path);
        }
        return 0;
    };
    auto parse_float = [&](const std::string& field) -> float
    {
        const std::string token = read_token(field);
        double parsed = 0.0;
        try
        {
            std::size_t consumed = 0;
            parsed = std::stod(token, &consumed);
            if (consumed != token.size())
            {
                throw std::invalid_argument("not a strict real token");
            }
        }
        catch (const std::exception&)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tcmap_in_file %s token '%s' is invalid or "
                "outside the finite double range\n\tInput file: %s\n",
                field.c_str(), token.c_str(), input_path);
        }
        const double float_max =
            static_cast<double>(std::numeric_limits<float>::max());
        if (!Double_Memory_Is_Finite(&parsed) || parsed > float_max ||
            parsed < -float_max)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tcmap_in_file %s token '%s' is outside the "
                "finite float range\n\tInput file: %s\n",
                field.c_str(), token.c_str(), input_path);
        }
        const float stored = static_cast<float>(parsed);
        if (!Float_Memory_Is_Finite(&stored) ||
            (parsed != 0.0 && stored == 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tcmap_in_file nonzero %s token '%s' cannot be "
                "represented as a finite float\n\tInput file: %s\n",
                field.c_str(), token.c_str(), input_path);
        }
        if (!Float_Memory_Is_Zero_Or_Normal(&stored))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tcmap_in_file %s token '%s' is a subnormal "
                "float; native float fields require a finite zero or normal "
                "value for consistent FTZ behavior\n\tInput file: %s\n",
                field.c_str(), token.c_str(), input_path);
        }
        return stored;
    };

    CMap parsed_cmap;
    try
    {
        const int total_cmap_numbers = parse_int("interaction count");
        if (total_cmap_numbers < 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "Xponge::Native_Load_CMap",
                "Reason:\n\t%s contains a negative CMAP interaction count\n",
                input_path);
        }
        parsed_cmap.unique_type_numbers = parse_int("unique type count");
        if (parsed_cmap.unique_type_numbers < 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "Xponge::Native_Load_CMap",
                "Reason:\n\t%s contains a negative CMAP type count\n",
                input_path);
        }

        // Counts in a native file are untrusted.  Read and append each record
        // instead of resizing from the declarations so a truncated file with
        // a huge count fails on its first missing record without attempting a
        // huge allocation.
        for (int i = 0; i < parsed_cmap.unique_type_numbers; i++)
        {
            const int resolution =
                parse_int("resolution for CMAP type " + std::to_string(i));
            if (resolution <= 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "Xponge::Native_Load_CMap",
                    "Reason:\n\t%s contains non-positive resolution %d for "
                    "CMAP type %d\n",
                    input_path, resolution, i);
            }
            const long long type_gridpoint_numbers =
                static_cast<long long>(resolution) * resolution;
            if (type_gridpoint_numbers > INT_MAX / 16 ||
                parsed_cmap.unique_gridpoint_numbers >
                    INT_MAX / 16 - type_gridpoint_numbers)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "Xponge::Native_Load_CMap",
                    "Reason:\n\t%s CMAP grid size overflows the supported "
                    "range\n",
                    input_path);
            }
            parsed_cmap.resolution.push_back(resolution);
            parsed_cmap.type_offset.push_back(
                16 * parsed_cmap.unique_gridpoint_numbers);
            parsed_cmap.unique_gridpoint_numbers +=
                static_cast<int>(type_gridpoint_numbers);
        }

        for (int i = 0; i < parsed_cmap.unique_gridpoint_numbers; i++)
        {
            parsed_cmap.grid_value.push_back(
                parse_float("grid value " + std::to_string(i)));
        }

        for (int i = 0; i < total_cmap_numbers; i++)
        {
            int values[6] = {};
            const char* names[6] = {"atom A", "atom B", "atom C",
                                    "atom D", "atom E", "type"};
            for (int field = 0; field < 6; field++)
            {
                values[field] =
                    parse_int("interaction " + std::to_string(i) + " " +
                              names[field]);
            }
            parsed_cmap.atom_a.push_back(values[0]);
            parsed_cmap.atom_b.push_back(values[1]);
            parsed_cmap.atom_c.push_back(values[2]);
            parsed_cmap.atom_d.push_back(values[3]);
            parsed_cmap.atom_e.push_back(values[4]);
            parsed_cmap.cmap_type.push_back(values[5]);
        }
        std::string trailing;
        if (input >> trailing)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "Xponge::Native_Load_CMap",
                "Reason:\n\t%s contains trailing data beginning with '%s' "
                "after the CMAP records\n",
                input_path, trailing.c_str());
        }
        if (input.bad())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tI/O error while checking the end of "
                "cmap_in_file\n\tInput file: %s\n",
                input_path);
        }
        input.clear();
        input.close();
        if (input.fail())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tI/O error while closing cmap_in_file\n"
                "\tInput file: %s\n",
                input_path);
        }
    }
    catch (const std::length_error&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tcmap_in_file declarations exceed the maximum "
            "supported container size\n\tInput file: %s\n",
            input_path);
        return;
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tcould not allocate storage while reading "
            "cmap_in_file\n\tInput file: %s\n",
            input_path);
        return;
    }

    // Publish only after all records, EOF, and checked close pass.
    *cmap = std::move(parsed_cmap);
}

static void Native_Load_CMap(System* system, CONTROLLER* controller)
{
    Native_Load_CMap(&system->classical_force_field.cmap, controller);
}

}  // namespace Xponge
