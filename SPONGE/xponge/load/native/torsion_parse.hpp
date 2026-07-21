#pragma once

#include "../common.hpp"

namespace Xponge
{

struct Native_Torsion_Input_Context
{
    const char* input_path;
    const char* input_name;
    const char* error_by;
    CONTROLLER* controller;
};

static void Native_Torsion_Throw_Storage_Error(
    const Native_Torsion_Input_Context& context, int interaction_count)
{
    context.controller->Throw_Formatted_SPONGE_Error(
        spongeErrorMallocFailed, context.error_by,
        "Reason:\n\tfailed to allocate transactional storage for %d %s "
        "interactions\n\tInput file: %s\n",
        interaction_count, context.input_name, context.input_path);
}

static std::string Native_Torsion_Read_Token(
    std::istream* input, const Native_Torsion_Input_Context& context,
    const char* field_name, int interaction)
{
    std::string token;
    bool read_succeeded = false;
    try
    {
        read_succeeded = static_cast<bool>(*input >> token);
    }
    catch (const std::ios_base::failure&)
    {
        if (interaction < 0)
        {
            context.controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, context.error_by,
                "Reason:\n\tI/O error while reading %s from %s\n"
                "\tInput file: %s\n",
                field_name, context.input_name, context.input_path);
        }
        else
        {
            context.controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, context.error_by,
                "Reason:\n\tI/O error while reading %s interaction %d %s\n"
                "\tInput file: %s\n",
                context.input_name, interaction, field_name,
                context.input_path);
        }
    }
    if (!read_succeeded)
    {
        if (input->bad())
        {
            if (interaction < 0)
            {
                context.controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, context.error_by,
                    "Reason:\n\tI/O error while reading %s from %s\n"
                    "\tInput file: %s\n",
                    field_name, context.input_name, context.input_path);
            }
            else
            {
                context.controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, context.error_by,
                    "Reason:\n\tI/O error while reading %s interaction %d "
                    "%s\n\tInput file: %s\n",
                    context.input_name, interaction, field_name,
                    context.input_path);
            }
        }
        if (interaction < 0)
        {
            context.controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, context.error_by,
                "Reason:\n\tthe format of %s is not right: missing %s\n"
                "\tInput file: %s\n",
                context.input_name, field_name, context.input_path);
        }
        else
        {
            context.controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, context.error_by,
                "Reason:\n\tthe format of %s is not right: interaction %d "
                "is missing %s\n\tInput file: %s\n",
                context.input_name, interaction, field_name,
                context.input_path);
        }
    }
    return token;
}

static int Native_Torsion_Parse_Int(const std::string& token,
                                    const Native_Torsion_Input_Context& context,
                                    const char* field_name, int interaction)
{
    try
    {
        std::size_t consumed = 0;
        long long value = std::stoll(token, &consumed, 10);
        if (consumed != token.size() ||
            value < static_cast<long long>(std::numeric_limits<int>::min()) ||
            value > static_cast<long long>(std::numeric_limits<int>::max()))
        {
            throw std::out_of_range("not a strict signed int");
        }
        return static_cast<int>(value);
    }
    catch (const std::exception&)
    {
        if (interaction < 0)
        {
            context.controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, context.error_by,
                "Reason:\n\t%s %s token '%s' is not a strict signed integer "
                "in range\n\tInput file: %s\n",
                context.input_name, field_name, token.c_str(),
                context.input_path);
        }
        else
        {
            context.controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, context.error_by,
                "Reason:\n\t%s interaction %d %s token '%s' is not a strict "
                "signed integer in range\n\tInput file: %s\n",
                context.input_name, interaction, field_name, token.c_str(),
                context.input_path);
        }
    }
    return 0;
}

static double Native_Torsion_Parse_Double(
    const std::string& token, const Native_Torsion_Input_Context& context,
    const char* field_name, int interaction)
{
    double value = 0.0;
    try
    {
        std::size_t consumed = 0;
        value = std::stod(token, &consumed);
        if (consumed != token.size())
        {
            throw std::invalid_argument("not a strict floating-point token");
        }
    }
    catch (const std::exception&)
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\t%s interaction %d %s token '%s' is invalid or "
            "outside the finite double range\n\tInput file: %s\n",
            context.input_name, interaction, field_name, token.c_str(),
            context.input_path);
    }
    if (!Double_Memory_Is_Finite(&value))
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\t%s interaction %d has a non-finite parameter (%s "
            "token '%s')\n\tInput file: %s\n",
            context.input_name, interaction, field_name, token.c_str(),
            context.input_path);
    }
    return value;
}

static float Native_Torsion_Checked_Float(
    double value, const Native_Torsion_Input_Context& context,
    const char* value_name, const char* value_kind, int interaction)
{
    const double float_max =
        static_cast<double>(std::numeric_limits<float>::max());
    if (!Double_Memory_Is_Finite(&value) || value > float_max ||
        value < -float_max)
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\t%s interaction %d %s %s is outside the finite float "
            "range\n\tInput file: %s\n",
            context.input_name, interaction, value_name, value_kind,
            context.input_path);
    }
    float narrowed = static_cast<float>(value);
    if (!Float_Memory_Is_Finite(&narrowed))
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\t%s interaction %d has a non-finite %s %s\n"
            "\tInput file: %s\n",
            context.input_name, interaction, value_name, value_kind,
            context.input_path);
    }
    if (value != 0.0 && narrowed == 0.0f)
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\t%s interaction %d nonzero %s %s underflows the "
            "finite float range\n\tInput file: %s\n",
            context.input_name, interaction, value_name, value_kind,
            context.input_path);
    }
    if (!Float_Memory_Is_Zero_Or_Normal(&narrowed))
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\t%s interaction %d %s %s is a subnormal float; "
            "torsion parameters require a finite zero or normal float for "
            "consistent FTZ behavior\n\tInput file: %s\n",
            context.input_name, interaction, value_name, value_kind,
            context.input_path);
    }
    return narrowed;
}

static void Native_Torsion_Ensure_End(
    std::istream* input, const Native_Torsion_Input_Context& context)
{
    std::string trailing;
    bool has_trailing_data = false;
    try
    {
        has_trailing_data = static_cast<bool>(*input >> trailing);
    }
    catch (const std::ios_base::failure&)
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\tI/O error while checking the end of %s\n"
            "\tInput file: %s\n",
            context.input_name, context.input_path);
    }
    if (has_trailing_data)
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\t%s has trailing data beginning with '%s'\n"
            "\tInput file: %s\n",
            context.input_name, trailing.c_str(), context.input_path);
    }
    if (input->bad())
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\tI/O error while checking the end of %s\n"
            "\tInput file: %s\n",
            context.input_name, context.input_path);
    }
}

static void Native_Torsion_Close(
    std::ifstream* input, const Native_Torsion_Input_Context& context)
{
    // Ensure_End consumes the normal EOF probe and therefore leaves eofbit
    // and failbit set.  Clear only after that state has been checked so close
    // failures remain distinguishable from a successful EOF check.
    input->clear();
    input->close();
    if (input->fail())
    {
        context.controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, context.error_by,
            "Reason:\n\tI/O error while closing %s\n"
            "\tInput file: %s\n",
            context.input_name, context.input_path);
    }
}

}  // namespace Xponge
