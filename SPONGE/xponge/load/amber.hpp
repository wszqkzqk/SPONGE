#pragma once

#include "../xponge.h"
#include "./common.hpp"

namespace Xponge
{

static std::string Amber_Trim(const std::string& value)
{
    std::size_t begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
    {
        return "";
    }
    std::size_t end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

static bool Amber_Is_Finite_Double(double value)
{
    return Double_Memory_Is_Finite(&value);
}

static bool Amber_Is_Finite_Float(float value)
{
    return Float_Memory_Is_Finite(&value);
}

static bool Amber_Is_Strict_Float_Lexeme(const std::string& value)
{
    if (value.empty()) return false;

    std::size_t position = 0;
    if (value[position] == '+' || value[position] == '-') position++;

    bool has_digit = false;
    while (position < value.size() && value[position] >= '0' &&
           value[position] <= '9')
    {
        has_digit = true;
        position++;
    }
    if (position < value.size() && value[position] == '.')
    {
        position++;
        while (position < value.size() && value[position] >= '0' &&
               value[position] <= '9')
        {
            has_digit = true;
            position++;
        }
    }
    if (!has_digit) return false;

    if (position < value.size() &&
        (value[position] == 'E' || value[position] == 'e'))
    {
        position++;
        if (position < value.size() &&
            (value[position] == '+' || value[position] == '-'))
        {
            position++;
        }
        std::size_t exponent_begin = position;
        while (position < value.size() && value[position] >= '0' &&
               value[position] <= '9')
        {
            position++;
        }
        if (position == exponent_begin) return false;
    }
    return position == value.size();
}

static int Amber_Parse_Int(const std::string& value, const char* section_name,
                           CONTROLLER* controller,
                           const char* error_by = "Xponge::Amber_Reader")
{
    std::string trimmed = Amber_Trim(value);
    try
    {
        std::size_t consumed = 0;
        long long parsed = std::stoll(trimmed, &consumed, 10);
        if (trimmed.empty() || consumed != trimmed.size() ||
            parsed < std::numeric_limits<int>::min() ||
            parsed > std::numeric_limits<int>::max())
        {
            throw std::invalid_argument("not a strict int");
        }
        return static_cast<int>(parsed);
    }
    catch (const std::exception&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s contains invalid integer field "
            "'%s'\n",
            section_name, trimmed.c_str());
    }
    return 0;
}

static double Amber_Parse_Double(const std::string& value,
                                 const char* section_name,
                                 CONTROLLER* controller,
                                 const char* error_by = "Xponge::Amber_Reader")
{
    std::string normalized = Amber_Trim(value);
    for (char& character : normalized)
    {
        if (character == 'D' || character == 'd')
        {
            character = 'E';
        }
    }
    if (!Amber_Is_Strict_Float_Lexeme(normalized))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s contains invalid or non-finite "
            "floating-point field '%s'\n",
            section_name, normalized.c_str());
    }
    try
    {
        std::size_t consumed = 0;
        double parsed = std::stod(normalized, &consumed);
        if (normalized.empty() || consumed != normalized.size() ||
            !Amber_Is_Finite_Double(parsed))
        {
            throw std::invalid_argument("not a strict finite double");
        }
        return parsed;
    }
    catch (const std::exception&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s contains invalid or non-finite "
            "floating-point field '%s'\n",
            section_name, normalized.c_str());
    }
    return 0.0;
}

static float Amber_Checked_Float(double value, const char* section_name,
                                 CONTROLLER* controller,
                                 const char* error_by = "Xponge::Amber_Reader",
                                 const char* source_context = "")
{
    double float_max = static_cast<double>(std::numeric_limits<float>::max());
    if (!Amber_Is_Finite_Double(value) || value > float_max ||
        value < -float_max)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s contains a value outside the finite "
            "float range\n%s",
            section_name, source_context);
    }
    float parsed = static_cast<float>(value);
    if (!Amber_Is_Finite_Float(parsed))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s overflows the finite float range\n%s",
            section_name, source_context);
    }
    if (!Float_Memory_Is_Zero_Or_Normal(&parsed))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s contains a subnormal value that "
            "is not representable as a finite zero or normal float\n%s",
            section_name, source_context);
    }
    if (value != 0.0 && parsed == 0.0f)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s contains a nonzero value that "
            "underflows the finite float range\n%s",
            section_name, source_context);
    }
    return parsed;
}

static float Amber_Parse_Float(const std::string& value,
                               const char* section_name, CONTROLLER* controller,
                               const char* error_by = "Xponge::Amber_Reader",
                               const char* source_context = "")
{
    return Amber_Checked_Float(
        Amber_Parse_Double(value, section_name, controller, error_by),
        section_name, controller, error_by, source_context);
}

static float Amber_Checked_Scale(float value, double multiplier,
                                 const char* section_name,
                                 CONTROLLER* controller, const char* error_by,
                                 const char* source_context = "")
{
    double scaled = static_cast<double>(value) * multiplier;
    return Amber_Checked_Float(scaled, section_name, controller, error_by,
                               source_context);
}

struct Amber_Format_Field
{
    char kind = 0;
    std::size_t width = 0;
};

struct Amber_Format_Instruction
{
    enum class Kind
    {
        kField,
        kGroupBegin,
        kGroupEnd
    };

    Kind kind = Kind::kField;
    std::size_t repeat = 1;
    Amber_Format_Field field;
};

using Amber_Format = std::vector<Amber_Format_Instruction>;

struct Amber_Section
{
    std::string name;
    std::vector<std::string> values;
    std::vector<std::string> raw_values;
    std::size_t flag_line = 0;
};

struct Amber_Parm7_Document
{
    std::string path;
    std::map<std::string, Amber_Section> sections;
    std::vector<int> pointers;
};

static std::string Amber_Source_Context(
    const Amber_Parm7_Document& document,
    const std::vector<const char*>& section_names)
{
    std::string context = "\tInput file: " + document.path + "\n";
    for (const char* name : section_names)
    {
        const auto found = document.sections.find(name);
        if (found == document.sections.end())
        {
            context += "\t%FLAG " + std::string(name) + " is absent\n";
        }
        else
        {
            context += "\t%FLAG " + std::string(name) + " begins on line " +
                       std::to_string(found->second.flag_line) + "\n";
        }
    }
    return context;
}

static std::size_t Amber_Parse_Format_Size(
    const std::string& sequence, std::size_t* position, bool* has_digits,
    const char* quantity, CONTROLLER* controller, const char* input_path,
    const char* section_name, std::size_t flag_line, std::size_t line_number)
{
    const char* error_by = "Xponge::Amber_Read_Parm7";
    std::size_t parsed = 0;
    *has_digits = false;
    while (*position < sequence.size() &&
           isdigit(static_cast<unsigned char>(sequence[*position])))
    {
        *has_digits = true;
        const std::size_t digit =
            static_cast<std::size_t>(sequence[*position] - '0');
        if (parsed > (std::numeric_limits<std::size_t>::max() - digit) / 10)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER %%FORMAT %s exceeds the platform size "
                "range on line %zu\n\tInput file: %s\n"
                "\t%%FLAG %s begins on line %zu\n",
                quantity, line_number, input_path, section_name, flag_line);
            return 0;
        }
        parsed = 10 * parsed + digit;
        (*position)++;
    }
    return parsed;
}

static Amber_Format Amber_Parse_Format(
    const std::string& format_line, CONTROLLER* controller,
    const char* input_path, const char* section_name, std::size_t flag_line,
    std::size_t line_number)
{
    const char* error_by = "Xponge::Amber_Read_Parm7";
    const std::string trimmed = Amber_Trim(format_line);
    const std::string prefix = "%FORMAT(";
    if (trimmed.rfind(prefix, 0) != 0 || trimmed.size() <= prefix.size() ||
        trimmed.back() != ')')
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tinvalid or missing AMBER %%FORMAT on line %zu\n",
            line_number);
    }
    const std::string sequence =
        trimmed.substr(prefix.size(), trimmed.size() - prefix.size() - 1);

    struct Group_Parse_State
    {
        bool has_field = false;
    };

    Amber_Format format;
    std::vector<Group_Parse_State> groups;
    bool top_level_has_field = false;
    std::size_t position = 0;
    while (position < sequence.size())
    {
        while (position < sequence.size() &&
               (sequence[position] == ' ' || sequence[position] == '\t' ||
                sequence[position] == ','))
        {
            position++;
        }
        if (position >= sequence.size()) break;

        if (sequence[position] == ')')
        {
            if (groups.empty())
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tunexpected ')' in AMBER %%FORMAT on line "
                    "%zu\n\tInput file: %s\n\t%%FLAG %s begins on line "
                    "%zu\n",
                    line_number, input_path, section_name, flag_line);
            }
            position++;
            const bool group_has_field = groups.back().has_field;
            groups.pop_back();
            if (!group_has_field)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tempty grouped AMBER %%FORMAT on line %zu\n"
                    "\tInput file: %s\n\t%%FLAG %s begins on line %zu\n",
                    line_number, input_path, section_name, flag_line);
            }
            Amber_Format_Instruction end;
            end.kind = Amber_Format_Instruction::Kind::kGroupEnd;
            format.push_back(end);
            if (groups.empty())
            {
                top_level_has_field = true;
            }
            else
            {
                groups.back().has_field = true;
            }
            continue;
        }

        bool negative = false;
        if (sequence[position] == '+' || sequence[position] == '-')
        {
            negative = sequence[position] == '-';
            position++;
            if (position >= sequence.size() ||
                !isdigit(static_cast<unsigned char>(sequence[position])))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tinvalid signed AMBER %%FORMAT descriptor on "
                    "line %zu\n\tInput file: %s\n\t%%FLAG %s begins on "
                    "line %zu\n",
                    line_number, input_path, section_name, flag_line);
            }
        }

        bool has_repeat = false;
        std::size_t repeat = Amber_Parse_Format_Size(
            sequence, &position, &has_repeat, "repeat", controller, input_path,
            section_name, flag_line, line_number);

        if (position < sequence.size() &&
            (sequence[position] == 'P' || sequence[position] == 'p'))
        {
            // A Fortran scale descriptor affects numeric rendering but not
            // the fixed-width record layout used to read a parm7 file.
            position++;
            continue;
        }
        if (negative)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER %%FORMAT field repeats cannot be negative "
                "on line %zu\n\tInput file: %s\n\t%%FLAG %s begins on line "
                "%zu\n",
                line_number, input_path, section_name, flag_line);
        }
        if (!has_repeat) repeat = 1;
        if (repeat == 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER %%FORMAT field repeats must be positive "
                "on line %zu\n\tInput file: %s\n\t%%FLAG %s begins on line "
                "%zu\n",
                line_number, input_path, section_name, flag_line);
        }

        if (position < sequence.size() && sequence[position] == '(')
        {
            position++;
            Amber_Format_Instruction begin;
            begin.kind = Amber_Format_Instruction::Kind::kGroupBegin;
            begin.repeat = repeat;
            format.push_back(begin);
            groups.push_back({});
            continue;
        }

        if (position >= sequence.size())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tincomplete AMBER %%FORMAT on line %zu\n",
                line_number);
        }
        const char kind = static_cast<char>(
            toupper(static_cast<unsigned char>(sequence[position])));
        position++;
        if (kind != 'I' && kind != 'E' && kind != 'D' && kind != 'F' &&
            kind != 'G' && kind != 'A')
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tunsupported AMBER %%FORMAT descriptor '%c' on "
                "line %zu\n\tInput file: %s\n\t%%FLAG %s begins on line "
                "%zu\n",
                kind, line_number, input_path, section_name, flag_line);
        }

        bool has_width = false;
        const std::size_t width = Amber_Parse_Format_Size(
            sequence, &position, &has_width, "width", controller, input_path,
            section_name, flag_line, line_number);
        if (!has_width || width == 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER %%FORMAT field widths must be positive on "
                "line %zu\n\tInput file: %s\n\t%%FLAG %s begins on line "
                "%zu\n",
                line_number, input_path, section_name, flag_line);
        }

        if (position < sequence.size() && sequence[position] == '.')
        {
            position++;
            if (position >= sequence.size() ||
                !isdigit(static_cast<unsigned char>(sequence[position])))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tinvalid AMBER %%FORMAT precision on line "
                    "%zu\n",
                    line_number);
            }
            while (position < sequence.size() &&
                   isdigit(static_cast<unsigned char>(sequence[position])))
            {
                position++;
            }
            // Fortran permits an optional exponent-width suffix (for example
            // E16.8E3). It does not change the fixed field width.
            if (position < sequence.size() &&
                (sequence[position] == 'E' || sequence[position] == 'e'))
            {
                position++;
                if (position >= sequence.size() ||
                    !isdigit(static_cast<unsigned char>(sequence[position])))
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, error_by,
                        "Reason:\n\tinvalid AMBER %%FORMAT exponent width on "
                        "line %zu\n\tInput file: %s\n\t%%FLAG %s begins on "
                        "line %zu\n",
                        line_number, input_path, section_name, flag_line);
                }
                while (position < sequence.size() &&
                       isdigit(static_cast<unsigned char>(sequence[position])))
                {
                    position++;
                }
            }
        }

        Amber_Format_Instruction field;
        field.kind = Amber_Format_Instruction::Kind::kField;
        field.repeat = repeat;
        field.field = {kind, width};
        format.push_back(field);
        if (groups.empty())
        {
            top_level_has_field = true;
        }
        else
        {
            groups.back().has_field = true;
        }
    }

    if (!groups.empty())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tunterminated grouped AMBER %%FORMAT on line %zu\n"
            "\tInput file: %s\n\t%%FLAG %s begins on line %zu\n",
            line_number, input_path, section_name, flag_line);
    }
    if (!top_level_has_field)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tempty AMBER %%FORMAT on line %zu\n", line_number);
    }
    return format;
}

template <typename Visitor>
static bool Amber_For_Each_Format_Field(const Amber_Format& format,
                                        Visitor&& visitor)
{
    struct Group_Frame
    {
        std::size_t begin_instruction;
        std::size_t repeats_remaining;
    };

    std::vector<Group_Frame> groups;
    std::size_t instruction_index = 0;
    while (instruction_index < format.size())
    {
        const Amber_Format_Instruction& instruction =
            format[instruction_index];
        if (instruction.kind == Amber_Format_Instruction::Kind::kField)
        {
            for (std::size_t repeats_remaining = instruction.repeat;
                 repeats_remaining > 0; repeats_remaining--)
            {
                if (!visitor(instruction.field)) return false;
            }
            instruction_index++;
        }
        else if (instruction.kind ==
                 Amber_Format_Instruction::Kind::kGroupBegin)
        {
            groups.push_back({instruction_index, instruction.repeat});
            instruction_index++;
        }
        else
        {
            Group_Frame& group = groups.back();
            group.repeats_remaining--;
            if (group.repeats_remaining > 0)
            {
                instruction_index = group.begin_instruction + 1;
            }
            else
            {
                groups.pop_back();
                instruction_index++;
            }
        }
    }
    return true;
}

static Amber_Parm7_Document Amber_Read_Parm7_Document(const char* path,
                                                      CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Read_Parm7";
    std::ifstream fin(path);
    if (!fin.is_open())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed to open amber_parm7\n");
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fin, line))
    {
        lines.push_back(line);
    }
    if (fin.bad())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed while reading amber_parm7 '%s'\n", path);
    }
    fin.clear();
    fin.close();
    if (fin.fail())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed while closing amber_parm7 '%s'\n", path);
    }

    Amber_Parm7_Document document;
    document.path = path;
    bool version_seen = false;
    bool flag_seen = false;
    for (std::size_t i = 0; i < lines.size();)
    {
        if (lines[i].rfind("%FLAG", 0) != 0)
        {
            const std::string trimmed = Amber_Trim(lines[i]);
            if (trimmed.empty() || lines[i].rfind("%COMMENT", 0) == 0)
            {
                i++;
                continue;
            }
            if (lines[i].rfind("%VERSION", 0) == 0)
            {
                if (version_seen || flag_seen)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, error_by,
                        "Reason:\n\tAMBER %%VERSION must occur at most once "
                        "before the first %%FLAG (line %zu)\n\tInput file: "
                        "%s\n",
                        i + 1, path);
                }
                version_seen = true;
                i++;
                continue;
            }
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tunexpected content outside an AMBER %%FLAG "
                "section on line %zu\n\tInput file: %s\n",
                i + 1, path);
        }
        flag_seen = true;
        std::string name = Amber_Trim(lines[i].substr(5));
        std::size_t flag_line = i + 1;
        if (name.empty())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tempty AMBER %%FLAG on line %zu\n", flag_line);
        }
        if (document.sections.count(name) != 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER singleton section %s occurs more than "
                "once\n",
                name.c_str());
        }
        i++;
        while (i < lines.size() && lines[i].rfind("%COMMENT", 0) == 0)
        {
            i++;
        }
        if (i >= lines.size() || lines[i].rfind("%FORMAT", 0) != 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER section %s has no %%FORMAT line\n",
                name.c_str());
        }
        const Amber_Format format = Amber_Parse_Format(
            lines[i], controller, path, name.c_str(), flag_line, i + 1);
        i++;

        Amber_Section section;
        section.name = name;
        section.flag_line = flag_line;
        while (i < lines.size() && lines[i].rfind("%FLAG", 0) != 0)
        {
            if (lines[i].rfind("%COMMENT", 0) == 0 || lines[i].empty())
            {
                i++;
                continue;
            }
            if (!lines[i].empty() && lines[i][0] == '%')
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tunexpected AMBER directive in section %s on "
                    "line %zu\n",
                    name.c_str(), i + 1);
            }
            std::size_t offset = 0;
            const bool format_completed = Amber_For_Each_Format_Field(
                format, [&](const Amber_Format_Field& field) {
                    if (offset >= lines[i].size()) return false;
                    const std::size_t available = lines[i].size() - offset;
                    const std::size_t field_characters =
                        std::min(field.width, available);
                    std::string raw =
                        lines[i].substr(offset, field_characters);
                    std::string value = Amber_Trim(raw);
                    if (!value.empty())
                    {
                        section.values.push_back(value);
                        if (name == "AMBER_ATOM_TYPE")
                        {
                            // Missing bytes at the end of a physical Fortran
                            // record are blank padding.  Keeping only the
                            // materialized prefix avoids allocating a
                            // user-declared width and preserves exact EP
                            // classification because the omitted suffix is
                            // entirely spaces.
                            section.raw_values.push_back(raw);
                        }
                    }
                    if (field.width >= available)
                    {
                        offset = lines[i].size();
                        return false;
                    }
                    offset += field.width;
                    return true;
                });
            if (format_completed && offset < lines[i].size() &&
                !Amber_Trim(lines[i].substr(offset)).empty())
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tAMBER section %s has data beyond its "
                    "%%FORMAT width on line %zu\n",
                    name.c_str(), i + 1);
            }
            i++;
        }
        document.sections.emplace(name, std::move(section));
    }
    return document;
}

static const Amber_Section* Amber_Find_Section(
    const Amber_Parm7_Document& document, const std::string& name)
{
    auto iter = document.sections.find(name);
    return iter == document.sections.end() ? nullptr : &iter->second;
}

static const Amber_Section& Amber_Require_Section(
    const Amber_Parm7_Document& document, const char* name,
    CONTROLLER* controller, const char* error_by)
{
    const Amber_Section* section = Amber_Find_Section(document, name);
    if (section == nullptr)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\trequired AMBER section %s is missing\n", name);
    }
    return *section;
}

static std::vector<int> Amber_Parse_Int_Section(const Amber_Section& section,
                                                CONTROLLER* controller,
                                                const char* error_by)
{
    std::vector<int> parsed;
    parsed.reserve(section.values.size());
    for (const std::string& value : section.values)
    {
        parsed.push_back(
            Amber_Parse_Int(value, section.name.c_str(), controller, error_by));
    }
    return parsed;
}

static std::vector<float> Amber_Parse_Float_Section(
    const Amber_Section& section, CONTROLLER* controller, const char* error_by,
    const char* source_context = "")
{
    std::vector<float> parsed;
    parsed.reserve(section.values.size());
    for (const std::string& value : section.values)
    {
        parsed.push_back(Amber_Parse_Float(
            value, section.name.c_str(), controller, error_by, source_context));
    }
    return parsed;
}

static std::vector<double> Amber_Parse_Double_Section(
    const Amber_Section& section, CONTROLLER* controller, const char* error_by)
{
    std::vector<double> parsed;
    parsed.reserve(section.values.size());
    for (const std::string& value : section.values)
    {
        parsed.push_back(Amber_Parse_Double(value, section.name.c_str(),
                                            controller, error_by));
    }
    return parsed;
}

template <typename T>
static void Amber_Require_Exact_Section_Size(const std::vector<T>& values,
                                             std::size_t expected_count,
                                             CONTROLLER* controller,
                                             const char* error_by,
                                             const char* section_name)
{
    if (values.size() != expected_count)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s has %zu values; expected %zu\n",
            section_name, values.size(), expected_count);
    }
}

static std::vector<int> Amber_Parse_Pointers(
    const Amber_Parm7_Document& document, CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Read_Parm7";
    const Amber_Section* section_pointer =
        Amber_Find_Section(document, "POINTERS");
    if (section_pointer == nullptr)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe AMBER POINTERS section is required\n");
    }
    const Amber_Section& section = *section_pointer;
    if (section.values.size() != 31 && section.values.size() != 32)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section POINTERS has %zu values; expected 31 "
            "or 32\n\tInput file: %s\n\t%%FLAG POINTERS begins on line "
            "%zu\n",
            section.values.size(), document.path.c_str(), section.flag_line);
    }
    std::vector<int> pointers =
        Amber_Parse_Int_Section(section, controller, error_by);
    for (std::size_t i = 0; i < pointers.size(); i++)
    {
        if (pointers[i] < 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER POINTERS entry %zu cannot be negative\n",
                i + 1);
        }
    }
    if (pointers[27] > 2)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tunsupported AMBER IFBOX value %d\n\tInput file: "
            "%s\n\t%%FLAG POINTERS begins on line %zu\n",
            pointers[27], document.path.c_str(), section.flag_line);
    }
    // The optional 32nd pointer is NCOPY, the number of PIMD/LES replicas.
    // NCOPY=1 is a single replica and is dynamically equivalent to an
    // ordinary topology; legacy writers also use 0.  Multiple replicas need
    // topology and dynamics support that the direct reader does not provide.
    if (pointers.size() == 32 && pointers[31] > 1)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER NCOPY=%d requests multiple PIMD/LES replicas, "
            "which are not supported by the direct reader\n"
            "\tInput file: %s\n\t%%FLAG POINTERS begins on line %zu\n",
            pointers[31], document.path.c_str(), section.flag_line);
    }

    const int unsupported_indices[] = {8, 9, 20, 21, 22, 23, 24, 25, 26, 29};
    const char* unsupported_names[] = {"NHPARM", "NPARM", "IFPERT", "NBPER",
                                       "NGPER",  "NDPER", "MBPER",  "MGPER",
                                       "MDPER",  "IFCAP"};
    for (std::size_t i = 0;
         i < sizeof(unsupported_indices) / sizeof(unsupported_indices[0]); i++)
    {
        int index = unsupported_indices[i];
        if (pointers[index] != 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER %s=%d is not supported by the direct "
                "reader\n\tInput file: %s\n\t%%FLAG POINTERS begins on line "
                "%zu\n",
                unsupported_names[i], pointers[index], document.path.c_str(),
                section.flag_line);
        }
    }

    const Amber_Section* ipol = Amber_Find_Section(document, "IPOL");
    if (ipol != nullptr)
    {
        Amber_Require_Exact_Section_Size(ipol->values, 1, controller, error_by,
                                         "IPOL");
        int value =
            Amber_Parse_Int(ipol->values[0], "IPOL", controller, error_by);
        if (value != 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tpolarizable AMBER topology IPOL=%d is not "
                "supported by the direct reader\n\tInput file: %s\n\t%%FLAG "
                "IPOL begins on line %zu\n",
                value, document.path.c_str(), ipol->flag_line);
        }
    }
    return pointers;
}

static std::vector<int> Amber_Build_Canonical_Nonbonded_Map(
    int atom_type_numbers, const std::vector<int>& nonbonded_parm_index,
    const std::vector<float>& hbond_pair_A,
    const std::vector<float>& hbond_pair_B, const char* input_path,
    const Amber_Section* nonbonded_index_section, CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Load_Classical_Force_Field";
    bool has_nonbonded_parm_index = nonbonded_index_section != nullptr;
    if (atom_type_numbers < 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe AMBER atom type count cannot be negative\n");
    }

    std::size_t type_numbers = static_cast<std::size_t>(atom_type_numbers);
    if (type_numbers > 0 &&
        type_numbers >
            std::numeric_limits<std::size_t>::max() / (type_numbers + 1))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER atom type pair table size overflows\n");
    }
    std::size_t pair_type_numbers = type_numbers * (type_numbers + 1) / 2;
    std::vector<int> source_index(pair_type_numbers, -1);
    if (!has_nonbonded_parm_index)
    {
        if (atom_type_numbers > 1)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tNONBONDED_PARM_INDEX is required when AMBER "
                "NTYPES is greater than 1\n");
        }
        // With zero or one type the mapping is provably unambiguous, preserving
        // compatibility with minimal single-type fixtures and legacy files.
        for (std::size_t i = 0; i < pair_type_numbers; i++)
        {
            source_index[i] = static_cast<int>(i);
        }
        return source_index;
    }

    if (type_numbers > 0 &&
        type_numbers > std::numeric_limits<std::size_t>::max() / type_numbers)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tNONBONDED_PARM_INDEX matrix size overflows\n");
    }
    Amber_Require_Exact_Section_Size(nonbonded_parm_index,
                                     type_numbers * type_numbers, controller,
                                     error_by, "NONBONDED_PARM_INDEX");
    for (int high = 0; high < atom_type_numbers; high++)
    {
        for (int low = 0; low <= high; low++)
        {
            int forward = nonbonded_parm_index[static_cast<std::size_t>(low) *
                                                   type_numbers +
                                               high];
            int reverse = nonbonded_parm_index[static_cast<std::size_t>(high) *
                                                   type_numbers +
                                               low];
            if (forward != reverse)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tNONBONDED_PARM_INDEX is asymmetric for atom "
                    "types %d and %d\n",
                    low + 1, high + 1);
            }

            std::size_t canonical_index =
                static_cast<std::size_t>(high) * (high + 1) / 2 + low;
            if (forward > 0)
            {
                if (static_cast<std::size_t>(forward) > pair_type_numbers)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, error_by,
                        "Reason:\n\tNONBONDED_PARM_INDEX value %d for atom "
                        "types %d and %d is outside [1, %zu]\n",
                        forward, low + 1, high + 1, pair_type_numbers);
                }
                source_index[canonical_index] = forward - 1;
                continue;
            }
            if (forward == 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tNONBONDED_PARM_INDEX contains zero for atom "
                    "types %d and %d\n",
                    low + 1, high + 1);
            }

            long long hbond_index_ll = -static_cast<long long>(forward) - 1;
            if (hbond_index_ll < 0 ||
                static_cast<std::size_t>(hbond_index_ll) >=
                    hbond_pair_A.size() ||
                static_cast<std::size_t>(hbond_index_ll) >= hbond_pair_B.size())
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tNONBONDED_PARM_INDEX value %d for atom "
                    "types %d and %d does not select a valid HBOND parameter\n",
                    forward, low + 1, high + 1);
            }
            std::size_t hbond_index = static_cast<std::size_t>(hbond_index_ll);
            if (hbond_pair_A[hbond_index] != 0.0f ||
                hbond_pair_B[hbond_index] != 0.0f)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tAMBER 10-12 nonbonded terms are not supported; "
                    "atom types %d and %d select nonzero HBOND parameter %zu"
                    "\n\tInput file: %s\n\t%%FLAG NONBONDED_PARM_INDEX "
                    "begins on line %zu\n",
                    low + 1, high + 1, hbond_index + 1, input_path,
                    nonbonded_index_section->flag_line);
            }
            // Negative indices with an all-zero HBOND entry are legacy
            // placeholders (commonly used by old water models), not a 12-6 LJ
            // interaction. Leave this canonical pair at zero.
            source_index[canonical_index] = -1;
        }
    }
    return source_index;
}

static void Amber_Remap_LJ_Pair_Matrix(
    const std::vector<float>& raw_pair_A, const std::vector<float>& raw_pair_B,
    bool has_pair_A, bool has_pair_B, const std::vector<int>& source_index,
    std::vector<float>* canonical_pair_A, std::vector<float>* canonical_pair_B,
    CONTROLLER* controller, const char* section_A, const char* section_B)
{
    const char* error_by = "Xponge::Amber_Load_Classical_Force_Field";
    if (has_pair_A != has_pair_B)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER sections %s and %s must either both be present "
            "or both be absent\n",
            section_A, section_B);
    }
    canonical_pair_A->clear();
    canonical_pair_B->clear();
    if (!has_pair_A)
    {
        return;
    }

    Amber_Require_Exact_Section_Size(raw_pair_A, source_index.size(),
                                     controller, error_by, section_A);
    Amber_Require_Exact_Section_Size(raw_pair_B, source_index.size(),
                                     controller, error_by, section_B);
    canonical_pair_A->assign(source_index.size(), 0.0f);
    canonical_pair_B->assign(source_index.size(), 0.0f);
    for (std::size_t canonical_index = 0; canonical_index < source_index.size();
         canonical_index++)
    {
        int source = source_index[canonical_index];
        if (source >= 0)
        {
            canonical_pair_A->at(canonical_index) = Amber_Checked_Scale(
                raw_pair_A.at(static_cast<std::size_t>(source)), 12.0,
                section_A, controller, error_by);
            canonical_pair_B->at(canonical_index) = Amber_Checked_Scale(
                raw_pair_B.at(static_cast<std::size_t>(source)), 6.0, section_B,
                controller, error_by);
        }
    }
}

static int Amber_Get_Atom_Numbers(const System* system);
static void Amber_Ensure_Atom_Numbers(System* system, int atom_numbers,
                                      CONTROLLER* controller,
                                      const char* error_by);

static void Amber_Load_Parm7(System* system,
                             const Amber_Parm7_Document& document,
                             const std::vector<int>& pointers,
                             CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Load_Parm7";
    int atom_numbers = pointers[0];
    int excluded_total = pointers[10];
    int residue_numbers = pointers[11];
    Amber_Ensure_Atom_Numbers(system, atom_numbers, controller, error_by);
    if (atom_numbers > 0 && residue_numbers <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER NRES must be positive when NATOM is positive\n");
    }

    const Amber_Section& mass_section =
        Amber_Require_Section(document, "MASS", controller, error_by);
    const Amber_Section& charge_section =
        Amber_Require_Section(document, "CHARGE", controller, error_by);
    const Amber_Section& residue_section = Amber_Require_Section(
        document, "RESIDUE_POINTER", controller, error_by);
    const Amber_Section& excluded_number_section = Amber_Require_Section(
        document, "NUMBER_EXCLUDED_ATOMS", controller, error_by);
    const Amber_Section& excluded_list_section = Amber_Require_Section(
        document, "EXCLUDED_ATOMS_LIST", controller, error_by);
    Amber_Require_Exact_Section_Size(mass_section.values,
                                     static_cast<std::size_t>(atom_numbers),
                                     controller, error_by, "MASS");
    Amber_Require_Exact_Section_Size(charge_section.values,
                                     static_cast<std::size_t>(atom_numbers),
                                     controller, error_by, "CHARGE");
    Amber_Require_Exact_Section_Size(residue_section.values,
                                     static_cast<std::size_t>(residue_numbers),
                                     controller, error_by, "RESIDUE_POINTER");
    Amber_Require_Exact_Section_Size(
        excluded_number_section.values, static_cast<std::size_t>(atom_numbers),
        controller, error_by, "NUMBER_EXCLUDED_ATOMS");
    Amber_Require_Exact_Section_Size(
        excluded_list_section.values, static_cast<std::size_t>(excluded_total),
        controller, error_by, "EXCLUDED_ATOMS_LIST");

    system->atoms.mass =
        Amber_Parse_Float_Section(mass_section, controller, error_by);
    system->atoms.charge =
        Amber_Parse_Float_Section(charge_section, controller, error_by);
    for (int i = 0; i < atom_numbers; i++)
    {
        if (system->atoms.mass[i] < 0.0f)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER MASS for atom %d cannot be negative\n",
                i + 1);
        }
    }

    std::vector<int> residue_pointer =
        Amber_Parse_Int_Section(residue_section, controller, error_by);
    system->residues.atom_numbers.clear();
    system->residues.atom_numbers.reserve(residue_numbers);
    for (int i = 0; i < residue_numbers; i++)
    {
        int pointer = residue_pointer[i];
        if (pointer < 1 || pointer > atom_numbers || (i == 0 && pointer != 1) ||
            (i > 0 && pointer <= residue_pointer[i - 1]))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tRESIDUE_POINTER entry %d is not a strictly "
                "increasing atom index beginning at 1\n",
                i + 1);
        }
        int next =
            i + 1 < residue_numbers ? residue_pointer[i + 1] : atom_numbers + 1;
        system->residues.atom_numbers.push_back(next - pointer);
    }

    const Amber_Section* radii_section = Amber_Find_Section(document, "RADII");
    system->generalized_born.radius.clear();
    if (radii_section != nullptr)
    {
        Amber_Require_Exact_Section_Size(radii_section->values,
                                         static_cast<std::size_t>(atom_numbers),
                                         controller, error_by, "RADII");
        system->generalized_born.radius =
            Amber_Parse_Float_Section(*radii_section, controller, error_by);
    }
    const Amber_Section* screen_section =
        Amber_Find_Section(document, "SCREEN");
    system->generalized_born.scale_factor.clear();
    if (screen_section != nullptr)
    {
        Amber_Require_Exact_Section_Size(screen_section->values,
                                         static_cast<std::size_t>(atom_numbers),
                                         controller, error_by, "SCREEN");
        system->generalized_born.scale_factor =
            Amber_Parse_Float_Section(*screen_section, controller, error_by);
    }

    std::vector<int> excluded_numbers =
        Amber_Parse_Int_Section(excluded_number_section, controller, error_by);
    std::vector<int> excluded_list =
        Amber_Parse_Int_Section(excluded_list_section, controller, error_by);
    const std::string exclusion_context = Amber_Source_Context(
        document, {"NUMBER_EXCLUDED_ATOMS", "EXCLUDED_ATOMS_LIST"});
    long long count_sum = 0;
    for (int i = 0; i < atom_numbers; i++)
    {
        if (excluded_numbers[i] < 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tNUMBER_EXCLUDED_ATOMS entry %d cannot be "
                "negative\n%s",
                i + 1, exclusion_context.c_str());
        }
        count_sum += excluded_numbers[i];
        if (count_sum > std::numeric_limits<int>::max())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER exclusion count sum overflows\n%s",
                exclusion_context.c_str());
        }
    }
    if (count_sum != excluded_total)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tNUMBER_EXCLUDED_ATOMS sums to %lld; POINTERS NNB is "
            "%d\n%s",
            count_sum, excluded_total, exclusion_context.c_str());
    }

    system->exclusions.excluded_atoms.assign(atom_numbers, {});
    std::size_t cursor = 0;
    for (int atom = 0; atom < atom_numbers; atom++)
    {
        int count = excluded_numbers[atom];
        if (count == 1 && excluded_list[cursor] == 0)
        {
            cursor++;
            continue;
        }
        int previous = 0;
        for (int j = 0; j < count; j++)
        {
            int excluded = excluded_list[cursor++];
            if (excluded <= 0 || excluded > atom_numbers ||
                excluded == atom + 1)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tEXCLUDED_ATOMS_LIST contains invalid atom %d "
                    "for source atom %d\n%s",
                    excluded, atom + 1, exclusion_context.c_str());
            }
            if (excluded < atom + 1)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tEXCLUDED_ATOMS_LIST atom %d for source atom "
                    "%d violates Amber triangular ordering; every excluded "
                    "atom must be greater than its source atom\n%s",
                    excluded, atom + 1, exclusion_context.c_str());
            }
            if (j > 0 && excluded <= previous)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tEXCLUDED_ATOMS_LIST for atom %d is not "
                    "strictly increasing\n%s",
                    atom + 1, exclusion_context.c_str());
            }
            previous = excluded;
            system->exclusions.excluded_atoms[atom].push_back(excluded - 1);
        }
    }
}

static std::vector<std::string> Amber_Tokenize_Restart_Header(
    const std::string& line, std::size_t line_number, CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Load_Rst7";
    std::vector<std::string> tokens;
    std::size_t position = 0;
    while (position < line.size())
    {
        while (position < line.size() &&
               isspace(static_cast<unsigned char>(line[position])))
        {
            position++;
        }
        if (position == line.size())
        {
            break;
        }
        std::size_t begin = position;
        if (line[position] == '+' || line[position] == '-')
        {
            position++;
        }
        bool has_digit = false;
        while (position < line.size() &&
               isdigit(static_cast<unsigned char>(line[position])))
        {
            has_digit = true;
            position++;
        }
        if (position < line.size() && line[position] == '.')
        {
            position++;
            while (position < line.size() &&
                   isdigit(static_cast<unsigned char>(line[position])))
            {
                has_digit = true;
                position++;
            }
        }
        if (!has_digit)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tinvalid numeric field in amber_rst7 on line "
                "%zu\n",
                line_number);
        }
        if (position < line.size() &&
            (line[position] == 'E' || line[position] == 'e' ||
             line[position] == 'D' || line[position] == 'd'))
        {
            position++;
            if (position < line.size() &&
                (line[position] == '+' || line[position] == '-'))
            {
                position++;
            }
            std::size_t exponent_begin = position;
            while (position < line.size() &&
                   isdigit(static_cast<unsigned char>(line[position])))
            {
                position++;
            }
            if (position == exponent_begin)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tinvalid exponent in amber_rst7 on line %zu\n",
                    line_number);
            }
        }
        tokens.push_back(line.substr(begin, position - begin));
        if (position < line.size() &&
            !isspace(static_cast<unsigned char>(line[position])) &&
            line[position] != '+' && line[position] != '-')
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\ttrailing characters in amber_rst7 numeric field "
                "on line %zu\n",
                line_number);
        }
    }
    return tokens;
}

static bool Amber_Is_Strict_Restart_Number(const std::string& raw)
{
    std::string value = Amber_Trim(raw);
    for (char& character : value)
    {
        if (character == 'D' || character == 'd')
        {
            character = 'E';
        }
    }
    return Amber_Is_Strict_Float_Lexeme(value);
}

static bool Amber_Try_Parse_Restart_Data_Line(
    const std::string& raw_line, std::size_t expected_fields,
    std::vector<std::string>* parsed_fields)
{
    std::string line = raw_line;
    if (!line.empty() && line.back() == '\r') line.pop_back();

    constexpr std::size_t kFieldWidth = 12;
    std::vector<std::string> fields;
    std::size_t fixed_width = 0;
    if (expected_fields <=
        std::numeric_limits<std::size_t>::max() / kFieldWidth)
    {
        fixed_width = expected_fields * kFieldWidth;
    }
    if (fixed_width != 0 && line.size() >= fixed_width &&
        Amber_Trim(line.substr(fixed_width)).empty())
    {
        bool fixed_width_valid = true;
        fields.reserve(expected_fields);
        for (std::size_t i = 0; i < expected_fields; i++)
        {
            std::string field = line.substr(i * kFieldWidth, kFieldWidth);
            if (!Amber_Is_Strict_Restart_Number(field))
            {
                fixed_width_valid = false;
                break;
            }
            fields.push_back(field);
        }
        if (fixed_width_valid)
        {
            *parsed_fields = std::move(fields);
            return true;
        }
        fields.clear();
    }

    bool looks_like_truncated_f12 = line.size() != fixed_width;
    for (std::size_t i = 0; looks_like_truncated_f12 && i < expected_fields;
         i++)
    {
        std::size_t decimal_position = i * kFieldWidth + 4;
        if (decimal_position >= line.size() || line[decimal_position] != '.')
        {
            looks_like_truncated_f12 = false;
        }
    }
    std::istringstream stream(line);
    std::string field;
    while (!looks_like_truncated_f12 && stream >> field)
    {
        if (!Amber_Is_Strict_Restart_Number(field))
        {
            fields.clear();
            break;
        }
        fields.push_back(field);
    }
    if (fields.size() == expected_fields)
    {
        *parsed_fields = std::move(fields);
        return true;
    }
    return false;
}

static std::vector<std::string> Amber_Parse_Restart_Data_Line(
    const std::string& raw_line, std::size_t expected_fields,
    const char* block_name, std::size_t line_number, CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Load_Rst7";
    std::vector<std::string> fields;
    if (Amber_Try_Parse_Restart_Data_Line(raw_line, expected_fields, &fields))
    {
        return fields;
    }

    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorBadFileFormat, error_by,
        "Reason:\n\tinvalid amber_rst7 %s record on line %zu; expected %zu "
        "fixed-width F12.7 fields or whitespace-delimited scientific fields"
        "\n",
        block_name, line_number, expected_fields);
    return {};
}

static std::vector<float> Amber_Parse_Restart_Block(
    const std::vector<std::string>& lines, std::size_t* line_index,
    std::size_t value_count, const char* block_name, CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Load_Rst7";
    std::vector<float> values;
    values.reserve(value_count);
    while (values.size() < value_count)
    {
        if (*line_index >= lines.size())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tamber_rst7 %s block is truncated\n", block_name);
        }
        std::size_t expected_fields =
            std::min<std::size_t>(6, value_count - values.size());
        std::vector<std::string> fields = Amber_Parse_Restart_Data_Line(
            lines[*line_index], expected_fields, block_name, *line_index + 1,
            controller);
        for (const std::string& field : fields)
        {
            values.push_back(
                Amber_Parse_Float(field, "amber_rst7", controller, error_by));
        }
        (*line_index)++;
    }
    return values;
}

static int Amber_Get_Atom_Numbers(const System* system)
{
    return Load_Get_Atom_Numbers(system);
}

static void Amber_Ensure_Atom_Numbers(System* system, int atom_numbers,
                                      CONTROLLER* controller,
                                      const char* error_by)
{
    int current_atom_numbers = Amber_Get_Atom_Numbers(system);
    if (current_atom_numbers < 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe retained atom arrays have inconsistent, "
            "misaligned, or unsupported sizes\n");
    }
    if (current_atom_numbers > 0 && current_atom_numbers != atom_numbers)
    {
        controller->Throw_SPONGE_Error(spongeErrorConflictingCommand, error_by,
                                       "Reason:\n\t'atom_numbers' is different "
                                       "in different input files\n");
    }
}

static std::size_t Amber_Checked_Product(std::size_t left, std::size_t right,
                                         const char* description,
                                         CONTROLLER* controller,
                                         const char* error_by)
{
    if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER %s size overflows\n", description);
    }
    return left * right;
}

static int Amber_Decode_Encoded_Atom(int raw_index, bool allow_sign,
                                     int atom_numbers,
                                     const char* interaction_name,
                                     const char* field_name,
                                     CONTROLLER* controller,
                                     const char* error_by,
                                     const char* source_context = "")
{
    if (!allow_sign && raw_index < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER %s field %s cannot be negative\n%s",
            interaction_name, field_name, source_context);
    }
    long long magnitude = static_cast<long long>(raw_index);
    if (magnitude < 0)
    {
        magnitude = -magnitude;
    }
    if (magnitude % 3 != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER %s field %s value %d is not a multiple of "
            "3\n%s",
            interaction_name, field_name, raw_index, source_context);
    }
    long long atom_index = magnitude / 3;
    if (atom_index < 0 || atom_index >= atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER %s field %s selects atom %lld outside [0, "
            "%d)\n%s",
            interaction_name, field_name, atom_index, atom_numbers,
            source_context);
    }
    return static_cast<int>(atom_index);
}

static int Amber_Decode_One_Based_Index(int raw_index, int upper_bound,
                                        const char* section_name,
                                        const char* index_kind,
                                        CONTROLLER* controller,
                                        const char* error_by)
{
    if (raw_index <= 0 || raw_index > upper_bound)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s contains an out-of-range %s index "
            "%d; expected [1, %d]\n",
            section_name, index_kind, raw_index, upper_bound);
    }
    return raw_index - 1;
}

static std::vector<float> Amber_Load_Optional_Float_Table(
    const Amber_Parm7_Document& document, const char* section_name,
    std::size_t expected_count, bool required, CONTROLLER* controller,
    const char* error_by)
{
    const Amber_Section* section = Amber_Find_Section(document, section_name);
    if (section == nullptr)
    {
        if (required)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\trequired AMBER section %s is missing\n",
                section_name);
        }
        return {};
    }
    Amber_Require_Exact_Section_Size(section->values, expected_count,
                                     controller, error_by, section_name);
    const std::string source_context =
        Amber_Source_Context(document, {section_name});
    return Amber_Parse_Float_Section(*section, controller, error_by,
                                     source_context.c_str());
}

static std::vector<double> Amber_Load_Optional_Double_Table(
    const Amber_Parm7_Document& document, const char* section_name,
    std::size_t expected_count, bool required, CONTROLLER* controller,
    const char* error_by)
{
    const Amber_Section* section = Amber_Find_Section(document, section_name);
    if (section == nullptr)
    {
        if (required)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\trequired AMBER section %s is missing\n",
                section_name);
        }
        return {};
    }
    Amber_Require_Exact_Section_Size(section->values, expected_count,
                                     controller, error_by, section_name);
    return Amber_Parse_Double_Section(*section, controller, error_by);
}

static std::vector<int> Amber_Load_Interaction_Tuples(
    const Amber_Parm7_Document& document, const char* section_name,
    int expected_count, std::size_t tuple_width, CONTROLLER* controller,
    const char* error_by)
{
    if (expected_count < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER pointer count for %s cannot be negative\n",
            section_name);
    }
    const Amber_Section* section = Amber_Find_Section(document, section_name);
    if (section == nullptr)
    {
        if (expected_count != 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER section %s is missing; POINTERS declares "
                "%d tuples\n",
                section_name, expected_count);
        }
        return {};
    }
    std::size_t expected_values =
        Amber_Checked_Product(static_cast<std::size_t>(expected_count),
                              tuple_width, section_name, controller, error_by);
    if (section->values.size() != expected_values)
    {
        if (section->values.size() % tuple_width == 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER section %s has %zu tuples; POINTERS "
                "declares %d\n",
                section_name, section->values.size() / tuple_width,
                expected_count);
        }
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER section %s has %zu values; expected a multiple "
            "of %zu and exactly %zu values\n",
            section_name, section->values.size(), tuple_width, expected_values);
    }
    return Amber_Parse_Int_Section(*section, controller, error_by);
}

struct Amber_Angle_Record
{
    int atom_a = 0;
    int atom_b = 0;
    int atom_c = 0;
    float k = 0.0f;
    float theta0 = 0.0f;
    bool moved_to_urey_bradley = false;
};

static void Amber_Append_Angle(const Amber_Angle_Record& angle, Angles* output)
{
    output->atom_a.push_back(angle.atom_a);
    output->atom_b.push_back(angle.atom_b);
    output->atom_c.push_back(angle.atom_c);
    output->k.push_back(angle.k);
    output->theta0.push_back(angle.theta0);
}

static bool Amber_Is_Exact_EP_Field(const std::string& raw_field)
{
    if (raw_field.size() < 2 || raw_field[0] != 'E' || raw_field[1] != 'P')
    {
        return false;
    }
    for (std::size_t i = 2; i < raw_field.size(); i++)
    {
        if (raw_field[i] != ' ')
        {
            return false;
        }
    }
    return true;
}

static void Amber_Load_Classical_Force_Field(
    System* system, const Amber_Parm7_Document& document,
    const std::vector<int>& pointers, CONTROLLER* controller)
{
    const char* error_by = "Xponge::Amber_Load_Classical_Force_Field";
    ClassicalForceField* ff = &system->classical_force_field;
    ff->bonds = Bonds{};
    ff->constraints = DistanceConstraints{};
    ff->angles = Angles{};
    ff->dihedrals = Torsions{};
    ff->impropers = Torsions{};
    ff->nb14 = NB14{};
    ff->lj = LennardJones{};
    ff->cmap = CMap{};
    ff->urey_bradley = UreyBradley{};

    int atom_numbers = pointers[0];
    int atom_type_numbers = pointers[1];
    int bond_with_hydrogen = pointers[2];
    int bond_without_hydrogen = pointers[3];  // MBONA, not NBONA
    int angle_with_hydrogen = pointers[4];
    int angle_without_hydrogen = pointers[5];  // MTHETA, not NTHETA
    int dihedral_with_hydrogen = pointers[6];
    int dihedral_without_hydrogen = pointers[7];  // MPHIA, not NPHIA
    int bond_type_numbers = pointers[15];
    int angle_type_numbers = pointers[16];
    int dihedral_type_numbers = pointers[17];
    int hbond_type_numbers = pointers[19];
    int extra_point_numbers = pointers[30];
    std::vector<bool> is_virtual_atom(static_cast<std::size_t>(atom_numbers),
                                      false);

    Amber_Ensure_Atom_Numbers(system, atom_numbers, controller, error_by);
    if (atom_numbers > 0 && atom_type_numbers <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER NTYPES must be positive when NATOM is "
            "positive\n");
    }

    bool has_bond_k =
        Amber_Find_Section(document, "BOND_FORCE_CONSTANT") != nullptr;
    bool has_bond_r0 =
        Amber_Find_Section(document, "BOND_EQUIL_VALUE") != nullptr;
    if (has_bond_k != has_bond_r0 || (bond_type_numbers > 0 && !has_bond_k))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tBOND_FORCE_CONSTANT and BOND_EQUIL_VALUE must both "
            "be present when AMBER bond parameter types are declared\n");
    }
    std::vector<float> bond_type_k = Amber_Load_Optional_Float_Table(
        document, "BOND_FORCE_CONSTANT",
        static_cast<std::size_t>(bond_type_numbers), false, controller,
        error_by);
    std::vector<float> bond_type_r0 = Amber_Load_Optional_Float_Table(
        document, "BOND_EQUIL_VALUE",
        static_cast<std::size_t>(bond_type_numbers), false, controller,
        error_by);
    for (std::size_t i = 0; i < bond_type_k.size(); i++)
    {
        if (bond_type_k[i] < 0.0f || bond_type_r0[i] < 0.0f)
        {
            std::string source_context = Amber_Source_Context(
                document, {"BOND_FORCE_CONSTANT", "BOND_EQUIL_VALUE"});
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER bond parameter type %zu has negative K or "
                "equilibrium distance\n%s",
                i + 1, source_context.c_str());
        }
    }

    bool has_angle_k =
        Amber_Find_Section(document, "ANGLE_FORCE_CONSTANT") != nullptr;
    bool has_angle_theta =
        Amber_Find_Section(document, "ANGLE_EQUIL_VALUE") != nullptr;
    if (has_angle_k != has_angle_theta ||
        (angle_type_numbers > 0 && !has_angle_k))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tANGLE_FORCE_CONSTANT and ANGLE_EQUIL_VALUE must both "
            "be present when AMBER angle parameter types are declared\n");
    }
    std::vector<float> angle_type_k = Amber_Load_Optional_Float_Table(
        document, "ANGLE_FORCE_CONSTANT",
        static_cast<std::size_t>(angle_type_numbers), false, controller,
        error_by);
    std::vector<float> angle_type_theta0 = Amber_Load_Optional_Float_Table(
        document, "ANGLE_EQUIL_VALUE",
        static_cast<std::size_t>(angle_type_numbers), false, controller,
        error_by);
    for (std::size_t i = 0; i < angle_type_k.size(); i++)
    {
        if (angle_type_k[i] < 0.0f || angle_type_theta0[i] < 0.0f ||
            angle_type_theta0[i] > CONSTANT_Pi)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER angle parameter type %zu is outside its "
                "physical domain\n",
                i + 1);
        }
    }

    std::vector<int> raw_bonds_h = Amber_Load_Interaction_Tuples(
        document, "BONDS_INC_HYDROGEN", bond_with_hydrogen, 3, controller,
        error_by);
    std::vector<int> raw_bonds_no_h = Amber_Load_Interaction_Tuples(
        document, "BONDS_WITHOUT_HYDROGEN", bond_without_hydrogen, 3,
        controller, error_by);
    std::vector<int> raw_bonds = raw_bonds_h;
    raw_bonds.insert(raw_bonds.end(), raw_bonds_no_h.begin(),
                     raw_bonds_no_h.end());
    std::vector<bool> bond_includes_hydrogen;
    bond_includes_hydrogen.reserve(raw_bonds.size() / 3);
    for (std::size_t i = 0; i < raw_bonds.size(); i += 3)
    {
        int atom_a =
            Amber_Decode_Encoded_Atom(raw_bonds[i], false, atom_numbers, "bond",
                                      "A", controller, error_by);
        int atom_b =
            Amber_Decode_Encoded_Atom(raw_bonds[i + 1], false, atom_numbers,
                                      "bond", "B", controller, error_by);
        int type_index = Amber_Decode_One_Based_Index(
            raw_bonds[i + 2], bond_type_numbers, "BONDS", "parameter",
            controller, error_by);
        if (atom_a == atom_b)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER bond cannot connect an atom to itself\n");
        }
        ff->bonds.atom_a.push_back(atom_a);
        ff->bonds.atom_b.push_back(atom_b);
        ff->bonds.k.push_back(bond_type_k[type_index]);
        ff->bonds.r0.push_back(bond_type_r0[type_index]);
        bond_includes_hydrogen.push_back(i < raw_bonds_h.size());
    }

    std::vector<int> raw_angles_h = Amber_Load_Interaction_Tuples(
        document, "ANGLES_INC_HYDROGEN", angle_with_hydrogen, 4, controller,
        error_by);
    std::vector<int> raw_angles_no_h = Amber_Load_Interaction_Tuples(
        document, "ANGLES_WITHOUT_HYDROGEN", angle_without_hydrogen, 4,
        controller, error_by);
    std::vector<int> raw_angles = raw_angles_h;
    raw_angles.insert(raw_angles.end(), raw_angles_no_h.begin(),
                      raw_angles_no_h.end());
    std::vector<Amber_Angle_Record> angle_records;
    angle_records.reserve(raw_angles.size() / 4);
    for (std::size_t i = 0; i < raw_angles.size(); i += 4)
    {
        Amber_Angle_Record angle;
        angle.atom_a =
            Amber_Decode_Encoded_Atom(raw_angles[i], false, atom_numbers,
                                      "angle", "A", controller, error_by);
        angle.atom_b =
            Amber_Decode_Encoded_Atom(raw_angles[i + 1], false, atom_numbers,
                                      "angle", "B", controller, error_by);
        angle.atom_c =
            Amber_Decode_Encoded_Atom(raw_angles[i + 2], false, atom_numbers,
                                      "angle", "C", controller, error_by);
        int type_index = Amber_Decode_One_Based_Index(
            raw_angles[i + 3], angle_type_numbers, "ANGLES", "parameter",
            controller, error_by);
        if (angle.atom_a == angle.atom_b || angle.atom_b == angle.atom_c ||
            angle.atom_a == angle.atom_c)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER angle contains repeated atoms\n");
        }
        angle.k = angle_type_k[type_index];
        angle.theta0 = angle_type_theta0[type_index];
        angle_records.push_back(angle);
    }

    system->virtual_atoms.records.clear();
    const Amber_Section* amber_atom_type_section =
        Amber_Find_Section(document, "AMBER_ATOM_TYPE");
    if (amber_atom_type_section != nullptr)
    {
        Amber_Require_Exact_Section_Size(amber_atom_type_section->values,
                                         static_cast<std::size_t>(atom_numbers),
                                         controller, error_by,
                                         "AMBER_ATOM_TYPE");
    }
    if (extra_point_numbers > 0 && amber_atom_type_section == nullptr)
    {
        std::string source_context =
            Amber_Source_Context(document, {"POINTERS", "AMBER_ATOM_TYPE"});
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER POINTERS NUMEXTRA=%d requires "
            "AMBER_ATOM_TYPE to identify Classic extra points\n%s",
            extra_point_numbers, source_context.c_str());
    }

    const Amber_Section* atomic_number_section =
        Amber_Find_Section(document, "ATOMIC_NUMBER");
    std::vector<int> atomic_numbers;
    if (atomic_number_section != nullptr)
    {
        Amber_Require_Exact_Section_Size(atomic_number_section->values,
                                         static_cast<std::size_t>(atom_numbers),
                                         controller, error_by, "ATOMIC_NUMBER");
        atomic_numbers = Amber_Parse_Int_Section(*atomic_number_section,
                                                 controller, error_by);
        for (int atom = 0; atom < atom_numbers; atom++)
        {
            if (atomic_numbers[atom] < 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tATOMIC_NUMBER for atom %d cannot be "
                    "negative\n",
                    atom + 1);
            }
        }
    }

    struct Amber_Bond_Edge
    {
        int other;
        float r0;
        std::size_t bond_index;
        bool includes_hydrogen;
    };
    std::vector<std::vector<Amber_Bond_Edge>> bond_adjacency(
        static_cast<std::size_t>(atom_numbers));
    for (std::size_t i = 0; i < ff->bonds.atom_a.size(); i++)
    {
        int atom_a = ff->bonds.atom_a[i];
        int atom_b = ff->bonds.atom_b[i];
        bond_adjacency[atom_a].push_back(
            {atom_b, ff->bonds.r0[i], i, bond_includes_hydrogen[i]});
        bond_adjacency[atom_b].push_back(
            {atom_a, ff->bonds.r0[i], i, bond_includes_hydrogen[i]});
    }

    int found_extra_points = 0;
    if (amber_atom_type_section != nullptr)
    {
        for (const std::string& raw_atom_type :
             amber_atom_type_section->raw_values)
        {
            found_extra_points += Amber_Is_Exact_EP_Field(raw_atom_type);
        }
    }
    if (found_extra_points != extra_point_numbers)
    {
        std::string source_context =
            Amber_Source_Context(document, {"POINTERS", "AMBER_ATOM_TYPE"});
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER POINTERS NUMEXTRA=%d but AMBER_ATOM_TYPE "
            "identifies %d exact EP atom(s)\n%s",
            extra_point_numbers, found_extra_points, source_context.c_str());
    }

    const std::string ep_edge_context = Amber_Source_Context(
        document,
        {"AMBER_ATOM_TYPE", "BONDS_WITHOUT_HYDROGEN", "BONDS_INC_HYDROGEN"});
    const std::string ep_parameter_context = Amber_Source_Context(
        document,
        {"AMBER_ATOM_TYPE", "BONDS_WITHOUT_HYDROGEN", "BOND_EQUIL_VALUE"});
    for (int atom = 0; atom < atom_numbers; atom++)
    {
        if (amber_atom_type_section == nullptr ||
            !Amber_Is_Exact_EP_Field(amber_atom_type_section->raw_values[atom]))
        {
            continue;
        }
        is_virtual_atom[atom] = true;

        const std::vector<Amber_Bond_Edge>& ep_edges = bond_adjacency[atom];
        if (ep_edges.size() != 1 || ep_edges[0].includes_hydrogen)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER NUMEXTRA atom %d must have exactly one "
                "parent edge in BONDS_WITHOUT_HYDROGEN and none in "
                "BONDS_INC_HYDROGEN; found %zu total edge(s)\n"
                "%s",
                atom + 1, ep_edges.size(), ep_edge_context.c_str());
        }
        int parent = ep_edges[0].other;
        const std::vector<Amber_Bond_Edge>& parent_edges =
            bond_adjacency[parent];
        std::vector<Amber_Bond_Edge> hydrogen_edges;
        int owner_without_hydrogen_edges = 0;
        for (const Amber_Bond_Edge& edge : parent_edges)
        {
            if (edge.includes_hydrogen)
            {
                hydrogen_edges.push_back(edge);
            }
            else
            {
                owner_without_hydrogen_edges++;
            }
        }
        if (owner_without_hydrogen_edges != 1 || hydrogen_edges.size() != 2)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER NUMEXTRA atom %d is not a Classic TIP4P "
                "frame: its owner must have exactly one "
                "BONDS_WITHOUT_HYDROGEN edge (the EP) and two "
                "BONDS_INC_HYDROGEN edges; found %d and %zu\n"
                "%s",
                atom + 1, owner_without_hydrogen_edges, hydrogen_edges.size(),
                ep_edge_context.c_str());
        }
        int hydrogen_a = hydrogen_edges[0].other;
        int hydrogen_b = hydrogen_edges[1].other;
        double r_ep = ep_edges[0].r0;
        if (!Amber_Is_Finite_Double(r_ep) || r_ep < 0.0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER NUMEXTRA atom %d O-EP equilibrium distance "
                "must be finite and non-negative\n%s",
                atom + 1, ep_parameter_context.c_str());
        }
        VirtualAtomRecord record;
        record.type = 5;
        record.virtual_atom = atom;
        record.from = {parent, hydrogen_a, hydrogen_b};
        record.parameter = {Amber_Checked_Float(
            r_ep, "AMBER NUMEXTRA TIP4P O-EP distance", controller, error_by)};
        system->virtual_atoms.records.push_back(record);
        // Amber Classic fix_masses canonicalizes every recognized extra
        // point to zero mass after frame assignment, regardless of the MASS
        // value stored in the topology.
        system->atoms.mass[atom] = 0.0f;
    }

    if (!ff->bonds.atom_a.empty())
    {
        Bonds chemical_bonds;
        for (std::size_t i = 0; i < ff->bonds.atom_a.size(); i++)
        {
            if (is_virtual_atom[ff->bonds.atom_a[i]] ||
                is_virtual_atom[ff->bonds.atom_b[i]])
            {
                continue;
            }
            chemical_bonds.atom_a.push_back(ff->bonds.atom_a[i]);
            chemical_bonds.atom_b.push_back(ff->bonds.atom_b[i]);
            chemical_bonds.k.push_back(ff->bonds.k[i]);
            chemical_bonds.r0.push_back(ff->bonds.r0[i]);
        }
        ff->bonds = std::move(chemical_bonds);
    }

    const char* ub_names[] = {"CHARMM_UREY_BRADLEY_COUNT",
                              "CHARMM_UREY_BRADLEY",
                              "CHARMM_UREY_BRADLEY_FORCE_CONSTANT",
                              "CHARMM_UREY_BRADLEY_EQUIL_VALUE"};
    int ub_flag_count = 0;
    for (const char* name : ub_names)
    {
        ub_flag_count += Amber_Find_Section(document, name) != nullptr;
    }
    if (ub_flag_count != 0 && ub_flag_count != 4)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe CHAMBER Urey-Bradley sections are incomplete\n");
    }
    if (ub_flag_count == 4)
    {
        const Amber_Section& count_section = Amber_Require_Section(
            document, "CHARMM_UREY_BRADLEY_COUNT", controller, error_by);
        Amber_Require_Exact_Section_Size(count_section.values, 2, controller,
                                         error_by, "CHARMM_UREY_BRADLEY_COUNT");
        std::vector<int> counts =
            Amber_Parse_Int_Section(count_section, controller, error_by);
        int ub_numbers = counts[0];
        int ub_type_numbers = counts[1];
        if (ub_numbers < 0 || ub_type_numbers < 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tCHAMBER Urey-Bradley counts cannot be negative\n");
        }
        const Amber_Section& tuple_section = Amber_Require_Section(
            document, "CHARMM_UREY_BRADLEY", controller, error_by);
        std::size_t expected_tuple_values =
            Amber_Checked_Product(static_cast<std::size_t>(ub_numbers), 3,
                                  "CHARMM_UREY_BRADLEY", controller, error_by);
        Amber_Require_Exact_Section_Size(tuple_section.values,
                                         expected_tuple_values, controller,
                                         error_by, "CHARMM_UREY_BRADLEY");
        std::vector<int> tuples =
            Amber_Parse_Int_Section(tuple_section, controller, error_by);
        std::vector<float> ub_k = Amber_Load_Optional_Float_Table(
            document, "CHARMM_UREY_BRADLEY_FORCE_CONSTANT",
            static_cast<std::size_t>(ub_type_numbers), true, controller,
            error_by);
        std::vector<float> ub_r0 = Amber_Load_Optional_Float_Table(
            document, "CHARMM_UREY_BRADLEY_EQUIL_VALUE",
            static_cast<std::size_t>(ub_type_numbers), true, controller,
            error_by);
        for (int i = 0; i < ub_type_numbers; i++)
        {
            if (ub_k[i] < 0.0f || ub_r0[i] < 0.0f)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tCHAMBER Urey-Bradley parameter type %d has "
                    "negative K or equilibrium distance\n",
                    i + 1);
            }
        }
        for (int i = 0; i < ub_numbers; i++)
        {
            int atom_a = Amber_Decode_One_Based_Index(
                tuples[3 * i], atom_numbers, "CHARMM_UREY_BRADLEY", "atom",
                controller, error_by);
            int atom_c = Amber_Decode_One_Based_Index(
                tuples[3 * i + 1], atom_numbers, "CHARMM_UREY_BRADLEY", "atom",
                controller, error_by);
            int type_index = Amber_Decode_One_Based_Index(
                tuples[3 * i + 2], ub_type_numbers, "CHARMM_UREY_BRADLEY",
                "parameter", controller, error_by);
            std::size_t matched_index = angle_records.size();
            int match_count = 0;
            for (std::size_t j = 0; j < angle_records.size(); j++)
            {
                const Amber_Angle_Record& angle = angle_records[j];
                bool endpoint_match =
                    (angle.atom_a == atom_a && angle.atom_c == atom_c) ||
                    (angle.atom_a == atom_c && angle.atom_c == atom_a);
                if (endpoint_match && !angle.moved_to_urey_bradley)
                {
                    matched_index = j;
                    match_count++;
                }
            }
            if (match_count != 1)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tCHAMBER Urey-Bradley term %d matches %d "
                    "standard angles; exactly one is required\n",
                    i + 1, match_count);
            }
            Amber_Angle_Record& angle = angle_records[matched_index];
            if (is_virtual_atom[angle.atom_a] ||
                is_virtual_atom[angle.atom_b] || is_virtual_atom[angle.atom_c])
            {
                std::string source_context = Amber_Source_Context(
                    document, {"CHARMM_UREY_BRADLEY", "ANGLES_INC_HYDROGEN",
                               "ANGLES_WITHOUT_HYDROGEN"});
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tAMBER NUMEXTRA virtual sites cannot "
                    "participate in CHARMM_UREY_BRADLEY terms\n%s",
                    source_context.c_str());
            }
            angle.moved_to_urey_bradley = true;
            ff->urey_bradley.atom_a.push_back(angle.atom_a);
            ff->urey_bradley.atom_b.push_back(angle.atom_b);
            ff->urey_bradley.atom_c.push_back(angle.atom_c);
            ff->urey_bradley.angle_k.push_back(angle.k);
            ff->urey_bradley.angle_theta0.push_back(angle.theta0);
            ff->urey_bradley.bond_k.push_back(ub_k[type_index]);
            ff->urey_bradley.bond_r0.push_back(ub_r0[type_index]);
        }
    }
    std::vector<Amber_Angle_Record> chemical_angle_records;
    chemical_angle_records.reserve(angle_records.size());
    for (const Amber_Angle_Record& angle : angle_records)
    {
        if (is_virtual_atom[angle.atom_a] || is_virtual_atom[angle.atom_b] ||
            is_virtual_atom[angle.atom_c])
        {
            // Amber Classic init_extra_pts(frameon=1) trims standard angles
            // involving extra points before force setup.
            continue;
        }
        chemical_angle_records.push_back(angle);
    }
    angle_records = std::move(chemical_angle_records);
    for (const Amber_Angle_Record& angle : angle_records)
    {
        if (!angle.moved_to_urey_bradley)
        {
            Amber_Append_Angle(angle, &ff->angles);
        }
    }

    bool has_plain_cmap = false;
    bool has_charmm_cmap = false;
    for (const auto& named_section : document.sections)
    {
        const std::string& name = named_section.first;
        has_plain_cmap = has_plain_cmap || name == "CMAP_COUNT" ||
                         name == "CMAP_RESOLUTION" || name == "CMAP_INDEX" ||
                         name.rfind("CMAP_PARAMETER_", 0) == 0;
        has_charmm_cmap = has_charmm_cmap || name == "CHARMM_CMAP_COUNT" ||
                          name == "CHARMM_CMAP_RESOLUTION" ||
                          name == "CHARMM_CMAP_INDEX" ||
                          name.rfind("CHARMM_CMAP_PARAMETER_", 0) == 0;
    }
    if (has_plain_cmap && has_charmm_cmap)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER CMAP and CHARMM_CMAP section families cannot be "
            "mixed\n");
    }
    if (has_plain_cmap || has_charmm_cmap)
    {
        std::string prefix = has_charmm_cmap ? "CHARMM_" : "";
        std::string count_name = prefix + "CMAP_COUNT";
        std::string resolution_name = prefix + "CMAP_RESOLUTION";
        std::string index_name = prefix + "CMAP_INDEX";
        const Amber_Section* count_section =
            Amber_Find_Section(document, count_name);
        const Amber_Section* resolution_section =
            Amber_Find_Section(document, resolution_name);
        const Amber_Section* index_section =
            Amber_Find_Section(document, index_name);
        if (count_section == nullptr || resolution_section == nullptr ||
            index_section == nullptr)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tthe AMBER CMAP count, resolution, parameter, and "
                "index sections are incomplete\n");
        }
        Amber_Require_Exact_Section_Size(count_section->values, 2, controller,
                                         error_by, count_name.c_str());
        std::vector<int> counts =
            Amber_Parse_Int_Section(*count_section, controller, error_by);
        int cmap_numbers = counts[0];
        int cmap_type_numbers = counts[1];
        if (cmap_numbers < 0 || cmap_type_numbers < 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER CMAP counts cannot be negative\n");
        }
        Amber_Require_Exact_Section_Size(
            resolution_section->values,
            static_cast<std::size_t>(cmap_type_numbers), controller, error_by,
            resolution_name.c_str());
        ff->cmap.resolution =
            Amber_Parse_Int_Section(*resolution_section, controller, error_by);
        ff->cmap.unique_type_numbers = cmap_type_numbers;
        ff->cmap.type_offset.resize(cmap_type_numbers);
        std::size_t total_gridpoints = 0;
        for (int i = 0; i < cmap_type_numbers; i++)
        {
            int resolution = ff->cmap.resolution[i];
            if (resolution <= 0)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tAMBER CMAP resolution %d must be positive\n",
                    i + 1);
            }
            std::size_t gridpoints =
                Amber_Checked_Product(static_cast<std::size_t>(resolution),
                                      static_cast<std::size_t>(resolution),
                                      "CMAP grid", controller, error_by);
            const std::size_t maximum_gridpoints =
                static_cast<std::size_t>(std::numeric_limits<int>::max()) / 16;
            if (total_gridpoints > maximum_gridpoints ||
                gridpoints > maximum_gridpoints - total_gridpoints)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tAMBER CMAP interpolation table exceeds the "
                    "supported kernel index range\n");
            }
            ff->cmap.type_offset[i] = static_cast<int>(16 * total_gridpoints);
            total_gridpoints += gridpoints;

            char suffix[32];
            snprintf(suffix, sizeof(suffix), "CMAP_PARAMETER_%02d", i + 1);
            std::string parameter_name = prefix + suffix;
            const Amber_Section* parameter_section =
                Amber_Find_Section(document, parameter_name);
            if (parameter_section == nullptr)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\trequired AMBER CMAP section %s is missing\n",
                    parameter_name.c_str());
            }
            Amber_Require_Exact_Section_Size(parameter_section->values,
                                             gridpoints, controller, error_by,
                                             parameter_name.c_str());
            std::vector<float> grid = Amber_Parse_Float_Section(
                *parameter_section, controller, error_by);
            ff->cmap.grid_value.insert(ff->cmap.grid_value.end(), grid.begin(),
                                       grid.end());
        }
        std::string parameter_prefix = prefix + "CMAP_PARAMETER_";
        int actual_parameter_sections = 0;
        for (const auto& named_section : document.sections)
        {
            if (named_section.first.rfind(parameter_prefix, 0) == 0)
            {
                actual_parameter_sections++;
            }
        }
        if (actual_parameter_sections != cmap_type_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER CMAP has %d parameter sections; count "
                "declares %d\n",
                actual_parameter_sections, cmap_type_numbers);
        }
        ff->cmap.unique_gridpoint_numbers = static_cast<int>(total_gridpoints);

        std::size_t expected_index_values =
            Amber_Checked_Product(static_cast<std::size_t>(cmap_numbers), 6,
                                  "CMAP index", controller, error_by);
        Amber_Require_Exact_Section_Size(index_section->values,
                                         expected_index_values, controller,
                                         error_by, index_name.c_str());
        std::vector<int> indices =
            Amber_Parse_Int_Section(*index_section, controller, error_by);
        ff->cmap.atom_a.resize(cmap_numbers);
        ff->cmap.atom_b.resize(cmap_numbers);
        ff->cmap.atom_c.resize(cmap_numbers);
        ff->cmap.atom_d.resize(cmap_numbers);
        ff->cmap.atom_e.resize(cmap_numbers);
        ff->cmap.cmap_type.resize(cmap_numbers);
        for (int i = 0; i < cmap_numbers; i++)
        {
            ff->cmap.atom_a[i] = Amber_Decode_One_Based_Index(
                indices[6 * i], atom_numbers, index_name.c_str(), "atom",
                controller, error_by);
            ff->cmap.atom_b[i] = Amber_Decode_One_Based_Index(
                indices[6 * i + 1], atom_numbers, index_name.c_str(), "atom",
                controller, error_by);
            ff->cmap.atom_c[i] = Amber_Decode_One_Based_Index(
                indices[6 * i + 2], atom_numbers, index_name.c_str(), "atom",
                controller, error_by);
            ff->cmap.atom_d[i] = Amber_Decode_One_Based_Index(
                indices[6 * i + 3], atom_numbers, index_name.c_str(), "atom",
                controller, error_by);
            ff->cmap.atom_e[i] = Amber_Decode_One_Based_Index(
                indices[6 * i + 4], atom_numbers, index_name.c_str(), "atom",
                controller, error_by);
            ff->cmap.cmap_type[i] = Amber_Decode_One_Based_Index(
                indices[6 * i + 5], cmap_type_numbers, index_name.c_str(),
                "parameter", controller, error_by);
            if (is_virtual_atom[ff->cmap.atom_a[i]] ||
                is_virtual_atom[ff->cmap.atom_b[i]] ||
                is_virtual_atom[ff->cmap.atom_c[i]] ||
                is_virtual_atom[ff->cmap.atom_d[i]] ||
                is_virtual_atom[ff->cmap.atom_e[i]])
            {
                std::string source_context =
                    Amber_Source_Context(document, {index_name.c_str()});
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tAMBER NUMEXTRA virtual sites cannot "
                    "participate in CMAP terms\n%s",
                    source_context.c_str());
            }
        }
    }

    bool has_dihedral_pk =
        Amber_Find_Section(document, "DIHEDRAL_FORCE_CONSTANT") != nullptr;
    bool has_dihedral_phase =
        Amber_Find_Section(document, "DIHEDRAL_PHASE") != nullptr;
    bool has_dihedral_periodicity =
        Amber_Find_Section(document, "DIHEDRAL_PERIODICITY") != nullptr;
    int dihedral_parameter_flags = static_cast<int>(has_dihedral_pk) +
                                   static_cast<int>(has_dihedral_phase) +
                                   static_cast<int>(has_dihedral_periodicity);
    if ((dihedral_parameter_flags != 0 && dihedral_parameter_flags != 3) ||
        (dihedral_type_numbers > 0 && dihedral_parameter_flags != 3))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tDIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY, and "
            "DIHEDRAL_PHASE must all be present when AMBER dihedral parameter "
            "types are declared\n");
    }
    // Amber Classic reads these arrays as double precision.  Keep that
    // precision through dihpar-compatible phase/coefficient preprocessing;
    // narrowing first can move values across its historical thresholds.
    std::vector<double> dihedral_type_pk = Amber_Load_Optional_Double_Table(
        document, "DIHEDRAL_FORCE_CONSTANT",
        static_cast<std::size_t>(dihedral_type_numbers), false, controller,
        error_by);
    std::vector<double> dihedral_type_phase = Amber_Load_Optional_Double_Table(
        document, "DIHEDRAL_PHASE",
        static_cast<std::size_t>(dihedral_type_numbers), false, controller,
        error_by);
    std::vector<double> dihedral_type_periodicity =
        Amber_Load_Optional_Double_Table(
            document, "DIHEDRAL_PERIODICITY",
            static_cast<std::size_t>(dihedral_type_numbers), false, controller,
            error_by);
    const Amber_Section* scee_section =
        Amber_Find_Section(document, "SCEE_SCALE_FACTOR");
    const Amber_Section* scnb_section =
        Amber_Find_Section(document, "SCNB_SCALE_FACTOR");
    std::vector<float> scee_scale_factor = Amber_Load_Optional_Float_Table(
        document, "SCEE_SCALE_FACTOR",
        static_cast<std::size_t>(dihedral_type_numbers), false, controller,
        error_by);
    std::vector<float> scnb_scale_factor = Amber_Load_Optional_Float_Table(
        document, "SCNB_SCALE_FACTOR",
        static_cast<std::size_t>(dihedral_type_numbers), false, controller,
        error_by);

    std::vector<int> raw_dihedrals_h = Amber_Load_Interaction_Tuples(
        document, "DIHEDRALS_INC_HYDROGEN", dihedral_with_hydrogen, 5,
        controller, error_by);
    std::vector<int> raw_dihedrals_no_h = Amber_Load_Interaction_Tuples(
        document, "DIHEDRALS_WITHOUT_HYDROGEN", dihedral_without_hydrogen, 5,
        controller, error_by);
    std::vector<int> raw_dihedrals = raw_dihedrals_h;
    raw_dihedrals.insert(raw_dihedrals.end(), raw_dihedrals_no_h.begin(),
                         raw_dihedrals_no_h.end());

    const Amber_Section* atom_type_section =
        Amber_Find_Section(document, "ATOM_TYPE_INDEX");
    if (atom_numbers > 0 && atom_type_section == nullptr)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tATOM_TYPE_INDEX is required for an AMBER topology "
            "with atoms\n");
    }
    if (atom_type_section != nullptr)
    {
        Amber_Require_Exact_Section_Size(
            atom_type_section->values, static_cast<std::size_t>(atom_numbers),
            controller, error_by, "ATOM_TYPE_INDEX");
        std::vector<int> atom_types =
            Amber_Parse_Int_Section(*atom_type_section, controller, error_by);
        ff->lj.atom_type.resize(atom_numbers);
        for (int i = 0; i < atom_numbers; i++)
        {
            ff->lj.atom_type[i] = Amber_Decode_One_Based_Index(
                atom_types[i], atom_type_numbers, "ATOM_TYPE_INDEX", "type",
                controller, error_by);
        }
    }
    ff->lj.atom_type_numbers = atom_type_numbers;

    std::size_t type_count = static_cast<std::size_t>(atom_type_numbers);
    if (type_count > 0 &&
        type_count > std::numeric_limits<std::size_t>::max() / (type_count + 1))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tAMBER Lennard-Jones pair table size overflows\n");
    }
    std::size_t pair_type_numbers = type_count * (type_count + 1) / 2;
    bool has_lj_a =
        Amber_Find_Section(document, "LENNARD_JONES_ACOEF") != nullptr;
    bool has_lj_b =
        Amber_Find_Section(document, "LENNARD_JONES_BCOEF") != nullptr;
    if (has_lj_a != has_lj_b || (atom_numbers > 0 && !has_lj_a))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tLENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF are "
            "required for an AMBER topology with atoms\n");
    }
    std::vector<float> raw_lj_a = Amber_Load_Optional_Float_Table(
        document, "LENNARD_JONES_ACOEF", pair_type_numbers, false, controller,
        error_by);
    std::vector<float> raw_lj_b = Amber_Load_Optional_Float_Table(
        document, "LENNARD_JONES_BCOEF", pair_type_numbers, false, controller,
        error_by);
    for (std::size_t i = 0; i < raw_lj_a.size(); i++)
    {
        if (raw_lj_a[i] < 0.0f || raw_lj_b[i] < 0.0f)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER Lennard-Jones pair parameter %zu is "
                "negative\n",
                i + 1);
        }
    }

    const Amber_Section* npi_section =
        Amber_Find_Section(document, "NONBONDED_PARM_INDEX");
    std::vector<int> nonbonded_parm_index;
    if (npi_section != nullptr)
    {
        nonbonded_parm_index =
            Amber_Parse_Int_Section(*npi_section, controller, error_by);
    }
    std::vector<float> hbond_a = Amber_Load_Optional_Float_Table(
        document, "HBOND_ACOEF", static_cast<std::size_t>(hbond_type_numbers),
        false, controller, error_by);
    std::vector<float> hbond_b = Amber_Load_Optional_Float_Table(
        document, "HBOND_BCOEF", static_cast<std::size_t>(hbond_type_numbers),
        false, controller, error_by);
    bool has_hbond_a = Amber_Find_Section(document, "HBOND_ACOEF") != nullptr;
    bool has_hbond_b = Amber_Find_Section(document, "HBOND_BCOEF") != nullptr;
    if (has_hbond_a != has_hbond_b || (hbond_type_numbers > 0 && !has_hbond_a))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tHBOND_ACOEF and HBOND_BCOEF must both be present when "
            "AMBER HBOND parameters are declared\n");
    }

    std::vector<int> canonical_nonbonded_map =
        Amber_Build_Canonical_Nonbonded_Map(
            atom_type_numbers, nonbonded_parm_index, hbond_a, hbond_b,
            document.path.c_str(), npi_section, controller);
    Amber_Remap_LJ_Pair_Matrix(raw_lj_a, raw_lj_b, has_lj_a, has_lj_b,
                               canonical_nonbonded_map, &ff->lj.pair_A,
                               &ff->lj.pair_B, controller,
                               "LENNARD_JONES_ACOEF", "LENNARD_JONES_BCOEF");

    bool has_lj14_a =
        Amber_Find_Section(document, "LENNARD_JONES_14_ACOEF") != nullptr;
    bool has_lj14_b =
        Amber_Find_Section(document, "LENNARD_JONES_14_BCOEF") != nullptr;
    std::vector<float> raw_lj14_a = Amber_Load_Optional_Float_Table(
        document, "LENNARD_JONES_14_ACOEF", pair_type_numbers, false,
        controller, error_by);
    std::vector<float> raw_lj14_b = Amber_Load_Optional_Float_Table(
        document, "LENNARD_JONES_14_BCOEF", pair_type_numbers, false,
        controller, error_by);
    std::vector<float> lj14_a;
    std::vector<float> lj14_b;
    Amber_Remap_LJ_Pair_Matrix(raw_lj14_a, raw_lj14_b, has_lj14_a, has_lj14_b,
                               canonical_nonbonded_map, &lj14_a, &lj14_b,
                               controller, "LENNARD_JONES_14_ACOEF",
                               "LENNARD_JONES_14_BCOEF");
    for (std::size_t i = 0; i < raw_lj14_a.size(); i++)
    {
        if (raw_lj14_a[i] < 0.0f || raw_lj14_b[i] < 0.0f)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER 1-4 Lennard-Jones pair parameter %zu is "
                "negative\n",
                i + 1);
        }
    }

    const Amber_Section* ccoef_section =
        Amber_Find_Section(document, "LENNARD_JONES_CCOEF");
    if (ccoef_section != nullptr)
    {
        Amber_Require_Exact_Section_Size(ccoef_section->values,
                                         pair_type_numbers, controller,
                                         error_by, "LENNARD_JONES_CCOEF");
        std::vector<float> ccoef =
            Amber_Parse_Float_Section(*ccoef_section, controller, error_by);
        for (std::size_t i = 0; i < ccoef.size(); i++)
        {
            if (ccoef[i] != 0.0f)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tnonzero AMBER LENNARD_JONES_CCOEF entry %zu "
                    "requires unsupported 12-6-4 nonbonded interactions\n"
                    "\tInput file: %s\n\t%%FLAG LENNARD_JONES_CCOEF begins on "
                    "line %zu\n",
                    i + 1, document.path.c_str(), ccoef_section->flag_line);
            }
        }
    }

    const char* improper_names[] = {
        "CHARMM_NUM_IMPROPERS", "CHARMM_IMPROPERS", "CHARMM_NUM_IMPR_TYPES",
        "CHARMM_IMPROPER_FORCE_CONSTANT", "CHARMM_IMPROPER_PHASE"};
    int improper_flag_count = 0;
    for (const char* name : improper_names)
    {
        improper_flag_count += Amber_Find_Section(document, name) != nullptr;
    }
    if (improper_flag_count != 0 && improper_flag_count != 5)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe CHAMBER improper sections are incomplete\n");
    }
    if (improper_flag_count == 5)
    {
        const Amber_Section& count_section = Amber_Require_Section(
            document, "CHARMM_NUM_IMPROPERS", controller, error_by);
        const Amber_Section& type_count_section = Amber_Require_Section(
            document, "CHARMM_NUM_IMPR_TYPES", controller, error_by);
        Amber_Require_Exact_Section_Size(count_section.values, 1, controller,
                                         error_by, "CHARMM_NUM_IMPROPERS");
        Amber_Require_Exact_Section_Size(type_count_section.values, 1,
                                         controller, error_by,
                                         "CHARMM_NUM_IMPR_TYPES");
        int improper_numbers =
            Amber_Parse_Int(count_section.values[0], "CHARMM_NUM_IMPROPERS",
                            controller, error_by);
        int improper_type_numbers =
            Amber_Parse_Int(type_count_section.values[0],
                            "CHARMM_NUM_IMPR_TYPES", controller, error_by);
        if (improper_numbers < 0 || improper_type_numbers < 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tCHAMBER improper counts cannot be negative\n");
        }
        const Amber_Section& tuple_section = Amber_Require_Section(
            document, "CHARMM_IMPROPERS", controller, error_by);
        std::size_t expected_values =
            Amber_Checked_Product(static_cast<std::size_t>(improper_numbers), 5,
                                  "CHARMM_IMPROPERS", controller, error_by);
        Amber_Require_Exact_Section_Size(tuple_section.values, expected_values,
                                         controller, error_by,
                                         "CHARMM_IMPROPERS");
        std::vector<int> tuples =
            Amber_Parse_Int_Section(tuple_section, controller, error_by);
        std::vector<float> improper_k = Amber_Load_Optional_Float_Table(
            document, "CHARMM_IMPROPER_FORCE_CONSTANT",
            static_cast<std::size_t>(improper_type_numbers), true, controller,
            error_by);
        std::vector<float> improper_phase = Amber_Load_Optional_Float_Table(
            document, "CHARMM_IMPROPER_PHASE",
            static_cast<std::size_t>(improper_type_numbers), true, controller,
            error_by);
        for (int i = 0; i < improper_type_numbers; i++)
        {
            if (improper_k[i] < 0.0f)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tCHAMBER improper force constant %d cannot be "
                    "negative\n",
                    i + 1);
            }
        }
        for (int i = 0; i < improper_numbers; i++)
        {
            int atom_a = Amber_Decode_One_Based_Index(
                tuples[5 * i], atom_numbers, "CHARMM_IMPROPERS", "atom",
                controller, error_by);
            int atom_b = Amber_Decode_One_Based_Index(
                tuples[5 * i + 1], atom_numbers, "CHARMM_IMPROPERS", "atom",
                controller, error_by);
            int atom_c = Amber_Decode_One_Based_Index(
                tuples[5 * i + 2], atom_numbers, "CHARMM_IMPROPERS", "atom",
                controller, error_by);
            int atom_d = Amber_Decode_One_Based_Index(
                tuples[5 * i + 3], atom_numbers, "CHARMM_IMPROPERS", "atom",
                controller, error_by);
            int type_index = Amber_Decode_One_Based_Index(
                tuples[5 * i + 4], improper_type_numbers, "CHARMM_IMPROPERS",
                "parameter", controller, error_by);
            if (is_virtual_atom[atom_a] || is_virtual_atom[atom_b] ||
                is_virtual_atom[atom_c] || is_virtual_atom[atom_d])
            {
                std::string source_context =
                    Amber_Source_Context(document, {"CHARMM_IMPROPERS"});
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tAMBER NUMEXTRA virtual sites cannot "
                    "participate in CHARMM_IMPROPERS terms\n%s",
                    source_context.c_str());
            }
            if (atom_a == atom_b || atom_a == atom_c || atom_a == atom_d ||
                atom_b == atom_c || atom_b == atom_d || atom_c == atom_d)
            {
                std::string source_context =
                    Amber_Source_Context(document, {"CHARMM_IMPROPERS"});
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tCHARMM_IMPROPERS term %d contains repeated "
                    "atoms\n%s",
                    i + 1, source_context.c_str());
            }
            ff->impropers.atom_a.push_back(atom_a);
            ff->impropers.atom_b.push_back(atom_b);
            ff->impropers.atom_c.push_back(atom_c);
            ff->impropers.atom_d.push_back(atom_d);
            ff->impropers.pk.push_back(improper_k[type_index]);
            ff->impropers.pn.push_back(0.0f);
            ff->impropers.ipn.push_back(0);
            ff->impropers.gamc.push_back(improper_phase[type_index]);
            ff->impropers.gams.push_back(0.0f);
        }
    }

    const std::vector<float>& nb14_pair_a = has_lj14_a ? lj14_a : ff->lj.pair_A;
    const std::vector<float>& nb14_pair_b = has_lj14_b ? lj14_b : ff->lj.pair_B;
    const char* nb14_pair_a_section =
        has_lj14_a ? "LENNARD_JONES_14_ACOEF" : "LENNARD_JONES_ACOEF";
    const char* nb14_pair_b_section =
        has_lj14_b ? "LENNARD_JONES_14_BCOEF" : "LENNARD_JONES_BCOEF";
    struct Amber_NB14_Source
    {
        std::string section;
        std::size_t term;
        int atom_a;
        int atom_b;
        int atom_c;
        int atom_d;
    };
    std::map<long long, Amber_NB14_Source> nb14_records;
    const double amber_pi = std::acos(-1.0);
    for (std::size_t i = 0; i < raw_dihedrals.size(); i += 5)
    {
        const char* source_section = i < raw_dihedrals_h.size()
                                         ? "DIHEDRALS_INC_HYDROGEN"
                                         : "DIHEDRALS_WITHOUT_HYDROGEN";
        std::size_t section_offset =
            i < raw_dihedrals_h.size() ? i : i - raw_dihedrals_h.size();
        const std::size_t source_term = section_offset / 5 + 1;
        const std::string term_source_context =
            Amber_Source_Context(document, {source_section});
        int atom_a = Amber_Decode_Encoded_Atom(
            raw_dihedrals[i], false, atom_numbers, "dihedral", "A", controller,
            error_by, term_source_context.c_str());
        int atom_b = Amber_Decode_Encoded_Atom(
            raw_dihedrals[i + 1], false, atom_numbers, "dihedral", "B",
            controller, error_by, term_source_context.c_str());
        int atom_c = Amber_Decode_Encoded_Atom(
            raw_dihedrals[i + 2], true, atom_numbers, "dihedral", "C",
            controller, error_by, term_source_context.c_str());
        int atom_d = Amber_Decode_Encoded_Atom(
            raw_dihedrals[i + 3], true, atom_numbers, "dihedral", "D",
            controller, error_by, term_source_context.c_str());
        const int raw_type_index = raw_dihedrals[i + 4];
        if (raw_type_index <= 0 || raw_type_index > dihedral_type_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER %s term %zu contains an out-of-range "
                "parameter index %d; expected [1, %d]\n%s",
                source_section, source_term, raw_type_index,
                dihedral_type_numbers, term_source_context.c_str());
        }
        int type_index = raw_type_index - 1;
        const std::string parameter_reference =
            "\tAMBER dihedral parameter type " +
            std::to_string(type_index + 1) + " is referenced by " +
            source_section + " term " + std::to_string(source_term) + "\n";
        double raw_pn = dihedral_type_periodicity[type_index];
        if (raw_pn < 0.0)
        {
            std::string source_context = Amber_Source_Context(
                document, {"DIHEDRAL_PERIODICITY", source_section});
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER dihedral parameter type %d has negative "
                "periodicity but is referenced by %s term %zu\n%s",
                type_index + 1, source_section, source_term,
                source_context.c_str());
        }
        if (is_virtual_atom[atom_a] || is_virtual_atom[atom_b] ||
            is_virtual_atom[atom_c] || is_virtual_atom[atom_d])
        {
            // Amber Classic init_extra_pts(frameon=1) trims every standard
            // bonded term involving an extra point before force setup.  This
            // also prevents a trimmed proper from generating a 1-4 pair.
            continue;
        }
        if (atom_a == atom_b || atom_a == atom_c || atom_a == atom_d ||
            atom_b == atom_c || atom_b == atom_d || atom_c == atom_d)
        {
            std::string source_context =
                Amber_Source_Context(document, {source_section});
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER %s term %zu contains repeated atoms\n%s",
                source_section, source_term, source_context.c_str());
        }
        bool ignore_end = raw_dihedrals[i + 2] < 0;
        bool periodic_improper = raw_dihedrals[i + 3] < 0;
        double pn_double = std::fabs(raw_pn);
        double legacy_ipn = pn_double + 1.0e-3;
        if (legacy_ipn > static_cast<double>(std::numeric_limits<int>::max()))
        {
            std::string source_context = Amber_Source_Context(
                document, {"DIHEDRAL_PERIODICITY", source_section});
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER dihedral parameter type %d periodicity "
                "is outside the legacy integer metadata range\n%s",
                type_index + 1, source_context.c_str());
        }
        const std::string pk_source_context =
            parameter_reference +
            Amber_Source_Context(document,
                                 {"DIHEDRAL_FORCE_CONSTANT", source_section});
        const std::string pn_source_context =
            parameter_reference +
            Amber_Source_Context(document,
                                 {"DIHEDRAL_PERIODICITY", source_section});
        const std::string phase_source_context =
            parameter_reference +
            Amber_Source_Context(document, {"DIHEDRAL_FORCE_CONSTANT",
                                            "DIHEDRAL_PHASE", source_section});
        float pk = Amber_Checked_Float(dihedral_type_pk[type_index],
                                       "DIHEDRAL_FORCE_CONSTANT", controller,
                                       error_by, pk_source_context.c_str());
        float pn =
            Amber_Checked_Float(pn_double, "DIHEDRAL_PERIODICITY", controller,
                                error_by, pn_source_context.c_str());

        // Match Amber Classic set.F90::dihpar exactly: phase snapping and
        // trig-component zeroing happen in double precision before the force
        // constant is applied.  In particular, a small force constant must
        // not cause an otherwise valid cosine/sine component to be erased.
        double phase = dihedral_type_phase[type_index];
        if (std::fabs(phase - amber_pi) <= 1.0e-3)
        {
            phase = std::copysign(amber_pi, phase);
        }
        double phase_cosine = std::cos(phase);
        double phase_sine = std::sin(phase);
        if (std::fabs(phase_cosine) <= 1.0e-6)
        {
            phase_cosine = 0.0;
        }
        if (std::fabs(phase_sine) <= 1.0e-6)
        {
            phase_sine = 0.0;
        }
        float gamc = Amber_Checked_Float(
            phase_cosine * dihedral_type_pk[type_index],
            "DIHEDRAL_FORCE_CONSTANT and DIHEDRAL_PHASE", controller, error_by,
            phase_source_context.c_str());
        float gams = Amber_Checked_Float(
            phase_sine * dihedral_type_pk[type_index],
            "DIHEDRAL_FORCE_CONSTANT and DIHEDRAL_PHASE", controller, error_by,
            phase_source_context.c_str());
        ff->dihedrals.atom_a.push_back(atom_a);
        ff->dihedrals.atom_b.push_back(atom_b);
        ff->dihedrals.atom_c.push_back(atom_c);
        ff->dihedrals.atom_d.push_back(atom_d);
        ff->dihedrals.pk.push_back(pk);
        ff->dihedrals.pn.push_back(pn);
        // ipn is retained only as Amber-compatible metadata.  The force and
        // energy kernels intentionally use the real-valued pn above.
        ff->dihedrals.ipn.push_back(static_cast<int>(legacy_ipn));
        ff->dihedrals.gamc.push_back(gamc);
        ff->dihedrals.gams.push_back(gams);

        if (ignore_end || periodic_improper)
        {
            continue;
        }
        float scee =
            scee_section == nullptr ? 1.2f : scee_scale_factor[type_index];
        float scnb =
            scnb_section == nullptr ? 2.0f : scnb_scale_factor[type_index];
        if (scee <= 0.0f || scnb <= 0.0f)
        {
            const std::string scale_source_context =
                parameter_reference +
                Amber_Source_Context(
                    document,
                    {"SCEE_SCALE_FACTOR", "SCNB_SCALE_FACTOR", source_section});
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tSCEE_SCALE_FACTOR and SCNB_SCALE_FACTOR for "
                "proper dihedral type %d must be finite and positive\n%s",
                type_index + 1, scale_source_context.c_str());
        }
        int type_a = ff->lj.atom_type[atom_a];
        int type_d = ff->lj.atom_type[atom_d];
        if (type_a > type_d)
        {
            std::swap(type_a, type_d);
        }
        std::size_t pair_type =
            static_cast<std::size_t>(type_d) * (type_d + 1) / 2 + type_a;
        if (pair_type >= nb14_pair_a.size() || pair_type >= nb14_pair_b.size())
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tAMBER 1-4 Lennard-Jones parameter lookup is out "
                "of range\n");
        }
        int endpoint_low = std::min(atom_a, atom_d);
        int endpoint_high = std::max(atom_a, atom_d);
        long long endpoint_key =
            static_cast<long long>(endpoint_low) * atom_numbers + endpoint_high;
        auto existing = nb14_records.find(endpoint_key);
        if (existing != nb14_records.end())
        {
            const Amber_NB14_Source& first = existing->second;
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tduplicate active AMBER 1-4 endpoint pair %d-%d\n"
                "\tFirst active proper: %s term %zu (atoms %d %d %d %d)\n"
                "\tSecond active proper: %s term %zu (atoms %d %d %d %d)\n"
                "\tOnly one active proper may generate a given 1-4 pair; "
                "continuation terms must encode atom C with a negative "
                "index\n\tInput file: %s\n",
                endpoint_low + 1, endpoint_high + 1, first.section.c_str(),
                first.term, first.atom_a + 1, first.atom_b + 1,
                first.atom_c + 1, first.atom_d + 1, source_section, source_term,
                atom_a + 1, atom_b + 1, atom_c + 1, atom_d + 1,
                document.path.c_str());
        }
        nb14_records.emplace(
            endpoint_key, Amber_NB14_Source{source_section, source_term, atom_a,
                                            atom_b, atom_c, atom_d});
        const std::string scee_source_context =
            parameter_reference +
            Amber_Source_Context(document,
                                 {"SCEE_SCALE_FACTOR", source_section});
        const std::string scnb_source_context =
            parameter_reference +
            Amber_Source_Context(document,
                                 {"SCNB_SCALE_FACTOR", source_section});
        const std::string nb14_a_source_context =
            parameter_reference +
            Amber_Source_Context(
                document,
                {nb14_pair_a_section, "SCNB_SCALE_FACTOR", source_section});
        const std::string nb14_b_source_context =
            parameter_reference +
            Amber_Source_Context(
                document,
                {nb14_pair_b_section, "SCNB_SCALE_FACTOR", source_section});
        float inverse_scnb = Amber_Checked_Float(
            1.0 / static_cast<double>(scnb), "SCNB_SCALE_FACTOR", controller,
            error_by, scnb_source_context.c_str());
        float inverse_scee = Amber_Checked_Float(
            1.0 / static_cast<double>(scee), "SCEE_SCALE_FACTOR", controller,
            error_by, scee_source_context.c_str());
        ff->nb14.atom_a.push_back(atom_a);
        ff->nb14.atom_b.push_back(atom_d);
        ff->nb14.A.push_back(Amber_Checked_Scale(
            nb14_pair_a[pair_type], inverse_scnb, nb14_pair_a_section,
            controller, error_by, nb14_a_source_context.c_str()));
        ff->nb14.B.push_back(Amber_Checked_Scale(
            nb14_pair_b[pair_type], inverse_scnb, nb14_pair_b_section,
            controller, error_by, nb14_b_source_context.c_str()));
        ff->nb14.cf_scale_factor.push_back(inverse_scee);
    }
}

static void Amber_Load_Rst7(System* system, const std::vector<int>* pointers,
                            CONTROLLER* controller,
                            bool* request_nopbc_default)
{
    if (!controller->Command_Exist("amber_rst7"))
    {
        return;
    }
    const char* error_by = "Xponge::Amber_Load_Rst7";
    std::ifstream fin(controller->Original_Command("amber_rst7"));
    if (!fin.is_open())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed to open amber_rst7\n");
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fin, line))
    {
        lines.push_back(line);
    }
    if (fin.bad())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed while reading amber_rst7\n");
    }
    fin.clear();
    fin.close();
    if (fin.fail())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed while closing amber_rst7\n");
    }
    while (lines.size() > 2 && Amber_Trim(lines.back()).empty())
    {
        lines.pop_back();
    }
    if (lines.size() < 2)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tamber_rst7 must contain a title and atom-count "
            "line\n");
    }
    std::vector<std::string> header_tokens =
        Amber_Tokenize_Restart_Header(lines[1], 2, controller);
    if (header_tokens.size() != 1 && header_tokens.size() != 2)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tamber_rst7 atom-count line must contain NATOM and "
            "optional time only\n");
    }
    int atom_numbers = Amber_Parse_Int(header_tokens[0], "amber_rst7 header",
                                       controller, error_by);
    if (atom_numbers <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tamber_rst7 NATOM must be positive\n");
    }
    Amber_Ensure_Atom_Numbers(system, atom_numbers, controller, error_by);
    system->start_time =
        header_tokens.size() == 2
            ? Amber_Parse_Double(header_tokens[1], "amber_rst7 time",
                                 controller, error_by)
            : 0.0;

    std::size_t coordinate_count =
        Amber_Checked_Product(static_cast<std::size_t>(atom_numbers), 3,
                              "restart coordinate", controller, error_by);
    std::size_t coordinate_line_count =
        coordinate_count / 6 + (coordinate_count % 6 != 0 ? 1 : 0);
    std::size_t data_line_count = lines.size() - 2;
    if (data_line_count < coordinate_line_count)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tamber_rst7 has %zu data lines after the header; at "
            "least %zu coordinate lines are required\n",
            data_line_count, coordinate_line_count);
    }
    std::size_t tail_line_count = data_line_count - coordinate_line_count;
    int ifbox = pointers == nullptr ? -1 : pointers->at(27);
    struct Restart_Layout
    {
        bool has_velocity;
        bool has_box;
    };
    std::vector<Restart_Layout> layouts;
    for (int has_velocity = 0; has_velocity <= 1; has_velocity++)
    {
        for (int has_box = 0; has_box <= 1; has_box++)
        {
            if ((ifbox == 0 && has_box != 0) || (ifbox > 0 && has_box == 0))
            {
                continue;
            }
            std::size_t expected_tail_lines =
                (has_velocity != 0 ? coordinate_line_count : 0) +
                (has_box != 0 ? 1 : 0);
            if (tail_line_count == expected_tail_lines)
            {
                if (ifbox < 0 && atom_numbers == 1 && tail_line_count == 1)
                {
                    // With one atom, a velocity block has three fields while
                    // a box has six, even though each occupies one line.
                    // rst7-only input has no IFBOX pointer, so distinguish the
                    // two layouts by the record itself instead of line count.
                    std::vector<std::string> candidate_fields;
                    std::size_t expected_fields =
                        has_velocity != 0 && has_box == 0 ? 3 : 6;
                    if (!Amber_Try_Parse_Restart_Data_Line(
                            lines[2 + coordinate_line_count], expected_fields,
                            &candidate_fields))
                    {
                        continue;
                    }
                }
                layouts.push_back({has_velocity != 0, has_box != 0});
            }
        }
    }
    if (layouts.size() != 1)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tamber_rst7 has %zu data lines after coordinates; this "
            "does not identify one unambiguous velocity/box layout for "
            "IFBOX=%d\n",
            tail_line_count, ifbox);
    }
    Restart_Layout layout = layouts[0];
    std::size_t line_index = 2;
    system->atoms.coordinate = Amber_Parse_Restart_Block(
        lines, &line_index, coordinate_count, "coordinate", controller);
    if (layout.has_velocity)
    {
        // Amber stores physical A/ps velocities divided by 20.455. SPONGE's
        // internal time unit is ps/20.455, so the file values are already in
        // SPONGE's internal velocity unit and must not be rescaled here.
        system->atoms.velocity = Amber_Parse_Restart_Block(
            lines, &line_index, coordinate_count, "velocity", controller);
    }
    else
    {
        system->atoms.velocity.assign(coordinate_count, 0.0f);
    }
    if (layout.has_box)
    {
        std::vector<float> box_values =
            Amber_Parse_Restart_Block(lines, &line_index, 6, "box", controller);
        system->box.box_length.assign(box_values.begin(),
                                      box_values.begin() + 3);
        system->box.box_angle.assign(box_values.begin() + 3, box_values.end());
        for (int i = 0; i < 3; i++)
        {
            if (system->box.box_length[i] <= 0.0f ||
                system->box.box_angle[i] <= 0.0f ||
                system->box.box_angle[i] >= 180.0f)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tamber_rst7 box lengths must be positive and "
                    "angles must lie in (0, 180) degrees\n");
            }
        }
        double alpha = static_cast<double>(system->box.box_angle[0]) *
                       CONSTANT_DEG_TO_RAD_DOUBLE;
        double beta = static_cast<double>(system->box.box_angle[1]) *
                      CONSTANT_DEG_TO_RAD_DOUBLE;
        double gamma = static_cast<double>(system->box.box_angle[2]) *
                       CONSTANT_DEG_TO_RAD_DOUBLE;
        double cos_alpha = std::cos(alpha);
        double cos_beta = std::cos(beta);
        double cos_gamma = std::cos(gamma);
        double determinant = 1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta -
                             cos_gamma * cos_gamma +
                             2.0 * cos_alpha * cos_beta * cos_gamma;
        if (!Amber_Is_Finite_Double(determinant) || determinant <= 0.0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tamber_rst7 box angles do not define a valid "
                "triclinic cell\n");
        }
    }
    else
    {
        if (controller->Command_Exist("pbc"))
        {
            if (controller->Get_Bool("pbc", error_by))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorConflictingCommand, error_by,
                    "Reason:\n\tAMBER IFBOX=0 restart is non-periodic but "
                    "pbc=true was requested\n");
            }
        }
        else
        {
            *request_nopbc_default = true;
        }
        system->box.box_length.resize(3);
        system->box.box_angle.assign(3, 90.0f);
        double axis_scale[3] = {};
        double largest_scale = 0.0;
        for (int axis = 0; axis < 3; axis++)
        {
            double min_coordinate = system->atoms.coordinate[axis];
            double max_coordinate = system->atoms.coordinate[axis];
            for (int atom = 1; atom < atom_numbers; atom++)
            {
                double coordinate = system->atoms.coordinate[3 * atom + axis];
                min_coordinate = std::min(min_coordinate, coordinate);
                max_coordinate = std::max(max_coordinate, coordinate);
            }
            double absolute_extent =
                std::max(std::fabs(min_coordinate), std::fabs(max_coordinate));
            double span = max_coordinate - min_coordinate;
            axis_scale[axis] = std::max(absolute_extent, span);
            largest_scale = std::max(largest_scale, axis_scale[axis]);
        }
        // NOPBC kernels now disable image translations explicitly, so this
        // box is metadata only.  Derive its scale from the coordinates rather
        // than fabricating a 100/1000-Angstrom safety margin.  A system whose
        // coordinates are all exactly at the origin has no intrinsic length;
        // use one Angstrom solely as a finite unit-scale placeholder.
        if (largest_scale == 0.0) largest_scale = 1.0;
        for (int axis = 0; axis < 3; axis++)
        {
            const double auxiliary_length =
                axis_scale[axis] == 0.0 ? largest_scale : axis_scale[axis];
            system->box.box_length[axis] = Amber_Checked_Float(
                auxiliary_length, "amber_rst7 NOPBC auxiliary cell", controller,
                error_by);
        }
    }
    if (line_index != lines.size())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tamber_rst7 contains trailing data lines\n");
    }
}

void Load_Amber_Inputs(System* system, CONTROLLER* controller)
{
    bool request_nopbc_default = false;
    Load_System_Transaction(
        system, controller, "Xponge::Load_Amber_Inputs",
        Load_System_Seed::kCurrent,
        [&](System* staged)
        {
            Amber_Parm7_Document document;
            std::vector<int> pointers;
            const bool has_parm7 =
                controller->Command_Exist("amber_parm7");
            const bool has_rst7 = controller->Command_Exist("amber_rst7");

            // A provided AMBER source replaces the fields owned by that
            // source.  Clear them before cross-file atom-count validation so
            // a previous load cannot either block a legitimate replacement
            // or leak stale force-field/coordinate data into the new System.
            // Omitted parm7/rst7 halves remain intentionally incremental.
            if (has_parm7)
            {
                staged->atoms.mass.clear();
                staged->atoms.charge.clear();
                staged->residues = Residues{};
                staged->exclusions = Exclusions{};
                staged->generalized_born = GeneralizedBorn{};
                staged->virtual_atoms = VirtualAtoms{};
                Load_Reset_Classical_Force_Field(
                    &staged->classical_force_field);
            }
            if (has_rst7)
            {
                staged->atoms.coordinate.clear();
                staged->atoms.velocity.clear();
                staged->box = Box{};
                staged->start_time = 0.0;
            }
            staged->source = InputSource::kAmber;
            if (has_parm7)
            {
                document = Amber_Read_Parm7_Document(
                    controller->Original_Command("amber_parm7"), controller);
                pointers = Amber_Parse_Pointers(document, controller);
                Amber_Load_Parm7(staged, document, pointers, controller);
            }
            Amber_Load_Rst7(staged, has_parm7 ? &pointers : nullptr,
                            controller, &request_nopbc_default);
            if (has_parm7)
            {
                Amber_Load_Classical_Force_Field(staged, document, pointers,
                                                 controller);
            }
            if (request_nopbc_default)
            {
                // The surrounding load transaction stages controller command
                // maps too, so this inferred default and the System publish
                // either commit together or both roll back on an exception.
                controller->Set_Command("pbc", "false");
            }
        });
}

}  // namespace Xponge
