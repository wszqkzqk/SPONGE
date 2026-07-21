#pragma once

#include "../../common.h"

#include <array>
#include <errno.h>
#include <new>

// ReaxFF force-field records have a fixed whitespace-delimited prefix followed
// by an optional !/# comment.  The helpers below deliberately do not use scanf:
// scanf accepts partial tokens ("1junk"), and a fixed fgets buffer turns one
// physical line into several logical records.
class REAXFF_INPUT_LINE
{
   public:
    explicit REAXFF_INPUT_LINE(const char* line) : cursor_(line) {}
    explicit REAXFF_INPUT_LINE(const std::string& line)
        : cursor_(line.c_str())
    {
    }

    bool Read_Int(int* value)
    {
        Skip_Whitespace();
        if (value == NULL || cursor_ == NULL || cursor_[0] == '\0')
            return false;

        errno = 0;
        char* end = NULL;
        const long parsed = strtol(cursor_, &end, 10);
        if (end == cursor_ || !Is_Token_End(end) || errno == ERANGE ||
            parsed < INT_MIN || parsed > INT_MAX)
        {
            return false;
        }
        *value = static_cast<int>(parsed);
        cursor_ = end;
        return true;
    }

    bool Read_Float(float* value)
    {
        Skip_Whitespace();
        if (value == NULL || cursor_ == NULL || cursor_[0] == '\0')
            return false;

        errno = 0;
        char* end = NULL;
        const float parsed = strtof(cursor_, &end);
        if (end == cursor_ || !Is_Token_End(end) || errno == ERANGE ||
            !Float_Memory_Is_Finite(&parsed))
        {
            return false;
        }
        *value = parsed;
        cursor_ = end;
        return true;
    }

    bool Read_Word(std::string* value)
    {
        Skip_Whitespace();
        if (value == NULL || cursor_ == NULL || cursor_[0] == '\0')
            return false;

        const char* end = cursor_;
        while (end[0] != '\0' &&
               !isspace(static_cast<unsigned char>(end[0])))
        {
            end++;
        }
        if (end == cursor_) return false;
        value->assign(cursor_, static_cast<std::size_t>(end - cursor_));
        cursor_ = end;
        return true;
    }

    bool Read_Word(char* value, std::size_t capacity)
    {
        Skip_Whitespace();
        if (value == NULL || capacity == 0 || cursor_ == NULL ||
            cursor_[0] == '\0')
            return false;

        const char* end = cursor_;
        while (end[0] != '\0' &&
               !isspace(static_cast<unsigned char>(end[0])))
        {
            end++;
        }
        const std::size_t length =
            static_cast<std::size_t>(end - cursor_);
        if (length == 0 || length >= capacity) return false;
        memcpy(value, cursor_, length);
        value[length] = '\0';
        cursor_ = end;
        return true;
    }

    bool At_End_Or_Comment()
    {
        Skip_Whitespace();
        return cursor_ == NULL || cursor_[0] == '\0' || cursor_[0] == '!' ||
               cursor_[0] == '#';
    }

   private:
    const char* cursor_ = NULL;

    void Skip_Whitespace()
    {
        if (cursor_ == NULL) return;
        while (cursor_[0] != '\0' &&
               isspace(static_cast<unsigned char>(cursor_[0])))
        {
            cursor_++;
        }
    }

    static bool Is_Token_End(const char* end)
    {
        return end != NULL &&
               (end[0] == '\0' ||
                isspace(static_cast<unsigned char>(end[0])));
    }
};

inline bool ReaxFF_Checked_Size_Multiply(std::size_t lhs, std::size_t rhs,
                                         std::size_t* result)
{
    if (result == NULL) return false;
    if (lhs != 0 && rhs > std::numeric_limits<std::size_t>::max() / lhs)
        return false;
    *result = lhs * rhs;
    return true;
}

inline bool ReaxFF_Checked_Int_Power(int base, int exponent, int* result)
{
    if (base < 0 || exponent < 0 || result == NULL) return false;
    int product = 1;
    for (int i = 0; i < exponent; i++)
    {
        if (base != 0 && product > INT_MAX / base) return false;
        product *= base;
    }
    *result = product;
    return true;
}

// Dense three-/four-body lookup tables are indexed with signed ints in the
// kernels.  In addition to arithmetic overflow, reject a table allocation set
// whose declared extent would exceed this explicit resource contract.  Callers
// include staging/host/device copies in element_size, so a malicious but
// non-overflowing type count fails before any large allocation.
constexpr std::size_t REAXFF_MAX_DENSE_PARAMETER_BYTES =
    static_cast<std::size_t>(256) * 1024 * 1024;

inline bool ReaxFF_Checked_Dense_Table_Count(int base, int exponent,
                                             std::size_t element_size,
                                             int* result)
{
    int count = 0;
    std::size_t bytes = 0;
    if (!ReaxFF_Checked_Int_Power(base, exponent, &count) ||
        !ReaxFF_Checked_Size_Multiply(static_cast<std::size_t>(count),
                                      element_size, &bytes) ||
        bytes > REAXFF_MAX_DENSE_PARAMETER_BYTES)
    {
        return false;
    }
    if (result != NULL) *result = count;
    return true;
}

enum REAXFF_INPUT_ERROR_KIND
{
    REAXFF_INPUT_NO_ERROR = 0,
    REAXFF_INPUT_OPEN_ERROR,
    REAXFF_INPUT_IO_ERROR,
    REAXFF_INPUT_FORMAT_ERROR,
    REAXFF_INPUT_RESOURCE_ERROR
};

struct REAXFF_INPUT_ERROR
{
    REAXFF_INPUT_ERROR_KIND kind = REAXFF_INPUT_NO_ERROR;
    std::string file_name;
    std::size_t line_number = 0;
    std::string reason;

    std::string Describe() const
    {
        std::ostringstream message;
        message << reason;
        if (line_number != 0) message << " at line " << line_number;
        if (!file_name.empty()) message << " in file " << file_name;
        return message.str();
    }
};

inline bool ReaxFF_Set_Input_Error(REAXFF_INPUT_ERROR* error,
                                   REAXFF_INPUT_ERROR_KIND kind,
                                   const std::string& file_name,
                                   std::size_t line_number,
                                   const std::string& reason)
{
    if (error != NULL)
    {
        error->kind = kind;
        error->file_name = file_name;
        error->line_number = line_number;
        error->reason = reason;
    }
    return false;
}

enum REAXFF_BYTE_READ_STATUS
{
    REAXFF_BYTE_READ_OK = 0,
    REAXFF_BYTE_READ_EOF,
    REAXFF_BYTE_READ_ERROR
};

// The source abstraction is intentionally tiny.  Production uses FILE*, while
// parser tests can inject a deterministic read error without non-portable
// fopencookie/funopen APIs.
class REAXFF_BYTE_SOURCE
{
   public:
    virtual ~REAXFF_BYTE_SOURCE() {}
    virtual REAXFF_BYTE_READ_STATUS Read_Byte(unsigned char* value) = 0;
};

class REAXFF_FILE_BYTE_SOURCE : public REAXFF_BYTE_SOURCE
{
   public:
    explicit REAXFF_FILE_BYTE_SOURCE(FILE* file) : file_(file) {}

    REAXFF_BYTE_READ_STATUS Read_Byte(unsigned char* value) override
    {
        if (file_ == NULL || value == NULL) return REAXFF_BYTE_READ_ERROR;
        const int next = fgetc(file_);
        if (next != EOF)
        {
            *value = static_cast<unsigned char>(next);
            return REAXFF_BYTE_READ_OK;
        }
        return ferror(file_) ? REAXFF_BYTE_READ_ERROR : REAXFF_BYTE_READ_EOF;
    }

   private:
    FILE* file_ = NULL;
};

class REAXFF_SCOPED_INPUT_FILE
{
   public:
    explicit REAXFF_SCOPED_INPUT_FILE(FILE* file) : file_(file) {}
    REAXFF_SCOPED_INPUT_FILE(const REAXFF_SCOPED_INPUT_FILE&) = delete;
    REAXFF_SCOPED_INPUT_FILE& operator=(const REAXFF_SCOPED_INPUT_FILE&) =
        delete;
    ~REAXFF_SCOPED_INPUT_FILE()
    {
        if (file_ != NULL) fclose(file_);
    }

    FILE* Get() const { return file_; }

    int Close()
    {
        if (file_ == NULL) return 0;
        FILE* closing = file_;
        file_ = NULL;
        return fclose(closing);
    }

   private:
    FILE* file_ = NULL;
};

struct REAXFF_SOURCE_LINE
{
    std::size_t number = 0;
    std::string text;
};

constexpr std::size_t REAXFF_MAX_LOGICAL_LINE_BYTES =
    static_cast<std::size_t>(16) * 1024 * 1024;
constexpr std::size_t REAXFF_MAX_INPUT_FILE_BYTES =
    static_cast<std::size_t>(64) * 1024 * 1024;
constexpr std::size_t REAXFF_MAX_INPUT_LOGICAL_LINES = 2000000;
constexpr int REAXFF_MAX_FORCE_FIELD_SECTION_RECORDS = 1000000;

inline bool ReaxFF_Read_All_Lines(REAXFF_BYTE_SOURCE* source,
                                  const std::string& file_name,
                                  std::vector<REAXFF_SOURCE_LINE>* lines,
                                  REAXFF_INPUT_ERROR* error)
{
    if (source == NULL || lines == NULL)
    {
        return ReaxFF_Set_Input_Error(
            error, REAXFF_INPUT_IO_ERROR, file_name, 0,
            "received a null input source or output line collection");
    }

    std::vector<REAXFF_SOURCE_LINE> staged_lines;
    std::string current;
    std::size_t current_line = 1;
    std::size_t total_bytes = 0;
    while (true)
    {
        unsigned char byte = 0;
        const REAXFF_BYTE_READ_STATUS status = source->Read_Byte(&byte);
        if (status == REAXFF_BYTE_READ_ERROR)
        {
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_IO_ERROR, file_name, current_line,
                "I/O error while reading input");
        }
        if (status == REAXFF_BYTE_READ_EOF)
        {
            if (!current.empty())
            {
                if (staged_lines.size() == REAXFF_MAX_INPUT_LOGICAL_LINES)
                {
                    return ReaxFF_Set_Input_Error(
                        error, REAXFF_INPUT_RESOURCE_ERROR, file_name,
                        current_line,
                        "input exceeds the supported logical-line count");
                }
                staged_lines.push_back({current_line, std::move(current)});
            }
            break;
        }

        if (total_bytes == REAXFF_MAX_INPUT_FILE_BYTES)
        {
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_RESOURCE_ERROR, file_name, current_line,
                "input exceeds the supported 64 MiB file-size limit");
        }
        total_bytes++;
        if (byte == 0)
        {
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, current_line,
                "NUL byte is not valid in a text input file");
        }
        if (byte == '\n')
        {
            if (staged_lines.size() == REAXFF_MAX_INPUT_LOGICAL_LINES)
            {
                return ReaxFF_Set_Input_Error(
                    error, REAXFF_INPUT_RESOURCE_ERROR, file_name,
                    current_line,
                    "input exceeds the supported logical-line count");
            }
            staged_lines.push_back({current_line, std::move(current)});
            current.clear();
            current_line++;
            continue;
        }
        if (current.size() == REAXFF_MAX_LOGICAL_LINE_BYTES)
        {
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_RESOURCE_ERROR, file_name, current_line,
                "logical line exceeds the supported 16 MiB limit");
        }
        current.push_back(static_cast<char>(byte));
    }

    *lines = std::move(staged_lines);
    if (error != NULL) *error = REAXFF_INPUT_ERROR();
    return true;
}

class REAXFF_LINE_CURSOR
{
   public:
    REAXFF_LINE_CURSOR(const std::vector<REAXFF_SOURCE_LINE>& lines,
                       const std::string& file_name)
        : lines_(lines), file_name_(file_name)
    {
    }

    bool Next(const char* stage, const REAXFF_SOURCE_LINE** line,
              REAXFF_INPUT_ERROR* error)
    {
        if (line == NULL)
        {
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_IO_ERROR, file_name_, 0,
                "received a null logical-line output pointer");
        }
        if (index_ >= lines_.size())
        {
            std::ostringstream reason;
            reason << "unexpected end of file while reading "
                   << (stage == NULL ? "input record" : stage);
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name_,
                lines_.empty() ? 1 : lines_.back().number + 1, reason.str());
        }
        *line = &lines_[index_++];
        return true;
    }

    std::size_t Remaining() const { return lines_.size() - index_; }

    bool Require_Only_Comments_Or_Blanks(REAXFF_INPUT_ERROR* error)
    {
        while (index_ < lines_.size())
        {
            const REAXFF_SOURCE_LINE& line = lines_[index_++];
            REAXFF_INPUT_LINE fields(line.text);
            if (!fields.At_End_Or_Comment())
            {
                return ReaxFF_Set_Input_Error(
                    error, REAXFF_INPUT_FORMAT_ERROR, file_name_, line.number,
                    "unexpected trailing data after the hydrogen-bond section");
            }
        }
        return true;
    }

   private:
    const std::vector<REAXFF_SOURCE_LINE>& lines_;
    const std::string& file_name_;
    std::size_t index_ = 0;
};

struct REAXFF_ATOM_TYPE_IR
{
    std::string name;
    std::array<std::array<float, 8>, 4> values{};
};

struct REAXFF_BOND_IR
{
    int type1 = 0;
    int type2 = 0;
    std::array<float, 8> line1{};
    std::array<float, 8> line2{};
};

struct REAXFF_OFF_DIAGONAL_IR
{
    int type1 = 0;
    int type2 = 0;
    std::array<float, 6> values{};
};

struct REAXFF_ANGLE_IR
{
    std::array<int, 3> types{};
    std::array<float, 7> values{};
};

struct REAXFF_TORSION_IR
{
    std::array<int, 4> types{};
    std::array<float, 7> values{};
};

struct REAXFF_HYDROGEN_BOND_IR
{
    std::array<int, 3> types{};
    std::array<float, 4> values{};
};

struct REAXFF_FORCE_FIELD_IR
{
    std::vector<float> general_parameters;
    std::vector<REAXFF_ATOM_TYPE_IR> atom_types;
    std::map<std::string, int> type_map;
    std::vector<REAXFF_BOND_IR> bonds;
    std::vector<REAXFF_OFF_DIAGONAL_IR> off_diagonal;
    std::vector<REAXFF_ANGLE_IR> angles;
    std::vector<REAXFF_TORSION_IR> torsions;
    std::vector<REAXFF_HYDROGEN_BOND_IR> hydrogen_bonds;
};

template <std::size_t Count>
inline bool ReaxFF_Read_Floats(REAXFF_INPUT_LINE* fields,
                               std::array<float, Count>* values)
{
    if (fields == NULL || values == NULL) return false;
    for (std::size_t i = 0; i < Count; i++)
        if (!fields->Read_Float(&(*values)[i])) return false;
    return true;
}

template <std::size_t Count>
inline bool ReaxFF_Read_Ints(REAXFF_INPUT_LINE* fields,
                             std::array<int, Count>* values)
{
    if (fields == NULL || values == NULL) return false;
    for (std::size_t i = 0; i < Count; i++)
        if (!fields->Read_Int(&(*values)[i])) return false;
    return true;
}

inline bool ReaxFF_Parse_Count_Record(const REAXFF_SOURCE_LINE& line,
                                      const std::string& file_name,
                                      const char* description, int minimum,
                                      int* count, REAXFF_INPUT_ERROR* error)
{
    REAXFF_INPUT_LINE fields(line.text);
    if (!fields.Read_Int(count) || !fields.At_End_Or_Comment() ||
        *count < minimum)
    {
        std::ostringstream reason;
        reason << "failed to parse " << description;
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_FORMAT_ERROR,
                                      file_name, line.number, reason.str());
    }
    return true;
}

inline bool ReaxFF_Check_Section_Count(
    int count, int physical_lines_per_record, std::size_t remaining_lines,
    const std::string& file_name, const REAXFF_SOURCE_LINE& count_line,
    const char* description, REAXFF_INPUT_ERROR* error)
{
    if (count > REAXFF_MAX_FORCE_FIELD_SECTION_RECORDS)
    {
        std::ostringstream reason;
        reason << description << " count " << count
               << " exceeds the supported limit "
               << REAXFF_MAX_FORCE_FIELD_SECTION_RECORDS;
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_RESOURCE_ERROR,
                                      file_name, count_line.number,
                                      reason.str());
    }
    const std::size_t per_record =
        static_cast<std::size_t>(physical_lines_per_record);
    if (per_record == 0 ||
        static_cast<std::size_t>(count) > remaining_lines / per_record)
    {
        std::ostringstream reason;
        reason << description << " declares " << count
               << " records, but the file does not contain enough complete "
                  "physical lines";
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_FORMAT_ERROR,
                                      file_name, count_line.number,
                                      reason.str());
    }
    return true;
}

inline bool ReaxFF_Type_Index_Is_Valid(int type, int atom_type_count)
{
    return type >= 1 && type <= atom_type_count;
}

inline bool ReaxFF_Parse_Force_Field(REAXFF_BYTE_SOURCE* source,
                                     const std::string& file_name,
                                     REAXFF_FORCE_FIELD_IR* output,
                                     REAXFF_INPUT_ERROR* error)
{
    if (output == NULL)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_IO_ERROR, file_name,
                                      0, "received a null force-field IR");
    }

    std::vector<REAXFF_SOURCE_LINE> lines;
    if (!ReaxFF_Read_All_Lines(source, file_name, &lines, error)) return false;
    REAXFF_LINE_CURSOR cursor(lines, file_name);
    REAXFF_FORCE_FIELD_IR staged;
    const REAXFF_SOURCE_LINE* line = NULL;

    if (!cursor.Next("force-field header", &line, error) ||
        !cursor.Next("general parameter count", &line, error))
        return false;
    int general_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name,
                                   "number of general parameters", 0,
                                   &general_count, error) ||
        !ReaxFF_Check_Section_Count(general_count, 1, cursor.Remaining(),
                                    file_name, *line, "general parameter",
                                    error))
        return false;
    for (int i = 0; i < general_count; i++)
    {
        if (!cursor.Next("general parameter", &line, error)) return false;
        float value = 0.0f;
        REAXFF_INPUT_LINE fields(line->text);
        if (!fields.Read_Float(&value) || !fields.At_End_Or_Comment())
        {
            std::ostringstream reason;
            reason << "failed to parse general parameter at index " << i + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        staged.general_parameters.push_back(value);
    }

    if (!cursor.Next("atom type count", &line, error)) return false;
    int atom_type_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name,
                                   "number of atom types", 1,
                                   &atom_type_count, error) ||
        !ReaxFF_Check_Section_Count(atom_type_count, 4, cursor.Remaining(),
                                    file_name, *line, "atom type", error))
        return false;
    if (!ReaxFF_Checked_Dense_Table_Count(
            atom_type_count, 4, 3 * 2 * sizeof(int), NULL))
    {
        std::ostringstream reason;
        reason << "atom type count " << atom_type_count
               << " cannot be represented within the ReaxFF dense "
                  "four-body table resource limit";
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_RESOURCE_ERROR,
                                      file_name, line->number, reason.str());
    }
    for (int header = 0; header < 3; header++)
        if (!cursor.Next("atom type section header", &line, error))
            return false;

    for (int atom = 0; atom < atom_type_count; atom++)
    {
        REAXFF_ATOM_TYPE_IR entry;
        for (int record = 0; record < 4; record++)
        {
            if (!cursor.Next("atom type parameter record", &line, error))
                return false;
            REAXFF_INPUT_LINE fields(line->text);
            bool valid = true;
            if (record == 0) valid = fields.Read_Word(&entry.name);
            valid = valid && ReaxFF_Read_Floats(&fields, &entry.values[record]);
            if (!valid || !fields.At_End_Or_Comment())
            {
                std::ostringstream reason;
                reason << "failed to parse atom type block line " << record + 1
                       << " for type index " << atom + 1;
                return ReaxFF_Set_Input_Error(
                    error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                    reason.str());
            }
        }
        if (entry.name.empty() || staged.type_map.count(entry.name) != 0)
        {
            std::ostringstream reason;
            reason << "duplicate or empty atom type name '" << entry.name
                   << "' at type index " << atom + 1;
            return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_FORMAT_ERROR,
                                          file_name, line->number,
                                          reason.str());
        }
        staged.type_map[entry.name] = atom;
        staged.atom_types.push_back(std::move(entry));
    }

    if (!cursor.Next("bond parameter count", &line, error)) return false;
    int bond_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name,
                                   "number of bond parameters", 0,
                                   &bond_count, error) ||
        !ReaxFF_Check_Section_Count(bond_count, 2, cursor.Remaining(),
                                    file_name, *line, "bond parameter", error))
        return false;
    if (!cursor.Next("bond parameter section header", &line, error))
        return false;
    for (int bond = 0; bond < bond_count; bond++)
    {
        REAXFF_BOND_IR entry;
        if (!cursor.Next("bond parameter block line 1", &line, error))
            return false;
        REAXFF_INPUT_LINE line1(line->text);
        bool valid = line1.Read_Int(&entry.type1) &&
                     line1.Read_Int(&entry.type2) &&
                     ReaxFF_Read_Floats(&line1, &entry.line1) &&
                     line1.At_End_Or_Comment();
        if (!valid)
        {
            std::ostringstream reason;
            reason << "failed to parse bond parameter block line 1 at index "
                   << bond + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        if (!ReaxFF_Type_Index_Is_Valid(entry.type1, atom_type_count) ||
            !ReaxFF_Type_Index_Is_Valid(entry.type2, atom_type_count))
        {
            std::ostringstream reason;
            reason << "bond type index out of range at bond parameter index "
                   << bond + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        if (!cursor.Next("bond parameter block line 2", &line, error))
            return false;
        REAXFF_INPUT_LINE line2(line->text);
        if (!ReaxFF_Read_Floats(&line2, &entry.line2) ||
            !line2.At_End_Or_Comment())
        {
            std::ostringstream reason;
            reason << "failed to parse bond parameter block line 2 at index "
                   << bond + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        staged.bonds.push_back(entry);
    }

    if (!cursor.Next("off-diagonal parameter count", &line, error))
        return false;
    int off_diagonal_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name,
                                   "number of off-diagonal parameters", 0,
                                   &off_diagonal_count, error) ||
        !ReaxFF_Check_Section_Count(off_diagonal_count, 1, cursor.Remaining(),
                                    file_name, *line,
                                    "off-diagonal parameter", error))
        return false;
    for (int off = 0; off < off_diagonal_count; off++)
    {
        if (!cursor.Next("off-diagonal parameter entry", &line, error))
            return false;
        REAXFF_OFF_DIAGONAL_IR entry;
        REAXFF_INPUT_LINE fields(line->text);
        if (!fields.Read_Int(&entry.type1) || !fields.Read_Int(&entry.type2) ||
            !ReaxFF_Read_Floats(&fields, &entry.values) ||
            !fields.At_End_Or_Comment())
        {
            std::ostringstream reason;
            reason << "failed to parse off-diagonal parameter entry at index "
                   << off + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        if (!ReaxFF_Type_Index_Is_Valid(entry.type1, atom_type_count) ||
            !ReaxFF_Type_Index_Is_Valid(entry.type2, atom_type_count))
        {
            std::ostringstream reason;
            reason << "off-diagonal atom type index out of range at index "
                   << off + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        staged.off_diagonal.push_back(entry);
    }

    if (!cursor.Next("three-body parameter count", &line, error))
        return false;
    int angle_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name,
                                   "number of three-body parameters", 0,
                                   &angle_count, error) ||
        !ReaxFF_Check_Section_Count(angle_count, 1, cursor.Remaining(),
                                    file_name, *line, "three-body parameter",
                                    error))
        return false;
    for (int angle = 0; angle < angle_count; angle++)
    {
        if (!cursor.Next("three-body parameter entry", &line, error))
            return false;
        REAXFF_ANGLE_IR entry;
        REAXFF_INPUT_LINE fields(line->text);
        if (!ReaxFF_Read_Ints(&fields, &entry.types) ||
            !ReaxFF_Read_Floats(&fields, &entry.values) ||
            !fields.At_End_Or_Comment())
        {
            std::ostringstream reason;
            reason << "failed to parse three-body parameter entry at index "
                   << angle + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        for (int type : entry.types)
        {
            if (!ReaxFF_Type_Index_Is_Valid(type, atom_type_count))
            {
                std::ostringstream reason;
                reason << "three-body atom type index out of range at index "
                       << angle + 1;
                return ReaxFF_Set_Input_Error(
                    error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                    reason.str());
            }
        }
        staged.angles.push_back(entry);
    }

    if (!cursor.Next("torsion parameter count", &line, error)) return false;
    int torsion_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name,
                                   "number of torsion parameters", 0,
                                   &torsion_count, error) ||
        !ReaxFF_Check_Section_Count(torsion_count, 1, cursor.Remaining(),
                                    file_name, *line, "torsion parameter",
                                    error))
        return false;
    for (int torsion = 0; torsion < torsion_count; torsion++)
    {
        if (!cursor.Next("torsion parameter entry", &line, error))
            return false;
        REAXFF_TORSION_IR entry;
        REAXFF_INPUT_LINE fields(line->text);
        if (!ReaxFF_Read_Ints(&fields, &entry.types) ||
            !ReaxFF_Read_Floats(&fields, &entry.values) ||
            !fields.At_End_Or_Comment())
        {
            std::ostringstream reason;
            reason << "failed to parse torsion parameter entry at index "
                   << torsion + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        const bool explicit_outer =
            ReaxFF_Type_Index_Is_Valid(entry.types[0], atom_type_count) &&
            ReaxFF_Type_Index_Is_Valid(entry.types[3], atom_type_count);
        const bool wildcard_outer =
            entry.types[0] == 0 && entry.types[3] == 0;
        if (!ReaxFF_Type_Index_Is_Valid(entry.types[1], atom_type_count) ||
            !ReaxFF_Type_Index_Is_Valid(entry.types[2], atom_type_count) ||
            (!explicit_outer && !wildcard_outer))
        {
            std::ostringstream reason;
            reason << "torsion atom type index out of range at index "
                   << torsion + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        staged.torsions.push_back(entry);
    }

    if (!cursor.Next("hydrogen-bond parameter count", &line, error))
        return false;
    int hydrogen_bond_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name,
                                   "number of hydrogen bond parameters", 0,
                                   &hydrogen_bond_count, error) ||
        !ReaxFF_Check_Section_Count(hydrogen_bond_count, 1,
                                    cursor.Remaining(), file_name, *line,
                                    "hydrogen bond parameter", error))
        return false;
    for (int hydrogen = 0; hydrogen < hydrogen_bond_count; hydrogen++)
    {
        if (!cursor.Next("hydrogen bond parameter entry", &line, error))
            return false;
        REAXFF_HYDROGEN_BOND_IR entry;
        REAXFF_INPUT_LINE fields(line->text);
        if (!ReaxFF_Read_Ints(&fields, &entry.types) ||
            !ReaxFF_Read_Floats(&fields, &entry.values) ||
            !fields.At_End_Or_Comment())
        {
            std::ostringstream reason;
            reason << "failed to parse hydrogen bond parameter entry at index "
                   << hydrogen + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        for (int type : entry.types)
        {
            if (!ReaxFF_Type_Index_Is_Valid(type, atom_type_count))
            {
                std::ostringstream reason;
                reason << "hydrogen bond atom type index out of range at index "
                       << hydrogen + 1;
                return ReaxFF_Set_Input_Error(
                    error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                    reason.str());
            }
        }
        staged.hydrogen_bonds.push_back(entry);
    }

    if (!cursor.Require_Only_Comments_Or_Blanks(error)) return false;
    *output = std::move(staged);
    if (error != NULL) *error = REAXFF_INPUT_ERROR();
    return true;
}

inline bool ReaxFF_Parse_Force_Field_File(const char* file_name,
                                          REAXFF_FORCE_FIELD_IR* output,
                                          REAXFF_INPUT_ERROR* error)
{
    if (output == NULL)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_IO_ERROR,
                                      file_name == NULL ? "" : file_name, 0,
                                      "received a null force-field IR");
    }
    if (file_name == NULL)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_OPEN_ERROR, "", 0,
                                      "force-field path is null");
    }
    FILE* opened_file = fopen(file_name, "rb");
    if (opened_file == NULL)
    {
        std::ostringstream reason;
        reason << "failed to open input: " << strerror(errno);
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_OPEN_ERROR,
                                      file_name, 0, reason.str());
    }
    REAXFF_SCOPED_INPUT_FILE file(opened_file);
    REAXFF_FILE_BYTE_SOURCE source(file.Get());
    REAXFF_FORCE_FIELD_IR staged;
    bool parsed = false;
    try
    {
        parsed = ReaxFF_Parse_Force_Field(&source, file_name, &staged, error);
    }
    catch (const std::bad_alloc&)
    {
        parsed = ReaxFF_Set_Input_Error(
            error, REAXFF_INPUT_RESOURCE_ERROR, file_name, 0,
            "memory allocation failed while staging input");
    }
    const int close_status = file.Close();
    if (parsed && close_status != 0)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_IO_ERROR, file_name,
                                      0, "I/O error while closing input");
    }
    if (parsed) *output = std::move(staged);
    return parsed;
}

inline bool ReaxFF_Parse_Type_File(REAXFF_BYTE_SOURCE* source,
                                   const std::string& file_name,
                                   int expected_atom_count,
                                   const REAXFF_FORCE_FIELD_IR& force_field,
                                   std::vector<int>* atom_types,
                                   std::vector<int>* is_hydrogen,
                                   REAXFF_INPUT_ERROR* error)
{
    if (atom_types == NULL)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_IO_ERROR, file_name,
                                      0, "received a null atom-type table");
    }
    std::vector<REAXFF_SOURCE_LINE> lines;
    if (!ReaxFF_Read_All_Lines(source, file_name, &lines, error)) return false;
    REAXFF_LINE_CURSOR cursor(lines, file_name);
    const REAXFF_SOURCE_LINE* line = NULL;
    if (!cursor.Next("atom count", &line, error)) return false;
    int declared_atom_count = 0;
    if (!ReaxFF_Parse_Count_Record(*line, file_name, "atom numbers", 0,
                                   &declared_atom_count, error))
        return false;
    if (declared_atom_count != expected_atom_count)
    {
        std::ostringstream reason;
        reason << "atom numbers (" << declared_atom_count
               << ") does not match system (" << expected_atom_count << ')';
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_FORMAT_ERROR,
                                      file_name, line->number, reason.str());
    }
    if (static_cast<std::size_t>(declared_atom_count) > cursor.Remaining())
    {
        std::ostringstream reason;
        reason << "atom type table declares " << declared_atom_count
               << " entries, but the file contains only " << cursor.Remaining()
               << " remaining physical lines";
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_FORMAT_ERROR,
                                      file_name, line->number, reason.str());
    }

    std::vector<int> staged_atom_types;
    std::vector<int> staged_is_hydrogen;
    for (int atom = 0; atom < declared_atom_count; atom++)
    {
        if (!cursor.Next("atom type entry", &line, error)) return false;
        std::string name;
        REAXFF_INPUT_LINE fields(line->text);
        if (!fields.Read_Word(&name) || !fields.At_End_Or_Comment())
        {
            std::ostringstream reason;
            reason << "failed to parse atom type at index " << atom + 1;
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        const auto found = force_field.type_map.find(name);
        if (found == force_field.type_map.end())
        {
            std::ostringstream reason;
            reason << "atom type " << name
                   << " not found in the force-field atom type table";
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                reason.str());
        }
        staged_atom_types.push_back(found->second);
        if (is_hydrogen != NULL)
            staged_is_hydrogen.push_back(name == "H" ? 1 : 0);
    }
    // Type files do not have named trailing sections.  Only blank or comment
    // lines are permitted after the declared atom records.
    while (cursor.Remaining() != 0)
    {
        if (!cursor.Next("trailing type-file line", &line, error)) return false;
        REAXFF_INPUT_LINE fields(line->text);
        if (!fields.At_End_Or_Comment())
        {
            return ReaxFF_Set_Input_Error(
                error, REAXFF_INPUT_FORMAT_ERROR, file_name, line->number,
                "unexpected trailing data after the atom type table");
        }
    }

    *atom_types = std::move(staged_atom_types);
    if (is_hydrogen != NULL) *is_hydrogen = std::move(staged_is_hydrogen);
    if (error != NULL) *error = REAXFF_INPUT_ERROR();
    return true;
}

inline bool ReaxFF_Parse_Type_File_Path(
    const char* file_name, int expected_atom_count,
    const REAXFF_FORCE_FIELD_IR& force_field, std::vector<int>* atom_types,
    std::vector<int>* is_hydrogen, REAXFF_INPUT_ERROR* error)
{
    if (atom_types == NULL)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_IO_ERROR,
                                      file_name == NULL ? "" : file_name, 0,
                                      "received a null atom-type table");
    }
    if (file_name == NULL)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_OPEN_ERROR, "", 0,
                                      "atom-type file path is null");
    }
    FILE* opened_file = fopen(file_name, "rb");
    if (opened_file == NULL)
    {
        std::ostringstream reason;
        reason << "failed to open input: " << strerror(errno);
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_OPEN_ERROR,
                                      file_name, 0, reason.str());
    }
    REAXFF_SCOPED_INPUT_FILE file(opened_file);
    REAXFF_FILE_BYTE_SOURCE source(file.Get());
    std::vector<int> staged_atom_types;
    std::vector<int> staged_is_hydrogen;
    bool parsed = false;
    try
    {
        parsed = ReaxFF_Parse_Type_File(
            &source, file_name, expected_atom_count, force_field,
            &staged_atom_types,
            is_hydrogen == NULL ? NULL : &staged_is_hydrogen, error);
    }
    catch (const std::bad_alloc&)
    {
        parsed = ReaxFF_Set_Input_Error(
            error, REAXFF_INPUT_RESOURCE_ERROR, file_name, 0,
            "memory allocation failed while staging input");
    }
    const int close_status = file.Close();
    if (parsed && close_status != 0)
    {
        return ReaxFF_Set_Input_Error(error, REAXFF_INPUT_IO_ERROR, file_name,
                                      0, "I/O error while closing input");
    }
    if (parsed)
    {
        *atom_types = std::move(staged_atom_types);
        if (is_hydrogen != NULL)
            *is_hydrogen = std::move(staged_is_hydrogen);
    }
    return parsed;
}
