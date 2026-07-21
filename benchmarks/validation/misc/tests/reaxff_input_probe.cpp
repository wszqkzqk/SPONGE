#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "manybody/reaxff/reaxff_input.h"

bool Float_Memory_Is_Finite(const void* address)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, address, sizeof(bits));
    return (bits & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000);
}

bool Float_Memory_Is_Normal(const void* address)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, address, sizeof(bits));
    const std::uint32_t exponent = bits & UINT32_C(0x7f800000);
    return exponent != 0 && exponent != UINT32_C(0x7f800000);
}

bool Float_Memory_Is_Zero_Or_Normal(const void* address)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, address, sizeof(bits));
    const std::uint32_t magnitude = bits & UINT32_C(0x7fffffff);
    const std::uint32_t exponent = bits & UINT32_C(0x7f800000);
    return magnitude == 0 ||
           (exponent != 0 && exponent != UINT32_C(0x7f800000));
}

bool Double_Memory_Is_Finite(const void* address)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, address, sizeof(bits));
    return (bits & UINT64_C(0x7ff0000000000000)) !=
           UINT64_C(0x7ff0000000000000);
}

namespace
{

class STRING_SOURCE : public REAXFF_BYTE_SOURCE
{
   public:
    explicit STRING_SOURCE(std::string input)
        : input_(std::move(input)), fail_at_(std::string::npos)
    {
    }

    STRING_SOURCE(std::string input, std::size_t fail_at)
        : input_(std::move(input)), fail_at_(fail_at)
    {
    }

    REAXFF_BYTE_READ_STATUS Read_Byte(unsigned char* value) override
    {
        if (position_ == fail_at_) return REAXFF_BYTE_READ_ERROR;
        if (position_ == input_.size()) return REAXFF_BYTE_READ_EOF;
        *value = static_cast<unsigned char>(input_[position_++]);
        return REAXFF_BYTE_READ_OK;
    }

   private:
    std::string input_;
    std::size_t position_ = 0;
    std::size_t fail_at_ = std::string::npos;
};

bool Check(bool condition, const char* message)
{
    if (condition) return true;
    std::fprintf(stderr, "%s\n", message);
    return false;
}

std::string Join(const std::vector<std::string>& lines)
{
    std::string result;
    for (const std::string& line : lines)
    {
        result += line;
        result += '\n';
    }
    return result;
}

std::vector<std::string> Complete_Force_Field_Lines(
    std::size_t header_length = 6)
{
    return {
        std::string(header_length, 'H'),
        "2 ! general parameter count",
        "1.0 ! general 1",
        "2.0 # general 2",
        "1 ! atom type count",
        "atom header 1",
        "atom header 2",
        "atom header 3",
        "A 1 1 1 1 1 1 1 1 ! atom line 1",
        "1 1 1 1 1 1 1 1 ! atom line 2",
        "1 1 1 1 1 1 1 1 ! atom line 3",
        "1 1 1 1 1 1 1 1 ! atom line 4",
        "1 ! bond count",
        "bond header",
        "1 1 1 1 1 1 1 1 1 1 ! bond line 1",
        "1 1 1 1 1 1 1 1 ! bond line 2",
        "1 ! off-diagonal count",
        "1 1 1 1 1 1 1 1 ! off-diagonal entry",
        "1 ! angle count",
        "1 1 1 1 1 1 1 1 1 1 ! angle entry",
        "1 ! torsion count",
        "1 1 1 1 1 1 1 1 1 1 1 ! torsion entry",
        "1 ! hydrogen-bond count",
        "1 1 1 1 1 1 1 ! hydrogen-bond entry"};
}

std::vector<std::string> Many_Atom_Type_Force_Field_Lines(int atom_types)
{
    std::vector<std::string> lines = {
        "header", "0", std::to_string(atom_types), "atom header 1",
        "atom header 2", "atom header 3"};
    for (int atom = 0; atom < atom_types; atom++)
    {
        lines.push_back("A" + std::to_string(atom) +
                        " 1 1 1 1 1 1 1 1");
        lines.push_back("1 1 1 1 1 1 1 1");
        lines.push_back("1 1 1 1 1 1 1 1");
        lines.push_back("1 1 1 1 1 1 1 1");
    }
    lines.push_back("0");
    lines.push_back("bond header");
    lines.push_back("0");
    lines.push_back("0");
    lines.push_back("0");
    lines.push_back("0");
    return lines;
}

bool Parse_Force_Field(const std::vector<std::string>& lines,
                       REAXFF_FORCE_FIELD_IR* output,
                       REAXFF_INPUT_ERROR* error)
{
    STRING_SOURCE source(Join(lines));
    return ReaxFF_Parse_Force_Field(&source, "memory.ffield", output, error);
}

bool Test_Long_Physical_Line_And_Complete_IR()
{
    std::vector<std::string> lines = Complete_Force_Field_Lines(4097);
    lines[2] += std::string(4097, 'c');
    REAXFF_FORCE_FIELD_IR force_field;
    REAXFF_INPUT_ERROR error;
    if (!Check(Parse_Force_Field(lines, &force_field, &error),
               error.Describe().c_str()))
        return false;
    return Check(force_field.general_parameters.size() == 2,
                 "general parameters were not parsed") &&
           Check(force_field.atom_types.size() == 1,
                 "atom type was not parsed") &&
           Check(force_field.bonds.size() == 1,
                 "bond entry was not parsed") &&
           Check(force_field.off_diagonal.size() == 1,
                 "off-diagonal entry was not parsed") &&
           Check(force_field.angles.size() == 1,
                 "angle entry was not parsed") &&
           Check(force_field.torsions.size() == 1,
                 "torsion entry was not parsed") &&
           Check(force_field.hydrogen_bonds.size() == 1,
                 "hydrogen-bond entry was not parsed");
}

bool Test_Every_Record_Rejects_Trailing_Tokens_Transactionally()
{
    const std::size_t strict_record_lines[] = {
        1,  2,  4,  8,  9,  10, 11, 12, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 23};
    for (std::size_t record_line : strict_record_lines)
    {
        std::vector<std::string> lines = Complete_Force_Field_Lines();
        const std::size_t comment =
            lines[record_line].find_first_of("!#");
        if (comment != std::string::npos)
            lines[record_line].erase(comment);
        lines[record_line] += " trailing_garbage";
        REAXFF_FORCE_FIELD_IR force_field;
        force_field.general_parameters.push_back(123.0f);
        REAXFF_INPUT_ERROR error;
        if (!Check(!Parse_Force_Field(lines, &force_field, &error),
                   "record with a trailing token was accepted") ||
            !Check(error.kind == REAXFF_INPUT_FORMAT_ERROR,
                   "trailing token did not produce a format error") ||
            !Check(force_field.general_parameters.size() == 1 &&
                       force_field.general_parameters[0] == 123.0f,
                   "failed parse modified the caller's existing IR"))
            return false;
    }
    return true;
}

bool Test_Truncation_And_Strict_EOF()
{
    std::vector<std::string> truncated = Complete_Force_Field_Lines();
    truncated.pop_back();
    REAXFF_FORCE_FIELD_IR force_field;
    REAXFF_INPUT_ERROR error;
    if (!Check(!Parse_Force_Field(truncated, &force_field, &error),
               "truncated hydrogen-bond section was accepted") ||
        !Check(error.kind == REAXFF_INPUT_FORMAT_ERROR,
               "clean truncation was not classified as a format error"))
        return false;

    std::vector<std::string> trailing = Complete_Force_Field_Lines();
    trailing.push_back("unexpected data");
    if (!Check(!Parse_Force_Field(trailing, &force_field, &error),
               "unexpected data after the final section was accepted"))
        return false;

    trailing.back() = "  # permitted trailing comment";
    trailing.push_back("");
    return Check(Parse_Force_Field(trailing, &force_field, &error),
                 "blank/comment-only tail was rejected");
}

bool Test_Count_And_Multiplication_Limits()
{
    std::vector<std::string> huge_count = Complete_Force_Field_Lines();
    huge_count[1] = "2147483647";
    REAXFF_FORCE_FIELD_IR force_field;
    REAXFF_INPUT_ERROR error;
    if (!Check(!Parse_Force_Field(huge_count, &force_field, &error),
               "INT_MAX section count was accepted") ||
        !Check(error.kind == REAXFF_INPUT_RESOURCE_ERROR,
                   "INT_MAX count did not hit the explicit resource limit"))
        return false;

    if (!Check(!Parse_Force_Field(Many_Atom_Type_Force_Field_Lines(100),
                                  &force_field, &error),
               "large non-overflowing atom type table was accepted") ||
        !Check(error.kind == REAXFF_INPUT_RESOURCE_ERROR,
               "large dense table did not fail at the parser resource "
               "boundary"))
        return false;

    int result = 0;
    std::size_t bytes = 0;
    return Check(!ReaxFF_Checked_Int_Power(50000, 2, &result),
                 "integer square overflow was not rejected") &&
           Check(!ReaxFF_Checked_Int_Power(2000, 3, &result),
                 "integer cube overflow was not rejected") &&
           Check(!ReaxFF_Checked_Int_Power(300, 4, &result),
                 "integer fourth-power overflow was not rejected") &&
           Check(!ReaxFF_Checked_Dense_Table_Count(100, 4, 2 * sizeof(int),
                                                   &result),
                 "large non-overflowing dense table was not rejected") &&
           Check(!ReaxFF_Checked_Size_Multiply(
                     std::numeric_limits<std::size_t>::max(), 2, &bytes),
                 "size_t byte multiplication overflow was not rejected");
}

bool Test_IO_Error_Is_Not_Clean_EOF()
{
    const std::string complete = Join(Complete_Force_Field_Lines());
    STRING_SOURCE failing_source(complete, 11);
    REAXFF_FORCE_FIELD_IR force_field;
    REAXFF_INPUT_ERROR error;
    if (!Check(!ReaxFF_Parse_Force_Field(&failing_source, "injected.ffield",
                                         &force_field, &error),
               "injected I/O error was accepted") ||
        !Check(error.kind == REAXFF_INPUT_IO_ERROR,
               "injected read failure was not classified as I/O"))
        return false;

    STRING_SOURCE clean_eof("header\n");
    return Check(!ReaxFF_Parse_Force_Field(&clean_eof, "short.ffield",
                                           &force_field, &error),
                 "clean premature EOF was accepted") &&
           Check(error.kind == REAXFF_INPUT_FORMAT_ERROR,
                 "clean premature EOF was confused with an I/O error");
}

bool Test_Type_File_Strictness()
{
    REAXFF_FORCE_FIELD_IR force_field;
    REAXFF_INPUT_ERROR error;
    if (!Parse_Force_Field(Complete_Force_Field_Lines(), &force_field, &error))
        return Check(false, error.Describe().c_str());

    std::vector<int> atom_types;
    std::vector<int> hydrogen;
    STRING_SOURCE valid("1 ! atoms\nA # type\n");
    if (!Check(ReaxFF_Parse_Type_File(&valid, "type.txt", 1, force_field,
                                      &atom_types, &hydrogen, &error),
               error.Describe().c_str()) ||
        !Check(atom_types.size() == 1 && atom_types[0] == 0,
               "valid type file produced the wrong mapping"))
        return false;

    STRING_SOURCE trailing_field("1\nA extra\n");
    atom_types.assign(1, 37);
    hydrogen.assign(1, 41);
    if (!Check(!ReaxFF_Parse_Type_File(&trailing_field, "type.txt", 1,
                                       force_field, &atom_types, &hydrogen,
                                       &error),
               "type record with an extra field was accepted") ||
        !Check(atom_types.size() == 1 && atom_types[0] == 37 &&
                   hydrogen.size() == 1 && hydrogen[0] == 41,
               "failed type parse modified the caller's existing tables"))
        return false;

    STRING_SOURCE extra_record("1\nA\nA\n");
    if (!Check(!ReaxFF_Parse_Type_File(&extra_record, "type.txt", 1,
                                       force_field, &atom_types, NULL, &error),
               "extra type record was accepted"))
        return false;

    STRING_SOURCE truncated("2\nA\n");
    return Check(!ReaxFF_Parse_Type_File(&truncated, "type.txt", 2,
                                         force_field, &atom_types, NULL,
                                         &error),
                 "truncated type table was accepted");
}

}  // namespace

int main()
{
    if (!Test_Long_Physical_Line_And_Complete_IR() ||
        !Test_Every_Record_Rejects_Trailing_Tokens_Transactionally() ||
        !Test_Truncation_And_Strict_EOF() ||
        !Test_Count_And_Multiplication_Limits() ||
        !Test_IO_Error_Is_Not_Clean_EOF() || !Test_Type_File_Strictness())
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
