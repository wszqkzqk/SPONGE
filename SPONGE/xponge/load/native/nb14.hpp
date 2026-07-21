#pragma once

#include "md_core_parse.hpp"

namespace Xponge
{

static int Native_NB14_Read_Count(Native_Core_Parser* parser,
                                  std::size_t current_count,
                                  const std::string& field)
{
    const int count = parser->Read_Int(field);
    if (count < 0)
    {
        parser->Fail(spongeErrorBadFileFormat,
                     "negative " + field + " " + std::to_string(count));
    }
    const std::size_t count_size = static_cast<std::size_t>(count);
    const std::vector<int> integer_storage;
    const std::vector<float> real_storage;
    if (current_count >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        current_count > integer_storage.max_size() ||
        current_count > real_storage.max_size() ||
        count_size >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                current_count ||
        count_size > integer_storage.max_size() - current_count ||
        count_size > real_storage.max_size() - current_count)
    {
        parser->Fail(spongeErrorBadFileFormat,
                     field + " cannot be represented by the NB14 kernels or "
                             "their host containers");
    }
    return count;
}

static int Native_NB14_Read_Atom(Native_Core_Parser* parser,
                                 const std::string& field, int atom_numbers)
{
    const int atom = parser->Read_Int(field);
    if (atom < 0 || atom >= atom_numbers)
    {
        parser->Fail(spongeErrorBadFileFormat,
                     field + " " + std::to_string(atom) + " is outside [0, " +
                         std::to_string(atom_numbers) + ")");
    }
    return atom;
}

static void Native_NB14_Append_Record(
    NB14* parsed, std::set<std::pair<int, int>>* seen_pairs,
    Native_Core_Parser* parser, int atom_a, int atom_b, float coefficient_a,
    float coefficient_b, float charge_scale, const std::string& record)
{
    if (atom_a == atom_b)
    {
        parser->Fail(spongeErrorBadFileFormat,
                     record + " is a self 1-4 interaction for atom " +
                         std::to_string(atom_a));
    }
    const std::pair<int, int> pair =
        std::minmax(static_cast<int>(atom_a), static_cast<int>(atom_b));
    try
    {
        if (!seen_pairs->insert(pair).second)
        {
            parser->Fail(spongeErrorBadFileFormat,
                         record + " duplicates atom pair (" +
                             std::to_string(pair.first) + ", " +
                             std::to_string(pair.second) + ")");
        }
    }
    catch (const std::length_error&)
    {
        parser->Fail(spongeErrorBadFileFormat,
                     "NB14 pair identity storage exceeds its supported size");
    }
    catch (const std::bad_alloc&)
    {
        parser->Fail(spongeErrorMallocFailed,
                     "could not allocate transactional NB14 pair identity "
                     "storage");
    }

    parser->Append(&parsed->atom_a, atom_a, record + " atom A");
    parser->Append(&parsed->atom_b, atom_b, record + " atom B");
    parser->Append(&parsed->A, coefficient_a, record + " LJ A");
    parser->Append(&parsed->B, coefficient_b, record + " LJ B");
    parser->Append(&parsed->cf_scale_factor, charge_scale,
                   record + " charge scale");
}

static void Native_Load_NB14(NB14* nb14, const int* atom_type,
                             const float* pair_A, const float* pair_B,
                             int atom_numbers, int atom_type_numbers,
                             CONTROLLER* controller,
                             const char* module_name = "nb14")
{
    const bool has_derived =
        controller->Command_Exist(module_name, "in_file");
    const bool has_explicit =
        controller->Command_Exist(module_name, "extra_in_file");
    if (!has_derived && !has_explicit)
    {
        return;
    }
    if (atom_numbers < 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "Xponge::Native_Load_NB14",
            "Reason:\n\tthe runtime atom count for NB14 parsing is negative\n");
        return;
    }

    NB14 parsed;
    std::set<std::pair<int, int>> seen_pairs;

    // Explicit coefficients are materialized first to preserve the native
    // format's historical ordering, but neither input is published until both
    // files have passed strict parsing, EOF, and checked-close validation.
    if (has_explicit)
    {
        const std::string input_name =
            std::string(module_name) + "_extra_in_file";
        Native_Core_Parser parser(
            controller->Original_Command(module_name, "extra_in_file"),
            input_name.c_str(), "Xponge::Native_Load_NB14", controller);
        const int count = Native_NB14_Read_Count(
            &parser, parsed.atom_a.size(), "explicit interaction count");
        for (int i = 0; i < count; i++)
        {
            const std::string record =
                Native_Core_Entry_Field("explicit interaction", i);
            const int atom_a = Native_NB14_Read_Atom(
                &parser, record + " atom A", atom_numbers);
            const int atom_b = Native_NB14_Read_Atom(
                &parser, record + " atom B", atom_numbers);
            const double raw_a = parser.Read_Double(record + " LJ A");
            const double raw_b = parser.Read_Double(record + " LJ B");
            const float coefficient_a = parser.Checked_Float(
                raw_a * 12.0, record + " LJ A scaled by 12");
            const float coefficient_b = parser.Checked_Float(
                raw_b * 6.0, record + " LJ B scaled by 6");
            const float charge_scale =
                parser.Read_Float(record + " charge scale");
            Native_NB14_Append_Record(
                &parsed, &seen_pairs, &parser, atom_a, atom_b, coefficient_a,
                coefficient_b, charge_scale, record);
        }
        parser.Ensure_End();
        parser.Close();
    }

    if (has_derived)
    {
        const std::string input_name = std::string(module_name) + "_in_file";
        Native_Core_Parser parser(
            controller->Original_Command(module_name, "in_file"),
            input_name.c_str(), "Xponge::Native_Load_NB14", controller);
        const int count = Native_NB14_Read_Count(
            &parser, parsed.atom_a.size(), "derived interaction count");
        if (count > 0)
        {
            if (atom_type == NULL || pair_A == NULL || pair_B == NULL ||
                atom_type_numbers <= 0)
            {
                parser.Fail(
                    spongeErrorConflictingCommand,
                    "derived NB14 interactions require initialized LJ atom "
                    "types and triangular pair coefficients");
            }
            const std::uint64_t type_count =
                static_cast<std::uint64_t>(atom_type_numbers);
            const std::uint64_t pair_count =
                type_count * (type_count + UINT64_C(1)) / UINT64_C(2);
            if (pair_count >
                static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
            {
                parser.Fail(spongeErrorConflictingCommand,
                            "the initialized LJ atom-type count has an "
                            "unsupported triangular pair count");
            }
        }
        for (int i = 0; i < count; i++)
        {
            const std::string record =
                Native_Core_Entry_Field("derived interaction", i);
            const int atom_a = Native_NB14_Read_Atom(
                &parser, record + " atom A", atom_numbers);
            const int atom_b = Native_NB14_Read_Atom(
                &parser, record + " atom B", atom_numbers);
            const float lj_scale = parser.Read_Float(record + " LJ scale");
            const float charge_scale =
                parser.Read_Float(record + " charge scale");
            const int type_a = atom_type[atom_a];
            const int type_b = atom_type[atom_b];
            if (type_a < 0 || type_a >= atom_type_numbers || type_b < 0 ||
                type_b >= atom_type_numbers)
            {
                parser.Fail(spongeErrorConflictingCommand,
                            record + " references an atom whose initialized "
                                     "LJ type is outside [0, " +
                                std::to_string(atom_type_numbers) + ")");
            }
            const int small_type = std::min(type_a, type_b);
            const int large_type = std::max(type_a, type_b);
            const std::uint64_t pair_index =
                static_cast<std::uint64_t>(large_type) *
                    static_cast<std::uint64_t>(large_type + 1) /
                    UINT64_C(2) +
                static_cast<std::uint64_t>(small_type);
            const float source_a = pair_A[pair_index];
            const float source_b = pair_B[pair_index];
            if (!Float_Memory_Is_Finite(&source_a) ||
                !Float_Memory_Is_Zero_Or_Normal(&source_a) ||
                !Float_Memory_Is_Finite(&source_b) ||
                !Float_Memory_Is_Zero_Or_Normal(&source_b))
            {
                parser.Fail(spongeErrorConflictingCommand,
                            record + " references a non-finite or subnormal "
                                     "initialized LJ coefficient");
            }
            const float coefficient_a = parser.Checked_Float(
                static_cast<double>(lj_scale) *
                    static_cast<double>(source_a),
                record + " scaled LJ A");
            const float coefficient_b = parser.Checked_Float(
                static_cast<double>(lj_scale) *
                    static_cast<double>(source_b),
                record + " scaled LJ B");
            Native_NB14_Append_Record(
                &parsed, &seen_pairs, &parser, atom_a, atom_b, coefficient_a,
                coefficient_b, charge_scale, record);
        }
        parser.Ensure_End();
        parser.Close();
    }

    *nb14 = std::move(parsed);
}

static void Native_Load_NB14(System* system, CONTROLLER* controller)
{
    const LennardJones& lj = system->classical_force_field.lj;
    const int atom_numbers = Load_Get_Atom_Numbers(system);
    if (!lj.atom_type.empty() &&
        lj.atom_type.size() != static_cast<std::size_t>(atom_numbers))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "Xponge::Native_Load_NB14",
            "Reason:\n\tthe initialized LJ atom-type array and runtime atom "
            "count differ\n");
        return;
    }
    if (controller->Command_Exist("nb14", "in_file"))
    {
        if (lj.atom_type_numbers < 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorConflictingCommand, "Xponge::Native_Load_NB14",
                "Reason:\n\tthe initialized LJ atom-type count is negative\n");
            return;
        }
        const std::uint64_t type_count =
            static_cast<std::uint64_t>(lj.atom_type_numbers);
        const std::uint64_t pair_count =
            type_count * (type_count + UINT64_C(1)) / UINT64_C(2);
        if (pair_count >
                static_cast<std::uint64_t>(
                    std::numeric_limits<std::size_t>::max()) ||
            lj.pair_A.size() != static_cast<std::size_t>(pair_count) ||
            lj.pair_B.size() != static_cast<std::size_t>(pair_count))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorConflictingCommand, "Xponge::Native_Load_NB14",
                "Reason:\n\tthe initialized LJ triangular coefficient "
                "arrays do not match the LJ atom-type count\n");
            return;
        }
    }
    const int* atom_type = lj.atom_type.empty() ? NULL : lj.atom_type.data();
    const float* pair_A = lj.pair_A.empty() ? NULL : lj.pair_A.data();
    const float* pair_B = lj.pair_B.empty() ? NULL : lj.pair_B.data();
    Native_Load_NB14(&system->classical_force_field.nb14, atom_type, pair_A,
                     pair_B, atom_numbers, lj.atom_type_numbers, controller);
}

}  // namespace Xponge
