#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <unordered_set>

#include "../ir/virtual_atoms.hpp"
#include "../xponge.h"
#include "./common.hpp"

namespace Xponge
{

namespace fs = std::filesystem;

static constexpr float Gromacs_Pi = 3.14159265358979323846f;

struct Gromacs_Source_Reference
{
    std::uint32_t file_id = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t line_number = 0;
};

struct Gromacs_Defaults
{
    int nbfunc = 1;
    int comb_rule = 2;
    bool gen_pairs = false;
    float fudge_lj = 1.0f;
    float fudge_qq = 1.0f;
};

struct Gromacs_Atom_Type
{
    std::string name;
    std::string bond_type;
    float mass = 0.0f;
    float charge = 0.0f;
    std::string ptype;
    float v = 0.0f;
    float w = 0.0f;
    Gromacs_Source_Reference source;
};

struct Gromacs_Bond_Type
{
    std::string ai;
    std::string aj;
    int funct = 0;
    float b0 = 0.0f;
    float kb = 0.0f;
    Gromacs_Source_Reference source;
};

struct Gromacs_Angle_Type
{
    std::string ai;
    std::string aj;
    std::string ak;
    int funct = 0;
    float theta0 = 0.0f;
    float k = 0.0f;
    float ub0 = 0.0f;
    float kub = 0.0f;
    Gromacs_Source_Reference source;
};

struct Gromacs_Dihedral_Type
{
    std::string ai;
    std::string aj;
    std::string ak;
    std::string al;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_Pair_Type
{
    std::string ai;
    std::string aj;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_CMap_Type
{
    std::string ai;
    std::string aj;
    std::string ak;
    std::string al;
    std::string am;
    int funct = 0;
    int resolution = 0;
    std::vector<float> grid;
    Gromacs_Source_Reference source;
};

struct Gromacs_Molecule_Atom
{
    int nr = 0;
    std::string type;
    int resnr = 0;
    std::string residue;
    std::string atom;
    int cgnr = 0;
    float charge = 0.0f;
    float mass = 0.0f;
    bool has_charge = false;
    bool has_mass = false;
    Gromacs_Source_Reference source;
};

struct Gromacs_Bond
{
    int ai = 0;
    int aj = 0;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_Pair
{
    int ai = 0;
    int aj = 0;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_Angle
{
    int ai = 0;
    int aj = 0;
    int ak = 0;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_Dihedral
{
    int ai = 0;
    int aj = 0;
    int ak = 0;
    int al = 0;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_Settle
{
    int ow = 0;
    int funct = 0;
    float doh = 0.0f;
    float dhh = 0.0f;
    Gromacs_Source_Reference source;
};

struct Gromacs_Constraint_Type
{
    std::string ai;
    std::string aj;
    int funct = 0;
    float distance = 0.0f;
    Gromacs_Source_Reference source;
};

struct Gromacs_Constraint
{
    int ai = 0;
    int aj = 0;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_CMap
{
    int ai = 0;
    int aj = 0;
    int ak = 0;
    int al = 0;
    int am = 0;
    int funct = 0;
    Gromacs_Source_Reference source;
};

struct Gromacs_Virtual_Site
{
    int site = 0;
    std::vector<int> from;
    int funct = 0;
    std::vector<float> parameters;
    Gromacs_Source_Reference source;
};

struct Gromacs_Molecule
{
    std::string name;
    int nrexcl = 0;
    std::vector<Gromacs_Molecule_Atom> atoms;
    std::vector<Gromacs_Bond> bonds;
    std::vector<Gromacs_Pair> pairs;
    std::vector<Gromacs_Angle> angles;
    std::vector<Gromacs_Dihedral> dihedrals;
    std::vector<Gromacs_Settle> settles;
    std::vector<Gromacs_Constraint> constraints;
    std::vector<std::pair<int, int>> exclusions;
    std::vector<Gromacs_CMap> cmaps;
    std::vector<Gromacs_Virtual_Site> virtual_sites;
};

struct Gromacs_Residue_Info
{
    int atom_numbers = 0;
};

struct Gromacs_System_Molecule
{
    std::string name;
    int count = 0;
    Gromacs_Source_Reference source;
};

template <std::size_t N>
struct Gromacs_Type_Key
{
    std::array<std::string, N> types;
    int funct = 0;

    bool operator==(const Gromacs_Type_Key& other) const
    {
        return funct == other.funct && types == other.types;
    }
};

template <std::size_t N>
struct Gromacs_Type_Key_Hash
{
    std::size_t operator()(const Gromacs_Type_Key<N>& key) const
    {
        std::size_t result = std::hash<int>{}(key.funct);
        for (const std::string& type : key.types)
        {
            std::size_t value = std::hash<std::string>{}(type);
            result ^= value + 0x9e3779b9U + (result << 6) + (result >> 2);
        }
        return result;
    }
};

using Gromacs_Pair_Type_Key = Gromacs_Type_Key<2>;
using Gromacs_Triple_Type_Key = Gromacs_Type_Key<3>;
using Gromacs_Quadruple_Type_Key = Gromacs_Type_Key<4>;
using Gromacs_Quintuple_Type_Key = Gromacs_Type_Key<5>;

template <std::size_t N>
static Gromacs_Type_Key<N> Gromacs_Make_Type_Key(
    std::array<std::string, N> types, int funct, bool reversible)
{
    if (reversible)
    {
        std::array<std::string, N> reversed = types;
        std::reverse(reversed.begin(), reversed.end());
        if (reversed < types)
        {
            types = std::move(reversed);
        }
    }
    return {std::move(types), funct};
}

struct Gromacs_Topology
{
    Gromacs_Defaults defaults;
    std::unordered_map<std::string, Gromacs_Atom_Type> atom_types;
    std::vector<Gromacs_Bond_Type> bond_types;
    std::vector<Gromacs_Angle_Type> angle_types;
    std::vector<Gromacs_Dihedral_Type> dihedral_types;
    std::vector<Gromacs_Pair_Type> pair_types;
    std::vector<Gromacs_Pair_Type> nonbond_params;
    std::vector<Gromacs_Constraint_Type> constraint_types;
    std::vector<Gromacs_CMap_Type> cmap_types;
    std::unordered_map<Gromacs_Pair_Type_Key, std::size_t,
                       Gromacs_Type_Key_Hash<2>>
        bond_type_index;
    std::unordered_map<Gromacs_Pair_Type_Key, std::size_t,
                       Gromacs_Type_Key_Hash<2>>
        constraint_type_index;
    std::unordered_map<Gromacs_Triple_Type_Key, std::size_t,
                       Gromacs_Type_Key_Hash<3>>
        angle_type_index;
    std::unordered_map<Gromacs_Quadruple_Type_Key, std::vector<std::size_t>,
                       Gromacs_Type_Key_Hash<4>>
        dihedral_type_index;
    std::unordered_map<Gromacs_Pair_Type_Key, std::size_t,
                       Gromacs_Type_Key_Hash<2>>
        pair_type_index;
    std::unordered_map<Gromacs_Pair_Type_Key, std::size_t,
                       Gromacs_Type_Key_Hash<2>>
        nonbond_parameter_index;
    std::unordered_map<Gromacs_Quintuple_Type_Key, std::size_t,
                       Gromacs_Type_Key_Hash<5>>
        cmap_type_index;
    std::vector<fs::path> source_files;
    std::unordered_map<std::string, Gromacs_Molecule> molecules;
    std::vector<Gromacs_System_Molecule> system_molecules;
};

struct Gromacs_Source_Line
{
    std::string text;
    fs::path file_path;
    std::size_t line_number = 0;
};

using Gromacs_Macros =
    std::unordered_map<std::string, std::vector<std::string>>;

static std::string Gromacs_Trim(const std::string& value)
{
    std::size_t begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos)
    {
        return "";
    }
    std::size_t end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

static std::string Gromacs_Strip_Comment(const std::string& line)
{
    std::string trimmed = Gromacs_Trim(line);
    // Legacy GROMACS force fields (including the distributed CHARMM port)
    // use leading '*' banner lines as comments outside topology sections.
    if (!trimmed.empty() && trimmed.front() == '*')
    {
        return "";
    }
    std::size_t comment = line.find(';');
    if (comment == std::string::npos)
    {
        return trimmed;
    }
    return Gromacs_Trim(line.substr(0, comment));
}

static std::vector<std::string> Gromacs_Split(const std::string& line)
{
    std::vector<std::string> tokens;
    std::istringstream stream(line);
    std::string token;
    while (stream >> token)
    {
        tokens.push_back(token);
    }
    return tokens;
}

static std::vector<std::string> Gromacs_Split_List(const std::string& value)
{
    std::vector<std::string> tokens;
    std::string current;
    for (char ch : value)
    {
        if (ch == ',' || ch == ';' || ch == ':' || ch == ' ' || ch == '\t')
        {
            if (!current.empty())
            {
                tokens.push_back(current);
                current.clear();
            }
        }
        else
        {
            current.push_back(ch);
        }
    }
    if (!current.empty())
    {
        tokens.push_back(current);
    }
    return tokens;
}

static bool Gromacs_Is_Macro_Name(const std::string& value)
{
    if (value.empty() ||
        !((value[0] >= 'A' && value[0] <= 'Z') ||
          (value[0] >= 'a' && value[0] <= 'z') || value[0] == '_'))
    {
        return false;
    }
    for (char ch : value)
    {
        if (!((ch >= 'A' && ch <= 'Z') || (ch >= 'a' && ch <= 'z') ||
              (ch >= '0' && ch <= '9') || ch == '_'))
        {
            return false;
        }
    }
    return true;
}

static bool Gromacs_Is_Integer_Token(const std::string& value)
{
    if (value.empty())
    {
        return false;
    }
    std::size_t i = (value[0] == '+' || value[0] == '-') ? 1 : 0;
    if (i == value.size())
    {
        return false;
    }
    for (; i < value.size(); i++)
    {
        if (value[i] < '0' || value[i] > '9')
        {
            return false;
        }
    }
    return true;
}

enum class Gromacs_Macro_Expansion_Status
{
    kSuccess,
    kCycle,
    kCapacityError,
};

// Expand object-like macros in two iterative passes.  The first pass performs
// a three-colour DFS and memoizes only the expanded token/byte counts, so
// shared or empty expansion graphs are evaluated once without retaining a
// full token vector for every alias.  Exact macro aliases are also collapsed.
// The second pass materializes only the final output.  There is no policy
// token limit: failure is tied only to the real representational/allocation
// capacity of the host.
static Gromacs_Macro_Expansion_Status Gromacs_Expand_Macro_Tokens(
    const std::vector<std::string>& input, const Gromacs_Macros& macros,
    std::vector<std::string>* output, const std::string** cycle_macro)
{
    enum class Expansion_State
    {
        kVisiting,
        kComplete,
    };
    struct Expansion_Info
    {
        Expansion_State state = Expansion_State::kVisiting;
        std::size_t token_count = 0;
        std::size_t byte_count = 0;
        const std::string* materialization_macro = NULL;
    };
    struct Count_Frame
    {
        const std::string* macro_name = NULL;
        const std::vector<std::string>* replacement = NULL;
        std::size_t next_token = 0;
        std::size_t token_count = 0;
        std::size_t byte_count = 0;
        bool exact_alias = true;
        const std::string* alias_target = NULL;
    };
    struct Materialization_Frame
    {
        const std::vector<std::string>* replacement = NULL;
        std::size_t next_token = 0;
    };

    output->clear();
    *cycle_macro = NULL;
    std::vector<Count_Frame> count_frames;
    std::unordered_map<const std::string*, Expansion_Info> expansion_info;
    try
    {
        const std::size_t maximum_token_count = output->max_size();
        const std::size_t maximum_byte_count = std::string().max_size();
        auto add_counts = [&](Count_Frame* frame, std::size_t token_count,
                              std::size_t byte_count)
        {
            if (token_count > maximum_token_count - frame->token_count ||
                byte_count > maximum_byte_count - frame->byte_count)
            {
                return false;
            }
            frame->token_count += token_count;
            frame->byte_count += byte_count;
            return true;
        };

        for (const std::string& root_token : input)
        {
            const auto root_macro = macros.find(root_token);
            if (root_macro == macros.end())
            {
                continue;
            }

            const std::string* root_name = &root_macro->first;
            auto root_info = expansion_info.find(root_name);
            if (root_info == expansion_info.end())
            {
                if (expansion_info.size() == expansion_info.max_size() ||
                    count_frames.size() == count_frames.max_size())
                {
                    return Gromacs_Macro_Expansion_Status::kCapacityError;
                }
                expansion_info.emplace(root_name, Expansion_Info{});
                count_frames.push_back(
                    {root_name, &root_macro->second, 0, 0, 0, true, NULL});
            }

            while (!count_frames.empty())
            {
                Count_Frame& frame = count_frames.back();
                if (frame.next_token == frame.replacement->size())
                {
                    const std::string* completed_name = frame.macro_name;
                    const std::size_t completed_tokens = frame.token_count;
                    const std::size_t completed_bytes = frame.byte_count;
                    const std::string* materialization_macro =
                        frame.exact_alias && frame.alias_target != NULL
                            ? frame.alias_target
                            : completed_name;
                    count_frames.pop_back();
                    auto completed = expansion_info.find(completed_name);
                    completed->second.state = Expansion_State::kComplete;
                    completed->second.token_count = completed_tokens;
                    completed->second.byte_count = completed_bytes;
                    completed->second.materialization_macro =
                        materialization_macro;
                    continue;
                }

                const std::string& token =
                    (*frame.replacement)[frame.next_token];
                const auto dependency = macros.find(token);
                if (dependency == macros.end())
                {
                    frame.next_token++;
                    frame.exact_alias = false;
                    if (!add_counts(&frame, 1, token.size()))
                    {
                        return Gromacs_Macro_Expansion_Status::kCapacityError;
                    }
                    continue;
                }

                const std::string* dependency_name = &dependency->first;
                auto dependency_info = expansion_info.find(dependency_name);
                if (dependency_info == expansion_info.end())
                {
                    if (expansion_info.size() == expansion_info.max_size() ||
                        count_frames.size() == count_frames.max_size())
                    {
                        return Gromacs_Macro_Expansion_Status::kCapacityError;
                    }
                    expansion_info.emplace(dependency_name, Expansion_Info{});
                    count_frames.push_back({dependency_name,
                                            &dependency->second, 0, 0, 0, true,
                                            NULL});
                    continue;
                }
                if (dependency_info->second.state ==
                    Expansion_State::kVisiting)
                {
                    *cycle_macro = dependency_name;
                    return Gromacs_Macro_Expansion_Status::kCycle;
                }

                frame.next_token++;
                if (dependency_info->second.token_count > 0)
                {
                    if (frame.exact_alias && frame.alias_target == NULL)
                    {
                        frame.alias_target = dependency_info->second.
                            materialization_macro;
                    }
                    else
                    {
                        frame.exact_alias = false;
                    }
                }
                if (!add_counts(&frame, dependency_info->second.token_count,
                                dependency_info->second.byte_count))
                {
                    return Gromacs_Macro_Expansion_Status::kCapacityError;
                }
            }
        }

        std::size_t total_token_count = 0;
        std::size_t total_byte_count = 0;
        for (const std::string& token : input)
        {
            const auto macro = macros.find(token);
            const std::size_t added_tokens =
                macro == macros.end()
                    ? 1
                    : expansion_info.find(&macro->first)->second.token_count;
            const std::size_t added_bytes =
                macro == macros.end()
                    ? token.size()
                    : expansion_info.find(&macro->first)->second.byte_count;
            if (added_tokens > maximum_token_count - total_token_count ||
                added_bytes > maximum_byte_count - total_byte_count)
            {
                return Gromacs_Macro_Expansion_Status::kCapacityError;
            }
            total_token_count += added_tokens;
            total_byte_count += added_bytes;
        }
        if (total_token_count > 0 &&
            total_token_count - 1 > maximum_byte_count - total_byte_count)
        {
            return Gromacs_Macro_Expansion_Status::kCapacityError;
        }

        output->reserve(total_token_count);
        std::vector<Materialization_Frame> materialization_frames;
        for (const std::string& root_token : input)
        {
            const std::string* current_token = &root_token;
            while (current_token != NULL)
            {
                const auto macro = macros.find(*current_token);
                if (macro == macros.end())
                {
                    output->push_back(*current_token);
                }
                else
                {
                    const Expansion_Info& info =
                        expansion_info.find(&macro->first)->second;
                    if (info.token_count > 0)
                    {
                        const auto materialization =
                            macros.find(*info.materialization_macro);
                        if (materialization_frames.size() ==
                            materialization_frames.max_size())
                        {
                            return Gromacs_Macro_Expansion_Status::
                                kCapacityError;
                        }
                        materialization_frames.push_back(
                            {&materialization->second, 0});
                    }
                }

                current_token = NULL;
                while (!materialization_frames.empty())
                {
                    Materialization_Frame& frame =
                        materialization_frames.back();
                    if (frame.next_token < frame.replacement->size())
                    {
                        current_token =
                            &(*frame.replacement)[frame.next_token++];
                        break;
                    }
                    materialization_frames.pop_back();
                }
            }
        }
        if (output->size() != total_token_count)
        {
            return Gromacs_Macro_Expansion_Status::kCapacityError;
        }
    }
    catch (const std::bad_alloc&)
    {
        return Gromacs_Macro_Expansion_Status::kCapacityError;
    }
    catch (const std::length_error&)
    {
        return Gromacs_Macro_Expansion_Status::kCapacityError;
    }
    return Gromacs_Macro_Expansion_Status::kSuccess;
}

static bool Gromacs_Join_Expanded_Tokens(
    const std::vector<std::string>& tokens, std::string* expanded_line)
{
    try
    {
        std::size_t required_size = 0;
        const std::size_t maximum_size = expanded_line->max_size();
        for (std::size_t i = 0; i < tokens.size(); i++)
        {
            if (i > 0)
            {
                if (required_size == maximum_size)
                {
                    return false;
                }
                required_size++;
            }
            if (tokens[i].size() > maximum_size - required_size)
            {
                return false;
            }
            required_size += tokens[i].size();
        }

        expanded_line->clear();
        expanded_line->reserve(required_size);
        for (std::size_t i = 0; i < tokens.size(); i++)
        {
            if (i > 0)
            {
                expanded_line->push_back(' ');
            }
            expanded_line->append(tokens[i]);
        }
    }
    catch (const std::bad_alloc&)
    {
        return false;
    }
    catch (const std::length_error&)
    {
        return false;
    }
    return true;
}

static bool Gromacs_Is_True(const std::string& value)
{
    return value == "yes" || value == "Yes" || value == "YES" || value == "1" ||
           value == "true" || value == "TRUE";
}

static bool Gromacs_Is_False(const std::string& value)
{
    return value == "no" || value == "No" || value == "NO" || value == "0" ||
           value == "false" || value == "FALSE";
}

static fs::path Gromacs_Resolve_Include(
    const fs::path& parent_dir, const std::string& include_name,
    const std::vector<fs::path>& include_dirs, CONTROLLER* controller,
    const char* error_by, const fs::path& source_file, std::size_t source_line)
{
    fs::path candidate = parent_dir / include_name;
    if (fs::exists(candidate))
    {
        return candidate;
    }
    for (const fs::path& include_dir : include_dirs)
    {
        candidate = include_dir / include_name;
        if (fs::exists(candidate))
        {
            return candidate;
        }
    }
    std::string reason = "Reason:\n\tfailed to resolve GROMACS include file '" +
                         include_name + "' at " + source_file.string() + ":" +
                         std::to_string(source_line) + "\n";
    controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                   reason.c_str());
    return {};
}

static void Gromacs_Preprocess_File(const fs::path& file_path,
                                    Gromacs_Macros* macros,
                                    const std::vector<fs::path>& include_dirs,
                                    std::vector<Gromacs_Source_Line>* lines,
                                    std::vector<fs::path>* include_stack,
                                    CONTROLLER* controller,
                                    const char* error_by)
{
    std::error_code canonical_error;
    fs::path canonical_path = fs::weakly_canonical(file_path, canonical_error);
    if (canonical_error)
    {
        canonical_path = fs::absolute(file_path).lexically_normal();
    }
    if (std::find(include_stack->begin(), include_stack->end(),
                  canonical_path) != include_stack->end())
    {
        std::string reason =
            "Reason:\n\tGROMACS topology include cycle detected at '" +
            canonical_path.string() + "'\n";
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                       reason.c_str());
    }
    include_stack->push_back(canonical_path);

    std::ifstream fin(file_path);
    if (!fin.is_open())
    {
        std::string reason =
            "Reason:\n\tfailed to open GROMACS topology file '" +
            file_path.string() + "'\n";
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                       reason.c_str());
    }

    struct Conditional_State
    {
        bool parent_active = true;
        bool branch_active = true;
        bool else_seen = false;
    };
    std::vector<Conditional_State> stack;
    auto is_active = [&stack]() -> bool
    {
        if (stack.empty())
        {
            return true;
        }
        return stack.back().parent_active && stack.back().branch_active;
    };

    std::string raw_line;
    std::size_t line_number = 0;
    while (std::getline(fin, raw_line))
    {
        line_number++;
        std::size_t logical_line_number = line_number;
        std::string line = Gromacs_Strip_Comment(Gromacs_Trim(raw_line));
        while (!line.empty() && line.back() == '\\')
        {
            line.pop_back();
            line = Gromacs_Trim(line);
            if (!std::getline(fin, raw_line))
            {
                std::string reason =
                    "Reason:\n\tunterminated GROMACS line continuation at " +
                    file_path.string() + ":" +
                    std::to_string(logical_line_number) + "\n";
                controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                               error_by, reason.c_str());
            }
            line_number++;
            std::string next_line =
                Gromacs_Strip_Comment(Gromacs_Trim(raw_line));
            if (!next_line.empty())
            {
                if (!line.empty())
                {
                    line += " ";
                }
                line += next_line;
            }
        }
        if (line.empty())
        {
            continue;
        }
        auto throw_preprocessor_storage_error = [&]()
        {
            std::string message =
                "Reason:\n\tGROMACS preprocessor data exceeds available "
                "host storage at " +
                file_path.string() + ":" +
                std::to_string(logical_line_number) + "\n";
            controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                           message.c_str());
        };
        if (!line.empty() && line[0] == '#')
        {
            auto throw_directive_error = [&](const std::string& reason)
            {
                std::string message =
                    "Reason:\n\t" + reason + " at " + file_path.string() + ":" +
                    std::to_string(logical_line_number) + "\n";
                controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                               error_by, message.c_str());
            };
            std::vector<std::string> tokens;
            try
            {
                tokens = Gromacs_Split(line);
            }
            catch (const std::bad_alloc&)
            {
                throw_preprocessor_storage_error();
            }
            catch (const std::length_error&)
            {
                throw_preprocessor_storage_error();
            }
            if (tokens.empty())
            {
                throw_directive_error("invalid empty GROMACS directive");
            }
            if (tokens[0] == "#ifdef" || tokens[0] == "#ifndef")
            {
                if (tokens.size() != 2 || !Gromacs_Is_Macro_Name(tokens[1]))
                {
                    throw_directive_error(
                        "invalid GROMACS conditional directive");
                }
                bool parent_active = is_active();
                bool defined = macros->count(tokens[1]) > 0;
                bool branch_active =
                    (tokens[0] == "#ifdef") ? defined : !defined;
                stack.push_back({parent_active, branch_active, false});
                continue;
            }
            if (tokens[0] == "#else")
            {
                if (tokens.size() != 1 || stack.empty() ||
                    stack.back().else_seen)
                {
                    throw_directive_error("invalid GROMACS #else directive");
                }
                stack.back().branch_active = !stack.back().branch_active;
                stack.back().else_seen = true;
                continue;
            }
            if (tokens[0] == "#endif")
            {
                if (tokens.size() != 1 || stack.empty())
                {
                    throw_directive_error("invalid GROMACS #endif directive");
                }
                stack.pop_back();
                continue;
            }
            // Non-conditional directives in an inactive branch have no
            // effect.  In particular, unsupported directives are only an
            // error when their branch is active.
            if (!is_active())
            {
                continue;
            }
            if (tokens[0] == "#include")
            {
                std::string argument =
                    Gromacs_Trim(line.substr(std::string("#include").size()));
                bool quoted = argument.size() >= 2 && argument.front() == '"' &&
                              argument.back() == '"';
                bool angled = argument.size() >= 2 && argument.front() == '<' &&
                              argument.back() == '>';
                if (!quoted && !angled)
                {
                    throw_directive_error("invalid GROMACS #include directive");
                }
                std::string include_name =
                    argument.substr(1, argument.size() - 2);
                if (include_name.empty())
                {
                    throw_directive_error("invalid GROMACS #include directive");
                }
                fs::path include_path = Gromacs_Resolve_Include(
                    file_path.parent_path(), include_name, include_dirs,
                    controller, error_by, file_path, logical_line_number);
                Gromacs_Preprocess_File(include_path, macros, include_dirs,
                                        lines, include_stack, controller,
                                        error_by);
                continue;
            }
            if (tokens[0] == "#define")
            {
                if (tokens.size() < 2 || !Gromacs_Is_Macro_Name(tokens[1]))
                {
                    throw_directive_error(
                        "unsupported or invalid GROMACS #define directive");
                }
                try
                {
                    (*macros)[tokens[1]] = std::vector<std::string>(
                        tokens.begin() + 2, tokens.end());
                }
                catch (const std::bad_alloc&)
                {
                    throw_preprocessor_storage_error();
                }
                catch (const std::length_error&)
                {
                    throw_preprocessor_storage_error();
                }
                continue;
            }
            if (tokens[0] == "#undef")
            {
                if (tokens.size() != 2 || !Gromacs_Is_Macro_Name(tokens[1]))
                {
                    throw_directive_error("invalid GROMACS #undef directive");
                }
                macros->erase(tokens[1]);
                continue;
            }
            if (tokens[0] == "#error")
            {
                throw_directive_error("active GROMACS #error directive");
            }
            throw_directive_error("unsupported GROMACS directive '" +
                                  tokens[0] + "'");
        }

        if (!is_active())
        {
            continue;
        }
        std::vector<std::string> input_tokens;
        try
        {
            input_tokens = Gromacs_Split(line);
        }
        catch (const std::bad_alloc&)
        {
            throw_preprocessor_storage_error();
        }
        catch (const std::length_error&)
        {
            throw_preprocessor_storage_error();
        }
        std::vector<std::string> expanded_tokens;
        const std::string* cycle_macro = NULL;
        const Gromacs_Macro_Expansion_Status expansion_status =
            Gromacs_Expand_Macro_Tokens(input_tokens, *macros,
                                        &expanded_tokens, &cycle_macro);
        if (expansion_status == Gromacs_Macro_Expansion_Status::kCycle)
        {
            std::string reason =
                "Reason:\n\tGROMACS macro expansion cycle involving '" +
                *cycle_macro + "' at " + file_path.string() + ":" +
                std::to_string(logical_line_number) + "\n";
            controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                           reason.c_str());
        }
        if (expansion_status ==
            Gromacs_Macro_Expansion_Status::kCapacityError)
        {
            throw_preprocessor_storage_error();
        }
        std::string expanded_line;
        if (!Gromacs_Join_Expanded_Tokens(expanded_tokens, &expanded_line))
        {
            throw_preprocessor_storage_error();
        }
        if (!expanded_line.empty())
        {
            if (lines->size() == lines->max_size())
            {
                throw_preprocessor_storage_error();
            }
            try
            {
                lines->push_back({std::move(expanded_line), file_path,
                                  logical_line_number});
            }
            catch (const std::bad_alloc&)
            {
                throw_preprocessor_storage_error();
            }
            catch (const std::length_error&)
            {
                throw_preprocessor_storage_error();
            }
        }
    }

    if (fin.bad())
    {
        std::string reason =
            "Reason:\n\tfailed while reading GROMACS topology file '" +
            file_path.string() + "'\n";
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                       reason.c_str());
    }
    fin.clear();
    fin.close();
    if (fin.fail())
    {
        std::string reason =
            "Reason:\n\tfailed while closing GROMACS topology file '" +
            file_path.string() + "'\n";
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                       reason.c_str());
    }

    if (!stack.empty())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tunterminated GROMACS preprocessor conditional\n");
    }
    include_stack->pop_back();
}

static bool Gromacs_Is_Recognized_Section(const std::string& section)
{
    static const std::set<std::string> recognized_sections = {
        "defaults",       "atomtypes",      "bondtypes",      "angletypes",
        "dihedraltypes",  "pairtypes",      "nonbond_params", "constrainttypes",
        "cmaptypes",      "moleculetype",   "atoms",          "bonds",
        "pairs",          "angles",         "dihedrals",      "settles",
        "constraints",    "exclusions",     "cmap",           "virtual_sites1",
        "virtual_sites2", "virtual_sites3", "virtual_sites4", "virtual_sitesn",
        "system",         "molecules"};
    return recognized_sections.count(section) > 0;
}

static void Gromacs_Throw_Topology_Error(CONTROLLER* controller,
                                         const char* error_by,
                                         const Gromacs_Source_Line& source_line,
                                         const std::string& reason)
{
    std::string message = "Reason:\n\t" + reason + " at " +
                          source_line.file_path.string() + ":" +
                          std::to_string(source_line.line_number) + "\n";
    controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                   message.c_str());
}

static bool Gromacs_Interaction_Parameters_Are_Supported(
    const std::vector<float>& parameters, std::size_t parameter_numbers)
{
    if (parameters.empty() || parameters.size() == parameter_numbers)
    {
        return true;
    }
    if (parameters.size() != 2 * parameter_numbers)
    {
        return false;
    }

    // Some standard force fields spell out identical A/B-state parameters
    // even in non-perturbed topologies.  Loading a genuinely different B
    // state as if it were A would be silent data loss, so only identical
    // duplicates are accepted until GROMACS free-energy states are supported.
    for (std::size_t i = 0; i < parameter_numbers; i++)
    {
        if (parameters[i] != parameters[i + parameter_numbers])
        {
            return false;
        }
    }
    return true;
}

// This target is built with -ffast-math, under which ordinary finite checks can
// be optimized away.  The opaque external helpers inspect IEEE-754 bits without
// exposing the floating-point value to that assumption.
static bool Gromacs_Is_Finite(float value)
{
    return Float_Memory_Is_Finite(&value);
}

static bool Gromacs_Is_Finite(double value)
{
    return Double_Memory_Is_Finite(&value);
}

static bool Gromacs_Is_Decimal_Number_Token(const std::string& token)
{
    if (token.empty())
    {
        return false;
    }
    std::size_t i = 0;
    if (token[i] == '+' || token[i] == '-')
    {
        i++;
    }
    bool mantissa_digit = false;
    while (i < token.size() && token[i] >= '0' && token[i] <= '9')
    {
        mantissa_digit = true;
        i++;
    }
    if (i < token.size() && token[i] == '.')
    {
        i++;
        while (i < token.size() && token[i] >= '0' && token[i] <= '9')
        {
            mantissa_digit = true;
            i++;
        }
    }
    if (!mantissa_digit)
    {
        return false;
    }
    if (i < token.size() && (token[i] == 'e' || token[i] == 'E'))
    {
        i++;
        if (i < token.size() && (token[i] == '+' || token[i] == '-'))
        {
            i++;
        }
        bool exponent_digit = false;
        while (i < token.size() && token[i] >= '0' && token[i] <= '9')
        {
            exponent_digit = true;
            i++;
        }
        if (!exponent_digit)
        {
            return false;
        }
    }
    return i == token.size();
}

static std::pair<double, double> Gromacs_Get_C6_C12_From_Sigma_Epsilon(
    double sigma_nm, double epsilon_kj)
{
    const double sigma = 10.0 * std::fabs(sigma_nm);
    const double epsilon = epsilon_kj / 4.184;
    const double sigma6 = std::pow(sigma, 6.0);
    const double c12 = 4.0 * epsilon * sigma6 * sigma6;
    // GROMACS uses a negative sigma as a sentinel for a purely repulsive
    // interaction: C6 is zero, while C12 is calculated from |sigma|.
    const double c6 = sigma_nm < 0.0f ? 0.0 : 4.0 * epsilon * sigma6;
    return {c6, c12};
}

static std::pair<double, double> Gromacs_Get_C6_C12(
    const Gromacs_Defaults& defaults, const Gromacs_Atom_Type& atom_i,
    const Gromacs_Atom_Type& atom_j)
{
    if (defaults.comb_rule == 1)
    {
        const double c6_product = static_cast<double>(atom_i.v) * atom_j.v;
        const double c12_product = static_cast<double>(atom_i.w) * atom_j.w;
        return {std::sqrt(c6_product) * 1000000.0 / 4.184,
                std::sqrt(c12_product) * 1000000000000.0 / 4.184};
    }

    if (defaults.comb_rule == 2)
    {
        const double epsilon_product = static_cast<double>(atom_i.w) * atom_j.w;
        double sigma_nm = 0.5 * (std::fabs(static_cast<double>(atom_i.v)) +
                                 std::fabs(static_cast<double>(atom_j.v)));
        if (atom_i.v < 0.0f || atom_j.v < 0.0f)
        {
            sigma_nm = -sigma_nm;
        }
        return Gromacs_Get_C6_C12_From_Sigma_Epsilon(
            sigma_nm, std::sqrt(epsilon_product));
    }

    // GROMACS combines the magnitudes geometrically and propagates a negative
    // sigma sentinel from either atom type.  Converting the two self pairs
    // first happens to work for ordinary positive epsilon, but loses the
    // official sign/error semantics for negative parameters.
    const double epsilon_product = static_cast<double>(atom_i.w) * atom_j.w;
    double sigma_nm = std::sqrt(std::fabs(static_cast<double>(atom_i.v) *
                                          static_cast<double>(atom_j.v)));
    if (atom_i.v < 0.0f || atom_j.v < 0.0f)
    {
        sigma_nm = -sigma_nm;
    }
    return Gromacs_Get_C6_C12_From_Sigma_Epsilon(sigma_nm,
                                                 std::sqrt(epsilon_product));
}

static std::pair<double, double> Gromacs_Get_C6_C12_From_Pair_Parameters(
    const Gromacs_Defaults& defaults, const std::vector<float>& parameters)
{
    if (parameters.size() < 2)
    {
        return {0.0f, 0.0f};
    }
    if (defaults.comb_rule == 1)
    {
        return {static_cast<double>(parameters[0]) * 1000000.0 / 4.184,
                static_cast<double>(parameters[1]) * 1000000000000.0 / 4.184};
    }
    return Gromacs_Get_C6_C12_From_Sigma_Epsilon(parameters[0], parameters[1]);
}

static int Gromacs_Count_Wildcards(const std::vector<std::string>& values)
{
    int count = 0;
    for (const std::string& value : values)
    {
        if (value == "X")
        {
            count++;
        }
    }
    return count;
}

static void Gromacs_Build_Type_Indices(Gromacs_Topology* topology)
{
    topology->bond_type_index.clear();
    topology->bond_type_index.reserve(topology->bond_types.size());
    for (std::size_t i = 0; i < topology->bond_types.size(); i++)
    {
        const Gromacs_Bond_Type& type = topology->bond_types[i];
        topology->bond_type_index[Gromacs_Make_Type_Key<2>(
            {type.ai, type.aj}, type.funct, true)] = i;
    }

    topology->constraint_type_index.clear();
    topology->constraint_type_index.reserve(topology->constraint_types.size());
    for (std::size_t i = 0; i < topology->constraint_types.size(); i++)
    {
        const Gromacs_Constraint_Type& type = topology->constraint_types[i];
        topology->constraint_type_index[Gromacs_Make_Type_Key<2>(
            {type.ai, type.aj}, type.funct, true)] = i;
    }

    topology->angle_type_index.clear();
    topology->angle_type_index.reserve(topology->angle_types.size());
    for (std::size_t i = 0; i < topology->angle_types.size(); i++)
    {
        const Gromacs_Angle_Type& type = topology->angle_types[i];
        topology->angle_type_index[Gromacs_Make_Type_Key<3>(
            {type.ai, type.aj, type.ak}, type.funct, true)] = i;
    }

    topology->dihedral_type_index.clear();
    topology->dihedral_type_index.reserve(topology->dihedral_types.size());
    for (std::size_t i = 0; i < topology->dihedral_types.size(); i++)
    {
        const Gromacs_Dihedral_Type& type = topology->dihedral_types[i];
        Gromacs_Quadruple_Type_Key key = Gromacs_Make_Type_Key<4>(
            {type.ai, type.aj, type.ak, type.al}, type.funct, true);
        std::vector<std::size_t>& group = topology->dihedral_type_index[key];
        bool continues_previous_group = false;
        if (type.funct == 9 && i > 0)
        {
            const Gromacs_Dihedral_Type& previous =
                topology->dihedral_types[i - 1];
            continues_previous_group =
                Gromacs_Make_Type_Key<4>(
                    {previous.ai, previous.aj, previous.ak, previous.al},
                    previous.funct, true) == key;
        }
        if (!continues_previous_group)
        {
            group.clear();
        }
        group.push_back(i);
    }

    topology->pair_type_index.clear();
    topology->pair_type_index.reserve(topology->pair_types.size());
    for (std::size_t i = 0; i < topology->pair_types.size(); i++)
    {
        const Gromacs_Pair_Type& type = topology->pair_types[i];
        topology->pair_type_index[Gromacs_Make_Type_Key<2>(
            {type.ai, type.aj}, type.funct, true)] = i;
    }

    topology->nonbond_parameter_index.clear();
    topology->nonbond_parameter_index.reserve(topology->nonbond_params.size());
    for (std::size_t i = 0; i < topology->nonbond_params.size(); i++)
    {
        const Gromacs_Pair_Type& type = topology->nonbond_params[i];
        topology->nonbond_parameter_index[Gromacs_Make_Type_Key<2>(
            {type.ai, type.aj}, type.funct, true)] = i;
    }

    topology->cmap_type_index.clear();
    topology->cmap_type_index.reserve(topology->cmap_types.size());
    for (std::size_t i = 0; i < topology->cmap_types.size(); i++)
    {
        const Gromacs_CMap_Type& type = topology->cmap_types[i];
        // Unlike the other bonded types, reversing a CMAP swaps its two grid
        // axes.  Keep its five-type key directional unless the grid is also
        // transposed.
        topology->cmap_type_index[Gromacs_Make_Type_Key<5>(
            {type.ai, type.aj, type.ak, type.al, type.am}, type.funct, false)] =
            i;
    }
}

static const Gromacs_Bond_Type* Gromacs_Find_Bond_Type(
    const Gromacs_Topology& topology, const std::string& ai,
    const std::string& aj, int funct)
{
    auto iter = topology.bond_type_index.find(
        Gromacs_Make_Type_Key<2>({ai, aj}, funct, true));
    if (iter == topology.bond_type_index.end())
    {
        return NULL;
    }
    return &topology.bond_types[iter->second];
}

static const Gromacs_Constraint_Type* Gromacs_Find_Constraint_Type(
    const Gromacs_Topology& topology, const std::string& ai,
    const std::string& aj, int funct)
{
    auto iter = topology.constraint_type_index.find(
        Gromacs_Make_Type_Key<2>({ai, aj}, funct, true));
    if (iter == topology.constraint_type_index.end())
    {
        return NULL;
    }
    return &topology.constraint_types[iter->second];
}

static const Gromacs_Angle_Type* Gromacs_Find_Angle_Type(
    const Gromacs_Topology& topology, const std::string& ai,
    const std::string& aj, const std::string& ak, int funct)
{
    auto iter = topology.angle_type_index.find(
        Gromacs_Make_Type_Key<3>({ai, aj, ak}, funct, true));
    if (iter == topology.angle_type_index.end())
    {
        return NULL;
    }
    return &topology.angle_types[iter->second];
}

static std::vector<const Gromacs_Dihedral_Type*> Gromacs_Find_Dihedral_Types(
    const Gromacs_Topology& topology, const std::string& ai,
    const std::string& aj, const std::string& ak, const std::string& al,
    int funct)
{
    const std::vector<std::size_t>* selected_group = NULL;
    std::size_t selected = 0;
    int best_wildcards = 1000;
    std::array<std::string, 4> values{ai, aj, ak, al};
    for (unsigned int mask = 0; mask < 16; mask++)
    {
        std::array<std::string, 4> pattern = values;
        for (std::size_t i = 0; i < pattern.size(); i++)
        {
            if ((mask & (1U << i)) != 0)
            {
                pattern[i] = "X";
            }
        }
        auto iter = topology.dihedral_type_index.find(
            Gromacs_Make_Type_Key<4>(std::move(pattern), funct, true));
        if (iter == topology.dihedral_type_index.end() || iter->second.empty())
        {
            continue;
        }
        std::size_t candidate = iter->second.back();
        const Gromacs_Dihedral_Type& type = topology.dihedral_types[candidate];
        int wildcards =
            Gromacs_Count_Wildcards({type.ai, type.aj, type.ak, type.al});
        if (wildcards < best_wildcards ||
            (wildcards == best_wildcards &&
             (selected_group == NULL || candidate > selected)))
        {
            best_wildcards = wildcards;
            selected = candidate;
            selected_group = &iter->second;
        }
    }

    std::vector<const Gromacs_Dihedral_Type*> matches;
    if (selected_group == NULL)
    {
        return matches;
    }
    if (funct != 9)
    {
        matches.push_back(&topology.dihedral_types[selected]);
        return matches;
    }
    matches.reserve(selected_group->size());
    for (std::size_t index : *selected_group)
    {
        matches.push_back(&topology.dihedral_types[index]);
    }
    return matches;
}

static const Gromacs_Pair_Type* Gromacs_Find_Pair_Type(
    const Gromacs_Topology& topology, const std::string& ai,
    const std::string& aj, int funct)
{
    auto iter = topology.pair_type_index.find(
        Gromacs_Make_Type_Key<2>({ai, aj}, funct, true));
    if (iter == topology.pair_type_index.end())
    {
        return NULL;
    }
    return &topology.pair_types[iter->second];
}

static const Gromacs_Pair_Type* Gromacs_Find_Nonbond_Parameter(
    const Gromacs_Topology& topology, const std::string& ai,
    const std::string& aj, int funct)
{
    auto iter = topology.nonbond_parameter_index.find(
        Gromacs_Make_Type_Key<2>({ai, aj}, funct, true));
    if (iter == topology.nonbond_parameter_index.end())
    {
        return NULL;
    }
    return &topology.nonbond_params[iter->second];
}

static int Gromacs_Find_CMap_Type(const Gromacs_Topology& topology,
                                  const std::string& ai, const std::string& aj,
                                  const std::string& ak, const std::string& al,
                                  const std::string& am, int funct)
{
    int match = -1;
    int best_wildcards = 1000;
    std::array<std::string, 5> values{ai, aj, ak, al, am};
    for (unsigned int mask = 0; mask < 32; mask++)
    {
        std::array<std::string, 5> pattern = values;
        for (std::size_t i = 0; i < pattern.size(); i++)
        {
            if ((mask & (1U << i)) != 0)
            {
                pattern[i] = "X";
            }
        }
        auto iter = topology.cmap_type_index.find(
            Gromacs_Make_Type_Key<5>(std::move(pattern), funct, false));
        if (iter == topology.cmap_type_index.end())
        {
            continue;
        }
        std::size_t candidate = iter->second;
        const Gromacs_CMap_Type& type = topology.cmap_types[candidate];
        int wildcards = Gromacs_Count_Wildcards(
            {type.ai, type.aj, type.ak, type.al, type.am});
        if (wildcards < best_wildcards ||
            (wildcards == best_wildcards &&
             (match < 0 || candidate > static_cast<std::size_t>(match))))
        {
            best_wildcards = wildcards;
            match = static_cast<int>(candidate);
        }
    }
    return match;
}

static Gromacs_Topology Gromacs_Parse_Topology(CONTROLLER* controller)
{
    const char* error_by = "Xponge::Load_Gromacs_Inputs";
    if (!controller->Command_Exist("gromacs_top"))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand, error_by,
            "Reason:\n\tgromacs_top is required for GROMACS input\n");
    }
    if (!controller->Command_Exist("gromacs_gro"))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand, error_by,
            "Reason:\n\tgromacs_gro is required for GROMACS input\n");
    }

    fs::path top_path = controller->Original_Command("gromacs_top");
    std::vector<fs::path> include_dirs;
    include_dirs.push_back(top_path.parent_path());
    if (controller->Command_Exist("gromacs_include_dir"))
    {
        for (const std::string& token :
             Gromacs_Split_List(
                 controller->Original_Command("gromacs_include_dir")))
        {
            include_dirs.push_back(token);
        }
    }
    Gromacs_Macros macros;
    if (controller->Command_Exist("gromacs_define"))
    {
        for (const std::string& token :
             Gromacs_Split_List(
                 controller->Original_Command("gromacs_define")))
        {
            std::string definition = token;
            if (definition.rfind("-D", 0) == 0)
            {
                definition.erase(0, 2);
            }
            std::size_t equal = definition.find('=');
            std::string name = definition.substr(0, equal);
            if (!Gromacs_Is_Macro_Name(name))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tinvalid macro name in gromacs_define\n");
            }
            if (equal == std::string::npos)
            {
                // Match the conventional command-line -DNAME expansion.
                macros[name] = {"1"};
            }
            else
            {
                std::string replacement = definition.substr(equal + 1);
                macros[name] = replacement.empty() ? std::vector<std::string>{}
                                                   : Gromacs_Split(replacement);
            }
        }
    }

    std::vector<Gromacs_Source_Line> lines;
    std::vector<fs::path> include_stack;
    Gromacs_Preprocess_File(top_path, &macros, include_dirs, &lines,
                            &include_stack, controller, error_by);

    Gromacs_Topology topology;
    std::unordered_map<std::string, std::uint32_t> source_file_ids;
    std::string current_section;
    Gromacs_Molecule* current_molecule = NULL;
    bool defaults_section_seen = false;
    bool defaults_seen = false;

    for (const Gromacs_Source_Line& source_line : lines)
    {
        std::string source_file_key = source_line.file_path.string();
        auto source_file_iter = source_file_ids.find(source_file_key);
        if (source_file_iter == source_file_ids.end())
        {
            if (topology.source_files.size() >=
                std::numeric_limits<std::uint32_t>::max())
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "too many GROMACS topology source files");
            }
            std::uint32_t file_id =
                static_cast<std::uint32_t>(topology.source_files.size());
            topology.source_files.push_back(source_line.file_path);
            source_file_iter =
                source_file_ids.emplace(std::move(source_file_key), file_id)
                    .first;
        }
        if (source_line.line_number > std::numeric_limits<std::uint32_t>::max())
        {
            Gromacs_Throw_Topology_Error(
                controller, error_by, source_line,
                "GROMACS topology source line number is out of range");
        }
        Gromacs_Source_Reference source_reference{
            source_file_iter->second,
            static_cast<std::uint32_t>(source_line.line_number)};
        const std::string& line = source_line.text;
        if (line.front() == '[' && line.back() == ']')
        {
            current_section = Gromacs_Trim(line.substr(1, line.size() - 2));
            if (!Gromacs_Is_Recognized_Section(current_section))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS topology section [ " +
                        current_section + " ]");
            }
            if (current_section == "defaults")
            {
                if (defaults_section_seen)
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "multiple GROMACS [ defaults ] sections are not "
                        "allowed");
                }
                defaults_section_seen = true;
            }
            continue;
        }

        std::vector<std::string> tokens = Gromacs_Split(line);
        if (tokens.empty())
        {
            continue;
        }
        auto require_molecule_atom =
            [&](int atom_index, const std::string& section)
        {
            if (current_molecule == NULL || atom_index < 1 ||
                static_cast<std::size_t>(atom_index) >
                    current_molecule->atoms.size())
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "atom index out of range in GROMACS [ " + section + " ]");
            }
        };
        auto parse_integer = [&](const std::string& token) -> int
        {
            if (!Gromacs_Is_Integer_Token(token))
            {
                Gromacs_Throw_Topology_Error(controller, error_by, source_line,
                                             "invalid integer field '" + token +
                                                 "' in GROMACS topology");
            }
            try
            {
                return std::stoi(token);
            }
            catch (const std::exception&)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "integer field out of range in GROMACS topology");
            }
            return 0;
        };
        auto parse_real = [&](const std::string& token) -> float
        {
            if (!Gromacs_Is_Decimal_Number_Token(token))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid real field '" + token + "' in GROMACS topology");
            }
            try
            {
                std::size_t consumed = 0;
                const double parsed = std::stod(token, &consumed);
                const double maximum_float =
                    static_cast<double>(std::numeric_limits<float>::max());
                if (consumed != token.size() || !Gromacs_Is_Finite(parsed))
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "non-finite or invalid real field '" + token +
                            "' in GROMACS topology");
                }
                if (std::fabs(parsed) > maximum_float)
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "real field '" + token +
                            "' is outside the finite float range used by "
                            "SPONGE");
                }
                const float stored = static_cast<float>(parsed);
                if (!Gromacs_Is_Finite(stored) ||
                    (parsed != 0.0 && stored == 0.0f))
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "nonzero real field '" + token +
                            "' cannot be represented by the SPONGE force "
                            "kernel");
                }
                if (!Float_Memory_Is_Zero_Or_Normal(&stored))
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "real field '" + token +
                            "' is a subnormal float; SPONGE requires finite "
                            "zero or normal topology values for consistent "
                            "FTZ behavior");
                }
                return stored;
            }
            catch (const std::invalid_argument&)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid real field in GROMACS topology");
            }
            catch (const std::out_of_range&)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "real field out of range in GROMACS topology");
            }
            return 0.0f;
        };
        auto parse_multiplicity = [&](const std::string& token) -> float
        {
            int value = parse_integer(token);
            constexpr int max_exact_float_integer = 1 << 24;
            if (value < -99)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "GROMACS dihedral multiplicity must not be less than "
                    "-99");
            }
            if (value > max_exact_float_integer)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "GROMACS dihedral multiplicity cannot be represented "
                    "exactly");
            }
            return static_cast<float>(value);
        };
        auto append_virtual_site =
            [&](Gromacs_Virtual_Site virtual_site, const std::string& section)
        {
            require_molecule_atom(virtual_site.site, section);
            for (int from : virtual_site.from)
            {
                require_molecule_atom(from, section);
                if (from == virtual_site.site)
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "GROMACS [ " + section +
                            " ] target cannot also be a constructing atom");
                }
            }
            for (const Gromacs_Virtual_Site& previous :
                 current_molecule->virtual_sites)
            {
                if (previous.site == virtual_site.site)
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "GROMACS virtual-site target is defined more than "
                        "once");
                }
            }
            virtual_site.source = source_reference;
            current_molecule->virtual_sites.push_back(std::move(virtual_site));
        };

        if (current_section == "defaults")
        {
            if (tokens.size() < 2 || tokens.size() > 6)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ defaults ] entry in GROMACS topology");
            }
            if (defaults_seen)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "multiple GROMACS [ defaults ] entries are not allowed");
            }
            defaults_seen = true;
            topology.defaults.nbfunc = parse_integer(tokens[0]);
            topology.defaults.comb_rule = parse_integer(tokens[1]);
            if (tokens.size() >= 3 && !Gromacs_Is_True(tokens[2]) &&
                !Gromacs_Is_False(tokens[2]))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid gen-pairs value in GROMACS [ defaults ]");
            }
            if (tokens.size() >= 3)
            {
                topology.defaults.gen_pairs = Gromacs_Is_True(tokens[2]);
            }
            if (tokens.size() >= 4)
            {
                topology.defaults.fudge_lj = parse_real(tokens[3]);
            }
            if (tokens.size() >= 5)
            {
                topology.defaults.fudge_qq = parse_real(tokens[4]);
            }
            if (topology.defaults.nbfunc != 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ defaults ] nonbonded function " +
                        std::to_string(topology.defaults.nbfunc));
            }
            if (topology.defaults.comb_rule < 1 ||
                topology.defaults.comb_rule > 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ defaults ] combination rule " +
                        std::to_string(topology.defaults.comb_rule));
            }
            if (!Gromacs_Is_Finite(topology.defaults.fudge_lj) ||
                !Gromacs_Is_Finite(topology.defaults.fudge_qq) ||
                topology.defaults.fudge_lj < 0.0f ||
                topology.defaults.fudge_qq < 0.0f)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "GROMACS [ defaults ] fudge factors must be finite and "
                    "non-negative");
            }
            if (tokens.size() == 6)
            {
                float repulsion_power = parse_real(tokens[5]);
                if (!Gromacs_Is_Finite(repulsion_power) ||
                    repulsion_power != 12.0f)
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "only repulsion power 12 is supported in GROMACS [ "
                        "defaults ]");
                }
            }
        }
        else if (current_section == "atomtypes")
        {
            if (tokens.size() < 6 || tokens.size() > 8)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ atomtypes ] entry in GROMACS topology");
            }
            Gromacs_Atom_Type atom_type;
            atom_type.source = source_reference;
            atom_type.name = tokens[0];
            atom_type.bond_type = atom_type.name;
            if (tokens.size() == 7 && !Gromacs_Is_Integer_Token(tokens[1]))
            {
                atom_type.bond_type = tokens[1];
            }
            else if (tokens.size() == 8)
            {
                atom_type.bond_type = tokens[1];
                if (!Gromacs_Is_Integer_Token(tokens[2]))
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "invalid atomic number in GROMACS [ atomtypes ]");
                }
            }
            atom_type.mass = parse_real(tokens[tokens.size() - 5]);
            atom_type.charge = parse_real(tokens[tokens.size() - 4]);
            atom_type.ptype = tokens[tokens.size() - 3];
            atom_type.v = parse_real(tokens[tokens.size() - 2]);
            atom_type.w = parse_real(tokens[tokens.size() - 1]);
            if (atom_type.ptype != "A" && atom_type.ptype != "S" &&
                atom_type.ptype != "V" && atom_type.ptype != "D")
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid GROMACS particle type '" + atom_type.ptype + "'");
            }
            topology.atom_types[atom_type.name] = atom_type;
        }
        else if (current_section == "bondtypes")
        {
            if (tokens.size() < 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ bondtypes ] entry in GROMACS topology");
            }
            Gromacs_Bond_Type bond_type;
            bond_type.source = source_reference;
            bond_type.ai = tokens[0];
            bond_type.aj = tokens[1];
            bond_type.funct = parse_integer(tokens[2]);
            std::size_t parameter_numbers = 0;
            if (bond_type.funct == 1 || bond_type.funct == 6)
            {
                parameter_numbers = 2;
            }
            else if (bond_type.funct != 5)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ bondtypes ] function " +
                        std::to_string(bond_type.funct));
            }
            if (tokens.size() != 3 + parameter_numbers)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ bondtypes ] "
                    "function " +
                        std::to_string(bond_type.funct));
            }
            if (parameter_numbers == 2)
            {
                bond_type.b0 = parse_real(tokens[3]);
                bond_type.kb = parse_real(tokens[4]);
            }
            topology.bond_types.push_back(bond_type);
        }
        else if (current_section == "angletypes")
        {
            if (tokens.size() < 4)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ angletypes ] entry in GROMACS topology");
            }
            Gromacs_Angle_Type angle_type;
            angle_type.source = source_reference;
            angle_type.ai = tokens[0];
            angle_type.aj = tokens[1];
            angle_type.ak = tokens[2];
            angle_type.funct = parse_integer(tokens[3]);
            std::size_t expected_tokens = 0;
            if (angle_type.funct == 1)
            {
                expected_tokens = 6;
            }
            else if (angle_type.funct == 5)
            {
                expected_tokens = 8;
            }
            else
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ angletypes ] function " +
                        std::to_string(angle_type.funct));
            }
            if (tokens.size() != expected_tokens)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ angletypes ] "
                    "function " +
                        std::to_string(angle_type.funct));
            }
            angle_type.theta0 = parse_real(tokens[4]);
            angle_type.k = parse_real(tokens[5]);
            if (angle_type.funct == 5)
            {
                angle_type.ub0 = parse_real(tokens[6]);
                angle_type.kub = parse_real(tokens[7]);
            }
            topology.angle_types.push_back(angle_type);
        }
        else if (current_section == "dihedraltypes")
        {
            if (tokens.size() < 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ dihedraltypes ] entry in GROMACS topology");
            }
            Gromacs_Dihedral_Type dihedral_type;
            dihedral_type.source = source_reference;
            std::size_t function_index = 4;
            if (Gromacs_Is_Integer_Token(tokens[2]))
            {
                // The official two-type shorthand specifies only the two
                // central bonded types; the terminal types are wildcards.
                dihedral_type.ai = "X";
                dihedral_type.aj = tokens[0];
                dihedral_type.ak = tokens[1];
                dihedral_type.al = "X";
                function_index = 2;
            }
            else
            {
                if (tokens.size() < 5)
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "invalid [ dihedraltypes ] entry in GROMACS "
                        "topology");
                }
                dihedral_type.ai = tokens[0];
                dihedral_type.aj = tokens[1];
                dihedral_type.ak = tokens[2];
                dihedral_type.al = tokens[3];
            }
            dihedral_type.funct = parse_integer(tokens[function_index]);
            std::size_t parameter_numbers = 0;
            if (dihedral_type.funct == 1 || dihedral_type.funct == 4 ||
                dihedral_type.funct == 9)
            {
                parameter_numbers = 3;
            }
            else if (dihedral_type.funct == 2)
            {
                parameter_numbers = 2;
            }
            else if (dihedral_type.funct == 3)
            {
                parameter_numbers = 6;
            }
            else if (dihedral_type.funct == 5)
            {
                parameter_numbers = 4;
            }
            else
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ dihedraltypes ] function " +
                        std::to_string(dihedral_type.funct));
            }
            if (tokens.size() != function_index + 1 + parameter_numbers)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ dihedraltypes ] "
                    "function " +
                        std::to_string(dihedral_type.funct));
            }
            for (std::size_t i = function_index + 1; i < tokens.size(); i++)
            {
                bool is_multiplicity =
                    parameter_numbers == 3 && i == function_index + 3;
                dihedral_type.parameters.push_back(
                    is_multiplicity ? parse_multiplicity(tokens[i])
                                    : parse_real(tokens[i]));
            }
            topology.dihedral_types.push_back(dihedral_type);
        }
        else if (current_section == "pairtypes")
        {
            if (tokens.size() < 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ pairtypes ] entry in GROMACS topology");
            }
            Gromacs_Pair_Type pair_type;
            pair_type.source = source_reference;
            pair_type.ai = tokens[0];
            pair_type.aj = tokens[1];
            pair_type.funct = parse_integer(tokens[2]);
            if (pair_type.funct != 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ pairtypes ] function " +
                        std::to_string(pair_type.funct));
            }
            if (tokens.size() != 5)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ pairtypes ] "
                    "function 1");
            }
            for (std::size_t i = 3; i < tokens.size(); i++)
            {
                pair_type.parameters.push_back(parse_real(tokens[i]));
            }
            topology.pair_types.push_back(pair_type);
        }
        else if (current_section == "nonbond_params")
        {
            if (tokens.size() != 5)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ nonbond_params ] entry in GROMACS topology");
            }
            Gromacs_Pair_Type nonbond_parameter;
            nonbond_parameter.source = source_reference;
            nonbond_parameter.ai = tokens[0];
            nonbond_parameter.aj = tokens[1];
            nonbond_parameter.funct = parse_integer(tokens[2]);
            if (nonbond_parameter.funct != 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ nonbond_params ] function " +
                        std::to_string(nonbond_parameter.funct));
            }
            for (std::size_t i = 3; i < tokens.size(); i++)
            {
                nonbond_parameter.parameters.push_back(parse_real(tokens[i]));
            }
            topology.nonbond_params.push_back(nonbond_parameter);
        }
        else if (current_section == "constrainttypes")
        {
            if (tokens.size() < 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ constrainttypes ] entry in GROMACS topology");
            }
            Gromacs_Constraint_Type constraint_type;
            constraint_type.source = source_reference;
            constraint_type.ai = tokens[0];
            constraint_type.aj = tokens[1];
            constraint_type.funct = parse_integer(tokens[2]);
            if (constraint_type.funct != 1 && constraint_type.funct != 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ constrainttypes ] function " +
                        std::to_string(constraint_type.funct));
            }
            if (tokens.size() != 4)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ constrainttypes ] "
                    "function " +
                        std::to_string(constraint_type.funct));
            }
            constraint_type.distance = parse_real(tokens[3]);
            topology.constraint_types.push_back(constraint_type);
        }
        else if (current_section == "cmaptypes")
        {
            if (tokens.size() < 8)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ cmaptypes ] entry in GROMACS topology");
            }
            Gromacs_CMap_Type cmap_type;
            cmap_type.source = source_reference;
            cmap_type.ai = tokens[0];
            cmap_type.aj = tokens[1];
            cmap_type.ak = tokens[2];
            cmap_type.al = tokens[3];
            cmap_type.am = tokens[4];
            cmap_type.funct = parse_integer(tokens[5]);
            if (cmap_type.funct != 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ cmaptypes ] function " +
                        std::to_string(cmap_type.funct));
            }
            int resolution_phi = parse_integer(tokens[6]);
            int resolution_psi = parse_integer(tokens[7]);
            if (resolution_phi != resolution_psi)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "non-square GROMACS CMAP grids are not supported yet");
            }
            if (resolution_phi <= 0)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "GROMACS CMAP resolution must be positive");
            }
            cmap_type.resolution = resolution_phi;
            std::size_t resolution =
                static_cast<std::size_t>(cmap_type.resolution);
            if (resolution >
                (std::numeric_limits<std::size_t>::max() - 8) / resolution)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "GROMACS CMAP grid size overflows the host size type");
            }
            std::size_t expected_grid_size = resolution * resolution;
            if (expected_grid_size >
                static_cast<std::size_t>(std::numeric_limits<int>::max()) / 16)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "GROMACS CMAP interpolation table exceeds the supported "
                    "kernel index range");
            }
            if (tokens.size() != 8 + expected_grid_size)
            {
                Gromacs_Throw_Topology_Error(controller, error_by, source_line,
                                             "invalid GROMACS CMAP grid size");
            }
            cmap_type.grid.reserve(expected_grid_size);
            for (std::size_t i = 8; i < tokens.size(); i++)
            {
                cmap_type.grid.push_back(parse_real(tokens[i]));
            }
            topology.cmap_types.push_back(cmap_type);
        }
        else if (current_section == "moleculetype")
        {
            if (tokens.size() != 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ moleculetype ] entry in GROMACS topology");
            }
            Gromacs_Molecule molecule;
            molecule.name = tokens[0];
            molecule.nrexcl = parse_integer(tokens[1]);
            if (molecule.nrexcl < 0)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "negative nrexcl in GROMACS [ moleculetype ]");
            }
            if (topology.molecules.count(molecule.name) > 0)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "duplicate GROMACS molecule type '" + molecule.name + "'");
            }
            topology.molecules[molecule.name] = molecule;
            current_molecule = &topology.molecules[molecule.name];
        }
        else if (current_section == "atoms")
        {
            if (current_molecule == NULL || tokens.size() < 6 ||
                tokens.size() > 8)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid or unsupported [ atoms ] entry in GROMACS "
                    "topology");
            }
            Gromacs_Molecule_Atom atom;
            atom.source = source_reference;
            atom.nr = parse_integer(tokens[0]);
            atom.type = tokens[1];
            atom.resnr = parse_integer(tokens[2]);
            atom.residue = tokens[3];
            atom.atom = tokens[4];
            atom.cgnr = parse_integer(tokens[5]);
            if (tokens.size() >= 7)
            {
                atom.charge = parse_real(tokens[6]);
                atom.has_charge = true;
            }
            if (atom.nr != static_cast<int>(current_molecule->atoms.size()) + 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "GROMACS [ atoms ] indices must be contiguous and start "
                    "at 1");
            }
            if (tokens.size() >= 8)
            {
                atom.mass = parse_real(tokens[7]);
                atom.has_mass = true;
            }
            current_molecule->atoms.push_back(atom);
        }
        else if (current_section == "virtual_sites1")
        {
            if (current_molecule == NULL || tokens.size() != 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ virtual_sites1 ] entry in GROMACS topology");
            }
            Gromacs_Virtual_Site virtual_site;
            virtual_site.site = parse_integer(tokens[0]);
            virtual_site.from = {parse_integer(tokens[1])};
            append_virtual_site(std::move(virtual_site), "virtual_sites1");
        }
        else if (current_section == "virtual_sites2")
        {
            if (current_molecule == NULL || tokens.size() < 4)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ virtual_sites2 ] entry in GROMACS topology");
            }
            Gromacs_Virtual_Site virtual_site;
            virtual_site.site = parse_integer(tokens[0]);
            virtual_site.from = {parse_integer(tokens[1]),
                                 parse_integer(tokens[2])};
            virtual_site.funct = parse_integer(tokens[3]);
            if (virtual_site.funct != 1 && virtual_site.funct != 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ virtual_sites2 ] function " +
                        std::to_string(virtual_site.funct));
            }
            if (tokens.size() != 5)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ virtual_sites2 ] "
                    "function " +
                        std::to_string(virtual_site.funct));
            }
            virtual_site.parameters = {parse_real(tokens[4])};
            append_virtual_site(std::move(virtual_site), "virtual_sites2");
        }
        else if (current_section == "virtual_sites3")
        {
            if (current_molecule == NULL || tokens.size() < 5)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ virtual_sites3 ] entry in GROMACS topology");
            }
            Gromacs_Virtual_Site virtual_site;
            virtual_site.site = parse_integer(tokens[0]);
            virtual_site.from = {parse_integer(tokens[1]),
                                 parse_integer(tokens[2]),
                                 parse_integer(tokens[3])};
            virtual_site.funct = parse_integer(tokens[4]);
            if (virtual_site.funct != 1 && virtual_site.funct != 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ virtual_sites3 ] function " +
                        std::to_string(virtual_site.funct));
            }
            if (tokens.size() != 7)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ virtual_sites3 ] "
                    "function " +
                        std::to_string(virtual_site.funct));
            }
            virtual_site.parameters = {parse_real(tokens[5]),
                                       parse_real(tokens[6])};
            append_virtual_site(std::move(virtual_site), "virtual_sites3");
        }
        else if (current_section == "virtual_sites4" ||
                 current_section == "virtual_sitesn")
        {
            Gromacs_Throw_Topology_Error(
                controller, error_by, source_line,
                "unsupported GROMACS [ " + current_section + " ] entry");
        }
        else if (current_section == "bonds")
        {
            if (current_molecule == NULL || tokens.size() < 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ bonds ] entry in GROMACS topology");
            }
            Gromacs_Bond bond;
            bond.source = source_reference;
            bond.ai = parse_integer(tokens[0]);
            bond.aj = parse_integer(tokens[1]);
            bond.funct = parse_integer(tokens[2]);
            require_molecule_atom(bond.ai, "bonds");
            require_molecule_atom(bond.aj, "bonds");
            std::size_t parameter_numbers = 0;
            if (bond.funct == 1 || bond.funct == 6)
            {
                parameter_numbers = 2;
            }
            else if (bond.funct != 5)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ bonds ] function " +
                        std::to_string(bond.funct));
            }
            for (std::size_t i = 3; i < tokens.size(); i++)
            {
                bond.parameters.push_back(parse_real(tokens[i]));
            }
            if (!Gromacs_Interaction_Parameters_Are_Supported(
                    bond.parameters, parameter_numbers))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid or perturbed parameter list for GROMACS [ bonds "
                    "] function " +
                        std::to_string(bond.funct));
            }
            current_molecule->bonds.push_back(bond);
        }
        else if (current_section == "pairs")
        {
            if (current_molecule == NULL || tokens.size() < 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ pairs ] entry in GROMACS topology");
            }
            Gromacs_Pair pair;
            pair.source = source_reference;
            pair.ai = parse_integer(tokens[0]);
            pair.aj = parse_integer(tokens[1]);
            pair.funct = parse_integer(tokens[2]);
            require_molecule_atom(pair.ai, "pairs");
            require_molecule_atom(pair.aj, "pairs");
            if (pair.funct != 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ pairs ] function " +
                        std::to_string(pair.funct));
            }
            for (std::size_t i = 3; i < tokens.size(); i++)
            {
                pair.parameters.push_back(parse_real(tokens[i]));
            }
            if (!Gromacs_Interaction_Parameters_Are_Supported(pair.parameters,
                                                              2))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid or perturbed parameter list for GROMACS [ pairs "
                    "] function 1");
            }
            current_molecule->pairs.push_back(pair);
        }
        else if (current_section == "angles")
        {
            if (current_molecule == NULL || tokens.size() < 4)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ angles ] entry in GROMACS topology");
            }
            Gromacs_Angle angle;
            angle.source = source_reference;
            angle.ai = parse_integer(tokens[0]);
            angle.aj = parse_integer(tokens[1]);
            angle.ak = parse_integer(tokens[2]);
            angle.funct = parse_integer(tokens[3]);
            require_molecule_atom(angle.ai, "angles");
            require_molecule_atom(angle.aj, "angles");
            require_molecule_atom(angle.ak, "angles");
            std::size_t parameter_numbers = 0;
            if (angle.funct == 1)
            {
                parameter_numbers = 2;
            }
            else if (angle.funct == 5)
            {
                parameter_numbers = 4;
            }
            else
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ angles ] function " +
                        std::to_string(angle.funct));
            }
            for (std::size_t i = 4; i < tokens.size(); i++)
            {
                angle.parameters.push_back(parse_real(tokens[i]));
            }
            if (!Gromacs_Interaction_Parameters_Are_Supported(
                    angle.parameters, parameter_numbers))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid or perturbed parameter list for GROMACS [ angles "
                    "] function " +
                        std::to_string(angle.funct));
            }
            current_molecule->angles.push_back(angle);
        }
        else if (current_section == "dihedrals")
        {
            if (current_molecule == NULL || tokens.size() < 5)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ dihedrals ] entry in GROMACS topology");
            }
            Gromacs_Dihedral dihedral;
            dihedral.source = source_reference;
            dihedral.ai = parse_integer(tokens[0]);
            dihedral.aj = parse_integer(tokens[1]);
            dihedral.ak = parse_integer(tokens[2]);
            dihedral.al = parse_integer(tokens[3]);
            dihedral.funct = parse_integer(tokens[4]);
            require_molecule_atom(dihedral.ai, "dihedrals");
            require_molecule_atom(dihedral.aj, "dihedrals");
            require_molecule_atom(dihedral.ak, "dihedrals");
            require_molecule_atom(dihedral.al, "dihedrals");
            std::size_t parameter_numbers = 0;
            if (dihedral.funct == 1 || dihedral.funct == 4 ||
                dihedral.funct == 9)
            {
                parameter_numbers = 3;
            }
            else if (dihedral.funct == 2)
            {
                parameter_numbers = 2;
            }
            else if (dihedral.funct == 3)
            {
                parameter_numbers = 6;
            }
            else if (dihedral.funct == 5)
            {
                parameter_numbers = 4;
            }
            else
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ dihedrals ] function " +
                        std::to_string(dihedral.funct));
            }
            for (std::size_t i = 5; i < tokens.size(); i++)
            {
                bool is_multiplicity =
                    parameter_numbers == 3 && (i - 5) % 3 == 2;
                dihedral.parameters.push_back(
                    is_multiplicity ? parse_multiplicity(tokens[i])
                                    : parse_real(tokens[i]));
            }
            if (!Gromacs_Interaction_Parameters_Are_Supported(
                    dihedral.parameters, parameter_numbers))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid or perturbed parameter list for GROMACS [ "
                    "dihedrals ] function " +
                        std::to_string(dihedral.funct));
            }
            current_molecule->dihedrals.push_back(dihedral);
        }
        else if (current_section == "settles")
        {
            if (current_molecule == NULL || tokens.size() != 4)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ settles ] entry");
            }
            Gromacs_Settle settle;
            settle.source = source_reference;
            settle.ow = parse_integer(tokens[0]);
            settle.funct = parse_integer(tokens[1]);
            if (settle.ow > std::numeric_limits<int>::max() - 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "atom index out of range in GROMACS [ settles ]");
            }
            require_molecule_atom(settle.ow, "settles");
            require_molecule_atom(settle.ow + 1, "settles");
            require_molecule_atom(settle.ow + 2, "settles");
            if (settle.funct != 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ settles ] function " +
                        std::to_string(settle.funct));
            }
            settle.doh = parse_real(tokens[2]);
            settle.dhh = parse_real(tokens[3]);
            current_molecule->settles.push_back(settle);
        }
        else if (current_section == "constraints")
        {
            if (current_molecule == NULL || tokens.size() < 3)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ constraints ] entry in GROMACS topology");
            }
            Gromacs_Constraint constraint;
            constraint.source = source_reference;
            constraint.ai = parse_integer(tokens[0]);
            constraint.aj = parse_integer(tokens[1]);
            constraint.funct = parse_integer(tokens[2]);
            require_molecule_atom(constraint.ai, "constraints");
            require_molecule_atom(constraint.aj, "constraints");
            if (constraint.funct != 1 && constraint.funct != 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ constraints ] function " +
                        std::to_string(constraint.funct));
            }
            for (std::size_t i = 3; i < tokens.size(); i++)
            {
                constraint.parameters.push_back(parse_real(tokens[i]));
            }
            if (!Gromacs_Interaction_Parameters_Are_Supported(
                    constraint.parameters, 1))
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid or perturbed parameter list for GROMACS [ "
                    "constraints ] function " +
                        std::to_string(constraint.funct));
            }
            current_molecule->constraints.push_back(constraint);
        }
        else if (current_section == "exclusions")
        {
            if (current_molecule == NULL || tokens.size() < 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ exclusions ] entry in GROMACS topology");
            }
            int atom_i = parse_integer(tokens[0]);
            require_molecule_atom(atom_i, "exclusions");
            for (std::size_t i = 1; i < tokens.size(); i++)
            {
                int atom_j = parse_integer(tokens[i]);
                if (atom_j < 1 ||
                    static_cast<std::size_t>(atom_j) >
                        current_molecule->atoms.size() ||
                    atom_j == atom_i)
                {
                    Gromacs_Throw_Topology_Error(
                        controller, error_by, source_line,
                        "invalid atom pair in GROMACS [ exclusions ]");
                }
                current_molecule->exclusions.push_back({atom_i, atom_j});
            }
        }
        else if (current_section == "cmap")
        {
            if (current_molecule == NULL || tokens.size() != 6)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid parameter count for GROMACS [ cmap ] entry");
            }
            Gromacs_CMap cmap;
            cmap.source = source_reference;
            cmap.ai = parse_integer(tokens[0]);
            cmap.aj = parse_integer(tokens[1]);
            cmap.ak = parse_integer(tokens[2]);
            cmap.al = parse_integer(tokens[3]);
            cmap.am = parse_integer(tokens[4]);
            cmap.funct = parse_integer(tokens[5]);
            require_molecule_atom(cmap.ai, "cmap");
            require_molecule_atom(cmap.aj, "cmap");
            require_molecule_atom(cmap.ak, "cmap");
            require_molecule_atom(cmap.al, "cmap");
            require_molecule_atom(cmap.am, "cmap");
            if (cmap.funct != 1)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "unsupported GROMACS [ cmap ] function " +
                        std::to_string(cmap.funct));
            }
            current_molecule->cmaps.push_back(cmap);
        }
        else if (current_section == "molecules")
        {
            if (tokens.size() != 2)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "invalid [ molecules ] entry in GROMACS topology");
            }
            int molecule_count = parse_integer(tokens[1]);
            if (molecule_count < 0)
            {
                Gromacs_Throw_Topology_Error(
                    controller, error_by, source_line,
                    "negative molecule count in GROMACS [ molecules ]");
            }
            topology.system_molecules.push_back(
                {tokens[0], molecule_count, source_reference});
        }
        else if (current_section != "system")
        {
            Gromacs_Throw_Topology_Error(
                controller, error_by, source_line,
                "data outside a supported GROMACS topology section");
        }
    }

    if (!defaults_seen)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tGROMACS topology is missing [ defaults ]\n");
    }
    for (auto& molecule_item : topology.molecules)
    {
        std::vector<std::pair<int, int>>& exclusions =
            molecule_item.second.exclusions;
        for (std::pair<int, int>& exclusion : exclusions)
        {
            if (exclusion.first > exclusion.second)
            {
                std::swap(exclusion.first, exclusion.second);
            }
        }
        std::sort(exclusions.begin(), exclusions.end());
        exclusions.erase(std::unique(exclusions.begin(), exclusions.end()),
                         exclusions.end());
    }
    Gromacs_Build_Type_Indices(&topology);
    return topology;
}

static void Gromacs_Throw_Reference_Error(
    CONTROLLER* controller, const char* error_by,
    const Gromacs_Topology& topology,
    const Gromacs_Source_Reference& source_reference, const std::string& reason)
{
    std::string message = "Reason:\n\t" + reason;
    if (source_reference.file_id < topology.source_files.size())
    {
        message += " at " +
                   topology.source_files[source_reference.file_id].string() +
                   ":" + std::to_string(source_reference.line_number);
    }
    message += "\n";
    controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                   message.c_str());
}

static std::vector<int> Gromacs_Validate_Virtual_Sites(
    const Gromacs_Topology& topology, const Gromacs_Molecule& molecule,
    CONTROLLER* controller, const char* error_by)
{
    std::vector<int> target_record(molecule.atoms.size(), -1);
    for (std::size_t record_index = 0;
         record_index < molecule.virtual_sites.size(); record_index++)
    {
        const Gromacs_Virtual_Site& virtual_site =
            molecule.virtual_sites[record_index];
        int target = virtual_site.site - 1;
        if (target < 0 || target >= static_cast<int>(molecule.atoms.size()))
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, virtual_site.source,
                "GROMACS virtual-site target atom index is out of range");
        }
        if (target_record[target] >= 0)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, virtual_site.source,
                "GROMACS virtual-site target is defined more than once");
        }
        target_record[target] = static_cast<int>(record_index);
        std::set<int> unique_sources;
        for (int source_atom : virtual_site.from)
        {
            int source = source_atom - 1;
            if (source < 0 || source >= static_cast<int>(molecule.atoms.size()))
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, virtual_site.source,
                    "GROMACS virtual-site constructing atom index is out of "
                    "range");
            }
            if (source == target)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, virtual_site.source,
                    "GROMACS virtual-site target cannot also be a "
                    "constructing atom");
            }
            if (!unique_sources.insert(source_atom).second)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, virtual_site.source,
                    "GROMACS virtual-site constructing atom is listed more "
                    "than once");
            }
        }
    }

    std::vector<int> indegree(molecule.virtual_sites.size(), 0);
    std::vector<std::vector<int>> consumers(molecule.virtual_sites.size());
    for (std::size_t record_index = 0;
         record_index < molecule.virtual_sites.size(); record_index++)
    {
        std::vector<int> dependencies;
        for (int source_atom : molecule.virtual_sites[record_index].from)
        {
            int dependency = target_record[source_atom - 1];
            if (dependency < 0 ||
                std::find(dependencies.begin(), dependencies.end(),
                          dependency) != dependencies.end())
            {
                continue;
            }
            dependencies.push_back(dependency);
            indegree[record_index]++;
            consumers[dependency].push_back(static_cast<int>(record_index));
        }
    }
    std::vector<int> ready;
    for (std::size_t record_index = 0; record_index < indegree.size();
         record_index++)
    {
        if (indegree[record_index] == 0)
        {
            ready.push_back(static_cast<int>(record_index));
        }
    }
    std::size_t processed = 0;
    for (std::size_t head = 0; head < ready.size(); head++)
    {
        int record_index = ready[head];
        processed++;
        for (int consumer : consumers[record_index])
        {
            indegree[consumer]--;
            if (indegree[consumer] == 0)
            {
                ready.push_back(consumer);
            }
        }
    }
    if (processed != molecule.virtual_sites.size())
    {
        for (std::size_t record_index = 0; record_index < indegree.size();
             record_index++)
        {
            if (indegree[record_index] > 0)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology,
                    molecule.virtual_sites[record_index].source,
                    "GROMACS virtual-site source graph contains a dependency "
                    "cycle");
            }
        }
    }
    return target_record;
}

static void Gromacs_Clean_Virtual_Site_Bondeds(
    const Gromacs_Topology& topology, Gromacs_Molecule* molecule,
    const std::vector<int>& target_record, CONTROLLER* controller,
    const char* error_by)
{
    if (molecule->virtual_sites.empty())
    {
        return;
    }

    auto virtual_site_for_atom =
        [&](int one_based_atom) -> const Gromacs_Virtual_Site*
    {
        if (one_based_atom < 1 ||
            one_based_atom > static_cast<int>(target_record.size()))
        {
            return NULL;
        }
        int record_index = target_record[one_based_atom - 1];
        return record_index < 0 ? NULL : &molecule->virtual_sites[record_index];
    };
    auto same_sources = [](const Gromacs_Virtual_Site& first,
                           const Gromacs_Virtual_Site& second)
    {
        if (first.from.size() != second.from.size())
        {
            return false;
        }
        for (int atom : second.from)
        {
            if (std::find(first.from.begin(), first.from.end(), atom) ==
                first.from.end())
            {
                return false;
            }
        }
        return true;
    };
    auto source_position =
        [](const Gromacs_Virtual_Site& virtual_site, int atom)
    {
        auto iter =
            std::find(virtual_site.from.begin(), virtual_site.from.end(), atom);
        return iter == virtual_site.from.end()
                   ? -1
                   : static_cast<int>(iter - virtual_site.from.begin());
    };

    // GROMACS decides whether an angle is constant from the original
    // constraint graph.  SETTLE is represented there by its first listed
    // atom pair; ordinary constraints contribute their explicit pair.
    std::vector<std::vector<int>> constraint_neighbors(molecule->atoms.size());
    auto add_constraint_edge = [&](int first, int second)
    {
        constraint_neighbors[first - 1].push_back(second);
        constraint_neighbors[second - 1].push_back(first);
    };
    for (const Gromacs_Constraint& constraint : molecule->constraints)
    {
        add_constraint_edge(constraint.ai, constraint.aj);
    }
    for (const Gromacs_Settle& settle : molecule->settles)
    {
        add_constraint_edge(settle.ow, settle.ow + 1);
    }
    auto sources_are_cyclically_constrained =
        [&](const Gromacs_Virtual_Site& virtual_site)
    {
        for (std::size_t i = 0; i < virtual_site.from.size(); i++)
        {
            int first = virtual_site.from[i];
            int second = virtual_site.from[(i + 1) % virtual_site.from.size()];
            const std::vector<int>& neighbors = constraint_neighbors[first - 1];
            if (std::find(neighbors.begin(), neighbors.end(), second) ==
                neighbors.end())
            {
                return false;
            }
        }
        return true;
    };
    auto is_fixed_distance_construction =
        [](const Gromacs_Virtual_Site& virtual_site)
    {
        // This matches clean_vsite_bonds in GROMACS: the immediate fixed-
        // distance removal applies to VSITE3FD (and the unsupported 3FAD/4FD
        // families), but not to VSITE2FD.
        return virtual_site.from.size() == 3 && virtual_site.funct == 2;
    };
    auto two_atom_term_is_constant = [&](int first_atom, int second_atom)
    {
        const int atoms[2] = {first_atom, second_atom};
        const Gromacs_Virtual_Site* first_virtual_site = NULL;
        int virtual_site_count = 0;
        bool all_fixed_distance = true;
        for (int i = 0; i < 2; i++)
        {
            const Gromacs_Virtual_Site* virtual_site =
                virtual_site_for_atom(atoms[i]);
            if (virtual_site == NULL)
            {
                continue;
            }
            virtual_site_count++;
            const bool fixed_distance =
                is_fixed_distance_construction(*virtual_site);
            all_fixed_distance = all_fixed_distance && fixed_distance;
            const int other_atom = atoms[1 - i];
            if (fixed_distance && virtual_site_for_atom(other_atom) == NULL &&
                other_atom == virtual_site->from.front())
            {
                return true;
            }
            if (first_virtual_site == NULL)
            {
                first_virtual_site = virtual_site;
            }
            else if (!same_sources(*first_virtual_site, *virtual_site))
            {
                return false;
            }
        }
        if (virtual_site_count == 0)
        {
            return false;
        }

        bool nonvirtual_atoms_are_first_two_sources = true;
        for (int atom : atoms)
        {
            if (virtual_site_for_atom(atom) != NULL)
            {
                continue;
            }
            int position = source_position(*first_virtual_site, atom);
            if (position < 0)
            {
                return false;
            }
            nonvirtual_atoms_are_first_two_sources =
                nonvirtual_atoms_are_first_two_sources && position < 2;
        }
        if (all_fixed_distance && nonvirtual_atoms_are_first_two_sources)
        {
            return true;
        }
        return sources_are_cyclically_constrained(*first_virtual_site);
    };
    auto angle_is_constant = [&](const Gromacs_Angle& angle)
    {
        const int atoms[3] = {angle.ai, angle.aj, angle.ak};
        const Gromacs_Virtual_Site* first_virtual_site = NULL;
        int virtual_site_count = 0;
        bool all_three_atom_fixed_angle_distance = true;
        for (int atom : atoms)
        {
            const Gromacs_Virtual_Site* virtual_site =
                virtual_site_for_atom(atom);
            if (virtual_site == NULL)
            {
                continue;
            }
            virtual_site_count++;
            all_three_atom_fixed_angle_distance =
                all_three_atom_fixed_angle_distance &&
                virtual_site->from.size() == 3 && virtual_site->funct == 3;
            if (first_virtual_site == NULL)
            {
                first_virtual_site = virtual_site;
            }
            else if (!same_sources(*first_virtual_site, *virtual_site))
            {
                return false;
            }
        }
        if (virtual_site_count == 0 || first_virtual_site->from.size() > 3)
        {
            return false;
        }

        bool nonvirtual_atoms_are_first_two_sources = true;
        for (int atom : atoms)
        {
            if (virtual_site_for_atom(atom) != NULL)
            {
                continue;
            }
            int position = source_position(*first_virtual_site, atom);
            if (position < 0)
            {
                return false;
            }
            nonvirtual_atoms_are_first_two_sources =
                nonvirtual_atoms_are_first_two_sources && position < 2;
        }
        if (all_three_atom_fixed_angle_distance &&
            nonvirtual_atoms_are_first_two_sources)
        {
            return true;
        }
        return sources_are_cyclically_constrained(*first_virtual_site);
    };
    auto dihedral_is_constant = [&](const Gromacs_Dihedral& dihedral)
    {
        const int atoms[4] = {dihedral.ai, dihedral.aj, dihedral.ak,
                              dihedral.al};
        const Gromacs_Virtual_Site* first_virtual_site = NULL;
        int virtual_site_count = 0;
        for (int atom : atoms)
        {
            const Gromacs_Virtual_Site* virtual_site =
                virtual_site_for_atom(atom);
            if (virtual_site == NULL)
            {
                continue;
            }
            virtual_site_count++;
            if (first_virtual_site == NULL)
            {
                first_virtual_site = virtual_site;
            }
            else if (!same_sources(*first_virtual_site, *virtual_site))
            {
                return false;
            }
        }
        if (virtual_site_count == 0)
        {
            return false;
        }
        for (int atom : atoms)
        {
            if (virtual_site_for_atom(atom) == NULL &&
                source_position(*first_virtual_site, atom) < 0)
            {
                return false;
            }
        }
        return true;
    };

    std::vector<Gromacs_Bond> cleaned_bonds;
    cleaned_bonds.reserve(molecule->bonds.size());
    for (Gromacs_Bond bond : molecule->bonds)
    {
        if (bond.funct == 1 && two_atom_term_is_constant(bond.ai, bond.aj))
        {
            // Chemical connectivity must survive cleanup because GROMACS
            // generates nrexcl exclusions before removing constant forces.
            bond.funct = 5;
            bond.parameters.clear();
        }
        cleaned_bonds.push_back(std::move(bond));
    }
    molecule->bonds = std::move(cleaned_bonds);

    std::vector<Gromacs_Angle> cleaned_angles;
    cleaned_angles.reserve(molecule->angles.size());
    for (const Gromacs_Angle& angle : molecule->angles)
    {
        if ((angle.funct == 1 || angle.funct == 5) && angle_is_constant(angle))
        {
            continue;
        }
        cleaned_angles.push_back(angle);
    }
    molecule->angles = std::move(cleaned_angles);

    std::vector<Gromacs_Dihedral> cleaned_dihedrals;
    cleaned_dihedrals.reserve(molecule->dihedrals.size());
    for (const Gromacs_Dihedral& dihedral : molecule->dihedrals)
    {
        // Input functions 1 and 9 both become F_PDIHS in GROMACS; function 2
        // becomes F_IDIHS.  Other supported dihedral families are retained.
        const bool cleanup_function =
            dihedral.funct == 1 || dihedral.funct == 2 || dihedral.funct == 9;
        if (cleanup_function && dihedral_is_constant(dihedral))
        {
            continue;
        }
        cleaned_dihedrals.push_back(dihedral);
    }
    molecule->dihedrals = std::move(cleaned_dihedrals);

    std::vector<Gromacs_Constraint> cleaned_constraints;
    cleaned_constraints.reserve(molecule->constraints.size());
    for (const Gromacs_Constraint& constraint : molecule->constraints)
    {
        if (two_atom_term_is_constant(constraint.ai, constraint.aj))
        {
            if (constraint.funct == 1)
            {
                Gromacs_Bond connection;
                connection.ai = constraint.ai;
                connection.aj = constraint.aj;
                connection.funct = 5;
                connection.source = constraint.source;
                molecule->bonds.push_back(std::move(connection));
            }
            continue;
        }
        if (virtual_site_for_atom(constraint.ai) != NULL ||
            virtual_site_for_atom(constraint.aj) != NULL)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, constraint.source,
                "GROMACS constraint involving a virtual-site target remains "
                "after constant virtual-site bonded cleanup");
        }
        cleaned_constraints.push_back(constraint);
    }
    molecule->constraints = std::move(cleaned_constraints);

    for (const Gromacs_Settle& settle : molecule->settles)
    {
        if (virtual_site_for_atom(settle.ow) != NULL ||
            virtual_site_for_atom(settle.ow + 1) != NULL ||
            virtual_site_for_atom(settle.ow + 2) != NULL)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, settle.source,
                "GROMACS SETTLE constraint involving a virtual-site target "
                "cannot be retained after virtual-site bonded cleanup");
        }
    }
}

static void Gromacs_Throw_Gro_Error(CONTROLLER* controller,
                                    const char* error_by,
                                    const fs::path& gro_path,
                                    std::size_t line_number,
                                    const std::string& reason)
{
    std::string message = "Reason:\n\t" + reason + " at " + gro_path.string();
    if (line_number > 0)
    {
        message += ":" + std::to_string(line_number);
    }
    message += "\n";
    controller->Throw_SPONGE_Error(spongeErrorBadFileFormat, error_by,
                                   message.c_str());
}

struct Gromacs_Expanded_System_Size
{
    std::size_t atom_count = 0;
    std::size_t molecule_copy_count = 0;
    std::size_t residue_count = 0;
    std::size_t virtual_site_count = 0;
    std::size_t bond_count = 0;
    std::size_t constraint_count = 0;
    std::size_t exclusion_count = 0;
    std::size_t cmap_count = 0;
    std::size_t angle_count = 0;
    std::size_t proper_dihedral_count = 0;
    std::size_t improper_dihedral_count = 0;
    std::size_t pair_count = 0;
    std::size_t lj_type_count = 0;
    std::size_t lj_pair_count = 0;
    std::size_t cmap_type_count = 0;
    std::size_t cmap_gridpoint_count = 0;
    std::unordered_map<std::string, int> atom_type_lj_id;
    std::vector<std::string> ordered_lj_types;
    std::unordered_map<const Gromacs_Molecule*,
                       std::vector<std::pair<int, int>>>
        molecule_exclusion_pairs;
};

static std::vector<std::pair<int, int>> Gromacs_Build_Molecule_Exclusions(
    const Gromacs_Topology& topology, const Gromacs_Molecule& molecule,
    std::size_t maximum_pair_count,
    const Gromacs_Source_Reference& system_source, CONTROLLER* controller,
    const char* error_by)
{
    const std::size_t atom_count = molecule.atoms.size();
    std::vector<std::vector<int>> adjacency(atom_count);
    auto require_atom = [&](int one_based_atom,
                            const Gromacs_Source_Reference& source,
                            const std::string& section)
    {
        if (one_based_atom < 1 ||
            static_cast<std::size_t>(one_based_atom) > atom_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                "GROMACS [ " + section + " ] atom index is out of range");
        }
    };
    auto add_edge = [&](int first, int second)
    {
        adjacency[first - 1].push_back(second - 1);
        adjacency[second - 1].push_back(first - 1);
    };
    for (const Gromacs_Bond& bond : molecule.bonds)
    {
        require_atom(bond.ai, bond.source, "bonds");
        require_atom(bond.aj, bond.source, "bonds");
        if (bond.funct == 1 || bond.funct == 5)
        {
            add_edge(bond.ai, bond.aj);
        }
    }
    for (const Gromacs_Constraint& constraint : molecule.constraints)
    {
        require_atom(constraint.ai, constraint.source, "constraints");
        require_atom(constraint.aj, constraint.source, "constraints");
        if (constraint.funct == 1)
        {
            add_edge(constraint.ai, constraint.aj);
        }
    }
    for (std::vector<int>& neighbors : adjacency)
    {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()),
                        neighbors.end());
    }

    std::vector<std::pair<int, int>> exclusion_pairs;
    if (molecule.nrexcl > 0)
    {
        std::vector<unsigned int> visit_generation(atom_count, 0);
        std::vector<int> depth(atom_count, 0);
        std::vector<int> queue;
        unsigned int generation = 0;
        for (std::size_t atom = 0; atom < atom_count; atom++)
        {
            generation++;
            if (generation == 0)
            {
                std::fill(visit_generation.begin(), visit_generation.end(), 0);
                generation = 1;
            }
            queue.clear();
            queue.push_back(static_cast<int>(atom));
            visit_generation[atom] = generation;
            depth[atom] = 0;
            for (std::size_t head = 0; head < queue.size(); head++)
            {
                const int current = queue[head];
                if (depth[current] >= molecule.nrexcl)
                {
                    continue;
                }
                for (int next : adjacency[current])
                {
                    if (visit_generation[next] == generation)
                    {
                        continue;
                    }
                    visit_generation[next] = generation;
                    depth[next] = depth[current] + 1;
                    queue.push_back(next);
                    if (static_cast<std::size_t>(next) > atom)
                    {
                        if (exclusion_pairs.size() == maximum_pair_count)
                        {
                            Gromacs_Throw_Reference_Error(
                                controller, error_by, topology, system_source,
                                "expanded GROMACS exclusion count exceeds the "
                                "supported kernel int range");
                        }
                        exclusion_pairs.emplace_back(
                            static_cast<int>(atom), next);
                    }
                }
            }
        }
    }
    std::sort(exclusion_pairs.begin(), exclusion_pairs.end());
    exclusion_pairs.erase(
        std::unique(exclusion_pairs.begin(), exclusion_pairs.end()),
        exclusion_pairs.end());
    if (molecule.exclusions.empty())
    {
        return exclusion_pairs;
    }

    // Explicit exclusions were normalized and deduplicated after parsing.
    // Merge them with the generated list instead of appending all input terms
    // and deduplicating afterwards: the latter can exceed container capacity
    // even when the final unique exclusion set is representable.
    std::vector<std::pair<int, int>> merged_pairs;
    std::size_t generated_index = 0;
    std::size_t explicit_index = 0;
    auto append_pair = [&](const std::pair<int, int>& pair)
    {
        if (!merged_pairs.empty() && merged_pairs.back() == pair)
        {
            return;
        }
        if (merged_pairs.size() == maximum_pair_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, system_source,
                "expanded GROMACS exclusion count exceeds the supported "
                "kernel int range");
        }
        merged_pairs.push_back(pair);
    };
    while (generated_index < exclusion_pairs.size() ||
           explicit_index < molecule.exclusions.size())
    {
        std::pair<int, int> explicit_pair;
        if (explicit_index < molecule.exclusions.size())
        {
            const std::pair<int, int>& exclusion =
                molecule.exclusions[explicit_index];
            require_atom(exclusion.first, Gromacs_Source_Reference{},
                         "exclusions");
            require_atom(exclusion.second, Gromacs_Source_Reference{},
                         "exclusions");
            if (exclusion.first == exclusion.second)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tGROMACS exclusion contains a self pair\n");
            }
            explicit_pair = {exclusion.first - 1, exclusion.second - 1};
        }

        if (explicit_index == molecule.exclusions.size() ||
            (generated_index < exclusion_pairs.size() &&
             exclusion_pairs[generated_index] < explicit_pair))
        {
            append_pair(exclusion_pairs[generated_index++]);
        }
        else if (generated_index == exclusion_pairs.size() ||
                 explicit_pair < exclusion_pairs[generated_index])
        {
            append_pair(explicit_pair);
            explicit_index++;
        }
        else
        {
            append_pair(explicit_pair);
            explicit_index++;
            generated_index++;
        }
    }
    return merged_pairs;
}

// Validate every [ molecules ] reference and compute the complete expanded
// size before any count-driven copy/allocation begins.  Comparing that checked
// sum with the already-open GRO stream prevents a hostile count from growing
// the System merely to discover a small coordinate header afterwards.
static Gromacs_Expanded_System_Size Gromacs_Check_Expanded_System_Size(
    Gromacs_Topology& topology, int gro_atom_count,
    const fs::path& gro_path, CONTROLLER* controller, const char* error_by)
{
    for (const Gromacs_System_Molecule& item : topology.system_molecules)
    {
        if (topology.molecules.find(item.name) == topology.molecules.end())
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "molecule '" + item.name +
                    "' referenced in [ molecules ] is not defined");
        }
    }

    Gromacs_Expanded_System_Size expanded;
    const std::size_t maximum_size =
        std::numeric_limits<std::size_t>::max();
    for (const Gromacs_System_Molecule& item : topology.system_molecules)
    {
        const Gromacs_Molecule& molecule = topology.molecules.at(item.name);
        if (item.count < 0)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "negative molecule count in GROMACS [ molecules ]");
        }
        if (item.count > 0 && molecule.atoms.empty())
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "molecule '" + item.name +
                    "' has no atoms and cannot be instantiated");
        }

        const std::size_t copy_count =
            static_cast<std::size_t>(item.count);
        if (copy_count > maximum_size - expanded.molecule_copy_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "expanded GROMACS molecule copy count overflows the host "
                "size type");
        }
        expanded.molecule_copy_count += copy_count;

        if (!molecule.atoms.empty() &&
            copy_count > maximum_size / molecule.atoms.size())
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "expanded GROMACS atom count overflows the host size type");
        }
        const std::size_t added_atoms = copy_count * molecule.atoms.size();
        if (added_atoms > maximum_size - expanded.atom_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "expanded GROMACS atom count overflows the host size type");
        }
        expanded.atom_count += added_atoms;
    }

    if (expanded.atom_count != static_cast<std::size_t>(gro_atom_count))
    {
        Gromacs_Throw_Gro_Error(
            controller, error_by, gro_path, 2,
            "GROMACS gro atom count does not match the checked expanded "
            "topology atom count (gro " + std::to_string(gro_atom_count) +
                ", topology " + std::to_string(expanded.atom_count) + ")");
    }
    const std::vector<float> component_storage;
    const std::vector<std::vector<int>> nested_index_storage;
    const std::size_t maximum_component_atom_count = std::min(
        {static_cast<std::size_t>(std::numeric_limits<int>::max() / 3),
         maximum_size / 3, component_storage.max_size() / 3,
         nested_index_storage.max_size()});
    if (expanded.atom_count > maximum_component_atom_count)
    {
        Gromacs_Throw_Gro_Error(
            controller, error_by, gro_path, 2,
            "checked expanded GROMACS atom count " +
                std::to_string(expanded.atom_count) +
                " cannot safely represent all 3 * atom count coordinate "
                "values and kernel indices");
    }
    try
    {
    // Virtual-site cleanup changes the exact bonded counts.  Perform it once
    // on the molecule templates before deriving any downstream expansion
    // sizes; instantiation revalidates the resulting virtual-site map.
    for (auto& molecule_item : topology.molecules)
    {
        const std::vector<int> target_record = Gromacs_Validate_Virtual_Sites(
            topology, molecule_item.second, controller, error_by);
        Gromacs_Clean_Virtual_Site_Bondeds(
            topology, &molecule_item.second, target_record, controller,
            error_by);
    }

    struct Molecule_Counts
    {
        std::size_t residues = 0;
        std::size_t virtual_sites = 0;
        std::size_t bonds = 0;
        std::size_t constraints = 0;
        std::size_t exclusions = 0;
        std::size_t cmaps = 0;
        std::size_t angles = 0;
        std::size_t proper_dihedrals = 0;
        std::size_t improper_dihedrals = 0;
        std::size_t pairs = 0;
    };
    std::unordered_map<const Gromacs_Molecule*, Molecule_Counts> count_cache;
    std::unordered_map<const Gromacs_Molecule*, std::size_t>
        molecule_copy_counts;
    for (const Gromacs_System_Molecule& item : topology.system_molecules)
    {
        if (item.count == 0)
        {
            continue;
        }
        const Gromacs_Molecule* molecule = &topology.molecules.at(item.name);
        std::size_t& copy_count = molecule_copy_counts[molecule];
        const std::size_t addition = static_cast<std::size_t>(item.count);
        if (addition > maximum_size - copy_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "expanded GROMACS molecule copy count overflows the host "
                "size type");
        }
        copy_count += addition;
    }
    const std::size_t maximum_kernel_count =
        static_cast<std::size_t>(std::numeric_limits<int>::max());
    const std::vector<int> integer_storage;
    const std::vector<float> real_storage;
    const std::vector<VirtualAtomRecord> virtual_site_storage;
    const std::vector<std::pair<int, int>> exclusion_pair_storage;
    const std::size_t maximum_force_count =
        std::min({maximum_kernel_count, integer_storage.max_size(),
                  real_storage.max_size()});
    const std::size_t maximum_residue_count =
        std::min(maximum_kernel_count, integer_storage.max_size());
    const std::size_t maximum_virtual_site_count =
        std::min(maximum_kernel_count, virtual_site_storage.max_size());
    const std::size_t maximum_exclusion_count =
        std::min(maximum_kernel_count, integer_storage.max_size());
    const std::size_t maximum_lj_pair_count =
        std::min(maximum_kernel_count, real_storage.max_size());
    auto add_local_count = [&](std::size_t addition, std::size_t* count,
                               const Gromacs_Source_Reference& source,
                               const std::string& description)
    {
        if (addition > maximum_size - *count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                "GROMACS molecule " + description +
                    " count overflows the host size type");
        }
        *count += addition;
    };
    auto add_expanded_count = [&](std::size_t molecule_count,
                                  std::size_t copy_count,
                                  std::size_t* expanded_count,
                                  const Gromacs_System_Molecule& item,
                                  const std::string& description,
                                  std::size_t maximum_count)
    {
        if (molecule_count != 0 && copy_count > maximum_size / molecule_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "expanded GROMACS " + description +
                    " count overflows the host size type");
        }
        const std::size_t addition = molecule_count * copy_count;
        if (addition > maximum_size - *expanded_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "expanded GROMACS " + description +
                    " count overflows the host size type");
        }
        if (*expanded_count > maximum_count ||
            addition > maximum_count - *expanded_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "expanded GROMACS " + description +
                    " count exceeds the supported kernel or host container "
                    "range");
        }
        *expanded_count += addition;
    };

    for (const Gromacs_System_Molecule& item : topology.system_molecules)
    {
        if (item.count == 0)
        {
            continue;
        }
        Gromacs_Molecule& molecule = topology.molecules.at(item.name);
        auto cached = count_cache.find(&molecule);
        if (cached == count_cache.end())
        {
            Molecule_Counts counts;
            if (!molecule.atoms.empty())
            {
                counts.residues = 1;
                for (std::size_t atom = 1; atom < molecule.atoms.size(); atom++)
                {
                    if (molecule.atoms[atom].resnr !=
                            molecule.atoms[atom - 1].resnr ||
                        molecule.atoms[atom].residue !=
                            molecule.atoms[atom - 1].residue)
                    {
                        counts.residues++;
                    }
                }
            }
            counts.virtual_sites = molecule.virtual_sites.size();
            auto require_interaction_atom =
                [&](int one_based_atom,
                    const Gromacs_Source_Reference& source,
                    const std::string& section)
                    -> const Gromacs_Molecule_Atom&
            {
                if (one_based_atom < 1 ||
                    static_cast<std::size_t>(one_based_atom) >
                        molecule.atoms.size())
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, source,
                        "GROMACS [ " + section +
                            " ] atom index is out of range");
                }
                const Gromacs_Molecule_Atom& atom =
                    molecule.atoms[one_based_atom - 1];
                if (topology.atom_types.find(atom.type) ==
                    topology.atom_types.end())
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, atom.source,
                        "undefined GROMACS atom type '" + atom.type + "'");
                }
                return atom;
            };
            for (const Gromacs_Bond& bond : molecule.bonds)
            {
                const Gromacs_Molecule_Atom& atom_i =
                    require_interaction_atom(bond.ai, bond.source, "bonds");
                const Gromacs_Molecule_Atom& atom_j =
                    require_interaction_atom(bond.aj, bond.source, "bonds");
                if (bond.funct != 5 && bond.parameters.size() < 2 &&
                    Gromacs_Find_Bond_Type(
                        topology,
                        topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        bond.funct) == NULL)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, bond.source,
                        "failed to find GROMACS bond type");
                }
                if (bond.funct != 5)
                {
                    add_local_count(1, &counts.bonds, bond.source, "bond");
                }
            }
            for (const Gromacs_Constraint& constraint : molecule.constraints)
            {
                const Gromacs_Molecule_Atom& atom_i = require_interaction_atom(
                    constraint.ai, constraint.source, "constraints");
                const Gromacs_Molecule_Atom& atom_j = require_interaction_atom(
                    constraint.aj, constraint.source, "constraints");
                if (constraint.parameters.empty() &&
                    Gromacs_Find_Constraint_Type(
                        topology,
                        topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        constraint.funct) == NULL)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, constraint.source,
                        "failed to resolve GROMACS constraint distance");
                }
            }
            add_local_count(molecule.constraints.size(), &counts.bonds,
                            item.source, "bond");
            add_local_count(molecule.constraints.size(), &counts.constraints,
                            item.source, "constraint");
            if (molecule.settles.size() > maximum_size / 3)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, item.source,
                    "GROMACS molecule SETTLE expansion count overflows the "
                    "host size type");
            }
            const std::size_t settle_terms = 3 * molecule.settles.size();
            for (const Gromacs_Settle& settle : molecule.settles)
            {
                require_interaction_atom(settle.ow, settle.source, "settles");
                require_interaction_atom(settle.ow + 1, settle.source,
                                         "settles");
                require_interaction_atom(settle.ow + 2, settle.source,
                                         "settles");
            }
            add_local_count(settle_terms, &counts.bonds, item.source, "bond");
            add_local_count(settle_terms, &counts.constraints, item.source,
                            "constraint");
            const std::size_t aggregate_copy_count =
                molecule_copy_counts.at(&molecule);
            const std::size_t exclusion_limit = std::min(
                (maximum_exclusion_count - expanded.exclusion_count) /
                    aggregate_copy_count,
                exclusion_pair_storage.max_size());
            std::vector<std::pair<int, int>> exclusion_pairs =
                Gromacs_Build_Molecule_Exclusions(
                    topology, molecule, exclusion_limit, item.source,
                    controller, error_by);
            counts.exclusions = exclusion_pairs.size();
            expanded.molecule_exclusion_pairs.emplace(
                &molecule, std::move(exclusion_pairs));
            add_expanded_count(counts.exclusions, aggregate_copy_count,
                               &expanded.exclusion_count, item, "exclusion",
                               maximum_exclusion_count);
            for (const Gromacs_CMap& cmap : molecule.cmaps)
            {
                const Gromacs_Molecule_Atom& atom_i = require_interaction_atom(
                    cmap.ai, cmap.source, "cmap");
                const Gromacs_Molecule_Atom& atom_j = require_interaction_atom(
                    cmap.aj, cmap.source, "cmap");
                const Gromacs_Molecule_Atom& atom_k = require_interaction_atom(
                    cmap.ak, cmap.source, "cmap");
                const Gromacs_Molecule_Atom& atom_l = require_interaction_atom(
                    cmap.al, cmap.source, "cmap");
                const Gromacs_Molecule_Atom& atom_m = require_interaction_atom(
                    cmap.am, cmap.source, "cmap");
                if (Gromacs_Find_CMap_Type(
                        topology,
                        topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        topology.atom_types.at(atom_k.type).bond_type,
                        topology.atom_types.at(atom_l.type).bond_type,
                        topology.atom_types.at(atom_m.type).bond_type,
                        cmap.funct) < 0)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, cmap.source,
                        "failed to find GROMACS CMAP type");
                }
            }
            counts.cmaps = molecule.cmaps.size();
            for (const Gromacs_Angle& angle : molecule.angles)
            {
                const Gromacs_Molecule_Atom& atom_i = require_interaction_atom(
                    angle.ai, angle.source, "angles");
                const Gromacs_Molecule_Atom& atom_j = require_interaction_atom(
                    angle.aj, angle.source, "angles");
                const Gromacs_Molecule_Atom& atom_k = require_interaction_atom(
                    angle.ak, angle.source, "angles");
                if (angle.parameters.size() < 2 &&
                    Gromacs_Find_Angle_Type(
                        topology,
                        topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        topology.atom_types.at(atom_k.type).bond_type,
                        angle.funct) == NULL)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, angle.source,
                        "failed to find GROMACS angle type");
                }
            }
            counts.angles = molecule.angles.size();
            for (const Gromacs_Pair& pair : molecule.pairs)
            {
                const Gromacs_Molecule_Atom& atom_i = require_interaction_atom(
                    pair.ai, pair.source, "pairs");
                const Gromacs_Molecule_Atom& atom_j = require_interaction_atom(
                    pair.aj, pair.source, "pairs");
                if (pair.parameters.size() < 2 &&
                    Gromacs_Find_Pair_Type(topology, atom_i.type, atom_j.type,
                                           pair.funct) == NULL &&
                    !topology.defaults.gen_pairs)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, pair.source,
                        "failed to resolve GROMACS pair interaction");
                }
            }
            counts.pairs = molecule.pairs.size();
            auto add_dihedral_terms =
                [&](int funct, std::size_t parameter_set_count,
                    const Gromacs_Source_Reference& source)
            {
                const std::size_t terms_per_set =
                    funct == 3 ? 5 : (funct == 5 ? 4 : 1);
                if (parameter_set_count > maximum_size / terms_per_set)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, source,
                        "GROMACS molecule dihedral expansion count overflows "
                        "the host size type");
                }
                const std::size_t terms =
                    parameter_set_count * terms_per_set;
                if (funct == 2)
                {
                    add_local_count(terms, &counts.improper_dihedrals, source,
                                    "improper dihedral");
                }
                else
                {
                    add_local_count(terms, &counts.proper_dihedrals, source,
                                    "proper dihedral");
                }
            };
            for (const Gromacs_Dihedral& dihedral : molecule.dihedrals)
            {
                const Gromacs_Molecule_Atom& atom_i =
                    require_interaction_atom(dihedral.ai, dihedral.source,
                                             "dihedrals");
                const Gromacs_Molecule_Atom& atom_j =
                    require_interaction_atom(dihedral.aj, dihedral.source,
                                             "dihedrals");
                const Gromacs_Molecule_Atom& atom_k =
                    require_interaction_atom(dihedral.ak, dihedral.source,
                                             "dihedrals");
                const Gromacs_Molecule_Atom& atom_l =
                    require_interaction_atom(dihedral.al, dihedral.source,
                                             "dihedrals");
                if (!dihedral.parameters.empty())
                {
                    add_dihedral_terms(dihedral.funct, 1, dihedral.source);
                    continue;
                }
                const std::vector<const Gromacs_Dihedral_Type*> types =
                    Gromacs_Find_Dihedral_Types(
                        topology,
                        topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        topology.atom_types.at(atom_k.type).bond_type,
                        topology.atom_types.at(atom_l.type).bond_type,
                        dihedral.funct);
                if (types.empty())
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, dihedral.source,
                        "failed to find GROMACS dihedral type");
                }
                add_dihedral_terms(dihedral.funct, types.size(),
                                    dihedral.source);
            }
            cached = count_cache.emplace(&molecule, counts).first;
        }

        const std::size_t copy_count = static_cast<std::size_t>(item.count);
        const Molecule_Counts& counts = cached->second;
        add_expanded_count(counts.residues, copy_count,
                           &expanded.residue_count, item, "residue",
                           maximum_residue_count);
        add_expanded_count(counts.virtual_sites, copy_count,
                           &expanded.virtual_site_count, item, "virtual-site",
                           maximum_virtual_site_count);
        add_expanded_count(counts.bonds, copy_count, &expanded.bond_count, item,
                           "bond", maximum_force_count);
        add_expanded_count(counts.constraints, copy_count,
                           &expanded.constraint_count, item, "constraint",
                           maximum_force_count);
        add_expanded_count(counts.cmaps, copy_count, &expanded.cmap_count, item,
                           "CMAP interaction", maximum_force_count);
        add_expanded_count(counts.angles, copy_count, &expanded.angle_count,
                           item, "angle", maximum_force_count);
        add_expanded_count(counts.proper_dihedrals, copy_count,
                           &expanded.proper_dihedral_count, item,
                           "proper dihedral", maximum_force_count);
        add_expanded_count(counts.improper_dihedrals, copy_count,
                           &expanded.improper_dihedral_count, item,
                           "improper dihedral", maximum_force_count);
        add_expanded_count(counts.pairs, copy_count, &expanded.pair_count, item,
                           "pair interaction", maximum_force_count);
    }

    std::unordered_map<std::string, int> lj_key_ids;
    auto triangular_fits = [&](std::size_t type_count)
    {
        std::size_t first = type_count;
        std::size_t second = type_count + 1;
        if (first % 2 == 0)
        {
            first /= 2;
        }
        else
        {
            second /= 2;
        }
        return second == 0 || first <= maximum_lj_pair_count / second;
    };
    std::size_t maximum_lj_type_count =
        std::min({maximum_kernel_count, integer_storage.max_size(),
                  lj_key_ids.max_size(),
                  expanded.ordered_lj_types.max_size()});
    std::size_t lower_type_count = 0;
    std::size_t upper_type_count = maximum_lj_type_count;
    while (lower_type_count < upper_type_count)
    {
        const std::size_t middle =
            lower_type_count + (upper_type_count - lower_type_count + 1) / 2;
        if (triangular_fits(middle))
        {
            lower_type_count = middle;
        }
        else
        {
            upper_type_count = middle - 1;
        }
    }
    maximum_lj_type_count = lower_type_count;
    auto add_lj_key = [&](const std::string& key,
                          const std::string& representative_type)
    {
        const auto existing = lj_key_ids.find(key);
        if (existing != lj_key_ids.end())
        {
            return existing->second;
        }
        if (lj_key_ids.size() == maximum_lj_type_count)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tGROMACS Lennard-Jones type count would exceed "
                "the supported pair-table or host container range\n");
        }
        const int id = static_cast<int>(lj_key_ids.size());
        lj_key_ids.emplace(key, id);
        expanded.ordered_lj_types.push_back(representative_type);
        return id;
    };
    // A nonbond override is keyed by the atom-type names, so those names must
    // remain distinct.  Without overrides, types can share a kernel ID only
    // when their original mixing-rule inputs are identical.  Comparing the
    // rounded self-pair A/B values is insufficient: it would, for example,
    // merge positive and negative epsilon values even though their cross pair
    // is not the same interaction.
    const bool preserve_named_atom_types = !topology.nonbond_params.empty();
    for (const Gromacs_System_Molecule& item : topology.system_molecules)
    {
        if (item.count == 0)
        {
            continue;
        }
        const Gromacs_Molecule& molecule = topology.molecules.at(item.name);
        for (const Gromacs_Molecule_Atom& atom : molecule.atoms)
        {
            const auto atom_type_iter = topology.atom_types.find(atom.type);
            if (atom_type_iter == topology.atom_types.end())
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, atom.source,
                    "undefined GROMACS atom type '" + atom.type + "'");
            }
            const Gromacs_Atom_Type& atom_type = atom_type_iter->second;
            const std::pair<double, double> self_c6_c12 =
                Gromacs_Get_C6_C12(topology.defaults, atom_type, atom_type);
            const double pair_a = 12.0 * self_c6_c12.second;
            const double pair_b = 6.0 * self_c6_c12.first;
            const double maximum_float =
                static_cast<double>(std::numeric_limits<float>::max());
            if (!Gromacs_Is_Finite(pair_a) || !Gromacs_Is_Finite(pair_b) ||
                std::fabs(pair_a) > maximum_float ||
                std::fabs(pair_b) > maximum_float)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, atom_type.source,
                    "GROMACS atom type '" + atom_type.name +
                        "' Lennard-Jones parameters cannot be represented by "
                        "SPONGE");
            }
            const float stored_a = static_cast<float>(pair_a);
            const float stored_b = static_cast<float>(pair_b);
            if (!Gromacs_Is_Finite(stored_a) ||
                !Gromacs_Is_Finite(stored_b) ||
                (pair_a != 0.0 && stored_a == 0.0f) ||
                (pair_b != 0.0 && stored_b == 0.0f) ||
                !Float_Memory_Is_Zero_Or_Normal(&stored_a) ||
                !Float_Memory_Is_Zero_Or_Normal(&stored_b))
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, atom_type.source,
                    "GROMACS atom type '" + atom_type.name +
                        "' Lennard-Jones parameters cannot be represented by "
                        "SPONGE as finite zero or normal floats");
            }
            int lj_type_id = 0;
            if (preserve_named_atom_types)
            {
                lj_type_id = add_lj_key(atom.type, atom.type);
            }
            else
            {
                char parameter_key[2 * sizeof(float)];
                std::memcpy(parameter_key, &atom_type.v, sizeof(float));
                std::memcpy(parameter_key + sizeof(float), &atom_type.w,
                            sizeof(float));
                lj_type_id =
                    add_lj_key(std::string(parameter_key, sizeof(parameter_key)),
                               atom.type);
            }
            const auto mapped =
                expanded.atom_type_lj_id.emplace(atom.type, lj_type_id);
            if (!mapped.second && mapped.first->second != lj_type_id)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorBadFileFormat, error_by,
                    "Reason:\n\tinternal GROMACS atom-type/Lennard-Jones ID "
                    "mapping mismatch\n");
            }
        }
    }
    expanded.lj_type_count = lj_key_ids.size();
    if (expanded.lj_type_count == maximum_size ||
        (expanded.lj_type_count > 0 &&
         expanded.lj_type_count >
             maximum_size / (expanded.lj_type_count + 1)))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tGROMACS Lennard-Jones pair table size overflows the "
            "host size type\n");
    }
    expanded.lj_pair_count =
        expanded.lj_type_count * (expanded.lj_type_count + 1) / 2;
    if (expanded.lj_type_count > maximum_kernel_count ||
        expanded.lj_type_count > integer_storage.max_size() ||
        expanded.lj_pair_count > maximum_lj_pair_count)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tGROMACS Lennard-Jones type or pair-table count "
            "exceeds the supported kernel or host container range\n");
    }

    expanded.cmap_type_count = topology.cmap_types.size();
    if (expanded.cmap_type_count > maximum_kernel_count ||
        expanded.cmap_type_count > integer_storage.max_size())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tGROMACS CMAP type count exceeds the supported kernel "
            "or host container range\n");
    }
    const std::size_t maximum_gridpoint_count =
        std::min(maximum_kernel_count / 16, real_storage.max_size() / 16);
    for (const Gromacs_CMap_Type& cmap_type : topology.cmap_types)
    {
        if (cmap_type.resolution <= 0)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, cmap_type.source,
                "GROMACS CMAP resolution must be positive");
        }
        const std::size_t resolution =
            static_cast<std::size_t>(cmap_type.resolution);
        if (resolution > maximum_size / resolution)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, cmap_type.source,
                "GROMACS CMAP resolution overflows the grid size");
        }
        const std::size_t type_gridpoint_count = resolution * resolution;
        if (cmap_type.grid.size() != type_gridpoint_count ||
            type_gridpoint_count > maximum_gridpoint_count ||
            expanded.cmap_gridpoint_count >
                maximum_gridpoint_count - type_gridpoint_count)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, cmap_type.source,
                "GROMACS CMAP interpolation table exceeds the supported "
                "kernel or host container range");
        }
        expanded.cmap_gridpoint_count += type_gridpoint_count;
    }
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed to allocate storage while preparing and "
            "counting the expanded GROMACS topology\n");
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe prepared GROMACS topology exceeds host "
            "container capacity\n");
    }
    return expanded;
}

static void Gromacs_Instantiate_System(
    Gromacs_Topology& topology, System* system, CONTROLLER* controller,
    const Gromacs_Expanded_System_Size& expanded_size);

static void Gromacs_Refresh_Initial_Virtual_Sites(
    const Gromacs_Topology& topology, System* system, const LTMatrix3& cell,
    const LTMatrix3& rcell, const fs::path& gro_path, CONTROLLER* controller,
    const char* error_by)
{
    const std::vector<VirtualAtomRecord>& records =
        system->virtual_atoms.records;
    if (records.empty())
    {
        return;
    }
    std::vector<Gromacs_Source_Reference> record_sources;
    record_sources.reserve(records.size());
    for (const Gromacs_System_Molecule& system_molecule :
         topology.system_molecules)
    {
        const Gromacs_Molecule& molecule =
            topology.molecules.at(system_molecule.name);
        for (int copy = 0; copy < system_molecule.count; copy++)
        {
            for (const Gromacs_Virtual_Site& virtual_site :
                 molecule.virtual_sites)
            {
                record_sources.push_back(virtual_site.source);
            }
        }
    }
    if (record_sources.size() != records.size())
    {
        Gromacs_Throw_Gro_Error(
            controller, error_by, gro_path, 0,
            "internal GROMACS virtual-site/source mapping mismatch");
    }

    VirtualAtomLayout layout;
    std::string validation_error;
    int atom_numbers = static_cast<int>(system->atoms.coordinate.size() / 3);
    if (!Validate_And_Build_Virtual_Atom_Layout(records, atom_numbers, &layout,
                                                &validation_error))
    {
        Gromacs_Throw_Gro_Error(
            controller, error_by, gro_path, 0,
            "invalid expanded GROMACS virtual-site graph: " + validation_error);
    }
    auto coordinate = [&](int atom)
    {
        return VECTOR{system->atoms.coordinate[3 * atom],
                      system->atoms.coordinate[3 * atom + 1],
                      system->atoms.coordinate[3 * atom + 2]};
    };
    auto finite_vector = [](const VECTOR& value)
    {
        return Gromacs_Is_Finite(value.x) && Gromacs_Is_Finite(value.y) &&
               Gromacs_Is_Finite(value.z) &&
               Float_Memory_Is_Zero_Or_Normal(&value.x) &&
               Float_Memory_Is_Zero_Or_Normal(&value.y) &&
               Float_Memory_Is_Zero_Or_Normal(&value.z);
    };
    auto store_coordinate =
        [&](int atom, const VECTOR& value, std::size_t record_index)
    {
        if (!Gromacs_Is_Finite(static_cast<double>(value.x)) ||
            !Gromacs_Is_Finite(static_cast<double>(value.y)) ||
            !Gromacs_Is_Finite(static_cast<double>(value.z)) ||
            !Float_Memory_Is_Zero_Or_Normal(&value.x) ||
            !Float_Memory_Is_Zero_Or_Normal(&value.y) ||
            !Float_Memory_Is_Zero_Or_Normal(&value.z))
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, record_sources[record_index],
                "GROMACS virtual-site coordinate is not a finite zero "
                "or normal float representable by SPONGE for initial "
                "GRO coordinates from '" +
                    gro_path.string() + "'");
        }
        system->atoms.coordinate[3 * atom] = value.x;
        system->atoms.coordinate[3 * atom + 1] = value.y;
        system->atoms.coordinate[3 * atom + 2] = value.z;
    };

    for (int level = 1; level <= layout.max_level; level++)
    {
        for (std::size_t record_index = 0; record_index < records.size();
             record_index++)
        {
            const VirtualAtomRecord& record = records[record_index];
            if (layout.atom_levels[record.virtual_atom] != level)
            {
                continue;
            }
            const VECTOR r1 = coordinate(record.from[0]);
            VECTOR virtual_coordinate = r1;
            if (record.type == 1)
            {
                VECTOR r21 = Get_Periodic_Displacement(
                    coordinate(record.from[1]), r1, cell, rcell);
                virtual_coordinate = r1 + record.parameter[0] * r21;
            }
            else if (record.type == 2)
            {
                VECTOR r21 = Get_Periodic_Displacement(
                    coordinate(record.from[1]), r1, cell, rcell);
                VECTOR r31 = Get_Periodic_Displacement(
                    coordinate(record.from[2]), r1, cell, rcell);
                virtual_coordinate =
                    r1 + record.parameter[0] * r21 + record.parameter[1] * r31;
            }
            else if (record.type == 3 && record.parameter[0] != 0.0f)
            {
                VECTOR r2 = coordinate(record.from[1]);
                VECTOR r3 = coordinate(record.from[2]);
                VECTOR r21 = Get_Periodic_Displacement(r2, r1, cell, rcell);
                VECTOR r32 = Get_Periodic_Displacement(r3, r2, cell, rcell);
                VECTOR direction = r21 + record.parameter[1] * r32;
                float direction_squared = direction * direction;
                if (!finite_vector(r1) || !finite_vector(r2) ||
                    !finite_vector(r3) || !finite_vector(r21) ||
                    !finite_vector(r32) || !finite_vector(direction) ||
                    !Gromacs_Is_Finite(direction_squared) ||
                    !Float_Memory_Is_Zero_Or_Normal(&direction_squared))
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology,
                        record_sources[record_index],
                        "GROMACS distance virtual site has non-finite "
                        "construction geometry in initial GRO coordinates "
                        "from '" +
                            gro_path.string() + "'");
                }
                if (direction_squared == 0.0f)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology,
                        record_sources[record_index],
                        "GROMACS distance virtual site has zero construction "
                        "direction in initial GRO coordinates from '" +
                            gro_path.string() + "'");
                }
                const float inverse_direction =
                    1.0f / std::sqrt(direction_squared);
                const VECTOR unit_direction = inverse_direction * direction;
                const VECTOR displacement =
                    record.parameter[0] * unit_direction;
                if (!Gromacs_Is_Finite(inverse_direction) ||
                    !Float_Memory_Is_Zero_Or_Normal(&inverse_direction) ||
                    !finite_vector(unit_direction) ||
                    !finite_vector(displacement))
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology,
                        record_sources[record_index],
                        "GROMACS distance virtual-site construction cannot "
                        "be represented for initial GRO coordinates from '" +
                            gro_path.string() + "'");
                }
                virtual_coordinate = r1 + displacement;
            }
            else if (record.type != 3)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology,
                    record_sources[record_index],
                    "unsupported expanded GROMACS virtual-site record type");
            }
            store_coordinate(record.virtual_atom, virtual_coordinate,
                             record_index);
        }
    }
}

static void Gromacs_Load_Gro(Gromacs_Topology& topology, System* system,
                             CONTROLLER* controller)
{
    const char* error_by = "Xponge::Load_Gromacs_Inputs";
    fs::path gro_path = controller->Original_Command("gromacs_gro");
    std::ifstream fin(gro_path);
    if (!fin.is_open())
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, 0,
                                "failed to open GROMACS gro file");
    }

    std::string line;
    if (!std::getline(fin, line))
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, 1,
                                "missing GROMACS gro title line");
    }
    auto is_horizontal_space = [](char ch) { return ch == ' ' || ch == '\t'; };
    double gro_start_time = 0.0;
    bool title_time_found = false;
    for (std::size_t i = 0; i < line.size(); i++)
    {
        if (line[i] != 't' || (i > 0 && !is_horizontal_space(line[i - 1])))
        {
            continue;
        }
        std::size_t equal = i + 1;
        while (equal < line.size() && is_horizontal_space(line[equal]))
        {
            equal++;
        }
        if (equal >= line.size() || line[equal] != '=')
        {
            continue;
        }
        std::size_t begin = equal + 1;
        while (begin < line.size() && is_horizontal_space(line[begin]))
        {
            begin++;
        }
        std::size_t end = begin;
        while (end < line.size() && !is_horizontal_space(line[end]) &&
               line[end] != '\r')
        {
            end++;
        }
        std::string token = line.substr(begin, end - begin);
        if (title_time_found || !Gromacs_Is_Decimal_Number_Token(token))
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, 1,
                "invalid or repeated time field in GROMACS gro title");
        }
        try
        {
            std::size_t parsed_characters = 0;
            double start_time = std::stod(token, &parsed_characters);
            if (parsed_characters != token.size() ||
                !Gromacs_Is_Finite(start_time))
            {
                Gromacs_Throw_Gro_Error(
                    controller, error_by, gro_path, 1,
                    "invalid time field in GROMACS gro title");
            }
            gro_start_time = start_time;
        }
        catch (const std::exception&)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, 1,
                "time field is out of range in GROMACS gro title");
        }
        title_time_found = true;
        i = end;
    }

    if (!std::getline(fin, line))
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, 2,
                                "missing atom count in GROMACS gro file");
    }
    std::string atom_count_token = Gromacs_Trim(line);
    int atom_numbers = 0;
    if (!Gromacs_Is_Integer_Token(atom_count_token))
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, 2,
                                "invalid atom count in GROMACS gro file");
    }
    try
    {
        atom_numbers = std::stoi(atom_count_token);
    }
    catch (const std::exception&)
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, 2,
                                "atom count is out of range in GROMACS gro "
                                "file");
    }
    if (atom_numbers < 0)
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, 2,
                                "negative atom count in GROMACS gro file");
    }

    const Gromacs_Expanded_System_Size expanded_size =
        Gromacs_Check_Expanded_System_Size(topology, atom_numbers, gro_path,
                                           controller, error_by);

    struct Expected_Gro_Atom
    {
        const Gromacs_Molecule_Atom* atom = NULL;
        bool require_atom_name = true;
    };
    std::unordered_map<std::string, std::set<std::string>>
        residue_topology_identities;
    for (const Gromacs_System_Molecule& system_molecule :
         topology.system_molecules)
    {
        const Gromacs_Molecule& molecule =
            topology.molecules.at(system_molecule.name);
        for (std::size_t atom_index = 0; atom_index < molecule.atoms.size();
             atom_index++)
        {
            const Gromacs_Molecule_Atom& atom = molecule.atoms[atom_index];
            residue_topology_identities[atom.residue.substr(0, 5)].insert(
                system_molecule.name + "#" + std::to_string(atom_index));
        }
    }
    const std::size_t atom_count = static_cast<std::size_t>(atom_numbers);
    const std::size_t coordinate_count = 3 * atom_count;
    std::vector<float> staged_coordinates;
    std::vector<float> staged_velocities;

    std::size_t expected_system_index = 0;
    int expected_copy_index = 0;
    std::size_t expected_atom_index = 0;
    auto next_expected_atom = [&](std::size_t line_number)
    {
        while (expected_system_index < topology.system_molecules.size() &&
               expected_copy_index >=
                   topology.system_molecules[expected_system_index].count)
        {
            expected_system_index++;
            expected_copy_index = 0;
            expected_atom_index = 0;
        }
        if (expected_system_index == topology.system_molecules.size())
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "internal GROMACS topology/GRO atom cursor mismatch");
            return Expected_Gro_Atom{};
        }

        const Gromacs_System_Molecule& system_molecule =
            topology.system_molecules[expected_system_index];
        const Gromacs_Molecule& molecule =
            topology.molecules.at(system_molecule.name);
        const Gromacs_Molecule_Atom& atom =
            molecule.atoms[expected_atom_index];
        // Single-atom molecule types commonly use a chemical atom name in the
        // topology (SOD/CLA) and an element-style display name in GRO (NA/CL).
        // Permit that cosmetic difference only when the truncated residue
        // name identifies one topology atom identity.  Repeated copies of the
        // same molecule type remain unambiguous, while two ion types both
        // named ION do not silently accept swapped GRO records.
        const std::string residue_name = atom.residue.substr(0, 5);
        const bool unambiguous_single_atom =
            molecule.atoms.size() == 1 &&
            residue_topology_identities.at(residue_name).size() == 1;
        Expected_Gro_Atom result{&atom, !unambiguous_single_atom};

        expected_atom_index++;
        if (expected_atom_index == molecule.atoms.size())
        {
            expected_atom_index = 0;
            expected_copy_index++;
        }
        return result;
    };

    auto parse_serial = [&](const std::string& field,
                            const std::string& field_name,
                            std::size_t line_number)
    {
        std::string token = Gromacs_Trim(field);
        if (!Gromacs_Is_Integer_Token(token))
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "invalid " + field_name + " in GROMACS gro atom record");
        }
        int value = -1;
        try
        {
            std::size_t parsed_characters = 0;
            value = std::stoi(token, &parsed_characters);
            if (parsed_characters != token.size())
            {
                Gromacs_Throw_Gro_Error(
                    controller, error_by, gro_path, line_number,
                    "invalid " + field_name + " in GROMACS gro atom record");
            }
        }
        catch (const std::exception&)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                field_name + " is out of range in GROMACS gro atom record");
        }
        // GROMACS writes these five-column serials modulo 100000.  Zero is
        // therefore a valid wrapped serial, while negative values are not.
        if (value < 0 || value > 99999)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                field_name + " is out of range in GROMACS gro atom record");
        }
        return value;
    };
    auto parse_fixed_real =
        [&](const std::string& field, std::size_t decimal_offset,
            const std::string& field_name, std::size_t line_number)
    {
        if (decimal_offset >= field.size() || field[decimal_offset] != '.')
        {
            Gromacs_Throw_Gro_Error(controller, error_by, gro_path, line_number,
                                    "inconsistent fixed-width " + field_name +
                                        " in GROMACS gro atom record");
        }
        std::string token = Gromacs_Trim(field);
        if (!Gromacs_Is_Decimal_Number_Token(token))
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "invalid " + field_name + " in GROMACS gro atom record");
        }
        try
        {
            std::size_t parsed_characters = 0;
            double value = std::stod(token, &parsed_characters);
            if (parsed_characters != token.size() || !Gromacs_Is_Finite(value))
            {
                Gromacs_Throw_Gro_Error(
                    controller, error_by, gro_path, line_number,
                    "non-finite " + field_name + " in GROMACS gro atom record");
            }
            return value;
        }
        catch (const std::exception&)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                field_name + " is out of range in GROMACS gro atom record");
        }
        return 0.0;
    };
    auto checked_float = [&](double value, double scale,
                             const std::string& field_name,
                             std::size_t line_number)
    {
        double converted = value * scale;
        if (!Gromacs_Is_Finite(converted) ||
            std::fabs(converted) >
                static_cast<double>(std::numeric_limits<float>::max()))
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                field_name + " cannot be represented by SPONGE");
        }
        float stored = static_cast<float>(converted);
        if (!Gromacs_Is_Finite(stored) || (converted != 0.0 && stored == 0.0f))
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                field_name + " cannot be represented by SPONGE");
        }
        if (!Float_Memory_Is_Zero_Or_Normal(&stored))
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                field_name +
                    " converts to a subnormal float; SPONGE requires a "
                    "finite zero or normal value for consistent FTZ "
                    "behavior");
        }
        return stored;
    };

    std::size_t coordinate_width = 0;
    bool velocity_layout_known = false;
    bool records_have_velocity = false;
    for (int i = 0; i < atom_numbers; i++)
    {
        std::size_t line_number = static_cast<std::size_t>(i) + 3;
        if (!std::getline(fin, line) || line.size() < 20)
        {
            Gromacs_Throw_Gro_Error(controller, error_by, gro_path, line_number,
                                    "invalid atom record in GROMACS gro file");
        }
        if (coordinate_width == 0)
        {
            std::size_t first_decimal = line.find('.', 20);
            std::size_t second_decimal =
                first_decimal == std::string::npos
                    ? std::string::npos
                    : line.find('.', first_decimal + 1);
            std::size_t third_decimal =
                second_decimal == std::string::npos
                    ? std::string::npos
                    : line.find('.', second_decimal + 1);
            if (first_decimal != 24 || second_decimal == std::string::npos ||
                third_decimal == std::string::npos ||
                second_decimal <= first_decimal ||
                third_decimal - second_decimal !=
                    second_decimal - first_decimal ||
                second_decimal - first_decimal < 6)
            {
                Gromacs_Throw_Gro_Error(
                    controller, error_by, gro_path, line_number,
                    "cannot infer a valid coordinate field width from the "
                    "first GROMACS gro atom record");
            }
            coordinate_width = second_decimal - first_decimal;
        }

        if (coordinate_width >
            (std::numeric_limits<std::size_t>::max() - 20) / 3)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "fixed-width coordinate layout overflows the host size type");
        }
        std::size_t coordinate_end = 20 + 3 * coordinate_width;
        if (line.size() < coordinate_end)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "truncated fixed-width coordinate record in GROMACS gro "
                "file");
        }

        parse_serial(line.substr(0, 5), "residue serial", line_number);
        parse_serial(line.substr(15, 5), "atom serial", line_number);
        std::string residue_name = Gromacs_Trim(line.substr(5, 5));
        std::string atom_name = Gromacs_Trim(line.substr(10, 5));
        const Expected_Gro_Atom expected_record =
            next_expected_atom(line_number);
        const Gromacs_Molecule_Atom& expected = *expected_record.atom;
        std::string expected_residue = expected.residue.substr(0, 5);
        std::string expected_atom = expected.atom.substr(0, 5);
        if (residue_name != expected_residue ||
            (expected_record.require_atom_name && atom_name != expected_atom))
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "GROMACS gro residue/atom name order does not match the "
                "expanded topology (expected '" +
                    expected_residue + "/" + expected_atom + "', found '" +
                    residue_name + "/" + atom_name + "')");
        }

        std::array<float, 3> coordinate_values = {};
        for (int axis = 0; axis < 3; axis++)
        {
            std::string field = line.substr(
                20 + static_cast<std::size_t>(axis) * coordinate_width,
                coordinate_width);
            double value =
                parse_fixed_real(field, 4, "coordinate field", line_number);
            coordinate_values[axis] =
                checked_float(value, 10.0, "coordinate field", line_number);
        }

        bool has_velocity = false;
        std::array<float, 3> velocity_values = {};
        if (line.size() > coordinate_end &&
            !Gromacs_Trim(line.substr(coordinate_end)).empty())
        {
            if (coordinate_width >
                (std::numeric_limits<std::size_t>::max() - coordinate_end) / 3)
            {
                Gromacs_Throw_Gro_Error(
                    controller, error_by, gro_path, line_number,
                    "fixed-width velocity layout overflows the host size "
                    "type");
            }
            std::size_t velocity_end = coordinate_end + 3 * coordinate_width;
            if (line.size() < velocity_end)
            {
                Gromacs_Throw_Gro_Error(
                    controller, error_by, gro_path, line_number,
                    "incomplete velocity fields in GROMACS gro atom record");
            }
            if (!Gromacs_Trim(line.substr(velocity_end)).empty())
            {
                Gromacs_Throw_Gro_Error(
                    controller, error_by, gro_path, line_number,
                    "unexpected trailing data in GROMACS gro atom record");
            }
            has_velocity = true;
            for (int axis = 0; axis < 3; axis++)
            {
                std::string field = line.substr(
                    coordinate_end +
                        static_cast<std::size_t>(axis) * coordinate_width,
                    coordinate_width);
                double value =
                    parse_fixed_real(field, 3, "velocity field", line_number);
                velocity_values[axis] = checked_float(
                    value, 10.0 / static_cast<double>(CONSTANT_TIME_CONVERTION),
                    "velocity field", line_number);
            }
        }
        if (!velocity_layout_known)
        {
            velocity_layout_known = true;
            records_have_velocity = has_velocity;
        }
        else if (records_have_velocity != has_velocity)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "mixed atom records with and without velocities in GROMACS "
                "gro file");
        }
        if (staged_coordinates.size() >
                staged_coordinates.max_size() - coordinate_values.size() ||
            staged_velocities.size() >
                staged_velocities.max_size() - velocity_values.size())
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "GROMACS gro coordinate staging exceeds host container "
                "capacity");
        }
        try
        {
            staged_coordinates.insert(staged_coordinates.end(),
                                      coordinate_values.begin(),
                                      coordinate_values.end());
            staged_velocities.insert(staged_velocities.end(),
                                     velocity_values.begin(),
                                     velocity_values.end());
        }
        catch (const std::bad_alloc&)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "GROMACS gro coordinate staging exceeds available host "
                "storage");
        }
        catch (const std::length_error&)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, line_number,
                "GROMACS gro coordinate staging exceeds host container "
                "capacity");
        }
    }

    if (staged_coordinates.size() != coordinate_count ||
        staged_velocities.size() != coordinate_count)
    {
        Gromacs_Throw_Gro_Error(
            controller, error_by, gro_path, 2,
            "internal checked GROMACS coordinate staging size mismatch");
    }

    std::size_t box_line_number = static_cast<std::size_t>(atom_numbers) + 3;
    if (!std::getline(fin, line))
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, box_line_number,
                                "missing box line in GROMACS gro file");
    }
    std::vector<std::string> tokens = Gromacs_Split(line);
    if (tokens.size() != 3 && tokens.size() != 9)
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, box_line_number,
                                "unsupported box line in GROMACS gro file");
    }

    double values[9] = {0.0};
    bool invalid_numeric_field = false;
    try
    {
        for (std::size_t i = 0; i < tokens.size(); i++)
        {
            if (!Gromacs_Is_Decimal_Number_Token(tokens[i]))
            {
                invalid_numeric_field = true;
                continue;
            }
            std::size_t parsed_characters = 0;
            values[i] = std::stod(tokens[i], &parsed_characters);
            if (parsed_characters != tokens[i].size() ||
                !Gromacs_Is_Finite(values[i]))
            {
                invalid_numeric_field = true;
            }
        }
    }
    catch (const std::exception&)
    {
        invalid_numeric_field = true;
    }
    if (invalid_numeric_field)
    {
        Gromacs_Throw_Gro_Error(
            controller, error_by, gro_path, box_line_number,
            "invalid numeric field in GROMACS gro box line");
    }

    // GRO stores triclinic boxes as
    // v1x v2y v3z v1y v1z v2x v2z v3x v3y.
    double vectors[3][3] = {
        {values[0], 0.0, 0.0}, {0.0, values[1], 0.0}, {0.0, 0.0, values[2]}};
    if (tokens.size() == 9)
    {
        vectors[0][1] = values[3];
        vectors[0][2] = values[4];
        vectors[1][0] = values[5];
        vectors[1][2] = values[6];
        vectors[2][0] = values[7];
        vectors[2][1] = values[8];

        // GROMACS keeps GRO triclinic cells in a canonical orientation:
        // v1y=v1z=v2z=0.  SPONGE stores lengths and angles and reconstructs
        // that same orientation.  Accepting an arbitrarily rotated cell here
        // without rotating all coordinates would silently change PBC.
        if (values[3] != 0.0 || values[4] != 0.0 || values[6] != 0.0)
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, box_line_number,
                "invalid non-canonical GROMACS triclinic box "
                "orientation");
        }
    }

    auto dot = [&](int a, int b)
    {
        return vectors[a][0] * vectors[b][0] + vectors[a][1] * vectors[b][1] +
               vectors[a][2] * vectors[b][2];
    };
    double lengths[3] = {std::sqrt(dot(0, 0)), std::sqrt(dot(1, 1)),
                         std::sqrt(dot(2, 2))};
    double determinant =
        vectors[0][0] *
            (vectors[1][1] * vectors[2][2] - vectors[1][2] * vectors[2][1]) -
        vectors[0][1] *
            (vectors[1][0] * vectors[2][2] - vectors[1][2] * vectors[2][0]) +
        vectors[0][2] *
            (vectors[1][0] * vectors[2][1] - vectors[1][1] * vectors[2][0]);
    if (values[0] <= 0.0 || values[1] <= 0.0 || values[2] <= 0.0 ||
        !Gromacs_Is_Finite(lengths[0]) || lengths[0] <= 0.0 ||
        !Gromacs_Is_Finite(lengths[1]) || lengths[1] <= 0.0 ||
        !Gromacs_Is_Finite(lengths[2]) || lengths[2] <= 0.0 ||
        !Gromacs_Is_Finite(determinant) || determinant <= 0.0)
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path, box_line_number,
                                "invalid GROMACS gro box geometry");
    }

    auto angle_degrees = [&](int a, int b)
    {
        const double cross_x =
            vectors[a][1] * vectors[b][2] - vectors[a][2] * vectors[b][1];
        const double cross_y =
            vectors[a][2] * vectors[b][0] - vectors[a][0] * vectors[b][2];
        const double cross_z =
            vectors[a][0] * vectors[b][1] - vectors[a][1] * vectors[b][0];
        const double cross_length = std::sqrt(
            cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
        return std::atan2(cross_length, dot(a, b)) * 180.0 /
               static_cast<double>(Gromacs_Pi);
    };
    float cell_vectors[3][3] = {};
    for (int vector_index = 0; vector_index < 3; vector_index++)
    {
        for (int component = 0; component < 3; component++)
        {
            cell_vectors[vector_index][component] =
                checked_float(vectors[vector_index][component], 10.0,
                              "box vector component", box_line_number);
        }
    }
    const std::array<float, 3> staged_box_length = {
        checked_float(lengths[0], 10.0, "box length", box_line_number),
        checked_float(lengths[1], 10.0, "box length", box_line_number),
        checked_float(lengths[2], 10.0, "box length", box_line_number)};
    const std::array<float, 3> staged_box_angle = {
        checked_float(angle_degrees(1, 2), 1.0, "box angle", box_line_number),
        checked_float(angle_degrees(0, 2), 1.0, "box angle", box_line_number),
        checked_float(angle_degrees(0, 1), 1.0, "box angle", box_line_number)};

    LTMatrix3 cell(cell_vectors[0][0], cell_vectors[1][0], cell_vectors[1][1],
                   cell_vectors[2][0], cell_vectors[2][1], cell_vectors[2][2]);
    LTMatrix3 inverse_cell = inv(cell);
    auto finite_matrix = [&](const LTMatrix3& matrix)
    {
        return Gromacs_Is_Finite(matrix.a11) && Gromacs_Is_Finite(matrix.a21) &&
               Gromacs_Is_Finite(matrix.a22) && Gromacs_Is_Finite(matrix.a31) &&
               Gromacs_Is_Finite(matrix.a32) && Gromacs_Is_Finite(matrix.a33) &&
               Float_Memory_Is_Zero_Or_Normal(&matrix.a11) &&
               Float_Memory_Is_Zero_Or_Normal(&matrix.a21) &&
               Float_Memory_Is_Zero_Or_Normal(&matrix.a22) &&
               Float_Memory_Is_Zero_Or_Normal(&matrix.a31) &&
               Float_Memory_Is_Zero_Or_Normal(&matrix.a32) &&
               Float_Memory_Is_Zero_Or_Normal(&matrix.a33);
    };
    if (cell.a11 <= 0.0f || cell.a22 <= 0.0f || cell.a33 <= 0.0f ||
        !finite_matrix(cell) || !finite_matrix(inverse_cell))
    {
        Gromacs_Throw_Gro_Error(
            controller, error_by, gro_path, box_line_number,
            "GROMACS gro box cannot be represented as a finite reversible "
            "SPONGE cell");
    }
    std::size_t trailing_line_number = box_line_number;
    while (std::getline(fin, line))
    {
        trailing_line_number++;
        if (!Gromacs_Trim(line).empty())
        {
            Gromacs_Throw_Gro_Error(
                controller, error_by, gro_path, trailing_line_number,
                "unexpected trailing data or an additional frame in GROMACS "
                "gro file");
        }
    }
    if (fin.bad())
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path,
                                trailing_line_number,
                                "failed while reading GROMACS gro file");
    }
    fin.clear();
    fin.close();
    if (fin.fail())
    {
        Gromacs_Throw_Gro_Error(controller, error_by, gro_path,
                                trailing_line_number,
                                "failed while closing GROMACS gro file");
    }

    // Publish only after the complete GRO frame (including the box and EOF)
    // has been parsed and validated.  A hostile count paired with a truncated
    // file therefore cannot trigger count-driven topology instantiation.
    Gromacs_Instantiate_System(topology, system, controller, expanded_size);
    Load_Ensure_Atom_Numbers(system, atom_numbers, controller, error_by);
    system->start_time = gro_start_time;
    system->atoms.coordinate = std::move(staged_coordinates);
    system->atoms.velocity = std::move(staged_velocities);
    system->box.box_length = {staged_box_length[0], staged_box_length[1],
                              staged_box_length[2]};
    system->box.box_angle = {staged_box_angle[0], staged_box_angle[1],
                             staged_box_angle[2]};
    Gromacs_Refresh_Initial_Virtual_Sites(topology, system, cell, inverse_cell,
                                          gro_path, controller, error_by);
}

static void Gromacs_Instantiate_System_Impl(
    Gromacs_Topology& topology, System* system, CONTROLLER* controller,
    const Gromacs_Expanded_System_Size& expanded_size)
{
    const char* error_by = "Xponge::Load_Gromacs_Inputs";
    system->source = InputSource::kGromacs;
    system->start_time = 0.0;
    system->atoms.mass.clear();
    system->atoms.charge.clear();
    system->residues.atom_numbers.clear();
    system->exclusions.excluded_atoms.clear();
    system->generalized_born = GeneralizedBorn{};
    system->virtual_atoms = VirtualAtoms{};
    Load_Reset_Classical_Force_Field(&system->classical_force_field);

    auto checked_parameter = [&](double value, double scale,
                                 const Gromacs_Source_Reference& source,
                                 const std::string& description)
    {
        const double converted = value * scale;
        const double maximum_float =
            static_cast<double>(std::numeric_limits<float>::max());
        if (!Gromacs_Is_Finite(value) || !Gromacs_Is_Finite(scale) ||
            !Gromacs_Is_Finite(converted) ||
            std::fabs(converted) > maximum_float)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                description + " cannot be represented by SPONGE");
        }
        const float stored = static_cast<float>(converted);
        if (!Gromacs_Is_Finite(stored) || (converted != 0.0 && stored == 0.0f))
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                description + " cannot be represented by SPONGE");
        }
        if (!Float_Memory_Is_Zero_Or_Normal(&stored))
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                description +
                    " converts to a subnormal float; SPONGE requires a "
                    "finite zero or normal value for consistent FTZ "
                    "behavior");
        }
        return stored;
    };

    struct Local_To_Global
    {
        int base = 0;
        int operator[](std::size_t local_index) const
        {
            return base + static_cast<int>(local_index);
        }
    };
    std::unordered_map<std::string, std::vector<int>> virtual_target_records;
    try
    {
        system->atoms.mass.reserve(expanded_size.atom_count);
        system->atoms.charge.reserve(expanded_size.atom_count);
        system->residues.atom_numbers.reserve(expanded_size.residue_count);
        system->virtual_atoms.records.reserve(
            expanded_size.virtual_site_count);
        system->classical_force_field.lj.atom_type.reserve(
            expanded_size.atom_count);
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed to allocate the checked expanded GROMACS "
            "system (%zu atoms, %zu molecule copies)\n",
            expanded_size.atom_count, expanded_size.molecule_copy_count);
    }
    catch (const std::length_error&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe checked expanded GROMACS system exceeds host "
            "container capacity (%zu atoms, %zu molecule copies)\n",
            expanded_size.atom_count, expanded_size.molecule_copy_count);
    }
    for (auto& molecule_item : topology.molecules)
    {
        std::vector<int> target_record = Gromacs_Validate_Virtual_Sites(
            topology, molecule_item.second, controller, error_by);
        Gromacs_Clean_Virtual_Site_Bondeds(topology, &molecule_item.second,
                                           target_record, controller, error_by);
        virtual_target_records.emplace(molecule_item.first,
                                       std::move(target_record));
    }

    std::size_t instantiated_copy_count = 0;
    for (const auto& item : topology.system_molecules)
    {
        auto iter = topology.molecules.find(item.name);
        if (iter == topology.molecules.end())
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, item.source,
                "molecule '" + item.name +
                    "' referenced in [ molecules ] is not defined");
        }
        const Gromacs_Molecule& molecule = iter->second;
        const std::vector<int>& virtual_targets =
            virtual_target_records.at(item.name);
        for (int copy = 0; copy < item.count; copy++)
        {
            const Local_To_Global local_to_global{
                static_cast<int>(system->atoms.mass.size())};
            int current_resnr = std::numeric_limits<int>::min();
            std::string current_residue;
            for (std::size_t i = 0; i < molecule.atoms.size(); i++)
            {
                const Gromacs_Molecule_Atom& atom = molecule.atoms[i];
                auto atom_type_iter = topology.atom_types.find(atom.type);
                if (atom_type_iter == topology.atom_types.end())
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, atom.source,
                        "undefined GROMACS atom type '" + atom.type + "'");
                }
                const Gromacs_Atom_Type& atom_type = atom_type_iter->second;
                if (atom_type.ptype != "A" && virtual_targets[i] < 0)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, atom.source,
                        "unsupported instantiated GROMACS particle type '" +
                            atom_type.ptype +
                            "' without a supported virtual-site definition");
                }
                float mass = atom.has_mass ? atom.mass : atom_type.mass;
                const Gromacs_Source_Reference& mass_source =
                    atom.has_mass ? atom.source : atom_type.source;
                if (virtual_targets[i] >= 0 && mass != 0.0f)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, mass_source,
                        "GROMACS virtual-site target must have zero mass");
                }
                if (virtual_targets[i] < 0 && mass <= 0.0f)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, mass_source,
                        "GROMACS atom mass must be positive");
                }
                if (system->atoms.mass.size() >=
                    static_cast<std::size_t>(std::numeric_limits<int>::max()))
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, atom.source,
                        "expanded GROMACS atom count exceeds the supported "
                        "kernel index range");
                }
                system->atoms.mass.push_back(mass);
                float charge = atom.has_charge ? atom.charge : atom_type.charge;
                system->atoms.charge.push_back(checked_parameter(
                    charge, static_cast<double>(CONSTANT_SPONGE_CHARGE_SCALE),
                    atom.source, "GROMACS atom charge"));
                system->classical_force_field.lj.atom_type.push_back(
                    expanded_size.atom_type_lj_id.at(atom.type));
                if (atom.resnr != current_resnr ||
                    atom.residue != current_residue)
                {
                    system->residues.atom_numbers.push_back(0);
                    current_resnr = atom.resnr;
                    current_residue = atom.residue;
                }
                system->residues.atom_numbers.back() += 1;
            }
            for (const Gromacs_Virtual_Site& virtual_site :
                 molecule.virtual_sites)
            {
                VirtualAtomRecord record;
                record.virtual_atom = local_to_global[virtual_site.site - 1];
                auto global_source = [&](std::size_t index)
                { return local_to_global[virtual_site.from[index] - 1]; };
                auto distance_angstrom = [&](float distance_nm)
                {
                    return checked_parameter(distance_nm, 10.0,
                                             virtual_site.source,
                                             "GROMACS virtual-site distance");
                };

                if (virtual_site.from.size() == 1)
                {
                    // A one-constructing-atom site is an exact coordinate copy.
                    record.type = 1;
                    record.from = {global_source(0), global_source(0)};
                    record.parameter = {0.0f};
                }
                else if (virtual_site.from.size() == 2 &&
                         virtual_site.funct == 1)
                {
                    record.type = 1;
                    record.from = {global_source(0), global_source(1)};
                    record.parameter = {virtual_site.parameters[0]};
                }
                else if (virtual_site.from.size() == 2 &&
                         virtual_site.funct == 2)
                {
                    // SPONGE type 3 uses r1 + d * normalize(r21 + k*r32).
                    // Repeating the second parent and setting k=0 reproduces
                    // GROMACS virtual_sites2 function 2 exactly.
                    record.type = 3;
                    record.from = {global_source(0), global_source(1),
                                   global_source(1)};
                    record.parameter = {
                        distance_angstrom(virtual_site.parameters[0]), 0.0f};
                }
                else if (virtual_site.from.size() == 3 &&
                         virtual_site.funct == 1)
                {
                    record.type = 2;
                    record.from = {global_source(0), global_source(1),
                                   global_source(2)};
                    record.parameter = {virtual_site.parameters[0],
                                        virtual_site.parameters[1]};
                }
                else if (virtual_site.from.size() == 3 &&
                         virtual_site.funct == 2)
                {
                    // GROMACS spells these parameters as (a, d), while the
                    // SPONGE type-3 record stores (d, k).
                    record.type = 3;
                    record.from = {global_source(0), global_source(1),
                                   global_source(2)};
                    record.parameter = {
                        distance_angstrom(virtual_site.parameters[1]),
                        virtual_site.parameters[0]};
                }
                else
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, virtual_site.source,
                        "unsupported GROMACS virtual-site representation");
                }
                system->virtual_atoms.records.push_back(std::move(record));
            }
            instantiated_copy_count++;
        }
    }

    if (system->atoms.mass.size() != expanded_size.atom_count ||
        system->atoms.charge.size() != expanded_size.atom_count ||
        system->residues.atom_numbers.size() != expanded_size.residue_count ||
        system->virtual_atoms.records.size() !=
            expanded_size.virtual_site_count ||
        system->classical_force_field.lj.atom_type.size() !=
            expanded_size.atom_count ||
        instantiated_copy_count != expanded_size.molecule_copy_count)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tinternal checked GROMACS expansion size mismatch\n");
    }

    const int atom_numbers = static_cast<int>(system->atoms.mass.size());
    Xponge::Bonds& bonds = system->classical_force_field.bonds;
    Xponge::DistanceConstraints& constraints =
        system->classical_force_field.constraints;
    Xponge::Angles& angles = system->classical_force_field.angles;
    Xponge::UreyBradley& urey = system->classical_force_field.urey_bradley;
    Xponge::Torsions& dihedrals = system->classical_force_field.dihedrals;
    Xponge::Torsions& impropers = system->classical_force_field.impropers;
    Xponge::NB14& nb14 = system->classical_force_field.nb14;
    Xponge::CMap& cmap = system->classical_force_field.cmap;
    try
    {
        system->exclusions.excluded_atoms.assign(atom_numbers, {});
        bonds.atom_a.reserve(expanded_size.bond_count);
        bonds.atom_b.reserve(expanded_size.bond_count);
        bonds.k.reserve(expanded_size.bond_count);
        bonds.r0.reserve(expanded_size.bond_count);
        constraints.atom_a.reserve(expanded_size.constraint_count);
        constraints.atom_b.reserve(expanded_size.constraint_count);
        constraints.r0.reserve(expanded_size.constraint_count);
        urey.atom_a.reserve(expanded_size.angle_count);
        urey.atom_b.reserve(expanded_size.angle_count);
        urey.atom_c.reserve(expanded_size.angle_count);
        urey.angle_k.reserve(expanded_size.angle_count);
        urey.angle_theta0.reserve(expanded_size.angle_count);
        urey.bond_k.reserve(expanded_size.angle_count);
        urey.bond_r0.reserve(expanded_size.angle_count);
        auto reserve_torsions = [](Xponge::Torsions* torsions,
                                   std::size_t count)
        {
            torsions->atom_a.reserve(count);
            torsions->atom_b.reserve(count);
            torsions->atom_c.reserve(count);
            torsions->atom_d.reserve(count);
            torsions->pk.reserve(count);
            torsions->pn.reserve(count);
            torsions->ipn.reserve(count);
            torsions->gamc.reserve(count);
            torsions->gams.reserve(count);
        };
        reserve_torsions(&dihedrals, expanded_size.proper_dihedral_count);
        reserve_torsions(&impropers, expanded_size.improper_dihedral_count);
        nb14.atom_a.reserve(expanded_size.pair_count);
        nb14.atom_b.reserve(expanded_size.pair_count);
        nb14.A.reserve(expanded_size.pair_count);
        nb14.B.reserve(expanded_size.pair_count);
        nb14.cf_scale_factor.reserve(expanded_size.pair_count);
        cmap.atom_a.reserve(expanded_size.cmap_count);
        cmap.atom_b.reserve(expanded_size.cmap_count);
        cmap.atom_c.reserve(expanded_size.cmap_count);
        cmap.atom_d.reserve(expanded_size.cmap_count);
        cmap.atom_e.reserve(expanded_size.cmap_count);
        cmap.cmap_type.reserve(expanded_size.cmap_count);
        cmap.resolution.resize(expanded_size.cmap_type_count);
        cmap.type_offset.resize(expanded_size.cmap_type_count);
        cmap.grid_value.reserve(expanded_size.cmap_gridpoint_count);
        system->classical_force_field.lj.pair_A.resize(
            expanded_size.lj_pair_count);
        system->classical_force_field.lj.pair_B.resize(
            expanded_size.lj_pair_count);
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tfailed to allocate the checked expanded GROMACS "
            "interaction storage\n");
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe checked expanded GROMACS interaction storage "
            "exceeds host container capacity\n");
    }

    if (topology.cmap_types.size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tGROMACS CMAP type count exceeds the supported "
            "kernel index range\n");
    }
    cmap.unique_type_numbers = static_cast<int>(topology.cmap_types.size());
    cmap.resolution.resize(cmap.unique_type_numbers);
    cmap.type_offset.resize(cmap.unique_type_numbers);
    std::size_t cumulative_gridpoint_numbers = 0;
    const std::size_t maximum_gridpoint_numbers =
        static_cast<std::size_t>(std::numeric_limits<int>::max()) / 16;
    for (int i = 0; i < cmap.unique_type_numbers; i++)
    {
        const Gromacs_CMap_Type& cmap_type = topology.cmap_types[i];
        const std::size_t resolution =
            static_cast<std::size_t>(cmap_type.resolution);
        if (cmap_type.resolution <= 0 ||
            resolution > std::numeric_limits<std::size_t>::max() / resolution)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, cmap_type.source,
                "GROMACS CMAP resolution overflows the grid size");
        }
        const std::size_t type_gridpoint_numbers = resolution * resolution;
        if (type_gridpoint_numbers > maximum_gridpoint_numbers ||
            cumulative_gridpoint_numbers >
                maximum_gridpoint_numbers - type_gridpoint_numbers)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, cmap_type.source,
                "GROMACS CMAP interpolation table exceeds the supported "
                "kernel index range");
        }
        cmap.resolution[i] = cmap_type.resolution;
        cmap.type_offset[i] =
            static_cast<int>(16 * cumulative_gridpoint_numbers);
        cumulative_gridpoint_numbers += type_gridpoint_numbers;
        for (float value : cmap_type.grid)
        {
            cmap.grid_value.push_back(
                checked_parameter(value, 1.0 / 4.184, cmap_type.source,
                                  "GROMACS CMAP grid value"));
        }
    }
    cmap.unique_gridpoint_numbers =
        static_cast<int>(cumulative_gridpoint_numbers);
    if (topology.cmap_types.size() != expanded_size.cmap_type_count ||
        cumulative_gridpoint_numbers != expanded_size.cmap_gridpoint_count)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tinternal checked GROMACS CMAP topology size "
            "mismatch\n");
    }

    auto checked_lj_ab =
        [&](const std::pair<double, double>& c6_c12, double scale,
            const Gromacs_Source_Reference& source, const std::string& context)
    {
        const double pair_a = 12.0 * c6_c12.second * scale;
        const double pair_b = 6.0 * c6_c12.first * scale;
        const double maximum_float =
            static_cast<double>(std::numeric_limits<float>::max());
        if (!Gromacs_Is_Finite(c6_c12.first) ||
            !Gromacs_Is_Finite(c6_c12.second) || !Gromacs_Is_Finite(scale) ||
            !Gromacs_Is_Finite(pair_a) || !Gromacs_Is_Finite(pair_b) ||
            std::fabs(pair_a) > maximum_float ||
            std::fabs(pair_b) > maximum_float)
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                context +
                    " Lennard-Jones parameters cannot be represented by "
                    "SPONGE");
        }
        const float stored_a = static_cast<float>(pair_a);
        const float stored_b = static_cast<float>(pair_b);
        if (!Gromacs_Is_Finite(stored_a) || !Gromacs_Is_Finite(stored_b) ||
            (pair_a != 0.0 && stored_a == 0.0f) ||
            (pair_b != 0.0 && stored_b == 0.0f))
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                context +
                    " Lennard-Jones parameters cannot be represented by "
                    "SPONGE");
        }
        if (!Float_Memory_Is_Zero_Or_Normal(&stored_a) ||
            !Float_Memory_Is_Zero_Or_Normal(&stored_b))
        {
            Gromacs_Throw_Reference_Error(
                controller, error_by, topology, source,
                context +
                    " Lennard-Jones conversion produces a subnormal float; "
                    "SPONGE requires finite zero or normal values for "
                    "consistent FTZ behavior");
        }
        return std::pair<float, float>{stored_a, stored_b};
    };

    const std::vector<std::string>& ordered_types =
        expanded_size.ordered_lj_types;
    const std::size_t atom_type_numbers = ordered_types.size();
    system->classical_force_field.lj.atom_type_numbers =
        static_cast<int>(atom_type_numbers);
    if (atom_type_numbers > 0 &&
        atom_type_numbers + 1 >
            std::numeric_limits<std::size_t>::max() / atom_type_numbers)
    {
        const Gromacs_Atom_Type& last_type =
            topology.atom_types.at(ordered_types.back());
        Gromacs_Throw_Reference_Error(
            controller, error_by, topology, last_type.source,
            "GROMACS Lennard-Jones pair table size overflows the host size "
            "type");
    }
    const std::size_t pair_type_numbers =
        atom_type_numbers * (atom_type_numbers + 1) / 2;
    if (atom_type_numbers != expanded_size.lj_type_count ||
        pair_type_numbers != expanded_size.lj_pair_count)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tinternal checked GROMACS Lennard-Jones table size "
            "mismatch\n");
    }
    if (pair_type_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        const Gromacs_Atom_Type& last_type =
            topology.atom_types.at(ordered_types.back());
        Gromacs_Throw_Reference_Error(
            controller, error_by, topology, last_type.source,
            "GROMACS Lennard-Jones pair table exceeds the supported kernel "
            "index range");
    }
    system->classical_force_field.lj.pair_A.resize(pair_type_numbers);
    system->classical_force_field.lj.pair_B.resize(pair_type_numbers);

    for (std::size_t type_j_index = 0; type_j_index < ordered_types.size();
         type_j_index++)
    {
        for (std::size_t type_i_index = 0; type_i_index <= type_j_index;
             type_i_index++)
        {
            const Gromacs_Atom_Type& type_i =
                topology.atom_types.at(ordered_types[type_i_index]);
            const Gromacs_Atom_Type& type_j =
                topology.atom_types.at(ordered_types[type_j_index]);
            std::pair<double, double> c6_c12 =
                Gromacs_Get_C6_C12(topology.defaults, type_i, type_j);
            const Gromacs_Pair_Type* nonbond_parameter =
                Gromacs_Find_Nonbond_Parameter(topology, type_i.name,
                                               type_j.name, 1);
            if (nonbond_parameter != NULL)
            {
                c6_c12 = Gromacs_Get_C6_C12_From_Pair_Parameters(
                    topology.defaults, nonbond_parameter->parameters);
            }
            const Gromacs_Source_Reference& parameter_source =
                nonbond_parameter == NULL ? type_j.source
                                          : nonbond_parameter->source;
            std::pair<float, float> pair_ab =
                checked_lj_ab(c6_c12, 1.0, parameter_source,
                              "GROMACS nonbonded type pair '" + type_i.name +
                                  "/" + type_j.name + "'");
            const std::size_t pair_id =
                type_j_index * (type_j_index + 1) / 2 + type_i_index;
            system->classical_force_field.lj.pair_A[pair_id] = pair_ab.first;
            system->classical_force_field.lj.pair_B[pair_id] = pair_ab.second;
        }
    }

    std::size_t molecule_atom_offset = 0;
    for (const auto& item : topology.system_molecules)
    {
        const Gromacs_Molecule& molecule = topology.molecules.at(item.name);
        for (int copy = 0; copy < item.count; copy++)
        {
            const Local_To_Global local_to_global{
                static_cast<int>(molecule_atom_offset)};

            auto require_local_atom =
                [&](int local_index, const Gromacs_Source_Reference& source,
                    const std::string& section)
            {
                if (local_index < 0 ||
                    local_index >= static_cast<int>(molecule.atoms.size()))
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, source,
                        "GROMACS [ " + section +
                            " ] atom index is out of range");
                }
            };

            auto append_bond = [&](int ai_local, int aj_local, float k,
                                   float r0)
            {
                bonds.atom_a.push_back(local_to_global[ai_local]);
                bonds.atom_b.push_back(local_to_global[aj_local]);
                bonds.k.push_back(k);
                bonds.r0.push_back(r0);
            };

            auto append_constraint = [&](int ai_local, int aj_local, float r0)
            {
                constraints.atom_a.push_back(local_to_global[ai_local]);
                constraints.atom_b.push_back(local_to_global[aj_local]);
                constraints.r0.push_back(r0);
            };

            for (const Gromacs_Bond& bond : molecule.bonds)
            {
                int ai_local = bond.ai - 1;
                int aj_local = bond.aj - 1;
                require_local_atom(ai_local, bond.source, "bonds");
                require_local_atom(aj_local, bond.source, "bonds");
                const Gromacs_Molecule_Atom& atom_i = molecule.atoms[ai_local];
                const Gromacs_Molecule_Atom& atom_j = molecule.atoms[aj_local];
                if (bond.funct == 5)
                {
                    // A GROMACS connection bond creates the molecular graph
                    // used by nrexcl, but has no force-field interaction.  The
                    // checked per-template exclusion plan already retained
                    // that connectivity.
                    continue;
                }
                const Gromacs_Bond_Type* type = NULL;
                if (bond.parameters.size() >= 2)
                {
                    append_bond(
                        ai_local, aj_local,
                        checked_parameter(bond.parameters[1],
                                          1.0 / (4.184 * 200.0), bond.source,
                                          "GROMACS bond force constant"),
                        checked_parameter(bond.parameters[0], 10.0, bond.source,
                                          "GROMACS bond distance"));
                }
                else
                {
                    type = Gromacs_Find_Bond_Type(
                        topology, topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        bond.funct);
                    if (type == NULL)
                    {
                        Gromacs_Throw_Reference_Error(
                            controller, error_by, topology, bond.source,
                            "failed to find GROMACS bond type");
                    }
                    append_bond(
                        ai_local, aj_local,
                        checked_parameter(type->kb, 1.0 / (4.184 * 200.0),
                                          type->source,
                                          "GROMACS bond force constant"),
                        checked_parameter(type->b0, 10.0, type->source,
                                          "GROMACS bond distance"));
                }
            }

            for (const Gromacs_Settle& settle : molecule.settles)
            {
                int oxygen_local = settle.ow - 1;
                int hydrogen_1_local = oxygen_local + 1;
                int hydrogen_2_local = oxygen_local + 2;
                require_local_atom(oxygen_local, settle.source, "settles");
                require_local_atom(hydrogen_1_local, settle.source, "settles");
                require_local_atom(hydrogen_2_local, settle.source, "settles");
                const float doh = checked_parameter(
                    settle.doh, 10.0, settle.source,
                    "GROMACS SETTLE oxygen-hydrogen distance");
                const float dhh = checked_parameter(
                    settle.dhh, 10.0, settle.source,
                    "GROMACS SETTLE hydrogen-hydrogen distance");
                append_bond(oxygen_local, hydrogen_1_local, 0.0f, doh);
                append_constraint(oxygen_local, hydrogen_1_local, doh);
                append_bond(oxygen_local, hydrogen_2_local, 0.0f, doh);
                append_constraint(oxygen_local, hydrogen_2_local, doh);
                append_bond(hydrogen_1_local, hydrogen_2_local, 0.0f, dhh);
                append_constraint(hydrogen_1_local, hydrogen_2_local, dhh);
            }

            for (const Gromacs_Constraint& constraint : molecule.constraints)
            {
                int ai_local = constraint.ai - 1;
                int aj_local = constraint.aj - 1;
                require_local_atom(ai_local, constraint.source, "constraints");
                require_local_atom(aj_local, constraint.source, "constraints");
                const Gromacs_Molecule_Atom& atom_i = molecule.atoms[ai_local];
                const Gromacs_Molecule_Atom& atom_j = molecule.atoms[aj_local];
                float distance = 0.0f;
                Gromacs_Source_Reference parameter_source = constraint.source;
                if (!constraint.parameters.empty())
                {
                    distance = constraint.parameters[0];
                }
                else
                {
                    const Gromacs_Constraint_Type* type =
                        Gromacs_Find_Constraint_Type(
                            topology,
                            topology.atom_types.at(atom_i.type).bond_type,
                            topology.atom_types.at(atom_j.type).bond_type,
                            constraint.funct);
                    if (type == NULL)
                    {
                        Gromacs_Throw_Reference_Error(
                            controller, error_by, topology, constraint.source,
                            "failed to resolve GROMACS constraint distance");
                    }
                    distance = type->distance;
                    parameter_source = type->source;
                }
                float distance_angstrom =
                    checked_parameter(distance, 10.0, parameter_source,
                                      "GROMACS constraint distance");
                append_bond(ai_local, aj_local, 0.0f, distance_angstrom);
                append_constraint(ai_local, aj_local, distance_angstrom);
            }
            const std::vector<std::pair<int, int>>& exclusion_pairs =
                expanded_size.molecule_exclusion_pairs.at(&molecule);
            try
            {
                for (std::size_t begin = 0; begin < exclusion_pairs.size();)
                {
                    std::size_t end = begin + 1;
                    while (end < exclusion_pairs.size() &&
                           exclusion_pairs[end].first ==
                               exclusion_pairs[begin].first)
                    {
                        end++;
                    }
                    system->exclusions
                        .excluded_atoms[local_to_global
                                            [exclusion_pairs[begin].first]]
                        .reserve(end - begin);
                    begin = end;
                }
                for (const std::pair<int, int>& exclusion : exclusion_pairs)
                {
                    system->exclusions
                        .excluded_atoms[local_to_global[exclusion.first]]
                        .push_back(local_to_global[exclusion.second]);
                }
            }
            catch (const std::bad_alloc&)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, item.source,
                    "failed to allocate checked expanded GROMACS exclusion "
                    "storage");
            }
            catch (const std::length_error&)
            {
                Gromacs_Throw_Reference_Error(
                    controller, error_by, topology, item.source,
                    "checked expanded GROMACS exclusion storage exceeds host "
                    "container capacity");
            }

            for (const Gromacs_CMap& cmap_item : molecule.cmaps)
            {
                int ai_local = cmap_item.ai - 1;
                int aj_local = cmap_item.aj - 1;
                int ak_local = cmap_item.ak - 1;
                int al_local = cmap_item.al - 1;
                int am_local = cmap_item.am - 1;
                require_local_atom(ai_local, cmap_item.source, "cmap");
                require_local_atom(aj_local, cmap_item.source, "cmap");
                require_local_atom(ak_local, cmap_item.source, "cmap");
                require_local_atom(al_local, cmap_item.source, "cmap");
                require_local_atom(am_local, cmap_item.source, "cmap");
                const Gromacs_Molecule_Atom& atom_i = molecule.atoms[ai_local];
                const Gromacs_Molecule_Atom& atom_j = molecule.atoms[aj_local];
                const Gromacs_Molecule_Atom& atom_k = molecule.atoms[ak_local];
                const Gromacs_Molecule_Atom& atom_l = molecule.atoms[al_local];
                const Gromacs_Molecule_Atom& atom_m = molecule.atoms[am_local];
                int cmap_type = Gromacs_Find_CMap_Type(
                    topology, topology.atom_types.at(atom_i.type).bond_type,
                    topology.atom_types.at(atom_j.type).bond_type,
                    topology.atom_types.at(atom_k.type).bond_type,
                    topology.atom_types.at(atom_l.type).bond_type,
                    topology.atom_types.at(atom_m.type).bond_type,
                    cmap_item.funct);
                if (cmap_type < 0)
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, cmap_item.source,
                        "failed to find GROMACS CMAP type");
                }
                cmap.atom_a.push_back(local_to_global[ai_local]);
                cmap.atom_b.push_back(local_to_global[aj_local]);
                cmap.atom_c.push_back(local_to_global[ak_local]);
                cmap.atom_d.push_back(local_to_global[al_local]);
                cmap.atom_e.push_back(local_to_global[am_local]);
                cmap.cmap_type.push_back(cmap_type);
            }

            for (const Gromacs_Angle& angle : molecule.angles)
            {
                int ai_local = angle.ai - 1;
                int aj_local = angle.aj - 1;
                int ak_local = angle.ak - 1;
                require_local_atom(ai_local, angle.source, "angles");
                require_local_atom(aj_local, angle.source, "angles");
                require_local_atom(ak_local, angle.source, "angles");
                const Gromacs_Molecule_Atom& atom_i = molecule.atoms[ai_local];
                const Gromacs_Molecule_Atom& atom_j = molecule.atoms[aj_local];
                const Gromacs_Molecule_Atom& atom_k = molecule.atoms[ak_local];
                float theta0 = 0.0f;
                float angle_k = 0.0f;
                float ub0 = 0.0f;
                float kub = 0.0f;
                if (angle.parameters.size() >= 2)
                {
                    theta0 = checked_parameter(
                        angle.parameters[0],
                        static_cast<double>(Gromacs_Pi) / 180.0, angle.source,
                        "GROMACS angle equilibrium value");
                    angle_k = checked_parameter(
                        angle.parameters[1], 1.0 / (4.184 * 2.0), angle.source,
                        "GROMACS angle force constant");
                    if (angle.funct == 5 && angle.parameters.size() >= 4)
                    {
                        ub0 = checked_parameter(
                            angle.parameters[2], 10.0, angle.source,
                            "GROMACS Urey-Bradley equilibrium distance");
                        kub = checked_parameter(
                            angle.parameters[3], 1.0 / (4.184 * 200.0),
                            angle.source,
                            "GROMACS Urey-Bradley force constant");
                    }
                }
                else
                {
                    const Gromacs_Angle_Type* type = Gromacs_Find_Angle_Type(
                        topology, topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        topology.atom_types.at(atom_k.type).bond_type,
                        angle.funct);
                    if (type == NULL)
                    {
                        Gromacs_Throw_Reference_Error(
                            controller, error_by, topology, angle.source,
                            "failed to find GROMACS angle type");
                    }
                    theta0 = checked_parameter(
                        type->theta0, static_cast<double>(Gromacs_Pi) / 180.0,
                        type->source, "GROMACS angle equilibrium value");
                    angle_k = checked_parameter(type->k, 1.0 / (4.184 * 2.0),
                                                type->source,
                                                "GROMACS angle force constant");
                    if (angle.funct == 5)
                    {
                        ub0 = checked_parameter(
                            type->ub0, 10.0, type->source,
                            "GROMACS Urey-Bradley equilibrium distance");
                        kub = checked_parameter(
                            type->kub, 1.0 / (4.184 * 200.0), type->source,
                            "GROMACS Urey-Bradley force constant");
                    }
                }
                urey.atom_a.push_back(local_to_global[ai_local]);
                urey.atom_b.push_back(local_to_global[aj_local]);
                urey.atom_c.push_back(local_to_global[ak_local]);
                urey.angle_k.push_back(angle_k);
                urey.angle_theta0.push_back(theta0);
                urey.bond_k.push_back(kub);
                urey.bond_r0.push_back(ub0);
            }

            for (const Gromacs_Dihedral& dihedral : molecule.dihedrals)
            {
                int ai_local = dihedral.ai - 1;
                int aj_local = dihedral.aj - 1;
                int ak_local = dihedral.ak - 1;
                int al_local = dihedral.al - 1;
                require_local_atom(ai_local, dihedral.source, "dihedrals");
                require_local_atom(aj_local, dihedral.source, "dihedrals");
                require_local_atom(ak_local, dihedral.source, "dihedrals");
                require_local_atom(al_local, dihedral.source, "dihedrals");
                const Gromacs_Molecule_Atom& atom_i = molecule.atoms[ai_local];
                const Gromacs_Molecule_Atom& atom_j = molecule.atoms[aj_local];
                const Gromacs_Molecule_Atom& atom_k = molecule.atoms[ak_local];
                const Gromacs_Molecule_Atom& atom_l = molecule.atoms[al_local];

                auto append_proper =
                    [&](float phase_deg, float k_kj, int multiplicity,
                        const Gromacs_Source_Reference& parameter_source)
                {
                    if (multiplicity < 0)
                    {
                        multiplicity = -multiplicity;
                        phase_deg = -phase_deg;
                    }
                    dihedrals.atom_a.push_back(local_to_global[ai_local]);
                    dihedrals.atom_b.push_back(local_to_global[aj_local]);
                    dihedrals.atom_c.push_back(local_to_global[ak_local]);
                    dihedrals.atom_d.push_back(local_to_global[al_local]);
                    dihedrals.ipn.push_back(multiplicity);
                    dihedrals.pn.push_back(static_cast<float>(multiplicity));
                    dihedrals.pk.push_back(checked_parameter(
                        k_kj, 1.0 / 4.184, parameter_source,
                        "GROMACS proper dihedral force constant"));
                    float phase = checked_parameter(
                        phase_deg, static_cast<double>(Gromacs_Pi) / 180.0,
                        parameter_source, "GROMACS proper dihedral phase");
                    dihedrals.gamc.push_back(checked_parameter(
                        static_cast<double>(std::cos(phase)) *
                            dihedrals.pk.back(),
                        1.0, parameter_source,
                        "GROMACS proper dihedral cosine coefficient"));
                    dihedrals.gams.push_back(checked_parameter(
                        static_cast<double>(std::sin(phase)) *
                            dihedrals.pk.back(),
                        1.0, parameter_source,
                        "GROMACS proper dihedral sine coefficient"));
                };
                auto append_improper =
                    [&](float phase_deg, float k_kj,
                        const Gromacs_Source_Reference& parameter_source)
                {
                    if (k_kj < 0.0f)
                    {
                        Gromacs_Throw_Reference_Error(
                            controller, error_by, topology, parameter_source,
                            "GROMACS harmonic improper force constant must "
                            "be non-negative");
                    }
                    impropers.atom_a.push_back(local_to_global[ai_local]);
                    impropers.atom_b.push_back(local_to_global[aj_local]);
                    impropers.atom_c.push_back(local_to_global[ak_local]);
                    impropers.atom_d.push_back(local_to_global[al_local]);
                    // GROMACS harmonic impropers use 1/2*k*(phi-phi0)^2,
                    // whereas SPONGE stores the coefficient of
                    // (phi-phi0)^2 directly.
                    impropers.pk.push_back(checked_parameter(
                        k_kj, 1.0 / (4.184 * 2.0), parameter_source,
                        "GROMACS improper dihedral force constant"));
                    impropers.pn.push_back(0.0f);
                    impropers.ipn.push_back(0);
                    impropers.gamc.push_back(checked_parameter(
                        phase_deg, static_cast<double>(Gromacs_Pi) / 180.0,
                        parameter_source, "GROMACS improper dihedral phase"));
                    impropers.gams.push_back(0.0f);
                };
                auto append_cosine_series =
                    [&](const std::array<double, 6>& coefficients,
                        int highest_multiplicity,
                        const std::string& description,
                        const Gromacs_Source_Reference& parameter_source)
                {
                    for (int multiplicity = 1;
                         multiplicity <= highest_multiplicity; multiplicity++)
                    {
                        dihedrals.atom_a.push_back(local_to_global[ai_local]);
                        dihedrals.atom_b.push_back(local_to_global[aj_local]);
                        dihedrals.atom_c.push_back(local_to_global[ak_local]);
                        dihedrals.atom_d.push_back(local_to_global[al_local]);
                        dihedrals.ipn.push_back(multiplicity);
                        dihedrals.pn.push_back(
                            static_cast<float>(multiplicity));
                        dihedrals.pk.push_back(
                            multiplicity == 1
                                ? checked_parameter(
                                      coefficients[0], 1.0 / 4.184,
                                      parameter_source,
                                      description + " constant coefficient")
                                : 0.0f);
                        dihedrals.gamc.push_back(checked_parameter(
                            coefficients[multiplicity], 1.0 / 4.184,
                            parameter_source,
                            description + " cosine coefficient"));
                        dihedrals.gams.push_back(0.0f);
                    }
                };
                auto append_ryckaert_bellemans =
                    [&](const std::vector<float>& parameters,
                        const Gromacs_Source_Reference& parameter_source)
                {
                    const double c0 = parameters[0];
                    const double c1 = parameters[1];
                    const double c2 = parameters[2];
                    const double c3 = parameters[3];
                    const double c4 = parameters[4];
                    const double c5 = parameters[5];
                    // GROMACS defines RB in the polymer convention
                    // psi=phi-pi.  Expanding powers of cos(psi) into the
                    // SPONGE cosine series gives these exact coefficients.
                    std::array<double, 6> coefficients = {
                        c0 + 0.5 * c2 + 0.375 * c4,
                        -c1 - 0.75 * c3 - 0.625 * c5,
                        0.5 * c2 + 0.5 * c4,
                        -0.25 * c3 - 0.3125 * c5,
                        0.125 * c4,
                        -0.0625 * c5};
                    append_cosine_series(coefficients, 5,
                                         "GROMACS Ryckaert-Bellemans",
                                         parameter_source);
                };
                auto append_fourier =
                    [&](const std::vector<float>& parameters,
                        const Gromacs_Source_Reference& parameter_source)
                {
                    std::array<double, 6> coefficients = {};
                    for (std::size_t term = 0; term < 4; term++)
                    {
                        const double half_coefficient =
                            0.5 * static_cast<double>(parameters[term]);
                        coefficients[0] += half_coefficient;
                        coefficients[term + 1] = term % 2 == 0
                                                     ? half_coefficient
                                                     : -half_coefficient;
                    }
                    append_cosine_series(coefficients, 4,
                                         "GROMACS Fourier dihedral",
                                         parameter_source);
                };

                if (!dihedral.parameters.empty())
                {
                    if (dihedral.funct == 2)
                    {
                        append_improper(dihedral.parameters[0],
                                        dihedral.parameters[1],
                                        dihedral.source);
                    }
                    else if (dihedral.funct == 3)
                    {
                        append_ryckaert_bellemans(dihedral.parameters,
                                                  dihedral.source);
                    }
                    else if (dihedral.funct == 5)
                    {
                        append_fourier(dihedral.parameters, dihedral.source);
                    }
                    else
                    {
                        append_proper(
                            dihedral.parameters[0], dihedral.parameters[1],
                            static_cast<int>(
                                std::lround(dihedral.parameters.size() >= 3
                                                ? dihedral.parameters[2]
                                                : 1.0f)),
                            dihedral.source);
                    }
                    continue;
                }

                std::vector<const Gromacs_Dihedral_Type*> types =
                    Gromacs_Find_Dihedral_Types(
                        topology, topology.atom_types.at(atom_i.type).bond_type,
                        topology.atom_types.at(atom_j.type).bond_type,
                        topology.atom_types.at(atom_k.type).bond_type,
                        topology.atom_types.at(atom_l.type).bond_type,
                        dihedral.funct);
                if (types.empty())
                {
                    Gromacs_Throw_Reference_Error(
                        controller, error_by, topology, dihedral.source,
                        "failed to find GROMACS dihedral type");
                }
                for (const Gromacs_Dihedral_Type* type : types)
                {
                    if (dihedral.funct == 2)
                    {
                        append_improper(type->parameters[0],
                                        type->parameters[1], type->source);
                    }
                    else if (dihedral.funct == 3)
                    {
                        append_ryckaert_bellemans(type->parameters,
                                                  type->source);
                    }
                    else if (dihedral.funct == 5)
                    {
                        append_fourier(type->parameters, type->source);
                    }
                    else
                    {
                        append_proper(type->parameters[0], type->parameters[1],
                                      static_cast<int>(std::lround(
                                          type->parameters.size() >= 3
                                              ? type->parameters[2]
                                              : 1.0f)),
                                      type->source);
                    }
                }
            }

            for (const Gromacs_Pair& pair : molecule.pairs)
            {
                int ai_local = pair.ai - 1;
                int aj_local = pair.aj - 1;
                require_local_atom(ai_local, pair.source, "pairs");
                require_local_atom(aj_local, pair.source, "pairs");
                const Gromacs_Molecule_Atom& atom_i = molecule.atoms[ai_local];
                const Gromacs_Molecule_Atom& atom_j = molecule.atoms[aj_local];
                std::pair<double, double> c6_c12{0.0, 0.0};
                double lj_scale = 1.0;
                Gromacs_Source_Reference parameter_source = pair.source;
                if (pair.parameters.size() >= 2)
                {
                    c6_c12 = Gromacs_Get_C6_C12_From_Pair_Parameters(
                        topology.defaults, pair.parameters);
                }
                else
                {
                    const Gromacs_Pair_Type* pair_type = Gromacs_Find_Pair_Type(
                        topology, atom_i.type, atom_j.type, pair.funct);
                    if (pair_type != NULL)
                    {
                        c6_c12 = Gromacs_Get_C6_C12_From_Pair_Parameters(
                            topology.defaults, pair_type->parameters);
                        parameter_source = pair_type->source;
                    }
                    else if (topology.defaults.gen_pairs)
                    {
                        const Gromacs_Atom_Type& type_i =
                            topology.atom_types.at(atom_i.type);
                        const Gromacs_Atom_Type& type_j =
                            topology.atom_types.at(atom_j.type);
                        c6_c12 = Gromacs_Get_C6_C12(topology.defaults, type_i,
                                                   type_j);
                        parameter_source = type_j.source;
                        const Gromacs_Pair_Type* nonbond_parameter =
                            Gromacs_Find_Nonbond_Parameter(
                                topology, atom_i.type, atom_j.type, 1);
                        if (nonbond_parameter != NULL)
                        {
                            c6_c12 = Gromacs_Get_C6_C12_From_Pair_Parameters(
                                topology.defaults,
                                nonbond_parameter->parameters);
                            parameter_source = nonbond_parameter->source;
                        }
                        lj_scale = topology.defaults.fudge_lj;
                    }
                    else
                    {
                        Gromacs_Throw_Reference_Error(
                            controller, error_by, topology, pair.source,
                            "failed to resolve GROMACS pair interaction");
                    }
                }
                std::pair<float, float> pair_ab =
                    checked_lj_ab(c6_c12, lj_scale, parameter_source,
                                  "GROMACS [ pairs ] interaction");
                nb14.atom_a.push_back(local_to_global[ai_local]);
                nb14.atom_b.push_back(local_to_global[aj_local]);
                nb14.A.push_back(pair_ab.first);
                nb14.B.push_back(pair_ab.second);
                nb14.cf_scale_factor.push_back(topology.defaults.fudge_qq);
            }
            molecule_atom_offset += molecule.atoms.size();
        }
    }

    std::size_t actual_exclusion_count = 0;
    for (const std::vector<int>& row :
         system->exclusions.excluded_atoms)
    {
        if (row.size() >
            std::numeric_limits<std::size_t>::max() - actual_exclusion_count)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorBadFileFormat, error_by,
                "Reason:\n\tinternal expanded GROMACS exclusion count "
                "overflow\n");
        }
        actual_exclusion_count += row.size();
    }
    if (molecule_atom_offset != expanded_size.atom_count ||
        bonds.atom_a.size() != expanded_size.bond_count ||
        bonds.atom_b.size() != expanded_size.bond_count ||
        bonds.k.size() != expanded_size.bond_count ||
        bonds.r0.size() != expanded_size.bond_count ||
        constraints.atom_a.size() != expanded_size.constraint_count ||
        constraints.atom_b.size() != expanded_size.constraint_count ||
        constraints.r0.size() != expanded_size.constraint_count ||
        actual_exclusion_count != expanded_size.exclusion_count ||
        cmap.atom_a.size() != expanded_size.cmap_count ||
        cmap.atom_b.size() != expanded_size.cmap_count ||
        cmap.atom_c.size() != expanded_size.cmap_count ||
        cmap.atom_d.size() != expanded_size.cmap_count ||
        cmap.atom_e.size() != expanded_size.cmap_count ||
        cmap.cmap_type.size() != expanded_size.cmap_count ||
        cmap.resolution.size() != expanded_size.cmap_type_count ||
        cmap.type_offset.size() != expanded_size.cmap_type_count ||
        cmap.grid_value.size() != expanded_size.cmap_gridpoint_count ||
        urey.atom_a.size() != expanded_size.angle_count ||
        urey.atom_b.size() != expanded_size.angle_count ||
        urey.atom_c.size() != expanded_size.angle_count ||
        urey.angle_k.size() != expanded_size.angle_count ||
        urey.angle_theta0.size() != expanded_size.angle_count ||
        urey.bond_k.size() != expanded_size.angle_count ||
        urey.bond_r0.size() != expanded_size.angle_count ||
        dihedrals.atom_a.size() != expanded_size.proper_dihedral_count ||
        dihedrals.atom_b.size() != expanded_size.proper_dihedral_count ||
        dihedrals.atom_c.size() != expanded_size.proper_dihedral_count ||
        dihedrals.atom_d.size() != expanded_size.proper_dihedral_count ||
        dihedrals.pk.size() != expanded_size.proper_dihedral_count ||
        dihedrals.pn.size() != expanded_size.proper_dihedral_count ||
        dihedrals.ipn.size() != expanded_size.proper_dihedral_count ||
        dihedrals.gamc.size() != expanded_size.proper_dihedral_count ||
        dihedrals.gams.size() != expanded_size.proper_dihedral_count ||
        impropers.atom_a.size() != expanded_size.improper_dihedral_count ||
        impropers.atom_b.size() != expanded_size.improper_dihedral_count ||
        impropers.atom_c.size() != expanded_size.improper_dihedral_count ||
        impropers.atom_d.size() != expanded_size.improper_dihedral_count ||
        impropers.pk.size() != expanded_size.improper_dihedral_count ||
        impropers.pn.size() != expanded_size.improper_dihedral_count ||
        impropers.ipn.size() != expanded_size.improper_dihedral_count ||
        impropers.gamc.size() != expanded_size.improper_dihedral_count ||
        impropers.gams.size() != expanded_size.improper_dihedral_count ||
        nb14.atom_a.size() != expanded_size.pair_count ||
        nb14.atom_b.size() != expanded_size.pair_count ||
        nb14.A.size() != expanded_size.pair_count ||
        nb14.B.size() != expanded_size.pair_count ||
        nb14.cf_scale_factor.size() != expanded_size.pair_count ||
        system->classical_force_field.lj.atom_type.size() !=
            expanded_size.atom_count ||
        system->classical_force_field.lj.pair_A.size() !=
            expanded_size.lj_pair_count ||
        system->classical_force_field.lj.pair_B.size() !=
            expanded_size.lj_pair_count ||
        system->classical_force_field.lj.atom_type_numbers !=
            static_cast<int>(expanded_size.lj_type_count))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tinternal checked GROMACS interaction expansion size "
            "mismatch\n");
    }
}

static void Gromacs_Instantiate_System(
    Gromacs_Topology& topology, System* system, CONTROLLER* controller,
    const Gromacs_Expanded_System_Size& expanded_size)
{
    try
    {
        Gromacs_Instantiate_System_Impl(topology, system, controller,
                                        expanded_size);
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Load_Gromacs_Inputs",
            "Reason:\n\tfailed to allocate storage while materializing the "
            "checked GROMACS system\n");
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, "Xponge::Load_Gromacs_Inputs",
            "Reason:\n\tthe checked GROMACS system exceeds host container "
            "capacity while being materialized\n");
    }
}

void Load_Gromacs_Inputs(System* system, CONTROLLER* controller)
{
    Load_System_Transaction(
        system, controller, "Xponge::Load_Gromacs_Inputs",
        Load_System_Seed::kEmpty,
        [&](System* staged)
        {
            Gromacs_Topology topology = Gromacs_Parse_Topology(controller);
            Gromacs_Load_Gro(topology, staged, controller);
        });
}

}  // namespace Xponge
