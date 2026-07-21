#pragma once

#include <cstdint>

#include "../common.hpp"

namespace Xponge
{

static bool Native_Core_Is_Strict_Decimal(const std::string& token)
{
    if (token.empty()) return false;

    std::size_t position = 0;
    if (token[position] == '+' || token[position] == '-') position++;

    bool has_digit = false;
    while (position < token.size() && token[position] >= '0' &&
           token[position] <= '9')
    {
        has_digit = true;
        position++;
    }
    if (position < token.size() && token[position] == '.')
    {
        position++;
        while (position < token.size() && token[position] >= '0' &&
               token[position] <= '9')
        {
            has_digit = true;
            position++;
        }
    }
    if (!has_digit) return false;

    if (position < token.size() &&
        (token[position] == 'e' || token[position] == 'E'))
    {
        position++;
        if (position < token.size() &&
            (token[position] == '+' || token[position] == '-'))
        {
            position++;
        }
        const std::size_t exponent_begin = position;
        while (position < token.size() && token[position] >= '0' &&
               token[position] <= '9')
        {
            position++;
        }
        if (position == exponent_begin) return false;
    }
    return position == token.size();
}

static bool Native_Core_Decimal_Significand_Is_Zero(const std::string& token)
{
    for (char character : token)
    {
        if (character == 'e' || character == 'E') break;
        if (character >= '1' && character <= '9') return false;
    }
    return true;
}

static bool Native_Core_Double_Is_Zero(double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value));
    memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT64_C(0x7fffffffffffffff)) == 0;
}

static bool Native_Core_Float_Is_Zero(float value)
{
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value));
    memcpy(&bits, &value, sizeof(bits));
    return (bits & UINT32_C(0x7fffffff)) == 0;
}

class Native_Core_Parser
{
   public:
    Native_Core_Parser(const char* input_path, const char* input_name,
                       const char* error_by, CONTROLLER* controller)
        : input_path_(input_path),
          input_name_(input_name),
          error_by_(error_by),
          controller_(controller),
          input_(input_path)
    {
        if (!input_.is_open())
        {
            Fail(spongeErrorOpenFileFailed,
                 "failed to open " + input_name_ + " for reading");
        }
    }

    std::string Read_Token(const std::string& field)
    {
        try
        {
            std::string token;
            if (!(input_ >> token))
            {
                if (input_.bad())
                {
                    Fail(spongeErrorBadFileFormat,
                         "I/O error while reading " + field + " from " +
                             input_name_);
                }
                Fail(spongeErrorBadFileFormat,
                     input_name_ + " is truncated while reading " + field);
            }
            return token;
        }
        catch (const std::length_error&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field +
                     " exceeds the maximum supported token length");
        }
        catch (const std::bad_alloc&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " could not allocate storage while reading " +
                     field);
        }
        return {};
    }

    std::vector<std::string> Read_Line_Tokens(const std::string& field)
    {
        try
        {
            std::string line;
            if (!std::getline(input_, line))
            {
                if (input_.bad())
                {
                    Fail(spongeErrorBadFileFormat,
                         "I/O error while reading " + field + " from " +
                             input_name_);
                }
                Fail(spongeErrorBadFileFormat,
                     input_name_ + " is truncated while reading " + field);
            }

            std::istringstream line_input(line);
            std::vector<std::string> tokens;
            std::string token;
            while (line_input >> token)
            {
                tokens.push_back(token);
            }
            if (line_input.bad())
            {
                Fail(spongeErrorBadFileFormat,
                     "I/O error while parsing " + field + " from " +
                         input_name_);
            }
            return tokens;
        }
        catch (const std::length_error&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field +
                     " exceeds the maximum supported line or token count");
        }
        catch (const std::bad_alloc&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " could not allocate storage while reading " +
                     field);
        }
        return {};
    }

    template <typename T, typename Value>
    void Append(std::vector<T>* values, Value&& value,
                const std::string& field) const
    {
        try
        {
            values->push_back(std::forward<Value>(value));
        }
        catch (const std::length_error&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " exceeds the maximum supported container "
                               "size while storing " +
                     field);
        }
        catch (const std::bad_alloc&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " could not allocate storage for " + field);
        }
    }

    template <typename T>
    void Assign(std::vector<T>* values, std::size_t count, const T& value,
                const std::string& field) const
    {
        try
        {
            values->assign(count, value);
        }
        catch (const std::length_error&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " exceeds the maximum supported container "
                               "size while storing " +
                     field);
        }
        catch (const std::bad_alloc&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " could not allocate storage for " + field);
        }
    }

    int Parse_Int(const std::string& token, const std::string& field) const
    {
        try
        {
            std::size_t consumed = 0;
            const long long parsed = std::stoll(token, &consumed, 10);
            if (consumed != token.size() ||
                parsed <
                    static_cast<long long>(std::numeric_limits<int>::min()) ||
                parsed >
                    static_cast<long long>(std::numeric_limits<int>::max()))
            {
                throw std::out_of_range("not a strict signed integer");
            }
            return static_cast<int>(parsed);
        }
        catch (const std::exception&)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " token '" + token +
                     "' is not a strict signed integer in range");
        }
        return 0;
    }

    int Read_Int(const std::string& field)
    {
        return Parse_Int(Read_Token(field), field);
    }

    double Parse_Double(const std::string& token,
                        const std::string& field) const
    {
        if (!Native_Core_Is_Strict_Decimal(token))
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " token '" + token +
                     "' is not a strict finite decimal");
        }

        char* parsed_end = nullptr;
        const double parsed = strtod(token.c_str(), &parsed_end);
        if (parsed_end != token.c_str() + token.size() ||
            !Double_Memory_Is_Finite(&parsed) ||
            (Native_Core_Double_Is_Zero(parsed) &&
             !Native_Core_Decimal_Significand_Is_Zero(token)))
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " token '" + token +
                     "' is outside the finite double range");
        }
        return parsed;
    }

    double Read_Double(const std::string& field)
    {
        return Parse_Double(Read_Token(field), field);
    }

    float Parse_Float(const std::string& token, const std::string& field) const
    {
        const double parsed = Parse_Double(token, field);
        const double float_max =
            static_cast<double>(std::numeric_limits<float>::max());
        if (parsed > float_max || parsed < -float_max)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " token '" + token +
                     "' is outside the finite float range");
        }

        const float stored = static_cast<float>(parsed);
        if (!Float_Memory_Is_Finite(&stored))
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " token '" + token +
                     "' is outside the finite float range");
        }
        if (!Native_Core_Double_Is_Zero(parsed) &&
            Native_Core_Float_Is_Zero(stored))
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " nonzero " + field + " token '" + token +
                     "' underflows the finite float range");
        }
        if (!Float_Memory_Is_Zero_Or_Normal(&stored))
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " token '" + token +
                     "' is a subnormal float; native float fields require "
                     "a finite zero or normal float for consistent FTZ "
                     "behavior");
        }
        return stored;
    }

    float Read_Float(const std::string& field)
    {
        return Parse_Float(Read_Token(field), field);
    }

    float Checked_Float(double value, const std::string& field) const
    {
        const double float_max =
            static_cast<double>(std::numeric_limits<float>::max());
        if (!Double_Memory_Is_Finite(&value) || value > float_max ||
            value < -float_max)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field +
                     " is outside the finite float range");
        }
        const float stored = static_cast<float>(value);
        if (!Float_Memory_Is_Finite(&stored) ||
            !Float_Memory_Is_Zero_Or_Normal(&stored) ||
            (!Native_Core_Double_Is_Zero(value) &&
             Native_Core_Float_Is_Zero(stored)))
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field +
                     " is not representable as a finite zero or normal "
                     "float on FTZ backends");
        }
        return stored;
    }

    int Validate_Triangular_Type_Count(int type_count,
                                       const std::string& field) const
    {
        if (type_count < 0)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " has a negative " + field);
        }
        const std::uint64_t count = static_cast<std::uint64_t>(type_count);
        const std::uint64_t pair_count = count * (count + 1) / 2;
        const std::vector<float> pair_storage;
        if (pair_count >
                static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
            pair_count > static_cast<std::uint64_t>(pair_storage.max_size()))
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " " + std::to_string(type_count) +
                     " has an unsupported triangular pair count");
        }
        return static_cast<int>(pair_count);
    }

    int Validate_Atom_Count(int atom_count, const std::string& field) const
    {
        if (atom_count < 0)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " has a negative " + field);
        }
        const std::size_t atom_count_size =
            static_cast<std::size_t>(atom_count);
        const std::vector<float> component_storage;
        if (atom_count > std::numeric_limits<int>::max() / 3 ||
            atom_count_size > std::numeric_limits<std::size_t>::max() / 3 ||
            atom_count_size > component_storage.max_size() / 3)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " " + field + " " + std::to_string(atom_count) +
                     " cannot safely represent all 3 * atom count values and "
                     "indices");
        }
        return atom_count;
    }

    void Ensure_Atom_Count_Matches(const System* system, int atom_count,
                                   const std::string& field) const
    {
        const int current_atom_count = Load_Get_Atom_Numbers(system);
        if (current_atom_count < 0)
        {
            Fail(spongeErrorConflictingCommand,
                 "the retained atom arrays have inconsistent, misaligned, "
                 "or unsupported sizes before reading " +
                     field);
        }
        if (current_atom_count > 0 && current_atom_count != atom_count)
        {
            Fail(spongeErrorConflictingCommand,
                 input_name_ + " " + field + " " + std::to_string(atom_count) +
                     " differs from the previously loaded atom count " +
                     std::to_string(current_atom_count));
        }
    }

    void Ensure_End()
    {
        std::string trailing;
        if (input_ >> trailing)
        {
            Fail(spongeErrorBadFileFormat,
                 input_name_ + " has trailing data beginning with '" +
                     trailing + "'");
        }
        if (input_.bad())
        {
            Fail(spongeErrorBadFileFormat,
                 "I/O error while checking the end of " + input_name_);
        }
    }

    void Close()
    {
        // Ensure_End necessarily leaves eofbit/failbit set after its final
        // extraction attempt.  Clear those checked state bits before asking
        // the stream buffer to report an actual close failure.
        input_.clear();
        input_.close();
        if (input_.fail())
        {
            Fail(spongeErrorBadFileFormat,
                 "I/O error while closing " + input_name_);
        }
    }

    void Fail(int error_number, const std::string& reason) const
    {
        const std::string message =
            "Reason:\n\t" + reason + "\n\tInput file: " + input_path_ + "\n";
        controller_->Throw_SPONGE_Error(error_number, error_by_.c_str(),
                                        message.c_str());
    }

   private:
    std::string input_path_;
    std::string input_name_;
    std::string error_by_;
    CONTROLLER* controller_;
    std::ifstream input_;
};

static std::string Native_Core_Entry_Field(const char* field, int entry)
{
    return std::string(field) + " entry " + std::to_string(entry);
}

}  // namespace Xponge
