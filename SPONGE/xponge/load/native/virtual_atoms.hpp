#pragma once

#include <sstream>

#include "../../ir/virtual_atoms.hpp"
#include "md_core_parse.hpp"

namespace Xponge
{

static void Native_Load_Virtual_Atoms(VirtualAtoms* virtual_atoms,
                                      CONTROLLER* controller,
                                      const char* module_name = "virtual_atom")
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    std::string file_path =
        controller->Original_Command(module_name, "in_file");
    std::size_t line_number = 0;
    auto fail = [&](const std::string& reason)
    {
        std::string message = "Reason:\n\t" + reason + " at " + file_path +
                              ":" + std::to_string(line_number) + "\n";
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                       "Xponge::Native_Load_Virtual_Atoms",
                                       message.c_str());
    };
    std::ifstream input(file_path);
    if (!input.is_open())
    {
        fail("failed to open virtual-atom input file");
    }

    VirtualAtoms parsed_virtual_atoms;
    try
    {
        std::string line;
        while (std::getline(input, line))
        {
            line_number++;
            std::stringstream stream(line);
            VirtualAtomRecord record;
            stream >> std::ws;
            if (stream.eof())
            {
                continue;
            }
            auto read_int = [&](int* value)
            {
                std::string token;
                if (!(stream >> token))
                {
                    return false;
                }
                long long parsed = 0;
                try
                {
                    std::size_t consumed = 0;
                    parsed = std::stoll(token, &consumed, 10);
                    if (consumed != token.size())
                    {
                        return false;
                    }
                }
                catch (const std::invalid_argument&)
                {
                    return false;
                }
                catch (const std::out_of_range&)
                {
                    fail("virtual-atom integer is outside the supported int "
                         "range");
                }
                if (parsed < std::numeric_limits<int>::min() ||
                    parsed > std::numeric_limits<int>::max())
                {
                    fail("virtual-atom integer is outside the supported int "
                         "range");
                }
                *value = static_cast<int>(parsed);
                return true;
            };
            auto read_float = [&](float* value)
            {
                std::string token;
                if (!(stream >> token))
                {
                    return false;
                }
                std::string lowercase = token;
                std::transform(
                    lowercase.begin(), lowercase.end(), lowercase.begin(),
                    [](unsigned char character)
                    { return static_cast<char>(std::tolower(character)); });
                if (lowercase.find("nan") != std::string::npos ||
                    lowercase.find("inf") != std::string::npos)
                {
                    fail("non-finite virtual-atom parameter");
                }
                if (!Native_Core_Is_Strict_Decimal(token))
                {
                    return false;
                }
                double parsed = 0.0;
                try
                {
                    std::size_t consumed = 0;
                    parsed = std::stod(token, &consumed);
                    if (consumed != token.size() ||
                        !Double_Memory_Is_Finite(&parsed))
                    {
                        return false;
                    }
                }
                catch (const std::invalid_argument&)
                {
                    return false;
                }
                catch (const std::out_of_range&)
                {
                    fail("virtual-atom parameter is outside the supported "
                         "finite range");
                }
                const double float_max =
                    static_cast<double>(std::numeric_limits<float>::max());
                if (parsed > float_max || parsed < -float_max)
                {
                    fail("virtual-atom parameter is outside the supported "
                         "finite float range");
                }
                const float stored = static_cast<float>(parsed);
                if (!Float_Memory_Is_Finite(&stored) ||
                    (parsed != 0.0 && stored == 0.0f))
                {
                    fail("virtual-atom parameter is outside the supported "
                         "finite float range");
                }
                if (!Float_Memory_Is_Zero_Or_Normal(&stored))
                {
                    fail("virtual-atom parameter is a subnormal float; native "
                         "float fields require a finite zero or normal value "
                         "for consistent FTZ behavior");
                }
                *value = stored;
                return true;
            };
            if (!read_int(&record.type) || !read_int(&record.virtual_atom))
            {
                fail("invalid virtual-atom record");
            }
            switch (record.type)
            {
                case 0:
                {
                    int from = 0;
                    float h = 0.0f;
                    if (!read_int(&from) || !read_float(&h))
                    {
                        fail("invalid type-0 virtual-atom record");
                    }
                    record.from.push_back(from);
                    record.parameter.push_back(h);
                    break;
                }
                case 1:
                {
                    int from1 = 0, from2 = 0;
                    float a = 0.0f;
                    if (!read_int(&from1) || !read_int(&from2) ||
                        !read_float(&a))
                    {
                        fail("invalid type-1 virtual-atom record");
                    }
                    record.from.push_back(from1);
                    record.from.push_back(from2);
                    record.parameter.push_back(a);
                    break;
                }
                case 2:
                case 3:
                {
                    int from1 = 0, from2 = 0, from3 = 0;
                    float a = 0.0f, b = 0.0f;
                    if (!read_int(&from1) || !read_int(&from2) ||
                        !read_int(&from3) || !read_float(&a) ||
                        !read_float(&b))
                    {
                        fail("invalid type-2/type-3 virtual-atom record");
                    }
                    record.from.push_back(from1);
                    record.from.push_back(from2);
                    record.from.push_back(from3);
                    record.parameter.push_back(a);
                    record.parameter.push_back(b);
                    break;
                }
                case 5:
                {
                    int from1 = 0, from2 = 0, from3 = 0;
                    float d = 0.0f;
                    if (!read_int(&from1) || !read_int(&from2) ||
                        !read_int(&from3) || !read_float(&d))
                    {
                        fail("invalid type-5 virtual-atom record");
                    }
                    record.from.push_back(from1);
                    record.from.push_back(from2);
                    record.from.push_back(from3);
                    record.parameter.push_back(d);
                    break;
                }
                default:
                    fail("unsupported virtual-atom record type " +
                         std::to_string(record.type));
            }
            stream >> std::ws;
            if (!stream.eof())
            {
                fail("unexpected trailing data in virtual-atom record");
            }
            for (float parameter : record.parameter)
            {
                if (!Virtual_Atom_Parameter_Is_Finite(parameter))
                {
                    fail("non-finite virtual-atom parameter");
                }
            }
            if (record.type == 5 && record.parameter[0] < 0.0f)
            {
                fail("negative type-5 virtual-atom distance");
            }
            if (parsed_virtual_atoms.records.size() >=
                static_cast<std::size_t>(std::numeric_limits<int>::max()))
            {
                fail("virtual-atom record count exceeds the supported int "
                     "range");
            }
            parsed_virtual_atoms.records.push_back(std::move(record));
        }
        if (input.bad())
        {
            fail("I/O error while reading virtual-atom records");
        }
        // A successful getline loop terminates with eofbit/failbit set.  Clear
        // that normal EOF state so a subsequent close failure can be tested
        // independently.
        input.clear();
        input.close();
        if (input.fail())
        {
            fail("I/O error while closing the virtual-atom input file");
        }
    }
    catch (const std::length_error&)
    {
        fail("virtual-atom input exceeds the maximum supported container "
             "size");
    }
    catch (const std::bad_alloc&)
    {
        fail("could not allocate storage while reading virtual-atom records");
    }
    *virtual_atoms = std::move(parsed_virtual_atoms);
}

static void Native_Load_Virtual_Atoms(System* system, CONTROLLER* controller)
{
    Native_Load_Virtual_Atoms(&system->virtual_atoms, controller);
}

}  // namespace Xponge
