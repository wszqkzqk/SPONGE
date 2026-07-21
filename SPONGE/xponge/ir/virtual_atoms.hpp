#pragma once

#include <queue>
#include <string>
#include <vector>

#include "../../common.h"
#include "forcefield.h"

namespace Xponge
{

struct VirtualAtomLayout
{
    std::vector<int> atom_levels;
    int max_level = 0;
};

static bool Virtual_Atom_Parameter_Is_Finite(float value)
{
    return Float_Memory_Is_Finite(&value);
}

static bool Validate_And_Build_Virtual_Atom_Layout(
    const std::vector<VirtualAtomRecord>& records, int atom_numbers,
    VirtualAtomLayout* layout, std::string* error)
{
    auto fail = [&](const std::string& reason)
    {
        if (error != nullptr)
        {
            *error = reason;
        }
        return false;
    };

    if (layout == nullptr)
    {
        return fail("internal error: virtual-atom layout output is null");
    }
    layout->atom_levels.clear();
    layout->max_level = 0;
    if (atom_numbers < 0)
    {
        return fail("atom count is negative");
    }
    if (records.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
        records.size() > static_cast<std::size_t>(atom_numbers))
    {
        return fail(
            "virtual-atom record count exceeds the supported atom or "
            "int range");
    }
    if (static_cast<std::size_t>(atom_numbers) >
            std::vector<int>().max_size() ||
        records.size() > std::vector<std::vector<int>>().max_size())
    {
        return fail(
            "virtual-atom validation allocation count exceeds the "
            "platform range");
    }
    layout->atom_levels.assign(static_cast<std::size_t>(atom_numbers), 0);

    std::vector<int> target_record(static_cast<std::size_t>(atom_numbers), -1);
    for (std::size_t record_index = 0; record_index < records.size();
         record_index++)
    {
        const VirtualAtomRecord& record = records[record_index];
        std::size_t expected_from = 0;
        std::size_t expected_parameters = 0;
        switch (record.type)
        {
            case 0:
                expected_from = 1;
                expected_parameters = 1;
                break;
            case 1:
                expected_from = 2;
                expected_parameters = 1;
                break;
            case 2:
            case 3:
                expected_from = 3;
                expected_parameters = 2;
                break;
            case 5:
                expected_from = 3;
                expected_parameters = 1;
                break;
            default:
                return fail("record " + std::to_string(record_index) +
                            " has unsupported type " +
                            std::to_string(record.type));
        }

        if (record.from.size() != expected_from ||
            record.parameter.size() != expected_parameters)
        {
            return fail("record " + std::to_string(record_index) +
                        " has the wrong number of sources or parameters for "
                        "type " +
                        std::to_string(record.type));
        }
        if (record.virtual_atom < 0 || record.virtual_atom >= atom_numbers)
        {
            return fail("record " + std::to_string(record_index) +
                        " has target atom index " +
                        std::to_string(record.virtual_atom) + " outside [0, " +
                        std::to_string(atom_numbers) + ")");
        }
        if (target_record[record.virtual_atom] >= 0)
        {
            return fail("atom " + std::to_string(record.virtual_atom) +
                        " is the target of more than one virtual-atom record");
        }
        target_record[record.virtual_atom] = static_cast<int>(record_index);

        for (int source : record.from)
        {
            if (source < 0 || source >= atom_numbers)
            {
                return fail("record " + std::to_string(record_index) +
                            " has source atom index " + std::to_string(source) +
                            " outside [0, " + std::to_string(atom_numbers) +
                            ")");
            }
            if (source == record.virtual_atom)
            {
                return fail("record " + std::to_string(record_index) +
                            " uses its target atom as a source");
            }
        }
        for (float parameter : record.parameter)
        {
            if (!Virtual_Atom_Parameter_Is_Finite(parameter))
            {
                return fail("record " + std::to_string(record_index) +
                            " has a non-finite parameter");
            }
        }
        if (record.type == 5)
        {
            if (record.parameter[0] < 0.0f)
            {
                return fail("record " + std::to_string(record_index) +
                            " has a negative type-5 distance");
            }
            if (record.from[0] == record.from[1] ||
                record.from[0] == record.from[2] ||
                record.from[1] == record.from[2])
            {
                return fail("record " + std::to_string(record_index) +
                            " has repeated type-5 source atoms");
            }
        }
        if (record.type == 0 &&
            (record.parameter[0] > std::numeric_limits<float>::max() * 0.5f ||
             record.parameter[0] < -std::numeric_limits<float>::max() * 0.5f))
        {
            return fail("record " + std::to_string(record_index) +
                        " has a type-0 reflection height whose doubled value "
                        "is outside the supported finite float range");
        }
    }

    std::vector<int> indegree(records.size(), 0);
    std::vector<int> record_level(records.size(), 1);
    std::vector<std::vector<int>> consumers(records.size());
    for (std::size_t record_index = 0; record_index < records.size();
         record_index++)
    {
        std::vector<int> dependencies;
        for (int source : records[record_index].from)
        {
            int dependency = target_record[source];
            if (dependency < 0)
            {
                continue;
            }
            bool duplicate = false;
            for (int previous : dependencies)
            {
                if (previous == dependency)
                {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate)
            {
                dependencies.push_back(dependency);
                indegree[record_index]++;
                consumers[dependency].push_back(static_cast<int>(record_index));
            }
        }
    }

    std::queue<int> ready;
    for (std::size_t record_index = 0; record_index < records.size();
         record_index++)
    {
        if (indegree[record_index] == 0)
        {
            ready.push(static_cast<int>(record_index));
        }
    }

    std::size_t processed = 0;
    while (!ready.empty())
    {
        int record_index = ready.front();
        ready.pop();
        processed++;

        int target = records[record_index].virtual_atom;
        int level = record_level[record_index];
        layout->atom_levels[target] = level;
        if (level > layout->max_level)
        {
            layout->max_level = level;
        }

        for (int consumer : consumers[record_index])
        {
            if (level == std::numeric_limits<int>::max())
            {
                return fail(
                    "virtual-atom dependency depth exceeds the supported "
                    "int range");
            }
            if (record_level[consumer] < level + 1)
            {
                record_level[consumer] = level + 1;
            }
            indegree[consumer]--;
            if (indegree[consumer] == 0)
            {
                ready.push(consumer);
            }
        }
    }

    if (processed != records.size())
    {
        return fail("virtual-atom source graph contains a dependency cycle");
    }
    return true;
}

}  // namespace Xponge
