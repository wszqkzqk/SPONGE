#pragma once

void MD_INFORMATION::non_bond_information::Initial(CONTROLLER* controller,
                                                   MD_INFORMATION* md_info)
{
    if (controller[0].Command_Exist("skin"))
    {
        controller->Check_Float(
            "skin", "MD_INFORMATION::non_bond_information::Initial");
        skin = atof(controller[0].Command("skin"));
    }
    else
    {
        skin = 2.0;
    }
    controller->printf("    skin set to %.2f Angstrom\n", skin);

    if (controller[0].Command_Exist("cutoff"))
    {
        controller->Check_Float(
            "cutoff", "MD_INFORMATION::non_bond_information::Initial");
        cutoff = atof(controller[0].Command("cutoff"));
    }
    else
    {
        cutoff = 10.0;
    }
    const double neighbor_radius =
        static_cast<double>(cutoff) + static_cast<double>(skin);
    if (!Float_Memory_Is_Finite(&cutoff) || !(cutoff > 0.0f) ||
        !Float_Memory_Is_Finite(&skin) || skin < 0.0f ||
        !Float_Memory_Is_Normal(&cutoff) ||
        !Float_Memory_Is_Zero_Or_Normal(&skin) ||
        !Double_Memory_Is_Finite(&neighbor_radius) || !(neighbor_radius > 0.0))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "MD_INFORMATION::non_bond_information::Initial",
            "Reason:\n\tcutoff must be finite and positive, skin must be "
            "finite and nonnegative, both must be normal when nonzero, and "
            "their sum must be finite; got "
            "cutoff %.9g and skin %.9g\n",
            cutoff, skin);
        return;
    }
    controller->printf("    cutoff set to %.2f Angstrom\n", cutoff);
    /*===========================
    读取排除表相关信息
    ============================*/
    const int atom_numbers = md_info->atom_numbers;
    const auto& excluded_rows = Xponge::system.exclusions.excluded_atoms;
    if (atom_numbers <= 0 ||
        excluded_rows.size() != static_cast<std::size_t>(atom_numbers))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat,
            "MD_INFORMATION::non_bond_information::Initial",
            "Reason:\n\tthe in-memory exclusion topology must contain "
            "exactly one row for each of %d atoms; got %zu rows\n",
            atom_numbers, excluded_rows.size());
        return;
    }

    if (!excluded_rows.empty())
    {
        controller->printf("    Start reading excluded list from Xponge:\n");
        std::size_t total = 0;
        for (int atom = 0; atom < atom_numbers; atom++)
        {
            const auto& row = excluded_rows[atom];
            if (row.size() >
                    static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
                row.size() >
                    static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                        total)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow,
                    "MD_INFORMATION::non_bond_information::Initial",
                    "Reason:\n\texclusion count overflows the supported int "
                    "range at atom %d\n",
                    atom);
                return;
            }
            int previous = atom;
            for (std::size_t entry = 0; entry < row.size(); entry++)
            {
                const int partner = row[entry];
                if (partner <= atom || partner >= atom_numbers)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat,
                        "MD_INFORMATION::non_bond_information::Initial",
                        "Reason:\n\texclusion row %d entry %zu has partner "
                        "%d; central exclusion topology must be a strict "
                        "upper-triangular CSR over [0, %d)\n",
                        atom, entry, partner, atom_numbers);
                    return;
                }
                if (partner <= previous)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat,
                        "MD_INFORMATION::non_bond_information::Initial",
                        "Reason:\n\texclusion row %d is not strictly "
                        "increasing at entry %zu (previous %d, current %d)\n",
                        atom, entry, previous, partner);
                    return;
                }
                previous = partner;
            }
            total += row.size();
        }
        excluded_atom_numbers = static_cast<int>(total);

        const std::size_t list_capacity = std::max<std::size_t>(1, total);
        if (list_capacity >
            std::numeric_limits<std::size_t>::max() / sizeof(int))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow,
                "MD_INFORMATION::non_bond_information::Initial",
                "Reason:\n\texclusion-list byte size overflows size_t\n");
            return;
        }

        Malloc_Safely((void**)&h_excluded_list_start,
                      sizeof(int) * atom_numbers);
        Malloc_Safely((void**)&h_excluded_numbers, sizeof(int) * atom_numbers);
        Malloc_Safely((void**)&h_excluded_list, sizeof(int) * list_capacity);
        int count = 0;
        for (int i = 0; i < atom_numbers; i++)
        {
            h_excluded_list_start[i] = count;
            h_excluded_numbers[i] = static_cast<int>(excluded_rows[i].size());
            for (int excluded_atom : excluded_rows[i])
            {
                h_excluded_list[count] = excluded_atom;
                count++;
            }
        }
        Device_Malloc_And_Copy_Safely((void**)&d_excluded_list_start,
                                      h_excluded_list_start,
                                      sizeof(int) * atom_numbers);
        Device_Malloc_And_Copy_Safely((void**)&d_excluded_numbers,
                                      h_excluded_numbers,
                                      sizeof(int) * atom_numbers);
        Device_Malloc_And_Copy_Safely((void**)&d_excluded_list, h_excluded_list,
                                      sizeof(int) * list_capacity);
        controller->printf("    End reading excluded list from Xponge\n\n");
    }
    else
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand,
            "MD_INFORMATION::non_bond_information::Initial",
            "Reason:\n\tno exclusion information found in Xponge::system\n");
    }
}

void MD_INFORMATION::non_bond_information::Excluded_List_Reform(
    CONTROLLER* controller, int atom_numbers)
{
    if (atom_numbers <= 0 || excluded_atom_numbers < 0 ||
        h_excluded_list_start == NULL || h_excluded_numbers == NULL ||
        h_excluded_list == NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "MD_INFORMATION::non_bond_information::Excluded_List_Reform",
            "Reason:\n\tthe triangular exclusion CSR is not initialized\n");
        return;
    }

    if (excluded_atom_numbers > std::numeric_limits<int>::max() / 2)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow,
            "MD_INFORMATION::non_bond_information::Excluded_List_Reform",
            "Reason:\n\tsymmetrizing %d triangular exclusions would "
            "overflow the supported int count\n",
            excluded_atom_numbers);
        return;
    }
    const int old_total = excluded_atom_numbers;
    const int new_total = old_total * 2;
    const std::size_t list_capacity =
        std::max<std::size_t>(1, static_cast<std::size_t>(new_total));
    if (list_capacity > std::numeric_limits<std::size_t>::max() / sizeof(int))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow,
            "MD_INFORMATION::non_bond_information::Excluded_List_Reform",
            "Reason:\n\tthe symmetric exclusion-list byte size overflows "
            "size_t\n");
        return;
    }

    std::vector<int> new_counts(static_cast<std::size_t>(atom_numbers), 0);
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
    {
        const int start = h_excluded_list_start[atom_i];
        const int count = h_excluded_numbers[atom_i];
        if (start < 0 || count < 0 || start > old_total - count)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MD_INFORMATION::non_bond_information::Excluded_List_Reform",
                "Reason:\n\texclusion row %d has invalid start/count %d/%d "
                "for total %d\n",
                atom_i, start, count, old_total);
            return;
        }
        int previous = atom_i;
        for (int entry = 0; entry < count; entry++)
        {
            const int atom_j = h_excluded_list[start + entry];
            if (atom_j <= atom_i || atom_j >= atom_numbers ||
                atom_j <= previous)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "MD_INFORMATION::non_bond_information::Excluded_List_"
                    "Reform",
                    "Reason:\n\texclusion row %d entry %d violates the "
                    "validated strict upper-triangular invariant (partner "
                    "%d, previous %d)\n",
                    atom_i, entry, atom_j, previous);
                return;
            }
            previous = atom_j;
            if (new_counts[atom_i] == std::numeric_limits<int>::max() ||
                new_counts[atom_j] == std::numeric_limits<int>::max())
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow,
                    "MD_INFORMATION::non_bond_information::Excluded_List_"
                    "Reform",
                    "Reason:\n\ta symmetric exclusion row overflows int\n");
                return;
            }
            new_counts[atom_i]++;
            new_counts[atom_j]++;
        }
    }

    std::vector<int> new_starts(static_cast<std::size_t>(atom_numbers), 0);
    int prefix = 0;
    for (int atom = 0; atom < atom_numbers; atom++)
    {
        new_starts[atom] = prefix;
        if (new_counts[atom] > new_total - prefix)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MD_INFORMATION::non_bond_information::Excluded_List_Reform",
                "Reason:\n\tsymmetric exclusion row counts do not sum to "
                "the proven doubled total\n");
            return;
        }
        prefix += new_counts[atom];
    }
    if (prefix != new_total)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "MD_INFORMATION::non_bond_information::Excluded_List_Reform",
            "Reason:\n\tsymmetric exclusion row counts sum to %d, expected "
            "%d\n",
            prefix, new_total);
        return;
    }

    int* new_list = NULL;
    Malloc_Safely((void**)&new_list, sizeof(int) * list_capacity);
    std::vector<int> cursors = new_starts;
    for (int atom_i = 0; atom_i < atom_numbers; atom_i++)
    {
        const int start = h_excluded_list_start[atom_i];
        const int count = h_excluded_numbers[atom_i];
        for (int entry = 0; entry < count; entry++)
        {
            const int atom_j = h_excluded_list[start + entry];
            new_list[cursors[atom_i]++] = atom_j;
            new_list[cursors[atom_j]++] = atom_i;
        }
    }

    // 释放host上的旧排除表
    free(h_excluded_list);
    h_excluded_list = new_list;
    memcpy(h_excluded_list_start, new_starts.data(),
           sizeof(int) * static_cast<std::size_t>(atom_numbers));
    memcpy(h_excluded_numbers, new_counts.data(),
           sizeof(int) * static_cast<std::size_t>(atom_numbers));
    excluded_atom_numbers = new_total;

#ifndef USE_CPU
    // 释放device上的旧排除表
    if (d_excluded_list != NULL) deviceFree(d_excluded_list);
    if (d_excluded_list_start != NULL) deviceFree(d_excluded_list_start);
    if (d_excluded_numbers != NULL) deviceFree(d_excluded_numbers);
#endif
    // 重新分配device上的排除表
    Device_Malloc_And_Copy_Safely((void**)&d_excluded_list, h_excluded_list,
                                  sizeof(int) * list_capacity);
    // 重新分配device上的排除表起点
    Device_Malloc_And_Copy_Safely((void**)&d_excluded_list_start,
                                  h_excluded_list_start,
                                  sizeof(int) * atom_numbers);
    // 重新分配device上的排除表每个原子排除数
    Device_Malloc_And_Copy_Safely((void**)&d_excluded_numbers,
                                  h_excluded_numbers,
                                  sizeof(int) * atom_numbers);
}
