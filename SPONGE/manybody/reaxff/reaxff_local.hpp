#pragma once

#include "atom_identity.h"

inline void REAXFF::Validate_Parallel_Layout(CONTROLLER* controller) const
{
    if (!is_initialized) return;
    if (CONTROLLER::PP_MPI_size != 1)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorNotImplemented, "REAXFF::Validate_Parallel_Layout",
            "REAXFF requires one PP rank, but PP_MPI_size is %d.  EEQ is a "
            "globally coupled linear system with a whole-system charge "
            "constraint; solving the current rank-local neighbor matrix and "
            "reducing only rank-local dot products gives different charges "
            "and a different Hamiltonian.  ReaxFF bond-order state is also "
            "not distributed across PP rank boundaries.  Refusing this "
            "configuration instead of running independent local systems.",
            CONTROLLER::PP_MPI_size);
    }
}

inline void REAXFF::Get_Local(CONTROLLER* controller, const int* atom_local,
                              int local_atom_numbers, int ghost_numbers)
{
    if (!is_initialized) return;
    Validate_Parallel_Layout(controller);
    if (CONTROLLER::PP_MPI_size != 1) return;
    if (local_atom_numbers != atom_numbers || ghost_numbers != 0 ||
        atom_local == NULL)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF::Get_Local",
            "The supported single-PP ReaxFF layout requires exactly %d owned "
            "atoms, no ghosts, and a non-null local-to-global map; received "
            "%d owned atoms, %d ghosts, and map %p.",
            atom_numbers, local_atom_numbers, ghost_numbers,
            static_cast<const void*>(atom_local));
        return;
    }

    std::vector<int> host_atom_local(static_cast<std::size_t>(atom_numbers));
    deviceMemcpy(host_atom_local.data(), atom_local,
                 sizeof(int) * static_cast<std::size_t>(atom_numbers),
                 deviceMemcpyDeviceToHost);
    std::vector<unsigned char> seen(static_cast<std::size_t>(atom_numbers), 0);
    for (int local_id = 0; local_id < atom_numbers; local_id++)
    {
        const int global_id = host_atom_local[local_id];
        if (global_id < 0 || global_id >= atom_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "REAXFF::Get_Local",
                "Local atom %d maps to invalid global atom %d for a %d-atom "
                "ReaxFF system.",
                local_id, global_id, atom_numbers);
            return;
        }
        if (seen[global_id])
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "REAXFF::Get_Local",
                "Global atom %d occurs more than once in the single-PP "
                "ReaxFF local-to-global map (second local index %d).",
                global_id, local_id);
            return;
        }
        seen[global_id] = 1;
    }

    auto gather =
        [&](int initialized, const int* global_values, int* local_values)
    {
        if (initialized)
        {
            ReaxFFAtomIdentity::Gather_Int_By_Global_Id(
                local_atom_numbers, atom_local, global_values, local_values);
        }
    };
    gather(eeq.is_initialized, eeq.d_atom_type_global, eeq.d_atom_type);
    gather(bond_order.is_initialized, bond_order.d_atom_type_global,
           bond_order.d_atom_type);
    gather(bond.is_initialized, bond.d_atom_type_global, bond.d_atom_type);
    gather(vdw.is_initialized, vdw.d_atom_type_global, vdw.d_atom_type);
    gather(ovun.is_initialized, ovun.d_atom_type_global, ovun.d_atom_type);
    gather(angle.is_initialized, angle.d_atom_type_global, angle.d_atom_type);
    gather(torsion.is_initialized, torsion.d_atom_type_global,
           torsion.d_atom_type);
    gather(hb.is_initialized, hb.d_atom_type_global, hb.d_atom_type);
    gather(hb.is_initialized, hb.d_is_hydrogen_global, hb.d_is_hydrogen);
}
