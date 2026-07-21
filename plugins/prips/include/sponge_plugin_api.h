#pragma once

#include <stdint.h>

extern "C"
{
    enum
    {
        SPONGE_PRIPS_API_VERSION = 5,
    };

    enum SPONGE_PLUGIN_FORCE_CAPABILITY
    {
        // Calculate_Force contributes its complete potential-energy term to
        // the local energy buffer whenever the evaluation requests energy.
        SPONGE_PLUGIN_FORCE_ENERGY_COMPLETE = 1u << 0,
        // Calculate_Force contributes its complete virial term to the local
        // virial buffer whenever the evaluation requests virial.
        SPONGE_PLUGIN_FORCE_VIRIAL_COMPLETE = 1u << 1,
        // The callback has no mutable state that can affect force results.
        SPONGE_PLUGIN_FORCE_PURE = 1u << 2,
        // Mutable evaluation state is isolated by the transaction hooks.
        SPONGE_PLUGIN_FORCE_TRANSACTIONAL = 1u << 3,
    };

    typedef struct SPONGE_PLUGIN_API
    {
        uint32_t api_version;
        int device_type;

        const char* (*get_command)(const char* key);
        void (*log_message)(const char* message);

        int (*get_mpi_rank)();
        int (*get_atom_numbers)();
        int (*get_steps)();
        void* (*get_coordinate_ptr)();
        void* (*get_force_ptr)();

        int (*get_neighbor_list_max_numbers)();
        int (*get_neighbor_list_count)(int atom_index);
        void* (*get_neighbor_list_index_ptr)();

        int (*get_local_atom_numbers)();
        int (*get_local_ghost_numbers)();
        int (*get_local_pp_rank)();
        int (*get_local_max_atom_numbers)();
        void* (*get_atom_local_ptr)();
        void* (*get_atom_local_label_ptr)();
        void* (*get_atom_local_id_ptr)();
        void* (*get_local_coordinate_ptr)();
        void* (*get_local_force_ptr)();

        // API v3 begins here.  Stable API revisions are append-only: never
        // insert, remove, or reorder members above or within an older block.
        // Force callbacks can be invoked for a committed dynamics sample or
        // for a transactional old/trial-state evaluation.  Plugins that own
        // adaptive/history state must only advance it for committed samples.
        int (*get_force_evaluation_commits_sampling_state)();
        int (*get_force_evaluation_is_exact)();

        // API v4 begins here.
        // These are the authoritative buffers for atoms owned by the current
        // PP rank.  Their shapes are (get_local_atom_numbers(),) and
        // (get_local_atom_numbers(), 6), respectively.  Virial elements use
        // LTMatrix3 memory order: a11, a21, a22, a31, a32, a33.
        void* (*get_local_energy_ptr)();
        void* (*get_local_virial_ptr)();

        // Nonzero only when the current Calculate_Force invocation requires
        // the corresponding complete Hamiltonian contribution.
        int (*get_force_evaluation_needs_energy)();
        int (*get_force_evaluation_needs_virial)();

        // Report an unrecoverable plugin contract/runtime error through the
        // core controller so MPI jobs terminate collectively instead of
        // leaving peer ranks blocked in a collective.
        void (*report_fatal_error)(const char* source, const char* message);

        // API v5 begins here.
        // Synchronously copy a contiguous buffer between pointers on the
        // active SPONGE device.  The call does not return until the source is
        // safe for the producer to release and subsequent force kernels can
        // consume the destination.  Returns zero on success.
        int (*copy_device_buffer)(void* destination, const void* source,
                                  uint64_t byte_count);

        // Logical device ordinal used by this SPONGE process.  DLPack views
        // and functional writeback results must declare this exact id.
        int (*get_device_id)();
    } SPONGE_PLUGIN_API;
}
