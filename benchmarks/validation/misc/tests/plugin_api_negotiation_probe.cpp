#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

#ifndef PROBE_API_VERSION
#error "PROBE_API_VERSION must be defined"
#endif

#ifdef _WIN32
#define PROBE_EXPORT extern "C" __declspec(dllexport)
#else
#define PROBE_EXPORT extern "C" __attribute__((visibility("default")))
#endif

#if PROBE_API_VERSION > 0
struct ProbeApi
{
    std::uint32_t api_version;
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

#if PROBE_API_VERSION >= 3
    int (*get_force_evaluation_commits_sampling_state)();
    int (*get_force_evaluation_is_exact)();
#endif

#if PROBE_API_VERSION >= 4
    void* (*get_local_energy_ptr)();
    void* (*get_local_virial_ptr)();
    int (*get_force_evaluation_needs_energy)();
    int (*get_force_evaluation_needs_virial)();
    void (*report_fatal_error)(const char* source, const char* message);
#endif

#if PROBE_API_VERSION >= 5
    int (*copy_device_buffer)(void* destination, const void* source,
                              std::uint64_t byte_count);
    int (*get_device_id)();
#endif
};

static const ProbeApi* saved_api = nullptr;
static std::uint32_t initial_version = 0;
#else
static int legacy_version_argument = 0;
static bool legacy_initialized = false;
#endif

PROBE_EXPORT std::string Name()
{
#if PROBE_API_VERSION > 0
    return "stable API negotiation probe v" + std::to_string(PROBE_API_VERSION);
#else
    return "legacy API negotiation probe";
#endif
}

PROBE_EXPORT std::string Version() { return "1"; }

PROBE_EXPORT std::string Version_Check(int version)
{
#if PROBE_API_VERSION > 0
    if (version == PROBE_API_VERSION) return {};
    return "probe rejected API version " + std::to_string(version);
#else
    legacy_version_argument = version;
    return {};
#endif
}

#if PROBE_API_VERSION > 0
PROBE_EXPORT void Initial_Stable(const ProbeApi* api)
{
    saved_api = api;
    initial_version = api == nullptr ? 0 : api->api_version;
}
#else
PROBE_EXPORT void Initial(void*, void*, void*, void*, void*, void*)
{
    legacy_initialized = true;
}
#endif

PROBE_EXPORT void After_Initial()
{
#if PROBE_API_VERSION > 0
    bool callbacks_valid =
        saved_api != nullptr && saved_api->get_command != nullptr &&
        saved_api->log_message != nullptr &&
        saved_api->get_mpi_rank != nullptr &&
        saved_api->get_atom_numbers != nullptr &&
        saved_api->get_steps != nullptr &&
        saved_api->get_coordinate_ptr != nullptr &&
        saved_api->get_force_ptr != nullptr &&
        saved_api->get_neighbor_list_max_numbers != nullptr &&
        saved_api->get_neighbor_list_count != nullptr &&
        saved_api->get_neighbor_list_index_ptr != nullptr &&
        saved_api->get_local_atom_numbers != nullptr &&
        saved_api->get_local_ghost_numbers != nullptr &&
        saved_api->get_local_pp_rank != nullptr &&
        saved_api->get_local_max_atom_numbers != nullptr &&
        saved_api->get_atom_local_ptr != nullptr &&
        saved_api->get_atom_local_label_ptr != nullptr &&
        saved_api->get_atom_local_id_ptr != nullptr &&
        saved_api->get_local_coordinate_ptr != nullptr &&
        saved_api->get_local_force_ptr != nullptr;

#if PROBE_API_VERSION >= 3
    callbacks_valid =
        callbacks_valid &&
        saved_api->get_force_evaluation_commits_sampling_state != nullptr &&
        saved_api->get_force_evaluation_is_exact != nullptr;
#endif
#if PROBE_API_VERSION >= 4
    callbacks_valid = callbacks_valid &&
                      saved_api->get_local_energy_ptr != nullptr &&
                      saved_api->get_local_virial_ptr != nullptr &&
                      saved_api->get_force_evaluation_needs_energy != nullptr &&
                      saved_api->get_force_evaluation_needs_virial != nullptr &&
                      saved_api->report_fatal_error != nullptr;
#endif

    int device_id = -1;
    bool copy_valid = true;
#if PROBE_API_VERSION >= 5
    callbacks_valid = callbacks_valid &&
                      saved_api->copy_device_buffer != nullptr &&
                      saved_api->get_device_id != nullptr;
    if (callbacks_valid)
    {
        device_id = saved_api->get_device_id();
        if (saved_api->device_type == 1)
        {
            char bytes[] = "abcdef";
            copy_valid =
                saved_api->copy_device_buffer(bytes + 1, bytes, 5) == 0 &&
                std::strcmp(bytes, "aabcde") == 0 &&
                saved_api->copy_device_buffer(bytes, bytes, sizeof(bytes)) == 0;
        }
    }
#endif

    char filename[64] = {};
    std::snprintf(filename, sizeof(filename), "stable_api_v%d.log",
                  PROBE_API_VERSION);
    if (FILE* output = std::fopen(filename, "w"))
    {
        std::fprintf(
            output,
            "initial=%u after=%u callbacks=%d copy=%d type=%d "
            "device=%d\n",
            initial_version, saved_api == nullptr ? 0 : saved_api->api_version,
            callbacks_valid ? 1 : 0, copy_valid ? 1 : 0,
            saved_api == nullptr ? -1 : saved_api->device_type, device_id);
        std::fclose(output);
    }
#else
    if (FILE* output = std::fopen("legacy_api.log", "w"))
    {
        std::fprintf(output, "argument=%d initialized=%d\n",
                     legacy_version_argument, legacy_initialized ? 1 : 0);
        std::fclose(output);
    }
#endif
}
