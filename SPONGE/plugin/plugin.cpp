#include "plugin.h"

#include <cstring>
#include <limits>
#include <memory>

namespace
{
MD_INFORMATION* g_plugin_md_info = NULL;
CONTROLLER* g_plugin_controller = NULL;
NEIGHBOR_LIST* g_plugin_neighbor_list = NULL;
DOMAIN_INFORMATION* g_plugin_domain_info = NULL;
std::vector<std::unique_ptr<SPONGE_PLUGIN_API>> g_prips_apis;
bool g_force_evaluation_commits_sampling_state = true;
bool g_force_evaluation_is_exact = false;
bool g_force_evaluation_needs_energy = false;
bool g_force_evaluation_needs_virial = false;
}  // namespace

std::map<std::string,
         std::function<void(COLLECTIVE_VARIABLE_CONTROLLER*, int, const char*)>>
    SPONGE_PLUGIN::cv_init_functions;
std::map<std::string, std::function<void(int, UNSIGNED_INT_VECTOR*, VECTOR,
                                         VECTOR*, VECTOR, int, int)>>
    SPONGE_PLUGIN::cv_compute_functions;

static std::string DlErrorString()
{
#ifdef _WIN32
    return std::to_string(static_cast<unsigned long>(dlerror()));
#else
    const char* err = dlerror();
    return err == NULL ? std::string() : std::string(err);
#endif
}

static int PluginBackendDeviceType()
{
#ifdef USE_HIP
    return 10;  // kDLROCM
#elif defined(USE_CUDA)
    return 2;  // kDLCUDA
#else
    return 1;  // kDLCPU
#endif
}

namespace
{
int PluginGetDeviceId()
{
#if defined(USE_CUDA) || defined(USE_HIP)
    return g_plugin_controller == NULL ? 0
                                       : g_plugin_controller->working_device;
#else
    // On CPU, CONTROLLER::working_device is an OpenMP thread count rather
    // than a DLPack device ordinal.
    return 0;
#endif
}

const char* PluginGetCommand(const char* key)
{
    if (g_plugin_controller == NULL || !g_plugin_controller->Command_Exist(key))
    {
        return NULL;
    }
    // Plugins commonly use controller values as paths (for example PRIPS'
    // Python script).  Return the complete scalar value rather than the first
    // whitespace-delimited compatibility token.
    return g_plugin_controller->Original_Command(key);
}

void PluginLogMessage(const char* message)
{
    if (g_plugin_controller == NULL || message == NULL) return;
    g_plugin_controller->printf("%s", message);
}

int PluginGetMPIRank() { return CONTROLLER::MPI_rank; }

int PluginGetAtomNumbers()
{
    return g_plugin_md_info == NULL ? 0 : g_plugin_md_info->atom_numbers;
}

int PluginGetSteps()
{
    return g_plugin_md_info == NULL ? 0 : g_plugin_md_info->sys.steps;
}

void* PluginGetCoordinatePtr()
{
    return g_plugin_md_info == NULL ? NULL : g_plugin_md_info->crd;
}

void* PluginGetForcePtr()
{
    return g_plugin_md_info == NULL ? NULL : g_plugin_md_info->frc;
}

int PluginGetNeighborListMaxNumbers()
{
    return g_plugin_neighbor_list == NULL
               ? 0
               : g_plugin_neighbor_list->max_neighbor_numbers;
}

int PluginGetNeighborListCount(int atom_index)
{
    if (g_plugin_neighbor_list == NULL ||
        g_plugin_neighbor_list->h_nl == NULL || g_plugin_md_info == NULL ||
        atom_index < 0 || atom_index >= g_plugin_md_info->atom_numbers)
    {
        return 0;
    }
    return g_plugin_neighbor_list->h_nl[atom_index].atom_numbers;
}

void* PluginGetNeighborListIndexPtr()
{
    if (g_plugin_neighbor_list == NULL || g_plugin_neighbor_list->h_nl == NULL)
    {
        return NULL;
    }
    return g_plugin_neighbor_list->h_nl->atom_serial;
}

int PluginGetLocalAtomNumbers()
{
    return g_plugin_domain_info == NULL ? 0
                                        : g_plugin_domain_info->atom_numbers;
}

int PluginGetLocalGhostNumbers()
{
    return g_plugin_domain_info == NULL ? 0
                                        : g_plugin_domain_info->ghost_numbers;
}

int PluginGetLocalPPRank()
{
    return g_plugin_domain_info == NULL ? 0 : g_plugin_domain_info->pp_rank;
}

int PluginGetLocalMaxAtomNumbers()
{
    return g_plugin_domain_info == NULL
               ? 0
               : g_plugin_domain_info->max_atom_numbers;
}

void* PluginGetAtomLocalPtr()
{
    return g_plugin_domain_info == NULL ? NULL
                                        : g_plugin_domain_info->atom_local;
}

void* PluginGetAtomLocalLabelPtr()
{
    return g_plugin_domain_info == NULL
               ? NULL
               : g_plugin_domain_info->atom_local_label;
}

void* PluginGetAtomLocalIdPtr()
{
    return g_plugin_domain_info == NULL ? NULL
                                        : g_plugin_domain_info->atom_local_id;
}

void* PluginGetLocalCoordinatePtr()
{
    return g_plugin_domain_info == NULL ? NULL : g_plugin_domain_info->crd;
}

void* PluginGetLocalForcePtr()
{
    return g_plugin_domain_info == NULL ? NULL : g_plugin_domain_info->frc;
}

void* PluginGetLocalEnergyPtr()
{
    return g_plugin_domain_info == NULL ? NULL : g_plugin_domain_info->d_energy;
}

void* PluginGetLocalVirialPtr()
{
    return g_plugin_domain_info == NULL ? NULL : g_plugin_domain_info->d_virial;
}

int PluginForceEvaluationCommitsSamplingState()
{
    return g_force_evaluation_commits_sampling_state ? 1 : 0;
}

int PluginForceEvaluationIsExact()
{
    return g_force_evaluation_is_exact ? 1 : 0;
}

int PluginForceEvaluationNeedsEnergy()
{
    return g_force_evaluation_needs_energy ? 1 : 0;
}

int PluginForceEvaluationNeedsVirial()
{
    return g_force_evaluation_needs_virial ? 1 : 0;
}

void PluginReportFatalError(const char* source, const char* message)
{
    if (g_plugin_controller != NULL)
    {
        g_plugin_controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            source == NULL ? "SPONGE plugin" : source,
            message == NULL ? "Reason:\n\tplugin reported a fatal error\n"
                            : message);
    }
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, spongeErrorSimulationBreakDown);
#else
    exit(spongeErrorSimulationBreakDown);
#endif
}

int PluginCopyDeviceBuffer(void* destination, const void* source,
                           uint64_t byte_count)
{
    if (byte_count == 0) return 0;
    if (destination == NULL || source == NULL ||
        byte_count >
            static_cast<uint64_t>(std::numeric_limits<std::size_t>::max()))
    {
        return -1;
    }
    const std::size_t bytes = static_cast<std::size_t>(byte_count);
#ifdef USE_CUDA
    cudaError_t status = cudaSetDevice(PluginGetDeviceId());
    if (status != cudaSuccess) return static_cast<int>(status);
    // Functional-return producers may have used any CUDA stream.  A
    // device-wide barrier before the copy is the conservative protocol that
    // makes their writes visible without importing framework-specific stream
    // handles.  The post-copy barrier keeps the producer object alive until
    // the destination is complete.  The pre-barrier is still required when
    // source == destination (for example an in-place CuPy/Torch update).
    status = cudaDeviceSynchronize();
    if (status != cudaSuccess) return static_cast<int>(status);
    if (destination != source)
    {
        status =
            cudaMemcpy(destination, source, bytes, cudaMemcpyDeviceToDevice);
        if (status != cudaSuccess) return static_cast<int>(status);
    }
    status = cudaDeviceSynchronize();
    return status == cudaSuccess ? 0 : static_cast<int>(status);
#elif defined(USE_HIP)
    hipError_t status = hipSetDevice(PluginGetDeviceId());
    if (status != hipSuccess) return static_cast<int>(status);
    status = hipDeviceSynchronize();
    if (status != hipSuccess) return static_cast<int>(status);
    if (destination != source)
    {
        status = hipMemcpy(destination, source, bytes, hipMemcpyDeviceToDevice);
        if (status != hipSuccess) return static_cast<int>(status);
    }
    status = hipDeviceSynchronize();
    return status == hipSuccess ? 0 : static_cast<int>(status);
#else
    if (destination != source) std::memmove(destination, source, bytes);
    return 0;
#endif
}

const SPONGE_PLUGIN_API* BuildPripsApi(uint32_t negotiated_version)
{
    std::unique_ptr<SPONGE_PLUGIN_API> api(new SPONGE_PLUGIN_API{});
    api->api_version = negotiated_version;
    api->device_type = PluginBackendDeviceType();
    api->get_command = PluginGetCommand;
    api->log_message = PluginLogMessage;
    api->get_mpi_rank = PluginGetMPIRank;
    api->get_atom_numbers = PluginGetAtomNumbers;
    api->get_steps = PluginGetSteps;
    api->get_coordinate_ptr = PluginGetCoordinatePtr;
    api->get_force_ptr = PluginGetForcePtr;
    api->get_neighbor_list_max_numbers = PluginGetNeighborListMaxNumbers;
    api->get_neighbor_list_count = PluginGetNeighborListCount;
    api->get_neighbor_list_index_ptr = PluginGetNeighborListIndexPtr;
    api->get_local_atom_numbers = PluginGetLocalAtomNumbers;
    api->get_local_ghost_numbers = PluginGetLocalGhostNumbers;
    api->get_local_pp_rank = PluginGetLocalPPRank;
    api->get_local_max_atom_numbers = PluginGetLocalMaxAtomNumbers;
    api->get_atom_local_ptr = PluginGetAtomLocalPtr;
    api->get_atom_local_label_ptr = PluginGetAtomLocalLabelPtr;
    api->get_atom_local_id_ptr = PluginGetAtomLocalIdPtr;
    api->get_local_coordinate_ptr = PluginGetLocalCoordinatePtr;
    api->get_local_force_ptr = PluginGetLocalForcePtr;
    api->get_force_evaluation_commits_sampling_state =
        PluginForceEvaluationCommitsSamplingState;
    api->get_force_evaluation_is_exact = PluginForceEvaluationIsExact;
    api->get_local_energy_ptr = PluginGetLocalEnergyPtr;
    api->get_local_virial_ptr = PluginGetLocalVirialPtr;
    api->get_force_evaluation_needs_energy = PluginForceEvaluationNeedsEnergy;
    api->get_force_evaluation_needs_virial = PluginForceEvaluationNeedsVirial;
    api->report_fatal_error = PluginReportFatalError;
    api->copy_device_buffer = PluginCopyDeviceBuffer;
    api->get_device_id = PluginGetDeviceId;
    const SPONGE_PLUGIN_API* result = api.get();
    g_prips_apis.push_back(std::move(api));
    return result;
}

bool UsesMonteCarloBarostat(MD_INFORMATION* md_info, CONTROLLER* controller)
{
    return md_info != NULL && controller != NULL &&
           md_info->mode == md_info->NPT &&
           (controller->Command_Choice("barostat", "monte_carlo_barostat") ||
            controller->Command_Choice("barostat_mode",
                                       "monte_carlo_barostat"));
}

void ValidateForcePluginContract(
    CONTROLLER* controller, const std::string& plugin_name,
    bool has_capability_getter, uint32_t capabilities,
    RuntimeFunction begin_transaction, RuntimeFunction commit_transaction,
    RuntimeFunction rollback_transaction, bool monte_carlo_barostat)
{
    const uint32_t known_capabilities = SPONGE_PLUGIN_FORCE_ENERGY_COMPLETE |
                                        SPONGE_PLUGIN_FORCE_VIRIAL_COMPLETE |
                                        SPONGE_PLUGIN_FORCE_PURE |
                                        SPONGE_PLUGIN_FORCE_TRANSACTIONAL;
    const uint32_t unknown_capabilities = capabilities & ~known_capabilities;
    const bool is_pure = (capabilities & SPONGE_PLUGIN_FORCE_PURE) != 0;
    const bool is_transactional =
        (capabilities & SPONGE_PLUGIN_FORCE_TRANSACTIONAL) != 0;

    if (has_capability_getter && unknown_capabilities != 0)
    {
        std::string reason = "Reason:\n\tCalculate_Force plugin '" +
                             plugin_name +
                             "' declares unknown force capability bits " +
                             std::to_string(unknown_capabilities) + "\n";
        controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                       "SPONGE_PLUGIN::Initial",
                                       reason.c_str());
    }
    if (is_pure && is_transactional)
    {
        std::string reason =
            "Reason:\n\tCalculate_Force plugin '" + plugin_name +
            "' declares both PURE and TRANSACTIONAL; these capabilities are "
            "mutually exclusive\n";
        controller->Throw_SPONGE_Error(spongeErrorConflictingCommand,
                                       "SPONGE_PLUGIN::Initial",
                                       reason.c_str());
    }
    if (is_transactional &&
        (begin_transaction == NULL || commit_transaction == NULL ||
         rollback_transaction == NULL))
    {
        std::string reason =
            "Reason:\n\ttransactional Calculate_Force plugin '" + plugin_name +
            "' must export Begin_Force_Transaction, "
            "Commit_Force_Transaction, and Rollback_Force_Transaction\n";
        controller->Throw_SPONGE_Error(spongeErrorNotImplemented,
                                       "SPONGE_PLUGIN::Initial",
                                       reason.c_str());
    }

    if (!monte_carlo_barostat) return;

    if (!has_capability_getter)
    {
        std::string reason =
            "Reason:\n\tCalculate_Force plugin '" + plugin_name +
            "' cannot be used by the Monte Carlo barostat because it does "
            "not explicitly export Get_Force_Capabilities\n";
        controller->Throw_SPONGE_Error(spongeErrorNotImplemented,
                                       "SPONGE_PLUGIN::Initial",
                                       reason.c_str());
    }
    if ((capabilities & SPONGE_PLUGIN_FORCE_ENERGY_COMPLETE) == 0)
    {
        std::string reason =
            "Reason:\n\tCalculate_Force plugin '" + plugin_name +
            "' cannot be used by the Monte Carlo barostat without the "
            "ENERGY_COMPLETE capability\n";
        controller->Throw_SPONGE_Error(spongeErrorNotImplemented,
                                       "SPONGE_PLUGIN::Initial",
                                       reason.c_str());
    }
    if ((capabilities & SPONGE_PLUGIN_FORCE_VIRIAL_COMPLETE) == 0)
    {
        std::string reason =
            "Reason:\n\tCalculate_Force plugin '" + plugin_name +
            "' cannot be used by the Monte Carlo barostat without the "
            "VIRIAL_COMPLETE capability\n";
        controller->Throw_SPONGE_Error(spongeErrorNotImplemented,
                                       "SPONGE_PLUGIN::Initial",
                                       reason.c_str());
    }
    if (!is_pure && !is_transactional)
    {
        std::string reason =
            "Reason:\n\tCalculate_Force plugin '" + plugin_name +
            "' cannot be used by the Monte Carlo barostat without exactly "
            "one of PURE or TRANSACTIONAL\n";
        controller->Throw_SPONGE_Error(spongeErrorNotImplemented,
                                       "SPONGE_PLUGIN::Initial",
                                       reason.c_str());
    }
}
}  // namespace

void SPONGE_PLUGIN::Initial(MD_INFORMATION* md_info, CONTROLLER* controller,
                            COLLECTIVE_VARIABLE_CONTROLLER* cv_controller,
                            NEIGHBOR_LIST* neighbor_list)
{
    if (!controller->Command_Exist("plugin"))
    {
        return;
    }

    controller->printf("START INITIALIZING SPONGE PLUGIN:\n");
    plugin_numbers = 0;
    g_plugin_md_info = md_info;
    g_plugin_controller = controller;
    g_plugin_neighbor_list = neighbor_list;
    g_plugin_domain_info = NULL;

    const std::string command(controller->Original_Command("plugin"));
    std::istringstream plugin_path_stream(command);
    std::vector<std::string> plugin_paths;
    std::string parsed_plugin_path;
    while (plugin_path_stream >> parsed_plugin_path)
    {
        plugin_paths.push_back(parsed_plugin_path);
    }
    if (plugin_paths.empty() ||
        plugin_paths.size() > static_cast<std::size_t>(INT_MAX))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "SPONGE_PLUGIN::Initial",
            "Reason:\n\tplugin must contain a nonempty list whose count "
            "fits in int\n");
    }
    plugin_numbers = static_cast<int>(plugin_paths.size());

    controller->printf("%d plugin(s) to load\n", plugin_numbers);
    Malloc_Safely((void**)&plugin_handles, sizeof(HMODULE) * plugin_numbers);
    Malloc_Safely((void**)&after_init_funcs,
                  sizeof(RuntimeFunction) * plugin_numbers);
    Malloc_Safely((void**)&force_funcs,
                  sizeof(RuntimeFunction) * plugin_numbers);
    Malloc_Safely((void**)&force_capabilities,
                  sizeof(uint32_t) * plugin_numbers);
    Malloc_Safely((void**)&begin_force_transaction_funcs,
                  sizeof(RuntimeFunction) * plugin_numbers);
    Malloc_Safely((void**)&commit_force_transaction_funcs,
                  sizeof(RuntimeFunction) * plugin_numbers);
    Malloc_Safely((void**)&rollback_force_transaction_funcs,
                  sizeof(RuntimeFunction) * plugin_numbers);
    Malloc_Safely((void**)&print_funcs,
                  sizeof(RuntimeFunction) * plugin_numbers);
    Malloc_Safely((void**)&set_domain_info_funcs,
                  sizeof(SetDomainInformationFunction) * plugin_numbers);

    int count = 0;
    std::string plugin_name, plugin_version, version_check_error;
    NameFunction name_func, version_func;
    VersionCheckFunction version_check_func;
    SetBackendDeviceTypeFunction set_backend_device_type_func;
    InitialStableFunction stable_initial_func;

    for (const std::string& plugin_path : plugin_paths)
    {
        int funcs_loaded = 1;
        int negotiated_stable_api_version = 0;

#ifdef _WIN32
        constexpr int dlopen_mode = 0;
#else
        constexpr int dlopen_mode = RTLD_LAZY | RTLD_GLOBAL;
#endif
        plugin_handles[count] = dlopen(plugin_path.c_str(), dlopen_mode);
        if (plugin_handles[count] == NULL)
        {
            std::string error_reason = "Reason:\n\tOpen Dynamic Library from ";
            error_reason += plugin_path;
            error_reason += " failed\n";
            error_reason += DlErrorString();
            controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                           "SPONGE_PLUGIN::Initial",
                                           error_reason.c_str());
        }

        name_func = (NameFunction)dlsym(plugin_handles[count], "Name");
        if (name_func == NULL)
        {
            std::string error_reason =
                "Reason:\n\tFind the name of the plugin from ";
            error_reason += plugin_path;
            error_reason += " failed\n";
            error_reason += DlErrorString();
            controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                           "SPONGE_PLUGIN::Initial",
                                           error_reason.c_str());
        }

        plugin_name = name_func();
        version_func = (NameFunction)dlsym(plugin_handles[count], "Version");
        if (version_func == NULL)
        {
            std::string error_reason =
                "Reason:\n\tFind the version of the plugin from ";
            error_reason += plugin_path;
            error_reason += " (" + plugin_name + ") failed\n";
            error_reason += DlErrorString();
            controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                           "SPONGE_PLUGIN::Initial",
                                           error_reason.c_str());
        }

        plugin_version = version_func();
        version_check_func =
            (VersionCheckFunction)dlsym(plugin_handles[count], "Version_Check");
        if (version_check_func == NULL)
        {
            std::string error_reason =
                "Reason:\n\tFind the version check function of the plugin "
                "from ";
            error_reason += plugin_path;
            error_reason += " (" + plugin_name + " version: " + plugin_version +
                            ") failed\n";
            error_reason += DlErrorString();
            controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                           "SPONGE_PLUGIN::Initial",
                                           error_reason.c_str());
        }

        stable_initial_func = (InitialStableFunction)dlsym(
            plugin_handles[count], "Initial_Stable");
        if (stable_initial_func != NULL)
        {
            // Stable API revisions are append-only.  Negotiate the newest
            // revision accepted by the plugin so v2-v4 binaries continue to
            // consume the prefix they were compiled against, while v5
            // plugins can require synchronous functional writeback.
            constexpr int oldest_stable_api_version = 2;
            for (int candidate = SPONGE_PRIPS_API_VERSION;
                 candidate >= oldest_stable_api_version; --candidate)
            {
                version_check_error = version_check_func(candidate);
                if (version_check_error.empty())
                {
                    negotiated_stable_api_version = candidate;
                    break;
                }
            }
        }
        else
        {
            version_check_error =
                version_check_func(controller->last_modify_date);
        }
        if (!version_check_error.empty())
        {
            std::string error_reason =
                "Reason:\n\tThe version check of the plugin from ";
            error_reason += plugin_path;
            error_reason += " (" + plugin_name + " version: " + plugin_version +
                            ") failed\n" + version_check_error;
            error_reason += DlErrorString();
            controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                           "SPONGE_PLUGIN::Initial",
                                           error_reason.c_str());
        }

        controller->printf(
            "Plugin %d:\n    name: %s\n    version: %s\n    path: %s\n    "
            "functions loaded: ",
            plugin_numbers, plugin_name.c_str(), plugin_version.c_str(),
            plugin_path.c_str());

        InitialFunction func =
            (InitialFunction)dlsym(plugin_handles[count], "Initial");
        if (func == NULL && stable_initial_func == NULL)
        {
            std::string error_reason =
                "Reason:\n\tFind the initial function of the plugin from ";
            error_reason += plugin_path;
            error_reason += " (" + plugin_name + " version: " + plugin_version +
                            ") failed\n";
            error_reason += DlErrorString();
            controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                           "SPONGE_PLUGIN::Initial",
                                           error_reason.c_str());
        }

        controller->printf(" Initial");

        set_backend_device_type_func = (SetBackendDeviceTypeFunction)dlsym(
            plugin_handles[count], "Set_Backend_Device_Type");
        if (set_backend_device_type_func != NULL)
        {
            funcs_loaded += 1;
            controller->printf(" Set_Backend_Device_Type");
            set_backend_device_type_func(PluginBackendDeviceType());
        }

        after_init_funcs[after_init_func_numbers] =
            (RuntimeFunction)dlsym(plugin_handles[count], "After_Initial");
        if (after_init_funcs[after_init_func_numbers] != NULL)
        {
            funcs_loaded += 1;
            after_init_func_numbers += 1;
            controller->printf(" After_Initial");
        }

        const int current_force_index = force_func_numbers;
        RuntimeFunction current_force_func =
            (RuntimeFunction)dlsym(plugin_handles[count], "Calculate_Force");
        if (current_force_func != NULL)
        {
            force_funcs[current_force_index] = current_force_func;
            funcs_loaded += 1;
            controller->printf(" Calculate_Force");
        }

        print_funcs[print_func_numbers] =
            (RuntimeFunction)dlsym(plugin_handles[count], "Mdout_Print");
        if (print_funcs[print_func_numbers] != NULL)
        {
            funcs_loaded += 1;
            print_func_numbers += 1;
            controller->printf(" Mdout_Print");
        }

        set_domain_info_funcs[set_domain_info_func_numbers] =
            (SetDomainInformationFunction)dlsym(plugin_handles[count],
                                                "Set_Domain_Information");
        if (set_domain_info_funcs[set_domain_info_func_numbers] != NULL)
        {
            funcs_loaded += 1;
            set_domain_info_func_numbers += 1;
            controller->printf(" Set_Domain_Information");
        }

        controller->printf(" (%d in total)\n", funcs_loaded);
        if (stable_initial_func != NULL)
        {
            stable_initial_func(BuildPripsApi(
                static_cast<uint32_t>(negotiated_stable_api_version)));
        }
        else
        {
            func(md_info, controller, neighbor_list, cv_controller, CV_MAP,
                 CV_INSTANCE_MAP);
        }

        // Stable plugins such as PRIPS finish loading their user module in
        // Initial_Stable.  Capabilities must therefore be queried only after
        // initialization, and are still required independently of which
        // loader entry point the plugin uses.
        if (current_force_func != NULL)
        {
            GetForceCapabilitiesFunction get_force_capabilities =
                (GetForceCapabilitiesFunction)dlsym(plugin_handles[count],
                                                    "Get_Force_Capabilities");
            begin_force_transaction_funcs[current_force_index] =
                (RuntimeFunction)dlsym(plugin_handles[count],
                                       "Begin_Force_Transaction");
            commit_force_transaction_funcs[current_force_index] =
                (RuntimeFunction)dlsym(plugin_handles[count],
                                       "Commit_Force_Transaction");
            rollback_force_transaction_funcs[current_force_index] =
                (RuntimeFunction)dlsym(plugin_handles[count],
                                       "Rollback_Force_Transaction");
            const uint32_t capabilities =
                get_force_capabilities == NULL ? 0 : get_force_capabilities();
            ValidateForcePluginContract(
                controller, plugin_name, get_force_capabilities != NULL,
                capabilities,
                begin_force_transaction_funcs[current_force_index],
                commit_force_transaction_funcs[current_force_index],
                rollback_force_transaction_funcs[current_force_index],
                UsesMonteCarloBarostat(md_info, controller));
            force_capabilities[current_force_index] = capabilities;
            force_func_numbers += 1;
        }

        count += 1;
    }

    controller->printf("END INITIALIZING SPONGE PLUGIN\n\n");
}

void SPONGE_PLUGIN::Set_Domain_Information(DOMAIN_INFORMATION* dd)
{
    g_plugin_domain_info = dd;
    for (int i = 0; i < set_domain_info_func_numbers; i++)
    {
        set_domain_info_funcs[i](dd);
    }
}

void SPONGE_PLUGIN::After_Initial()
{
    for (int i = 0; i < after_init_func_numbers; i++)
    {
        after_init_funcs[i]();
    }
}

void SPONGE_PLUGIN::Calculate_Force(bool commit_sampling_state,
                                    bool exact_state, bool needs_energy,
                                    bool needs_virial)
{
    g_force_evaluation_commits_sampling_state = commit_sampling_state;
    g_force_evaluation_is_exact = exact_state;
    g_force_evaluation_needs_energy = needs_energy;
    g_force_evaluation_needs_virial = needs_virial;
    // Every callback gets an isolated transaction.  Old/trial MC probes are
    // always rolled back; the one final committed evaluation advances both
    // transient plugin state and sampling/history state exactly once.
    Begin_Force_Transaction();
    for (int i = 0; i < force_func_numbers; i++)
    {
        force_funcs[i]();
    }
    if (commit_sampling_state)
        Commit_Force_Transaction();
    else
        Rollback_Force_Transaction();
}

void SPONGE_PLUGIN::Begin_Force_Transaction()
{
    if (force_transaction_active)
    {
        g_plugin_controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "SPONGE_PLUGIN::Begin_Force_Transaction",
            "Reason:\n\ta plugin force transaction is already active\n");
    }
    force_transaction_active = true;
    for (int i = 0; i < force_func_numbers; i++)
    {
        if ((force_capabilities[i] & SPONGE_PLUGIN_FORCE_TRANSACTIONAL) != 0)
        {
            begin_force_transaction_funcs[i]();
        }
    }
}

void SPONGE_PLUGIN::Commit_Force_Transaction()
{
    if (!force_transaction_active)
    {
        g_plugin_controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "SPONGE_PLUGIN::Commit_Force_Transaction",
            "Reason:\n\tno plugin force transaction is active\n");
    }
    for (int i = 0; i < force_func_numbers; i++)
    {
        if ((force_capabilities[i] & SPONGE_PLUGIN_FORCE_TRANSACTIONAL) != 0)
        {
            commit_force_transaction_funcs[i]();
        }
    }
    force_transaction_active = false;
}

void SPONGE_PLUGIN::Rollback_Force_Transaction()
{
    if (!force_transaction_active)
    {
        g_plugin_controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "SPONGE_PLUGIN::Rollback_Force_Transaction",
            "Reason:\n\tno plugin force transaction is active\n");
    }
    for (int i = force_func_numbers - 1; i >= 0; i--)
    {
        if ((force_capabilities[i] & SPONGE_PLUGIN_FORCE_TRANSACTIONAL) != 0)
        {
            rollback_force_transaction_funcs[i]();
        }
    }
    force_transaction_active = false;
}

void SPONGE_PLUGIN::Mdout_Print()
{
    for (int i = 0; i < print_func_numbers; i++)
    {
        print_funcs[i]();
    }
}
