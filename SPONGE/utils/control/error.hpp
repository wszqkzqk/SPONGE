#pragma once

// SPONGE错误类型
enum spongeError
{
    spongeSuccess = 0,
    // 1000以下的错误留给 deviceError
    // 未实现的功能
    spongeErrorNotImplemented = 1001,
    // 文件格式（编码、换行符）问题 或 数据格式不正确
    spongeErrorBadFileFormat,
    // 冲突的命令
    spongeErrorConflictingCommand,
    // 缺失的命令
    spongeErrorMissingCommand,
    // 类型错误的命令
    spongeErrorTypeErrorCommand,
    // 值错误的命令
    spongeErrorValueErrorCommand,
    // 模拟崩溃
    spongeErrorSimulationBreakDown,
    // 内存分配失败
    spongeErrorMallocFailed,
    // 越界
    spongeErrorOverflow,
    // 打开文件失败
    spongeErrorOpenFileFailed,
};

inline bool CONTROLLER::Register_Fatal_Cleanup(
    const Fatal_Cleanup_Callback callback, void* context) noexcept
{
    if (callback == NULL) return false;

    Lock_Fatal_Cleanup_Registry();
    if (fatal_sequence_claimed_.load(std::memory_order_acquire))
    {
        Unlock_Fatal_Cleanup_Registry();
        return false;
    }
    for (std::size_t i = 0; i < fatal_cleanup_count_; ++i)
    {
        if (fatal_cleanup_entries_[i].callback == callback &&
            fatal_cleanup_entries_[i].context == context)
        {
            Unlock_Fatal_Cleanup_Registry();
            return false;
        }
    }
    if (fatal_cleanup_count_ >= CONTROLLER::FATAL_CLEANUP_CAPACITY)
    {
        Unlock_Fatal_Cleanup_Registry();
        return false;
    }
    fatal_cleanup_entries_[fatal_cleanup_count_].callback = callback;
    fatal_cleanup_entries_[fatal_cleanup_count_].context = context;
    fatal_cleanup_entries_[fatal_cleanup_count_].state =
        FATAL_CLEANUP_REGISTERED;
    ++fatal_cleanup_count_;
    Unlock_Fatal_Cleanup_Registry();
    return true;
}

inline bool CONTROLLER::Unregister_Fatal_Cleanup(
    const Fatal_Cleanup_Callback callback, void* context) noexcept
{
    if (callback == NULL) return false;

    Lock_Fatal_Cleanup_Registry();
    for (std::size_t upper = fatal_cleanup_count_; upper != 0; --upper)
    {
        const std::size_t index = upper - 1;
        if (fatal_cleanup_entries_[index].callback != callback ||
            fatal_cleanup_entries_[index].context != context)
        {
            continue;
        }

        Fatal_Cleanup_Entry& entry = fatal_cleanup_entries_[index];
        if (entry.state == FATAL_CLEANUP_DONE)
        {
            Unlock_Fatal_Cleanup_Registry();
            return true;
        }
        if (entry.state == FATAL_CLEANUP_RUNNING)
        {
            // No callback can wait for itself (or an outer callback in a
            // same-thread nested fatal path) to return. False preserves the
            // lifetime contract: the context is not released to this caller.
            if (fatal_sequence_owner_)
            {
                Unlock_Fatal_Cleanup_Registry();
                return false;
            }

            Unlock_Fatal_Cleanup_Registry();
            for (;;)
            {
                Lock_Fatal_Cleanup_Registry();
                const bool callback_returned =
                    fatal_cleanup_entries_[index].state == FATAL_CLEANUP_DONE;
                Unlock_Fatal_Cleanup_Registry();
                if (callback_returned) return true;
                std::this_thread::yield();
            }
        }
        if (entry.state == FATAL_CLEANUP_REGISTERED)
        {
            if (fatal_sequence_claimed_.load(std::memory_order_acquire))
            {
                // The fatal owner may already be scanning the stack. Marking
                // the entry done under the same lock makes removal atomic with
                // claiming it for execution.
                entry.state = FATAL_CLEANUP_DONE;
            }
            else
            {
                // Before fatal cleanup starts, compacting preserves exact
                // registration order when the freed capacity is reused.
                for (std::size_t move = index + 1; move < fatal_cleanup_count_;
                     ++move)
                {
                    fatal_cleanup_entries_[move - 1] =
                        fatal_cleanup_entries_[move];
                }
                --fatal_cleanup_count_;
                fatal_cleanup_entries_[fatal_cleanup_count_] = {
                    NULL, NULL, FATAL_CLEANUP_EMPTY};
            }
            Unlock_Fatal_Cleanup_Registry();
            return true;
        }
    }
    Unlock_Fatal_Cleanup_Registry();
    return false;
}

inline void CONTROLLER::Lock_Fatal_Cleanup_Registry() noexcept
{
    while (fatal_cleanup_lock_.test_and_set(std::memory_order_acquire))
    {
    }
}

inline void CONTROLLER::Unlock_Fatal_Cleanup_Registry() noexcept
{
    fatal_cleanup_lock_.clear(std::memory_order_release);
}

inline void CONTROLLER::Begin_Fatal_Sequence() noexcept
{
    if (fatal_sequence_owner_) return;

    bool expected = false;
    if (fatal_sequence_claimed_.compare_exchange_strong(
            expected, true, std::memory_order_acq_rel,
            std::memory_order_acquire))
    {
        fatal_sequence_owner_ = true;
        return;
    }

    // Only the first fatal caller may clean up, diagnose, or terminate. A
    // concurrent caller cannot safely return to code that has already entered
    // a process-fatal state, so it waits for the owner to terminate the
    // process.
    Wait_For_Fatal_Termination();
}

[[noreturn]] inline void CONTROLLER::Wait_For_Fatal_Termination() noexcept
{
    fatal_sequence_waiters_.fetch_add(1, std::memory_order_release);
    for (;;)
    {
        std::this_thread::yield();
    }
}

inline void CONTROLLER::Run_Fatal_Cleanups() noexcept
{
    for (;;)
    {
        Lock_Fatal_Cleanup_Registry();
        std::size_t index = FATAL_CLEANUP_CAPACITY;
        for (std::size_t upper = fatal_cleanup_count_; upper != 0; --upper)
        {
            const std::size_t candidate = upper - 1;
            if (fatal_cleanup_entries_[candidate].state ==
                FATAL_CLEANUP_REGISTERED)
            {
                index = candidate;
                break;
            }
        }
        if (index == FATAL_CLEANUP_CAPACITY)
        {
            Unlock_Fatal_Cleanup_Registry();
            return;
        }

        Fatal_Cleanup_Entry& entry = fatal_cleanup_entries_[index];
        const Fatal_Cleanup_Callback callback = entry.callback;
        void* const context = entry.context;
        entry.state = FATAL_CLEANUP_RUNNING;
        Unlock_Fatal_Cleanup_Registry();

        try
        {
            callback(context);
        }
        catch (...)
        {
        }

        Lock_Fatal_Cleanup_Registry();
        if (fatal_cleanup_entries_[index].state == FATAL_CLEANUP_RUNNING)
        {
            fatal_cleanup_entries_[index].state = FATAL_CLEANUP_DONE;
        }
        Unlock_Fatal_Cleanup_Registry();
    }
}

inline void CONTROLLER::Resolve_Fatal_Base(const int error_number,
                                           const char** error_name,
                                           const char** error_reason) noexcept
{
    switch (error_number)
    {
        case spongeErrorNotImplemented:
            *error_name = "spongeErrorNotImplemented";
            *error_reason =
                "The function has not been implemented in SPONGE yet";
            return;
        case spongeErrorBadFileFormat:
            *error_name = "spongeErrorBadFileFormat";
            *error_reason = "The format of the file is bad";
            return;
        case spongeErrorConflictingCommand:
            *error_name = "spongeErrorConflictingCommand";
            *error_reason = "Some commands are conflicting";
            return;
        case spongeErrorMissingCommand:
            *error_name = "spongeErrorMissingCommand";
            *error_reason = "Missing required command(s)";
            return;
        case spongeErrorTypeErrorCommand:
            *error_name = "spongeErrorTypeErrorCommand";
            *error_reason = "The type of the command is wrong";
            return;
        case spongeErrorValueErrorCommand:
            *error_name = "spongeErrorValueErrorCommand";
            *error_reason = "The value of the command is wrong";
            return;
        case spongeErrorSimulationBreakDown:
            *error_name = "spongeErrorSimulationBrokenDown";
            *error_reason = "The system was broken down";
            return;
        case spongeErrorMallocFailed:
            *error_name = "spongeErrorMallocFailed";
            *error_reason = "Fail to allocate memory";
            return;
        case spongeErrorOverflow:
            *error_name = "spongeErrorOverflow";
            *error_reason = "Boundary was overflowed";
            return;
        case spongeErrorOpenFileFailed:
            *error_name = "spongeErrorOpenFileFailed";
            *error_reason = "Fail to open file";
            return;
        default:
            *error_name = "spongeErrorUnclassified";
            *error_reason = "Unclassified Error";
            return;
    }
}

inline bool CONTROLLER::Capture_Fatal_Diagnostic_From_Base(
    const char* error_name, const char* error_reason, const char* error_by,
    const char* extra_error_string, std::vector<char>* diagnostic) noexcept
{
    const std::size_t name_length = strlen(error_name);
    const std::size_t reason_length = strlen(error_reason);
    const std::size_t by_length = error_by == NULL ? 0 : strlen(error_by);
    const std::size_t extra_length =
        extra_error_string == NULL ? 0 : strlen(extra_error_string);
    const bool extra_needs_newline =
        extra_length != 0 && extra_error_string[extra_length - 1] != '\n';

    std::size_t required = 0;
    const auto add_size = [&required](const std::size_t addition) noexcept
    {
        if (addition > std::numeric_limits<std::size_t>::max() - required)
            return false;
        required += addition;
        return true;
    };
    if (!add_size(1) || !add_size(name_length) ||
        (error_by != NULL &&
         (!add_size(sizeof(" raised by ") - 1) || !add_size(by_length))) ||
        !add_size(1) || !add_size(reason_length) || !add_size(1) ||
        !add_size(extra_length) || (extra_needs_newline && !add_size(1)))
    {
        return false;
    }

    try
    {
        diagnostic->resize(required);
    }
    catch (...)
    {
        diagnostic->clear();
        return false;
    }

    char* cursor = diagnostic->data();
    const auto append =
        [&cursor](const char* source, const std::size_t length) noexcept
    {
        if (length != 0)
        {
            memcpy(cursor, source, length);
            cursor += length;
        }
    };
    append("\n", 1);
    append(error_name, name_length);
    if (error_by != NULL)
    {
        append(" raised by ", sizeof(" raised by ") - 1);
        append(error_by, by_length);
    }
    append("\n", 1);
    append(error_reason, reason_length);
    append("\n", 1);
    append(extra_error_string, extra_length);
    if (extra_needs_newline) append("\n", 1);
    return true;
}

inline std::size_t CONTROLLER::Capture_Fatal_Diagnostic_Fixed_From_Base(
    const char* error_name, const char* error_reason, const char* error_by,
    const char* extra_error_string, char* diagnostic,
    const std::size_t capacity) noexcept
{
    if (diagnostic == NULL || capacity == 0) return 0;

    static constexpr char TRUNCATION_MARKER[] =
        "\n[diagnostic truncated: dynamic allocation unavailable]\n";
    constexpr std::size_t marker_length = sizeof(TRUNCATION_MARKER) - 1;

    const std::size_t name_length = strlen(error_name);
    const std::size_t reason_length = strlen(error_reason);
    const std::size_t by_length = error_by == NULL ? 0 : strlen(error_by);
    const std::size_t extra_length =
        extra_error_string == NULL ? 0 : strlen(extra_error_string);
    const bool extra_needs_newline =
        extra_length != 0 && extra_error_string[extra_length - 1] != '\n';

    std::size_t required = 0;
    const auto add_size = [&required](const std::size_t addition) noexcept
    {
        if (addition > std::numeric_limits<std::size_t>::max() - required)
            return false;
        required += addition;
        return true;
    };
    const bool length_representable =
        add_size(1) && add_size(name_length) &&
        (error_by == NULL ||
         (add_size(sizeof(" raised by ") - 1) && add_size(by_length))) &&
        add_size(1) && add_size(reason_length) && add_size(1) &&
        add_size(extra_length) && (!extra_needs_newline || add_size(1));
    const bool truncated = !length_representable || required > capacity;
    const std::size_t copied_marker_length =
        truncated ? std::min(capacity, marker_length) : 0;
    const std::size_t content_capacity = capacity - copied_marker_length;

    std::size_t cursor = 0;
    const auto append =
        [&](const char* source, const std::size_t length) noexcept
    {
        const std::size_t available = content_capacity - cursor;
        const std::size_t copied = std::min(available, length);
        if (copied != 0)
        {
            memcpy(diagnostic + cursor, source, copied);
            cursor += copied;
        }
    };
    append("\n", 1);
    append(error_name, name_length);
    if (error_by != NULL)
    {
        append(" raised by ", sizeof(" raised by ") - 1);
        append(error_by, by_length);
    }
    append("\n", 1);
    append(error_reason, reason_length);
    append("\n", 1);
    append(extra_error_string, extra_length);
    if (extra_needs_newline) append("\n", 1);

    if (truncated)
    {
        memcpy(diagnostic + cursor, TRUNCATION_MARKER, copied_marker_length);
        cursor += copied_marker_length;
    }
    return cursor;
}

inline bool CONTROLLER::Capture_Fatal_Diagnostic(
    const int error_number, const char* error_by,
    const char* extra_error_string, std::vector<char>* diagnostic) noexcept
{
    const char* error_name = NULL;
    const char* error_reason = NULL;
    Resolve_Fatal_Base(error_number, &error_name, &error_reason);
    return Capture_Fatal_Diagnostic_From_Base(
        error_name, error_reason, error_by, extra_error_string, diagnostic);
}

inline void CONTROLLER::Write_Fatal_Bytes(const char* text,
                                          std::size_t remaining) noexcept
{
    if (text == NULL) return;

    while (remaining != 0)
    {
#ifdef _WIN32
        const unsigned int chunk = remaining > static_cast<std::size_t>(INT_MAX)
                                       ? static_cast<unsigned int>(INT_MAX)
                                       : static_cast<unsigned int>(remaining);
        const int written = _write(2, text, chunk);
#else
        const std::size_t chunk =
            remaining > static_cast<std::size_t>(SSIZE_MAX)
                ? static_cast<std::size_t>(SSIZE_MAX)
                : remaining;
        const ssize_t written = write(2, text, chunk);
#endif
        if (written > 0)
        {
            text += static_cast<std::size_t>(written);
            remaining -= static_cast<std::size_t>(written);
            continue;
        }
        if (written < 0 && errno == EINTR) continue;
        return;
    }
}

inline void CONTROLLER::Write_Fatal_Text(const char* text) noexcept
{
    if (text != NULL) Write_Fatal_Bytes(text, strlen(text));
}

[[noreturn]] inline void CONTROLLER::Finish_Fatal_Sequence(
    const int error_number, const char* captured_diagnostic,
    const std::size_t captured_length) noexcept
{
    Run_Fatal_Cleanups();
    if (captured_diagnostic != NULL)
    {
        Write_Fatal_Bytes(captured_diagnostic, captured_length);
    }
    else
    {
        // Allocation failure and spongeErrorMallocFailed must not consult
        // caller-owned pointers after cleanup.  A static base diagnostic is
        // always available, even if device-specific lookup is no longer safe.
        const char* error_name = NULL;
        const char* error_reason = NULL;
        Resolve_Fatal_Base(error_number, &error_name, &error_reason);
        Write_Fatal_Text("\n");
        Write_Fatal_Text(error_name);
        Write_Fatal_Text("\n");
        Write_Fatal_Text(error_reason);
        Write_Fatal_Text("\n");
    }
    Terminate_Fatal_Sequence(error_number);
}

[[noreturn]] inline void CONTROLLER::Finish_Uncaptured_Fatal_Sequence(
    const int error_number, const char* error_name, const char* error_reason,
    const char* error_by, const char* extra_error_string) noexcept
{
    // Heap exhaustion must not force diagnostics to consume caller-owned
    // storage or touch a potentially blocking stderr before cleanup. Snapshot
    // a bounded diagnostic on this stack, then use the same cleanup-first
    // finish sequence as the dynamically captured path.
    char diagnostic[FATAL_DIAGNOSTIC_FALLBACK_CAPACITY];
    const std::size_t diagnostic_length =
        Capture_Fatal_Diagnostic_Fixed_From_Base(
            error_name, error_reason, error_by, extra_error_string, diagnostic,
            sizeof(diagnostic));
    Finish_Fatal_Sequence(error_number, diagnostic, diagnostic_length);
}

[[noreturn]] inline void CONTROLLER::Terminate_Fatal_Sequence(
    const int error_number) noexcept
{
#ifdef USE_MPI
    // MPI_Initialized is explicitly valid before MPI_Init. Never invoke
    // MPI_Abort after a failed MPI_Init_thread call left MPI unavailable.
    int initialized = 0;
    int finalized = 0;
    if (MPI_Initialized(&initialized) == MPI_SUCCESS && initialized != 0 &&
        MPI_Finalized(&finalized) == MPI_SUCCESS && finalized == 0)
    {
        // Some MPI shims and faulty implementations return from MPI_Abort.
        // The local hard-termination fallback is therefore unconditional.
        (void)MPI_Abort(MPI_COMM_WORLD, error_number);
    }
#endif
    std::_Exit(error_number);
}

inline void CONTROLLER::Throw_SPONGE_Error(const int error_number,
                                           const char* error_by,
                                           const char* extra_error_string)
{
    if (error_number == 0) return;
    Begin_Fatal_Sequence();
    // A cleanup callback is allowed to destroy the state that supplied
    // error_by or extra_error_string.  Snapshot every byte before running any
    // callback, then never consult caller-owned diagnostic storage again. Try
    // to preserve the complete diagnostic even for allocation errors; a fixed
    // stack snapshot is the allocation-free fallback.

    std::vector<char> diagnostic;
    if (!Capture_Fatal_Diagnostic(error_number, error_by, extra_error_string,
                                  &diagnostic))
    {
        const char* error_name = NULL;
        const char* error_reason = NULL;
        Resolve_Fatal_Base(error_number, &error_name, &error_reason);
        Finish_Uncaptured_Fatal_Sequence(error_number, error_name, error_reason,
                                         error_by, extra_error_string);
    }
    Finish_Fatal_Sequence(error_number, diagnostic.data(), diagnostic.size());
}

#ifdef GPU_ARCH_NAME
inline void CONTROLLER::Throw_Device_Error(const deviceError_t error_number,
                                           const char* error_by,
                                           const char* extra_error_string)
{
    const int numeric_error = static_cast<int>(error_number);
    if (numeric_error == 0) return;
    Begin_Fatal_Sequence();
    const char* error_name = deviceGetErrorName(error_number);
    const char* error_reason = deviceGetErrorString(error_number);
    if (error_name == NULL) error_name = "deviceErrorUnclassified";
    if (error_reason == NULL) error_reason = "Unclassified device error";

    std::vector<char> diagnostic;
    if (!Capture_Fatal_Diagnostic_From_Base(error_name, error_reason, error_by,
                                            extra_error_string, &diagnostic))
    {
        // Device-provided strings may cease to be safe after device cleanup.
        // Consume them before callbacks if the owned snapshot cannot be made.
        Finish_Uncaptured_Fatal_Sequence(numeric_error, error_name,
                                         error_reason, error_by,
                                         extra_error_string);
    }
    Finish_Fatal_Sequence(numeric_error, diagnostic.data(), diagnostic.size());
}
#endif

inline void CONTROLLER::Throw_Formatted_SPONGE_Error(const int error_number,
                                                     const char* error_by,
                                                     const char* format, ...)
{
    if (error_number == 0) return;
    Begin_Fatal_Sequence();
    // Formatting is best-effort even for spongeErrorMallocFailed.  If its
    // vector allocation fails, the catch path below delegates to the direct,
    // allocation-free diagnostic path with a static explanation.
    if (format == NULL)
    {
        Throw_SPONGE_Error(
            error_number, error_by,
            "Reason:\n\tSPONGE received a null error-message format\n");
        return;
    }

    va_list args;
    va_start(args, format);
    va_list length_args;
    va_copy(length_args, args);
    const int required_characters = vsnprintf(NULL, 0, format, length_args);
    va_end(length_args);
    if (required_characters < 0)
    {
        va_end(args);
        Throw_SPONGE_Error(
            error_number, error_by,
            "Reason:\n\tfailed to format a SPONGE error diagnostic\n");
        return;
    }

    const std::size_t required_size =
        static_cast<std::size_t>(required_characters);
    if (required_size == std::numeric_limits<std::size_t>::max())
    {
        va_end(args);
        Throw_SPONGE_Error(
            error_number, error_by,
            "Reason:\n\ta SPONGE error diagnostic exceeds the supported "
            "length\n");
        return;
    }
    std::vector<char> error_reason;
    try
    {
        error_reason.resize(required_size + 1);
    }
    catch (const std::bad_alloc&)
    {
        va_end(args);
        Throw_SPONGE_Error(
            error_number, error_by,
            "Reason:\n\tfailed to allocate the complete SPONGE error "
            "diagnostic\n");
        return;
    }
    catch (const std::length_error&)
    {
        va_end(args);
        Throw_SPONGE_Error(
            error_number, error_by,
            "Reason:\n\tthe SPONGE error diagnostic exceeds the supported "
            "length\n");
        return;
    }
    const int written =
        vsnprintf(error_reason.data(), error_reason.size(), format, args);
    va_end(args);
    if (written < 0 || written != required_characters)
    {
        Throw_SPONGE_Error(
            error_number, error_by,
            "Reason:\n\tfailed to format a SPONGE error diagnostic\n");
        return;
    }
    Throw_SPONGE_Error(error_number, error_by, error_reason.data());
}

inline void CONTROLLER::Check_Error(float energy)
{
#ifdef GPU_ARCH_NAME
    deviceError_t device_error = deviceGetLastError();
    if (device_error == deviceErrorInvalidConfiguration ||
        device_error == deviceErrorInvalidValue ||
        device_error == deviceErrorLaunchOutOfResources)
    {
        Throw_Device_Error(device_error, "CONTROLLER::Check_Error",
                           "Reasons:\n\tA device kernel function is launched "
                           "with wrong parameters, and this should be a bug. "
                           "Please report the issue to the developers.");
    }
    else if (device_error != 0)
    {
        Throw_Device_Error(
            device_error, "CONTROLLER::Check_Error",
            "Possible reasons:\n\t1. the energy of the system is not fully "
            "minimized\n\t2. bad dt (too large)\n\t3. bad thermostat/barostat "
            "parameters\n\t4. bad force field parameters\n");
    }
#endif
    if (isnan(energy) || isinf(energy) || isnan(printf_sum) ||
        isinf(printf_sum))
    {
        Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "CONTROLLER::Check_Error",
            "Possible reasons:\n\t1. the energy of the system is not fully "
            "minimized\n\t2. bad dt (too large)\n\t3. bad thermostat/barostat "
            "parameters\n\t4. bad force field parameters\n");
    }
    printf_sum = 0;
}

inline void CONTROLLER::Check_Int(const char* command, const char* error_by)
{
    const char* str = this->Command(command);
    if (!is_str_int(str))
    {
        char error_reason[CHAR_LENGTH_MAX];
        sprintf(
            error_reason,
            "Reason:\n\t the value '%s' of the command '%s' is not an int\n",
            str, command);
        this->Throw_SPONGE_Error(spongeErrorTypeErrorCommand, error_by,
                                 error_reason);
    }
}

inline void CONTROLLER::Check_Int(const char* prefix, const char* command,
                                  const char* error_by)
{
    const char* str = this->Command(prefix, command);
    if (!is_str_int(str))
    {
        char error_reason[CHAR_LENGTH_MAX];
        sprintf(
            error_reason,
            "Reason:\n\t the value '%s' of the command '%s' is not an int\n",
            str, command);
        this->Throw_SPONGE_Error(spongeErrorTypeErrorCommand, error_by,
                                 error_reason);
    }
}

inline void CONTROLLER::Check_Float(const char* command, const char* error_by)
{
    const char* str = this->Command(command);
    if (!is_str_float(str))
    {
        char error_reason[CHAR_LENGTH_MAX];
        sprintf(
            error_reason,
            "Reason:\n\t the value '%s' of the command '%s' is not a float\n",
            str, command);
        this->Throw_SPONGE_Error(spongeErrorTypeErrorCommand, error_by,
                                 error_reason);
    }
}

inline void CONTROLLER::Check_Float(const char* prefix, const char* command,
                                    const char* error_by)
{
    const char* str = this->Command(prefix, command);
    if (!is_str_float(str))
    {
        char error_reason[CHAR_LENGTH_MAX];
        sprintf(
            error_reason,
            "Reason:\n\t the value '%s' of the command '%s' is not a float\n",
            str, command);
        this->Throw_SPONGE_Error(spongeErrorTypeErrorCommand, error_by,
                                 error_reason);
    }
}

inline bool CONTROLLER::Get_Bool(const char* command, const char* error_by)
{
    const char* str = this->Command(command);
    if (is_str_equal(str, "true"))
    {
        return true;
    }
    else if (is_str_equal(str, "false"))
    {
        return false;
    }
    else
    {
        Check_Int(command, error_by);
        return atoi(str);
    }
}

inline bool CONTROLLER::Get_Bool(const char* prefix, const char* command,
                                 const char* error_by)
{
    const char* str = this->Command(prefix, command);
    if (is_str_equal(str, "true"))
    {
        return true;
    }
    else if (is_str_equal(str, "false"))
    {
        return false;
    }
    else
    {
        Check_Int(prefix, command, error_by);
        return atoi(str);
    }
}
