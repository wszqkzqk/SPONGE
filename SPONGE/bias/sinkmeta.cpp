#include "sinkmeta.h"

#include <sys/stat.h>

#include <atomic>
#include <cerrno>
#include <cstdint>
#include <iomanip>
#include <iterator>

#include "meta_io_transaction.h"

#ifdef _WIN32
#include <share.h>
#else
#include <fcntl.h>
#include <sys/types.h>
#endif

static void Copy_Meta_String(CONTROLLER* controller, char* destination,
                             std::size_t capacity, const char* source,
                             const char* field_name)
{
    if (destination == NULL || capacity == 0 || source == NULL)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
            "Metadynamics %s contains a null string or destination.",
            field_name);
    }
    const std::size_t length = strlen(source);
    if (length >= capacity)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "METADYNAMICS::Initial",
            "Metadynamics %s is %zu bytes, but its storage permits at most "
            "%zu bytes.",
            field_name, length, capacity - 1);
    }
    memcpy(destination, source, length + 1);
}

static int Meta_Effective_IO_Error(int error_number)
{
    return error_number == 0 ? EIO : error_number;
}

static fs::path Meta_Native_Path(CONTROLLER* controller,
                                 const std::string& utf8_name,
                                 const char* error_by)
{
    try
    {
        return fs::u8path(utf8_name);
    }
    catch (const std::length_error&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, error_by,
            "Metadynamics path '%s' exceeds the host path limit.",
            utf8_name.c_str());
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorMallocFailed, error_by,
            "Unable to allocate native path storage for '%s'.",
            utf8_name.c_str());
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to interpret metadynamics path '%s': %s.",
            utf8_name.c_str(), error.what());
    }
    return fs::path();
}

static std::string Meta_Random_Temporary_Basename(CONTROLLER* controller,
                                                  const char* error_by,
                                                  const std::string& final_name)
{
    std::uint8_t random_bytes[16] = {};
#ifdef _WIN32
    using BCryptGenRandomFunction =
        long(WINAPI*)(void*, unsigned char*, unsigned long, unsigned long);
    HMODULE bcrypt = LoadLibraryW(L"bcrypt.dll");
    BCryptGenRandomFunction generate =
        bcrypt == NULL ? NULL
                       : reinterpret_cast<BCryptGenRandomFunction>(
                             GetProcAddress(bcrypt, "BCryptGenRandom"));
    const long random_status =
        generate == NULL
            ? -1073741823L
            : generate(NULL, random_bytes, sizeof(random_bytes), 0x00000002UL);
    if (bcrypt != NULL) FreeLibrary(bcrypt);
    if (random_status < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to obtain system randomness for a temporary output next "
            "to '%s' (Windows status 0x%08lx).",
            final_name.c_str(), static_cast<unsigned long>(random_status));
    }
#else
    errno = 0;
    const int random_descriptor = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    int random_error =
        random_descriptor < 0 ? Meta_Effective_IO_Error(errno) : 0;
    std::size_t completed = 0;
    while (random_error == 0 && completed < sizeof(random_bytes))
    {
        errno = 0;
        const ssize_t count = read(random_descriptor, random_bytes + completed,
                                   sizeof(random_bytes) - completed);
        if (count > 0)
        {
            completed += static_cast<std::size_t>(count);
        }
        else if (count < 0 && errno == EINTR)
        {
            continue;
        }
        else
        {
            random_error = count < 0 ? Meta_Effective_IO_Error(errno) : EIO;
        }
    }
    if (random_descriptor >= 0 && close(random_descriptor) != 0 &&
        random_error == 0)
    {
        random_error = Meta_Effective_IO_Error(errno);
    }
    if (random_error != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to obtain system randomness for a temporary output next "
            "to '%s': %s.",
            final_name.c_str(), strerror(random_error));
    }
#endif

    static const char hexadecimal[] = "0123456789abcdef";
    char suffix[2 * sizeof(random_bytes) + 1] = {};
    for (std::size_t index = 0; index < sizeof(random_bytes); ++index)
    {
        suffix[2 * index] = hexadecimal[random_bytes[index] >> 4];
        suffix[2 * index + 1] = hexadecimal[random_bytes[index] & 0x0f];
    }
    try
    {
        return std::string(".sponge-tmp.") + suffix;
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorMallocFailed, error_by,
            "Unable to allocate a random temporary-output name next to '%s'.",
            final_name.c_str());
    }
    return std::string();
}

class Meta_Atomic_Output
{
   private:
    enum FatalCleanupState
    {
        fatal_cleanup_active = 0,
        fatal_cleanup_claimed = 1,
        fatal_cleanup_normal_busy = 2,
        fatal_cleanup_disarmed = 3,
    };

    enum NormalCleanupState
    {
        normal_cleanup_idle = 0,
        normal_cleanup_busy = 1,
        normal_cleanup_done = 2,
    };

    struct Fatal_Cleanup_Record
    {
        std::atomic<int> state{fatal_cleanup_disarmed};
        std::atomic<int> cleanup_error{0};
#ifndef _WIN32
        std::string temporary_basename;
        int directory_descriptor = -1;
        int temporary_descriptor = -1;
        dev_t temporary_device = 0;
        ino_t temporary_inode = 0;
#endif
    };

   public:
    Meta_Atomic_Output(CONTROLLER* controller, const std::string& final_name,
                       const char* error_by)
        : controller_(controller), final_name_(final_name), error_by_(error_by)
    {
        Open();
    }

    Meta_Atomic_Output(const Meta_Atomic_Output&) = delete;
    Meta_Atomic_Output& operator=(const Meta_Atomic_Output&) = delete;

    ~Meta_Atomic_Output() noexcept { Close_And_Remove(true); }

    FILE* Stream() const noexcept { return stream_; }

#ifdef SPONGE_META_ATOMIC_OUTPUT_PROBE
    using Probe_Hook = void (*)(Meta_Atomic_Output*);
    static void Set_Probe_Post_Rename_Hook(Probe_Hook hook) noexcept
    {
        probe_post_rename_hook_ = hook;
    }
    static void Set_Probe_Normal_Cleanup_Hook(Probe_Hook hook) noexcept
    {
        probe_normal_cleanup_hook_ = hook;
    }
    static void Set_Probe_Competing_Cleanup_Hook(Probe_Hook hook) noexcept
    {
        probe_competing_cleanup_hook_ = hook;
    }
    void Probe_Invoke_Fatal_Cleanup() noexcept
    {
        Fatal_Cleanup(&fatal_cleanup_record_);
    }
    void Probe_Close() noexcept { Close_And_Remove(true); }
#endif

    void Commit()
    {
        const char* failed_operation = NULL;
        int io_error = 0;
        if (stream_ == NULL)
        {
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Cannot commit metadynamics output '%s' from a null "
                "temporary stream.",
                final_name_.c_str());
            return;
        }
        if (ferror(stream_) != 0)
        {
            failed_operation = "write";
            io_error = EIO;
        }
        if (failed_operation == NULL)
        {
            errno = 0;
            if (fflush(stream_) != 0)
            {
                failed_operation = "flush";
                io_error = Meta_Effective_IO_Error(errno);
            }
        }
        if (failed_operation == NULL && ferror(stream_) != 0)
        {
            failed_operation = "flush";
            io_error = EIO;
        }
#ifndef _WIN32
        if (failed_operation == NULL && preserve_replaced_mode_)
        {
            errno = 0;
            if (fchmod(fileno(stream_), replaced_mode_) != 0)
            {
                failed_operation = "set permissions on";
                io_error = Meta_Effective_IO_Error(errno);
            }
        }
#endif
        if (failed_operation == NULL)
        {
            errno = 0;
#ifdef _WIN32
            if (_commit(_fileno(stream_)) != 0)
#else
            if (fsync(fileno(stream_)) != 0)
#endif
            {
                failed_operation = "synchronize";
                io_error = Meta_Effective_IO_Error(errno);
            }
        }
#ifndef _WIN32
        int identity_error = 0;
        // POSIX provides no rename-by-open-fd primitive.  This fresh check
        // removes stale/cached-identity cleanup bugs, while the unpredictable
        // basename and parent-directory write permissions remain the security
        // boundary until renameat executes.
        if (!Verify_Temporary_Identity(&identity_error) &&
            failed_operation == NULL)
        {
            failed_operation = "verify the identity of";
            io_error = identity_error;
        }
#endif
#ifdef _WIN32
        errno = 0;
        const int close_status = fclose(stream_);
        const int close_error = errno;
        stream_ = NULL;
        if (failed_operation == NULL && close_status != 0)
        {
            failed_operation = "close";
            io_error = Meta_Effective_IO_Error(close_error);
        }
#endif
        if (failed_operation != NULL)
        {
            const int cleanup_error = Close_And_Remove(true);
            if (cleanup_error != 0)
            {
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorOpenFileFailed, error_by_,
                    "Failed to %s the complete temporary output for '%s': "
                    "%s; the previous output was preserved, but removing the "
                    "temporary file also failed with error %d.",
                    failed_operation, final_name_.c_str(), strerror(io_error),
                    cleanup_error);
            }
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Failed to %s the complete temporary output for '%s': %s; "
                "the previous output was preserved.",
                failed_operation, final_name_.c_str(), strerror(io_error));
            return;
        }

#ifdef _WIN32
        bool replace_failed =
            ReplaceFileW(final_path_.c_str(), temporary_path_.c_str(), NULL, 0,
                         NULL, NULL) == 0;
        unsigned long replace_error =
            replace_failed ? GetLastError() : ERROR_SUCCESS;
        if (replace_failed && (replace_error == ERROR_FILE_NOT_FOUND ||
                               replace_error == ERROR_PATH_NOT_FOUND))
        {
            // There is no replaced file whose metadata must be retained.  Do
            // not request replacement here: if another writer creates the name
            // concurrently, failing is safer than discarding its metadata.
            replace_failed =
                MoveFileExW(temporary_path_.c_str(), final_path_.c_str(),
                            MOVEFILE_WRITE_THROUGH) == 0;
            replace_error = replace_failed ? GetLastError() : ERROR_SUCCESS;
        }
#else
        errno = 0;
        const bool replace_failed =
            renameat(directory_descriptor_, temporary_basename_.c_str(),
                     directory_descriptor_, final_basename_.c_str()) != 0;
        const int replace_error = Meta_Effective_IO_Error(errno);
#endif
        if (replace_failed)
        {
            const int cleanup_error = Close_And_Remove(true);
#ifdef _WIN32
            if (cleanup_error != 0)
            {
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorOpenFileFailed, error_by_,
                    "Unable to atomically replace '%s' (Windows error %lu); "
                    "the previous output was preserved, but removing the "
                    "temporary file also failed: %s.",
                    final_name_.c_str(), replace_error,
                    strerror(cleanup_error));
            }
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Unable to atomically replace '%s' with the completed "
                "temporary output (Windows error %lu); the previous output "
                "was preserved.",
                final_name_.c_str(), replace_error);
#else
            if (cleanup_error != 0)
            {
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorOpenFileFailed, error_by_,
                    "Unable to atomically replace '%s': %s; the previous "
                    "output was preserved, but removing the temporary file "
                    "also failed with error %d.",
                    final_name_.c_str(), strerror(replace_error),
                    cleanup_error);
            }
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Unable to atomically replace '%s' with the completed "
                "temporary output: %s; the previous output was preserved.",
                final_name_.c_str(), strerror(replace_error));
#endif
            return;
        }
        temporary_exists_ = false;
#ifdef SPONGE_META_ATOMIC_OUTPUT_PROBE
        if (probe_post_rename_hook_ != NULL) probe_post_rename_hook_(this);
#endif

#ifndef _WIN32
        errno = 0;
        const int committed_close_status = fclose(stream_);
        const int committed_close_error = Meta_Effective_IO_Error(errno);
        stream_ = NULL;

        errno = 0;
        const int directory_sync_status = fsync(directory_descriptor_);
        const int directory_error = Meta_Effective_IO_Error(errno);
        Close_And_Remove(true);
        if (committed_close_status != 0)
        {
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Metadynamics output '%s' was atomically replaced after its "
                "contents were synchronized, but closing the committed file "
                "failed: %s.%s",
                final_name_.c_str(), strerror(committed_close_error),
                directory_sync_status == 0
                    ? " The parent directory was synchronized."
                    : " Synchronizing the parent directory also failed.");
            return;
        }
        if (directory_sync_status != 0)
        {
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Metadynamics output '%s' was atomically replaced, but "
                "synchronizing its parent directory failed: %s.  The new "
                "output is visible, but crash durability is not guaranteed.",
                final_name_.c_str(), strerror(directory_error));
            return;
        }
#else
        Close_And_Remove(true);
#endif
    }

   private:
    static void Fatal_Cleanup(void* context) noexcept
    {
        Fatal_Cleanup_Record* record =
            static_cast<Fatal_Cleanup_Record*>(context);
        int expected = fatal_cleanup_active;
        if (!record->state.compare_exchange_strong(
                expected, fatal_cleanup_claimed, std::memory_order_acq_rel,
                std::memory_order_acquire))
        {
            while (record->state.load(std::memory_order_acquire) ==
                   fatal_cleanup_normal_busy)
            {
            }
            return;
        }

        int cleanup_error = 0;
#ifdef _WIN32
        // A CRT FILE opened without an explicit FILE_SHARE_DELETE contract has
        // no race-free conditional delete-by-file-ID operation.  Do not close
        // the writer's FILE from this thread and do not delete a pathname that
        // may have been replaced.  Process termination closes the stream; a
        // random temporary name can remain after a Windows fatal exit.
        cleanup_error = EBUSY;
#else
        cleanup_error = Remove_Matching_Path(
            record->directory_descriptor, record->temporary_basename,
            record->temporary_descriptor, record->temporary_device,
            record->temporary_inode);
        if (cleanup_error == ENOENT) cleanup_error = 0;
        errno = 0;
        if (close(record->temporary_descriptor) != 0 && cleanup_error == 0)
        {
            cleanup_error = Meta_Effective_IO_Error(errno);
        }
        errno = 0;
        if (close(record->directory_descriptor) != 0 && cleanup_error == 0)
        {
            cleanup_error = Meta_Effective_IO_Error(errno);
        }
#endif
        record->cleanup_error.store(cleanup_error, std::memory_order_release);
    }

#ifndef _WIN32
    static int Duplicate_Cleanup_Descriptor(int descriptor,
                                            int* error_number) noexcept
    {
        for (;;)
        {
            errno = 0;
#ifdef F_DUPFD_CLOEXEC
            const int duplicate = fcntl(descriptor, F_DUPFD_CLOEXEC, 0);
#else
            const int duplicate = dup(descriptor);
#endif
            if (duplicate >= 0)
            {
#ifndef F_DUPFD_CLOEXEC
                errno = 0;
                if (fcntl(duplicate, F_SETFD, FD_CLOEXEC) != 0)
                {
                    *error_number = Meta_Effective_IO_Error(errno);
                    close(duplicate);
                    return -1;
                }
#endif
                return duplicate;
            }
            if (errno == EINTR) continue;
            *error_number = Meta_Effective_IO_Error(errno);
            return -1;
        }
    }
#endif

    bool Prepare_Fatal_Cleanup(int* error_number) noexcept
    {
#ifndef _WIN32
        try
        {
            fatal_cleanup_record_.temporary_basename = temporary_basename_;
        }
        catch (const std::length_error&)
        {
            *error_number = EOVERFLOW;
            return false;
        }
        catch (const std::bad_alloc&)
        {
            *error_number = ENOMEM;
            return false;
        }

        int duplicate_error = 0;
        const int cleanup_directory = Duplicate_Cleanup_Descriptor(
            directory_descriptor_, &duplicate_error);
        if (cleanup_directory < 0)
        {
            *error_number = duplicate_error;
            return false;
        }
        const int cleanup_temporary =
            Duplicate_Cleanup_Descriptor(fileno(stream_), &duplicate_error);
        if (cleanup_temporary < 0)
        {
            close(cleanup_directory);
            *error_number = duplicate_error;
            return false;
        }
        fatal_cleanup_record_.directory_descriptor = cleanup_directory;
        fatal_cleanup_record_.temporary_descriptor = cleanup_temporary;
        fatal_cleanup_record_.temporary_device = temporary_device_;
        fatal_cleanup_record_.temporary_inode = temporary_inode_;
#endif
        fatal_cleanup_record_.cleanup_error.store(0, std::memory_order_relaxed);
        fatal_cleanup_record_.state.store(fatal_cleanup_active,
                                          std::memory_order_release);
        fatal_cleanup_prepared_ = true;
        return true;
    }

    void Open()
    {
        if (final_name_.empty() || final_name_.find('\0') != std::string::npos)
        {
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, error_by_,
                "Metadynamics output path must be nonempty and contain no "
                "embedded null byte.");
            return;
        }
        try
        {
            final_path_ = Meta_Native_Path(controller_, final_name_, error_by_);
            parent_path_ = final_path_.parent_path();
            if (parent_path_.empty()) parent_path_ = fs::path(".");
            const fs::path filename = final_path_.filename();
            if (filename.empty() || filename == fs::path(".") ||
                filename == fs::path(".."))
            {
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand, error_by_,
                    "Metadynamics output path '%s' does not name a file.",
                    final_name_.c_str());
                return;
            }
#ifndef _WIN32
            final_basename_ = filename.native();
#endif
        }
        catch (const std::length_error&)
        {
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, error_by_,
                "Metadynamics output path '%s' exceeds the host path limit.",
                final_name_.c_str());
            return;
        }
        catch (const std::bad_alloc&)
        {
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorMallocFailed, error_by_,
                "Unable to allocate native path storage for metadynamics "
                "output '%s'.",
                final_name_.c_str());
            return;
        }
        catch (const std::exception& error)
        {
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Unable to prepare metadynamics output path '%s': %s.",
                final_name_.c_str(), error.what());
            return;
        }

#ifndef _WIN32
        errno = 0;
        directory_descriptor_ = open(parent_path_.native().c_str(),
                                     O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (directory_descriptor_ < 0)
        {
            const int open_error = Meta_Effective_IO_Error(errno);
            controller_->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by_,
                "Unable to open the parent directory of metadynamics output "
                "'%s': %s.",
                final_name_.c_str(), strerror(open_error));
            return;
        }
#endif

        for (std::uint64_t attempt = 0;; ++attempt)
        {
            int candidate_error = 0;
            try
            {
                temporary_basename_ = Meta_Random_Temporary_Basename(
                    controller_, error_by_, final_name_);
#ifdef _WIN32
                temporary_path_ =
                    parent_path_ / Meta_Native_Path(controller_,
                                                    temporary_basename_,
                                                    error_by_);
#endif
            }
            catch (const std::length_error&)
            {
                Close_And_Remove(false);
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow, error_by_,
                    "Unable to construct a temporary output name next to "
                    "'%s': the path exceeds the host string limit.",
                    final_name_.c_str());
                return;
            }
            catch (const std::bad_alloc&)
            {
                Close_And_Remove(false);
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorMallocFailed, error_by_,
                    "Unable to allocate a temporary output name next to "
                    "'%s'.",
                    final_name_.c_str());
                return;
            }

            errno = 0;
#ifdef _WIN32
            FILE* candidate_file = NULL;
            const errno_t status =
                _wfopen_s(&candidate_file, temporary_path_.c_str(), L"wx");
            candidate_error = status == 0 ? 0 : static_cast<int>(status);
#else
            const int candidate_descriptor =
                openat(directory_descriptor_, temporary_basename_.c_str(),
                       O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0666);
            candidate_error = Meta_Effective_IO_Error(errno);
            FILE* candidate_file = NULL;
            if (candidate_descriptor >= 0)
            {
                struct stat temporary_status = {};
                struct stat final_status = {};
                bool candidate_ready = true;
                errno = 0;
                if (fstat(candidate_descriptor, &temporary_status) != 0 ||
                    !S_ISREG(temporary_status.st_mode))
                {
                    candidate_error = Meta_Effective_IO_Error(errno);
                    candidate_ready = false;
                }
                errno = 0;
                if (candidate_ready &&
                    fstatat(directory_descriptor_, final_basename_.c_str(),
                            &final_status, AT_SYMLINK_NOFOLLOW) == 0)
                {
                    if (S_ISREG(final_status.st_mode))
                    {
                        replaced_mode_ = final_status.st_mode & 07777;
                        preserve_replaced_mode_ = true;
                        if (fchmod(candidate_descriptor, replaced_mode_) != 0)
                        {
                            candidate_error = Meta_Effective_IO_Error(errno);
                            candidate_ready = false;
                        }
                    }
                }
                else if (candidate_ready && errno != ENOENT)
                {
                    candidate_error = Meta_Effective_IO_Error(errno);
                    candidate_ready = false;
                }
                if (candidate_ready)
                {
                    temporary_device_ = temporary_status.st_dev;
                    temporary_inode_ = temporary_status.st_ino;
                    candidate_file = fdopen(candidate_descriptor, "w");
                }
                if (!candidate_ready || candidate_file == NULL)
                {
                    if (candidate_file == NULL && candidate_ready)
                    {
                        candidate_error = Meta_Effective_IO_Error(errno);
                    }
                    const int cleanup_error = Remove_Matching_Path(
                        directory_descriptor_, temporary_basename_,
                        candidate_descriptor, temporary_status.st_dev,
                        temporary_status.st_ino);
                    close(candidate_descriptor);
                    if (cleanup_error != 0 && candidate_error == 0)
                    {
                        candidate_error = cleanup_error;
                    }
                }
            }
#endif
            if (candidate_file != NULL)
            {
                stream_ = candidate_file;
                temporary_exists_ = true;
                int protection_error = 0;
                if (!Prepare_Fatal_Cleanup(&protection_error))
                {
                    const int cleanup_error = Close_And_Remove(false);
                    controller_->Throw_Formatted_SPONGE_Error(
                        spongeErrorOpenFileFailed, error_by_,
                        "Cannot prepare fatal cleanup for metadynamics output "
                        "'%s': %s%s%s.",
                        final_name_.c_str(), strerror(protection_error),
                        cleanup_error == 0
                            ? ""
                            : "; temporary cleanup also failed: ",
                        cleanup_error == 0 ? "" : strerror(cleanup_error));
                    return;
                }
                if (!CONTROLLER::Register_Fatal_Cleanup(Fatal_Cleanup,
                                                        &fatal_cleanup_record_))
                {
                    const int cleanup_error = Close_And_Remove(false);
                    controller_->Throw_Formatted_SPONGE_Error(
                        spongeErrorOverflow, error_by_,
                        "Cannot protect metadynamics output '%s': "
                        "fatal-cleanup "
                        "registration is unavailable (duplicate, capacity, or "
                        "fatal-shutdown state)%s%s.",
                        final_name_.c_str(),
                        cleanup_error == 0 ? ""
                                           : ", and temporary cleanup "
                                             "failed: ",
                        cleanup_error == 0 ? "" : strerror(cleanup_error));
                    return;
                }
                registered_ = true;
                return;
            }
            if (candidate_error != EEXIST)
            {
                Close_And_Remove(false);
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorOpenFileFailed, error_by_,
                    "Unable to create an exclusive temporary output next to "
                    "'%s': %s.",
                    final_name_.c_str(),
                    strerror(Meta_Effective_IO_Error(candidate_error)));
                return;
            }
            if (attempt == std::numeric_limits<std::uint64_t>::max())
            {
                Close_And_Remove(false);
                controller_->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow, error_by_,
                    "Exhausted the representable temporary-output name space "
                    "next to '%s'.",
                    final_name_.c_str());
                return;
            }
        }
    }

    int Remove_Temporary() noexcept
    {
        if (!temporary_exists_) return 0;
        errno = 0;
#ifdef _WIN32
        const int remove_status = _wremove(temporary_path_.c_str());
        const int remove_error = remove_status == 0 || errno == ENOENT
                                     ? 0
                                     : Meta_Effective_IO_Error(errno);
        if (remove_error == 0) temporary_exists_ = false;
        return remove_error;
#else
        if (stream_ == NULL)
        {
            temporary_exists_ = false;
            return ESTALE;
        }
        const int descriptor = fileno(stream_);
        const int remove_error = Remove_Matching_Path(
            directory_descriptor_, temporary_basename_, descriptor,
            temporary_device_, temporary_inode_);
        if (remove_error == 0 || remove_error == ENOENT)
        {
            temporary_exists_ = false;
            return 0;
        }
        temporary_exists_ = false;
        return remove_error;
#endif
    }

#ifndef _WIN32
    static bool Descriptor_Matches_Path(int directory_descriptor,
                                        const std::string& basename,
                                        int descriptor, dev_t expected_device,
                                        ino_t expected_inode,
                                        int* error_number) noexcept
    {
        if (directory_descriptor < 0 || descriptor < 0)
        {
            *error_number = EIO;
            return false;
        }
        struct stat descriptor_status = {};
        struct stat path_status = {};
        errno = 0;
        if (fstat(descriptor, &descriptor_status) != 0)
        {
            *error_number = Meta_Effective_IO_Error(errno);
            return false;
        }
        errno = 0;
        if (fstatat(directory_descriptor, basename.c_str(), &path_status,
                    AT_SYMLINK_NOFOLLOW) != 0)
        {
            *error_number = Meta_Effective_IO_Error(errno);
            return false;
        }
        if (descriptor_status.st_dev != expected_device ||
            descriptor_status.st_ino != expected_inode ||
            path_status.st_dev != expected_device ||
            path_status.st_ino != expected_inode ||
            !S_ISREG(descriptor_status.st_mode) ||
            !S_ISREG(path_status.st_mode))
        {
            *error_number = ESTALE;
            return false;
        }
        return true;
    }

    static int Remove_Matching_Path(int directory_descriptor,
                                    const std::string& basename, int descriptor,
                                    dev_t expected_device,
                                    ino_t expected_inode) noexcept
    {
        int identity_error = 0;
        if (!Descriptor_Matches_Path(directory_descriptor, basename, descriptor,
                                     expected_device, expected_inode,
                                     &identity_error))
        {
            return identity_error;
        }
        // POSIX has no conditional unlink-by-inode operation.  The fresh
        // comparison avoids acting on a stale cached identity; the random
        // basename and write permission on the parent directory remain the
        // security boundary for the comparison-to-unlink interval.
        errno = 0;
        if (unlinkat(directory_descriptor, basename.c_str(), 0) != 0)
        {
            return Meta_Effective_IO_Error(errno);
        }
        return 0;
    }

    bool Verify_Temporary_Identity(int* error_number) noexcept
    {
        if (!temporary_exists_ || stream_ == NULL || directory_descriptor_ < 0)
        {
            *error_number = EIO;
            return false;
        }
        const int descriptor = fileno(stream_);
        return Descriptor_Matches_Path(
            directory_descriptor_, temporary_basename_, descriptor,
            temporary_device_, temporary_inode_, error_number);
    }
#endif

    int Close_And_Remove(bool unregister) noexcept
    {
        int normal_expected = normal_cleanup_idle;
        if (!normal_cleanup_state_.compare_exchange_strong(
                normal_expected, normal_cleanup_busy, std::memory_order_acq_rel,
                std::memory_order_acquire))
        {
#ifdef SPONGE_META_ATOMIC_OUTPUT_PROBE
            if (probe_competing_cleanup_hook_ != NULL)
            {
                probe_competing_cleanup_hook_(this);
            }
#endif
            while (normal_cleanup_state_.load(std::memory_order_acquire) ==
                   normal_cleanup_busy)
            {
                std::this_thread::yield();
            }
            return normal_cleanup_result_.load(std::memory_order_acquire);
        }

        // Exactly one normal caller may ever touch registered_, stream_, or
        // the ordinary/cleanup descriptors.  In particular, fatal state
        // `claimed` coordinates the fatal callback but does not by itself
        // elect a unique normal follower.
        const auto finish_normal_cleanup = [this](const int result) noexcept
        {
            normal_cleanup_result_.store(result, std::memory_order_release);
            // This release publication is the owner's final object access.
            // A waiting caller may return (and permit object destruction) as
            // soon as it observes `done`.
            normal_cleanup_state_.store(normal_cleanup_done,
                                        std::memory_order_release);
            return result;
        };

#ifdef SPONGE_META_ATOMIC_OUTPUT_PROBE
        if (fatal_cleanup_prepared_ && probe_normal_cleanup_hook_ != NULL)
        {
            probe_normal_cleanup_hook_(this);
        }
#endif

        bool normal_cleanup_owner = !fatal_cleanup_prepared_;
        int observed_state = fatal_cleanup_disarmed;
        if (fatal_cleanup_prepared_)
        {
            int expected = fatal_cleanup_active;
            normal_cleanup_owner =
                fatal_cleanup_record_.state.compare_exchange_strong(
                    expected, fatal_cleanup_normal_busy,
                    std::memory_order_acq_rel, std::memory_order_acquire);
            observed_state =
                normal_cleanup_owner ? fatal_cleanup_normal_busy : expected;
        }

        if (!normal_cleanup_owner &&
            observed_state == fatal_cleanup_normal_busy)
        {
            while (fatal_cleanup_record_.state.load(
                       std::memory_order_acquire) == fatal_cleanup_normal_busy)
            {
                std::this_thread::yield();
            }
            const int cleanup_error = fatal_cleanup_record_.cleanup_error.load(
                std::memory_order_acquire);
            return finish_normal_cleanup(cleanup_error);
        }

        if (!normal_cleanup_owner && observed_state == fatal_cleanup_claimed)
        {
            if (unregister && registered_)
            {
                if (!CONTROLLER::Unregister_Fatal_Cleanup(
                        Fatal_Cleanup, &fatal_cleanup_record_))
                {
                    const int cleanup_error =
                        fatal_cleanup_record_.cleanup_error.load(
                            std::memory_order_acquire);
                    return finish_normal_cleanup(cleanup_error);
                }
                registered_ = false;
            }
            if (stream_ != NULL)
            {
                fclose(stream_);
                stream_ = NULL;
            }
#ifndef _WIN32
            if (directory_descriptor_ >= 0)
            {
                close(directory_descriptor_);
                directory_descriptor_ = -1;
            }
#endif
            temporary_exists_ = false;
            const int cleanup_error = fatal_cleanup_record_.cleanup_error.load(
                std::memory_order_acquire);
            return finish_normal_cleanup(cleanup_error);
        }

        int remove_error = 0;
#ifndef _WIN32
        if (normal_cleanup_owner) remove_error = Remove_Temporary();
#endif
        if (stream_ != NULL)
        {
            fclose(stream_);
            stream_ = NULL;
        }
#ifdef _WIN32
        if (normal_cleanup_owner) remove_error = Remove_Temporary();
#endif
#ifndef _WIN32
        if (directory_descriptor_ >= 0)
        {
            close(directory_descriptor_);
            directory_descriptor_ = -1;
        }
#endif
        if (normal_cleanup_owner && fatal_cleanup_prepared_)
        {
#ifndef _WIN32
            errno = 0;
            if (close(fatal_cleanup_record_.temporary_descriptor) != 0 &&
                remove_error == 0)
            {
                remove_error = Meta_Effective_IO_Error(errno);
            }
            errno = 0;
            if (close(fatal_cleanup_record_.directory_descriptor) != 0 &&
                remove_error == 0)
            {
                remove_error = Meta_Effective_IO_Error(errno);
            }
#endif
            fatal_cleanup_record_.cleanup_error.store(
                remove_error, std::memory_order_release);
            fatal_cleanup_record_.state.store(fatal_cleanup_disarmed,
                                              std::memory_order_release);
        }
        if (unregister && registered_)
        {
            if (CONTROLLER::Unregister_Fatal_Cleanup(Fatal_Cleanup,
                                                     &fatal_cleanup_record_))
            {
                registered_ = false;
            }
        }
        return finish_normal_cleanup(remove_error);
    }

    CONTROLLER* controller_ = NULL;
    std::string final_name_;
    const char* error_by_ = NULL;
    fs::path final_path_;
    fs::path parent_path_;
    std::string temporary_basename_;
    FILE* stream_ = NULL;
    bool registered_ = false;
    bool temporary_exists_ = false;
    bool fatal_cleanup_prepared_ = false;
    Fatal_Cleanup_Record fatal_cleanup_record_;
    std::atomic<int> normal_cleanup_state_{normal_cleanup_idle};
    std::atomic<int> normal_cleanup_result_{0};
#ifdef SPONGE_META_ATOMIC_OUTPUT_PROBE
    inline static Probe_Hook probe_post_rename_hook_ = NULL;
    inline static Probe_Hook probe_normal_cleanup_hook_ = NULL;
    inline static Probe_Hook probe_competing_cleanup_hook_ = NULL;
#endif
#ifdef _WIN32
    fs::path temporary_path_;
#else
    int directory_descriptor_ = -1;
    std::string final_basename_;
    dev_t temporary_device_ = 0;
    ino_t temporary_inode_ = 0;
    mode_t replaced_mode_ = 0;
    bool preserve_replaced_mode_ = false;
#endif
};

#ifndef SPONGE_META_ATOMIC_OUTPUT_PROBE
static std::vector<std::string> Tokenize_Meta_Line(const std::string& line)
{
    std::istringstream input(line);
    std::vector<std::string> tokens;
    std::string token;
    while (input >> token) tokens.push_back(token);
    return tokens;
}

static float Parse_Meta_File_Float(
    CONTROLLER* controller, const std::string& token, std::size_t line_number,
    const char* field_name, const char* error_by = "META::Read_Potential",
    const char* input_name = "Potential input")
{
    errno = 0;
    char* end = NULL;
    const float parsed = strtof(token.c_str(), &end);
    if (end == token.c_str() || end == NULL || *end != '\0' ||
        errno == ERANGE || !Float_Memory_Is_Zero_Or_Normal(&parsed))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "%s line %zu has invalid %s token '%s'; expected a "
            "finite zero-or-normal float.",
            input_name, line_number, field_name, token.c_str());
    }
    return parsed;
}

static int Parse_Meta_File_Int(CONTROLLER* controller, const std::string& token,
                               std::size_t line_number, const char* field_name,
                               const char* error_by = "META::Read_Potential",
                               const char* input_name = "Potential input")
{
    errno = 0;
    char* end = NULL;
    const long parsed = strtol(token.c_str(), &end, 10);
    if (end == token.c_str() || end == NULL || *end != '\0' ||
        errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "%s line %zu has invalid %s token '%s'; expected a "
            "base-10 int.",
            input_name, line_number, field_name, token.c_str());
    }
    return static_cast<int>(parsed);
}

static std::uint64_t Parse_Meta_File_Hex_U64(
    CONTROLLER* controller, const std::string& token, std::size_t line_number,
    const char* field_name, const char* error_by, const char* input_name)
{
    errno = 0;
    char* end = NULL;
    const unsigned long long parsed = strtoull(token.c_str(), &end, 16);
    if (end == token.c_str() || end == NULL || *end != '\0' ||
        errno == ERANGE || token[0] == '+' || token[0] == '-')
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "%s line %zu has invalid %s token '%s'; expected an unsigned "
            "64-bit hexadecimal integer.",
            input_name, line_number, field_name, token.c_str());
    }
    return static_cast<std::uint64_t>(parsed);
}

enum class MetaCoordinateFormat
{
    round_trip,
    legacy_potential,
    legacy_edge,
};

static bool Meta_Coordinate_Token_Matches(const std::string& token,
                                          float observed, float expected,
                                          MetaCoordinateFormat format)
{
    if (format == MetaCoordinateFormat::round_trip) return observed == expected;
    if (format == MetaCoordinateFormat::legacy_potential)
    {
        std::ostringstream canonical;
        canonical << expected;
        return token == canonical.str();
    }
    char canonical[64];
    const int written = snprintf(canonical, sizeof(canonical), "%f",
                                 static_cast<double>(expected));
    return written >= 0 &&
           static_cast<std::size_t>(written) < sizeof(canonical) &&
           token == canonical;
}

static bool Read_Nonempty_Meta_Line(std::ifstream& input, std::string* line,
                                    std::size_t* line_number)
{
    while (std::getline(input, *line))
    {
        ++*line_number;
        if (!Tokenize_Meta_Line(*line).empty()) return true;
    }
    return false;
}

static float Store_Checked_Meta_Float(CONTROLLER* controller, int error_number,
                                      const char* error_by,
                                      const char* quantity, double value);

static float Evaluate_Gaussian_Switch(const float rij, const float center,
                                      const float inv_w, const float period,
                                      float& df)
{
    float distant = rij - center;
    if (period != 0.0f)
    {
        distant -= roundf(distant / period) * period;
    }
    const float dx = distant * inv_w;
    const float f = expf(-dx * dx / 2.0f);
    df = -dx * inv_w * f;
    return f;
}

void MetaGrid::Initial(const std::vector<int>& npts,
                       const std::vector<float>& lo,
                       const std::vector<float>& up,
                       const std::vector<bool>& periodic,
                       const std::vector<float>& validated_spacing,
                       const std::vector<double>& validated_inverse_spacing)
{
    ndim = static_cast<int>(npts.size());
    num_points = npts;
    is_periodic = periodic;
    lower = lo;
    upper = up;
    spacing = validated_spacing;
    inv_spacing = validated_inverse_spacing;
    total_size = 1;
    for (int d = 0; d < ndim; ++d)
    {
        total_size *= num_points[d];
    }
}

void MetaGrid::Alloc_Device()
{
    if (!potential.empty() && d_potential == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_potential, potential.data(),
                                      sizeof(float) * potential.size());
    else if (!potential.empty())
        deviceMemcpy(d_potential, potential.data(),
                     sizeof(float) * potential.size(),
                     deviceMemcpyHostToDevice);
    if (!force.empty() && d_force == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_force, force.data(),
                                      sizeof(float) * force.size());
    else if (!force.empty())
        deviceMemcpy(d_force, force.data(), sizeof(float) * force.size(),
                     deviceMemcpyHostToDevice);
    if (!normal_lse.empty() && d_normal_lse == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_normal_lse, normal_lse.data(),
                                      sizeof(float) * normal_lse.size());
    else if (!normal_lse.empty())
        deviceMemcpy(d_normal_lse, normal_lse.data(),
                     sizeof(float) * normal_lse.size(),
                     deviceMemcpyHostToDevice);
    if (!normal_force.empty() && d_normal_force == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_normal_force,
                                      normal_force.data(),
                                      sizeof(float) * normal_force.size());
    else if (!normal_force.empty())
        deviceMemcpy(d_normal_force, normal_force.data(),
                     sizeof(float) * normal_force.size(),
                     deviceMemcpyHostToDevice);
    if (ndim > 0)
    {
        if (d_num_points == NULL)
            Device_Malloc_And_Copy_Safely(
                (void**)&d_num_points, num_points.data(), sizeof(int) * ndim);
        if (d_lower == NULL)
            Device_Malloc_And_Copy_Safely((void**)&d_lower, lower.data(),
                                          sizeof(float) * ndim);
        if (d_spacing == NULL)
            Device_Malloc_And_Copy_Safely((void**)&d_spacing, spacing.data(),
                                          sizeof(float) * ndim);
    }
}

void MetaGrid::Sync_To_Device()
{
    if (d_potential && !potential.empty())
        deviceMemcpy(d_potential, potential.data(),
                     sizeof(float) * potential.size(),
                     deviceMemcpyHostToDevice);
    if (d_force && !force.empty())
        deviceMemcpy(d_force, force.data(), sizeof(float) * force.size(),
                     deviceMemcpyHostToDevice);
    if (d_normal_lse && !normal_lse.empty())
        deviceMemcpy(d_normal_lse, normal_lse.data(),
                     sizeof(float) * normal_lse.size(),
                     deviceMemcpyHostToDevice);
    if (d_normal_force && !normal_force.empty())
        deviceMemcpy(d_normal_force, normal_force.data(),
                     sizeof(float) * normal_force.size(),
                     deviceMemcpyHostToDevice);
}

void MetaGrid::Sync_To_Host()
{
    if (d_potential && !potential.empty())
        deviceMemcpy(potential.data(), d_potential,
                     sizeof(float) * potential.size(),
                     deviceMemcpyDeviceToHost);
    if (d_force && !force.empty())
        deviceMemcpy(force.data(), d_force, sizeof(float) * force.size(),
                     deviceMemcpyDeviceToHost);
    if (d_normal_lse && !normal_lse.empty())
        deviceMemcpy(normal_lse.data(), d_normal_lse,
                     sizeof(float) * normal_lse.size(),
                     deviceMemcpyDeviceToHost);
    if (d_normal_force && !normal_force.empty())
        deviceMemcpy(normal_force.data(), d_normal_force,
                     sizeof(float) * normal_force.size(),
                     deviceMemcpyDeviceToHost);
}

int MetaGrid::Get_Flat_Index(const std::vector<float>& values) const
{
    int idx = 0;
    int fac = 1;
    for (int d = 0; d < ndim; ++d)
    {
        const double scaled =
            (static_cast<double>(values[d]) - lower[d]) * inv_spacing[d];
        int i = 0;
        if (is_periodic[d])
        {
            double wrapped = std::fmod(scaled, num_points[d]);
            if (wrapped < 0.0) wrapped += num_points[d];
            i = static_cast<int>(std::floor(wrapped));
            if (i >= num_points[d]) i = 0;
        }
        else if (scaled >= num_points[d])
        {
            i = num_points[d] - 1;
        }
        else if (scaled > 0.0)
        {
            i = static_cast<int>(std::floor(scaled));
        }
        idx += i * fac;
        fac *= num_points[d];
    }
    return idx;
}

std::vector<float> MetaGrid::Get_Coordinates(int flat_index) const
{
    std::vector<float> coords(ndim);
    Get_Coordinates(flat_index, coords.data());
    return coords;
}

void MetaGrid::Get_Coordinates(int flat_index, float* coordinates) const
{
    for (int d = 0; d < ndim; ++d)
    {
        int i = flat_index % num_points[d];
        flat_index /= num_points[d];
        coordinates[d] =
            static_cast<float>(static_cast<double>(lower[d]) +
                               (static_cast<double>(i) + 0.5) * spacing[d]);
    }
}

void MetaScatter::Initial(const std::vector<int>& npts,
                          const std::vector<float>& period,
                          const std::vector<std::vector<float>>& coor)
{
    ndim = static_cast<int>(npts.size());
    num_points = static_cast<int>(coor.size());
    coordinates = coor;
    periods = period;
    coordinates_flat.resize(num_points * ndim);
    for (int index = 0; index < num_points; ++index)
    {
        for (int d = 0; d < ndim; ++d)
        {
            coordinates_flat[index * ndim + d] = coordinates[index][d];
        }
    }
}

void MetaScatter::Alloc_Device()
{
    if (!coordinates_flat.empty() && d_coordinates == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_coordinates,
                                      coordinates_flat.data(),
                                      sizeof(float) * coordinates_flat.size());
    else if (!coordinates_flat.empty())
        deviceMemcpy(d_coordinates, coordinates_flat.data(),
                     sizeof(float) * coordinates_flat.size(),
                     deviceMemcpyHostToDevice);
    if (!periods.empty() && d_periods == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_periods, periods.data(),
                                      sizeof(float) * periods.size());
    else if (!periods.empty())
        deviceMemcpy(d_periods, periods.data(), sizeof(float) * periods.size(),
                     deviceMemcpyHostToDevice);
    if (!potential.empty() && d_potential == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_potential, potential.data(),
                                      sizeof(float) * potential.size());
    else if (!potential.empty())
        deviceMemcpy(d_potential, potential.data(),
                     sizeof(float) * potential.size(),
                     deviceMemcpyHostToDevice);
    if (!force.empty() && d_force == NULL)
        Device_Malloc_And_Copy_Safely((void**)&d_force, force.data(),
                                      sizeof(float) * force.size());
    else if (!force.empty())
        deviceMemcpy(d_force, force.data(), sizeof(float) * force.size(),
                     deviceMemcpyHostToDevice);
}

void MetaScatter::Sync_To_Device()
{
    if (d_potential && !potential.empty())
        deviceMemcpy(d_potential, potential.data(),
                     sizeof(float) * potential.size(),
                     deviceMemcpyHostToDevice);
    if (d_force && !force.empty())
        deviceMemcpy(d_force, force.data(), sizeof(float) * force.size(),
                     deviceMemcpyHostToDevice);
}

void MetaScatter::Sync_To_Host()
{
    if (d_potential && !potential.empty())
        deviceMemcpy(potential.data(), d_potential,
                     sizeof(float) * potential.size(),
                     deviceMemcpyDeviceToHost);
    if (d_force && !force.empty())
        deviceMemcpy(force.data(), d_force, sizeof(float) * force.size(),
                     deviceMemcpyDeviceToHost);
}

int MetaScatter::Get_Index(const std::vector<float>& values) const
{
    float min_dist = std::numeric_limits<float>::max();
    int min_idx = 0;
    for (int i = 0; i < num_points; ++i)
    {
        float dist = 0;
        for (int d = 0; d < ndim; ++d)
        {
            float diff = values[d] - coordinates[i][d];
            if (periods[d] > 0)
            {
                diff -= std::round(diff / periods[d]) * periods[d];
            }
            dist += diff * diff;
        }
        if (dist < min_dist)
        {
            min_dist = dist;
            min_idx = i;
        }
    }
    return min_idx;
}

std::vector<int> MetaScatter::Get_Neighbor(const std::vector<float>& values,
                                           const float* cutoff) const
{
    std::vector<int> neighbors;
    for (int i = 0; i < num_points; ++i)
    {
        bool within = true;
        for (int d = 0; d < ndim; ++d)
        {
            float diff = values[d] - coordinates[i][d];
            if (periods[d] > 0)
            {
                diff -= std::round(diff / periods[d]) * periods[d];
            }
            if (std::fabs(diff) > cutoff[d])
            {
                within = false;
                break;
            }
        }
        if (within)
        {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

const std::vector<float>& MetaScatter::Get_Coordinate(int index) const
{
    return coordinates[index];
}

static std::vector<float> normalize(const std::vector<float>& v)
{
    float norm = 0.;
    for (auto vi : v)
    {
        norm += vi * vi;
    }
    if (norm == 0.0)
    {
        throw std::runtime_error("Zero-length vector cannot be normalized.");
    }
    std::vector<float> new_v;
    for (int i = 0; i < v.size(); ++i)
    {
        new_v.push_back(v[i] / sqrt(norm));
    }
    return new_v;
}

static std::vector<float> cross_product(const std::vector<float>& a,
                                        const std::vector<float>& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

static float determinant(const std::vector<std::vector<float>>& matrix)
{
    int n = matrix.size();
    if (n == 1)
    {
        return matrix[0][0];
    }

    float det = 0;
    for (int i = 0; i < n; ++i)
    {
        std::vector<std::vector<float>> submatrix(n - 1,
                                                  std::vector<float>(n - 1));

        for (int j = 1; j < n; ++j)
        {
            int subCol = 0;
            for (int k = 0; k < n; ++k)
            {
                if (k == i) continue;
                submatrix[j - 1][subCol] = matrix[j][k];
                subCol++;
            }
        }

        float subDet = determinant(submatrix);
        det += (i % 2 == 0 ? 1 : -1) * matrix[0][i] * subDet;
    }

    return det;
}

META::Axis META::Rotate_Vector(const Axis& tang_vector)
{
    std::vector<float> normal_vector;
    int reference_axis = 0;
    if (fabs(tang_vector[reference_axis]) > 0.99)
    {
        ++reference_axis;
    }
    for (int i = 0; i < ndim; ++i)
    {
        if (i == reference_axis)
        {
            normal_vector.push_back(1.);
        }
        else
        {
            normal_vector.push_back(0.);
        }
    }
    Axis jb;
    float i_min = tang_vector[reference_axis];
    float e1 = sqrtf(1 - i_min * i_min);
    float e2 = -i_min / e1;
    for (int i = 0; i < ndim; ++i)
    {
        if (i == reference_axis)
        {
            jb.push_back(e1);
        }
        else
        {
            jb.push_back(tang_vector[i] * e2);
        }
    }
    if (ndim == 2)
    {
        std::vector<std::vector<float>> determinant_v =
            std::vector<Axis>{tang_vector, jb};
        float sign = determinant(determinant_v);
        return Axis{jb[0] * sign, jb[1] * sign};
    }
    return jb;
}

void META::Cartesian_To_Path(const Axis& Cartesian_values, Axis& Path_values)
{
    double cumulative_s = 0.0;
    Axis tang_vector(ndim, 0.);
    int index = mscatter->Get_Index(Cartesian_values);
    const Axis& values = (index < scatter_size - 1)
                             ? mscatter->Get_Coordinate(index)
                             : mscatter->Get_Coordinate(index - 1);
    const Axis& neighbor = (index < scatter_size - 1)
                               ? mscatter->Get_Coordinate(index + 1)
                               : mscatter->Get_Coordinate(index);
    Tang_Vector(tang_vector, values, neighbor);
    double projected_last =
        Project_To_Path(tang_vector, neighbor, Cartesian_values);
    Path_values.push_back(cumulative_s + projected_last);
    Axis normal_vector = Rotate_Vector(tang_vector);
    Path_values.push_back(
        Project_To_Path(normal_vector, values, Cartesian_values));
    if (ndim == 3)
    {
        Axis binormal_vector =
            normalize(cross_product(tang_vector, normal_vector));
        Path_values.push_back(
            Project_To_Path(binormal_vector, values, Cartesian_values));
    }
}

float META::Project_To_Path(const Gdata& tang_vector, const Axis& values,
                            const Axis& Cartesian)
{
    float projected_s = 0.;
    for (int i = 0; i < ndim; ++i)
    {
        projected_s += (Cartesian[i] - values[i]) * tang_vector[i];
    }
    return projected_s;
}

double META::Tang_Vector(Gdata& tang_vector, const Axis& values,
                         const Axis& neighbor)
{
    double square = 0;
    for (int i = 0; i < ndim; ++i)
    {
        double distance = neighbor[i] - values[i];
        tang_vector[i] = distance;
        square += distance * distance;
    }
    double segment_s = sqrt(square);
    for (int i = 0; i < ndim; ++i)
    {
        tang_vector[i] /= segment_s;
    }
    return segment_s;
}

void META::Set_Grid(CONTROLLER* controller)  //
{
    std::vector<int> ngrid;
    std::vector<float> lower, upper, periodic, spacing;
    std::vector<double> inverse_spacing;
    std::vector<bool> isperiodic;
    if (ndim <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "META::Set_Grid",
            "Metadynamics requires at least one CV dimension.");
    }
    border_upper.resize(ndim);
    border_lower.resize(ndim);
    est_values_.resize(ndim);
    est_sum_force_.resize(ndim);
    if (n_grids == NULL || cv_mins == NULL || cv_maxs == NULL ||
        cv_periods == NULL || cutoff == NULL ||
        cv_deltas.size() != static_cast<std::size_t>(ndim) ||
        sigmas.size() != static_cast<std::size_t>(ndim) ||
        periods.size() != static_cast<std::size_t>(ndim))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "META::Set_Grid",
            "Metadynamics grid storage is incomplete before construction.");
    }
    std::size_t checked_grid_size = 1;
    for (int d = 0; d < ndim; ++d)
    {
        if (n_grids[d] <= 0 ||
            checked_grid_size > static_cast<std::size_t>(INT_MAX) /
                                    static_cast<std::size_t>(n_grids[d]))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "META::Set_Grid",
                "Metadynamics grid shape overflows its int index space at "
                "dimension %d (current product %zu, extent %d).",
                d, checked_grid_size, n_grids[d]);
        }
        checked_grid_size *= static_cast<std::size_t>(n_grids[d]);

        if (!Float_Memory_Is_Normal(&cv_deltas[d]) || !(cv_deltas[d] > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "META::Set_Grid",
                "Metadynamics grid spacing for dimension %d is not a finite "
                "positive normal float.",
                d);
        }
        const double inverse = 1.0 / static_cast<double>(cv_deltas[d]);
        if (!Double_Memory_Is_Finite(&inverse) || !(inverse > 0.0))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "META::Set_Grid",
                "The inverse grid spacing for metadynamics dimension %d is "
                "not finite and positive.",
                d);
        }
        spacing.push_back(cv_deltas[d]);
        inverse_spacing.push_back(inverse);
    }
    if (checked_grid_size >
        static_cast<std::size_t>(INT_MAX) / static_cast<std::size_t>(ndim))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "META::Set_Grid",
            "Metadynamics grid has %zu points and %d dimensions; the "
            "flattened force index exceeds INT_MAX.",
            checked_grid_size, ndim);
    }
    const std::size_t checked_force_size =
        checked_grid_size * static_cast<std::size_t>(ndim);

    float stored_log_normalization = 0.0f;
    if (usegrid)
    {
        double normalization = 1.0;
        const double sqrtpi = std::sqrt(static_cast<double>(CONSTANT_Pi));
        for (int i = 0; i < ndim; i++)
        {
            normalization /=
                static_cast<double>(cv_deltas[i]) * sigmas[i] / sqrtpi;
        }
        if (!Double_Memory_Is_Finite(&normalization) || !(normalization > 0.0))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "META::Set_Grid",
                "Metadynamics grid normalization is not finite and positive.");
        }
        const double log_normalization = std::log(normalization);
        if (!Double_Memory_Is_Finite(&log_normalization) ||
            log_normalization > std::numeric_limits<float>::max() ||
            log_normalization < -std::numeric_limits<float>::max())
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "META::Set_Grid",
                "The logarithm of the metadynamics grid normalization is "
                "outside the finite float range.");
        }
        stored_log_normalization = static_cast<float>(log_normalization);
        if (!Float_Memory_Is_Zero_Or_Normal(&stored_log_normalization))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "META::Set_Grid",
                "The logarithm of the metadynamics grid normalization is not "
                "representable as a zero-or-normal float.");
        }
    }
    else if (use_scatter)
    {
        if (scatter_size <= 0 ||
            static_cast<std::size_t>(scatter_size) > checked_grid_size ||
            static_cast<std::size_t>(scatter_size) >
                static_cast<std::size_t>(INT_MAX) /
                    static_cast<std::size_t>(ndim) ||
            tcoor.size() != static_cast<std::size_t>(ndim))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "META::Set_Grid",
                "Metadynamics scatter shape is invalid: %d points, %d "
                "dimensions, grid capacity %zu.",
                scatter_size, ndim, checked_grid_size);
        }
        for (int d = 0; d < ndim; ++d)
        {
            if (tcoor[d] == NULL)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorValueErrorCommand, "META::Set_Grid",
                    "Metadynamics scatter coordinate dimension %d is null.", d);
            }
            for (int index = 0; index < scatter_size; ++index)
            {
                if (!Float_Memory_Is_Zero_Or_Normal(tcoor[d] + index))
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorValueErrorCommand, "META::Set_Grid",
                        "Metadynamics scatter coordinate (%d, %d) is not a "
                        "finite zero-or-normal float.",
                        index, d);
                }
            }
        }
        if (static_cast<std::size_t>(scatter_size) == checked_grid_size)
        {
            double normalization = 1.0;
            const double sqrt_two_pi =
                std::sqrt(static_cast<double>(CONSTANT_Pi) * 2.0);
            for (int d = 0; d < ndim; ++d)
            {
                normalization /=
                    static_cast<double>(cv_deltas[d]) * sigmas[d] / sqrt_two_pi;
            }
            if (!Double_Memory_Is_Finite(&normalization) ||
                !(normalization > 0.0))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow, "META::Set_Grid",
                    "Metadynamics edge normalization is not finite and "
                    "positive.");
            }
            const double log_normalization = std::log(normalization);
            if (!Double_Memory_Is_Finite(&log_normalization) ||
                log_normalization > std::numeric_limits<float>::max() ||
                log_normalization < -std::numeric_limits<float>::max())
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow, "META::Set_Grid",
                    "The logarithm of the metadynamics edge normalization is "
                    "outside the finite float range.");
            }
            const float stored_log = static_cast<float>(log_normalization);
            if (!Float_Memory_Is_Zero_Or_Normal(&stored_log))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow, "META::Set_Grid",
                    "The logarithm of the metadynamics edge normalization is "
                    "not representable as a zero-or-normal float.");
            }
        }
    }
    else
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "META::Set_Grid",
            "Metadynamics has neither a grid nor a scatter representation.");
    }

    Device_Malloc_Safely((void**)&d_hill_centers, sizeof(float) * ndim);
    Device_Malloc_Safely((void**)&d_hill_inv_w, sizeof(float) * ndim);
    Device_Malloc_Safely((void**)&d_hill_periods, sizeof(float) * ndim);
    Device_Malloc_And_Copy_Safely((void**)&d_cutoff, cutoff,
                                  sizeof(float) * ndim);
    for (size_t i = 0; i < ndim; ++i)
    {
        ngrid.push_back(n_grids[i]);
        lower.push_back(cv_mins[i]);
        upper.push_back(cv_maxs[i]);
        periodic.push_back(cv_periods[i]);
        isperiodic.push_back(cv_periods[i] > 0 ? true : false);
    }
    mgrid = new MetaGrid();
    mgrid->Initial(ngrid, lower, upper, isperiodic, spacing, inverse_spacing);
    if (mgrid->total_size != static_cast<int>(checked_grid_size))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "META::Set_Grid",
            "Metadynamics grid size changed after validated construction.");
    }
    mgrid->normal_force.assign(checked_force_size, 0.0f);
    mgrid->normal_lse.assign(mgrid->total_size, 0.0f);
    mgrid->potential.assign(mgrid->total_size, 0.0f);
    if (usegrid)
    {
        mgrid->force.assign(checked_force_size, 0.0f);
        mgrid->normal_lse.assign(mgrid->total_size, stored_log_normalization);
        mscatter = NULL;
        Sum_Hills(history_freq);
        mgrid->Alloc_Device();
    }
    else if (use_scatter)
    {
        if (mask > 0)
        {
            mgrid->force.assign(checked_force_size, 0.0f);
        }
        std::vector<int> nscatter;
        for (size_t i = 0; i < ndim; ++i)
        {
            nscatter.push_back(n_grids[i]);
        }
        max_index = scatter_size / 2;
        std::vector<std::vector<float>> coor;
        for (size_t j = 0; j < scatter_size; ++j)
        {
            std::vector<float> p;
            for (size_t i = 0; i < ndim; ++i)
            {
                p.push_back(tcoor[i][j]);
            }
            coor.push_back(p);
        }
        mscatter = new MetaScatter();
        mscatter->Initial(nscatter, periodic, coor);
        mscatter->force.assign(static_cast<std::size_t>(scatter_size) * ndim,
                               0.0f);
        mscatter->potential.assign(scatter_size, 0.0f);
        Edge_Effect(1, scatter_size);
        Sum_Hills(history_freq);
        mgrid->Alloc_Device();
        mscatter->Alloc_Device();
    }
    else
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "META::Set_Grid",
            "Metadynamics has neither a grid nor a scatter representation.");
    }
}
void META::Estimate(const Axis& values, const bool need_potential,
                    const bool need_force)
{
    potential_local = 0;
    potential_backup = 0;

    double derived_shift = static_cast<double>(potential_max) +
                           static_cast<double>(dip) * CONSTANT_kB * temperature;
    if (do_negative)
    {
        if (grw)
        {
            derived_shift = (static_cast<double>(welltemp_factor) + dip) *
                            CONSTANT_kB * temperature;
        }
    }
    const float shift = Store_Checked_Meta_Float(
        controller, spongeErrorSimulationBreakDown, "META::Estimate",
        "sink normalization factor", derived_shift);
    const int nf_idx = mgrid->Get_Flat_Index(values);
    if (do_negative)
    {
        normalization_base_factor = shift;
        int reference_index = nf_idx;
        if (usegrid)
        {
            reference_index = 0;
        }
        else if (convmeta)
        {
            if (mscatter == NULL || max_index < 0 || max_index >= scatter_size)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "META::Estimate",
                    "Sink normalization has an invalid maximum scatter index.");
            }
            reference_index =
                mgrid->Get_Flat_Index(mscatter->Get_Coordinate(max_index));
        }
        normalization_reference_lse = mgrid->normal_lse[reference_index];
        if (!Float_Memory_Is_Zero_Or_Normal(&normalization_reference_lse))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Estimate",
                "Sink normalization references an invalid log normalization.");
        }
        new_max = subhill ? Normalization(values, shift, true) : 0.0f;
    }
    double force_max = 0.0;
    float normalforce_sum = 0.0;
    for (size_t i = 0; i < ndim; ++i)
    {
        est_sum_force_[i] = 0.0f;
    }
    for (size_t i = 0; i < ndim; ++i)
    {
        Dpotential_local[i] = 0.0;
        force_max += fabs(mgrid->normal_force[nf_idx * ndim + i]);
    }
    if (force_max > max_force && need_force && mask)
    {
        exit_tag += 1.0;
    }
    if (use_scatter)
    {
        if (subhill)
        {
            Hill hill = Hill(values, sigmas, periods, 1.0);
            std::vector<int> indices;
            if (do_cutoff)
            {
                indices = mscatter->Get_Neighbor(values, cutoff);
            }
            else
            {
                indices = std::vector<int>(scatter_size);
                std::iota(indices.begin(), indices.end(), 0);
            }
            for (auto index : indices)
            {
                const Axis& neighbor = mscatter->Get_Coordinate(index);
                const Gdata& tder = hill.Calc_Hill(neighbor);
                normalforce_sum += hill.potential;
                float factor =
                    (mask > 0)
                        ? mgrid->potential[mgrid->Get_Flat_Index(neighbor)]
                        : mscatter->potential[index];
                if (need_force)
                {
                    for (size_t i = 0; i < ndim; ++i)
                    {
                        est_sum_force_[i] += tder[i];
                        Dpotential_local[i] -= (factor)*tder[i];
                    }
                }
                potential_backup += factor * hill.potential;
            }
        }
        else
        {
            int sidx = mscatter->Get_Index(values);
            potential_backup =
                (mask > 0) ? mgrid->potential[mgrid->Get_Flat_Index(values)]
                           : mscatter->potential[sidx];
            potential_local = potential_backup - Calc_V_Shift(values);
            if (need_force)
            {
                int fidx = (mask > 0) ? mgrid->Get_Flat_Index(values) : sidx;
                for (int i = 0; i < cvs.size(); ++i)
                {
                    Dpotential_local[i] +=
                        (mask > 0) ? mgrid->force[fidx * ndim + i]
                                   : mscatter->force[fidx * ndim + i];
                }
            }
        }
    }
    else if (usegrid)
    {
        if (subhill)
        {
            Hill hill = Hill(values, sigmas, periods, 1.0);
            Axis vminus(ndim), vplus(ndim);
            for (size_t i = 0; i < ndim; ++i)
            {
                float lower = values[i] - cutoff[i];
                float upper = values[i] + cutoff[i] + 0.000001;
                if (periods[i] > 0)
                {
                    vminus[i] = lower;
                    vplus[i] = upper;
                }
                else
                {
                    vminus[i] = std::fmax(lower, cv_mins[i]);
                    vplus[i] = std::fmin(upper, cv_maxs[i]);
                }
            }
            Axis loop_flag = vminus;
            int index = 0;
            while (index >= 0)
            {
                const Gdata& tder = hill.Calc_Hill(loop_flag);
                float factor =
                    mgrid->potential[mgrid->Get_Flat_Index(loop_flag)];
                potential_backup += factor * hill.potential;
                if (need_force)
                {
                    for (size_t i = 0; i < ndim; ++i)
                    {
                        Dpotential_local[i] -= (factor - new_max) * tder[i];
                    }
                }
                index = ndim - 1;
                while (index >= 0)
                {
                    loop_flag[index] += cv_deltas[index];
                    if (loop_flag[index] > vplus[index])
                    {
                        loop_flag[index] = vminus[index];
                        --index;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
        else
        {
            int gidx = mgrid->Get_Flat_Index(values);
            potential_backup = mgrid->potential[gidx];
            if (need_force)
            {
                for (int i = 0; i < cvs.size(); ++i)
                {
                    Dpotential_local[i] += mgrid->force[gidx * ndim + i];
                }
            }
        }
        if (do_borderwall)
        {
            for (size_t i = 0; i < ndim; ++i)
            {
                border_upper[i] = cv_maxs[i] - cutoff[i];
                border_lower[i] = cv_mins[i] + cutoff[i];
            }
        }
    }
    if (need_potential)
    {
        potential_local = potential_backup - Calc_V_Shift(values);
    }
    if (need_force)
    {
        if (subhill)
        {
            if (!convmeta)
            {
                if (!Float_Memory_Is_Normal(&normalforce_sum) ||
                    !(normalforce_sum > 0.0f))
                {
                    controller->Throw_SPONGE_Error(
                        spongeErrorSimulationBreakDown, "META::Estimate",
                        "Subhill normalization sum is not finite and "
                        "positive.");
                }
                new_max = Store_Checked_Meta_Float(
                    controller, spongeErrorSimulationBreakDown,
                    "META::Estimate", "subhill normalization factor",
                    static_cast<double>(shift) / normalforce_sum);
            }
            for (int i = 0; i < cvs.size(); ++i)
            {
                Dpotential_local[i] = Store_Checked_Meta_Float(
                    controller, spongeErrorSimulationBreakDown,
                    "META::Estimate", "subhill potential derivative",
                    static_cast<double>(Dpotential_local[i]) +
                        static_cast<double>(new_max) * est_sum_force_[i]);
            }
        }
        else if (do_negative)
        {
            const double exponent =
                static_cast<double>(mgrid->normal_lse[nf_idx]) -
                normalization_reference_lse;
            const double normalization_scale = std::exp(exponent);
            if (!Double_Memory_Is_Finite(&exponent) ||
                !Double_Memory_Is_Finite(&normalization_scale) ||
                !(normalization_scale > 0.0))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "META::Estimate",
                    "Sink force normalization is not representable in stable "
                    "log space.");
            }
            const double force_factor =
                normalization_base_factor * normalization_scale;
            for (int i = 0; i < cvs.size(); ++i)
            {
                Dpotential_local[i] = Store_Checked_Meta_Float(
                    controller, spongeErrorSimulationBreakDown,
                    "META::Estimate", "sink potential derivative",
                    static_cast<double>(Dpotential_local[i]) +
                        force_factor * mgrid->normal_force[nf_idx * ndim + i]);
            }
        }
    }
    return;
}

static void Write_CV_Header(FILE* output, int ndim, const CV_LIST& cvs)
{
    for (int i = 0; i < ndim; ++i)
    {
        const char* name = NULL;
        if (i < static_cast<int>(cvs.size()) && cvs[i] != NULL &&
            !cvs[i]->module_name.empty())
        {
            name = cvs[i]->module_name.c_str();
        }
        if (name != NULL)
            fprintf(output, "%s\t", name);
        else
            fprintf(output, "cv%d\t", i + 1);
    }
}

void META::Write_Potential(void)
{
    if (!is_initialized) return;
    if (CONTROLLER::MPI_rank != CONTROLLER::MPI_size - 1) return;
    try
    {
        Axis point;
        Axis lower;
        if (subhill || (!usegrid && !use_scatter))
        {
            point.assign(ndim, 0.0f);
            lower.assign(ndim, 0.0f);
        }
        else if (mgrid != NULL)
        {
            point.resize(ndim);
        }
        Meta_Atomic_Output atomic_output(controller, write_potential_file_name,
                                         "META::Write_Potential");
        FILE* output = atomic_output.Stream();
        const int precision = std::numeric_limits<float>::max_digits10;
        if (subhill || (!usegrid && !use_scatter))
        {
            fprintf(output, "# ");
            Write_CV_Header(output, ndim, cvs);
            fprintf(output, "potential_local\tpotential_backup%s\n",
                    kde ? "" : "\tpotential_raw");
            for (int d = 0; d < ndim; ++d)
            {
                lower[d] = cv_mins[d] + 0.5f * cv_deltas[d];
                point[d] = lower[d];
            }
            int d = 0;
            while (d >= 0)
            {
                Estimate(point, true, false);
                for (float value : point)
                    fprintf(output, "%.*g\t", precision,
                            static_cast<double>(value));
                fprintf(output, "%.*g\t%.*g", precision,
                        static_cast<double>(potential_local), precision,
                        static_cast<double>(potential_backup));
                if (!kde && mgrid != NULL)
                    fprintf(
                        output, "\t%.*g", precision,
                        static_cast<double>(
                            mgrid->potential[mgrid->Get_Flat_Index(point)]));
                else if (!kde && mscatter != NULL)
                    fprintf(
                        output, "\t%.*g", precision,
                        static_cast<double>(
                            mscatter->potential[mscatter->Get_Index(point)]));
                fprintf(output, "\n");
                d = ndim - 1;
                while (d >= 0)
                {
                    point[d] += cv_deltas[d];
                    if (point[d] > cv_maxs[d])
                    {
                        point[d] = lower[d];
                        --d;
                    }
                    else
                        break;
                }
            }
        }
        else if (mscatter != NULL)
        {
            fprintf(output, "# ");
            Write_CV_Header(output, ndim, cvs);
            fprintf(output, "potential_raw\tpotential_shifted\n");
            for (int index = 0; index < scatter_size; ++index)
            {
                const Axis& point = mscatter->Get_Coordinate(index);
                const float shift = Calc_V_Shift(point);
                for (float value : point)
                    fprintf(output, "%.*g\t", precision,
                            static_cast<double>(value));
                const double shifted =
                    static_cast<double>(mscatter->potential[index]) - shift;
                if (!Double_Memory_Is_Finite(&shifted))
                    controller->Throw_SPONGE_Error(
                        spongeErrorSimulationBreakDown, "META::Write_Potential",
                        "Shifted scatter potential is non-finite.");
                fprintf(output, "%.*g\t%.*g\n", precision,
                        static_cast<double>(mscatter->potential[index]),
                        std::numeric_limits<double>::max_digits10, shifted);
            }
        }
        else if (mgrid != NULL)
        {
            fprintf(output, "# ");
            Write_CV_Header(output, ndim, cvs);
            fprintf(output, "potential_raw\tpotential_shifted\tvshift\n");
            for (int index = 0; index < mgrid->total_size; ++index)
            {
                mgrid->Get_Coordinates(index, point.data());
                const float shift = Calc_V_Shift(point);
                for (float value : point)
                    fprintf(output, "%.*g\t", precision,
                            static_cast<double>(value));
                const double shifted =
                    static_cast<double>(mgrid->potential[index]) - shift;
                if (!Double_Memory_Is_Finite(&shifted))
                    controller->Throw_SPONGE_Error(
                        spongeErrorSimulationBreakDown, "META::Write_Potential",
                        "Shifted grid potential is non-finite.");
                fprintf(output, "%.*g\t%.*g\t%.*g\n", precision,
                        static_cast<double>(mgrid->potential[index]),
                        std::numeric_limits<double>::max_digits10, shifted,
                        precision, static_cast<double>(shift));
            }
        }
        atomic_output.Commit();
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMallocFailed, "META::Write_Potential",
            "Unable to allocate bounded potential-output working storage.");
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "META::Write_Potential",
            "Potential-output working storage exceeds the host limit.");
    }
}

void META::Write_Directly(void)
{
    if (!is_initialized || !(use_scatter || usegrid)) return;
    if (CONTROLLER::MPI_rank != CONTROLLER::MPI_size - 1) return;
    try
    {
        Axis grid_point;
        if (mscatter == NULL) grid_point.resize(ndim);
        Meta_Atomic_Output atomic_output(controller, write_directly_file_name,
                                         "META::Write_Directly");
        FILE* output = atomic_output.Stream();
        const int precision = std::numeric_limits<float>::max_digits10;
        auto write_value = [&](float value, const char* field)
        {
            if (!Float_Memory_Is_Zero_Or_Normal(&value))
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "META::Write_Directly",
                    "Cannot persist a non-finite or subnormal %s.", field);
            fprintf(output, "%.*g", precision, static_cast<double>(value));
        };
        fprintf(output, "SPONGE_META_POTENTIAL_V1 %d %s %s\n", ndim,
                mscatter != NULL ? "scatter" : "grid",
                subhill ? "subhill" : "d_force");
        for (int d = 0; d < ndim; ++d)
        {
            write_value(cv_mins[d], "CV minimum");
            fprintf(output, "\t");
            write_value(cv_maxs[d], "CV maximum");
            fprintf(output, "\t");
            write_value(cv_deltas[d], "CV spacing");
            fprintf(output, "\n");
        }
        for (int d = 0; d < ndim; ++d) fprintf(output, "%d\t", n_grids[d]);
        const int count = mscatter != NULL ? scatter_size : mgrid->total_size;
        fprintf(output, "%d\n", count);
        for (int index = 0; index < count; ++index)
        {
            const Axis* point = NULL;
            if (mscatter != NULL)
                point = &mscatter->Get_Coordinate(index);
            else
            {
                mgrid->Get_Coordinates(index, grid_point.data());
                point = &grid_point;
            }
            Estimate(*point, true, false);
            for (float value : *point)
            {
                write_value(value, "potential coordinate");
                fprintf(output, "\t");
            }
            write_value(potential_local, "evaluated potential");
            fprintf(output, "\t");
            if (subhill)
            {
                write_value(potential_backup, "backup potential");
                fprintf(output, "\t");
            }
            else
            {
                const float* force = mscatter != NULL
                                         ? &mscatter->force[index * ndim]
                                         : &mgrid->force[index * ndim];
                for (int d = 0; d < ndim; ++d)
                {
                    write_value(force[d], "potential derivative");
                    fprintf(output, "\t");
                }
            }
            write_value(mscatter != NULL ? mscatter->potential[index]
                                         : mgrid->potential[index],
                        "raw potential");
            fprintf(output, "\n");
        }
        atomic_output.Commit();
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMallocFailed, "META::Write_Directly",
            "Unable to allocate bounded restart-output working storage.");
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "META::Write_Directly",
            "Restart-output working storage exceeds the host limit.");
    }
}
void META::Read_Potential(CONTROLLER* controller)
{
    std::ifstream input(Meta_Native_Path(controller, read_potential_file_name,
                                         "META::Read_Potential"),
                        std::ios::in | std::ios::binary);
    if (!input.is_open())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, "META::Read_Potential",
            "Unable to open metadynamics potential input '%s'.",
            read_potential_file_name.c_str());
    }

    std::size_t line_number = 0;
    std::string line;
    if (!std::getline(input, line))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Read_Potential",
            "Metadynamics potential input is empty.");
    }
    ++line_number;
    const std::vector<std::string> header_tokens = Tokenize_Meta_Line(line);
    const bool versioned_format =
        !header_tokens.empty() &&
        header_tokens[0] == "SPONGE_META_POTENTIAL_V1";
    if (versioned_format)
    {
        if (header_tokens.size() != 4)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input line 1 has %zu header fields; the V1 format "
                "requires exactly version, dimension, representation, and "
                "payload.",
                header_tokens.size());
        }
        const int file_ndim = Parse_Meta_File_Int(controller, header_tokens[1],
                                                  line_number, "dimension");
        const char* expected_representation = use_scatter ? "scatter" : "grid";
        const char* expected_payload = subhill ? "subhill" : "d_force";
        if (file_ndim != ndim || header_tokens[2] != expected_representation ||
            header_tokens[3] != expected_payload)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input header declares %dD %s %s data, but the "
                "active configuration requires %dD %s %s data.",
                file_ndim, header_tokens[2].c_str(), header_tokens[3].c_str(),
                ndim, expected_representation, expected_payload);
        }
    }
    else
    {
        const std::string expected_dimension = std::to_string(ndim) + "D-Meta";
        const bool declares_subhill =
            std::find(header_tokens.begin(), header_tokens.end(), "subhill") !=
            header_tokens.end();
        const bool declares_force =
            std::find(header_tokens.begin(), header_tokens.end(), "d_force") !=
            header_tokens.end();
        if (header_tokens.empty() || header_tokens[0] != expected_dimension ||
            declares_subhill == declares_force || declares_subhill != subhill)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Legacy potential input header must declare dimension '%s' "
                "and exactly the configured '%s' payload.",
                expected_dimension.c_str(), subhill ? "subhill" : "d_force");
        }
    }

    std::vector<float> parsed_mins(ndim);
    std::vector<float> parsed_maxs(ndim);
    std::vector<float> parsed_deltas(ndim);
    std::vector<int> parsed_grids(ndim);
    for (int d = 0; d < ndim; ++d)
    {
        if (!Read_Nonempty_Meta_Line(input, &line, &line_number))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input ended before grid metadata for dimension %d.",
                d);
        }
        const std::vector<std::string> tokens = Tokenize_Meta_Line(line);
        if (tokens.size() != 3)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input line %zu must contain exactly minimum, "
                "maximum, and spacing for dimension %d; got %zu tokens.",
                line_number, d, tokens.size());
        }
        parsed_mins[d] = Parse_Meta_File_Float(controller, tokens[0],
                                               line_number, "CV minimum");
        parsed_maxs[d] = Parse_Meta_File_Float(controller, tokens[1],
                                               line_number, "CV maximum");
        parsed_deltas[d] = Parse_Meta_File_Float(controller, tokens[2],
                                                 line_number, "CV spacing");
        if (!(parsed_maxs[d] > parsed_mins[d]) ||
            !Float_Memory_Is_Normal(&parsed_deltas[d]) ||
            !(parsed_deltas[d] > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input dimension %d requires maximum > minimum and "
                "a positive normal spacing.",
                d);
        }
    }

    if (!Read_Nonempty_Meta_Line(input, &line, &line_number))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Read_Potential",
            "Potential input ended before its grid shape and record count.");
    }
    const std::vector<std::string> shape_tokens = Tokenize_Meta_Line(line);
    if (shape_tokens.size() != static_cast<std::size_t>(ndim) + 1)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Read_Potential",
            "Potential input line %zu must contain %d grid extents followed "
            "by one record count; got %zu tokens.",
            line_number, ndim, shape_tokens.size());
    }
    std::size_t checked_grid_size = 1;
    for (int d = 0; d < ndim; ++d)
    {
        parsed_grids[d] = Parse_Meta_File_Int(controller, shape_tokens[d],
                                              line_number, "grid extent");
        if (parsed_grids[d] <= 1 ||
            checked_grid_size > static_cast<std::size_t>(INT_MAX) /
                                    static_cast<std::size_t>(parsed_grids[d]))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input grid extent %d is invalid or overflows the "
                "flat grid size at dimension %d.",
                parsed_grids[d], d);
        }
        checked_grid_size *= static_cast<std::size_t>(parsed_grids[d]);

        if (versioned_format)
        {
            const double derived_delta =
                (static_cast<double>(parsed_maxs[d]) - parsed_mins[d]) /
                parsed_grids[d];
            const float stored_derived_delta =
                static_cast<float>(derived_delta);
            if (!Double_Memory_Is_Finite(&derived_delta) ||
                !(derived_delta > 0.0) ||
                !Float_Memory_Is_Normal(&stored_derived_delta) ||
                parsed_deltas[d] != stored_derived_delta)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "META::Read_Potential",
                    "Potential input dimension %d declares spacing %.9g, but "
                    "its bounds and %d-point extent require %.9g.",
                    d, static_cast<double>(parsed_deltas[d]), parsed_grids[d],
                    static_cast<double>(stored_derived_delta));
            }
        }
    }
    const int record_count = Parse_Meta_File_Int(controller, shape_tokens[ndim],
                                                 line_number, "record count");
    if (record_count <= 0 ||
        (usegrid &&
         static_cast<std::size_t>(record_count) != checked_grid_size) ||
        (use_scatter &&
         static_cast<std::size_t>(record_count) > checked_grid_size))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Read_Potential",
            "Potential input record count %d is incompatible with its "
            "%zu-point "
            "%s representation.",
            record_count, checked_grid_size, usegrid ? "grid" : "scatter");
    }
    if (checked_grid_size >
        static_cast<std::size_t>(INT_MAX) / static_cast<std::size_t>(ndim))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Read_Potential",
            "Potential input has %zu grid points and %d dimensions; its "
            "flattened force storage exceeds INT_MAX.",
            checked_grid_size, ndim);
    }

    std::vector<std::vector<float>> parsed_coordinates(
        ndim, std::vector<float>(record_count));
    std::vector<float> potential_from_file;
    std::vector<Gdata> force_from_file;
    potential_from_file.reserve(record_count);
    force_from_file.reserve(record_count);
    const std::size_t compact_subhill_tokens =
        static_cast<std::size_t>(ndim) + 3;
    const std::size_t derivative_record_tokens =
        2 * static_cast<std::size_t>(ndim) + 2;
    bool legacy_scatter_matches_grid =
        !versioned_format && use_scatter &&
        static_cast<std::size_t>(record_count) == checked_grid_size;
    for (int record = 0; record < record_count; ++record)
    {
        if (!Read_Nonempty_Meta_Line(input, &line, &line_number))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input ended after %d of %d declared records.",
                record, record_count);
        }
        const std::vector<std::string> tokens = Tokenize_Meta_Line(line);
        const bool valid_compact_subhill =
            subhill && tokens.size() == compact_subhill_tokens;
        const bool valid_legacy_grid_subhill =
            subhill && !versioned_format && usegrid &&
            tokens.size() == derivative_record_tokens;
        const bool valid_derivative_record =
            !subhill && tokens.size() == derivative_record_tokens;
        if (!(valid_compact_subhill || valid_legacy_grid_subhill ||
              valid_derivative_record))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input record %d on line %zu has %zu tokens; the "
                "configured %s payload requires %zu%s.",
                record, line_number, tokens.size(),
                subhill ? "subhill" : "d_force",
                subhill ? compact_subhill_tokens : derivative_record_tokens,
                valid_legacy_grid_subhill ? "" : " tokens");
        }
        Gdata force(ndim, 0.0f);
        int flat_coordinate_index = record;
        for (int d = 0; d < ndim; ++d)
        {
            parsed_coordinates[d][record] = Parse_Meta_File_Float(
                controller, tokens[d], line_number, "CV coordinate");
            if (usegrid || legacy_scatter_matches_grid)
            {
                const int coordinate_index =
                    flat_coordinate_index % parsed_grids[d];
                flat_coordinate_index /= parsed_grids[d];
                // Legacy writers serialized the metadata spacing with %f,
                // but generated coordinates from the full-precision spacing
                // implied by the bounds and extent before formatting each
                // coordinate with the default six significant digits.  Using
                // the rounded metadata spacing here rejects files produced by
                // SPONGE itself (for example a three-point [0, 1] grid).
                const double coordinate_spacing =
                    versioned_format ? static_cast<double>(parsed_deltas[d])
                                     : (static_cast<double>(parsed_maxs[d]) -
                                        parsed_mins[d]) /
                                           parsed_grids[d];
                const float expected_coordinate = static_cast<float>(
                    static_cast<double>(parsed_mins[d]) +
                    (static_cast<double>(coordinate_index) + 0.5) *
                        coordinate_spacing);
                const bool coordinate_matches = Meta_Coordinate_Token_Matches(
                    tokens[d], parsed_coordinates[d][record],
                    expected_coordinate,
                    versioned_format ? MetaCoordinateFormat::round_trip
                                     : MetaCoordinateFormat::legacy_potential);
                if (usegrid && !coordinate_matches)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Read_Potential",
                        "Potential input record %d line %zu has coordinate "
                        "%.9g in dimension %d; grid order requires %.9g.",
                        record, line_number,
                        static_cast<double>(parsed_coordinates[d][record]), d,
                        static_cast<double>(expected_coordinate));
                }
                if (use_scatter && !coordinate_matches)
                {
                    legacy_scatter_matches_grid = false;
                }
            }
        }
        // Validate every serialized diagnostic even when only the raw
        // potential and force are restored.
        for (std::size_t token = ndim; token < tokens.size(); ++token)
        {
            Parse_Meta_File_Float(controller, tokens[token], line_number,
                                  "potential/force value");
        }
        potential_from_file.push_back(Parse_Meta_File_Float(
            controller, tokens.back(), line_number, "raw potential"));
        if (!subhill)
        {
            for (int d = 0; d < ndim; ++d)
            {
                force[d] = Parse_Meta_File_Float(
                    controller,
                    tokens[static_cast<std::size_t>(ndim) + 1 +
                           static_cast<std::size_t>(d)],
                    line_number, "potential derivative");
            }
        }
        force_from_file.push_back(force);
    }
    if (legacy_scatter_matches_grid)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Read_Potential",
            "Legacy potential input is structurally valid as both grid and "
            "scatter data; use the V1 header to declare its representation.");
    }
    while (std::getline(input, line))
    {
        ++line_number;
        if (!Tokenize_Meta_Line(line).empty())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Potential",
                "Potential input contains undeclared data on line %zu.",
                line_number);
        }
    }
    input.close();

    Malloc_Safely((void**)&cv_mins, sizeof(float) * ndim);
    Malloc_Safely((void**)&cv_maxs, sizeof(float) * ndim);
    Malloc_Safely((void**)&n_grids, sizeof(int) * ndim);
    memcpy(cv_mins, parsed_mins.data(), sizeof(float) * ndim);
    memcpy(cv_maxs, parsed_maxs.data(), sizeof(float) * ndim);
    memcpy(n_grids, parsed_grids.data(), sizeof(int) * ndim);
    cv_deltas = parsed_deltas;
    scatter_size = record_count;
    if (use_scatter)
    {
        tcoor.clear();
        tcoor.reserve(ndim);
        for (int d = 0; d < ndim; ++d)
        {
            float* coordinates = NULL;
            Malloc_Safely((void**)&coordinates, sizeof(float) * record_count);
            memcpy(coordinates, parsed_coordinates[d].data(),
                   sizeof(float) * record_count);
            tcoor.push_back(coordinates);
        }
    }

    Set_Grid(controller);
    auto max_it = std::max_element(potential_from_file.begin(),
                                   potential_from_file.end());
    potential_max = *max_it;
    if (usegrid)
    {
        mgrid->potential = potential_from_file;  // potential
        // calculate derivative force dpotential
        if (!subhill)
        {
            for (int idx = 0; idx < mgrid->total_size; ++idx)
            {
                for (int d = 0; d < ndim; ++d)
                {
                    mgrid->force[idx * ndim + d] = force_from_file[idx][d];
                }
            }
        }
    }
    else if (use_scatter)
    {
        mscatter->potential = potential_from_file;
        if (convmeta)
        {
            max_index = std::distance(potential_from_file.begin(), max_it);
        }
        if (!subhill)
        {
            mscatter->force.resize(static_cast<std::size_t>(scatter_size) *
                                   ndim);
            for (int idx = 0; idx < scatter_size; ++idx)
            {
                for (int d = 0; d < ndim; ++d)
                {
                    mscatter->force[idx * ndim + d] = force_from_file[idx][d];
                }
            }
        }
        if (mask)
        {
            for (int index = 0; index < mscatter->num_points; ++index)
            {
                const Axis& coor = mscatter->Get_Coordinate(index);
                int gidx = mgrid->Get_Flat_Index(coor);
                mgrid->potential[gidx] = potential_from_file[index];

                for (int d = 0; d < ndim; ++d)
                {
                    mgrid->force[gidx * ndim + d] = force_from_file[index][d];
                }
            }
        }
    }
    if (mgrid != NULL) mgrid->Sync_To_Device();
    if (mscatter != NULL) mscatter->Sync_To_Device();
}

void META::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized)
    {
        return;
    }
    if (CONTROLLER::MPI_size == 1)
    {
        controller->Step_Print(this->module_name, potential_local);
        controller->Step_Print("rbias", rbias);
        controller->Step_Print("rct", rct);
        return;
    }
#ifdef USE_MPI
    if (CONTROLLER::MPI_rank == CONTROLLER::MPI_size - 1)
    {
        MPI_Send(&potential_local, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&rbias, 1, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&rct, 1, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
    }
    if (CONTROLLER::MPI_rank == 0)
    {
        MPI_Recv(&potential_local, 1, MPI_FLOAT, CONTROLLER::MPI_size - 1, 0,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rbias, 1, MPI_FLOAT, CONTROLLER::MPI_size - 1, 1,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&rct, 1, MPI_FLOAT, CONTROLLER::MPI_size - 1, 2,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        controller->Step_Print(this->module_name, potential_local);
        controller->Step_Print("rbias", rbias);
        controller->Step_Print("rct", rct);
    }
#endif
}

using TimePoint = std::chrono::time_point<std::chrono::steady_clock>;

std::string GetTime(TimePoint& local_time)
{
    local_time = std::chrono::steady_clock::now();
    return std::string();
}

std::string GetDuration(const TimePoint& late_time, const TimePoint& early_time,
                        float& duration)
{
    // Some constants.
    const auto elapsed = late_time - early_time;
    const double elapsed_seconds =
        std::chrono::duration<double>(elapsed).count();
    const std::uint64_t elapsed_microseconds =
        static_cast<std::uint64_t>(std::max(0.0, elapsed_seconds) * 1000000.0);
    const std::uint64_t seconds = elapsed_microseconds / 1000000ULL;
    const std::uint64_t microseconds = elapsed_microseconds % 1000000ULL;
    const size_t Second2Day = 86400L;
    const size_t Second2Hour = 3600L;
    const size_t Second2Minute = 60L;
    const std::uint64_t day = seconds / Second2Day;
    const std::uint64_t hour = seconds % Second2Day / Second2Hour;
    const std::uint64_t min =
        seconds % Second2Day % Second2Hour / Second2Minute;
    const std::uint64_t second =
        seconds % Second2Day % Second2Hour % Second2Minute;
    // Calculate duration in second.
    const int BufferSize = 2048;
    char buffer[BufferSize];
    snprintf(buffer, BufferSize,
             "%lu days %lu hours %lu minutes %lu seconds %.1f milliseconds",
             static_cast<unsigned long>(day), static_cast<unsigned long>(hour),
             static_cast<unsigned long>(min),
             static_cast<unsigned long>(second), microseconds * 0.001);
    duration = static_cast<float>(std::max(0.0, elapsed_seconds));
    return std::string(buffer);
}

static __global__ void Update_Edge_Effect_Grid(
    int total_size, int ndim, int scatter_size, const int* num_points,
    const float* lower, const float* spacing, const float* sigmas,
    const float* periods, const float* scatter_coordinates, int calculate_force,
    float* normal_lse, float* normal_force, int* error_flag)
{
    SIMPLE_DEVICE_FOR(gidx, total_size)
    {
        if (scatter_size <= 0)
        {
            atomicExch(error_flag, 1);
            continue;
        }
        double max_log = 0.0;
        double sum_exp = 0.0;
        bool has_log = false;
        for (int index = 0; index < scatter_size; ++index)
        {
            double pregauss = 0.0;
            int flat_index = gidx;
            for (int d = 0; d < ndim; ++d)
            {
                const int coordinate_index = flat_index % num_points[d];
                flat_index /= num_points[d];
                const double value =
                    static_cast<double>(lower[d]) +
                    (static_cast<double>(coordinate_index) + 0.5) * spacing[d];
                double diff = value - scatter_coordinates[index * ndim + d];
                if (periods[d] != 0.0f)
                {
                    diff -= round(diff / periods[d]) * periods[d];
                }
                const double scaled = diff * sigmas[d];
                pregauss -= 0.5 * scaled * scaled;
            }
            if (!isfinite(pregauss))
            {
                atomicExch(error_flag, 1);
                has_log = false;
                break;
            }
            if (!has_log)
            {
                max_log = pregauss;
                sum_exp = 1.0;
                has_log = true;
            }
            else if (pregauss > max_log)
            {
                sum_exp = sum_exp * exp(max_log - pregauss) + 1.0;
                max_log = pregauss;
            }
            else
            {
                sum_exp += exp(pregauss - max_log);
            }
        }
        if (!has_log) continue;
        const double lse = max_log + log(sum_exp);
        const float stored_lse = static_cast<float>(lse);
        if (!isfinite(lse) || fabs(lse) > FLT_MAX ||
            (stored_lse != 0.0f && fabsf(stored_lse) < FLT_MIN))
        {
            atomicExch(error_flag, 1);
            continue;
        }
        normal_lse[gidx] = stored_lse;

        if (calculate_force)
        {
            for (int force_dim = 0; force_dim < ndim; ++force_dim)
            {
                double log_derivative = 0.0;
                for (int index = 0; index < scatter_size; ++index)
                {
                    double pregauss = 0.0;
                    double force_diff = 0.0;
                    int flat_index = gidx;
                    for (int d = 0; d < ndim; ++d)
                    {
                        const int coordinate_index = flat_index % num_points[d];
                        flat_index /= num_points[d];
                        const double value =
                            static_cast<double>(lower[d]) +
                            (static_cast<double>(coordinate_index) + 0.5) *
                                spacing[d];
                        double diff =
                            value - scatter_coordinates[index * ndim + d];
                        if (periods[d] != 0.0f)
                        {
                            diff -= round(diff / periods[d]) * periods[d];
                        }
                        if (d == force_dim) force_diff = diff;
                        const double scaled = diff * sigmas[d];
                        pregauss -= 0.5 * scaled * scaled;
                    }
                    const double derivative =
                        -force_diff * sigmas[force_dim] * sigmas[force_dim];
                    log_derivative += derivative * exp(pregauss - lse);
                }
                const float stored_derivative =
                    static_cast<float>(log_derivative);
                if (!isfinite(log_derivative) ||
                    fabs(log_derivative) > FLT_MAX ||
                    (stored_derivative != 0.0f &&
                     fabsf(stored_derivative) < FLT_MIN) ||
                    (log_derivative != 0.0 && stored_derivative == 0.0f))
                {
                    atomicExch(error_flag, 1);
                    continue;
                }
                normal_force[gidx * ndim + force_dim] = stored_derivative;
            }
        }
    }
}

static __global__ void Update_Grid_With_Hill(
    int total_size, int ndim, const int* num_points, const float* lower,
    const float* spacing, const float* hill_centers, const float* hill_inv_w,
    const float* hill_periods, float factor, int update_force, float* potential,
    float* force)
{
    SIMPLE_DEVICE_FOR(idx, total_size)
    {
        int flat = idx;
        float pot = 1.0f;
        for (int d = 0; d < ndim; ++d)
        {
            int i = flat % num_points[d];
            flat /= num_points[d];
            const float coord =
                static_cast<float>(static_cast<double>(lower[d]) +
                                   (static_cast<double>(i) + 0.5) * spacing[d]);
            float diff = coord - hill_centers[d];
            if (hill_periods[d] > 0.0f)
            {
                diff -= roundf(diff / hill_periods[d]) * hill_periods[d];
            }
            const float x = diff * hill_inv_w[d];
            pot *= expf(-0.5f * x * x);
        }
        potential[idx] += factor * pot;

        if (update_force)
        {
            flat = idx;
            for (int d = 0; d < ndim; ++d)
            {
                const int i = flat % num_points[d];
                flat /= num_points[d];
                const float coord = static_cast<float>(
                    static_cast<double>(lower[d]) +
                    (static_cast<double>(i) + 0.5) * spacing[d]);
                float diff = coord - hill_centers[d];
                if (hill_periods[d] > 0.0f)
                {
                    diff -= roundf(diff / hill_periods[d]) * hill_periods[d];
                }
                const float x = diff * hill_inv_w[d];
                force[idx * ndim + d] += factor * (-x * hill_inv_w[d] * pot);
            }
        }
    }
}

static __global__ void Update_Scatter_With_Hill(
    int num_points, int ndim, const float* coordinates, const float* periods,
    const float* hill_centers, const float* hill_inv_w, float factor,
    int update_force, int use_cutoff, const float* cutoff, float* potential,
    float* force)
{
    SIMPLE_DEVICE_FOR(index, num_points)
    {
        const float* coord = coordinates + index * ndim;
        bool within = true;
        float hill_potential = 1.0f;
        for (int d = 0; d < ndim; ++d)
        {
            float diff = coord[d] - hill_centers[d];
            if (periods[d] > 0.0f)
            {
                diff -= roundf(diff / periods[d]) * periods[d];
            }
            if (use_cutoff && fabsf(diff) > cutoff[d])
            {
                within = false;
                break;
            }
            const float x = diff * hill_inv_w[d];
            hill_potential *= expf(-0.5f * x * x);
        }
        if (within)
        {
            potential[index] += factor * hill_potential;

            if (update_force && force != NULL)
            {
                float* data = force + index * ndim;
                for (int i = 0; i < ndim; ++i)
                {
                    float diff = coord[i] - hill_centers[i];
                    if (periods[i] > 0.0f)
                    {
                        diff -= roundf(diff / periods[i]) * periods[i];
                    }
                    const float x = diff * hill_inv_w[i];
                    data[i] += factor * (-x * hill_inv_w[i] * hill_potential);
                }
            }
        }
    }
}

static float Store_Checked_Meta_Float(CONTROLLER* controller, int error_number,
                                      const char* error_by,
                                      const char* quantity, double value)
{
    if (!Double_Memory_Is_Finite(&value) ||
        value > std::numeric_limits<float>::max() ||
        value < -std::numeric_limits<float>::max())
    {
        controller->Throw_Formatted_SPONGE_Error(
            error_number, error_by,
            "Metadynamics %s is outside the finite float range (%.17g).",
            quantity, value);
    }
    const float stored = static_cast<float>(value);
    if (!Float_Memory_Is_Zero_Or_Normal(&stored) ||
        (value != 0.0 && stored == 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            error_number, error_by,
            "Metadynamics %s is not representable as a zero-or-normal float "
            "(%.17g).",
            quantity, value);
    }
    return stored;
}

static double Checked_Meta_Log_Add_Exp(CONTROLLER* controller,
                                       const char* error_by, double lhs,
                                       double rhs)
{
    if (!Double_Memory_Is_Finite(&lhs) || !Double_Memory_Is_Finite(&rhs))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Metadynamics log-sum-exp received a non-finite operand.");
    }
    const double maximum = std::max(lhs, rhs);
    const double result = maximum + std::log1p(std::exp(-std::fabs(lhs - rhs)));
    if (!Double_Memory_Is_Finite(&result))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Metadynamics log-sum-exp produced a non-finite result.");
    }
    return result;
}

static double Checked_Meta_Log_Partition(CONTROLLER* controller,
                                         const char* error_by, float factor,
                                         float* value_max,
                                         const std::vector<float>& values)
{
    if (values.empty())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, error_by,
            "Metadynamics cannot evaluate a partition function over an empty "
            "potential representation.");
    }
    if (!Float_Memory_Is_Normal(&factor) || !(factor > 0.0f))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, error_by,
            "Metadynamics partition factor must be a finite positive normal "
            "float.");
    }

    double maximum_scaled = -std::numeric_limits<double>::infinity();
    float maximum_value = -std::numeric_limits<float>::max();
    for (float value : values)
    {
        if (!Float_Memory_Is_Zero_Or_Normal(&value))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, error_by,
                "Metadynamics partition potentials must be finite "
                "zero-or-normal floats.");
        }
        maximum_value = std::max(maximum_value, value);
        maximum_scaled =
            std::max(maximum_scaled, static_cast<double>(factor) * value);
    }
    double sum = 0.0;
    for (float value : values)
    {
        sum += std::exp(static_cast<double>(factor) * value - maximum_scaled);
    }
    const double result = maximum_scaled + std::log(sum);
    if (!Double_Memory_Is_Finite(&sum) || !(sum > 0.0) ||
        !Double_Memory_Is_Finite(&result))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, error_by,
            "Metadynamics partition evaluation produced a non-finite or "
            "non-positive derived value.");
    }
    *value_max = maximum_value;
    return result;
}

static bool Lock_Meta_Append_File(int descriptor, bool exclusive,
                                  int* error_number)
{
#ifdef _WIN32
    const intptr_t native_handle = _get_osfhandle(descriptor);
    if (native_handle == -1)
    {
        *error_number = Meta_Effective_IO_Error(errno);
        return false;
    }
    OVERLAPPED overlap = {};
    if (LockFileEx(reinterpret_cast<HANDLE>(native_handle),
                   exclusive ? LOCKFILE_EXCLUSIVE_LOCK : 0, 0, MAXDWORD,
                   MAXDWORD, &overlap) == 0)
    {
        *error_number = EACCES;
        return false;
    }
    return true;
#else
    struct flock lock = {};
    lock.l_type = exclusive ? F_WRLCK : F_RDLCK;
    lock.l_whence = SEEK_SET;
    for (;;)
    {
        if (fcntl(descriptor, F_SETLKW, &lock) == 0) return true;
        if (errno != EINTR)
        {
            *error_number = Meta_Effective_IO_Error(errno);
            return false;
        }
    }
#endif
}

static void Unlock_Meta_Append_File(int descriptor)
{
#ifdef _WIN32
    const intptr_t native_handle = _get_osfhandle(descriptor);
    if (native_handle != -1)
    {
        OVERLAPPED overlap = {};
        UnlockFileEx(reinterpret_cast<HANDLE>(native_handle), 0, MAXDWORD,
                     MAXDWORD, &overlap);
    }
#else
    struct flock lock = {};
    lock.l_type = F_UNLCK;
    lock.l_whence = SEEK_SET;
    fcntl(descriptor, F_SETLK, &lock);
#endif
}

static std::int64_t Meta_Append_File_End(int descriptor)
{
#ifdef _WIN32
    return static_cast<std::int64_t>(_lseeki64(descriptor, 0, SEEK_END));
#else
    return static_cast<std::int64_t>(lseek(descriptor, 0, SEEK_END));
#endif
}

static bool Read_Meta_Append_At(int descriptor, std::int64_t offset,
                                char* buffer, std::size_t size,
                                int* error_number)
{
#ifdef _WIN32
    if (_lseeki64(descriptor, offset, SEEK_SET) < 0)
    {
        *error_number = Meta_Effective_IO_Error(errno);
        return false;
    }
    std::size_t completed = 0;
    while (completed < size)
    {
        const unsigned request = static_cast<unsigned>(std::min<std::size_t>(
            size - completed, static_cast<std::size_t>(INT_MAX)));
        const int read_size = _read(descriptor, buffer + completed, request);
        if (read_size > 0)
        {
            completed += static_cast<std::size_t>(read_size);
            continue;
        }
        if (read_size < 0 && errno == EINTR) continue;
        *error_number = read_size < 0 ? Meta_Effective_IO_Error(errno) : EIO;
        return false;
    }
#else
    std::size_t completed = 0;
    while (completed < size)
    {
        const ssize_t read_size =
            pread(descriptor, buffer + completed, size - completed,
                  static_cast<off_t>(offset + completed));
        if (read_size > 0)
        {
            completed += static_cast<std::size_t>(read_size);
            continue;
        }
        if (read_size < 0 && errno == EINTR) continue;
        *error_number = read_size < 0 ? Meta_Effective_IO_Error(errno) : EIO;
        return false;
    }
#endif
    return true;
}

static bool Find_Complete_Meta_Append_Prefix(int descriptor,
                                             std::int64_t file_size,
                                             std::int64_t* complete_size,
                                             int* error_number)
{
    char buffer[4096];
    std::int64_t cursor = file_size;
    while (cursor > 0)
    {
        const std::size_t chunk_size = static_cast<std::size_t>(
            std::min<std::int64_t>(cursor, sizeof(buffer)));
        const std::int64_t chunk_start =
            cursor - static_cast<std::int64_t>(chunk_size);
        if (!Read_Meta_Append_At(descriptor, chunk_start, buffer, chunk_size,
                                 error_number))
        {
            return false;
        }
        for (std::size_t index = chunk_size; index != 0; --index)
        {
            if (buffer[index - 1] == '\n')
            {
                *complete_size = chunk_start + static_cast<std::int64_t>(index);
                return true;
            }
        }
        cursor = chunk_start;
    }
    *complete_size = 0;
    return true;
}

enum class MetaLineReadStatus
{
    line,
    end,
    too_long,
    io_error,
    allocation_error,
};

class Meta_Bounded_Line_Reader
{
   public:
    Meta_Bounded_Line_Reader(int descriptor, std::int64_t committed_size,
                             std::size_t maximum_line_size)
        : descriptor_(descriptor),
          committed_size_(committed_size),
          maximum_line_size_(maximum_line_size)
    {
    }

    MetaLineReadStatus Next(std::string* line, int* error_number)
    {
        line->clear();
        for (;;)
        {
            if (buffer_cursor_ == buffer_size_)
            {
                if (next_read_offset_ == committed_size_)
                {
                    return line->empty() ? MetaLineReadStatus::end
                                         : MetaLineReadStatus::io_error;
                }
                const std::size_t request =
                    static_cast<std::size_t>(std::min<std::int64_t>(
                        committed_size_ - next_read_offset_, sizeof(buffer_)));
                if (!Read_Meta_Append_At(descriptor_, next_read_offset_,
                                         buffer_, request, error_number))
                {
                    return MetaLineReadStatus::io_error;
                }
                next_read_offset_ += static_cast<std::int64_t>(request);
                buffer_cursor_ = 0;
                buffer_size_ = request;
            }

            const void* newline_address = memchr(buffer_ + buffer_cursor_, '\n',
                                                 buffer_size_ - buffer_cursor_);
            const std::size_t segment_size =
                newline_address == NULL
                    ? buffer_size_ - buffer_cursor_
                    : static_cast<const char*>(newline_address) -
                          (buffer_ + buffer_cursor_);
            if (line->size() > maximum_line_size_ ||
                segment_size > maximum_line_size_ - line->size())
            {
                return MetaLineReadStatus::too_long;
            }
            try
            {
                line->append(buffer_ + buffer_cursor_, segment_size);
            }
            catch (const std::length_error&)
            {
                return MetaLineReadStatus::too_long;
            }
            catch (const std::bad_alloc&)
            {
                *error_number = ENOMEM;
                return MetaLineReadStatus::allocation_error;
            }
            buffer_cursor_ += segment_size;
            if (newline_address != NULL)
            {
                ++buffer_cursor_;
                return MetaLineReadStatus::line;
            }
        }
    }

   private:
    int descriptor_ = -1;
    std::int64_t committed_size_ = 0;
    std::int64_t next_read_offset_ = 0;
    std::size_t maximum_line_size_ = 0;
    char buffer_[8192] = {};
    std::size_t buffer_cursor_ = 0;
    std::size_t buffer_size_ = 0;
};

static bool Resize_Meta_Append_File(int descriptor, std::int64_t size,
                                    int* error_number)
{
#ifdef _WIN32
    const errno_t status = _chsize_s(descriptor, size);
    if (status != 0)
    {
        *error_number = static_cast<int>(status);
        return false;
    }
#else
    if (ftruncate(descriptor, static_cast<off_t>(size)) != 0)
    {
        *error_number = Meta_Effective_IO_Error(errno);
        return false;
    }
#endif
    return true;
}

static bool Sync_Meta_Append_File(int descriptor, int* error_number)
{
#ifdef _WIN32
    if (_commit(descriptor) != 0)
#else
    if (fsync(descriptor) != 0)
#endif
    {
        *error_number = Meta_Effective_IO_Error(errno);
        return false;
    }
    return true;
}

static bool Write_Meta_Append_Part(int descriptor, const char* data,
                                   std::size_t size, int* error_number)
{
    errno = 0;
#ifdef _WIN32
    const int write_size =
        _write(descriptor, data, static_cast<unsigned>(size));
#else
    const ssize_t write_size = write(descriptor, data, size);
#endif
    if (write_size >= 0 && static_cast<std::size_t>(write_size) == size)
    {
        return true;
    }
    *error_number = write_size < 0 ? Meta_Effective_IO_Error(errno) : EIO;
    return false;
}

static void Append_Meta_Record(CONTROLLER* controller, const std::string& name,
                               const std::string& record, const char* error_by)
{
    if (record.size() < 2 || record.back() != '\n' ||
        record.find('\n') != record.size() - 1)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, error_by,
            "A metadynamics record for '%s' must contain one nonempty line "
            "terminated by exactly one newline commit marker.",
            name.c_str());
    }
#ifdef _WIN32
    if (record.size() > static_cast<std::size_t>(INT_MAX))
#else
    if (record.size() > static_cast<std::size_t>(SSIZE_MAX))
#endif
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, error_by,
            "A single metadynamics record for '%s' has %zu bytes and exceeds "
            "the platform's atomic append limit.",
            name.c_str(), record.size());
    }

    // POSIX record locks are owned by the process, not by a calling thread,
    // and closing any descriptor for the inode releases that process's locks.
    // Keep every local hills descriptor lifecycle disjoint; the kernel lock
    // below then extends the same critical section across processes.
    const std::lock_guard<std::mutex> process_lock(
        sponge_meta_io::Process_Mutex());
    const fs::path native_name = Meta_Native_Path(controller, name, error_by);

    int descriptor = -1;
#ifdef _WIN32
    const errno_t open_status =
        _wsopen_s(&descriptor, native_name.c_str(),
                  _O_RDWR | _O_CREAT | _O_APPEND | _O_BINARY, _SH_DENYNO,
                  _S_IREAD | _S_IWRITE);
    if (open_status != 0)
#else
    fs::path parent_path;
    std::string basename;
    try
    {
        parent_path = native_name.parent_path();
        if (parent_path.empty()) parent_path = fs::path(".");
        const fs::path filename = native_name.filename();
        if (filename.empty() || filename == fs::path(".") ||
            filename == fs::path(".."))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, error_by,
                "Metadynamics append path '%s' does not name a file.",
                name.c_str());
        }
        basename = filename.native();
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorMallocFailed, error_by,
            "Unable to allocate append-path storage for '%s'.", name.c_str());
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to prepare metadynamics append path '%s': %s.",
            name.c_str(), error.what());
    }
    errno = 0;
    const int parent_descriptor =
        open(parent_path.native().c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
    if (parent_descriptor < 0)
    {
        const int parent_error = Meta_Effective_IO_Error(errno);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to open the parent directory of metadynamics append "
            "output '%s': %s.",
            name.c_str(), strerror(parent_error));
    }
    struct stat existing_path_status = {};
    errno = 0;
    if (fstatat(parent_descriptor, basename.c_str(), &existing_path_status,
                AT_SYMLINK_NOFOLLOW) == 0)
    {
        if (!S_ISREG(existing_path_status.st_mode))
        {
            close(parent_descriptor);
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by,
                "Metadynamics append output '%s' is not a stable regular "
                "file: %s.",
                name.c_str(), strerror(EINVAL));
        }
    }
    else if (errno != ENOENT)
    {
        const int inspect_error = Meta_Effective_IO_Error(errno);
        close(parent_descriptor);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to inspect metadynamics append output '%s' before "
            "opening it: %s.",
            name.c_str(), strerror(inspect_error));
    }
    const int append_open_flags =
        O_RDWR | O_APPEND | O_CLOEXEC | O_NOFOLLOW | O_NONBLOCK;
    int append_open_error = 0;
    // Elect exactly one creator before opening the shared inode.  A single
    // nonexclusive O_CREAT open can transiently return ENOENT when several
    // processes create the same basename on Darwin/APFS; O_EXCL gives the
    // losers the unambiguous existing-file path instead.
    for (;;)
    {
        errno = 0;
        descriptor = openat(parent_descriptor, basename.c_str(),
                            append_open_flags | O_CREAT | O_EXCL, 0666);
        if (descriptor >= 0) break;
        append_open_error = Meta_Effective_IO_Error(errno);
        if (append_open_error != EINTR) break;
    }
    if (descriptor < 0 && append_open_error == EEXIST)
    {
        for (;;)
        {
            errno = 0;
            descriptor =
                openat(parent_descriptor, basename.c_str(), append_open_flags);
            if (descriptor >= 0) break;
            append_open_error = Meta_Effective_IO_Error(errno);
            if (append_open_error != EINTR) break;
        }
    }
    if (descriptor < 0)
#endif
    {
#ifdef _WIN32
        const int open_error = static_cast<int>(open_status);
#else
        const int open_error = append_open_error;
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to open metadynamics append output '%s': %s.", name.c_str(),
            strerror(Meta_Effective_IO_Error(open_error)));
    }

    int append_identity_error = 0;
#ifdef _WIN32
    struct _stat64 append_status = {};
    errno = 0;
    if (_fstat64(descriptor, &append_status) != 0)
    {
        append_identity_error = Meta_Effective_IO_Error(errno);
    }
    else if ((append_status.st_mode & _S_IFMT) != _S_IFREG)
    {
        append_identity_error = EINVAL;
    }
#else
    struct stat append_descriptor_status = {};
    struct stat append_path_status = {};
    errno = 0;
    if (fstat(descriptor, &append_descriptor_status) != 0)
    {
        append_identity_error = Meta_Effective_IO_Error(errno);
    }
    else if (!S_ISREG(append_descriptor_status.st_mode))
    {
        append_identity_error = EINVAL;
    }
    else
    {
        errno = 0;
        if (fstatat(parent_descriptor, basename.c_str(), &append_path_status,
                    AT_SYMLINK_NOFOLLOW) != 0)
        {
            append_identity_error = Meta_Effective_IO_Error(errno);
        }
        else if (!S_ISREG(append_path_status.st_mode) ||
                 append_path_status.st_dev != append_descriptor_status.st_dev ||
                 append_path_status.st_ino != append_descriptor_status.st_ino)
        {
            append_identity_error = ESTALE;
        }
    }
#endif
    if (append_identity_error != 0)
    {
#ifdef _WIN32
        _close(descriptor);
#else
        close(descriptor);
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Metadynamics append output '%s' is not a stable regular file: "
            "%s.",
            name.c_str(), strerror(append_identity_error));
    }

    int io_error = 0;
    if (!Lock_Meta_Append_File(descriptor, true, &io_error))
    {
#ifdef _WIN32
        _close(descriptor);
#else
        close(descriptor);
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to lock metadynamics append output '%s': %s.", name.c_str(),
            strerror(io_error));
    }

    std::int64_t committed_size = Meta_Append_File_End(descriptor);
    if (committed_size < 0)
    {
        io_error = Meta_Effective_IO_Error(errno);
        Unlock_Meta_Append_File(descriptor);
#ifdef _WIN32
        _close(descriptor);
#else
        close(descriptor);
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "Unable to inspect metadynamics append output '%s': %s.",
            name.c_str(), strerror(io_error));
    }

    if (committed_size > 0)
    {
        char final_byte = 0;
        if (!Read_Meta_Append_At(descriptor, committed_size - 1, &final_byte, 1,
                                 &io_error))
        {
            Unlock_Meta_Append_File(descriptor);
#ifdef _WIN32
            _close(descriptor);
#else
            close(descriptor);
            close(parent_descriptor);
#endif
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOpenFileFailed, error_by,
                "Unable to inspect the final metadynamics record in '%s': "
                "%s.",
                name.c_str(), strerror(io_error));
        }
        if (final_byte != '\n')
        {
            std::int64_t complete_size = 0;
            if (!Find_Complete_Meta_Append_Prefix(descriptor, committed_size,
                                                  &complete_size, &io_error) ||
                !Resize_Meta_Append_File(descriptor, complete_size,
                                         &io_error) ||
                !Sync_Meta_Append_File(descriptor, &io_error))
            {
                Unlock_Meta_Append_File(descriptor);
#ifdef _WIN32
                _close(descriptor);
#else
                close(descriptor);
                close(parent_descriptor);
#endif
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOpenFileFailed, error_by,
                    "Unable to remove an incomplete final metadynamics "
                    "record from '%s': %s.",
                    name.c_str(), strerror(io_error));
            }
            committed_size = complete_size;
        }
    }

    const std::size_t payload_size = record.size() - 1;
    int write_error = 0;
    if (!Write_Meta_Append_Part(descriptor, record.data(), payload_size,
                                &write_error))
    {
        int rollback_error = 0;
        const bool rollback_ok =
            Resize_Meta_Append_File(descriptor, committed_size,
                                    &rollback_error) &&
            Sync_Meta_Append_File(descriptor, &rollback_error);
        Unlock_Meta_Append_File(descriptor);
#ifdef _WIN32
        _close(descriptor);
#else
        close(descriptor);
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            rollback_ok
                ? "Failed to append a metadynamics record payload to '%s': %s; "
                  "the previous complete prefix was restored."
                : "Failed to append a metadynamics record payload to '%s': %s; "
                  "rollback also failed (error %d).  The unterminated tail has "
                  "no "
                  "commit marker and will be ignored on restart.",
            name.c_str(), strerror(write_error), rollback_error);
    }

    if (!Sync_Meta_Append_File(descriptor, &io_error))
    {
        const int sync_error = io_error;
        int rollback_error = 0;
        const bool rollback_ok =
            Resize_Meta_Append_File(descriptor, committed_size,
                                    &rollback_error) &&
            Sync_Meta_Append_File(descriptor, &rollback_error);
        Unlock_Meta_Append_File(descriptor);
#ifdef _WIN32
        _close(descriptor);
#else
        close(descriptor);
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            rollback_ok
                ? "Failed to synchronize a metadynamics record payload in "
                  "'%s': %s; the previous complete prefix was restored."
                : "Failed to synchronize a metadynamics record payload in "
                  "'%s': %s; rollback also failed (error %d).  The "
                  "unterminated "
                  "tail has no commit marker and will be ignored on restart.",
            name.c_str(), strerror(sync_error), rollback_error);
    }

    write_error = 0;
    if (!Write_Meta_Append_Part(descriptor, "\n", 1, &write_error))
    {
        int rollback_error = 0;
        const bool rollback_ok =
            Resize_Meta_Append_File(descriptor, committed_size,
                                    &rollback_error) &&
            Sync_Meta_Append_File(descriptor, &rollback_error);
        Unlock_Meta_Append_File(descriptor);
#ifdef _WIN32
        _close(descriptor);
#else
        close(descriptor);
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            rollback_ok
                ? "Failed to append the metadynamics record commit marker to "
                  "'%s': %s; the previous complete prefix was restored."
                : "Failed to append the metadynamics record commit marker to "
                  "'%s': %s; rollback also failed (error %d).  No complete new "
                  "record was observed.",
            name.c_str(), strerror(write_error), rollback_error);
    }

    if (!Sync_Meta_Append_File(descriptor, &io_error))
    {
        const int sync_error = io_error;
        int rollback_error = 0;
        const bool rollback_ok =
            Resize_Meta_Append_File(descriptor, committed_size,
                                    &rollback_error) &&
            Sync_Meta_Append_File(descriptor, &rollback_error);
        Unlock_Meta_Append_File(descriptor);
#ifdef _WIN32
        _close(descriptor);
#else
        close(descriptor);
        close(parent_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            rollback_ok
                ? "Failed to synchronize the metadynamics record commit marker "
                  "in '%s': %s; the previous complete prefix was restored."
                : "Failed to synchronize the metadynamics record commit marker "
                  "in '%s': %s; rollback also failed (error %d).  The record's "
                  "durable commit status is indeterminate and requires manual "
                  "inspection before restart.",
            name.c_str(), strerror(sync_error), rollback_error);
    }

#ifndef _WIN32
    errno = 0;
    if (fsync(parent_descriptor) != 0)
    {
        const int directory_error = Meta_Effective_IO_Error(errno);
        Unlock_Meta_Append_File(descriptor);
        close(descriptor);
        close(parent_descriptor);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "The metadynamics record in '%s' was synchronized, but its parent "
            "directory could not be synchronized: %s.  The record is visible, "
            "but crash durability of a newly created directory entry is not "
            "guaranteed.",
            name.c_str(), strerror(directory_error));
    }
#endif

    Unlock_Meta_Append_File(descriptor);
    errno = 0;
#ifdef _WIN32
    const int close_status = _close(descriptor);
#else
    const int close_status = close(descriptor);
#endif
    const int close_error =
        close_status == 0 ? 0 : Meta_Effective_IO_Error(errno);
#ifndef _WIN32
    errno = 0;
    const int parent_close_status = close(parent_descriptor);
    const int parent_close_error =
        parent_close_status == 0 ? 0 : Meta_Effective_IO_Error(errno);
#endif
    if (close_status != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "The complete metadynamics record in '%s' was synchronized, but "
            "closing the append descriptor failed: %s.",
            name.c_str(), strerror(close_error));
    }
#ifndef _WIN32
    if (parent_close_status != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, error_by,
            "The complete metadynamics record in '%s' was synchronized and "
            "its data descriptor was closed, but closing the parent directory "
            "descriptor failed: %s.",
            name.c_str(), strerror(parent_close_error));
    }
#endif
}

static void hilllog(CONTROLLER* controller, const std::string& fn,
                    const std::vector<float>& hillcenter, float hillheight,
                    bool do_negative, float potential_max, int max_index,
                    bool has_scatter, int scatter_index, bool has_mask,
                    float exit_tag)
{
    if (fn.empty()) return;
    std::ostringstream record;
    record.exceptions(std::ios::badbit | std::ios::failbit);
    record << std::setprecision(std::numeric_limits<float>::max_digits10);
    auto format_float = [&](float value, const char* field_name)
    {
        if (!Float_Memory_Is_Zero_Or_Normal(&value))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Add_Potential",
                "Cannot persist a non-finite or subnormal hills %s.",
                field_name);
        }
        record << value << "\t";
    };
    for (float center : hillcenter) format_float(center, "center");
    format_float(hillheight, "height");
    if (do_negative)
    {
        format_float(potential_max, "maximum potential");
        if (max_index < 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Add_Potential",
                "Cannot persist a negative hills projection index.");
        }
        record << max_index << "\t";
    }
    if (has_scatter)
    {
        if (scatter_index < 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Add_Potential",
                "Cannot persist a negative hills scatter index.");
        }
        record << scatter_index << "\t";
    }
    if (has_mask)
    {
        if (exit_tag < 0.0f)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Add_Potential",
                "Cannot persist a negative committed hills exit tag.");
        }
        format_float(exit_tag, "committed exit tag");
    }
    record << "\n";
    Append_Meta_Record(controller, fn, record.str(), "META::Add_Potential");
}

static void Hash_Meta_Edge_U32(std::uint64_t* hash, std::uint32_t value)
{
    for (int shift = 0; shift < 32; shift += 8)
    {
        *hash ^= static_cast<std::uint8_t>(value >> shift);
        *hash *= UINT64_C(1099511628211);
    }
}

static void Hash_Meta_Edge_Float(std::uint64_t* hash, float value)
{
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SinkMeta edge fingerprint assumes 32-bit floats");
    memcpy(&bits, &value, sizeof(bits));
    Hash_Meta_Edge_U32(hash, bits);
}

static std::uint64_t Meta_Edge_Config_Fingerprint(const META& meta,
                                                  bool includes_force)
{
    std::uint64_t hash = UINT64_C(1469598103934665603);
    Hash_Meta_Edge_U32(&hash, static_cast<std::uint32_t>(meta.ndim));
    Hash_Meta_Edge_U32(&hash,
                       static_cast<std::uint32_t>(meta.mgrid->total_size));
    Hash_Meta_Edge_U32(&hash, static_cast<std::uint32_t>(meta.scatter_size));
    Hash_Meta_Edge_U32(&hash, includes_force ? 1U : 0U);
    for (int d = 0; d < meta.ndim; ++d)
    {
        Hash_Meta_Edge_U32(
            &hash, static_cast<std::uint32_t>(meta.mgrid->num_points[d]));
        Hash_Meta_Edge_Float(&hash, meta.mgrid->lower[d]);
        Hash_Meta_Edge_Float(&hash, meta.mgrid->spacing[d]);
        Hash_Meta_Edge_Float(&hash, meta.sigmas[d]);
        Hash_Meta_Edge_Float(&hash, meta.periods[d]);
    }
    for (float coordinate : meta.mscatter->coordinates_flat)
    {
        Hash_Meta_Edge_Float(&hash, coordinate);
    }
    return hash;
}

bool META::Read_Edge_File(const char* file_name, std::vector<float>& potential)
{
    const std::string input_name = file_name == NULL ? "" : file_name;
    std::ifstream input(
        Meta_Native_Path(controller, input_name, "META::Read_Edge_File"),
        std::ios::in | std::ios::binary);
    if (!input.is_open()) return false;
    if (mgrid == NULL || mgrid->total_size <= 0 || ndim <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Read_Edge_File",
            "Edge cache was requested before a valid metadynamics grid was "
            "constructed.");
    }

    const int total = mgrid->total_size;
    const bool requires_force = do_negative || mask != 0;
    const std::size_t expected_fields =
        static_cast<std::size_t>(ndim) + 1 +
        (requires_force ? static_cast<std::size_t>(ndim) : 0);
    std::vector<float> parsed_lse;
    std::vector<Gdata> parsed_force;
    parsed_lse.reserve(total);
    parsed_force.reserve(total);
    std::size_t line_number = 0;
    bool format_decided = false;
    bool versioned_format = false;
    std::string line;
    while (std::getline(input, line))
    {
        ++line_number;
        const std::vector<std::string> tokens = Tokenize_Meta_Line(line);
        if (tokens.empty()) continue;
        if (!format_decided)
        {
            format_decided = true;
            if (tokens[0] == "SPONGE_META_EDGE_V1")
            {
                versioned_format = true;
                if (tokens.size() != 6)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Read_Edge_File",
                        "Edge cache line %zu has %zu header fields; V1 "
                        "requires version, dimension, grid size, scatter "
                        "size, payload, and configuration fingerprint.",
                        line_number, tokens.size());
                }
                const int file_ndim = Parse_Meta_File_Int(
                    controller, tokens[1], line_number, "dimension",
                    "META::Read_Edge_File", "Edge cache");
                const int file_total = Parse_Meta_File_Int(
                    controller, tokens[2], line_number, "grid size",
                    "META::Read_Edge_File", "Edge cache");
                const int file_scatter = Parse_Meta_File_Int(
                    controller, tokens[3], line_number, "scatter size",
                    "META::Read_Edge_File", "Edge cache");
                const char* expected_payload =
                    requires_force ? "normal_lse_force" : "normal_lse";
                const std::uint64_t file_fingerprint = Parse_Meta_File_Hex_U64(
                    controller, tokens[5], line_number,
                    "configuration fingerprint", "META::Read_Edge_File",
                    "Edge cache");
                const std::uint64_t expected_fingerprint =
                    Meta_Edge_Config_Fingerprint(*this, requires_force);
                const bool configuration_matches =
                    file_ndim == ndim && file_total == total &&
                    file_scatter == scatter_size &&
                    tokens[4] == expected_payload &&
                    file_fingerprint == expected_fingerprint;
                if (!configuration_matches && !has_edge_file_input)
                {
                    return false;
                }
                if (!configuration_matches)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Read_Edge_File",
                        "Edge cache header declares %dD grid=%d scatter=%d "
                        "%s fingerprint=%016llx, but the active configuration "
                        "requires %dD grid=%d scatter=%d %s "
                        "fingerprint=%016llx.",
                        file_ndim, file_total, file_scatter, tokens[4].c_str(),
                        static_cast<unsigned long long>(file_fingerprint), ndim,
                        total, scatter_size, expected_payload,
                        static_cast<unsigned long long>(expected_fingerprint));
                }
                continue;
            }
            if (has_edge_file_input)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorBadFileFormat, "META::Read_Edge_File",
                    "An explicit edge_in_file must use the versioned V1 "
                    "format so its configuration can be validated.");
            }
            return false;
        }

        const std::size_t record = parsed_lse.size();
        if (record >= static_cast<std::size_t>(total) ||
            tokens.size() != expected_fields)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Read_Edge_File",
                "Edge cache line %zu has %zu fields at record %zu; exactly "
                "%zu fields and %d records are required.",
                line_number, tokens.size(), record, expected_fields, total);
        }
        const Axis expected_coordinates =
            mgrid->Get_Coordinates(static_cast<int>(record));
        for (int d = 0; d < ndim; ++d)
        {
            const float coordinate = Parse_Meta_File_Float(
                controller, tokens[d], line_number, "grid coordinate",
                "META::Read_Edge_File", "Edge cache");
            if (!Meta_Coordinate_Token_Matches(
                    tokens[d], coordinate, expected_coordinates[d],
                    versioned_format ? MetaCoordinateFormat::round_trip
                                     : MetaCoordinateFormat::legacy_edge))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "META::Read_Edge_File",
                    "Edge cache line %zu record %zu has coordinate %.9g in "
                    "dimension %d; grid order requires %.9g.",
                    line_number, record, static_cast<double>(coordinate), d,
                    static_cast<double>(expected_coordinates[d]));
            }
        }

        const float stored_lse = Parse_Meta_File_Float(
            controller, tokens[ndim], line_number, "normal_lse",
            "META::Read_Edge_File", "Edge cache");
        parsed_lse.push_back(stored_lse);
        Gdata force(ndim, 0.0f);
        if (requires_force)
        {
            for (int d = 0; d < ndim; ++d)
            {
                force[d] = Parse_Meta_File_Float(
                    controller,
                    tokens[static_cast<std::size_t>(ndim) + 1 +
                           static_cast<std::size_t>(d)],
                    line_number, "normalization derivative",
                    "META::Read_Edge_File", "Edge cache");
            }
        }
        parsed_force.push_back(force);
    }
    if (!format_decided) return false;
    if (parsed_lse.size() != static_cast<std::size_t>(total))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Read_Edge_File",
            "Edge cache contains %zu records; the active grid requires %d.",
            parsed_lse.size(), total);
    }

    potential = parsed_lse;
    mgrid->normal_lse = parsed_lse;
    sum_max = *std::max_element(parsed_lse.begin(), parsed_lse.end());
    if (requires_force)
    {
        for (int idx = 0; idx < total; ++idx)
        {
            for (int d = 0; d < ndim; ++d)
            {
                mgrid->normal_force[static_cast<std::size_t>(idx) * ndim + d] =
                    parsed_force[idx][d];
            }
        }
    }
    return true;
}
// Load hills from output file.
int META::Load_Hills(const std::string& fn)
{
    if (ndim <= 0 || sigmas.size() != static_cast<std::size_t>(ndim) ||
        periods.size() != static_cast<std::size_t>(ndim) ||
        (do_negative && (mgrid == NULL ||
                         mgrid->normal_lse.size() !=
                             static_cast<std::size_t>(mgrid->total_size) ||
                         hills.size() != vsink.size())))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Load_Hills",
            "Metadynamics history replay configuration is incomplete.");
    }
    const std::lock_guard<std::mutex> process_lock(
        sponge_meta_io::Process_Mutex());
    const fs::path hills_path =
        Meta_Native_Path(controller, fn, "META::Load_Hills");
#ifdef _WIN32
    int hills_descriptor = -1;
    const errno_t hills_open_status =
        _wsopen_s(&hills_descriptor, hills_path.c_str(), _O_RDONLY | _O_BINARY,
                  _SH_DENYNO, _S_IREAD);
    const int hills_open_error = static_cast<int>(hills_open_status);
    if (hills_open_status != 0)
#else
    errno = 0;
    const int hills_descriptor =
        open(hills_path.native().c_str(),
             O_RDONLY | O_CLOEXEC | O_NOFOLLOW | O_NONBLOCK);
    const int hills_open_error = Meta_Effective_IO_Error(errno);
    if (hills_descriptor < 0)
#endif
    {
        if (hills_open_error == ENOENT)
        {
            controller->printf(
                "    No existing hills record was found; metadynamics history "
                "starts empty\n");
            return 0;
        }
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, "META::Load_Hills",
            "Unable to open metadynamics hills input '%s': %s.", fn.c_str(),
            strerror(Meta_Effective_IO_Error(hills_open_error)));
    }
    int hills_identity_error = 0;
#ifdef _WIN32
    struct _stat64 hills_status = {};
    errno = 0;
    if (_fstat64(hills_descriptor, &hills_status) != 0)
    {
        hills_identity_error = Meta_Effective_IO_Error(errno);
    }
    else if ((hills_status.st_mode & _S_IFMT) != _S_IFREG)
    {
        hills_identity_error = EINVAL;
    }
#else
    struct stat hills_status = {};
    errno = 0;
    if (fstat(hills_descriptor, &hills_status) != 0)
    {
        hills_identity_error = Meta_Effective_IO_Error(errno);
    }
    else if (!S_ISREG(hills_status.st_mode))
    {
        hills_identity_error = EINVAL;
    }
#endif
    if (hills_identity_error != 0)
    {
#ifdef _WIN32
        _close(hills_descriptor);
#else
        close(hills_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, "META::Load_Hills",
            "Metadynamics hills input '%s' is not a regular file: %s.",
            fn.c_str(), strerror(hills_identity_error));
    }
    int hills_io_error = 0;
    if (!Lock_Meta_Append_File(hills_descriptor, false, &hills_io_error))
    {
#ifdef _WIN32
        _close(hills_descriptor);
#else
        close(hills_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, "META::Load_Hills",
            "Unable to lock metadynamics hills input '%s': %s.", fn.c_str(),
            strerror(hills_io_error));
    }
    const std::int64_t hills_size = Meta_Append_File_End(hills_descriptor);
    if (hills_size < 0)
    {
        hills_io_error = Meta_Effective_IO_Error(errno);
        Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
        _close(hills_descriptor);
#else
        close(hills_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "META::Load_Hills",
            "Unable to determine the size of hills input '%s': %s.", fn.c_str(),
            strerror(hills_io_error));
    }
    bool has_incomplete_tail = false;
    std::int64_t committed_size = hills_size;
    if (hills_size > 0)
    {
        char final_byte = 0;
        if (!Read_Meta_Append_At(hills_descriptor, hills_size - 1, &final_byte,
                                 1, &hills_io_error))
        {
            Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
            _close(hills_descriptor);
#else
            close(hills_descriptor);
#endif
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Load_Hills",
                "Unable to inspect the final byte of hills input '%s': %s.",
                fn.c_str(), strerror(hills_io_error));
        }
        has_incomplete_tail = final_byte != '\n';
        if (has_incomplete_tail &&
            !Find_Complete_Meta_Append_Prefix(hills_descriptor, hills_size,
                                              &committed_size, &hills_io_error))
        {
            Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
            _close(hills_descriptor);
#else
            close(hills_descriptor);
#endif
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "META::Load_Hills",
                "Unable to locate the complete prefix of hills input '%s': "
                "%s.",
                fn.c_str(), strerror(hills_io_error));
        }
    }
    const int cvsize = ndim;
    std::size_t expected_words = static_cast<std::size_t>(cvsize) + 1;
    if (do_negative) expected_words += 2;
    if (mscatter != NULL) expected_words += 1;
    if (mask) expected_words += 1;
    if (expected_words > (std::numeric_limits<std::size_t>::max() - 4096) / 64)
    {
        Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
        _close(hills_descriptor);
#else
        close(hills_descriptor);
#endif
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "META::Load_Hills",
            "The active metadynamics dimension cannot be represented as a "
            "bounded hills-record size.");
    }
    const std::size_t maximum_line_size = 4096 + expected_words * 64;
    Meta_Bounded_Line_Reader hills_reader(hills_descriptor, committed_size,
                                          maximum_line_size);
    std::vector<Hill> parsed_hills;
    Axis parsed_vsink;
    std::string line;
    std::size_t line_number = 0;
    int num_hills = 0;
    int last_projection_index = 0;
    float last_potential_max = 0.0f;
    try
    {
        for (;;)
        {
            hills_io_error = 0;
            const MetaLineReadStatus read_status =
                hills_reader.Next(&line, &hills_io_error);
            if (read_status == MetaLineReadStatus::end) break;
            ++line_number;
            if (read_status != MetaLineReadStatus::line)
            {
                Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
                _close(hills_descriptor);
#else
                close(hills_descriptor);
#endif
                if (read_status == MetaLineReadStatus::too_long)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Load_Hills",
                        "Hills line %zu exceeds the configuration-derived "
                        "limit of %zu bytes.",
                        line_number, maximum_line_size);
                }
                if (read_status == MetaLineReadStatus::allocation_error)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorMallocFailed, "META::Load_Hills",
                        "Unable to allocate at most %zu bytes for hills line "
                        "%zu.",
                        maximum_line_size, line_number);
                }
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "META::Load_Hills",
                    "Failed while reading hills line %zu from '%s': %s.",
                    line_number, fn.c_str(),
                    strerror(Meta_Effective_IO_Error(hills_io_error)));
            }
            Axis values;
            const std::vector<std::string> words = Tokenize_Meta_Line(line);
            if (words.empty()) continue;
            if (words.size() != expected_words)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "META::Load_Hills",
                    "Hills line %zu record %d has %zu fields; exactly %zu are "
                    "required by the active metadynamics configuration.",
                    line_number, num_hills, words.size(), expected_words);
            }
            for (int i = 0; i < cvsize; ++i)
            {
                values.push_back(Parse_Meta_File_Float(
                    controller, words[i], line_number, "hill center",
                    "META::Load_Hills", "Hills input"));
            }
            std::size_t cursor = static_cast<std::size_t>(cvsize);
            const float theight = Parse_Meta_File_Float(
                controller, words[cursor++], line_number, "hill height",
                "META::Load_Hills", "Hills input");
            if (do_negative)
            {
                const float p_max = Parse_Meta_File_Float(
                    controller, words[cursor++], line_number,
                    "maximum potential", "META::Load_Hills", "Hills input");
                const int p_id = Parse_Meta_File_Int(
                    controller, words[cursor++], line_number,
                    "sink projection index", "META::Load_Hills", "Hills input");
                const int projection_count =
                    mscatter != NULL ? scatter_size : mgrid->total_size;
                if (p_id < 0 || p_id >= projection_count)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Load_Hills",
                        "Hills line %zu sink projection index %d is outside "
                        "[0, %d).",
                        line_number, p_id, projection_count);
                }
                const int source_grid_index = mgrid->Get_Flat_Index(values);
                const int target_grid_index =
                    mscatter != NULL
                        ? mgrid->Get_Flat_Index(mscatter->Get_Coordinate(p_id))
                        : p_id;
                const float source_lse = mgrid->normal_lse[source_grid_index];
                const float target_lse = mgrid->normal_lse[target_grid_index];
                if (!Float_Memory_Is_Zero_Or_Normal(&source_lse) ||
                    !Float_Memory_Is_Zero_Or_Normal(&target_lse))
                {
                    controller->Throw_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Load_Hills",
                        "Hills sink projection references a non-finite or "
                        "subnormal edge normalization.");
                }
                const double coefficient =
                    static_cast<double>(p_max) +
                    static_cast<double>(dip) * CONSTANT_kB * temperature;
                const double lse_difference =
                    static_cast<double>(source_lse) - target_lse;
                const double projection_scale = std::exp(lse_difference);
                if (!Double_Memory_Is_Finite(&coefficient) ||
                    !Double_Memory_Is_Finite(&lse_difference) ||
                    !Double_Memory_Is_Finite(&projection_scale) ||
                    !(projection_scale > 0.0))
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Load_Hills",
                        "Hills line %zu produces an unrepresentable sink "
                        "projection scale.",
                        line_number);
                }
                const double sink_projection = coefficient * projection_scale;
                parsed_vsink.push_back(Store_Checked_Meta_Float(
                    controller, spongeErrorBadFileFormat, "META::Load_Hills",
                    "history sink projection", sink_projection));
                last_projection_index = p_id;
                last_potential_max = p_max;
            }
            if (mscatter != NULL)
            {
                const int scatter_index = Parse_Meta_File_Int(
                    controller, words[cursor++], line_number, "scatter index",
                    "META::Load_Hills", "Hills input");
                if (scatter_index < 0 || scatter_index >= scatter_size)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Load_Hills",
                        "Hills line %zu scatter index %d is outside [0, %d).",
                        line_number, scatter_index, scatter_size);
                }
            }
            if (mask)
            {
                const float committed_exit_tag = Parse_Meta_File_Float(
                    controller, words[cursor++], line_number,
                    "committed exit tag", "META::Load_Hills", "Hills input");
                if (committed_exit_tag < 0.0f)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Load_Hills",
                        "Hills line %zu has a negative committed exit tag.",
                        line_number);
                }
            }
            parsed_hills.emplace_back(values, sigmas, periods, theight);
            if (num_hills == INT_MAX)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow, "META::Load_Hills",
                    "Hills input contains more than INT_MAX records.");
            }
            ++num_hills;
        }
    }
    catch (const std::length_error&)
    {
        Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
        _close(hills_descriptor);
#else
        close(hills_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "META::Load_Hills",
            "Hills input '%s' exceeds representable in-memory history "
            "storage near line %zu.",
            fn.c_str(), line_number);
    }
    catch (const std::bad_alloc&)
    {
        Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
        _close(hills_descriptor);
#else
        close(hills_descriptor);
#endif
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorMallocFailed, "META::Load_Hills",
            "Unable to allocate in-memory history storage while reading '%s' "
            "near line %zu.",
            fn.c_str(), line_number);
    }
    if (has_incomplete_tail)
    {
        controller->printf(
            "    Ignoring incomplete final hills record on line %zu; the next "
            "append will remove it\n",
            line_number + 1);
    }
    Unlock_Meta_Append_File(hills_descriptor);
#ifdef _WIN32
    const int hills_close_status = _close(hills_descriptor);
#else
    const int hills_close_status = close(hills_descriptor);
#endif
    if (hills_close_status != 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOpenFileFailed, "META::Load_Hills",
            "Closing metadynamics hills input '%s' failed: %s.", fn.c_str(),
            strerror(Meta_Effective_IO_Error(errno)));
    }
    if (parsed_hills.size() != static_cast<std::size_t>(num_hills) ||
        (do_negative && parsed_vsink.size() != parsed_hills.size()) ||
        (!do_negative && !parsed_vsink.empty()))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Load_Hills",
            "Parsed hills and sink-projection histories are inconsistent.");
    }
    if (parsed_hills.size() > hills.max_size() - hills.size() ||
        parsed_vsink.size() > vsink.max_size() - vsink.size())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "META::Load_Hills",
            "Appending the parsed hills history exceeds vector capacity.");
    }
    try
    {
        sponge_meta_io::Publish_History(&hills, &vsink, &parsed_hills,
                                        parsed_vsink);
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "META::Load_Hills",
            "Parsed hills history exceeds representable vector storage.");
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMallocFailed, "META::Load_Hills",
            "Unable to reserve transactional hills-history storage.");
    }
    if (do_negative && num_hills > 0)
    {
        max_index = last_projection_index;
        potential_max = last_potential_max;
    }
    return num_hills;
}
float META::Calc_Hill(const Axis& values, const int i)
{
    float potential = 0;
    for (int j = 0; j < i; ++j)
    {
        Hill& hill = hills[j];
        const Gdata& tder = hill.Calc_Hill(values);
        potential += hill.potential * hill.height;
    }
    return potential;
}
float META::Sum_Hills(int history_freq)
{
    if (history_freq == 0)
    {
        return 0.;
    }
    TimePoint start_time, end_time;
    float duration;
    GetTime(start_time);
    int nhills = Load_Hills("myhill.log");
    if (nhills == 0)
    {
        controller->printf(
            "    No prior hills require reweighting; history replay is "
            "empty\n");
        return 0.0f;
    }
    controller->printf("Load hills file successfully, now calculate RCT...\n");
    Axis grid_coordinate;
    if (!use_scatter) grid_coordinate.resize(ndim);
    const std::string history_name = "history.log";
    Meta_Atomic_Output history_output(controller, history_name,
                                      "META::Sum_Hills");
    FILE* history_file = history_output.Stream();
    const int history_precision = std::numeric_limits<float>::max_digits10;
    // first loop: history
    float old_potential = 0.0f;
    Update_Reweighting_Factors("META::Sum_Hills");
    float total_gputime = 0.;
    for (int hill_index = 0; hill_index < nhills; ++hill_index)
    {
        Hill& hill = hills[hill_index];
        const Axis& values = hill.centers_;
        old_potential = Calc_Hill(values, hill_index);
        if (history_freq != 0 && (hill_index % history_freq == 0))
        {
            mgrid->potential.assign(mgrid->total_size, 0.0f);
            TimePoint tstart, tend;
            float gputime;
            GetTime(tstart);
            if (use_scatter)
            {
                for (int iter = 0; iter < scatter_size; ++iter)
                {
                    mscatter->potential[iter] =
                        Calc_Hill(mscatter->Get_Coordinate(iter), hill_index);
                }
                potential_max = 0.0f;
                const double Z_0 = Checked_Meta_Log_Partition(
                    controller, "META::Sum_Hills", minus_beta_f, &potential_max,
                    mscatter->potential);
                const double Z_V = Checked_Meta_Log_Partition(
                    controller, "META::Sum_Hills", minus_beta_f_plus_v,
                    &potential_max, mscatter->potential);
                const double checked_rct = static_cast<double>(CONSTANT_kB) *
                                           temperature * (Z_0 - Z_V);
                rct = Store_Checked_Meta_Float(
                    controller, spongeErrorValueErrorCommand, "META::Sum_Hills",
                    "history reweighting correction", checked_rct);
            }
            else  // use grid
            {
                for (int idx = 0; idx < mgrid->total_size; ++idx)
                {
                    mgrid->Get_Coordinates(idx, grid_coordinate.data());
                    mgrid->potential[idx] =
                        Calc_Hill(grid_coordinate, hill_index);
                }
                potential_max = 0.0f;
                const double Z_0 = Checked_Meta_Log_Partition(
                    controller, "META::Sum_Hills", minus_beta_f, &potential_max,
                    mgrid->potential);
                const double Z_V = Checked_Meta_Log_Partition(
                    controller, "META::Sum_Hills", minus_beta_f_plus_v,
                    &potential_max, mgrid->potential);
                const double checked_rct = static_cast<double>(CONSTANT_kB) *
                                           temperature * (Z_0 - Z_V);
                rct = Store_Checked_Meta_Float(
                    controller, spongeErrorValueErrorCommand, "META::Sum_Hills",
                    "history reweighting correction", checked_rct);
            }
            GetTime(tend);
            GetDuration(tend, tstart, gputime);
            total_gputime += gputime;
            const float history_rbias = Store_Checked_Meta_Float(
                controller, spongeErrorValueErrorCommand, "META::Sum_Hills",
                "history reweighted bias",
                static_cast<double>(old_potential) - rct);
            if (do_negative)
            {
                if (static_cast<std::size_t>(hill_index) >= vsink.size())
                {
                    controller->Throw_SPONGE_Error(
                        spongeErrorBadFileFormat, "META::Sum_Hills",
                        "Hills history is missing a required sink projection.");
                }
                fprintf(history_file, "%.*g\t%.*g\t%.*g\t%.*g\n",
                        history_precision, static_cast<double>(old_potential),
                        history_precision, static_cast<double>(history_rbias),
                        history_precision, static_cast<double>(rct),
                        history_precision,
                        static_cast<double>(vsink[hill_index]));
            }
            else
            {
                fprintf(history_file, "%.*g\t%.*g\t%.*g\n", history_precision,
                        static_cast<double>(old_potential), history_precision,
                        static_cast<double>(history_rbias), history_precision,
                        static_cast<double>(rct));
            }
        }
    }
    history_output.Commit();
    GetTime(end_time);
    GetDuration(end_time, start_time, duration);
    int hours = floor(duration / 3600);
    float nohour = duration - 3600 * hours;
    int mins = floor(nohour / 60);
    float seconds = nohour - 60 * mins;
    controller->printf(
        "The RBIAS & RCT calculation cost %f of %f seconds: %d hour %d min %f "
        "second\n",
        duration > 0.0f ? total_gputime / duration : 0.0f, duration, hours,
        mins, seconds);
    return old_potential;
}
void META::Edge_Effect(const int dim, const int scatter_size)
{
    std::vector<float> potential_from_file;
    const char* file_name = edge_file_name.c_str();

    int total = mgrid->total_size;
    if (scatter_size == total)
    {
        double normalization = 1.0;
        const double sqrtpi = std::sqrt(static_cast<double>(CONSTANT_Pi) * 2.0);
        for (int i = 0; i < ndim; i++)
        {
            normalization /=
                static_cast<double>(cv_deltas[i]) * sigmas[i] / sqrtpi;
        }
        if (!Double_Memory_Is_Finite(&normalization) || !(normalization > 0.0))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "META::Edge_Effect",
                "Metadynamics edge normalization is not finite and positive.");
        }
        const double log_normalization = std::log(normalization);
        if (!Double_Memory_Is_Finite(&log_normalization) ||
            log_normalization > std::numeric_limits<float>::max() ||
            log_normalization < -std::numeric_limits<float>::max())
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "META::Edge_Effect",
                "The logarithm of the metadynamics edge normalization is "
                "outside the finite float range.");
        }
        const float stored_log_normalization =
            static_cast<float>(log_normalization);
        if (!Float_Memory_Is_Zero_Or_Normal(&stored_log_normalization))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "META::Edge_Effect",
                "The logarithm of the metadynamics edge normalization is not "
                "representable as a zero-or-normal float.");
        }
        mgrid->normal_lse.assign(mgrid->total_size, stored_log_normalization);
    }
    bool readsuccess = Read_Edge_File(file_name, potential_from_file);
    if (has_edge_file_input && !readsuccess)
    {
        controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                       "META::Edge_Effect",
                                       "Failed to read edge_in_file");
    }
    if (!readsuccess)
    {
        controller->printf("Calculation the %d grid of edge effect\n", total);
        const bool requires_force = do_negative || mask != 0;
        const int precision = std::numeric_limits<float>::max_digits10;
        const std::uint64_t config_fingerprint =
            Meta_Edge_Config_Fingerprint(*this, requires_force);
        Axis esigmas;
        float adjust_factor = 1.0;
        for (int i = 0; i < ndim; ++i)
        {
            esigmas.push_back(sigmas[i] * adjust_factor);
        }
        mgrid->normal_lse.assign(mgrid->total_size, 0.0f);
        mgrid->normal_force.assign(mgrid->total_size * ndim, 0.0f);
        mgrid->Alloc_Device();
        mscatter->Alloc_Device();
        deviceMemcpy(d_hill_inv_w, esigmas.data(), sizeof(float) * ndim,
                     deviceMemcpyHostToDevice);
        int edge_error = 0;
        int* d_edge_error = NULL;
        Device_Malloc_And_Copy_Safely((void**)&d_edge_error, &edge_error,
                                      sizeof(edge_error));
        Launch_Device_Kernel(
            Update_Edge_Effect_Grid, total, 32, 0, NULL, total, ndim,
            scatter_size, mgrid->d_num_points, mgrid->d_lower, mgrid->d_spacing,
            d_hill_inv_w, mscatter->d_periods, mscatter->d_coordinates,
            static_cast<int>(requires_force), mgrid->d_normal_lse,
            mgrid->d_normal_force, d_edge_error);
        deviceMemcpy(&edge_error, d_edge_error, sizeof(edge_error),
                     deviceMemcpyDeviceToHost);
#ifndef CPU_ARCH_NAME
        deviceFree(d_edge_error);
#endif
        if (edge_error != 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Edge_Effect",
                "Edge normalization is not representable in stable double "
                "log-space or zero-or-normal float storage.");
        }
        mgrid->Sync_To_Host();
        sum_max = -std::numeric_limits<float>::max();
        for (int gidx = 0; gidx < mgrid->total_size; ++gidx)
        {
            float logsumhills = mgrid->normal_lse[gidx];
            if (!Float_Memory_Is_Zero_Or_Normal(&logsumhills))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "META::Edge_Effect",
                    "Edge normalization at grid record %d is non-finite or "
                    "subnormal.",
                    gidx);
            }
            sum_max = fmaxf(logsumhills, sum_max);
            if (requires_force)
            {
                float* nf_data = &mgrid->normal_force[gidx * ndim];
                for (int i = 0; i < ndim; ++i)
                {
                    if (!Float_Memory_Is_Zero_Or_Normal(&nf_data[i]))
                    {
                        controller->Throw_Formatted_SPONGE_Error(
                            spongeErrorSimulationBreakDown, "META::Edge_Effect",
                            "Edge normalization derivative at grid record %d "
                            "dimension %d is non-finite or subnormal.",
                            gidx, i);
                    }
                }
            }
        }

        Axis values(ndim);
        Meta_Atomic_Output atomic_output(controller, edge_file_name,
                                         "META::Edge_Effect");
        FILE* temp_file = atomic_output.Stream();
        fprintf(temp_file, "SPONGE_META_EDGE_V1 %d %d %d %s %016llx\n", ndim,
                total, scatter_size,
                requires_force ? "normal_lse_force" : "normal_lse",
                static_cast<unsigned long long>(config_fingerprint));
        for (int gidx = 0; gidx < mgrid->total_size; ++gidx)
        {
            mgrid->Get_Coordinates(gidx, values.data());
            for (float value : values)
            {
                fprintf(temp_file, "%.*g\t", precision,
                        static_cast<double>(value));
            }
            fprintf(temp_file, "%.*g", precision,
                    static_cast<double>(mgrid->normal_lse[gidx]));
            if (requires_force)
            {
                const float* nf_data = &mgrid->normal_force[gidx * ndim];
                for (int i = 0; i < ndim; ++i)
                {
                    fprintf(temp_file, "\t%.*g", precision,
                            static_cast<double>(nf_data[i]));
                }
            }
            fprintf(temp_file, "\n");
        }
        atomic_output.Commit();
    }
    if (dim == 1)
    {
        Pick_Scatter("lnbias.dat");
    }
}

void META::Pick_Scatter(const std::string fn)
{
    Meta_Atomic_Output atomic_output(controller, fn, "META::Pick_Scatter");
    FILE* output = atomic_output.Stream();
    const int precision = std::numeric_limits<float>::max_digits10;
    fprintf(output, "SPONGE_META_SCATTER_LOG_V1 %d\n", scatter_size);
    for (int index = 0; index < scatter_size; ++index)
    {
        const Axis& neighbor = mscatter->Get_Coordinate(index);
        const float lnbias = mgrid->normal_lse[mgrid->Get_Flat_Index(neighbor)];
        if (!Float_Memory_Is_Zero_Or_Normal(&lnbias))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Pick_Scatter",
                "Scatter normalization at index %d is not representable.",
                index);
        }
        fprintf(output, "%d\t%.*g\n", index, precision,
                static_cast<double>(lnbias));
    }
    atomic_output.Commit();
}
float META::Normalization(const Axis& values, float factor, bool do_normalise)
{
    if (!Float_Memory_Is_Zero_Or_Normal(&factor))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Normalization",
            "Metadynamics normalization factor is non-finite or subnormal.");
    }
    if (!do_normalise || factor == 0.0f) return factor;
    int reference_index = mgrid->Get_Flat_Index(values);
    if (usegrid)
    {
        reference_index = 0;
    }
    else if (convmeta)
    {
        if (mscatter == NULL || max_index < 0 || max_index >= scatter_size)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Normalization",
                "Metadynamics normalization has an invalid scatter index.");
        }
        reference_index =
            mgrid->Get_Flat_Index(mscatter->Get_Coordinate(max_index));
    }
    const float reference_lse = mgrid->normal_lse[reference_index];
    if (!Float_Memory_Is_Zero_Or_Normal(&reference_lse))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Normalization",
            "Metadynamics normalization references an invalid log value.");
    }
    const double scale = std::exp(-static_cast<double>(reference_lse));
    if (!Double_Memory_Is_Finite(&scale) || !(scale > 0.0))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Normalization",
            "Metadynamics normalization scale is not representable in stable "
            "double log space.");
    }
    return Store_Checked_Meta_Float(controller, spongeErrorSimulationBreakDown,
                                    "META::Normalization", "normalized factor",
                                    static_cast<double>(factor) * scale);
}
void META::Validate_Runtime_Temperature(const char* error_by) const
{
    if (!Float_Memory_Is_Finite(&temperature) ||
        !Float_Memory_Is_Normal(&temperature) || !(temperature > 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            is_initialized ? spongeErrorSimulationBreakDown
                           : spongeErrorValueErrorCommand,
            error_by,
            "Metadynamics requires a finite, positive, normal runtime "
            "temperature; got %.9g K.",
            static_cast<double>(temperature));
    }
}

float META::Update_Reweighting_Factors(const char* error_by)
{
    Validate_Runtime_Temperature(error_by);
    const int error_number = is_initialized ? spongeErrorSimulationBreakDown
                                            : spongeErrorValueErrorCommand;
    const double thermal_energy =
        static_cast<double>(CONSTANT_kB) * static_cast<double>(temperature);
    const double factor_denominator =
        static_cast<double>(welltemp_factor) - 1.0;
    const double beta = 1.0 / thermal_energy;
    const double beta_f_plus_v = beta / factor_denominator;
    const double beta_f = static_cast<double>(welltemp_factor) * beta_f_plus_v;
    const double derived[] = {thermal_energy, factor_denominator, beta,
                              beta_f_plus_v, beta_f};
    for (double value : derived)
    {
        if (!Double_Memory_Is_Finite(&value) || !(value > 0.0))
        {
            controller->Throw_Formatted_SPONGE_Error(
                error_number, error_by,
                "Metadynamics cannot derive finite positive thermal and "
                "well-tempered factors from temperature %.9g K and bias "
                "factor %.9g.",
                static_cast<double>(temperature),
                static_cast<double>(welltemp_factor));
        }
    }
    if (beta > std::numeric_limits<float>::max() ||
        beta_f_plus_v > std::numeric_limits<float>::max() ||
        beta_f > std::numeric_limits<float>::max())
    {
        controller->Throw_Formatted_SPONGE_Error(
            error_number, error_by,
            "Metadynamics inverse thermal energy is not representable as a "
            "finite float at temperature %.9g K.",
            static_cast<double>(temperature));
    }
    const float checked_beta = static_cast<float>(beta);
    const float checked_beta_f_plus_v = static_cast<float>(beta_f_plus_v);
    const float checked_beta_f = static_cast<float>(beta_f);
    if (!Float_Memory_Is_Normal(&checked_beta) ||
        !Float_Memory_Is_Normal(&checked_beta_f_plus_v) ||
        !Float_Memory_Is_Normal(&checked_beta_f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            error_number, error_by,
            "Metadynamics inverse thermal and well-tempered factors must be "
            "representable as positive normal floats at temperature %.9g K "
            "and bias factor %.9g.",
            static_cast<double>(temperature),
            static_cast<double>(welltemp_factor));
    }
    minus_beta_f_plus_v = checked_beta_f_plus_v;
    minus_beta_f = checked_beta_f;
    return checked_beta;
}

void META::Get_Height(const Axis& values)
{
    Validate_Runtime_Temperature("META::Get_Height");
    Estimate(values, true, false);
    height = height_0;
    if (is_welltemp != 1)
    {
        return;
    }
    const double denominator = (static_cast<double>(welltemp_factor) - 1.0) *
                               static_cast<double>(CONSTANT_kB) *
                               static_cast<double>(temperature);
    const double exponent =
        -static_cast<double>(potential_backup) / denominator;
    const double tempered_height =
        static_cast<double>(height_0) * std::exp(exponent);
    if (!Double_Memory_Is_Finite(&denominator) || !(denominator > 0.0) ||
        !Double_Memory_Is_Finite(&exponent) ||
        !Double_Memory_Is_Finite(&tempered_height) ||
        (height_0 != 0.0f && tempered_height == 0.0) ||
        tempered_height > std::numeric_limits<float>::max() ||
        tempered_height < -std::numeric_limits<float>::max())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Get_Height",
            "Metadynamics produced an unrepresentable well-tempered hill "
            "height at temperature %.9g K.",
            static_cast<double>(temperature));
    }
    const float checked_height = static_cast<float>(tempered_height);
    if (!Float_Memory_Is_Zero_Or_Normal(&checked_height) ||
        (tempered_height != 0.0 && checked_height == 0.0f))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Get_Height",
            "Metadynamics well-tempered hill height underflowed to an "
            "unrepresentable float.");
    }
    height = checked_height;
}

float META::Calc_V_Shift(const Axis& values)
{
    if (!do_negative)
    {
        return 0.;
    }
    const int nidx = mgrid->Get_Flat_Index(values);
    const float point_lse = mgrid->normal_lse[nidx];
    if (!Float_Memory_Is_Zero_Or_Normal(&point_lse) ||
        !Float_Memory_Is_Zero_Or_Normal(&normalization_reference_lse) ||
        !Double_Memory_Is_Finite(&normalization_base_factor))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Calc_V_Shift",
            "Sink shift references a non-finite or subnormal normalization.");
    }
    const double exponent =
        static_cast<double>(point_lse) - normalization_reference_lse;
    const double scale = std::exp(exponent);
    if (!Double_Memory_Is_Finite(&exponent) ||
        !Double_Memory_Is_Finite(&scale) || !(scale > 0.0))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Calc_V_Shift",
            "Sink shift scale is not representable in stable double log "
            "space.");
    }
    double result = normalization_base_factor * scale;
    if (convmeta)
    {
        return Store_Checked_Meta_Float(
            controller, spongeErrorSimulationBreakDown, "META::Calc_V_Shift",
            "sink shift", result);
    }
    result *= static_cast<double>(point_lse) - sum_max;
    return Store_Checked_Meta_Float(controller, spongeErrorSimulationBreakDown,
                                    "META::Calc_V_Shift",
                                    "growth-reweighted sink shift", result);
}
void META::Get_Reweighting_Bias(float temp)
{
    const float beta = Update_Reweighting_Factors("META::Get_Reweighting_Bias");
    const float sampled_potential =
        Store_Checked_Meta_Float(controller, spongeErrorSimulationBreakDown,
                                 "META::Get_Reweighting_Bias",
                                 "sampled bias potential", potential_local);
    bias = sampled_potential;
    const double sampled_backup = potential_backup;
    const float sampled_new_max = new_max;
    const double sampled_normalization_base = normalization_base_factor;
    const float sampled_normalization_reference = normalization_reference_lse;
    const float checked_shift = Store_Checked_Meta_Float(
        controller, spongeErrorSimulationBreakDown,
        "META::Get_Reweighting_Bias", "sampled bias shift", temp);
    double Z_0_sink = 0.0;
    double Z_V_sink = 0.0;
    bool has_partition_term = false;
    float next_potential_max = -std::numeric_limits<float>::max();
    int next_max_index = 0;
    auto accumulate_point = [&](const Axis& coor, int point_index,
                                double log_weight, const char* representation)
    {
        Estimate(coor, true, false);
        const float checked_potential = Store_Checked_Meta_Float(
            controller, spongeErrorSimulationBreakDown,
            "META::Get_Reweighting_Bias", "bias potential", potential_backup);
        const float checked_vshift = Store_Checked_Meta_Float(
            controller, spongeErrorSimulationBreakDown,
            "META::Get_Reweighting_Bias", "bias shift", Calc_V_Shift(coor));
        const double exponent_0 =
            log_weight + static_cast<double>(minus_beta_f) * checked_potential;
        const double exponent_v =
            log_weight +
            static_cast<double>(minus_beta_f_plus_v) * checked_potential +
            static_cast<double>(beta) * checked_vshift;
        if (!Double_Memory_Is_Finite(&exponent_0) ||
            !Double_Memory_Is_Finite(&exponent_v))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Get_Reweighting_Bias",
                "The %s partition exponent at point %d is non-finite.",
                representation, point_index);
        }
        if (!has_partition_term)
        {
            Z_0_sink = exponent_0;
            Z_V_sink = exponent_v;
            has_partition_term = true;
        }
        else
        {
            Z_0_sink = Checked_Meta_Log_Add_Exp(
                controller, "META::Get_Reweighting_Bias", Z_0_sink, exponent_0);
            Z_V_sink = Checked_Meta_Log_Add_Exp(
                controller, "META::Get_Reweighting_Bias", Z_V_sink, exponent_v);
        }
        if (checked_potential > next_potential_max)
        {
            next_potential_max = checked_potential;
            next_max_index = point_index;
        }
    };
    if (mscatter != NULL)
    {
        for (int iter = 0; iter < scatter_size; ++iter)
        {
            accumulate_point(mscatter->Get_Coordinate(iter), iter, 0.0,
                             "scatter");
        }
    }
    else if (mgrid != NULL)
    {
        double log_cell_volume = 0.0;
        for (int d = 0; d < ndim; ++d)
        {
            const double spacing = mgrid->spacing[d];
            if (!Double_Memory_Is_Finite(&spacing) || !(spacing > 0.0))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "META::Get_Reweighting_Bias",
                    "Grid spacing in dimension %d is not finite and positive.",
                    d);
            }
            log_cell_volume += std::log(spacing);
        }
        if (!Double_Memory_Is_Finite(&log_cell_volume))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Get_Reweighting_Bias",
                "The metadynamics grid cell volume is not representable in "
                "log space.");
        }
        for (int idx = 0; idx < mgrid->total_size; ++idx)
        {
            accumulate_point(mgrid->Get_Coordinates(idx), idx, log_cell_volume,
                             "grid");
        }
    }
    if (!has_partition_term)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Get_Reweighting_Bias",
            "Metadynamics reweighting has no grid or scatter partition terms.");
    }
    potential_max =
        Store_Checked_Meta_Float(controller, spongeErrorSimulationBreakDown,
                                 "META::Get_Reweighting_Bias",
                                 "maximum bias potential", next_potential_max);
    max_index = next_max_index;
    const double checked_rct =
        static_cast<double>(CONSTANT_kB) * temperature * (Z_0_sink - Z_V_sink);
    rct = Store_Checked_Meta_Float(controller, spongeErrorSimulationBreakDown,
                                   "META::Get_Reweighting_Bias",
                                   "reweighting correction", checked_rct);
    rbias = Store_Checked_Meta_Float(
        controller, spongeErrorSimulationBreakDown,
        "META::Get_Reweighting_Bias", "reweighted bias",
        sampled_backup - static_cast<double>(rct) - checked_shift);
    potential_local = sampled_potential;
    potential_backup = static_cast<float>(sampled_backup);
    new_max = sampled_new_max;
    normalization_base_factor = sampled_normalization_base;
    normalization_reference_lse = sampled_normalization_reference;
}

void META::Add_Potential(float temp, int steps)
{
    if (Will_Add_Potential(steps))
    {
        Axis values;
        for (int i = 0; i < cvs.size(); ++i)
        {
            values.push_back(cvs[i]->value);
        }
        Get_Height(values);
        float vshift = Calc_V_Shift(values);
        Get_Reweighting_Bias(vshift);
        Hill hill = Hill(values, sigmas, periods, height);
        hills.push_back(hill);
        const bool has_scatter = mscatter != NULL;
        const int scatter_index =
            has_scatter ? mscatter->Get_Index(values) : -1;
        // Record the committed sampling-history count.  Transactional force
        // evaluations restore exit_tag and therefore cannot inflate it.
        hilllog(controller, "myhill.log", values, height, do_negative,
                potential_max, max_index, has_scatter, scatter_index, mask != 0,
                exit_tag);
        exit_tag = 0.0;
        if (!kde && subhill)
        {
            const Gdata& tder = hill.Calc_Hill(values);
            if (mgrid != NULL)
            {
                mgrid->potential[mgrid->Get_Flat_Index(values)] +=
                    height * hill.potential;
            }
            else if (mscatter != NULL)
            {
                int sidx = mscatter->Get_Index(values);
                mscatter->potential[sidx] += height * hill.potential;
                deviceMemcpy(mscatter->d_potential + sidx,
                             &mscatter->potential[sidx], sizeof(float),
                             deviceMemcpyHostToDevice);
            }
            return;
        }
        float factor = Normalization(values, height,
                                     kde);  // height with normalized factor
        if (use_scatter)
        {
            deviceMemcpy(d_hill_centers, hill.centers_.data(),
                         sizeof(float) * ndim, deviceMemcpyHostToDevice);
            deviceMemcpy(d_hill_inv_w, hill.inv_w_.data(), sizeof(float) * ndim,
                         deviceMemcpyHostToDevice);
            int update_force = (!mscatter->force.empty()) ? 1 : 0;
            Launch_Device_Kernel(
                Update_Scatter_With_Hill,
                (static_cast<std::size_t>(scatter_size) + 255) / 256, 256, 0,
                NULL, scatter_size, ndim, mscatter->d_coordinates,
                mscatter->d_periods, d_hill_centers, d_hill_inv_w, factor,
                update_force, do_cutoff ? 1 : 0, d_cutoff,
                mscatter->d_potential, mscatter->d_force);
            mscatter->Sync_To_Host();
        }
        // Update grid potential and force with hill on device
        if (mgrid != NULL)
        {
            deviceMemcpy(d_hill_centers, hill.centers_.data(),
                         sizeof(float) * ndim, deviceMemcpyHostToDevice);
            deviceMemcpy(d_hill_inv_w, hill.inv_w_.data(), sizeof(float) * ndim,
                         deviceMemcpyHostToDevice);
            deviceMemcpy(d_hill_periods, periods.data(), sizeof(float) * ndim,
                         deviceMemcpyHostToDevice);
            int update_force = (!subhill && !mgrid->force.empty()) ? 1 : 0;
            Launch_Device_Kernel(
                Update_Grid_With_Hill,
                (static_cast<std::size_t>(mgrid->total_size) + 255) / 256, 256,
                0, NULL, mgrid->total_size, ndim, mgrid->d_num_points,
                mgrid->d_lower, mgrid->d_spacing, d_hill_centers, d_hill_inv_w,
                d_hill_periods, factor, update_force, mgrid->d_potential,
                mgrid->d_force);
            mgrid->Sync_To_Host();
        }
    }
}

bool META::Will_Add_Potential(int step) const
{
    if (!is_initialized) return false;
    if (potential_update_interval <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Will_Add_Potential",
            "Metadynamics potential_update_interval became %d after "
            "initialization; refusing to silently disable hill deposition.",
            potential_update_interval);
        return false;
    }
    if (step < 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "META::Will_Add_Potential",
            "Metadynamics received invalid negative simulation step %d.", step);
        return false;
    }
    return step % potential_update_interval == 0;
}

void META::Initial(CONTROLLER* controller,
                   COLLECTIVE_VARIABLE_CONTROLLER* cv_controller,
                   char* module_name, float sys_temperature)
{
    this->controller = controller;
    if (module_name == NULL)
    {
        Copy_Meta_String(controller, this->module_name,
                         sizeof(this->module_name), "meta", "module name");
    }
    else
    {
        Copy_Meta_String(controller, this->module_name,
                         sizeof(this->module_name), module_name, "module name");
    }
    if (!cv_controller->Command_Exist(this->module_name, "CV"))
    {
        controller->printf("META IS NOT INITIALIZED\n\n");
        return;
    }
    else
    {
        std::vector<std::string> cv_str =
            cv_controller->Ask_For_String_Parameter(this->module_name, "CV",
                                                    ndim);
        std::string cvv =
            std::accumulate(cv_str.begin(), cv_str.end(), std::string(""));
        controller->printf("%s contains %d dimension META\n", cvv.c_str(),
                           ndim);
    }
    temperature = sys_temperature;
    Validate_Runtime_Temperature("METADYNAMICS::Initial");
    if (cv_controller->Command_Exist(this->module_name, "dip"))
    {
        dip = cv_controller->Ask_For_Float_Parameter(this->module_name, "dip",
                                                     1)[0];
        if (!Float_Memory_Is_Zero_Or_Normal(&dip))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "dip must be a finite zero-or-normal float");
        }
    }
    if (cv_controller->Command_Exist(this->module_name, "welltemp_factor"))
    {
        float* temp_value = cv_controller->Ask_For_Float_Parameter(
            this->module_name, "welltemp_factor");
        welltemp_factor = temp_value[0];
        free(temp_value);
        if (Float_Memory_Is_Finite(&welltemp_factor) &&
            Float_Memory_Is_Normal(&welltemp_factor) && welltemp_factor > 1.0f)
        {
            is_welltemp = 1;
        }
        else
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "welltemp_factor must be finite, normal, and greater than "
                "1");
        }
    }
    Update_Reweighting_Factors("METADYNAMICS::Initial");
    cvs = cv_controller->Ask_For_CV(this->module_name, -1);
    if (cvs.empty() || cvs.size() > static_cast<std::size_t>(INT_MAX))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
            "Metadynamics requires between 1 and %d CV dimensions; got %zu.",
            INT_MAX, cvs.size());
    }
    for (std::size_t i = 0; i < cvs.size(); ++i)
    {
        if (cvs[i] == NULL)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "Metadynamics CV dimension %zu resolved to a null object.", i);
        }
    }
    if (cv_controller->Command_Exist(this->module_name, "Ndim"))
    {
        ndim = *cv_controller->Ask_For_Int_Parameter(this->module_name, "Ndim");
        if (ndim != cvs.size())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand, "METADYNAMICS::Initial",
                "Configured Ndim is %d, but the CV list contains %zu "
                "dimensions.",
                ndim, cvs.size());
        }
    }
    else
    {
        ndim = cvs.size();
    }
    controller->printf("START INITIALIZING %dD-META:\n", ndim);
    read_potential_file_name = "Meta_Potential.txt";
    write_potential_file_name = "Meta_Potential.txt";
    if (controller->Command_Exist("default_in_file_prefix"))
    {
        read_potential_file_name = std::string(controller->Original_Command(
                                       "default_in_file_prefix")) +
                                   "_Meta_Potential.txt";
    }
    if (controller->Command_Exist("default_out_file_prefix"))
    {
        write_potential_file_name = std::string(controller->Original_Command(
                                        "default_out_file_prefix")) +
                                    "_Meta_Potential.txt";
    }
    write_directly_file_name = "Meta_directly.txt";
    edge_file_name = "sumhill.log";
    has_edge_file_input = false;
    if (cv_controller->Command_Exist(this->module_name, "edge_in_file"))
    {
        has_edge_file_input = true;
        edge_file_name = cv_controller->Ask_For_String_Parameter(
            this->module_name, "edge_in_file")[0];
        if (edge_file_name.empty())
        {
            controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                           "METADYNAMICS::Initial",
                                           "edge_in_file must not be empty");
        }
    }

    height_0 = 1.0f;
    if (cv_controller->Command_Exist(this->module_name, "height"))
    {
        float* temp_value =
            cv_controller->Ask_For_Float_Parameter(this->module_name, "height");
        height_0 = temp_value[0];
        free(temp_value);
    }
    if (!Float_Memory_Is_Zero_Or_Normal(&height_0))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
            "height must be a finite zero-or-normal float");
    }
    if (cv_controller->Command_Exist(this->module_name, "wall_height"))
    {
        do_borderwall = true;
        float* temp_value = cv_controller->Ask_For_Float_Parameter(
            this->module_name, "wall_height");
        border_potential_height = temp_value[0];
        free(temp_value);
        if (!Float_Memory_Is_Normal(&border_potential_height) ||
            !(border_potential_height > 0.0f))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "wall_height must be a finite positive normal float");
        }
    }
    if (cv_controller->Command_Exist(this->module_name, "potential_out_file"))
    {
        write_potential_file_name = cv_controller->Ask_For_String_Parameter(
            this->module_name, "potential_out_file")[0];
        if (write_potential_file_name.empty())
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "potential_out_file must not be empty");
        }
    }

    bool has_potential_update_interval = false;
    if (cv_controller->Command_Exist(this->module_name,
                                     "potential_update_interval"))
    {
        int* temp_value = cv_controller->Ask_For_Int_Parameter(
            this->module_name, "potential_update_interval");
        potential_update_interval = temp_value[0];
        has_potential_update_interval = true;
        free(temp_value);
    }
    if (controller->Command_Exist("write_information_interval"))
    {
        controller->Check_Int("write_information_interval",
                              "METADYNAMICS::Initial");
        write_information_interval =
            atoi(controller->Command("write_information_interval"));
    }
    else
    {
        write_information_interval = 1000;
    }
    if (write_information_interval <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
            "write_information_interval must be positive when "
            "metadynamics is enabled; got %d.",
            write_information_interval);
    }
    if (!has_potential_update_interval)
    {
        potential_update_interval = write_information_interval;
        controller->printf(
            "    Potential update interval defaults to the validated "
            "write_information_interval (%d)\n",
            potential_update_interval);
    }
    if (potential_update_interval <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
            "potential_update_interval must be positive; got %d.",
            potential_update_interval);
    }

    if (cv_controller->Command_Exist(this->module_name, "subhill"))
    {
        subhill = true;
        controller->printf("    reading subhill for meta: 1\n");
    }
    if (cv_controller->Command_Exist(this->module_name, "kde"))
    {
        int kde_dim = cv_controller->Ask_For_Int_Parameter(this->module_name,
                                                           "kde", 1)[0];
        if (kde_dim)
        {
            kde = true;
            subhill = true;
            controller->printf("    reading kde's subhill for meta: %d\n",
                               kde_dim);
        }
    }
    if (cv_controller->Command_Exist(this->module_name, "mask"))
    {
        mask = cv_controller->Ask_For_Int_Parameter(this->module_name, "mask",
                                                    1)[0];
        if (mask)
        {
            controller->printf("    reading mask dimension meta: %d\n", mask);
            if (cv_controller->Command_Exist(this->module_name, "max_force"))
            {
                max_force = cv_controller->Ask_For_Float_Parameter(
                    this->module_name, "max_force", 1)[0];
                if (!Float_Memory_Is_Normal(&max_force) || !(max_force > 0.0f))
                {
                    controller->Throw_SPONGE_Error(
                        spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                        "max_force must be a finite positive normal float");
                }
            }
        }
    }
    if (cv_controller->Command_Exist(this->module_name, "sink"))
    {
        int sub_dim = cv_controller->Ask_For_Int_Parameter(this->module_name,
                                                           "sink", 1)[0];
        if (sub_dim > 0)
        {
            do_negative = true;
            controller->printf(
                "    reading sink/submarine dimension for meta: %d\n", sub_dim);
        }
    }
    if (cv_controller->Command_Exist(this->module_name, "sumhill_freq"))
    {
        history_freq = cv_controller->Ask_For_Int_Parameter(
            this->module_name, "sumhill_freq", 1)[0];
        if (history_freq < 0)
        {
            controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                           "METADYNAMICS::Initial",
                                           "sumhill_freq must be nonnegative");
        }
    }
    if (cv_controller->Command_Exist(this->module_name, "convmeta"))
    {
        do_negative = true;
        convmeta = cv_controller->Ask_For_Int_Parameter(this->module_name,
                                                        "convmeta", 1)[0];
    }
    if (cv_controller->Command_Exist(this->module_name, "grw"))
    {
        do_negative = true;
        grw = cv_controller->Ask_For_Int_Parameter(this->module_name, "grw",
                                                   1)[0];
    }
    cv_periods = cv_controller->Ask_For_Float_Parameter(
        this->module_name, "CV_period", cvs.size(), 1, false);
    cv_sigmas = cv_controller->Ask_For_Float_Parameter(this->module_name,
                                                       "CV_sigma", cvs.size());
    cutoff = cv_controller->Ask_For_Float_Parameter(
        this->module_name, "CV_sigma", cvs.size(), 1, false, 0., -3);
    if (cv_controller->Command_Exist(this->module_name, "cutoff"))
    {
        do_cutoff = true;
        cutoff = cv_controller->Ask_For_Float_Parameter(this->module_name,
                                                        "cutoff", cvs.size());
    }
    for (int i = 0; i < cvs.size(); i++)
    {
        if (!Float_Memory_Is_Normal(&cv_sigmas[i]) || !(cv_sigmas[i] > 0.0f))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "CV_sigma must be a finite positive normal float");
        }
        if (!Float_Memory_Is_Zero_Or_Normal(&cv_periods[i]) ||
            cv_periods[i] < 0.0f)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "CV_period must be zero or a finite positive normal float");
        }
        if (!do_cutoff)
        {
            const double derived_cutoff =
                3.0 * static_cast<double>(cv_sigmas[i]);
            if (!Double_Memory_Is_Finite(&derived_cutoff) ||
                derived_cutoff > std::numeric_limits<float>::max())
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow, "METADYNAMICS::Initial",
                    "Three times CV_sigma is outside the finite float range");
            }
            cutoff[i] = static_cast<float>(derived_cutoff);
        }
        else if (!Float_Memory_Is_Normal(&cutoff[i]) || !(cutoff[i] > 0.0f))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "cutoff must contain finite positive normal floats");
        }
        if (!Float_Memory_Is_Normal(&cutoff[i]))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "METADYNAMICS::Initial",
                "The derived metadynamics cutoff is not a normal float");
        }
        const double inverse_sigma =
            (kde ? 1.414 : 1.0) / static_cast<double>(cv_sigmas[i]);
        if (!Double_Memory_Is_Finite(&inverse_sigma) ||
            inverse_sigma > std::numeric_limits<float>::max())
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "METADYNAMICS::Initial",
                "The inverse CV_sigma is outside the finite float range");
        }
        const float stored_inverse_sigma = static_cast<float>(inverse_sigma);
        if (!Float_Memory_Is_Normal(&stored_inverse_sigma))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "METADYNAMICS::Initial",
                "The inverse CV_sigma is not a positive normal float");
        }
        cv_sigmas[i] = stored_inverse_sigma;
    }
    for (int i = 0; i < ndim; i++)
    {
        sigmas.push_back(cv_sigmas[i]);
        periods.push_back(cv_periods[i]);
    }
    const bool has_potential_input =
        cv_controller->Command_Exist(this->module_name, "potential_in_file");
    const bool has_scatter_input =
        cv_controller->Command_Exist(this->module_name, "scatter_in_file");
    if (has_potential_input && has_scatter_input)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "METADYNAMICS::Initial",
            "potential_in_file and scatter_in_file are mutually exclusive");
    }
    if (has_potential_input)
    {
        read_potential_file_name = cv_controller->Ask_For_String_Parameter(
            this->module_name, "potential_in_file")[0];
        if (read_potential_file_name.empty())
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                "potential_in_file must not be empty");
        }
        if (usegrid || use_scatter)
        {
            Read_Potential(controller);
        }
    }
    else if (has_scatter_input)
    {
        usegrid = false;
        use_scatter = true;
        controller->printf("    Use %d scatter point for CV!\n", scatter_size);
        read_potential_file_name = cv_controller->Ask_For_String_Parameter(
            this->module_name, "scatter_in_file")[0];
        if (read_potential_file_name.empty())
        {
            controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                           "METADYNAMICS::Initial",
                                           "scatter_in_file must not be empty");
        }
        if (usegrid || use_scatter)
        {
            Read_Potential(controller);
        }
    }
    else
    {
        if (cv_controller->Command_Exist(this->module_name, "scatter"))
        {
            scatter_size = *(cv_controller->Ask_For_Int_Parameter(
                this->module_name, "scatter", 1));
            if (scatter_size < 0)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                    "scatter must be a nonnegative point count");
            }
            if (scatter_size > 0)
            {
                usegrid = false;
                use_scatter = true;
                controller->printf("    Use %d scatter point for CV!\n",
                                   scatter_size);
                for (int i = 0; i < cvs.size(); i++)
                {
                    tcoor.push_back(cv_controller->Ask_For_Float_Parameter(
                        cvs[i]->module_name.c_str(), "CV_point", scatter_size,
                        1, false));
                }
            }
            else
            {
                controller->printf("    Not using scatter point for CV\n");
                use_scatter = false;
            }
        }

        cv_mins = cv_controller->Ask_For_Float_Parameter(
            this->module_name, "CV_minimal", cvs.size());
        cv_maxs = cv_controller->Ask_For_Float_Parameter(
            this->module_name, "CV_maximum", cvs.size());
        n_grids = cv_controller->Ask_For_Int_Parameter(this->module_name,
                                                       "CV_grid", cvs.size());
        for (int i = 0; i < cvs.size(); ++i)
        {
            if (!Float_Memory_Is_Zero_Or_Normal(&cv_mins[i]) ||
                !Float_Memory_Is_Zero_Or_Normal(&cv_maxs[i]) ||
                !(cv_maxs[i] > cv_mins[i]))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                    "CV bounds must be finite zero-or-normal floats with "
                    "CV_maximum greater than CV_minimal");
            }
            if (n_grids[i] <= 1)
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorValueErrorCommand, "METADYNAMICS::Initial",
                    "CV_grid should always be greater than 1");
            }
            const double delta =
                (static_cast<double>(cv_maxs[i]) - cv_mins[i]) / n_grids[i];
            if (!Double_Memory_Is_Finite(&delta) || !(delta > 0.0) ||
                delta > std::numeric_limits<float>::max())
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow, "METADYNAMICS::Initial",
                    "A metadynamics CV grid spacing is outside the finite "
                    "positive float range");
            }
            const float stored_delta = static_cast<float>(delta);
            if (!Float_Memory_Is_Normal(&stored_delta))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorOverflow, "METADYNAMICS::Initial",
                    "A metadynamics CV grid spacing is not a positive normal "
                    "float");
            }
            cv_deltas.push_back(stored_delta);
        }
        Set_Grid(controller);
    }
    Malloc_Safely((void**)&Dpotential_local, sizeof(float) * ndim);
    memset(Dpotential_local, 0, sizeof(float) * ndim);
    controller->Step_Print_Initial("meta", "%f");
    controller->Step_Print_Initial("rbias", "%f");
    controller->Step_Print_Initial("rct", "%f");
    controller->printf("    potential output file: %s\n",
                       write_potential_file_name.c_str());
    controller->printf("    edge effect file: %s\n", edge_file_name.c_str());
    is_initialized = 1;
    controller->printf("END INITIALIZING META\n\n");
}

META::Hill::Hill(const Axis& centers, const Axis& inv_w, const Axis& period,
                 const float& theight)
    : centers_(centers), inv_w_(inv_w), periods_(period), height(theight)
{
    const int n = static_cast<int>(centers_.size());
    dx_.resize(n);
    df_.resize(n);
    tder_.resize(n);
}

const META::Gdata& META::Hill::Calc_Hill(const Axis& values)
{
    const int n = static_cast<int>(values.size());
    for (int i = 0; i < n; ++i)
    {
        dx_[i] = Evaluate_Gaussian_Switch(values[i], centers_[i], inv_w_[i],
                                          periods_[i], df_[i]);
    }
    potential = 1.0;
    for (int i = 0; i < n; ++i)
    {
        tder_[i] = 1.0;
        potential *= dx_[i];
        for (int j = 0; j < n; ++j)
        {
            if (j != i)
            {
                tder_[i] *= dx_[j];
            }
            else
            {
                tder_[i] *= df_[j];
            }
        }
    }
    return tder_;
}

static __global__ void Add_Frc(const int atom_numbers, VECTOR* frc,
                               VECTOR* cv_grad, float dheight_dcv)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        frc[i] = frc[i] - dheight_dcv * cv_grad[i];
    }
}

static __global__ void Add_Potential_Kernel(float* d_potential,
                                            const float to_add)
{
    d_potential[0] += to_add;
}

static __global__ void Add_Virial(LTMatrix3* d_virial, const float dU_dCV,
                                  const LTMatrix3* cv_virial)
{
    d_virial[0] = d_virial[0] - dU_dCV * cv_virial[0];
}

void META::Meta_Force_With_Energy_And_Virial(int atom_numbers, VECTOR* frc,
                                             int need_potential,
                                             int need_pressure,
                                             float* d_potential,
                                             LTMatrix3* d_virial)
{
    if (!is_initialized)
    {
        return;
    }
    Potential_And_Derivative(need_potential);
    if (do_borderwall)
    {
        Border_Derivative(border_upper.data(), border_lower.data(), cutoff,
                          Dpotential_local);
    }

    for (int i = 0; i < cvs.size(); ++i)
    {
        Launch_Device_Kernel(Add_Frc, (atom_numbers + 31 / 32), 32, 0, NULL,
                             atom_numbers, frc, cvs[i]->crd_grads,
                             Dpotential_local[i]);
        if (need_pressure)
        {
            Launch_Device_Kernel(Add_Virial, 1, 1, 0, NULL, d_virial,
                                 Dpotential_local[i], cvs[i]->virial);
        }
    }
    if (need_potential)
    {
        Launch_Device_Kernel(Add_Potential_Kernel, 1, 1, 0, NULL, d_potential,
                             potential_local);
    }
}

void META::Potential_And_Derivative(const int need_potential)
{
    if (!is_initialized)
    {
        return;
    }
    for (int i = 0; i < cvs.size(); ++i)
    {
        est_values_[i] = cvs[i]->value;
        Dpotential_local[i] = 0.f;
    }
    Estimate(est_values_, need_potential, true);
}

void META::Border_Derivative(float* border_upper, float* border_lower,
                             float* cutoff, float* Dpotential_local)
{
    for (int i = 0; i < cvs.size(); ++i)
    {
        float h_cv = cvs[i]->value;
        if (h_cv - border_lower[i] < cutoff[i])
        {
            float distance = border_lower[i] - h_cv;
            if (periods[i] > 0)
            {
                distance -= roundf(distance / cv_periods[i]) * cv_periods[i];
            }
            Dpotential_local[i] =
                Dpotential_local[i] - border_potential_height * expf(distance);
        }
        else if (border_upper[i] - h_cv < cutoff[i])
        {
            float distance = h_cv - border_upper[i];
            if (periods[i] > 0)
            {
                distance -= roundf(distance / cv_periods[i]) * cv_periods[i];
            }
            Dpotential_local[i] =
                Dpotential_local[i] + border_potential_height * expf(distance);
        }
    }
}

void META::Do_Metadynamics(int atom_numbers, VECTOR* crd, LTMatrix3 cell,
                           LTMatrix3 rcell, LTMatrix3 reference_cell, int step,
                           int need_potential, int need_pressure, VECTOR* frc,
                           float* d_potential, LTMatrix3* d_virial,
                           float sys_temp, bool commit_history)
{
    if (this->is_initialized)
    {
        temperature = sys_temp;
        Validate_Runtime_Temperature("META::Do_Metadynamics");
        if (ndim <= 0 || cvs.size() != static_cast<std::size_t>(ndim))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "META::Do_Metadynamics",
                "Metadynamics runtime dimensionality is inconsistent: "
                "Ndim=%d, CV count=%zu.",
                ndim, cvs.size());
        }
        const float saved_exit_tag = exit_tag;
        int need = CV_NEED_GPU_VALUE | CV_NEED_CRD_GRADS;
        if (need_pressure)
        {
            need |= CV_NEED_VIRIAL;
        }

        for (int i = 0; i < cvs.size(); i = i + 1)
        {
            if (cvs[i] == NULL)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "META::Do_Metadynamics",
                    "Metadynamics CV dimension %d became null at runtime.", i);
            }
            this->cvs[i]->Compute(atom_numbers, crd, cell, rcell,
                                  reference_cell, need, step);
            if (!Float_Memory_Is_Zero_Or_Normal(&this->cvs[i]->value))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "META::Do_Metadynamics",
                    "Metadynamics CV dimension %d produced a non-finite or "
                    "subnormal value.",
                    i);
            }
        }
        Meta_Force_With_Energy_And_Virial(atom_numbers, frc, need_potential,
                                          need_pressure, d_potential, d_virial);
        if (commit_history)
        {
            Add_Potential(sys_temp, step);
        }
        else
        {
            // Estimate() uses exit_tag as persistent sampling history for the
            // mask path.  Trial/replay force evaluations may update transient
            // diagnostics, but must not advance that history.
            exit_tag = saved_exit_tag;
        }
    }
}
#endif  // SPONGE_META_ATOMIC_OUTPUT_PROBE
