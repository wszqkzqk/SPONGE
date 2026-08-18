#pragma once

#include <cerrno>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <string>
#include <thread>
#include <utility>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

namespace SpongeH5MD
{
inline bool Remove_File_If_Exists(const std::string& path,
                                  std::string* error_message = nullptr)
{
    if (std::remove(path.c_str()) == 0 || errno == ENOENT) return true;
    if (error_message != nullptr)
    {
        *error_message = "failed to remove temporary file " + path + ": " +
                         std::strerror(errno);
    }
    return false;
}

inline bool Atomic_Replace_File(const std::string& temporary_path,
                                const std::string& destination_path,
                                std::string* error_message = nullptr)
{
#ifdef _WIN32
    constexpr int kSharingRetryCount = 50;
    constexpr auto kSharingRetryDelay = std::chrono::milliseconds(10);
    DWORD last_error = ERROR_SUCCESS;
    for (int attempt = 0; attempt <= kSharingRetryCount; ++attempt)
    {
        if (MoveFileExA(temporary_path.c_str(), destination_path.c_str(),
                        MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH))
        {
            return true;
        }
        last_error = GetLastError();
        const bool retryable = last_error == ERROR_SHARING_VIOLATION ||
                               last_error == ERROR_LOCK_VIOLATION ||
                               last_error == ERROR_ACCESS_DENIED;
        if (!retryable || attempt == kSharingRetryCount) break;
        std::this_thread::sleep_for(kSharingRetryDelay);
    }
    if (error_message != nullptr)
    {
        *error_message = "failed to atomically replace " + destination_path +
                         " with " + temporary_path + ": Windows error " +
                         std::to_string(last_error);
    }
    return false;
#else
    if (std::rename(temporary_path.c_str(), destination_path.c_str()) == 0)
    {
        return true;
    }
    if (error_message != nullptr)
    {
        *error_message = "failed to atomically replace " + destination_path +
                         " with " + temporary_path + ": " +
                         std::strerror(errno);
    }
    return false;
#endif
}

class TemporaryFileGuard
{
   public:
    explicit TemporaryFileGuard(std::string path) : path_(std::move(path)) {}
    ~TemporaryFileGuard()
    {
        if (armed_) Remove_File_If_Exists(path_);
    }

    TemporaryFileGuard(const TemporaryFileGuard&) = delete;
    TemporaryFileGuard& operator=(const TemporaryFileGuard&) = delete;

    void Release() { armed_ = false; }

   private:
    std::string path_;
    bool armed_ = true;
};
}  // namespace SpongeH5MD
