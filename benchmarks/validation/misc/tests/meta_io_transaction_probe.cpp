#include <fcntl.h>
#include <unistd.h>

#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <new>
#include <string>
#include <thread>
#include <vector>

#include "bias/meta_io_transaction.h"

namespace
{

struct Allocation_Control
{
    int allocation_count = 0;
    int fail_on = 0;
};

Allocation_Control allocation_control;

template <typename T>
struct Failing_Allocator
{
    using value_type = T;

    Failing_Allocator() noexcept = default;

    template <typename U>
    Failing_Allocator(const Failing_Allocator<U>&) noexcept
    {
    }

    T* allocate(std::size_t count)
    {
        ++allocation_control.allocation_count;
        if (allocation_control.fail_on != 0 &&
            allocation_control.allocation_count == allocation_control.fail_on)
        {
            throw std::bad_alloc();
        }
        return static_cast<T*>(::operator new(count * sizeof(T)));
    }

    void deallocate(T* pointer, std::size_t) noexcept
    {
        ::operator delete(pointer);
    }

    template <typename U>
    bool operator==(const Failing_Allocator<U>&) const noexcept
    {
        return true;
    }

    template <typename U>
    bool operator!=(const Failing_Allocator<U>&) const noexcept
    {
        return false;
    }
};

struct Probe_Hill
{
    int value = 0;
    Probe_Hill() = default;
    explicit Probe_Hill(int input) : value(input) {}
    Probe_Hill(const Probe_Hill&) = default;
    Probe_Hill& operator=(const Probe_Hill&) = default;
    Probe_Hill(Probe_Hill&&) noexcept = default;
    Probe_Hill& operator=(Probe_Hill&&) noexcept = default;
};

int Test_Second_Reserve_Failure()
{
    using Hill_Vector = std::vector<Probe_Hill, Failing_Allocator<Probe_Hill>>;
    using Sink_Vector = std::vector<float, Failing_Allocator<float>>;

    allocation_control = {};
    Hill_Vector hills;
    Sink_Vector sink;
    Hill_Vector parsed_hills;
    Sink_Vector parsed_sink;
    hills.reserve(1);
    sink.reserve(1);
    parsed_hills.reserve(1);
    parsed_sink.reserve(1);
    hills.emplace_back(11);
    sink.emplace_back(12.0f);
    parsed_hills.emplace_back(21);
    parsed_sink.emplace_back(22.0f);

    allocation_control.allocation_count = 0;
    allocation_control.fail_on = 2;
    bool failed = false;
    try
    {
        sponge_meta_io::Publish_History(&hills, &sink, &parsed_hills,
                                        parsed_sink);
    }
    catch (const std::bad_alloc&)
    {
        failed = true;
    }
    allocation_control.fail_on = 0;

    if (!failed || hills.size() != 1 || sink.size() != 1 ||
        hills[0].value != 11 || sink[0] != 12.0f)
    {
        std::fprintf(stderr,
                     "second reserve failure published partial history\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool Lock_File(int descriptor, short type)
{
    struct flock lock = {};
    lock.l_type = type;
    lock.l_whence = SEEK_SET;
    for (;;)
    {
        if (fcntl(descriptor, F_SETLKW, &lock) == 0) return true;
        if (errno != EINTR) return false;
    }
}

bool Write_All(int descriptor, const char* data, std::size_t size)
{
    std::size_t completed = 0;
    while (completed < size)
    {
        const ssize_t written =
            write(descriptor, data + completed, size - completed);
        if (written > 0)
        {
            completed += static_cast<std::size_t>(written);
            continue;
        }
        if (written < 0 && errno == EINTR) continue;
        return false;
    }
    return true;
}

bool Snapshot_Is_Complete(const char* path)
{
    const std::lock_guard<std::mutex> process_lock(
        sponge_meta_io::Process_Mutex());
    const int descriptor = open(path, O_RDONLY | O_CLOEXEC);
    if (descriptor < 0) return errno == ENOENT;
    if (!Lock_File(descriptor, F_RDLCK))
    {
        close(descriptor);
        return false;
    }
    const off_t size = lseek(descriptor, 0, SEEK_END);
    std::string contents;
    if (size >= 0) contents.resize(static_cast<std::size_t>(size));
    std::size_t completed = 0;
    while (size >= 0 && completed < contents.size())
    {
        const ssize_t count = pread(descriptor, &contents[completed],
                                    contents.size() - completed, completed);
        if (count > 0)
        {
            completed += static_cast<std::size_t>(count);
            continue;
        }
        if (count < 0 && errno == EINTR) continue;
        completed = contents.size() + 1;
        break;
    }
    Lock_File(descriptor, F_UNLCK);
    close(descriptor);
    if (size < 0 || completed != contents.size()) return false;
    if (!contents.empty() && contents.back() != '\n') return false;
    std::size_t cursor = 0;
    while (cursor < contents.size())
    {
        const std::size_t newline = contents.find('\n', cursor);
        if (newline == std::string::npos || newline - cursor != 7 ||
            contents.compare(cursor, 7, "payload") != 0)
        {
            return false;
        }
        cursor = newline + 1;
    }
    return true;
}

bool Append_Record(const char* path)
{
    const std::lock_guard<std::mutex> process_lock(
        sponge_meta_io::Process_Mutex());
    const int descriptor =
        open(path, O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC, 0666);
    if (descriptor < 0) return false;
    if (!Lock_File(descriptor, F_WRLCK))
    {
        close(descriptor);
        return false;
    }
    const bool payload_ok = Write_All(descriptor, "payload", 7);
    std::this_thread::yield();
    const bool marker_ok = payload_ok && Write_All(descriptor, "\n", 1);
    Lock_File(descriptor, F_UNLCK);
    close(descriptor);
    return marker_ok;
}

int Test_Threaded_Append_And_Load(const char* path)
{
    unlink(path);
    constexpr int writer_count = 4;
    constexpr int reader_count = 4;
    constexpr int records_per_writer = 200;
    std::atomic<bool> start{false};
    std::atomic<bool> valid{true};
    std::vector<std::thread> threads;
    for (int writer = 0; writer < writer_count; ++writer)
    {
        threads.emplace_back(
            [&]()
            {
                while (!start.load(std::memory_order_acquire))
                {
                }
                for (int record = 0; record < records_per_writer; ++record)
                {
                    if (!Append_Record(path))
                    {
                        valid.store(false, std::memory_order_release);
                        return;
                    }
                }
            });
    }
    for (int reader = 0; reader < reader_count; ++reader)
    {
        threads.emplace_back(
            [&]()
            {
                while (!start.load(std::memory_order_acquire))
                {
                }
                for (int snapshot = 0; snapshot < records_per_writer;
                     ++snapshot)
                {
                    if (!Snapshot_Is_Complete(path))
                    {
                        valid.store(false, std::memory_order_release);
                        return;
                    }
                }
            });
    }
    start.store(true, std::memory_order_release);
    for (std::thread& thread : threads) thread.join();
    if (!valid.load(std::memory_order_acquire) || !Snapshot_Is_Complete(path))
    {
        std::fprintf(stderr,
                     "threaded append/load observed a partial record\n");
        return EXIT_FAILURE;
    }

    FILE* input = std::fopen(path, "rb");
    if (input == NULL) return EXIT_FAILURE;
    int newline_count = 0;
    for (int byte = std::fgetc(input); byte != EOF; byte = std::fgetc(input))
    {
        if (byte == '\n') ++newline_count;
    }
    std::fclose(input);
    return newline_count == writer_count * records_per_writer ? EXIT_SUCCESS
                                                              : EXIT_FAILURE;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], "publish_failure") == 0)
    {
        return Test_Second_Reserve_Failure();
    }
    if (argc == 3 && std::strcmp(argv[1], "thread_io") == 0)
    {
        return Test_Threaded_Append_And_Load(argv[2]);
    }
    return EXIT_FAILURE;
}
