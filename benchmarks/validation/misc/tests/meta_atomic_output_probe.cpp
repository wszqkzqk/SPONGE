#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>

#define SPONGE_META_ATOMIC_OUTPUT_PROBE
#include "bias/sinkmeta.cpp"

int CONTROLLER::MPI_rank = 0;

namespace
{

std::atomic<bool> post_rename_entered{false};
std::atomic<bool> fatal_callback_done{false};
std::atomic<bool> normal_cleanup_entered{false};
std::atomic<bool> competing_cleanup_observed{false};

void Wait_For_Fatal_Callback(Meta_Atomic_Output*)
{
    post_rename_entered.store(true, std::memory_order_release);
    while (!fatal_callback_done.load(std::memory_order_acquire))
    {
        std::this_thread::yield();
    }
}

void Wait_For_Competing_Cleanup(Meta_Atomic_Output*)
{
    normal_cleanup_entered.store(true, std::memory_order_release);
    while (!competing_cleanup_observed.load(std::memory_order_acquire))
    {
        std::this_thread::yield();
    }
}

void Mark_Competing_Cleanup_Observed(Meta_Atomic_Output*)
{
    competing_cleanup_observed.store(true, std::memory_order_release);
}

int Test_Commit_Interleaving(const char* path)
{
    std::remove(path);
    CONTROLLER controller{};
    Meta_Atomic_Output output(&controller, path, "meta_atomic_output_probe");
    if (std::fputs("complete output\n", output.Stream()) < 0)
    {
        return EXIT_FAILURE;
    }

    post_rename_entered.store(false, std::memory_order_relaxed);
    fatal_callback_done.store(false, std::memory_order_relaxed);
    Meta_Atomic_Output::Set_Probe_Post_Rename_Hook(Wait_For_Fatal_Callback);
    std::thread fatal_thread(
        [&]()
        {
            while (!post_rename_entered.load(std::memory_order_acquire))
            {
                std::this_thread::yield();
            }
            output.Probe_Invoke_Fatal_Cleanup();
            fatal_callback_done.store(true, std::memory_order_release);
        });
    output.Commit();
    fatal_thread.join();
    Meta_Atomic_Output::Set_Probe_Post_Rename_Hook(NULL);

    FILE* input = std::fopen(path, "rb");
    if (input == NULL) return EXIT_FAILURE;
    char contents[64] = {};
    const std::size_t size = std::fread(contents, 1, sizeof(contents), input);
    std::fclose(input);
    return size == std::strlen("complete output\n") &&
                   std::memcmp(contents, "complete output\n", size) == 0
               ? EXIT_SUCCESS
               : EXIT_FAILURE;
}

int Test_Close_Interleaving(const char* path)
{
    constexpr int repetitions = 300;
    CONTROLLER controller{};
    for (int repetition = 0; repetition < repetitions; ++repetition)
    {
        std::remove(path);
        Meta_Atomic_Output output(&controller, path,
                                  "meta_atomic_output_probe");
        if (std::fputs("temporary output\n", output.Stream()) < 0)
        {
            return EXIT_FAILURE;
        }
        std::atomic<bool> start{false};
        std::thread normal_thread(
            [&]()
            {
                while (!start.load(std::memory_order_acquire))
                {
                    std::this_thread::yield();
                }
                output.Probe_Close();
            });
        std::thread fatal_thread(
            [&]()
            {
                while (!start.load(std::memory_order_acquire))
                {
                    std::this_thread::yield();
                }
                output.Probe_Invoke_Fatal_Cleanup();
            });
        start.store(true, std::memory_order_release);
        normal_thread.join();
        fatal_thread.join();
    }
    return EXIT_SUCCESS;
}

int Test_Concurrent_Normal_Close(const char* path, bool fatal_already_claimed)
{
    constexpr int repetitions = 300;
    CONTROLLER controller{};
    for (int repetition = 0; repetition < repetitions; ++repetition)
    {
        std::remove(path);
        Meta_Atomic_Output output(&controller, path,
                                  "meta_atomic_output_probe");
        if (std::fputs("temporary output\n", output.Stream()) < 0)
        {
            return EXIT_FAILURE;
        }
        if (fatal_already_claimed) output.Probe_Invoke_Fatal_Cleanup();

        normal_cleanup_entered.store(false, std::memory_order_relaxed);
        competing_cleanup_observed.store(false, std::memory_order_relaxed);
        Meta_Atomic_Output::Set_Probe_Normal_Cleanup_Hook(
            Wait_For_Competing_Cleanup);
        Meta_Atomic_Output::Set_Probe_Competing_Cleanup_Hook(
            Mark_Competing_Cleanup_Observed);
        std::thread owner_thread([&]() { output.Probe_Close(); });
        std::thread competing_thread(
            [&]()
            {
                while (!normal_cleanup_entered.load(std::memory_order_acquire))
                {
                    std::this_thread::yield();
                }
                output.Probe_Close();
            });
        owner_thread.join();
        competing_thread.join();
        Meta_Atomic_Output::Set_Probe_Normal_Cleanup_Hook(NULL);
        Meta_Atomic_Output::Set_Probe_Competing_Cleanup_Hook(NULL);
    }
    return EXIT_SUCCESS;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 3) return EXIT_FAILURE;
    if (std::strcmp(argv[1], "commit_fatal") == 0)
    {
        return Test_Commit_Interleaving(argv[2]);
    }
    if (std::strcmp(argv[1], "close_fatal") == 0)
    {
        return Test_Close_Interleaving(argv[2]);
    }
    if (std::strcmp(argv[1], "close_close") == 0)
    {
        return Test_Concurrent_Normal_Close(argv[2], false);
    }
    if (std::strcmp(argv[1], "claimed_close_close") == 0)
    {
        return Test_Concurrent_Normal_Close(argv[2], true);
    }
    return EXIT_FAILURE;
}
