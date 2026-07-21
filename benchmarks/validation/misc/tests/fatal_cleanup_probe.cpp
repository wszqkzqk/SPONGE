#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>
#include <thread>
#include <vector>

#include "control.h"

int CONTROLLER::MPI_rank = 0;

namespace
{

bool reject_cpp_allocations = false;

}  // namespace

void* operator new(std::size_t size)
{
    if (reject_cpp_allocations) throw std::bad_alloc();
    void* allocation = std::malloc(size == 0 ? 1 : size);
    if (allocation == NULL) throw std::bad_alloc();
    return allocation;
}

void* operator new[](std::size_t size) { return ::operator new(size); }

void operator delete(void* allocation) noexcept { std::free(allocation); }

void operator delete[](void* allocation) noexcept { std::free(allocation); }

void operator delete(void* allocation, std::size_t) noexcept
{
    std::free(allocation);
}

void operator delete[](void* allocation, std::size_t) noexcept
{
    std::free(allocation);
}

namespace
{

struct Marker_Context
{
    const char* path = NULL;
    char marker = '\0';
    CONTROLLER* reentry_controller = NULL;
};

struct Callback_Operations_Context
{
    CONTROLLER* controller = NULL;
    Marker_Context* older = NULL;
    const char* path = NULL;
};

struct Fatal_Race_Context
{
    CONTROLLER* controller = NULL;
    const char* path = NULL;
    std::atomic<bool> target_started{false};
    std::atomic<bool> registration_finished{false};
    std::atomic<bool> unregister_entered{false};
    std::atomic<bool> unregister_returned{false};
    std::atomic<bool> registration_succeeded{false};
    std::atomic<bool> unregistration_succeeded{false};
};

struct Dual_Fatal_Context
{
    const char* path = NULL;
    std::atomic<bool> owner_cleanup_started{false};
    std::atomic<bool> secondary_calling{false};
};

struct Diagnostic_Source_Context
{
    Marker_Context marker;
    char* error_by = NULL;
    char* extra = NULL;
    char* formatted = NULL;
    bool release_sources = false;
};

const char* atexit_marker_path = NULL;

void Append_Marker(void* opaque_context)
{
    Marker_Context* context = static_cast<Marker_Context*>(opaque_context);
    FILE* output = std::fopen(context->path, "ab");
    if (output != NULL)
    {
        std::fputc(static_cast<unsigned char>(context->marker), output);
        std::fclose(output);
    }
    if (context->reentry_controller != NULL)
    {
        context->reentry_controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "fatal cleanup reentry");
    }
}

void Ignore_Context(void*) {}

void Append_Then_Throw(void* opaque_context)
{
    Append_Marker(opaque_context);
    throw 7;
}

void Exercise_Callback_Registry_Operations(void* opaque_context)
{
    Callback_Operations_Context* context =
        static_cast<Callback_Operations_Context*>(opaque_context);
    const bool self_unregistered =
        context->controller->Unregister_Fatal_Cleanup(
            Exercise_Callback_Registry_Operations, context);
    const bool older_unregistered =
        context->controller->Unregister_Fatal_Cleanup(Append_Marker,
                                                      context->older);
    const bool registration_accepted =
        context->controller->Register_Fatal_Cleanup(Ignore_Context, context);
    Marker_Context result{
        context->path,
        !self_unregistered && older_unregistered && !registration_accepted
            ? 'O'
            : 'E',
        NULL};
    Append_Marker(&result);
}

void Run_Race_Target(void* opaque_context)
{
    Fatal_Race_Context* context =
        static_cast<Fatal_Race_Context*>(opaque_context);
    Marker_Context started{context->path, 'T', NULL};
    Append_Marker(&started);
    context->target_started.store(true, std::memory_order_release);
    while (!context->registration_finished.load(std::memory_order_acquire) ||
           !context->unregister_entered.load(std::memory_order_acquire))
    {
        std::this_thread::yield();
    }
    if (context->unregister_returned.load(std::memory_order_acquire))
    {
        Marker_Context premature{context->path, 'E', NULL};
        Append_Marker(&premature);
    }
}

void Run_Race_Gate(void* opaque_context)
{
    Fatal_Race_Context* context =
        static_cast<Fatal_Race_Context*>(opaque_context);
    while (!context->unregister_returned.load(std::memory_order_acquire))
    {
        std::this_thread::yield();
    }
    const bool valid =
        !context->registration_succeeded.load(std::memory_order_acquire) &&
        context->unregistration_succeeded.load(std::memory_order_acquire);
    Marker_Context result{context->path, valid ? 'G' : 'E', NULL};
    Append_Marker(&result);
}

void Hold_First_Fatal(void* opaque_context)
{
    Dual_Fatal_Context* context =
        static_cast<Dual_Fatal_Context*>(opaque_context);
    Marker_Context marker{context->path, 'T', NULL};
    Append_Marker(&marker);
    context->owner_cleanup_started.store(true, std::memory_order_release);
    while (!context->secondary_calling.load(std::memory_order_acquire))
    {
        std::this_thread::yield();
    }
    while (CONTROLLER::Fatal_Sequence_Waiter_Count_For_Test() == 0)
    {
        std::this_thread::yield();
    }
}

void Append_Atexit_Marker()
{
    if (atexit_marker_path == NULL) return;
    Marker_Context marker{atexit_marker_path, 'X', NULL};
    Append_Marker(&marker);
}

void Destroy_Diagnostic_Sources(void* opaque_context)
{
    Diagnostic_Source_Context* context =
        static_cast<Diagnostic_Source_Context*>(opaque_context);
    Append_Marker(&context->marker);
    if (context->release_sources)
    {
        std::free(context->error_by);
        std::free(context->extra);
        std::free(context->formatted);
        context->error_by = NULL;
        context->extra = NULL;
        context->formatted = NULL;
        return;
    }
    std::strcpy(context->error_by, "DESTROYED_ERROR_BY");
    std::strcpy(context->extra, "DESTROYED_EXTRA");
    std::strcpy(context->formatted, "DESTROYED_FORMATTED");
}

char* Duplicate_Diagnostic_Source(const char* source)
{
    const std::size_t length = std::strlen(source);
    char* duplicate = static_cast<char*>(std::malloc(length + 1));
    if (duplicate != NULL) std::memcpy(duplicate, source, length + 1);
    return duplicate;
}

int Fail(const char* message)
{
    std::fprintf(stderr, "probe setup failed: %s\n", message);
    return EXIT_FAILURE;
}

[[noreturn]] void Throw_Fatal(CONTROLLER* controller)
{
    controller->Throw_SPONGE_Error(
        spongeErrorValueErrorCommand, "fatal_cleanup_probe",
        "The probe intentionally enters the fatal termination path.");
    std::abort();
}

int Test_Fatal(const char* marker_path)
{
    CONTROLLER registering_controller{};
    CONTROLLER throwing_controller{};
    Marker_Context marker{marker_path, 'F', NULL};
    if (!registering_controller.Register_Fatal_Cleanup(Append_Marker, &marker))
    {
        return Fail("fatal callback registration was rejected");
    }
    if (registering_controller.Register_Fatal_Cleanup(Append_Marker, &marker))
    {
        return Fail("duplicate callback/context registration was accepted");
    }
    Throw_Fatal(&throwing_controller);
}

int Test_Unregister(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context marker{marker_path, 'U', NULL};
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &marker))
    {
        return Fail("callback registration was rejected");
    }
    if (!controller.Unregister_Fatal_Cleanup(Append_Marker, &marker))
    {
        return Fail("registered callback was not found");
    }
    if (controller.Unregister_Fatal_Cleanup(Append_Marker, &marker))
    {
        return Fail("callback was unregistered more than once");
    }
    Throw_Fatal(&controller);
}

int Test_Reentrant_Lifo(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context first{marker_path, 'A', NULL};
    Marker_Context second{marker_path, 'B', &controller};
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &first) ||
        !controller.Register_Fatal_Cleanup(Append_Marker, &second))
    {
        return Fail("LIFO callback registration was rejected");
    }
    Throw_Fatal(&controller);
}

int Test_Unregister_Preserves_Order(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context first{marker_path, 'A', NULL};
    Marker_Context removed{marker_path, 'B', NULL};
    Marker_Context third{marker_path, 'C', NULL};
    Marker_Context newest{marker_path, 'D', NULL};
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &first) ||
        !controller.Register_Fatal_Cleanup(Append_Marker, &removed) ||
        !controller.Register_Fatal_Cleanup(Append_Marker, &third))
    {
        return Fail("ordered callback registration was rejected");
    }
    if (!controller.Unregister_Fatal_Cleanup(Append_Marker, &removed))
    {
        return Fail("middle callback was not found");
    }
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &newest))
    {
        return Fail("registration after unregistration was rejected");
    }
    Throw_Fatal(&controller);
}

int Test_Capacity(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context contexts[CONTROLLER::FATAL_CLEANUP_CAPACITY] = {};
    if (controller.Register_Fatal_Cleanup(NULL, NULL))
    {
        return Fail("null callback registration was accepted");
    }
    for (std::size_t i = 0; i < CONTROLLER::FATAL_CLEANUP_CAPACITY; ++i)
    {
        contexts[i] = {marker_path, 'X', NULL};
        if (!controller.Register_Fatal_Cleanup(Append_Marker, &contexts[i]))
        {
            return Fail(
                "registry reported full before its documented capacity");
        }
    }
    Marker_Context overflow{marker_path, 'Y', NULL};
    if (controller.Register_Fatal_Cleanup(Append_Marker, &overflow))
    {
        return Fail("full registry silently accepted another callback");
    }
    const std::size_t middle = CONTROLLER::FATAL_CLEANUP_CAPACITY / 2;
    if (!controller.Unregister_Fatal_Cleanup(Append_Marker, &contexts[middle]))
    {
        return Fail("callback could not be removed from a full registry");
    }
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &overflow))
    {
        return Fail("a freed registry slot could not be reused");
    }
    Marker_Context still_overflowing{marker_path, 'Z', NULL};
    if (controller.Register_Fatal_Cleanup(Append_Marker, &still_overflowing))
    {
        return Fail("re-filled registry silently overwrote a callback");
    }
    Throw_Fatal(&controller);
}

int Test_Concurrent_Registry(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context shared{marker_path, 'C', NULL};
    constexpr int thread_count = 8;
    constexpr int repetitions = 200;
    for (int repetition = 0; repetition < repetitions; ++repetition)
    {
        std::atomic<bool> start{false};
        std::atomic<int> accepted{0};
        std::vector<std::thread> workers;
        workers.reserve(thread_count);
        for (int i = 0; i < thread_count; ++i)
        {
            workers.emplace_back(
                [&]()
                {
                    while (!start.load(std::memory_order_acquire))
                    {
                    }
                    if (controller.Register_Fatal_Cleanup(Ignore_Context,
                                                          &shared))
                    {
                        accepted.fetch_add(1, std::memory_order_relaxed);
                    }
                });
        }
        start.store(true, std::memory_order_release);
        for (std::thread& worker : workers) worker.join();
        if (accepted.load(std::memory_order_relaxed) != 1)
        {
            return Fail("concurrent duplicate registration was not atomic");
        }
        if (!controller.Unregister_Fatal_Cleanup(Ignore_Context, &shared))
        {
            return Fail("concurrently registered callback was not present");
        }
    }
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &shared))
    {
        return Fail("registry was corrupted by concurrent registration");
    }
    Throw_Fatal(&controller);
}

int Test_Concurrent_Fatal(const char* marker_path)
{
    CONTROLLER controller{};
    controller.mdinfo = std::tmpfile();
    if (controller.mdinfo == NULL)
        return Fail("could not create non-null mdinfo");

    Dual_Fatal_Context context{marker_path};
    if (!controller.Register_Fatal_Cleanup(Hold_First_Fatal, &context))
    {
        return Fail("concurrent-fatal callback registration was rejected");
    }

    std::thread owner(
        [&]()
        {
            controller.Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                          "fatal sequence owner");
        });
    while (!context.owner_cleanup_started.load(std::memory_order_acquire))
    {
        std::this_thread::yield();
    }

    std::thread secondary(
        [&]()
        {
            context.secondary_calling.store(true, std::memory_order_release);
            controller.Throw_SPONGE_Error(spongeErrorMissingCommand,
                                          "secondary fatal caller");
        });
    owner.join();
    secondary.join();
    return Fail("concurrent fatal paths unexpectedly returned");
}

int Test_Callback_Operations(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context older{marker_path, 'A', NULL};
    Callback_Operations_Context operations{&controller, &older, marker_path};
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &older) ||
        !controller.Register_Fatal_Cleanup(
            Exercise_Callback_Registry_Operations, &operations))
    {
        return Fail("callback-operation registration was rejected");
    }
    Throw_Fatal(&controller);
}

int Test_Callback_Throw(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context older{marker_path, 'A', NULL};
    Marker_Context throwing{marker_path, 'B', NULL};
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &older) ||
        !controller.Register_Fatal_Cleanup(Append_Then_Throw, &throwing))
    {
        return Fail("throwing-callback registration was rejected");
    }
    Throw_Fatal(&controller);
}

int Test_Register_Unregister_Fatal_Race(const char* marker_path)
{
    CONTROLLER controller{};
    Fatal_Race_Context context{};
    context.controller = &controller;
    context.path = marker_path;
    if (!controller.Register_Fatal_Cleanup(Run_Race_Gate, &context) ||
        !controller.Register_Fatal_Cleanup(Run_Race_Target, &context))
    {
        return Fail("race callback registration was rejected");
    }

    std::thread registrar(
        [&]()
        {
            while (!context.target_started.load(std::memory_order_acquire))
            {
                std::this_thread::yield();
            }
            const bool accepted =
                controller.Register_Fatal_Cleanup(Ignore_Context, &context);
            context.registration_succeeded.store(accepted,
                                                 std::memory_order_release);
            context.registration_finished.store(true,
                                                std::memory_order_release);
        });
    std::thread unregistrar(
        [&]()
        {
            while (!context.target_started.load(std::memory_order_acquire))
            {
                std::this_thread::yield();
            }
            context.unregister_entered.store(true, std::memory_order_release);
            const bool removed =
                controller.Unregister_Fatal_Cleanup(Run_Race_Target, &context);
            context.unregistration_succeeded.store(removed,
                                                   std::memory_order_release);
            context.unregister_returned.store(true, std::memory_order_release);
        });

    Throw_Fatal(&controller);
    registrar.join();
    unregistrar.join();
    return Fail("fatal race unexpectedly returned");
}

int Test_Hard_Exit(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context cleanup{marker_path, 'A', NULL};
    atexit_marker_path = marker_path;
    if (std::atexit(Append_Atexit_Marker) != 0)
    {
        return Fail("atexit registration failed");
    }
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &cleanup))
    {
        return Fail("hard-exit callback registration was rejected");
    }
    Throw_Fatal(&controller);
}

int Test_Diagnostic_Allocation_Failure(const char* marker_path, bool formatted)
{
    CONTROLLER controller{};
    Marker_Context marker{marker_path, formatted ? 'M' : 'D', NULL};
    reject_cpp_allocations = true;
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &marker))
    {
        reject_cpp_allocations = false;
        return Fail("allocation-failure callback registration was rejected");
    }
    if (formatted)
    {
        controller.Throw_Formatted_SPONGE_Error(
            spongeErrorMallocFailed, "fatal_cleanup_probe", "%s",
            "The diagnostic allocation is intentionally rejected.");
    }
    else
    {
        controller.Throw_SPONGE_Error(spongeErrorMallocFailed,
                                      "fatal_cleanup_probe");
    }
    reject_cpp_allocations = false;
    return Fail("allocation-failure fatal path unexpectedly returned");
}

int Test_Diagnostic_Snapshot(const char* marker_path, bool formatted,
                             bool release_sources, bool allocation_error)
{
    CONTROLLER controller{};
    char error_by_storage[64] = "ORIGINAL_ERROR_BY";
    char extra_storage[64] = "ORIGINAL_EXTRA";
    char formatted_storage[64] = "ORIGINAL_FORMATTED";
    Diagnostic_Source_Context context{};
    context.marker = {marker_path, release_sources ? 'F' : 'S', NULL};
    context.release_sources = release_sources;
    if (release_sources)
    {
        context.error_by = Duplicate_Diagnostic_Source(error_by_storage);
        context.extra = Duplicate_Diagnostic_Source(extra_storage);
        context.formatted = Duplicate_Diagnostic_Source(formatted_storage);
        if (context.error_by == NULL || context.extra == NULL ||
            context.formatted == NULL)
        {
            std::free(context.error_by);
            std::free(context.extra);
            std::free(context.formatted);
            return Fail("could not allocate diagnostic sources");
        }
    }
    else
    {
        context.error_by = error_by_storage;
        context.extra = extra_storage;
        context.formatted = formatted_storage;
    }
    if (!controller.Register_Fatal_Cleanup(Destroy_Diagnostic_Sources,
                                           &context))
    {
        if (release_sources)
        {
            std::free(context.error_by);
            std::free(context.extra);
            std::free(context.formatted);
        }
        return Fail("diagnostic-source callback registration was rejected");
    }
    if (formatted)
    {
        controller.Throw_Formatted_SPONGE_Error(
            allocation_error ? spongeErrorMallocFailed
                             : spongeErrorValueErrorCommand,
            context.error_by, "Formatted source: %s", context.formatted);
    }
    else
    {
        controller.Throw_SPONGE_Error(allocation_error
                                          ? spongeErrorMallocFailed
                                          : spongeErrorValueErrorCommand,
                                      context.error_by, context.extra);
    }
    return Fail("diagnostic-snapshot fatal path unexpectedly returned");
}

int Test_Diagnostic_Capture_Allocation_Failure(const char* marker_path)
{
    CONTROLLER controller{};
    char error_by[64] = "UNAVAILABLE_ERROR_BY";
    char extra[64] = "UNAVAILABLE_EXTRA";
    char formatted[64] = "UNUSED_FORMATTED";
    Diagnostic_Source_Context context{
        {marker_path, 'C', NULL}, error_by, extra, formatted, false};
    if (!controller.Register_Fatal_Cleanup(Destroy_Diagnostic_Sources,
                                           &context))
    {
        return Fail("capture-failure callback registration was rejected");
    }
    reject_cpp_allocations = true;
    controller.Throw_SPONGE_Error(spongeErrorValueErrorCommand, error_by,
                                  extra);
    reject_cpp_allocations = false;
    return Fail("diagnostic capture failure unexpectedly returned");
}

int Test_Diagnostic_Capture_Truncation(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context marker{marker_path, 'T', NULL};
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &marker))
    {
        return Fail("truncated-diagnostic callback registration was rejected");
    }
    char extra[16384];
    std::memset(extra, 'q', sizeof(extra));
    constexpr char ending[] = "FALLBACK_TAIL_SHOULD_BE_TRUNCATED";
    std::memcpy(extra + sizeof(extra) - sizeof(ending), ending, sizeof(ending));
    reject_cpp_allocations = true;
    controller.Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                  "TRUNCATED_ERROR_BY", extra);
    reject_cpp_allocations = false;
    return Fail("truncated diagnostic fatal path unexpectedly returned");
}

int Test_Long_Formatted_Diagnostic(const char* marker_path)
{
    CONTROLLER controller{};
    Marker_Context marker{marker_path, 'L', NULL};
    if (!controller.Register_Fatal_Cleanup(Append_Marker, &marker))
    {
        return Fail("long-diagnostic callback registration was rejected");
    }
    char diagnostic[4097];
    std::memset(diagnostic, 'x', sizeof(diagnostic));
    constexpr char ending[] = "LONG_DIAGNOSTIC_END";
    const std::size_t ending_offset = sizeof(diagnostic) - sizeof(ending);
    std::memcpy(diagnostic + ending_offset, ending, sizeof(ending));
    controller.Throw_Formatted_SPONGE_Error(
        spongeErrorValueErrorCommand, "fatal_cleanup_probe", "%s", diagnostic);
    return Fail("long formatted fatal path unexpectedly returned");
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 3) return Fail("expected MODE and MARKER_PATH arguments");
    if (std::strcmp(argv[1], "fatal") == 0) return Test_Fatal(argv[2]);
    if (std::strcmp(argv[1], "unregister") == 0)
    {
        return Test_Unregister(argv[2]);
    }
    if (std::strcmp(argv[1], "reentrant_lifo") == 0)
    {
        return Test_Reentrant_Lifo(argv[2]);
    }
    if (std::strcmp(argv[1], "preserve_order") == 0)
    {
        return Test_Unregister_Preserves_Order(argv[2]);
    }
    if (std::strcmp(argv[1], "capacity") == 0)
    {
        return Test_Capacity(argv[2]);
    }
    if (std::strcmp(argv[1], "concurrent_registry") == 0)
    {
        return Test_Concurrent_Registry(argv[2]);
    }
    if (std::strcmp(argv[1], "concurrent_fatal") == 0)
    {
        return Test_Concurrent_Fatal(argv[2]);
    }
    if (std::strcmp(argv[1], "callback_operations") == 0)
    {
        return Test_Callback_Operations(argv[2]);
    }
    if (std::strcmp(argv[1], "callback_throw") == 0)
    {
        return Test_Callback_Throw(argv[2]);
    }
    if (std::strcmp(argv[1], "fatal_registry_race") == 0)
    {
        return Test_Register_Unregister_Fatal_Race(argv[2]);
    }
    if (std::strcmp(argv[1], "hard_exit") == 0)
    {
        return Test_Hard_Exit(argv[2]);
    }
    if (std::strcmp(argv[1], "direct_allocation_failure") == 0)
    {
        return Test_Diagnostic_Allocation_Failure(argv[2], false);
    }
    if (std::strcmp(argv[1], "formatted_allocation_failure") == 0)
    {
        return Test_Diagnostic_Allocation_Failure(argv[2], true);
    }
    if (std::strcmp(argv[1], "direct_snapshot_mutation") == 0)
    {
        return Test_Diagnostic_Snapshot(argv[2], false, false, false);
    }
    if (std::strcmp(argv[1], "formatted_snapshot_mutation") == 0)
    {
        return Test_Diagnostic_Snapshot(argv[2], true, false, false);
    }
    if (std::strcmp(argv[1], "direct_snapshot_free") == 0)
    {
        return Test_Diagnostic_Snapshot(argv[2], false, true, false);
    }
    if (std::strcmp(argv[1], "formatted_snapshot_free") == 0)
    {
        return Test_Diagnostic_Snapshot(argv[2], true, true, false);
    }
    if (std::strcmp(argv[1], "direct_malloc_snapshot_free") == 0)
    {
        return Test_Diagnostic_Snapshot(argv[2], false, true, true);
    }
    if (std::strcmp(argv[1], "formatted_malloc_snapshot_free") == 0)
    {
        return Test_Diagnostic_Snapshot(argv[2], true, true, true);
    }
    if (std::strcmp(argv[1], "capture_allocation_failure") == 0)
    {
        return Test_Diagnostic_Capture_Allocation_Failure(argv[2]);
    }
    if (std::strcmp(argv[1], "capture_allocation_truncation") == 0)
    {
        return Test_Diagnostic_Capture_Truncation(argv[2]);
    }
    if (std::strcmp(argv[1], "long_formatted") == 0)
    {
        return Test_Long_Formatted_Diagnostic(argv[2]);
    }
    return Fail("unknown mode");
}
