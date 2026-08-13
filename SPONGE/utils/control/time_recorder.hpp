#pragma once

#include <chrono>

namespace Chrono = std::chrono;
using Clock = Chrono::high_resolution_clock;
using time_recorder_t = Chrono::time_point<Clock>;

// GPU 构建下默认用 device event 计时：Start/Stop 只在默认流上记录事件，
// 不做全设备同步，CPU 可以跨段、跨步提前入队 kernel。段耗时在下一次
// Start 或 Flush 时结算，届时事件通常已完成，结算基本不阻塞；同时这把
// CPU 领先 GPU 的深度自然限制在一个计时段以内。
// use_chrono_only = true 的计时器（如主要耗在 CPU 侧等待的 MPI
// Communication 段）退回到全设备同步 + 墙钟计时。
struct TIME_RECORDER
{
   private:
    time_recorder_t start_timestamp;
    time_recorder_t end_timestamp;
#ifdef GPU_ARCH_NAME
    deviceEvent_t event_start = nullptr;
    deviceEvent_t event_end = nullptr;
    bool event_started = false;
    bool event_pending = false;

    void Ensure_Events()
    {
        if (event_start == nullptr)
        {
            deviceEventCreate(&event_start);
            deviceEventCreate(&event_end);
        }
    }

    void Settle_Pending_Event()
    {
        if (!event_pending) return;
        deviceEventSynchronize(event_end);
        float milliseconds = 0.0f;
        deviceEventElapsedTime(&milliseconds, event_start, event_end);
        time += milliseconds * 1e-3;
        event_pending = false;
    }
#endif

   public:
    double time = 0;
    // true 时退回“全设备同步 + 墙钟”计时（用于 MPI 通信等 CPU 侧等待段）
    bool use_chrono_only = false;

    TIME_RECORDER() = default;
    // 拷贝只携带计时数值与模式；事件句柄不随拷贝转移，由拷贝体按需自建
    TIME_RECORDER(const TIME_RECORDER& other)
        : time(other.time), use_chrono_only(other.use_chrono_only)
    {
    }
    TIME_RECORDER& operator=(const TIME_RECORDER& other)
    {
        if (this != &other)
        {
            time = other.time;
            use_chrono_only = other.use_chrono_only;
#ifdef GPU_ARCH_NAME
            event_started = false;
            event_pending = false;
#endif
        }
        return *this;
    }

#ifdef GPU_ARCH_NAME
    ~TIME_RECORDER()
    {
        if (event_start != nullptr)
        {
            deviceEventDestroy(event_start);
            deviceEventDestroy(event_end);
        }
    }
#endif

    void Start()
    {
#ifdef GPU_ARCH_NAME
        if (!use_chrono_only)
        {
            Ensure_Events();
            Settle_Pending_Event();
            deviceEventRecord(event_start, 0);
            event_started = true;
            return;
        }
        hostDeviceSynchronize();
#endif
        start_timestamp = Clock::now();
    }

    void Stop()
    {
#ifdef GPU_ARCH_NAME
        if (!use_chrono_only)
        {
            if (event_started)
            {
                deviceEventRecord(event_end, 0);
                event_pending = true;
                event_started = false;
            }
            return;
        }
        hostDeviceSynchronize();
#endif
        end_timestamp = Clock::now();
        time += Chrono::duration_cast<Chrono::duration<double>>(end_timestamp -
                                                                start_timestamp)
                    .count();
    }

    // 汇总打印前调用：把尚未结算的 event 段计入 time
    void Flush()
    {
#ifdef GPU_ARCH_NAME
        if (!use_chrono_only) Settle_Pending_Event();
#endif
    }

    void Clear()
    {
        time = 0;
        start_timestamp = time_recorder_t();
        end_timestamp = time_recorder_t();
#ifdef GPU_ARCH_NAME
        event_started = false;
        event_pending = false;
#endif
    }
};
