#pragma once

#include <atomic>
#include <cerrno>
#include <cstdlib>
#include <new>
#include <stdexcept>
#include <thread>

#include "common.h"
#include "utils/control/string.hpp"
#include "utils/control/time_recorder.hpp"

using StringMap = std::map<std::string, std::string>;
using CheckMap = std::map<std::string, int>;
using StringVector = std::vector<std::string>;
using IntPairMap = std::map<int, std::pair<int, int>>;

struct CONTROLLER
{
    static unsigned int device_optimized_block;
    static unsigned int device_warp;
    static unsigned int device_max_thread;
    static int force_replica_count;

    // CPU MPI 相关
    // PP, PM 可能有多个进程，需要构造通信器；
    // CC 只有单进程情况，不需要构造通信器。
    static int MPI_rank;
    static int MPI_size;
    static int PP_MPI_size;
    static int PM_MPI_size;
    static int CC_MPI_size;

    static int PP_MPI_rank;
    static int PM_MPI_rank;
    static int CC_MPI_rank;

    static MPI_Comm pp_comm;
    static MPI_Comm pm_comm;

    // DEVICE MPI 相关（NCCL)
    static D_MPI_Comm D_MPI_COMM_WORLD;
    static D_MPI_Comm d_pp_comm;
    static D_MPI_Comm d_pm_comm;

    // 自身信息
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260216;

    void Initial(int argc, char** argv, const char* subpackage_hint);
    void Clear();

    // 文件输出缓冲大小
    int buffer_frame;
    // 设置文件缓冲区大小
    void Set_File_Buffer(FILE* file, size_t one_frame_size);

    // 输入控制部分
    StringMap original_commands;      // 由文件读入的原始命令
    StringMap commands;               // 由文件读入的命令（去除空格）
    CheckMap command_check;           // 检查输入是否都使用了
    CheckMap choice_check;            // 检查选项是否都使用了
    bool workspace_from_cli = false;  // 记录workspace是否由命令行设置
    bool mdin_is_toml = false;        // 当前mdin是否为TOML输入
    std::string mdin_toml_source_path;
    std::string mdin_toml_content;
    void Get_Command(
        const std::string& line,
        const std::string&
            prefix);  // 内部解析argument时专用，设置命令，不外部调用
    void Set_Command(
        const char* Flag, const char* Value, int Check = 1,
        const char* prefix = NULL,
        bool preserve_full_value =
            false);  // 内部解析argument时专用，设置命令，不外部调用
    void Arguments_Parse(int argc, char** argv);  // 对终端输入进行分析
    void Commands_From_In_File(
        int argc, char** argv,
        const char* subpackage_hint);    // 对mdin输入进行分析并打印日志信息
    void Default_Set();                  // 对最基本的功能进行默认设置
    int working_device = 0;              // 使用的设备
    void Init_Host_MPI();                // 对主机MPI初始化
    void Init_Device();                  // 对设备初始化
    void Initial_Force_Replica_Count();  // 初始化力/能量/virial副本数量
    void Init_Device_MPI();  // 对设备MPI初始化 (优先初始化xccl，否则使用mpi)
    // 本部分的上面的内容最好不要外部调用

    void
    Input_Check();  // 检查所有输入是否都被使用了（防止错误的输入未被检查到）

    bool Command_Exist(const char* key);  // 判断文件读入的命令中是否有key
    bool Command_Exist(const char* prefix,
                       const char* key);  // 判断文件读入的命令中是否有key
    // 判断是否存在key且值为value。未设置时返回no_set_return_value，可控制比对是否对大小写敏感
    bool Command_Choice(const char* key, const char* value,
                        bool case_sensitive = 0);
    // 判断是否存在prefix_key且值为value。未设置时返回no_set_return_value，可控制比对是否对大小写敏感
    bool Command_Choice(const char* prefix, const char* key, const char* value,
                        bool case_sensitive = 0);
    const char* Command(const char* key);  // 获得文件读入的命令key对应的value
    const char* Command(const char* prefix,
                        const char* key);  // 获得文件读入的命令key对应的value
    const char* Original_Command(
        const char* key);  // 获得文件读入的命令key对应的value
    const char* Original_Command(
        const char* prefix,
        const char* key);  // 获得文件读入的命令key对应的value

    // 计时部分
    TIME_RECORDER core_time;                              // 总计时器
    std::map<std::string, TIME_RECORDER> time_recorders;  // 各分部计时器
    std::vector<std::string> time_recorder_names;  // 各分部计时器名字（按顺序）
    float simulation_speed;                        // 模拟运行速度（纳秒/天）
    TIME_RECORDER* Get_Time_Recorder(const char* name);  // 计时器名字
    void Final_Time_Summary(int steps, float time_factor, const char* unit_name,
                            const int mode);  // 最后总结总CPU时间

    // 输出控制部分
    FILE* mdinfo = NULL;  // 屏幕信息打印文件
    FILE* mdout = NULL;
    StringMap outputs_content;  // 记录每步输出数值
    StringMap outputs_format;   // 记录每步输出的格式
    StringVector outputs_key;   // 记录每部输出的表头
    // 本部分的上面的内容最好不要外部调用

    float printf_sum = 0;
    void printf(const char* fmt,
                ...);  // 重载printf，使得printf能够同时打印到mdinfo和屏幕
    void MPI_printf(const char* fmt,
                    ...);  // 上面printf的MPI版本，将会按MPI_rank顺序打印
    void Step_Print_Initial(
        const char* head,
        const char* format);  // 其他模块初始化时调用，获得对应的表头和格式
    void Step_Print(
        const char* head, const float* pointer,
        const bool add_to_total =
            false);  // 其他模块打印时调用，获得对应的表头和数值，使之能以不同的格式打印在屏幕和mdout
    void Step_Print(
        const char* head, const float pointer,
        const bool add_to_total =
            false);  // 其他模块打印时调用，获得对应的表头和数值，使之能以不同的格式打印在屏幕和mdout
    void Step_Print(
        const char* head, const double pointer,
        const bool add_to_total =
            false);  // 其他模块打印时调用，获得对应的表头和数值，使之能以不同的格式打印在屏幕和mdout
    void Step_Print(
        const char* head,
        const int
            pointer);  // 其他模块打印时调用，获得对应的表头和数值，使之能以不同的格式打印在屏幕和mdout
    void Step_Print(
        const char* head,
        const char*
            pointer);  // 其他模块打印时调用，获得对应的表头和数值，使之能以不同的格式打印在屏幕和mdout
    void Print_First_Line_To_Mdout(
        FILE* mdout =
            NULL);  // 模拟开始前的操作，将表头打印到mdout，并在屏幕打印一个分割线
    void Print_To_Screen_And_Mdout(
        FILE* mdout =
            NULL);  // 模拟开始每步的调用，使得其他部分的结果打印到屏幕和mdout

    // 错误处理
    void Check_Error(float potential);
    void Check_Int(const char* command, const char* error_by);
    void Check_Float(const char* command, const char* error_by);
    bool Get_Bool(const char* command, const char* error_by);
    void Check_Int(const char* prefix, const char* command,
                   const char* error_by);
    void Check_Float(const char* prefix, const char* command,
                     const char* error_by);
    bool Get_Bool(const char* prefix, const char* command,
                  const char* error_by);
    // Fixed-capacity fatal cleanup stack shared by one linked SPONGE core
    // image. A callback/context pair must remain valid until successful
    // unregistration or callback return.
    // Successful unregistration guarantees that the callback is no longer
    // running. The sole exception is a callback unregistering itself: that
    // operation returns false rather than deadlocking or granting a false
    // lifetime guarantee. Registration returns false for a null callback, a
    // duplicate pair, a full stack, or a registry whose fatal path has begun.
    using Fatal_Cleanup_Callback = void (*)(void*);
    static constexpr std::size_t FATAL_CLEANUP_CAPACITY = 32;
    static constexpr std::size_t FATAL_DIAGNOSTIC_FALLBACK_CAPACITY = 8192;
    static bool Register_Fatal_Cleanup(Fatal_Cleanup_Callback callback,
                                       void* context) noexcept;
    static bool Unregister_Fatal_Cleanup(Fatal_Cleanup_Callback callback,
                                         void* context) noexcept;
#ifdef SPONGE_FATAL_CLEANUP_TESTING
    static std::size_t Fatal_Sequence_Waiter_Count_For_Test() noexcept
    {
        return fatal_sequence_waiters_.load(std::memory_order_acquire);
    }
#endif
    void Throw_SPONGE_Error(const int error_number, const char* error_by = NULL,
                            const char* extra_error_string = NULL);
#ifdef GPU_ARCH_NAME
    void Throw_Device_Error(deviceError_t error_number,
                            const char* error_by = NULL,
                            const char* extra_error_string = NULL);
#endif
    void Throw_Formatted_SPONGE_Error(const int error_number,
                                      const char* error_by, const char* format,
                                      ...);

    // 警告
    StringVector warnings;
    bool warn_of_initialization;
    void Warn(const char* warning);
    void Deprecated(const char* deprecated_command,
                    const char* recommanded_command, const char* version,
                    const char* reason);

    // 输出文件
    FILE* Get_Output_File(bool binary, const char* command,
                          const char* default_suffix,
                          const char* default_filename);
    FILE* Get_Output_File(bool binary, const char* prefix, const char* command,
                          const char* default_suffix,
                          const char* default_filename);

   private:
    static void Begin_Fatal_Sequence() noexcept;
    static void Run_Fatal_Cleanups() noexcept;
    static void Resolve_Fatal_Base(int error_number, const char** error_name,
                                   const char** error_reason) noexcept;
    static bool Capture_Fatal_Diagnostic_From_Base(
        const char* error_name, const char* error_reason, const char* error_by,
        const char* extra_error_string, std::vector<char>* diagnostic) noexcept;
    static std::size_t Capture_Fatal_Diagnostic_Fixed_From_Base(
        const char* error_name, const char* error_reason, const char* error_by,
        const char* extra_error_string, char* diagnostic,
        std::size_t capacity) noexcept;
    static bool Capture_Fatal_Diagnostic(
        int error_number, const char* error_by, const char* extra_error_string,
        std::vector<char>* diagnostic) noexcept;
    static void Write_Fatal_Bytes(const char* text,
                                  std::size_t length) noexcept;
    static void Write_Fatal_Text(const char* text) noexcept;
    [[noreturn]] static void Finish_Fatal_Sequence(
        int error_number, const char* captured_diagnostic,
        std::size_t captured_length) noexcept;
    [[noreturn]] static void Finish_Uncaptured_Fatal_Sequence(
        int error_number, const char* error_name, const char* error_reason,
        const char* error_by, const char* extra_error_string) noexcept;
    [[noreturn]] static void Terminate_Fatal_Sequence(
        int error_number) noexcept;
    [[noreturn]] static void Wait_For_Fatal_Termination() noexcept;
    static void Lock_Fatal_Cleanup_Registry() noexcept;
    static void Unlock_Fatal_Cleanup_Registry() noexcept;

    enum Fatal_Cleanup_Entry_State : unsigned char
    {
        FATAL_CLEANUP_EMPTY = 0,
        FATAL_CLEANUP_REGISTERED = 1,
        FATAL_CLEANUP_RUNNING = 2,
        FATAL_CLEANUP_DONE = 3,
    };

    struct Fatal_Cleanup_Entry
    {
        Fatal_Cleanup_Callback callback;
        void* context;
        Fatal_Cleanup_Entry_State state;
    };

    inline static Fatal_Cleanup_Entry
        fatal_cleanup_entries_[FATAL_CLEANUP_CAPACITY] = {};
    inline static std::size_t fatal_cleanup_count_ = 0;
    // The spin lock protects only short registry state transitions. Cleanup
    // callbacks always execute with it released. The fatal-sequence claim is
    // separate so exactly one thread owns cleanup, diagnostics, and process
    // termination; all other fatal callers wait for that termination.
    inline static std::atomic_flag fatal_cleanup_lock_ = ATOMIC_FLAG_INIT;
    inline static std::atomic<bool> fatal_sequence_claimed_{false};
    inline static std::atomic<std::size_t> fatal_sequence_waiters_{0};
    inline static thread_local bool fatal_sequence_owner_ = false;
};

// 一些常用的控制函数
#include "utils/control/error.hpp"
#include "utils/control/malloc_and_file.hpp"
#include "utils/control/os.hpp"
#include "utils/control/print.hpp"
#include "utils/control/warning.hpp"

// 读入复杂参数
#include "utils/control/configuration_reader.hpp"
