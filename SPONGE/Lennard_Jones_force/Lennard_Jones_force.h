#pragma once
#include "../common.h"
#include "../control.h"

// 用于计算LJ_Force时使用的坐标和记录的原子LJ种类序号与原子电荷
#ifndef VECTOR_LJ_DEFINE
#define VECTOR_LJ_DEFINE
#define TWO_DIVIDED_BY_SQRT_PI 1.1283791670218446f
__host__ __device__ __forceinline__ int Get_LJ_Type(int a, int b);
__host__ __device__ __forceinline__ int Get_LJ_Type(unsigned int a,
                                                    unsigned int b);

struct VECTOR_LJ
{
    VECTOR crd;
    int LJ_type;
    float charge;
    int global_atom;
    friend __host__ __device__ __forceinline__ int Get_LJ_Type(int a, int b)
    {
        if (a < 0 || b < 0) return -1;
        return Get_LJ_Type(static_cast<unsigned int>(a),
                           static_cast<unsigned int>(b));
    }
    friend __host__ __device__ __forceinline__ int Get_LJ_Type(unsigned int a,
                                                               unsigned int b)
    {
        // 65,535 types are the largest triangular table whose last index fits
        // in a signed int.  Use a 64-bit product so every representable table
        // is supported without signed-overflow UB in host or device code.
        constexpr unsigned int max_type_index = 65534U;
        if (a > max_type_index || b > max_type_index) return -1;
        const unsigned int high = a > b ? a : b;
        const unsigned int low = a > b ? b : a;
        const unsigned long long pair =
            static_cast<unsigned long long>(high) * (high + 1ULL) / 2ULL + low;
        return static_cast<int>(pair);
    }
    friend __device__ __host__ __forceinline__ VECTOR Get_Periodic_Displacement(
        VECTOR_LJ uvec_a, VECTOR_LJ uvec_b, LTMatrix3 cell, LTMatrix3 rcell)
    {
        return Get_Periodic_Displacement(uvec_a.crd, uvec_b.crd, cell, rcell);
    }
    friend __device__ __host__ __forceinline__ float Get_LJ_Energy(
        VECTOR_LJ r1, VECTOR_LJ r2, float dr_abs, const float A, const float B)
    {
        float dr_6 = powf(dr_abs, -6.0f);
        return (0.083333333f * A * dr_6 - 0.166666667f * B) * dr_6;
    }
    friend __device__ __host__ __forceinline__ float Get_LJ_Force(
        VECTOR_LJ r1, VECTOR_LJ r2, float dr_abs, const float A, const float B)
    {
        return (B - A * powf(dr_abs, -6.0f)) * powf(dr_abs, -8.0f);
    }
    friend __device__ __host__ __forceinline__ float Get_LJ_Virial(
        VECTOR_LJ r1, VECTOR_LJ r2, float dr_abs, const float A, const float B)
    {
        float dr_6 = powf(dr_abs, -6.0f);
        return -(B - A * dr_6) * dr_6;
    }
    friend __device__ __host__ __forceinline__ float Get_Direct_Coulomb_Energy(
        VECTOR_LJ r1, VECTOR_LJ r2, float dr_abs, const float pme_beta)
    {
        return r1.charge * r2.charge * erfcf(pme_beta * dr_abs) / dr_abs;
    }
    friend __device__ __host__ __forceinline__ float Get_Direct_Coulomb_Force(
        VECTOR_LJ r1, VECTOR_LJ r2, float dr_abs, const float pme_beta)
    {
        float beta_dr = pme_beta * dr_abs;
        return r1.charge * r2.charge * powf(dr_abs, -3.0f) *
               (beta_dr * TWO_DIVIDED_BY_SQRT_PI * expf(-beta_dr * beta_dr) +
                erfcf(beta_dr));
    }
};
// 从 16B 对齐双数组（float4: crd.xyz+charge，int2: LJ_type+global_atom）
// 装载一个 VECTOR_LJ 值；消费 kernel 内的语义与旧 AoS 布局逐语句一致
__host__ __device__ __forceinline__ VECTOR_LJ Load_VECTOR_LJ(
    const float4* crd_q, const int2* type_g, const int idx)
{
    const float4 cq = crd_q[idx];
    const int2 tg = type_g[idx];
    VECTOR_LJ r = {{cq.x, cq.y, cq.z}, tg.x, cq.w, tg.y};
    return r;
}
__global__ void Copy_LJ_Type_To_New_Crd(const int atom_numbers,
                                        VECTOR_LJ* new_crd, const int* LJ_type);
// 每步重打包：只把 crd gather 进 lj_crd_q 的 xyz；w（charge）与
// LJ_type/global_atom 是静态的，由 Get_Local 在初始化/域刷新时填充
__global__ void Repack_LJ_Crd(const int atom_numbers, const VECTOR* crd,
                              float4* lj_crd_q);
#endif

// S3：tile 表视图（设计文档 _local/LJ_TILING_DESIGN.md §2）。LJ_TILE 定义在
// neighbor_list/neighbor_list.h，此处前向声明避免头文件互拉。指针在每次调用
// 时从 NEIGHBOR_LIST 现取（溢出恢复扩容会重绑定 device 数组），tile_numbers
// 取 host 侧回读的 h_lj_tile_numbers。
struct LJ_TILE;
struct LJ_TILE_SET
{
    const LJ_TILE* tiles = NULL;
    const int* cluster_atoms = NULL;
    const int* cluster_flags = NULL;
    // tile 按 cluster_i 分组排序的下标表（kernel 按行消费）
    const int* tile_sorted = NULL;
    int tile_numbers = 0;
    // cluster 原子槽总数（(local+ghost cluster)×8，含 padding），用于
    // cluster 序打包数组的定容与 gather 范围
    int cluster_atom_slots = 0;
};

// 用于记录与计算LJ相关的信息
struct LENNARD_JONES_INFORMATION
{
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260216;
    CONTROLLER* controller = NULL;

    // a = LJ_A between atom[i] and atom[j]
    // b = LJ_B between atom[i] and atom[j]
    // E_lj = a/12 * r^-12 - b/6 * r^-6;
    // F_lj = (a * r^-14 - b * r ^ -6) * dr
    int atom_numbers = 0;       // 原子数
    int atom_type_numbers = 0;  // 原子种类数
    int pair_type_numbers = 0;  // 原子对种类数

    int* h_atom_LJ_type = NULL;  // 原子对应的LJ种类
    int* d_atom_LJ_type = NULL;  // 原子对应的LJ种类

    float* h_LJ_A = NULL;  // LJ的A系数
    float* h_LJ_B = NULL;  // LJ的B系数
    float* d_LJ_A = NULL;  // LJ的A系数
    float* d_LJ_B = NULL;  // LJ的B系数

    float* h_LJ_energy_atom = NULL;  // 每个原子的LJ的能量
    float h_LJ_energy_sum = 0;       // 所有原子的LJ能量和
    float* d_LJ_energy_atom = NULL;  // 每个原子的LJ的能量
    float* d_LJ_energy_sum = NULL;
    ;                            // 所有原子的LJ能量和
    float h_LJ_long_energy = 0;  // 长程修正能量

    // 初始化
    void Initial(CONTROLLER* controller, float cutoff,
                 const char* module_name = NULL);
    // 分配内存
    void LJ_Malloc();
    // 参数传到GPU上
    void Parameter_Host_To_Device();

    float cutoff = 10.0;
    VECTOR_LJ* crd_with_LJ_parameters = NULL;

    // S3：warp-per-tile 力 kernel 开关（mdin 命令 "LJ use_tile"；CUDA 构建
    // 默认开，CPU/HIP 恒关走旧 warp-per-atom kernel）。开关只是意图，真正
    // 走 tile 路径还要求调用方传入非空的 LJ_TILE_SET。
    bool use_tile = false;

    /*
        以下用于区域分解
    */
    int local_atom_numbers = 0;
    int ghost_numbers = 0;
    // 局域原子的坐标/电荷与 LJ_type/global_atom 打包（atom 序，16B 对齐双数组）：
    // d_lj_crd_q 的 xyz 每步由 Repack_LJ_Crd 刷新，w（charge）与 d_lj_type_g
    // 是静态的，由 Get_Local 在初始化/域刷新时填充
    float4* d_lj_crd_q = NULL;
    int2* d_lj_type_g = NULL;
    // S3：cluster 序的 LJ 打包数据（slot = cluster*8 + 槽位，cluster 内 8
    // 槽连续；tile kernel 的 cluster 载入因此是全合并 LDG，无需经
    // cluster_atoms 的二级间接）。每步由 Gather_LJ_Cluster_Data 从
    // d_lj_crd_q/d_lj_type_g 重排；容量随邻居表重建的 cluster 槽数惰性增长
    float4* d_lj_cluster_crd_q = NULL;
    int2* d_lj_cluster_type_g = NULL;
    int lj_cluster_array_capacity = 0;
    int* d_pair_overlap_error = NULL;
    int* d_local_metadata_error = NULL;
    bool local_metadata_is_ready = false;
    void Get_Local(int* atom_local, int local_atom_numbers, int ghost_numbers,
                   const float* charge);  // 获取局域粒子信息
    bool Validate_Local_State(const char* error_by, int global_atom_numbers,
                              int local_atom_numbers, int ghost_numbers,
                              int solvent_numbers);
    void Reset_Pair_Overlap_Error();
    bool Check_Pair_Overlap_Error(const char* error_by);

    // 可以根据外界传入的need_atom_energy和need_virial，选择性计算能量和维里。其中的维里对PME直接部分计算的原子能量，在和PME其他部分加和后即维里。
    // tile_set 非空时走 S3 的 warp-per-tile kernel（覆盖全部 local i 的半表
    // 对，溶剂 dispatch 必须由调用方同时停用，否则水-水双计）；为空走旧的
    // warp-per-atom kernel（i 范围 [0, local_atom_numbers - solvent_numbers)）。
    void LJ_PME_Direct_Force_With_Atom_Energy_And_Virial(
        const int atom_numbers, const int local_atom_numbers,
        const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
        const float* charge, VECTOR* frc, const LTMatrix3 cell,
        const LTMatrix3 rcell, const ATOM_GROUP* nl, const float pme_beta,
        const int need_atom_energy, float* atom_energy, const int need_virial,
        LTMatrix3* atom_virial, float* atom_direct_pme_energy,
        const LJ_TILE_SET* tile_set);

#ifdef KPCCL_TASKLOOP
    void LJ_PME_Direct_Force_With_Atom_Energy_And_Virial_Init(
        const int atom_numbers, const int local_atom_numbers,
        const int ghost_numbers, const VECTOR* crd, const float* charge,
        VECTOR* frc, const LTMatrix3 cell, const LTMatrix3 rcell,
        const ATOM_GROUP* nl, const float pme_beta, const int need_atom_energy,
        float* atom_energy, const int need_virial, LTMatrix3* atom_virial,
        float* atom_direct_pme_energy);

    void LJ_PME_Direct_Force_With_Atom_Energy_And_Virial_KPCCL(
        const int atom_numbers, const int local_atom_numbers,
        const int ghost_numbers, const VECTOR* crd, const float* charge,
        VECTOR* frc, const LTMatrix3 cell, const LTMatrix3 rcell,
        const ATOM_GROUP* nl, const float pme_beta, const int need_atom_energy,
        float* atom_energy, const int need_virial, LTMatrix3* atom_virial,
        float* atom_direct_pme_energy);
#endif

    // 长程能量和维里修正
    float long_range_factor = 0;
    // 求力的时候对能量和维里的长程修正
    void Long_Range_Correction(int need_pressure, LTMatrix3* d_virial,
                               int need_potential, float* d_potential,
                               const float volume);

    void Step_Print(CONTROLLER* controller);
};
