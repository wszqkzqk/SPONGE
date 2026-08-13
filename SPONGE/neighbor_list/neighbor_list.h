#ifndef NEIGHBOR_LIST_H
#define NEIGHBOR_LIST_H
#include "../common.h"
#include "../control.h"
#include "full_neighbor_list.h"

// LJ cluster-pair tile 表（设计文档 _local/LJ_TILING_DESIGN.md §1）。
// cluster = 格桶顺序的 8 原子一组（编译期常量，4 留作调参对照）。
#define LJ_TILE_CLUSTER_SIZE 8
struct LJ_TILE
{
    int cluster_i;             // 非 ghost cluster（local 段编号）
    int cluster_j;             // 可为 ghost cluster（统一编号：local 段在前）
    unsigned long long mask;   // 8x8，bit = i_local*8 + j_local
};

// 每 global 原子的排除范围缓存条目（list_start, n, min, max）
struct LJ_EXCL_RANGE
{
    int list_start;
    int n;
    int min_atom;
    int max_atom;
};

struct NEIGHBOR_LIST
{
    enum UPDATE_ERROR
    {
        UPDATE_OK = 0,
        UPDATE_INVALID_GEOMETRY = 1
    };

    bool is_initialized = 0;
    bool throw_error_when_overflow = 0;
    int active_local_atom_numbers = 0;
    int last_update_error = UPDATE_OK;
    int last_error_atom = -1;
    int* d_update_error = NULL;

    // 构建状态块：溢出标志与错误码合并在同一块连续 device 内存里，
    // 一次 memset 复位、一次 DtoH 取回。下面四个指针是块内字段的别名，
    // 不单独分配/释放。
    enum UPDATE_STATUS_FIELD
    {
        UPDATE_STATUS_GRID_OVERFLOW = 0,
        UPDATE_STATUS_GRID_GHOST_OVERFLOW = 1,
        UPDATE_STATUS_LIST_OVERFLOW = 2,
        UPDATE_STATUS_TILE_OVERFLOW = 3,
        UPDATE_STATUS_ERROR = 4,  // [4]=错误码，[5]=出错原子
        UPDATE_STATUS_TILE_COUNT = 6,  // [6]=本次重建的 tile 数
        UPDATE_STATUS_INTS = 8
    };
    int* d_update_status = NULL;

    // 是否需要构建半近邻表（默认需要）
    bool is_needed_half = true;
    // 是否需要构建全近邻表（默认不需要，由具体势函数开启）
    bool is_needed_full = false;

    void Initial(CONTROLLER* controller, int atom_numbers, float cutoff,
                 float skin, LTMatrix3 cell, LTMatrix3 rcell);
    bool Update(int* atom_local, int local_atom_numbers, int ghost_numbers,
                VECTOR* crd, LTMatrix3 cell, LTMatrix3 rcell, int step,
                int update, int* excluded_list_start = NULL,
                int* excluded_list = NULL,
                int* excluded_numbers = NULL);  // 这里用NULL先不考虑排除表
    // Build a complete list before returning.  If a grid, half-list, or
    // full-list capacity is exceeded, grow the affected storage and retry the
    // same state instead of allowing force kernels to consume a truncated
    // list.  The return value reports whether any exported buffer was rebound.
    bool Update_With_Overflow_Recovery(
        CONTROLLER* controller, int* atom_local, int local_atom_numbers,
        int ghost_numbers, VECTOR* crd, LTMatrix3 cell, LTMatrix3 rcell,
        int step, int update, int* excluded_list_start = NULL,
        int* excluded_list = NULL, int* excluded_numbers = NULL);
    void Clear();

    float cutoff = 0.0f;
    float skin = 0.0f;

    // local的粒子数目
    int atom_numbers = 0;
    // ghost粒子数目
    int ghost_numbers = 0;
    int* neighbor_num = NULL;  // 当前区域不同方向的邻居数目
    // 当前区域的x,y,z坐标的最小的点和最大的点
    VECTOR min_corner;
    VECTOR max_corner;
    // 当前区域的box_length，dom_box_length=max_corner-min_corner
    VECTOR dom_box_length;
    // 近邻表
    int* d_temp = NULL;
    ATOM_GROUP *h_nl = NULL, *d_nl = NULL;
    // 每个原子的最大近邻数
    int max_neighbor_numbers = 0;
    // 每个格子的最大原子数
    int max_atom_in_grid_numbers = 0;
    // 每个原子的近邻数溢出
    // Zero means that the current build fit.  A positive value is the minimum
    // capacity required by at least one row/grid in the rejected build.
    int h_neighbor_list_overflow = 0, *d_neighbor_list_overflow = NULL;
    // 每个格子的原子数溢出
    int h_neighbor_grid_overflow = 0, *d_neighbor_grid_overflow = NULL;

    // 每个格子的最大ghost粒子数目
    int max_ghost_in_grid_numbers = 0;
    // 每个格子的ghost数溢出
    int h_neighbor_grid_ghost_overflow = 0,
        *d_neighbor_grid_ghost_overflow = NULL;

    // ---- LJ cluster-pair tile 表（S2，仅 USE_GPU 路径分配/构建，暂无消费者）----
    // cluster = 格桶顺序的 8 原子一组；local 段与 ghost 段统一编号
    //（ghost 段接在 local 段之后），与 d_lj_crd_q 的 local+ghost 布局一致
    int* d_grid_cluster_base = NULL;        // 每格 local cluster 起始编号（独占前缀）
    int* d_grid_ghost_cluster_base = NULL;  // 每格 ghost cluster 起始编号（含 local 段偏移）
    int* d_lj_cluster_numbers = NULL;       // [0]=local cluster 总数，[1]=ghost cluster 总数
    // 每 cluster 8 个统一坐标索引（crd/d_lj_crd_q 的下标），尾部 padding 为 -1
    int* d_lj_cluster_atoms = NULL;
    int* d_lj_cluster_flags = NULL;  // bit0 = 是否 ghost cluster
    int lj_cluster_capacity = 0;
    LJ_TILE* d_lj_tiles = NULL;
    int lj_tile_capacity = 0;
    // 每 global 原子的排除范围 (list_start, n, min, max) 缓存：排除表在
    // 运行期是静态的，首次重建时一次性预计算，把 tile 构建里的排除查询
    // 从 3-4 级依赖 gather 降为单次 16B 读
    LJ_EXCL_RANGE* d_lj_excl_range = NULL;
    int lj_excl_range_built = 0;
    // 溢出恢复时跨 Clear/Initial 保留的 tile 容量提示（Clear 不复位）
    int lj_tile_capacity_hint = 0;
    int h_lj_tile_numbers = 0;  // 上次成功重建的 tile 数（host 侧回读）
    // tile 表容量溢出（0=本次构建装得下；否则为被拒绝构建所需的最小容量）
    int h_neighbor_tile_overflow = 0, *d_neighbor_tile_overflow = NULL;
    int* d_lj_tile_count = NULL;  // 状态块槽位 6 的别名：构建期 atomicAdd 计数

    struct GRIDS
    {
        // 总的格点数
        int grid_numbers = 0;
        // 格点在三个方向的数量
        int Nx = 0, Ny = 0, Nz = 0;
        // 每个格点的周围格点数目
        int *h_neighbor_grid_numbers = NULL, *d_neighbor_grid_numbers = NULL;
        // 每个格点的周围格点
        int *h_neighbor_grids = NULL, *d_neighbor_grids = NULL;
        // 每个格点的邻居格原子数独占前缀和（每次重建时重算，
        // 行长 MAX_GRID_NEIGHBORS+1，仅 GPU 路径使用）
        int* d_grid_neighbor_prefix = NULL;
        // 每个格点内的原子数量
        int *h_grid_atom_numbers = NULL, *d_grid_atom_numbers = NULL;
        // 每个格点内的ghost数量
        int *h_grid_ghost_numbers = NULL, *d_grid_ghost_numbers = NULL;
        // 每个格点内的原子
        int *h_grid_atoms = NULL, *d_grid_atoms = NULL;
        // 每个格点内的ghost
        int *h_grid_ghosts = NULL, *d_grid_ghosts = NULL;
        // 每个格点内原子的坐标
        VECTOR* d_grid_atom_crd = NULL;
        // 每个格点内ghost的坐标
        VECTOR* d_grid_ghost_crd = NULL;
        // 初始化格点信息
        bool Initial(CONTROLLER* controller, int max_atom_in_grid_numbers,
                     int max_ghost_in_grid_numbers, LTMatrix3 cell,
                     LTMatrix3 rcell, float grid_length);
        // 释放内存
        void Clear();
    } grids;

    struct UPDATOR
    {
        TIME_RECORDER* time_recorder = NULL;
        // 更新间隔
        int refresh_interval = 0;
        // 是否需要更新
        int h_need_update = 0, *d_need_update = NULL;
        // 当某原子移动若干距离以后更新
        float skin_permit = 0.5;
        // 上次更新的坐标，用于判断是否需要更新
        VECTOR* old_crd = NULL;
        bool Initial(CONTROLLER* controller, int atom_numbers);
        void Check(int atom_numbers, float skin, VECTOR* crd, LTMatrix3 cell,
                   LTMatrix3 rcell);
        void Update(int* atom_local, int local_atom_numbers, int ghost_numbers,
                    int need_copy, VECTOR* crd, LTMatrix3 cell, LTMatrix3 rcell,
                    NEIGHBOR_LIST::GRIDS* grids, int max_atom_in_grid_numbers,
                    int max_ghost_in_grid_numbers, int max_neighbor_numbers,
                    float grid_length, int* d_neighbor_grid_overflow,
                    int* d_neighbor_grid_ghost_overflow,
                    int* d_neighbor_list_overflow, int* d_update_error,
                    ATOM_GROUP* d_nl, int* excluded_list_start = NULL,
                    int* excluded_list = NULL, int* excluded_numbers = NULL,
                    NEIGHBOR_LIST* owner = NULL);
        void Clear();
    } updator;

    enum NEIGHBOR_LIST_UPDATE_PARAMETER
    {
        CONDITIONAL_UPDATE = 0,
        FORCED_UPDATE = 1
    };

    // 全连接近邻表（用于需要全连接表的计算如SW）
    FULL_NEIGHBOR_LIST full_neighbor_list;

    // 静态变量：各模块可以注册它们需要的特殊截断距离
    float cutoff_full = 0.0f;  // 默认值为0表示不需要cutoff_full
};

#endif
