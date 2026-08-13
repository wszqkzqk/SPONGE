#include "Lennard_Jones_force.h"

#include <cstdint>
#include <cstdlib>

#include "../neighbor_list/neighbor_list.h"
#include "../xponge/load/native/lj.hpp"
#include "../xponge/xponge.h"
#include "pair_activity.h"
// #include "assert.h"

// 由LJ坐标和转化系数求距离
__global__ void Copy_LJ_Type_To_New_Crd(const int atom_numbers,
                                        VECTOR_LJ* new_crd, const int* LJ_type)
{
    SIMPLE_DEVICE_FOR(atom_i, atom_numbers)
    {
        new_crd[atom_i].LJ_type = LJ_type[atom_i];
        new_crd[atom_i].global_atom = atom_i;
    }
}

__global__ void Repack_LJ_Crd(const int atom_numbers, const VECTOR* crd,
                              float4* lj_crd_q)
{
    SIMPLE_DEVICE_FOR(atom_i, atom_numbers)
    {
        const VECTOR r = crd[atom_i];
        lj_crd_q[atom_i].x = r.x;
        lj_crd_q[atom_i].y = r.y;
        lj_crd_q[atom_i].z = r.z;
    }
}

// S3：把 atom 序的 LJ 打包数据按 cluster 序重排（每步一次；cluster 内 8 槽
// 连续，tile kernel 的 cluster 载入变为全合并 LDG，且不再需要经
// cluster_atoms 的二级间接 gather）。padding 槽（cluster_atoms 为 -1）写零，
// 构建保证对应掩码位恒 0，内容不会被消费。
__global__ void Gather_LJ_Cluster_Data(const int cluster_atom_slots,
                                       const int* cluster_atoms,
                                       const float4* crd_q, const int2* type_g,
                                       float4* cluster_crd_q,
                                       int2* cluster_type_g)
{
    SIMPLE_DEVICE_FOR(slot, cluster_atom_slots)
    {
        const int atom = cluster_atoms[slot];
        if (atom >= 0)
        {
            cluster_crd_q[slot] = crd_q[atom];
            cluster_type_g[slot] = type_g[atom];
        }
        else
        {
            const float4 zero4 = {0.0f, 0.0f, 0.0f, 0.0f};
            const int2 pad2 = {0, -1};
            cluster_crd_q[slot] = zero4;
            cluster_type_g[slot] = pad2;
        }
    }
}

static __global__ void device_add(float* variable, const float adder)
{
    variable[0] += adder;
}

template <bool need_force, bool need_energy, bool need_virial,
          bool need_coulomb>
static __global__ void Lennard_Jones_And_Direct_Coulomb_Device(
    const int local_atom_numbers, const int solvent_numbers,
    const ATOM_GROUP* nl, const float4* crd_q, const int2* type_g,
    const LTMatrix3 cell, const LTMatrix3 rcell, const float* LJ_type_A,
    const float* LJ_type_B, const float cutoff, VECTOR* frc,
    const float pme_beta, float* atom_energy, LTMatrix3* atom_virial,
    float* atom_direct_cf_energy, float* atom_LJ_ene, int* pair_overlap_error)
{
#ifdef USE_GPU
    int atom_i = 0 + blockDim.y * blockIdx.x + threadIdx.y;
    if (atom_i < local_atom_numbers - solvent_numbers)
#else
#pragma omp parallel for schedule(dynamic)
    for (int atom_i = 0; atom_i < local_atom_numbers - solvent_numbers;
         atom_i++)
#endif
    {
        VECTOR frc_record = {0.0f, 0.0f, 0.0f};
        LTMatrix3 virial = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        float energy_lj = 0.0f;
        float energy_coulomb = 0.0f;
        float energy_total = 0.0f;
        ATOM_GROUP nl_i = nl[atom_i];
        VECTOR_LJ r1 = Load_VECTOR_LJ(crd_q, type_g, atom_i);
#ifdef USE_GPU
        for (int j = threadIdx.x; j < nl_i.atom_numbers; j += blockDim.x)
#else
        for (int j = 0; j < nl_i.atom_numbers; j += 1)
#endif
        {
            int atom_j = nl_i.atom_serial[j];
            float ij_factor = atom_j < local_atom_numbers ? 1.0f : 0.5f;
            VECTOR_LJ r2 = Load_VECTOR_LJ(crd_q, type_g, atom_j);
            VECTOR dr = Get_Periodic_Displacement(r2, r1, cell, rcell);
            float dr_abs = norm3df(dr.x, dr.y, dr.z);
            if (dr_abs < cutoff)
            {
                int atom_pair_LJ_type = Get_LJ_Type(r1.LJ_type, r2.LJ_type);
                float A = LJ_type_A[atom_pair_LJ_type];
                float B = LJ_type_B[atom_pair_LJ_type];
                const PairwiseInteraction::Pair_Activity activity =
                    PairwiseInteraction::Classify(
                        A, B,
                        need_coulomb && PairwiseInteraction::Coulomb_Is_Active(
                                            r1.charge, r2.charge));
                if (!activity.Any())
                {
                    continue;
                }
                if (dr_abs == 0.0f)
                {
                    PairwiseInteraction::Fail_Exact_Overlap(
                        r1.global_atom, r2.global_atom,
                        PairwiseInteraction::Components(activity.lennard_jones,
                                                        activity.coulomb),
                        pair_overlap_error);
                    continue;
                }
                if (need_force)
                {
                    float frc_abs = 0.0f;
                    if (activity.lennard_jones)
                    {
                        frc_abs = Get_LJ_Force(r1, r2, dr_abs, A, B);
                    }
                    if (activity.coulomb)
                    {
                        float frc_cf_abs =
                            Get_Direct_Coulomb_Force(r1, r2, dr_abs, pme_beta);
                        frc_abs = frc_abs - frc_cf_abs;
                    }
                    VECTOR frc_lin = frc_abs * dr;
                    frc_record = frc_record + frc_lin;
                    if (atom_j < local_atom_numbers)
                    {
                        atomicAdd(frc + atom_j, -frc_lin);
                    }
                    if (need_virial)
                    {
                        virial = virial - ij_factor * Get_Virial_From_Force_Dis(
                                                          frc_lin, dr);
                    }
                }
                if (need_energy)
                {
                    if (activity.lennard_jones)
                    {
                        energy_lj +=
                            ij_factor * Get_LJ_Energy(r1, r2, dr_abs, A, B);
                    }
                    if (activity.coulomb)
                    {
                        energy_coulomb +=
                            ij_factor *
                            Get_Direct_Coulomb_Energy(r1, r2, dr_abs, pme_beta);
                    }
                }
            }
        }
        energy_total = energy_lj + energy_coulomb;
        if (need_force)
        {
            Warp_Sum_To(frc + atom_i, frc_record, warpSize);
        }
        if (need_energy)
        {
            Warp_Sum_To(atom_energy + atom_i, energy_total, warpSize);
            Warp_Sum_To(atom_LJ_ene + atom_i, energy_lj, warpSize);
            if (need_coulomb)
                Warp_Sum_To(atom_direct_cf_energy + atom_i, energy_coulomb,
                            warpSize);
        }
        if (need_virial)
        {
            Warp_Sum_To(atom_virial + atom_i, virial, warpSize);
        }
    }
}

#ifdef USE_CUDA
// ---- S3：tile LJ+直接库仑 kernel（设计文档 LJ_TILING_DESIGN.md §2）----
// tile = 8×8=64 对；lane l 负责 i=l/4、j∈{l%4, l%4+4} 两对。tile 表按
// cluster_i 分组排序（neighbor_list 重建时计数排序），一个 warp 连续消费
// LJ_TILE_WARP_SPAN 个 tile，同一 i-cluster 段内 f_i/能量/维里在寄存器
// 累计、段末一次归约写出（i 侧全局原子加是全 kernel 的吞吐瓶颈，实测占
// 首版 tile kernel 的 2/3）。活对的逐语句语义与
// Lennard_Jones_And_Direct_Coulomb_Device 完全一致（dr、dr_abs<cutoff
// 复判、Get_LJ_Type 查表、pair_activity 分类、零重叠 trap），差别只在
// 数据来源（tile 掩码 + cluster 序打包坐标）与归约写出方式。

// 4-lane 组内 / stride-4 跨组的 xor-shuffle 归约辅助
static __device__ __forceinline__ VECTOR LJ_Tile_Xor_Sum(VECTOR v,
                                                         const int offset)
{
    v.x += __shfl_xor_sync(FULL_MASK, v.x, offset);
    v.y += __shfl_xor_sync(FULL_MASK, v.y, offset);
    v.z += __shfl_xor_sync(FULL_MASK, v.z, offset);
    return v;
}

static __device__ __forceinline__ LTMatrix3 LJ_Tile_Xor_Sum(LTMatrix3 m,
                                                            const int offset)
{
    m.a11 += __shfl_xor_sync(FULL_MASK, m.a11, offset);
    m.a21 += __shfl_xor_sync(FULL_MASK, m.a21, offset);
    m.a22 += __shfl_xor_sync(FULL_MASK, m.a22, offset);
    m.a31 += __shfl_xor_sync(FULL_MASK, m.a31, offset);
    m.a32 += __shfl_xor_sync(FULL_MASK, m.a32, offset);
    m.a33 += __shfl_xor_sync(FULL_MASK, m.a33, offset);
    return m;
}

// SPONGE_LJ_FORCE_ABLATE：tile kernel 消融测时（仅诊断，结果不正确）：
// bit0 跳过对相互作用评估（保留 dr/掩码），bit1 跳过归约与写出，
// bit2 跳过 cluster 打包数据载入（r1/r2 用常量，atom 下标照常读）
static int LJ_Force_Ablate_Mode()
{
    static int mode = -1;
    if (mode < 0)
    {
        const char* v = std::getenv("SPONGE_LJ_FORCE_ABLATE");
        mode = v != NULL ? atoi(v) : 0;
    }
    return mode;
}

// 每 warp 连续消费的分组序 tile 数：tile 表已按 cluster_i 排序，warp 顺序
// 走过同一段（同一 i-cluster）时 f_i/能量/维里在寄存器累计，段末一次归约
// 写出——i 侧全局原子加从每 tile 8 次降为每段 8 次（段均 ~76 tile）
#define LJ_TILE_WARP_SPAN 64

template <bool need_force, bool need_energy, bool need_virial,
          bool need_coulomb, int ABLATE = 0>
static __global__ void __launch_bounds__(256)
Lennard_Jones_And_Direct_Coulomb_Tile_Device(
    const int tile_numbers, const LJ_TILE* tiles, const int* tile_sorted,
    const int* cluster_atoms, const int* cluster_flags,
    const float4* cluster_crd_q, const int2* cluster_type_g,
    const LTMatrix3 cell, const LTMatrix3 rcell, const float* LJ_type_A,
    const float* LJ_type_B, const float cutoff, VECTOR* frc,
    const float pme_beta, float* atom_energy, LTMatrix3* atom_virial,
    float* atom_direct_cf_energy, float* atom_LJ_ene, int* pair_overlap_error)
{
    // 要求 32 lane warp（CUDA 上恒成立）；tile 段界是 warp 一致的，越界的
    // warp 整体退出，循环内 shuffle 均在满 warp 下执行
    const int lane = threadIdx.x & (warpSize - 1);
    const int warp_global =
        blockIdx.x * (blockDim.x / warpSize) + threadIdx.x / warpSize;
    const int tile_begin = warp_global * LJ_TILE_WARP_SPAN;
    if (tile_begin >= tile_numbers) return;
    const int tile_end = min(tile_begin + LJ_TILE_WARP_SPAN, tile_numbers);

    const int i_local = lane >> 2;  // 0..7，4 lane 一组共享同一 i
    const int j_local[2] = {lane & 3, (lane & 3) + 4};

    // 当前 i 段的累计状态（cluster 序打包：slot = cluster*8 + 槽位，连续
    // 8 槽全合并读；同 4-lane 组读同一 i 槽、stride-4 的 8 个 lane 读同一
    // j 槽，L1 广播合并。atom 下标仅用于力回写与 padding 判定）
    int cur_ci = -1;
    int atom_i = -1;
    VECTOR_LJ r1 = {{0.11f, 0.23f, 0.37f}, 0, 0.5f, 0};
    VECTOR frc_i = {0.0f, 0.0f, 0.0f};
    float energy_lj = 0.0f;
    float energy_coulomb = 0.0f;
    LTMatrix3 virial = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    // 段末冲刷：4-lane 组内 xor 归约后组长一次原子加。组内无活对（或活对
    // 全被 cutoff 复判剔除）时部分和恒零，跳过原子加——语义等价（加 0
    // 不改变值），而全局原子加是本 kernel 的吞吐瓶颈
    auto flush_i = [&]()
    {
        if (ABLATE & 2) return;
        if (need_force)
        {
            frc_i = LJ_Tile_Xor_Sum(frc_i, 1);
            frc_i = LJ_Tile_Xor_Sum(frc_i, 2);
            if ((lane & 3) == 0 && atom_i >= 0 &&
                (frc_i.x != 0.0f || frc_i.y != 0.0f || frc_i.z != 0.0f))
            {
                atomicAdd(frc + atom_i, frc_i);
            }
        }
        if (need_energy)
        {
            // 能量按对记入 i 侧的三个 per-atom 数组（×ij_factor 已在循环内乘）
            float energy_total = energy_lj + energy_coulomb;
            for (int offset = 1; offset < 4; offset <<= 1)
            {
                energy_total += __shfl_xor_sync(FULL_MASK, energy_total, offset);
                energy_lj += __shfl_xor_sync(FULL_MASK, energy_lj, offset);
                energy_coulomb +=
                    __shfl_xor_sync(FULL_MASK, energy_coulomb, offset);
            }
            if ((lane & 3) == 0 && atom_i >= 0)
            {
                atomicAdd(atom_energy + atom_i, energy_total);
                atomicAdd(atom_LJ_ene + atom_i, energy_lj);
                if (need_coulomb)
                    atomicAdd(atom_direct_cf_energy + atom_i, energy_coulomb);
            }
        }
        if (need_virial)
        {
            virial = LJ_Tile_Xor_Sum(virial, 1);
            virial = LJ_Tile_Xor_Sum(virial, 2);
            if ((lane & 3) == 0 && atom_i >= 0)
            {
                atomicAdd(atom_virial + atom_i, virial);
            }
        }
    };

    for (int t = tile_begin; t < tile_end; ++t)
    {
        const LJ_TILE tile = tiles[tile_sorted[t]];
        if (tile.cluster_i != cur_ci)
        {
            if (cur_ci >= 0) flush_i();
            cur_ci = tile.cluster_i;
            const int slot_i = cur_ci * LJ_TILE_CLUSTER_SIZE + i_local;
            atom_i = cluster_atoms[slot_i];
            if (!(ABLATE & 4))
                r1 = Load_VECTOR_LJ(cluster_crd_q, cluster_type_g, slot_i);
            frc_i.x = 0.0f;
            frc_i.y = 0.0f;
            frc_i.z = 0.0f;
            energy_lj = 0.0f;
            energy_coulomb = 0.0f;
            virial.a11 = 0.0f;
            virial.a21 = 0.0f;
            virial.a22 = 0.0f;
            virial.a31 = 0.0f;
            virial.a32 = 0.0f;
            virial.a33 = 0.0f;
        }
        // ij_factor 是 tile 级常量：local×local tile 恒 1，local×ghost tile
        // 恒 0.5（ghost cluster_j 的能量/维里减半、j 力不写、i 力全额）
        const bool j_is_ghost = (cluster_flags[tile.cluster_j] & 1) != 0;
        const float ij_factor = j_is_ghost ? 0.5f : 1.0f;
        const int* atoms_j =
            cluster_atoms + tile.cluster_j * LJ_TILE_CLUSTER_SIZE;
        const int slot_j[2] = {tile.cluster_j * LJ_TILE_CLUSTER_SIZE +
                                   j_local[0],
                               tile.cluster_j * LJ_TILE_CLUSTER_SIZE +
                                   j_local[1]};
        const int atom_j[2] = {atoms_j[j_local[0]], atoms_j[j_local[1]]};
        VECTOR_LJ r2[2] = {{{0.42f, 0.53f, 0.61f}, 0, -0.5f, 1},
                           {{0.74f, 0.85f, 0.96f}, 0, 0.25f, 2}};
        if (!(ABLATE & 4))
        {
            r2[0] = Load_VECTOR_LJ(cluster_crd_q, cluster_type_g, slot_j[0]);
            r2[1] = Load_VECTOR_LJ(cluster_crd_q, cluster_type_g, slot_j[1]);
        }
        VECTOR frc_j[2] = {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 0.0f}};

#pragma unroll
        for (int k = 0; k < 2; ++k)
        {
            // 掩码 bit 在重建时刻判定（带 skin），运行期按 cutoff 复判
            if (((tile.mask >> (i_local * LJ_TILE_CLUSTER_SIZE + j_local[k])) &
                 1ULL) == 0ULL)
            {
                continue;
            }
            if (atom_i < 0 || atom_j[k] < 0) continue;  // padding 保险
            VECTOR dr =
                Get_Periodic_Displacement(r2[k].crd, r1.crd, cell, rcell);
            float dr_abs = norm3df(dr.x, dr.y, dr.z);
            if (ABLATE & 1)
            {
                energy_lj += dr_abs * 1.0e-30f;  // 消融：跳过相互作用评估
                continue;
            }
            if (dr_abs < cutoff)
            {
                int atom_pair_LJ_type = Get_LJ_Type(r1.LJ_type, r2[k].LJ_type);
                float A = LJ_type_A[atom_pair_LJ_type];
                float B = LJ_type_B[atom_pair_LJ_type];
                const PairwiseInteraction::Pair_Activity activity =
                    PairwiseInteraction::Classify(
                        A, B,
                        need_coulomb && PairwiseInteraction::Coulomb_Is_Active(
                                            r1.charge, r2[k].charge));
                if (!activity.Any())
                {
                    continue;
                }
                if (dr_abs == 0.0f)
                {
                    PairwiseInteraction::Fail_Exact_Overlap(
                        r1.global_atom, r2[k].global_atom,
                        PairwiseInteraction::Components(activity.lennard_jones,
                                                        activity.coulomb),
                        pair_overlap_error);
                    continue;
                }
                if (need_force)
                {
                    float frc_abs = 0.0f;
                    if (activity.lennard_jones)
                    {
                        frc_abs = Get_LJ_Force(r1, r2[k], dr_abs, A, B);
                    }
                    if (activity.coulomb)
                    {
                        float frc_cf_abs = Get_Direct_Coulomb_Force(
                            r1, r2[k], dr_abs, pme_beta);
                        frc_abs = frc_abs - frc_cf_abs;
                    }
                    VECTOR frc_lin = frc_abs * dr;
                    frc_i = frc_i + frc_lin;
                    frc_j[k] = frc_j[k] - frc_lin;
                    if (need_virial)
                    {
                        virial = virial - ij_factor * Get_Virial_From_Force_Dis(
                                                          frc_lin, dr);
                    }
                }
                if (need_energy)
                {
                    if (activity.lennard_jones)
                    {
                        energy_lj +=
                            ij_factor * Get_LJ_Energy(r1, r2[k], dr_abs, A, B);
                    }
                    if (activity.coulomb)
                    {
                        energy_coulomb +=
                            ij_factor * Get_Direct_Coulomb_Energy(r1, r2[k],
                                                                  dr_abs,
                                                                  pme_beta);
                    }
                }
            }
        }

        // f_j：同一 j 的部分和在 stride-4 的 8 个 lane 上，跨组 xor 归约后
        // lane 0..3 各落一个 j（零和跳过写出）；ghost tile 整体跳过
        if (!(ABLATE & 2) && need_force && !j_is_ghost)
        {
#pragma unroll
            for (int k = 0; k < 2; ++k)
            {
                frc_j[k] = LJ_Tile_Xor_Sum(frc_j[k], 4);
                frc_j[k] = LJ_Tile_Xor_Sum(frc_j[k], 8);
                frc_j[k] = LJ_Tile_Xor_Sum(frc_j[k], 16);
            }
            if (lane < 4)
            {
                if (atom_j[0] >= 0 &&
                    (frc_j[0].x != 0.0f || frc_j[0].y != 0.0f ||
                     frc_j[0].z != 0.0f))
                    atomicAdd(frc + atom_j[0], frc_j[0]);
                if (atom_j[1] >= 0 &&
                    (frc_j[1].x != 0.0f || frc_j[1].y != 0.0f ||
                     frc_j[1].z != 0.0f))
                    atomicAdd(frc + atom_j[1], frc_j[1]);
            }
        }
    }
    if (cur_ci >= 0) flush_i();
    if (ABLATE & 2)
    {
        // 消融：跳过全部归约与写出；以下条件运行时恒假，仅防死代码消除
        if (frc_i.x == 1.0e30f && energy_lj == 1.0e30f && virial.a11 == 1.0e30f)
        {
            pair_overlap_error[0] = atom_i;
        }
    }
}
#endif  // USE_CUDA

void LENNARD_JONES_INFORMATION::LJ_Malloc()
{
    Malloc_Safely((void**)&h_atom_LJ_type, sizeof(int) * atom_numbers);
    Malloc_Safely((void**)&h_LJ_A, sizeof(float) * pair_type_numbers);
    Malloc_Safely((void**)&h_LJ_B, sizeof(float) * pair_type_numbers);
    Malloc_Safely((void**)&h_LJ_energy_atom, sizeof(float) * atom_numbers);
}

void LENNARD_JONES_INFORMATION::Initial(CONTROLLER* controller, float cutoff,
                                        const char* module_name)
{
    this->controller = controller;
    if (module_name == NULL)
    {
        strcpy(this->module_name, "LJ");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }
#ifdef USE_CUDA
    // S3：tile kernel 的 lane 映射假设 32-lane warp（CUDA 上恒成立），
    // CUDA 构建默认启用，可用 mdin 命令 "LJ use_tile = 0" 回退旧 kernel
    use_tile = CONTROLLER::device_warp == 32;
#endif
    if (controller->Command_Exist(this->module_name, "use_tile"))
    {
        use_tile = controller->Get_Bool(this->module_name, "use_tile",
                                        "LENNARD_JONES_INFORMATION::Initial");
#ifndef USE_CUDA
        if (use_tile)
        {
            controller->printf(
                "    LJ tile force kernel is CUDA-only; falling back to the "
                "legacy kernel\n");
            use_tile = false;
        }
#endif
    }
    controller->printf("START INITIALIZING LENNADR JONES INFORMATION:\n");
    const auto& lj = Xponge::system.classical_force_field.lj;
    Xponge::LennardJones local_lj;
    const Xponge::LennardJones* lj_to_use = NULL;
    if (module_name == NULL)
    {
        lj_to_use = &lj;
    }
    else if (controller->Command_Exist(this->module_name, "in_file"))
    {
        Xponge::Native_Load_LJ(&local_lj, controller,
                               Xponge::Load_Get_Atom_Numbers(&Xponge::system),
                               this->module_name);
        lj_to_use = &local_lj;
    }
    if (lj_to_use != NULL)
    {
        if (lj_to_use->atom_type.size() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_INFORMATION::Initial",
                "Reason:\n\tLJ atom count cannot be represented by the "
                "runtime index type\n");
            return;
        }
        atom_numbers = static_cast<int>(lj_to_use->atom_type.size());
        atom_type_numbers = lj_to_use->atom_type_numbers;
        if (atom_type_numbers < 0 || atom_type_numbers > 65535)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_INFORMATION::Initial",
                "Reason:\n\tLJ atom type count %d is outside the "
                "representable range [0, 65535]\n",
                atom_type_numbers);
            return;
        }
        const std::uint64_t type_count =
            static_cast<std::uint64_t>(atom_type_numbers);
        const std::uint64_t expected_pair_count =
            type_count * (type_count + 1ULL) / 2ULL;
        if (expected_pair_count >
                static_cast<std::uint64_t>(std::numeric_limits<int>::max()) ||
            lj_to_use->pair_A.size() != expected_pair_count ||
            lj_to_use->pair_B.size() != expected_pair_count ||
            (atom_numbers > 0 && atom_type_numbers == 0))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "LENNARD_JONES_INFORMATION::Initial",
                "Reason:\n\tLJ table shape is inconsistent: atoms=%d, "
                "types=%d, expected pairs=%llu, A/B pairs=%llu/%llu\n",
                atom_numbers, atom_type_numbers,
                static_cast<unsigned long long>(expected_pair_count),
                static_cast<unsigned long long>(lj_to_use->pair_A.size()),
                static_cast<unsigned long long>(lj_to_use->pair_B.size()));
            return;
        }
        pair_type_numbers = static_cast<int>(expected_pair_count);
        for (int atom = 0; atom < atom_numbers; atom++)
        {
            const int type = lj_to_use->atom_type[atom];
            if (type < 0 || type >= atom_type_numbers)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand,
                    "LENNARD_JONES_INFORMATION::Initial",
                    "Reason:\n\tLJ atom %d has type %d outside [0, %d)\n", atom,
                    type, atom_type_numbers);
                return;
            }
        }
        for (int pair = 0; pair < pair_type_numbers; pair++)
        {
            if (!Float_Memory_Is_Finite(lj_to_use->pair_A.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_to_use->pair_A.data() +
                                                pair) ||
                !Float_Memory_Is_Finite(lj_to_use->pair_B.data() + pair) ||
                !Float_Memory_Is_Zero_Or_Normal(lj_to_use->pair_B.data() +
                                                pair))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorConflictingCommand,
                    "LENNARD_JONES_INFORMATION::Initial",
                    "Reason:\n\tLJ pair %d contains a non-finite or "
                    "subnormal coefficient\n",
                    pair);
                return;
            }
        }
    }
    if (atom_numbers > 0)
    {
        controller->printf("    atom_numbers is %d\n", atom_numbers);
        controller->printf("    atom_LJ_type_number is %d\n",
                           atom_type_numbers);
        LJ_Malloc();

        for (int i = 0; i < pair_type_numbers; i++)
        {
            h_LJ_A[i] = lj_to_use->pair_A[i];
            h_LJ_B[i] = lj_to_use->pair_B[i];
        }
        for (int i = 0; i < atom_numbers; i++)
        {
            h_atom_LJ_type[i] = lj_to_use->atom_type[i];
        }
        Parameter_Host_To_Device();
        is_initialized = 1;
    }
    if (is_initialized)
    {
        this->cutoff = cutoff;
        Device_Malloc_Safely((void**)&crd_with_LJ_parameters,
                             sizeof(VECTOR_LJ) * atom_numbers);
        Launch_Device_Kernel(
            Copy_LJ_Type_To_New_Crd,
            (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers,
            crd_with_LJ_parameters, d_atom_LJ_type);
        controller->printf("    Start initializing long range LJ correction\n");
        // 全对求和 Σ_i Σ_j B[type_i, type_j] 等于按类型直方图的
        // Σ_a count_a · Σ_b count_b · B[pair(a,b)]，后者按固定顺序双精度
        // 累加，结果确定；此前的全对 kernel 复杂度为 O(N²) 且依赖 float
        // 原子加顺序，本身就有运行间波动。
        std::vector<int64_t> type_count(atom_type_numbers, 0);
        for (int i = 0; i < atom_numbers; i++)
        {
            type_count[h_atom_LJ_type[i]] += 1;
        }
        double c6_sum = 0.0;
        for (int itype = 0; itype < atom_type_numbers; itype++)
        {
            if (type_count[itype] == 0) continue;
            double inner_sum = 0.0;
            for (int jtype = 0; jtype < atom_type_numbers; jtype++)
            {
                if (type_count[jtype] == 0) continue;
                inner_sum += static_cast<double>(type_count[jtype]) *
                             static_cast<double>(
                                 h_LJ_B[Get_LJ_Type(itype, jtype)]);
            }
            c6_sum += static_cast<double>(type_count[itype]) * inner_sum;
        }
        long_range_factor = static_cast<float>(c6_sum);
        printf("        Total C6 factor is %e\n", long_range_factor);

        long_range_factor *=
            -2.0f / 3.0f * CONSTANT_Pi / cutoff / cutoff / cutoff / 6.0f;
        controller->printf("        long range correction factor is: %e\n",
                           long_range_factor);
        controller->printf("    End initializing long range LJ correction\n");
    }
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial("LJ_short", "%.2f");
        controller->Step_Print_Initial("LJ_long", "%.2f");
        controller->Step_Print_Initial("LJ", "%.2f");
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
        controller->printf("    LJ tile force kernel (use_tile): %s\n",
                           use_tile ? "on" : "off");
    }
    controller->printf("END INITIALIZING LENNADR JONES INFORMATION\n\n");
}

static __global__ void get_local_device(int* atom_local, int local_atom_numbers,
                                        int ghost_numbers,
                                        int global_atom_numbers,
                                        int* d_atom_LJ_type,
                                        const float* charge, float4* d_lj_crd_q,
                                        int2* d_lj_type_g,
                                        int* invalid_local_index)
{
    SIMPLE_DEVICE_FOR(i, local_atom_numbers + ghost_numbers)
    {
        int atom_i = atom_local[i];
        if (atom_i < 0 || atom_i >= global_atom_numbers)
        {
            atomicExch(invalid_local_index, i);
        }
        else
        {
            d_lj_type_g[i].x = d_atom_LJ_type[atom_i];
            d_lj_type_g[i].y = atom_i;
            d_lj_crd_q[i].w = charge[i];
        }
    }
}

void LENNARD_JONES_INFORMATION::Get_Local(int* atom_local,
                                          int local_atom_numbers,
                                          int ghost_numbers,
                                          const float* charge)
{
    if (!is_initialized) return;
    local_metadata_is_ready = false;
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > atom_numbers ||
        ghost_numbers > atom_numbers - local_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "LENNARD_JONES_INFORMATION::Get_Local",
            "Reason:\n\t%s received invalid local/ghost atom counts %d/%d "
            "for %d global atoms\n",
            module_name, local_atom_numbers, ghost_numbers, atom_numbers);
        return;
    }
    this->local_atom_numbers = local_atom_numbers;
    this->ghost_numbers = ghost_numbers;
    const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
    if (local_coordinate_numbers == 0)
    {
        local_metadata_is_ready = true;
        return;
    }
    deviceMemset(d_local_metadata_error, -1, sizeof(int));
    Launch_Device_Kernel(
        get_local_device,
        (local_coordinate_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, atom_local, local_atom_numbers,
        ghost_numbers, atom_numbers, d_atom_LJ_type, charge, d_lj_crd_q,
        d_lj_type_g, d_local_metadata_error);
    int invalid_local_index = -1;
    deviceMemcpy(&invalid_local_index, d_local_metadata_error, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (invalid_local_index >= 0)
    {
        int invalid_global_atom = -1;
        deviceMemcpy(&invalid_global_atom, atom_local + invalid_local_index,
                     sizeof(int), deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "LENNARD_JONES_INFORMATION::Get_Local",
            "Reason:\n\t%s local coordinate %d maps to global atom %d "
            "outside [0, %d)\n",
            module_name, invalid_local_index, invalid_global_atom,
            atom_numbers);
        return;
    }
    local_metadata_is_ready = true;
}

bool LENNARD_JONES_INFORMATION::Validate_Local_State(const char* error_by,
                                                     int global_atom_numbers,
                                                     int local_atom_numbers,
                                                     int ghost_numbers,
                                                     int solvent_numbers)
{
    if (global_atom_numbers != atom_numbers ||
        local_atom_numbers != this->local_atom_numbers ||
        ghost_numbers != this->ghost_numbers || !local_metadata_is_ready ||
        solvent_numbers < 0 || solvent_numbers > local_atom_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s local nonbond metadata mismatch: call has "
            "global/local/ghost/solvent counts %d/%d/%d/%d, initialized "
            "state has %d/%d/%d and ready=%d\n",
            module_name, global_atom_numbers, local_atom_numbers, ghost_numbers,
            solvent_numbers, atom_numbers, this->local_atom_numbers,
            this->ghost_numbers, static_cast<int>(local_metadata_is_ready));
        return false;
    }
    return true;
}

void LENNARD_JONES_INFORMATION::Reset_Pair_Overlap_Error()
{
#ifndef GPU_ARCH_NAME
    deviceMemset(d_pair_overlap_error, 0, 3 * sizeof(int));
#endif
}

bool LENNARD_JONES_INFORMATION::Check_Pair_Overlap_Error(const char* error_by)
{
#ifndef GPU_ARCH_NAME
    int overlap_error[3] = {0, -1, -1};
    deviceMemcpy(overlap_error, d_pair_overlap_error, sizeof(overlap_error),
                 deviceMemcpyDeviceToHost);
    if (overlap_error[0] != PairwiseInteraction::PAIR_COMPONENT_NONE)
    {
        const char* component =
            overlap_error[0] ==
                    (PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES |
                     PairwiseInteraction::PAIR_COMPONENT_COULOMB)
                ? "LJ and Coulomb components"
            : overlap_error[0] &
                    PairwiseInteraction::PAIR_COMPONENT_LENNARD_JONES
                ? "LJ component"
                : "Coulomb component";
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, error_by,
            "Reason:\n\t%s global atoms %d %d overlap exactly with active "
            "%s; an active hard nonbond pair has undefined inverse-distance "
            "energy and force\n",
            module_name, overlap_error[1], overlap_error[2], component);
        return true;
    }
#else
    (void)error_by;
#endif
    return false;
}

static __global__ void Long_Range_Virial_Correction(LTMatrix3* d_virial,
                                                    const float factor)
{
    d_virial[0].a11 += factor;
    d_virial[0].a22 += factor;
    d_virial[0].a33 += factor;
}

void LENNARD_JONES_INFORMATION::Long_Range_Correction(int need_pressure,
                                                      LTMatrix3* d_virial,
                                                      int need_potential,
                                                      float* d_potential,
                                                      const float volume)
{
    if (is_initialized && CONTROLLER::PP_MPI_rank == 0)
    {
        if (need_pressure)
        {
            Launch_Device_Kernel(Long_Range_Virial_Correction, 1, 1, 0, 0,
                                 d_virial, 2 * long_range_factor / volume);
        }
        if (need_potential)
        {
            Launch_Device_Kernel(device_add, 1, 1, 0, 0, d_potential,
                                 long_range_factor / volume);

            h_LJ_long_energy = long_range_factor / volume;
        }
    }
}

void LENNARD_JONES_INFORMATION::Parameter_Host_To_Device()
{
    Device_Malloc_And_Copy_Safely((void**)&d_atom_LJ_type, h_atom_LJ_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_A, h_LJ_A,
                                  sizeof(float) * pair_type_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_B, h_LJ_B,
                                  sizeof(float) * pair_type_numbers);
    Device_Malloc_And_Copy_Safely((void**)&d_LJ_energy_sum, h_LJ_energy_atom,
                                  sizeof(float));
    Device_Malloc_Safely((void**)&d_LJ_energy_atom,
                         sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_lj_crd_q, sizeof(float4) * atom_numbers);
    Device_Malloc_Safely((void**)&d_lj_type_g, sizeof(int2) * atom_numbers);
    Device_Malloc_Safely((void**)&d_local_metadata_error, sizeof(int));
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_pair_overlap_error, 3 * sizeof(int));
#endif
}

void LENNARD_JONES_INFORMATION::LJ_PME_Direct_Force_With_Atom_Energy_And_Virial(
    const int atom_numbers, const int local_atom_numbers,
    const int solvent_numbers, const int ghost_numbers, const VECTOR* crd,
    const float* charge, VECTOR* frc, const LTMatrix3 cell,
    const LTMatrix3 rcell, const ATOM_GROUP* nl, const float pme_beta,
    const int need_atom_energy, float* atom_energy, const int need_virial,
    LTMatrix3* atom_virial, float* atom_direct_pme_energy,
    const LJ_TILE_SET* tile_set)
{
    if (is_initialized)
    {
        if (!Validate_Local_State(
                "LENNARD_JONES_INFORMATION::"
                "LJ_PME_Direct_Force_With_Atom_Energy_And_Virial",
                atom_numbers, local_atom_numbers, ghost_numbers,
                solvent_numbers))
        {
            return;
        }
        Launch_Device_Kernel(
            Repack_LJ_Crd,
            (this->atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL,
            this->local_atom_numbers + this->ghost_numbers, crd, d_lj_crd_q);
        if (need_atom_energy)
        {
            deviceMemset(atom_direct_pme_energy, 0,
                         sizeof(float) * this->atom_numbers);
            deviceMemset(d_LJ_energy_atom, 0,
                         sizeof(float) * this->atom_numbers);
        }

        if (atom_numbers == 0 || local_atom_numbers == 0) return;

        Reset_Pair_Overlap_Error();

#ifdef USE_CUDA
        // S3 tile 路径：tile 表覆盖全部 local i 的半表对（溶剂对也在内，
        // 调用方必须同时停用溶剂 dispatch，否则水-水双计），solvent_numbers
        // 在此不起作用。tile_set 为空（use_tile=0 / 邻居表未产出 tile 表 /
        // SITS 选择性施加）时走下方旧 kernel，主/溶剂分工保持不变。
        if (tile_set != NULL && use_tile && tile_set->tiles != NULL &&
            tile_set->cluster_atoms != NULL && tile_set->cluster_flags != NULL &&
            tile_set->tile_sorted != NULL)
        {
            // cluster 序打包数组：容量随 cluster 槽数惰性增长（溢出恢复扩容
            // 重建后槽数可能变大），每步从 atom 序数组重排一次
            if (tile_set->cluster_atom_slots > lj_cluster_array_capacity)
            {
                if (d_lj_cluster_crd_q != NULL)
                    Free_Single_Device_Pointer((void**)&d_lj_cluster_crd_q);
                if (d_lj_cluster_type_g != NULL)
                    Free_Single_Device_Pointer((void**)&d_lj_cluster_type_g);
                lj_cluster_array_capacity = 0;
                if (!Device_Malloc_Safely(
                        (void**)&d_lj_cluster_crd_q,
                        sizeof(float4) * tile_set->cluster_atom_slots) ||
                    !Device_Malloc_Safely(
                        (void**)&d_lj_cluster_type_g,
                        sizeof(int2) * tile_set->cluster_atom_slots))
                {
                    controller->Throw_SPONGE_Error(
                        spongeErrorMallocFailed,
                        "LENNARD_JONES_INFORMATION::"
                        "LJ_PME_Direct_Force_With_Atom_Energy_And_Virial");
                    return;
                }
                lj_cluster_array_capacity = tile_set->cluster_atom_slots;
            }
            if (tile_set->tile_numbers > 0 && tile_set->cluster_atom_slots > 0)
            {
                Launch_Device_Kernel(
                    Gather_LJ_Cluster_Data,
                    (tile_set->cluster_atom_slots +
                     CONTROLLER::device_max_thread - 1) /
                        CONTROLLER::device_max_thread,
                    CONTROLLER::device_max_thread, 0, NULL,
                    tile_set->cluster_atom_slots, tile_set->cluster_atoms,
                    d_lj_crd_q, d_lj_type_g, d_lj_cluster_crd_q,
                    d_lj_cluster_type_g);
                // tile kernel 固定 256 线程块（8 warp），每 warp 连续消费
                // LJ_TILE_WARP_SPAN 个分组序 tile
                dim3 tileBlock(256);
                dim3 tileGrid(
                    (tile_set->tile_numbers +
                     (256 / 32) * LJ_TILE_WARP_SPAN - 1) /
                    ((256 / 32) * LJ_TILE_WARP_SPAN));
                auto ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                    true, false, false, true>;
                if (!need_atom_energy && !need_virial)
                {
                    // 消融测时实例（仅诊断；生产 ABLATE=0 编译期消除）。
                    // 仅对主力实例 <true,false,false,true> 提供
                    switch (LJ_Force_Ablate_Mode())
                    {
                        case 1:
                            ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                                true, false, false, true, 1>;
                            break;
                        case 2:
                            ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                                true, false, false, true, 2>;
                            break;
                        case 3:
                            ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                                true, false, false, true, 3>;
                            break;
                        case 4:
                            ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                                true, false, false, true, 4>;
                            break;
                        case 5:
                            ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                                true, false, false, true, 5>;
                            break;
                        case 6:
                            ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                                true, false, false, true, 6>;
                            break;
                        case 7:
                            ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                                true, false, false, true, 7>;
                            break;
                        default:
                            break;
                    }
                }
                else if (need_atom_energy && !need_virial)
                {
                    ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                        true, true, false, true>;
                }
                else if (!need_atom_energy && need_virial)
                {
                    ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                        true, false, true, true>;
                }
                else
                {
                    ft = Lennard_Jones_And_Direct_Coulomb_Tile_Device<
                        true, true, true, true>;
                }
                Launch_Device_Kernel(
                    ft, tileGrid, tileBlock, 0, NULL, tile_set->tile_numbers,
                    tile_set->tiles, tile_set->tile_sorted,
                    tile_set->cluster_atoms, tile_set->cluster_flags,
                    d_lj_cluster_crd_q, d_lj_cluster_type_g, cell, rcell,
                    d_LJ_A, d_LJ_B, cutoff, frc, pme_beta, atom_energy,
                    atom_virial, atom_direct_pme_energy, d_LJ_energy_atom,
                    d_pair_overlap_error);
            }
            Check_Pair_Overlap_Error(
                "LENNARD_JONES_INFORMATION::"
                "LJ_PME_Direct_Force_With_Atom_Energy_And_Virial");
            return;
        }
#endif

        dim3 blockSize = {
            CONTROLLER::device_warp,
            CONTROLLER::device_max_thread / CONTROLLER::device_warp};
        dim3 gridSize = (local_atom_numbers + blockSize.y - 1) / blockSize.y;
        auto f =
            Lennard_Jones_And_Direct_Coulomb_Device<true, false, false, true>;
        if (!need_atom_energy && !need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, false, false,
                                                        true>;
        }
        else if (need_atom_energy && !need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, true, false,
                                                        true>;
        }
        else if (!need_atom_energy && need_virial)
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, false, true,
                                                        true>;
        }
        else
        {
            f = Lennard_Jones_And_Direct_Coulomb_Device<true, true, true, true>;
        }
        Launch_Device_Kernel(
            f, gridSize, blockSize, 0, NULL, local_atom_numbers,
            solvent_numbers, nl, d_lj_crd_q, d_lj_type_g, cell, rcell,
            d_LJ_A, d_LJ_B, cutoff, frc, pme_beta, atom_energy, atom_virial,
            atom_direct_pme_energy, d_LJ_energy_atom, d_pair_overlap_error);
        Check_Pair_Overlap_Error(
            "LENNARD_JONES_INFORMATION::"
            "LJ_PME_Direct_Force_With_Atom_Energy_And_Virial");
    }
}

void LENNARD_JONES_INFORMATION::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized || CONTROLLER::MPI_rank >= CONTROLLER::PP_MPI_size)
        return;
    Sum_Of_List(d_LJ_energy_atom, d_LJ_energy_sum, atom_numbers);
    deviceMemcpy(&h_LJ_energy_sum, d_LJ_energy_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &h_LJ_energy_sum, 1, MPI_FLOAT, MPI_SUM,
                  CONTROLLER::pp_comm);
#endif
    controller->Step_Print("LJ_short", h_LJ_energy_sum);
    controller->Step_Print("LJ_long", h_LJ_long_energy);
    controller->Step_Print("LJ", h_LJ_energy_sum + h_LJ_long_energy, true);
}
