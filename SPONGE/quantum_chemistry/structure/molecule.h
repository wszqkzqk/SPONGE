#pragma once

#include "../../common.h"
#include "elements.hpp"

// 量化壳层
struct QC_SHELL
{
    int l;
    std::vector<float> exps;
    std::vector<float> coeffs;
};

// 量化分子
struct QC_MOLECULE
{
    // 原子数
    int natm;
    // 电子总数
    int nelectron;
    // 笛卡尔 AO 总数
    int nao_cart;
    // 笛卡尔 AO 总数平方（初始化时已验证可用 int 索引）
    int nao_cart2 = 0;
    // 球谐 AO 总数
    int nao_sph;
    // 是否使用球谐基（1: 是，0: 否）
    int is_spherical;
    // 当前有效 AO 总数（由 is_spherical 选择球谐或笛卡尔）
    int nao = 0;
    // 当前有效 AO 总数平方
    int nao2 = 0;
    // 壳层总数
    int nbas;
    // 分子总电荷
    int charge;
    // 自旋多重度
    int multiplicity;
    // 壳层数据
    std::vector<QC_SHELL> shells;

    // Immutable chemical identity. ECP setup must never modify this array;
    // element-specific bases and DFT-grid metadata are selected from it.
    std::vector<int> h_atomic_numbers;
    int* d_atomic_numbers = NULL;

    // Effective nuclear charge used by Coulomb/integral kernels. For an ECP
    // atom this is atomic_number - n_core and therefore is not an element ID.
    std::vector<int> h_Z;
    int* d_Z = NULL;

    int Atomic_Number(int atom_index) const
    {
        return h_atomic_numbers.at(static_cast<std::size_t>(atom_index));
    }

    int Effective_Nuclear_Charge(int atom_index) const
    {
        return h_Z.at(static_cast<std::size_t>(atom_index));
    }

    // 每个壳层中心坐标
    std::vector<VECTOR> h_centers;
    VECTOR* d_centers = NULL;

    // 每个壳层的角量子数
    std::vector<int> h_l_list;
    int* d_l_list = NULL;

    // 所有壳层指数参数拼接
    std::vector<float> h_exps;
    float* d_exps = NULL;

    // 所有壳层收缩系数拼接
    std::vector<float> h_coeffs;
    float* d_coeffs = NULL;

    // 每个壳层的起始偏移
    std::vector<int> h_shell_offsets;
    int* d_shell_offsets = NULL;

    // 每个壳层的高斯函数数量
    std::vector<int> h_shell_sizes;
    int* d_shell_sizes = NULL;

    // 每个壳层的原子轨道起始偏移（笛卡尔）
    std::vector<int> h_ao_offsets;
    int* d_ao_offsets = NULL;

    // 每个壳层的原子轨道起始偏移（球谐/有效 AO）
    std::vector<int> h_ao_offsets_sph;
    int* d_ao_offsets_sph = NULL;

    // 原子参数数组（积分内核使用）
    std::vector<int> h_atm;
    int* d_atm = NULL;

    // 基函数参数数组（积分内核使用）
    std::vector<int> h_bas;
    int* d_bas = NULL;

    // 环境参数数组（积分内核使用）
    std::vector<float> h_env;
    float* d_env = NULL;

    // 原子坐标 (按原子索引, ECP/梯度内核使用)
    std::vector<VECTOR> h_atom_coords;  // [natm]
    VECTOR* d_atom_coords = NULL;

    // ECP 数据
    bool has_ecp = false;
    // 每个原子的核心电子数和 l_max (-1 = 无 ECP)
    std::vector<int> h_ecp_n_core;  // [natm]
    std::vector<int> h_ecp_l_max;   // [natm]
    int* d_ecp_l_max = NULL;
    // 扁平化通道数据 (device kernel 使用)
    int ecp_total_channels = 0;
    int ecp_total_terms = 0;
    std::vector<float> h_ecp_d;              // 所有 d_k 拼接
    std::vector<float> h_ecp_zeta;           // 所有 ζ_k 拼接
    std::vector<int> h_ecp_n;                // 所有 n_k 拼接
    std::vector<int> h_ecp_l;                // 每个 channel 的 l 值
    std::vector<int> h_ecp_channel_offsets;  // 每个 channel 在扁平数组中的偏移
    std::vector<int> h_ecp_channel_sizes;    // 每个 channel 的 term 数
    std::vector<int>
        h_ecp_atom_channel_range;  // [natm+1] 原子 i 的 channel 范围
    // Device 指针
    float* d_ecp_d = NULL;
    float* d_ecp_zeta = NULL;
    int* d_ecp_n = NULL;
    int* d_ecp_l = NULL;
    int* d_ecp_channel_offsets = NULL;
    int* d_ecp_channel_sizes = NULL;
    int* d_ecp_atom_channel_range = NULL;
};
