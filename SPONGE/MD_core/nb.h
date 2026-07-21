#pragma once

struct non_bond_information
{
    MD_INFORMATION* md_info = NULL;
    float cutoff = 10.0;
    float skin = 2.0;
    int excluded_atom_numbers;   // 排除表总长
    int* d_excluded_list_start;  // 记录每个原子的剔除表起点
    int* d_excluded_list;        // 剔除表
    int* d_excluded_numbers;     // 记录每个原子需要剔除的原子个数
    int* h_excluded_list_start;  // 记录每个原子的剔除表起点
    int* h_excluded_list;        // 剔除表
    int* h_excluded_numbers;     // 记录每个原子需要剔除的原子个数
    void Initial(CONTROLLER* controller, MD_INFORMATION* md_info);
    void Excluded_List_Reform(
        CONTROLLER* controller,
        int atom_numbers);  // 为了区域分解后近邻表计算正确性，补齐跨残基的近邻表
};
