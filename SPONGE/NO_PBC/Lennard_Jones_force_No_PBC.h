#pragma once
#include "../common.h"
#include "../control.h"

// 用于记录与计算LJ相关的信息
struct LENNARD_JONES_NO_PBC_INFORMATION
{
    CONTROLLER* controller = NULL;
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260216;

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
    float* d_LJ_energy_sum = NULL;   // 所有原子的LJ能量和
    int* d_pair_overlap_error = NULL;

    // 初始化
    void Initial(CONTROLLER* controller, const char* module_name = NULL);
    // 清除内存
    void Clear();
    // 分配内存
    void LJ_Malloc();
    // 参数传到GPU上
    void Parameter_Host_To_Device();
    void Reset_Pair_Overlap_Error();
    bool Check_Pair_Overlap_Error(const char* error_by);

    // 可以根据外界传入的need_atom_energ选择性计算能量
    void LJ_Force_With_Atom_Energy(const int atom_numbers, const VECTOR* crd,
                                   VECTOR* frc, const int need_atom_energy,
                                   float* atom_energy,
                                   const int* excluded_list_start,
                                   const int* excluded_list,
                                   const int* excluded_atom_numbers);

    void Step_Print(CONTROLLER* controller);
};
