#pragma once
#include "../common.h"
#include "../control.h"

// 用于记录与计算LJ相关的信息
struct COULOMB_FORCE_NO_PBC_INFORMATION
{
    CONTROLLER* controller = NULL;
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260216;

    // E_lj = qi * qj / r;
    // F_lj = qi * qj / r ^ 2;
    int atom_numbers = 0;  // 原子数

    float* h_Coulomb_energy_atom = NULL;  // 每个原子的Coulomb的能量
    float h_Coulomb_energy_sum = 0;       // 所有原子的Coulomb能量和
    float* d_Coulomb_energy_atom = NULL;  // 每个原子的Coulomb的能量
    float* d_Coulomb_energy_sum = NULL;   // 所有原子的Coulomb能量和
    int* d_pair_overlap_error = NULL;

    // 初始化
    void Initial(CONTROLLER* controller, int atom_numbers,
                 const char* module_name = NULL);
    // 分配内存
    void Malloc();
    // 清除内存
    void Clear();
    void Reset_Pair_Overlap_Error();
    bool Check_Pair_Overlap_Error(const char* error_by);

    // 可以根据外界传入的need_atom_energ选择性计算能量
    void Coulomb_Force_With_Atom_Energy(const int atom_numbers,
                                        const VECTOR* crd, const float* charge,
                                        VECTOR* frc, const int need_atom_energy,
                                        float* atom_energy,
                                        const int* excluded_list_start,
                                        const int* excluded_list,
                                        const int* excluded_atom_numbers);

    void Step_Print(CONTROLLER* controller);
};
