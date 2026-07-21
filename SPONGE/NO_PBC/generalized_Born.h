#pragma once
#include "../common.h"
#include "../control.h"

// 用于记录与计算GB相关的信息
// Reference:
// 1.Theory and application of the Generalized Born solvation model in
// macromolecules simulations
// DOI:10.1002/1097-0282(2000)56:4<275::AID-BIP10024>3.0.CO;2-E
// 2.Parametrized models of aqueous free energies of solvation based on pairwise
// descreening of solute atomic charges from a dielectric medium
// DOI:10.1021/jp961710n
// equation 9-11
struct GENERALIZED_BORN_INFORMATION
{
    CONTROLLER* controller = NULL;
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260216;

    int atom_numbers = 0;  // 原子数

    float* h_GB_energy_atom = NULL;  // 每个原子的GB的能量
    float h_GB_energy_sum = 0;       // 所有原子的GB能量和
    float* d_GB_energy_atom = NULL;  // 每个原子的GB的能量
    float* d_GB_energy_sum = NULL;   // 所有原子的GB能量和

    float* h_GB_self_radius = NULL;   // 自身半径(radii - offset)
    float* d_GB_self_radius = NULL;   // 自身半径(radii - offset)
    float* h_GB_other_radius = NULL;  // 屏蔽半径scaler * (radii - offset)
    float* d_GB_other_radius = NULL;  // 屏蔽半径scaler * (radii - offset)

    float* d_GB_effective_radius = NULL;  // 有效半径
    float* d_dE_da = NULL;                // 能量对effective_radius的导数
    int* d_pair_overlap_error = NULL;
    int* d_effective_radius_error = NULL;

    float relative_dielectric_constant = 78.5;
    float radii_offset = 0.09;
    float radii_cutoff = 25.0;
    float radii_cutoff_square = 625.0;

    // 初始化
    void Initial(CONTROLLER* controller, int expected_atom_numbers,
                 float default_radii_cutoff,
                 const char* module_name = NULL);
    // 清除内存
    void Clear();
    // 分配内存
    void Malloc();
    void Reset_Pair_Overlap_Error();
    bool Check_Pair_Overlap_Error(const char* error_by);
    void Reset_Effective_Radius_Error();
    bool Check_Effective_Radius_Error(const char* error_by);

    bool Get_Effective_Born_Radius(int atom_numbers, const VECTOR* crd);

    // 可以根据外界传入的need_atom_energ选择性计算能量
    void GB_Force_With_Atom_Energy(const int atom_numbers, const VECTOR* crd,
                                   const float* charge, VECTOR* frc,
                                   float* atom_energy);

    void Step_Print(CONTROLLER* controller);
};
