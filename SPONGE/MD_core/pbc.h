#pragma once

struct periodic_box_condition_information
{
    MD_INFORMATION* md_info = NULL;
    CONTROLLER* controller = NULL;
    bool is_initialized = false;
    bool pbc = false;

    // 初始化时的盒子，用于判断是否变化合适
    LTMatrix3 cell0;
    // Input/reference cell used for volume, density, stress normalization and
    // explicitly scaled coordinates even when interactions are nonperiodic.
    LTMatrix3 reference_cell;
    // 盒子的下三角矩阵
    LTMatrix3 cell;
    // 倒空间的盒子的下三角矩阵
    LTMatrix3 rcell;
    // 盒子的质量倒数
    LTMatrix3 cell_mass_inverse;
    // 盒子的动量
    LTMatrix3 cell_momentum;

    void Initial(CONTROLLER* controller, MD_INFORMATION* md_info);
    void No_PBC_Check(CONTROLLER* controller);
    void PBC_Check();
    LTMatrix3 Get_Cell(VECTOR box_length, VECTOR box_angle);
    void Update_Box(LTMatrix3 g);
    void Update_Box_From_Input(VECTOR box_length, VECTOR box_angle);
    bool Check_Change_Large();
};
