#pragma once
#include "../common.h"
#include "../control.h"

// 用于记录与计算MC控压相关的信息
struct MC_BAROSTAT_INFORMATION
{
    char module_name[CHAR_LENGTH_MAX];
    CONTROLLER* controller = NULL;
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260216;

    enum
    {
        NO = 0,
        XY = 1,
        XZ = 2,
        YZ = 3,
        XYZ = 4
    } couple_dimension = XYZ;
    int only_direction = 0;
    int surface_number = 0;
    float surface_tension = 0.0f;

    int xyz = 0;                // 这次更新的是xyz的哪一个
    int update_interval = 100;  // 每多少步进行一次MC尝试
    int check_interval =
        20;  // 每多少次MC尝试以后进行一次DeltaV_max取值的检查（使得MC接受概率在40%~50%
    int accept = 0;      // 接受与否标志，1表示接受，0表示拒绝
    VECTOR* crd_backup;  // 备份坐标，以便还原
    // 成功率记录相关

    int total_count[3] = {0, 0, 0};    // 总共进行的MC尝试次数
    int accep_count[3] = {0, 0, 0};    // 接受的MC尝试次数
    float accept_rate[3] = {0, 0, 0};  // 接受率
    float accept_rate_low;             // 接受率的低限
    float accept_rate_high;            // 接受率的高限

    // 每次允许变化最大的边长
    float Delta_Box_Length_Max[3];

    //// 温度常数， 等于1 / k_B / T
    // float beta;

    // 设定的压强
    float p0;

    // 每次尝试变化的体积
    VECTOR Delta_Box_Length;
    VECTOR New_Box_Length;
    double DeltaV;
    float newV;
    double VDevided;
    float DeltaS;

    // 坐标系数因子
    LTMatrix3 g;  // 盒子缩放矩阵

    // 概率项
    float energy_old;
    float energy_new;
    double extra_term;
    double accept_possibility;

    // scale coordinate atomically
    void Scale_Coordinate_Atomically(int atom_numbers, VECTOR* crd);

    // 程序中使用的单位与物理的换算
    //  单位: mN/m -> kcal/mol/A^2
    //  1 / ( (1 kcal/mol) * (4.184074e3 J/kcal) / (6.023e23 mol^-1) * (1e20
    //  m^2/A^2) * (1e3 N/mN) ) 程序/(kcal/mol/A^2 ) = 物理/(mN/m) *
    //  TENSION_UNIT_FACTOR
    const float TENSION_UNIT_FACTOR = 1.439506e-3f;

    // 初始化
    void Initial(CONTROLLER* controller, int atom_numbers,
                 float target_pressure, VECTOR boxlength, LTMatrix3 cell,
                 const char* module_name = NULL);

    // 输入盒子信息，得到坐标变化因子
    void Volume_Change_Attempt(VECTOR boxlength, float dt);

    bool Will_Attempt(int steps) const;
    LTMatrix3 Get_Exact_Reverse_G(float dt) const;

    // 判断是否接受
    int Check_MC_Barostat_Accept();

    // 对猜测值根据接受率更新
    void Delta_Box_Length_Max_Update();

    // 判断步数条件，决定是否让计算力的时候也计算能量
    void Ask_For_Calculate_Potential(int steps, int* need_potential);
};
