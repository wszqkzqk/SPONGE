#pragma once

#include <cstdint>

#include "../common.h"
#include "../control.h"
#include "gradient/grad_workspace.h"
#include "scf/eigensolver_policy.hpp"
#include "structure/cart2sph.h"
#include "structure/dft.h"
#include "structure/ecp.h"
#include "structure/integral_tasks.h"
#include "structure/matrix.h"
#include "structure/method.h"
#include "structure/molecule.h"
#include "structure/scf_workspace.h"

#define ONE_E_BATCH_SIZE 4096
#define QC_GRAD_ERI_THREADS 64
#define QC_GRAD_GAMMA_POOL_BLOCKS 64
#define QC_GRAD_GAMMA_POOL_SLOTS \
    (QC_GRAD_ERI_THREADS * QC_GRAD_GAMMA_POOL_BLOCKS)
#define PI_25 17.4934183276248628469f
#define HR_BASE_MAX 17
#define HR_SIZE_MAX 83521
#define ONEE_MD_BASE 9
#define ONEE_MD_IDX(t, u, v, n) \
    ((((t) * ONEE_MD_BASE + (u)) * ONEE_MD_BASE + (v)) * ONEE_MD_BASE + (n))
#define ERI_BATCH_SIZE 128
#define QC_BOUNDS_POOL_SLOTS 1024
#define MAX_CART_SHELL 15
#define MAX_SHELL_ERI \
    (MAX_CART_SHELL * MAX_CART_SHELL * MAX_CART_SHELL * MAX_CART_SHELL)

struct QUANTUM_CHEMISTRY
{
   public:
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    bool periodic_boundary = true;
    int last_modify_date = 20260216;
    int atom_numbers = 0;

    float scf_energy = 0.0f;
    FILE* scf_output_file = NULL;

    // 本地原子映射
    std::vector<int> atom_local;
    int* d_atom_local = NULL;

    // 计算方法
    QC_METHOD method = QC_METHOD::HF;
    // 初始猜测
    QC_INITIAL_GUESS initial_guess = QC_INITIAL_GUESS::SAP;
    bool need_initial_guess = true;
    // DFT信息
    QC_DFT dft;
    // 分子信息
    QC_MOLECULE mol;
    // 积分任务
    QC_INTEGRAL_TASKS task_ctx;
    // SCF计算内容
    QC_SCF_WORKSPACE scf_ws;
    // 梯度工作空间
    QC_GRAD_WORKSPACE grad_ws;

    BLAS_HANDLE blas_handle{};
    SOLVER_HANDLE solver_handle{};

    // 笛卡尔基组转球形基组
    QC_CARTESIAN_TO_SPHERICAL cart2sph;

    int need_gradient = 1;  // 是否计算梯度/力 (qc_need_gradient)

    // 外部入口
    void Initial(CONTROLLER* controller, const int atom_numbers,
                 const VECTOR* crd, const VECTOR box_length,
                 const char* module_name = NULL,
                 std::uint64_t coordinate_generation = 0);
    void Solve_SCF(const VECTOR* crd, const VECTOR box_length,
                   bool need_energy = true, int md_step = -1,
                   bool commit_sampling_state = true,
                   std::uint64_t coordinate_generation = 0);
    void Compute_Gradient(VECTOR* local_frc, const VECTOR* global_crd,
                          const int* global_to_local, int owned_atom_numbers,
                          const VECTOR box_length, int need_virial = 0,
                          LTMatrix3* atom_virial = NULL);
    void Accumulate_Energy(float* local_atom_energy, const int* global_to_local,
                           int owned_atom_numbers);

    // 外部查询与输出
    void Step_Print(CONTROLLER* controller);

   private:
    CONTROLLER* controller = NULL;

    // This cache is the last accepted SCF history. Force-evaluation scratch
    // is restored from it before every solve; only a committed evaluation may
    // replace it.
    float* d_accepted_alpha_density = NULL;
    float* d_accepted_beta_density = NULL;
    float* d_accepted_env = NULL;
    float* d_accepted_aux_env = NULL;
    int* d_scf_validation_failure = NULL;
    int* d_nuclear_overlap_pair = NULL;
    int* d_nuclear_geometry_failure = NULL;
    bool accepted_need_initial_guess = true;
    // True only when the accepted density passed the strict repeated
    // ensemble KKT certificate.  Trial force evaluations restore this marker
    // but may not publish a replacement for it.
    bool accepted_density_is_ensemble = false;
    std::uint64_t accepted_coordinate_generation = 0;

    // 轨道基组��称（RI 初始化时用于映射辅助基）
    std::string orbital_basis_name;
    // ECP 选择: AUTO=根据基组自动选择, NONE=禁用, DEF2_ECP/LANL2DZ=显式指定
    QC_ECP_TYPE ecp_type = QC_ECP_TYPE::AUTO;

    // 初始化内部流程
    bool Parsing_Arguments(CONTROLLER* controller, const int atom_numbers,
                           const char*& qc_type_file,
                           std::string& basis_set_name);
    void Initial_Molecule(CONTROLLER* controller, const char* qc_type_file,
                          const std::string& basis_set_name);
    void Initial_Integral_Tasks(CONTROLLER* controller);
    void Memory_Allocate(CONTROLLER* controller);
    void Build_SCF_Workspace();

    // 积分与基组变换内部流程
    void Build_Cart2Sph_Matrix();
    void Cart2Sph_OneE_Integrals();
    void Cart2Sph_Single_Matrix(float* d_cart, float* d_sph);

    // 坐标更新
    bool Initialize_Coordinates_From_MD(const VECTOR* crd,
                                        const VECTOR box_length);
    bool Update_Coordinates_From_MD(const VECTOR* crd,
                                    const VECTOR box_length);
    void Refresh_Coordinate_Derived_State_From_Env();
    bool Validate_Nuclear_Geometry(const VECTOR box_length);

    // 积分
    void Compute_OneE_Integrals();
    void Compute_ECP_Matrix();
    void Compute_Nuclear_Repulsion(const VECTOR box_length);
    void Compute_Analytical_Norms();
    void Build_Shell_Pair_Bounds();
    void Prepare_Integrals();

    // DFT VXC 构建
    void Update_DFT_Grid();
    void Build_DFT_VXC();

    // RI (Density Fitting) 内部流程
    void Initial_Auxiliary_Basis(CONTROLLER* controller);
    void RI_Memory_Allocate();
    void RI_Precompute();
    void Build_Fock_RI();
    bool Factor_RI_Spin_Density(const float* d_density, double density_scale,
                                float* d_factor, int* factor_rank,
                                QC_SCF_Eigensolver_Channel spin_channel);

    // SCF 循环内部流程
    void Build_Fock(int iter, bool force_full_rebuild);
    void Accumulate_SCF_Energy(int iter);
    void Apply_DIIS(int iter);
    bool Diagonalize_And_Build_Density();
    bool Check_Convergence(int iter, int md_step, double energy,
                           double delta_energy, bool physical_iteration);
    bool Start_Ensemble_Probe(double energy, double density_residual);
    int Advance_Ensemble_Search(int iter, int md_step, double energy,
                                double delta_energy);
    double Ensemble_Commutator_RMS();
    bool Build_Overlap_X();
    void Reset_SCF_State();
    bool Build_Initial_Guess();
    void Build_DFT_XC_Gradient();
    void Build_RI_Gradient();
    bool Diag_Guess_And_Build_P();
    void Compute_Spin_Square();
};
