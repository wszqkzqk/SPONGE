#pragma once

#include "../../common.h"
#include "../scf/convergence_policy.hpp"
#include "../scf/ensemble_policy.hpp"
#include "ri_workspace.h"

// 持久 AO 核心矩阵与能量结果缓存
struct QC_SCF_Core_Matrices
{
    float* d_S = NULL;
    float* d_T = NULL;
    float* d_V = NULL;
    float* d_V_ECP = NULL;
    float* d_H_core = NULL;
    double* d_scf_energy = NULL;
    double* d_nuc_energy_dev = NULL;
};

// 单个自旋通道的 Fock、密度与轨道系数工作区
struct QC_SCF_Spin_Channel
{
    std::vector<float> h_F;
    float* d_F = NULL;
    std::vector<float> h_P;
    float* d_P = NULL;
    std::vector<float> h_P_new;
    float* d_P_new = NULL;
    std::vector<float> h_C;
    float* d_C = NULL;
    double* d_F_double = NULL;
    double* d_F_for_grad = NULL;  // 梯度用: 缓存 Build_Fock 后、DIIS 前的 Fock
};

// 重叠正交化、本征分解与双精度临时缓冲
struct QC_SCF_Ortho_Workspace
{
    std::vector<float> h_X;
    double* d_X = NULL;
    std::vector<float> h_W;
    float* d_W = NULL;
    float* d_W_alpha = NULL;  // UHF: 保存 alpha 特征值（beta 会覆盖 d_W）
    float* d_norms = NULL;

    double* d_dwork_nao2_1 = NULL;
    double* d_dwork_nao2_2 = NULL;
    double* d_dwork_nao2_3 = NULL;
    double* d_dwork_nao2_4 = NULL;
    double* d_dW_double = NULL;
    double* d_solver_work_double = NULL;

    int lwork_double = 0;
    double lindep_threshold = 1e-6;
    int nao_eff = 0;
};

// DIIS/EDIIS/ADIIS/MESA 迭代加速相关缓冲与历史状态
//
// MESA 算法参考:
//   S. Lehtola, "OpenOrbitalOptimizer - a reusable open source library for
//   self-consistent field calculations", arXiv:2503.23034 (2025).
//
// EDIIS 参考:
//   K. N. Kudin, G. E. Scuseria, E. Cancès, J. Chem. Phys. 116, 8255 (2002).
//
// ADIIS 参考:
//   X. Hu, W. Yang, J. Chem. Phys. 132, 054109 (2010).
struct QC_SCF_DIIS_Workspace
{
    double* d_diis_err = NULL;
    float *d_diis_w1 = NULL, *d_diis_w2 = NULL, *d_diis_w3 = NULL,
          *d_diis_w4 = NULL;
    // CDIIS: Fock 和误差向量历史
    std::vector<double*> d_diis_f_hist;
    std::vector<double*> d_diis_e_hist;
    std::vector<double*> d_diis_f_hist_b;
    std::vector<double*> d_diis_e_hist_b;

    // EDIIS/ADIIS: 密度矩阵历史和 SCF 能量历史
    std::vector<double*> d_diis_d_hist;
    std::vector<double*> d_diis_d_hist_b;
    std::vector<double> energy_hist;

    int diis_hist_count = 0;
    int diis_hist_head = 0;
    // MESA → CDIIS 切换阈值
    double mesa_to_cdiis_threshold = 1e-1;

    // 最近一轮的 DIIS 误差范数（用于动态 level shift）
    double last_enorm = 1e10;

    double* d_diis_accum = NULL;

    // 批量 Tr 用连续缓冲: 两块 [diis_space × nao²]，存 gather 后的历史向量
    double* d_gather_a = NULL;
    double* d_gather_b = NULL;
    // 批量 Tr 输出: [diis_space × diis_space]
    double* d_dot_out = NULL;
};

// direct SCF 的 pair density、线程私有 Fock 与 ERI 工作池
struct QC_SCF_Direct_Workspace
{
    float* d_pair_density_coul = NULL;
    float* d_pair_density_exx = NULL;
    float* d_pair_density_exx_b = NULL;
    float* d_hr_pool = NULL;
    void* h_cpu_bra_terms = NULL;
    void* h_cpu_ket_terms = NULL;

    int fock_thread_count = 1;
    double* d_F_thread = NULL;
    double* d_F_b_thread = NULL;

    float* d_Ptot = NULL;
    float* d_P_coul = NULL;

    // Incremental Fock: previous densities and accumulated ERI Fock
    float* d_P_coul_prev = NULL;
    float* d_P_exx_prev = NULL;
    float* d_P_exx_b_prev = NULL;
    double* d_F_eri_accum = NULL;  // CPU: double-precision accumulator
    double* d_F_eri_b_accum = NULL;
    float* d_F_eri_accum_f = NULL;  // GPU: float accumulator
    float* d_F_eri_b_accum_f = NULL;
};

enum QC_SCF_Fock_Build_Mode
{
    QC_SCF_FOCK_BUILD_FULL = 0,
    QC_SCF_FOCK_BUILD_INCREMENTAL = 1,
    QC_SCF_FOCK_BUILD_RI = 2,
};

enum QC_SCF_Ensemble_Phase
{
    QC_SCF_ENSEMBLE_INACTIVE = 0,
    QC_SCF_ENSEMBLE_PROBE_UPPER = 1,
    QC_SCF_ENSEMBLE_EVALUATE_ROOT = 2,
    QC_SCF_ENSEMBLE_VERIFY_KKT_REPEAT = 3,
    QC_SCF_ENSEMBLE_VERIFY_LINE_BEST = 4,
    QC_SCF_ENSEMBLE_VERIFY_COMMITTED = 5,
    QC_SCF_ENSEMBLE_VERIFY_LINE_STATIONARITY = 6,
    QC_SCF_ENSEMBLE_RECOVER_COMMITTED = 7,
    QC_SCF_ENSEMBLE_PROBE_SPECTRAL_ORBITAL = 8,
};

struct QC_SCF_Ensemble_Atom
{
    float* d_alpha = NULL;
    float* d_beta = NULL;
    double weight = 0.0;
};

// Dedicated storage for the zero-temperature ODA/Frank-Wolfe fallback.  The
// line is D(x)=origin+x*direction for both spin channels; every point is
// evaluated with an unshifted, non-extrapolated Fock build.
struct QC_SCF_Ensemble_Workspace
{
    QC_SCF_Ensemble_Phase phase = QC_SCF_ENSEMBLE_INACTIVE;
    QC_SCF_Ensemble_Bracket bracket;
    double line_origin_derivative = 0.0;
    double line_origin_energy = 0.0;
    double direction_density_rms = 0.0;
    double line_derivative_tolerance = 1.0e-6;
    double line_energy_guard_multiplier = 2.0;
    double current_fraction = 0.0;
    double maximum_fraction = 1.0;
    // Energy of the current committed density, evaluated only after an
    // active-set reconstruction receives a fresh raw F[P] build.  A trial
    // energy from the preceding line search must never update this value.
    double committed_energy = 0.0;
    double verification_global_fw_gap = 0.0;
    double verification_commutator = 0.0;
    double verification_active_fw_gap = 0.0;
    double active_fw_gap = 0.0;
    double spectral_orbital_origin_commutator = 0.0;
    // Fixed-origin-F contribution of the nominal-particle and float-
    // representation repair, separate from the spectrum-preserving orbital
    // rotation.  Spectral raw-energy guards remove only this normal first-
    // order term; committed energies and final KKT certificates stay raw.
    double spectral_orbital_repair_linear_component = 0.0;
    int probe_evaluations = 0;
    bool confirmed = false;
    bool interior_minimum_confirmed = false;
    bool corrective_line_used_quasi_newton = false;
    bool corrective_line_used_pairwise = false;
    bool corrective_force_pairwise = false;
    std::vector<QC_SCF_Ensemble_Atom> active_atoms;
    std::vector<double> line_weight_direction;
    // Reduced-coordinate inverse Hessian and the preceding raw active-hull
    // observation used by the fully-corrective quasi-Newton solve.  Coordinate
    // i represents weight i+1, with weight 0 fixed by normalization.
    std::vector<double> corrective_inverse_hessian;
    std::vector<double> corrective_reference_weights;
    std::vector<double> corrective_reference_linear_values;

    float* d_origin_alpha = NULL;
    float* d_direction_alpha = NULL;
    float* d_origin_beta = NULL;
    float* d_direction_beta = NULL;
    // The physical, symmetric Fock operator at the origin of a spectral
    // orbital line.  Later raw Fock builds overwrite each spin channel's
    // regular d_F_double, so the line's rounding objective needs dedicated
    // persistent storage.
    double* d_spectral_origin_fock_alpha = NULL;
    double* d_spectral_origin_fock_beta = NULL;
    double* d_commutator = NULL;
    double* d_accum = NULL;
};

// SCF 配置、收敛状态与能量累计缓冲
struct QC_SCF_Runtime_State
{
    bool unrestricted = false;
    int n_alpha = 0;
    int n_beta = 0;
    float occ_factor = 2.0f;
    int max_scf_iter = 300;
    bool use_diis = true;
    int diis_start_iter = 8;
    int diis_space = 6;
    double diis_reg = 1e-10;
    double energy_tol = 1e-6;
    double density_tol = 1e-6;
    int convergence_streak = 0;
    QC_SCF_Convergence_State convergence;
    double configured_level_shift = 0.25;
    double level_shift = 0.25;
    bool print_iter = false;
    QC_SCF_Fock_Build_Mode fock_build_mode = QC_SCF_FOCK_BUILD_FULL;

    double spin_square = 0.0;
    double spin_square_exact = 0.0;

    double* d_e = NULL;
    double* d_e_b = NULL;
    double* d_pvxc = NULL;
    double* d_prev_energy = NULL;
    double* d_delta_e = NULL;
    double* d_density_residual = NULL;
    int* d_converged = NULL;
};

// SCF 工作空间（初始化分配，循环中复用）
struct QC_SCF_WORKSPACE
{
    QC_SCF_Core_Matrices core;
    QC_SCF_Spin_Channel alpha;
    QC_SCF_Spin_Channel beta;
    QC_SCF_Ortho_Workspace ortho;
    QC_SCF_DIIS_Workspace diis;
    QC_SCF_Direct_Workspace direct;
    QC_SCF_Runtime_State runtime;
    QC_SCF_Ensemble_Workspace ensemble;
    QC_RI_WORKSPACE ri;
};
