#pragma once

#include <vector>

#include "../structure/integral_tasks.h"
#include "../structure/molecule.h"

// 球谐→笛卡尔密度变换 (含归一化):
//   M_cart[nc×nc] = U (N · M_sph · N) U^T
// 在 CPU 上执行 (矩阵通常 <100×100, 远小于 PCIe 开销)
// h_norms: [ns], h_C/U: [nc×ns] row-major, h_M_sph: [ns×ns]
void QC_Sph2Cart_Density_Host(int ns, int nc, const std::vector<float>& h_norms,
                              const std::vector<float>& h_C,
                              const std::vector<float>& h_M_sph,
                              std::vector<float>& h_M_cart);

// 计算 ECP 矩阵 V_ECP (Cartesian 基)
// 输出到 d_V_ECP[nao_cart × nao_cart]，调用前需 memset 为 0
enum QC_ECP_EVALUATION_FAILURE_KIND
{
    QC_ECP_EVALUATION_OK = 0,
    QC_ECP_PRIMITIVE_QUADRATURE_NOT_CONVERGED = 1,
    QC_ECP_NONFINITE_VALUE_OR_ESTIMATE = 2,
    QC_ECP_MATRIX_ESTIMATE_EXCEEDS_BUDGET = 3,
    QC_ECP_MATRIX_STORAGE_EXCEEDS_BUDGET = 4,
    QC_ECP_GRADIENT_ESTIMATE_EXCEEDS_BUDGET = 5,
    QC_ECP_RESOURCE_ALLOCATION_FAILED = 6,
};

inline const char* QC_ECP_Evaluation_Failure_Kind_Name(
    QC_ECP_EVALUATION_FAILURE_KIND kind)
{
    switch (kind)
    {
        case QC_ECP_EVALUATION_OK:
            return "no failure";
        case QC_ECP_PRIMITIVE_QUADRATURE_NOT_CONVERGED:
            return "primitive quadrature did not converge";
        case QC_ECP_NONFINITE_VALUE_OR_ESTIMATE:
            return "non-finite input, accumulated value, or error estimate";
        case QC_ECP_MATRIX_ESTIMATE_EXCEEDS_BUDGET:
            return "contracted ECP matrix integration estimate exceeds budget";
        case QC_ECP_MATRIX_STORAGE_EXCEEDS_BUDGET:
            return "float ECP matrix storage error exceeds its representation budget";
        case QC_ECP_GRADIENT_ESTIMATE_EXCEEDS_BUDGET:
            return "final signed ECP gradient estimate exceeds budget";
        case QC_ECP_RESOURCE_ALLOCATION_FAILED:
            return "ECP evaluation scratch allocation failed";
    }
    return "unknown ECP evaluation failure";
}

struct QC_ECP_EVALUATION_FAILURE
{
    QC_ECP_EVALUATION_FAILURE_KIND kind = QC_ECP_EVALUATION_OK;
    int task_id = -1;
    int shell_i = -1;
    int shell_j = -1;
    int atom = -1;
    int observable_atom = -1;
    int direction = -1;
    int channel = -1;
    int term = -1;
    int primitive_i = -1;
    int primitive_j = -1;
    int cartesian_i = -1;
    int cartesian_j = -1;
    int channel_l = -1;
    int n_k = -1;
    double value = 0.0;
    double estimated_error = 0.0;
};

bool QC_Compute_V_ECP(const QC_MOLECULE& mol,
                      const QC_INTEGRAL_TASKS& task_ctx, float* d_V_ECP,
                      QC_ECP_EVALUATION_FAILURE* failure);

// 计算 ECP 梯度: Tr[P · dV_ECP/dR] 对各原子的贡献
// 累加到 d_grad[natm * 3] (double 精度, Hartree/Bohr)
// d_P_cart_eff: Cartesian 基密度矩阵 (含归一化, nao_cart × nao_cart)
//   球谐基时为 U(NPN)U^T, 笛卡尔基时为 P.*norms*norms'
bool QC_Compute_ECP_Gradient(const QC_MOLECULE& mol,
                             const QC_INTEGRAL_TASKS& task_ctx,
                             const int* d_shell_atom,
                             const float* d_P_cart_eff, double* d_grad,
                             QC_ECP_EVALUATION_FAILURE* failure);
