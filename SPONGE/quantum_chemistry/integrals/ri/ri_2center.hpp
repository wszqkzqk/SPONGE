#pragma once

// 二中心 Coulomb 积分 (P|Q) = ∫∫ P(r1) 1/r12 Q(r2) dr1 dr2
// 使用 McMurchie-Davidson 方案。

#include "ri_workspace.hpp"

// 二中心 Coulomb 积分内核
// 每个 worker 处理若干辅助 shell 对 (P_sh, Q_sh)，并独占一份动态
// E/Boys/R workspace。
static __global__ void QC_RI_2Center_Kernel(
    const int n_tasks, const QC_ONE_E_TASK* tasks, const VECTOR* centers,
    const int* l_list, const float* exps, const float* coeffs,
    const int* shell_offsets, const int* shell_sizes, const int* ao_offsets,
    int naux, QC_RI_DEVICE_WORKSPACE workspace, double* out_metric)
{
    SIMPLE_DEVICE_FOR(worker_id, workspace.n_workers)
    {
        float* worker_e = workspace.e + (size_t)worker_id * workspace.e_stride;
        const size_t e_tensor_size = QC_RI_E_Tensor_Size(
            workspace.e_dim_a, workspace.e_dim_b, workspace.e_dim_n);
        float* E_Px = worker_e;
        float* E_Py = E_Px + e_tensor_size;
        float* E_Pz = E_Py + e_tensor_size;
        float* E_Qx = E_Pz + e_tensor_size;
        float* E_Qy = E_Qx + e_tensor_size;
        float* E_Qz = E_Qy + e_tensor_size;
        float* R_vals = workspace.r + (size_t)worker_id * workspace.r_stride;
        double* F_vals =
            workspace.boys + (size_t)worker_id * workspace.boys_stride;

        for (size_t task_id = (size_t)worker_id; task_id < (size_t)n_tasks;
             task_id += (size_t)workspace.n_workers)
        {
            const QC_ONE_E_TASK sh = tasks[task_id];
            const int P_sh = sh.x;
            const int Q_sh = sh.y;

            const int lP = l_list[P_sh], lQ = l_list[Q_sh];
            const int nP = (int)QC_RI_Cartesian_Count(lP);
            const int nQ = (int)QC_RI_Cartesian_Count(lQ);
            const int offP = ao_offsets[P_sh];
            const int offQ = ao_offsets[Q_sh];

            const VECTOR A = centers[P_sh];
            const VECTOR B = centers[Q_sh];
            const float Ax = A.x, Ay = A.y, Az = A.z;
            const float Bx = B.x, By = B.y, Bz = B.z;
            const float dist_sq = (Ax - Bx) * (Ax - Bx) +
                                  (Ay - By) * (Ay - By) + (Az - Bz) * (Az - Bz);

            for (int idxP = 0; idxP < nP; idxP++)
            {
                for (int idxQ = 0; idxQ < nQ; idxQ++)
                {
                    int lxP, lyP, lzP, lxQ, lyQ, lzQ;
                    QC_Get_Lxyz_Device(lP, idxP, lxP, lyP, lzP);
                    QC_Get_Lxyz_Device(lQ, idxQ, lxQ, lyQ, lzQ);

                    double total = 0.0;

                    for (int pi = 0; pi < shell_sizes[P_sh]; pi++)
                    {
                        const float eP = exps[shell_offsets[P_sh] + pi];
                        const float cP = coeffs[shell_offsets[P_sh] + pi];

                        for (int pj = 0; pj < shell_sizes[Q_sh]; pj++)
                        {
                            const float eQ = exps[shell_offsets[Q_sh] + pj];
                            const float cQ = coeffs[shell_offsets[Q_sh] + pj];

                            QC_RI_Compute_MD_Coeffs(E_Px, workspace.e_dim_a,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n, lxP, 0,
                                                    0.0f, 0.0f, 0.5f / eP);
                            QC_RI_Compute_MD_Coeffs(E_Py, workspace.e_dim_a,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n, lyP, 0,
                                                    0.0f, 0.0f, 0.5f / eP);
                            QC_RI_Compute_MD_Coeffs(E_Pz, workspace.e_dim_a,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n, lzP, 0,
                                                    0.0f, 0.0f, 0.5f / eP);
                            QC_RI_Compute_MD_Coeffs(E_Qx, workspace.e_dim_a,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n, lxQ, 0,
                                                    0.0f, 0.0f, 0.5f / eQ);
                            QC_RI_Compute_MD_Coeffs(E_Qy, workspace.e_dim_a,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n, lyQ, 0,
                                                    0.0f, 0.0f, 0.5f / eQ);
                            QC_RI_Compute_MD_Coeffs(E_Qz, workspace.e_dim_a,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n, lzQ, 0,
                                                    0.0f, 0.0f, 0.5f / eQ);

                            const float alpha_pq = eP * eQ / (eP + eQ);
                            const float T_val = alpha_pq * dist_sq;
                            const int L_tot = lP + lQ;
                            const int r_n_stride = L_tot + 1;
                            QC_RI_Compute_Boys_Double(F_vals, T_val, L_tot);
                            float AB[3] = {Ax - Bx, Ay - By, Az - Bz};
                            QC_RI_Compute_R_Tensor(R_vals, F_vals, alpha_pq, AB,
                                                   L_tot);

                            const double prefactor =
                                (double)cP * (double)cQ *
                                (2.0 * CONSTANT_Pi * CONSTANT_Pi *
                                 sqrt(CONSTANT_Pi)) /
                                ((double)eP * (double)eQ *
                                 sqrt((double)(eP + eQ)));

                            double v_sum = 0.0;
                            for (int t = 0; t <= lxP; t++)
                            {
                                const double ePx = (double)E_Px[QC_RI_E_Index(
                                    lxP, 0, t, workspace.e_dim_b,
                                    workspace.e_dim_n)];
                                if (ePx == 0.0) continue;
                                for (int u = 0; u <= lyP; u++)
                                {
                                    const double ePy =
                                        (double)E_Py[QC_RI_E_Index(
                                            lyP, 0, u, workspace.e_dim_b,
                                            workspace.e_dim_n)];
                                    if (ePy == 0.0) continue;
                                    for (int v = 0; v <= lzP; v++)
                                    {
                                        const double ePz =
                                            (double)E_Pz[QC_RI_E_Index(
                                                lzP, 0, v, workspace.e_dim_b,
                                                workspace.e_dim_n)];
                                        if (ePz == 0.0) continue;
                                        for (int tt = 0; tt <= lxQ; tt++)
                                        {
                                            const double eQx =
                                                (double)E_Qx[QC_RI_E_Index(
                                                    lxQ, 0, tt,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n)];
                                            if (eQx == 0.0) continue;
                                            for (int uu = 0; uu <= lyQ; uu++)
                                            {
                                                const double eQy =
                                                    (double)E_Qy[QC_RI_E_Index(
                                                        lyQ, 0, uu,
                                                        workspace.e_dim_b,
                                                        workspace.e_dim_n)];
                                                if (eQy == 0.0) continue;
                                                for (int vv = 0; vv <= lzQ;
                                                     vv++)
                                                {
                                                    const double eQz = (double)
                                                        E_Qz[QC_RI_E_Index(
                                                            lzQ, 0, vv,
                                                            workspace.e_dim_b,
                                                            workspace.e_dim_n)];
                                                    if (eQz == 0.0) continue;
                                                    const double sign =
                                                        ((tt + uu + vv) & 1)
                                                            ? -1.0
                                                            : 1.0;
                                                    v_sum +=
                                                        ePx * ePy * ePz * eQx *
                                                        eQy * eQz * sign *
                                                        (double)R_vals
                                                            [QC_RI_R_Index(
                                                                t + tt, u + uu,
                                                                v + vv, 0,
                                                                r_n_stride)];
                                                }
                                            }
                                        }
                                    }
                                }
                            }

                            total += prefactor * v_sum;
                        }
                    }

                    const size_t P_idx = (size_t)(offP + idxP);
                    const size_t Q_idx = (size_t)(offQ + idxQ);
                    out_metric[P_idx * (size_t)naux + Q_idx] = total;
                    if (P_sh != Q_sh)
                        out_metric[Q_idx * (size_t)naux + P_idx] = total;
                }
            }
        }
    }
}

static inline void QC_Launch_RI_2Center_Kernel(
    int threads, int n_tasks, const QC_ONE_E_TASK* tasks, const VECTOR* centers,
    const int* l_list, const float* exps, const float* coeffs,
    const int* shell_offsets, const int* shell_sizes, const int* ao_offsets,
    int naux, const QC_RI_INTEGRAL_WORKSPACE& workspace, double* out_metric)
{
    if (n_tasks == 0) return;
    const int grid = (workspace.device.n_workers + threads - 1) / threads;
    Launch_Device_Kernel(QC_RI_2Center_Kernel, grid, threads, 0, 0, n_tasks,
                         tasks, centers, l_list, exps, coeffs, shell_offsets,
                         shell_sizes, ao_offsets, naux, workspace.device,
                         out_metric);
}
