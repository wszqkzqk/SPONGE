#pragma once

// 二中心 Coulomb 积分导数 d(P|Q)/dR 的内核
// 使用 McMurchie-Davidson 方案，与 ri_2center.hpp 对应
// 导数公式: d(P|Q)/dA_{P,x} = 2*alpha_P * ((P+1_x)|Q) - l_{P,x} * ((P-1_x)|Q)
// 平移不变性: d(P|Q)/dA_Q = -d(P|Q)/dA_P

#include "ri_workspace.hpp"

// 二中心积分导数内核
// 有界 worker 池按 stride 遍历辅助 shell P_sh，内部循环 Q_sh <= P_sh。
// 每个 worker 独占动态 Cartesian/E/R/Boys workspace。
static __global__ void QC_RI_2Center_Grad_Kernel(
    const int naux_bas, const VECTOR* aux_centers, const int* aux_l_list,
    const float* aux_exps, const float* aux_coeffs,
    const int* aux_shell_offsets, const int* aux_shell_sizes,
    const int* aux_ao_offsets_cart, const int* aux_ao_offsets_sph,
    const float* aux_norms, const float* U_aux, int naux_cart, int naux_sph,
    const double* D2_eff, const int* shell_atom_aux,
    QC_RI_DEVICE_WORKSPACE workspace, double* grad)
{
    SIMPLE_DEVICE_FOR(worker_id, workspace.n_workers)
    {
        const int e_dim_a = workspace.e_dim_a;
        const int e_dim_b = workspace.e_dim_b;
        const int e_dim_n = workspace.e_dim_n;
        double* d_cart =
            workspace.cart + (size_t)worker_id * workspace.cart_stride;
        float* worker_e = workspace.e + (size_t)worker_id * workspace.e_stride;
        const size_t e_tensor_size =
            QC_RI_E_Tensor_Size(e_dim_a, e_dim_b, e_dim_n);
        float* E_Px = worker_e;
        float* E_Py = E_Px + e_tensor_size;
        float* E_Pz = E_Py + e_tensor_size;
        float* E_Qx = E_Pz + e_tensor_size;
        float* E_Qy = E_Qx + e_tensor_size;
        float* E_Qz = E_Qy + e_tensor_size;
        float* R_vals = workspace.r + (size_t)worker_id * workspace.r_stride;
        double* F_vals =
            workspace.boys + (size_t)worker_id * workspace.boys_stride;

        for (size_t P_sh_index = (size_t)worker_id;
             P_sh_index < (size_t)naux_bas;
             P_sh_index += (size_t)workspace.n_workers)
        {
            const int P_sh = (int)P_sh_index;
            for (int Q_sh = 0; Q_sh <= P_sh; Q_sh++)
            {
                const int lP = aux_l_list[P_sh], lQ = aux_l_list[Q_sh];
                const int nP_cart = (int)QC_RI_Cartesian_Count(lP);
                const int nQ_cart = (int)QC_RI_Cartesian_Count(lQ);
                const int offP_cart = aux_ao_offsets_cart[P_sh];
                const int offQ_cart = aux_ao_offsets_cart[Q_sh];
                const int offP_sph = aux_ao_offsets_sph[P_sh];
                const int offQ_sph = aux_ao_offsets_sph[Q_sh];
                const int nP_sph = 2 * lP + 1;
                const int nQ_sph = 2 * lQ + 1;

                const VECTOR A = aux_centers[P_sh];
                const VECTOR B = aux_centers[Q_sh];
                const float Ax = A.x, Ay = A.y, Az = A.z;
                const float Bx = B.x, By = B.y, Bz = B.z;
                const float dist_sq = (Ax - Bx) * (Ax - Bx) +
                                      (Ay - By) * (Ay - By) +
                                      (Az - Bz) * (Az - Bz);

                const int atom_P = shell_atom_aux[P_sh];
                const int atom_Q = shell_atom_aux[Q_sh];

                // 清零 workspace 用于本 shell pair
                const size_t d_cart_size =
                    (size_t)nP_cart * (size_t)nQ_cart * 3;
                for (size_t k = 0; k < d_cart_size; k++) d_cart[k] = 0.0;

                for (int idxP = 0; idxP < nP_cart; idxP++)
                {
                    int lxP, lyP, lzP;
                    QC_Get_Lxyz_Device(lP, idxP, lxP, lyP, lzP);

                    for (int idxQ = 0; idxQ < nQ_cart; idxQ++)
                    {
                        int lxQ, lyQ, lzQ;
                        QC_Get_Lxyz_Device(lQ, idxQ, lxQ, lyQ, lzQ);

                        double dA[3] = {0.0, 0.0, 0.0};

                        for (int pi = 0; pi < aux_shell_sizes[P_sh]; pi++)
                        {
                            const float eP =
                                aux_exps[aux_shell_offsets[P_sh] + pi];
                            const float cP =
                                aux_coeffs[aux_shell_offsets[P_sh] + pi];

                            for (int pj = 0; pj < aux_shell_sizes[Q_sh]; pj++)
                            {
                                const float eQ =
                                    aux_exps[aux_shell_offsets[Q_sh] + pj];
                                const float cQ =
                                    aux_coeffs[aux_shell_offsets[Q_sh] + pj];

                                QC_RI_Compute_MD_Coeffs(E_Px, e_dim_a, e_dim_b,
                                                        e_dim_n, lxP + 1, 0,
                                                        0.0f, 0.0f, 0.5f / eP);
                                QC_RI_Compute_MD_Coeffs(E_Py, e_dim_a, e_dim_b,
                                                        e_dim_n, lyP + 1, 0,
                                                        0.0f, 0.0f, 0.5f / eP);
                                QC_RI_Compute_MD_Coeffs(E_Pz, e_dim_a, e_dim_b,
                                                        e_dim_n, lzP + 1, 0,
                                                        0.0f, 0.0f, 0.5f / eP);

                                QC_RI_Compute_MD_Coeffs(E_Qx, e_dim_a, e_dim_b,
                                                        e_dim_n, lxQ, 0, 0.0f,
                                                        0.0f, 0.5f / eQ);
                                QC_RI_Compute_MD_Coeffs(E_Qy, e_dim_a, e_dim_b,
                                                        e_dim_n, lyQ, 0, 0.0f,
                                                        0.0f, 0.5f / eQ);
                                QC_RI_Compute_MD_Coeffs(E_Qz, e_dim_a, e_dim_b,
                                                        e_dim_n, lzQ, 0, 0.0f,
                                                        0.0f, 0.5f / eQ);

                                const float alpha_pq = eP * eQ / (eP + eQ);
                                const float T_val = alpha_pq * dist_sq;
                                const int L_tot = lP + lQ;
                                const int R_order = L_tot + 1;
                                const int R_n_stride = R_order + 1;

                                QC_RI_Compute_Boys_Double(F_vals, T_val,
                                                          R_order);
                                float AB[3] = {Ax - Bx, Ay - By, Az - Bz};
                                QC_RI_Compute_R_Tensor(R_vals, F_vals, alpha_pq,
                                                       AB, R_order);

                                const double prefactor =
                                    (double)cP * (double)cQ *
                                    (2.0 * CONSTANT_Pi * CONSTANT_Pi *
                                     sqrt(CONSTANT_Pi)) /
                                    ((double)eP * (double)eQ *
                                     sqrt((double)(eP + eQ)));

                                auto contract_2c = [&](int axP, int ayP,
                                                       int azP) -> double
                                {
                                    if (axP < 0 || ayP < 0 || azP < 0)
                                        return 0.0;
                                    double v_sum = 0.0;
                                    for (int t = 0; t <= axP; t++)
                                    {
                                        double ePx = (double)E_Px[QC_RI_E_Index(
                                            axP, 0, t, e_dim_b, e_dim_n)];
                                        if (ePx == 0.0) continue;
                                        for (int u = 0; u <= ayP; u++)
                                        {
                                            double ePy =
                                                (double)E_Py[QC_RI_E_Index(
                                                    ayP, 0, u, e_dim_b,
                                                    e_dim_n)];
                                            if (ePy == 0.0) continue;
                                            for (int v = 0; v <= azP; v++)
                                            {
                                                double ePz =
                                                    (double)E_Pz[QC_RI_E_Index(
                                                        azP, 0, v, e_dim_b,
                                                        e_dim_n)];
                                                if (ePz == 0.0) continue;
                                                for (int tt = 0; tt <= lxQ;
                                                     tt++)
                                                {
                                                    double eQx = (double)
                                                        E_Qx[QC_RI_E_Index(
                                                            lxQ, 0, tt, e_dim_b,
                                                            e_dim_n)];
                                                    if (eQx == 0.0) continue;
                                                    for (int uu = 0; uu <= lyQ;
                                                         uu++)
                                                    {
                                                        double eQy = (double)
                                                            E_Qy[QC_RI_E_Index(
                                                                lyQ, 0, uu,
                                                                e_dim_b,
                                                                e_dim_n)];
                                                        if (eQy == 0.0)
                                                            continue;
                                                        for (int vv = 0;
                                                             vv <= lzQ; vv++)
                                                        {
                                                            double eQz = (double)
                                                                E_Qz[QC_RI_E_Index(
                                                                    lzQ, 0, vv,
                                                                    e_dim_b,
                                                                    e_dim_n)];
                                                            if (eQz == 0.0)
                                                                continue;
                                                            double sign =
                                                                ((tt + uu +
                                                                  vv) &
                                                                 1)
                                                                    ? -1.0
                                                                    : 1.0;
                                                            v_sum +=
                                                                ePx * ePy *
                                                                ePz * eQx *
                                                                eQy * eQz *
                                                                sign *
                                                                (double)R_vals
                                                                    [QC_RI_R_Index(
                                                                        t + tt,
                                                                        u + uu,
                                                                        v + vv,
                                                                        0,
                                                                        R_n_stride)];
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    return v_sum;
                                };

                                double dx = 2.0 * (double)eP *
                                            contract_2c(lxP + 1, lyP, lzP);
                                if (lxP > 0)
                                    dx -= (double)lxP *
                                          contract_2c(lxP - 1, lyP, lzP);
                                double dy = 2.0 * (double)eP *
                                            contract_2c(lxP, lyP + 1, lzP);
                                if (lyP > 0)
                                    dy -= (double)lyP *
                                          contract_2c(lxP, lyP - 1, lzP);
                                double dz = 2.0 * (double)eP *
                                            contract_2c(lxP, lyP, lzP + 1);
                                if (lzP > 0)
                                    dz -= (double)lzP *
                                          contract_2c(lxP, lyP, lzP - 1);

                                dA[0] += prefactor * dx;
                                dA[1] += prefactor * dy;
                                dA[2] += prefactor * dz;
                            }
                        }

                        const size_t cart_idx =
                            ((size_t)idxP * (size_t)nQ_cart + (size_t)idxQ) * 3;
                        d_cart[cart_idx + 0] = dA[0];
                        d_cart[cart_idx + 1] = dA[1];
                        d_cart[cart_idx + 2] = dA[2];
                    }
                }

                // Cart2Sph + 归一化 + D2_eff 收缩 → atomicAdd 到 grad
                for (int ps = 0; ps < nP_sph; ps++)
                {
                    const int P_sph = offP_sph + ps;
                    const double normP = (double)aux_norms[P_sph];

                    for (int qs = 0; qs < nQ_sph; qs++)
                    {
                        const int Q_sph = offQ_sph + qs;
                        const double normQ = (double)aux_norms[Q_sph];

                        double d_sph[3] = {0.0, 0.0, 0.0};
                        for (int pc = 0; pc < nP_cart; pc++)
                        {
                            double u_p = (double)
                                U_aux[(offP_cart + pc) * naux_sph + P_sph];
                            if (u_p == 0.0) continue;
                            for (int qc = 0; qc < nQ_cart; qc++)
                            {
                                double u_q = (double)
                                    U_aux[(offQ_cart + qc) * naux_sph + Q_sph];
                                if (u_q == 0.0) continue;
                                double w = u_p * u_q;
                                const size_t cidx =
                                    ((size_t)pc * (size_t)nQ_cart +
                                     (size_t)qc) *
                                    3;
                                d_sph[0] += w * d_cart[cidx + 0];
                                d_sph[1] += w * d_cart[cidx + 1];
                                d_sph[2] += w * d_cart[cidx + 2];
                            }
                        }

                        double norm_pq = normP * normQ;
                        d_sph[0] *= norm_pq;
                        d_sph[1] *= norm_pq;
                        d_sph[2] *= norm_pq;

                        double dens = D2_eff[(size_t)P_sph * naux_sph + Q_sph];
                        if (P_sh != Q_sh)
                            dens += D2_eff[(size_t)Q_sph * naux_sph + P_sph];

                        for (int d = 0; d < 3; d++)
                        {
                            double contrib = dens * d_sph[d];
                            atomicAdd(&grad[atom_P * 3 + d], contrib);
                            atomicAdd(&grad[atom_Q * 3 + d], -contrib);
                        }
                    }
                }
            }
        }
    }
}
