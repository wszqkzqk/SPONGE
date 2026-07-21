#pragma once

// 三中心 Coulomb 积分 (P|μν) = ∫∫ P(r1) 1/r12 μ(r2)ν(r2) dr1 dr2
// 使用 McMurchie-Davidson 方案。primitive 循环在外，Cartesian 分量
// 循环在内，E 系数和 Boys/R 张量只计算一次后查表收缩。

#include "ri_workspace.hpp"

struct QC_RI_3C_TASK
{
    int P_sh, mu_sh, nu_sh;
};

static __global__ void QC_RI_3Center_Kernel(
    const int n_tasks, const QC_RI_3C_TASK* tasks, const VECTOR* aux_centers,
    const int* aux_l_list, const float* aux_exps, const float* aux_coeffs,
    const int* aux_shell_offsets, const int* aux_shell_sizes,
    const int* aux_ao_offsets, const VECTOR* orb_centers, const int* orb_l_list,
    const float* orb_exps, const float* orb_coeffs,
    const int* orb_shell_offsets, const int* orb_shell_sizes,
    const int* orb_ao_offsets, int naux, int out_mu_dim, int out_nu_dim,
    int mu_offset_base, int nu_offset_base, bool fill_symmetric,
    QC_RI_DEVICE_WORKSPACE workspace, double* out_eri3c)
{
    SIMPLE_DEVICE_FOR(worker_id, workspace.n_workers)
    {
        float* worker_e = workspace.e + (size_t)worker_id * workspace.e_stride;
        const size_t e_tensor_size = QC_RI_E_Tensor_Size(
            workspace.e_dim_a, workspace.e_dim_b, workspace.e_dim_n);
        float* E_Px = worker_e;
        float* E_Py = E_Px + e_tensor_size;
        float* E_Pz = E_Py + e_tensor_size;
        float* E_Kx = E_Pz + e_tensor_size;
        float* E_Ky = E_Kx + e_tensor_size;
        float* E_Kz = E_Ky + e_tensor_size;
        float* R_vals = workspace.r + (size_t)worker_id * workspace.r_stride;
        double* F_vals =
            workspace.boys + (size_t)worker_id * workspace.boys_stride;
        double* buf =
            workspace.cart + (size_t)worker_id * workspace.cart_stride;

        for (size_t task_id = (size_t)worker_id; task_id < (size_t)n_tasks;
             task_id += (size_t)workspace.n_workers)
        {
            const QC_RI_3C_TASK task = tasks[task_id];
            const int P_sh = task.P_sh;
            const int mu_sh = task.mu_sh;
            const int nu_sh = task.nu_sh;

            const int lP = aux_l_list[P_sh];
            const int lmu = orb_l_list[mu_sh];
            const int lnu = orb_l_list[nu_sh];

            const int nP = (int)QC_RI_Cartesian_Count(lP);
            const int nmu = (int)QC_RI_Cartesian_Count(lmu);
            const int nnu = (int)QC_RI_Cartesian_Count(lnu);

            const int offP = aux_ao_offsets[P_sh];
            const int offmu = orb_ao_offsets[mu_sh];
            const int offnu = orb_ao_offsets[nu_sh];

            const VECTOR A = aux_centers[P_sh];
            const VECTOR B = orb_centers[mu_sh];
            const VECTOR C = orb_centers[nu_sh];

            const size_t n_cart = (size_t)nP * (size_t)nmu * (size_t)nnu;
            for (size_t i = 0; i < n_cart; i++) buf[i] = 0.0;

            for (int pP = 0; pP < aux_shell_sizes[P_sh]; pP++)
            {
                const float eP = aux_exps[aux_shell_offsets[P_sh] + pP];
                const float cP = aux_coeffs[aux_shell_offsets[P_sh] + pP];

                QC_RI_Compute_MD_Coeffs(E_Px, workspace.e_dim_a,
                                        workspace.e_dim_b, workspace.e_dim_n,
                                        lP, 0, 0.0f, 0.0f, 0.5f / eP);
                QC_RI_Compute_MD_Coeffs(E_Py, workspace.e_dim_a,
                                        workspace.e_dim_b, workspace.e_dim_n,
                                        lP, 0, 0.0f, 0.0f, 0.5f / eP);
                QC_RI_Compute_MD_Coeffs(E_Pz, workspace.e_dim_a,
                                        workspace.e_dim_b, workspace.e_dim_n,
                                        lP, 0, 0.0f, 0.0f, 0.5f / eP);

                for (int p_mu = 0; p_mu < orb_shell_sizes[mu_sh]; p_mu++)
                {
                    const float e_mu =
                        orb_exps[orb_shell_offsets[mu_sh] + p_mu];
                    const float c_mu =
                        orb_coeffs[orb_shell_offsets[mu_sh] + p_mu];

                    for (int p_nu = 0; p_nu < orb_shell_sizes[nu_sh]; p_nu++)
                    {
                        const float e_nu =
                            orb_exps[orb_shell_offsets[nu_sh] + p_nu];
                        const float c_nu =
                            orb_coeffs[orb_shell_offsets[nu_sh] + p_nu];

                        const float g_ket = e_mu + e_nu;
                        const float BC2 = (B.x - C.x) * (B.x - C.x) +
                                          (B.y - C.y) * (B.y - C.y) +
                                          (B.z - C.z) * (B.z - C.z);
                        const float K_ket = expf(-e_mu * e_nu / g_ket * BC2);

                        const float Qx = (e_mu * B.x + e_nu * C.x) / g_ket;
                        const float Qy = (e_mu * B.y + e_nu * C.y) / g_ket;
                        const float Qz = (e_mu * B.z + e_nu * C.z) / g_ket;

                        QC_RI_Compute_MD_Coeffs(
                            E_Kx, workspace.e_dim_a, workspace.e_dim_b,
                            workspace.e_dim_n, lmu, lnu, Qx - B.x, Qx - C.x,
                            0.5f / g_ket);
                        QC_RI_Compute_MD_Coeffs(
                            E_Ky, workspace.e_dim_a, workspace.e_dim_b,
                            workspace.e_dim_n, lmu, lnu, Qy - B.y, Qy - C.y,
                            0.5f / g_ket);
                        QC_RI_Compute_MD_Coeffs(
                            E_Kz, workspace.e_dim_a, workspace.e_dim_b,
                            workspace.e_dim_n, lmu, lnu, Qz - B.z, Qz - C.z,
                            0.5f / g_ket);

                        const float alpha_pq = eP * g_ket / (eP + g_ket);
                        const float AQ2 = (A.x - Qx) * (A.x - Qx) +
                                          (A.y - Qy) * (A.y - Qy) +
                                          (A.z - Qz) * (A.z - Qz);
                        const float T_val = alpha_pq * AQ2;
                        const int L_tot = lP + lmu + lnu;
                        const int r_n_stride = L_tot + 1;
                        QC_RI_Compute_Boys_Double(F_vals, T_val, L_tot);
                        float AQ[3] = {A.x - Qx, A.y - Qy, A.z - Qz};
                        QC_RI_Compute_R_Tensor(R_vals, F_vals, alpha_pq, AQ,
                                               L_tot);

                        const double prefactor =
                            (double)cP * (double)c_mu * (double)c_nu *
                            (double)K_ket *
                            (2.0 * CONSTANT_Pi * CONSTANT_Pi *
                             sqrt(CONSTANT_Pi)) /
                            ((double)eP * (double)g_ket *
                             sqrt((double)(eP + g_ket)));

                        for (int idxP = 0; idxP < nP; idxP++)
                        {
                            int lxP, lyP, lzP;
                            QC_Get_Lxyz_Device(lP, idxP, lxP, lyP, lzP);

                            for (int idx_mu = 0; idx_mu < nmu; idx_mu++)
                            {
                                int lx_mu, ly_mu, lz_mu;
                                QC_Get_Lxyz_Device(lmu, idx_mu, lx_mu, ly_mu,
                                                   lz_mu);

                                for (int idx_nu = 0; idx_nu < nnu; idx_nu++)
                                {
                                    int lx_nu, ly_nu, lz_nu;
                                    QC_Get_Lxyz_Device(lnu, idx_nu, lx_nu,
                                                       ly_nu, lz_nu);

                                    double v_sum = 0.0;
                                    for (int t = 0; t <= lxP; t++)
                                    {
                                        const double ePx =
                                            (double)E_Px[QC_RI_E_Index(
                                                lxP, 0, t, workspace.e_dim_b,
                                                workspace.e_dim_n)];
                                        if (ePx == 0.0) continue;
                                        for (int u = 0; u <= lyP; u++)
                                        {
                                            const double ePy =
                                                (double)E_Py[QC_RI_E_Index(
                                                    lyP, 0, u,
                                                    workspace.e_dim_b,
                                                    workspace.e_dim_n)];
                                            if (ePy == 0.0) continue;
                                            for (int v = 0; v <= lzP; v++)
                                            {
                                                const double ePz =
                                                    (double)E_Pz[QC_RI_E_Index(
                                                        lzP, 0, v,
                                                        workspace.e_dim_b,
                                                        workspace.e_dim_n)];
                                                if (ePz == 0.0) continue;
                                                for (int tt = 0;
                                                     tt <= lx_mu + lx_nu; tt++)
                                                {
                                                    const double eKx = (double)
                                                        E_Kx[QC_RI_E_Index(
                                                            lx_mu, lx_nu, tt,
                                                            workspace.e_dim_b,
                                                            workspace.e_dim_n)];
                                                    if (eKx == 0.0) continue;
                                                    for (int uu = 0;
                                                         uu <= ly_mu + ly_nu;
                                                         uu++)
                                                    {
                                                        const double eKy = (double)
                                                            E_Ky[QC_RI_E_Index(
                                                                ly_mu, ly_nu,
                                                                uu,
                                                                workspace
                                                                    .e_dim_b,
                                                                workspace
                                                                    .e_dim_n)];
                                                        if (eKy == 0.0)
                                                            continue;
                                                        for (int vv = 0;
                                                             vv <=
                                                             lz_mu + lz_nu;
                                                             vv++)
                                                        {
                                                            const double eKz = (double)
                                                                E_Kz[QC_RI_E_Index(
                                                                    lz_mu,
                                                                    lz_nu, vv,
                                                                    workspace
                                                                        .e_dim_b,
                                                                    workspace
                                                                        .e_dim_n)];
                                                            if (eKz == 0.0)
                                                                continue;
                                                            const double sign =
                                                                ((tt + uu +
                                                                  vv) &
                                                                 1)
                                                                    ? -1.0
                                                                    : 1.0;
                                                            v_sum +=
                                                                ePx * ePy *
                                                                ePz * eKx *
                                                                eKy * eKz *
                                                                sign *
                                                                (double)R_vals
                                                                    [QC_RI_R_Index(
                                                                        t + tt,
                                                                        u + uu,
                                                                        v + vv,
                                                                        0,
                                                                        r_n_stride)];
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    const size_t cart_idx =
                                        ((size_t)idxP * (size_t)nmu +
                                         (size_t)idx_mu) *
                                            (size_t)nnu +
                                        (size_t)idx_nu;
                                    buf[cart_idx] += prefactor * v_sum;
                                }
                            }
                        }
                    }
                }
            }

            for (int idxP = 0; idxP < nP; idxP++)
            {
                const size_t P_idx = (size_t)(offP + idxP);
                for (int idx_mu = 0; idx_mu < nmu; idx_mu++)
                {
                    const int mu_idx = offmu + idx_mu;
                    const int mu_local = mu_idx - mu_offset_base;
                    for (int idx_nu = 0; idx_nu < nnu; idx_nu++)
                    {
                        const int nu_idx = offnu + idx_nu;
                        const int nu_local = nu_idx - nu_offset_base;
                        const size_t cart_idx =
                            ((size_t)idxP * (size_t)nmu + (size_t)idx_mu) *
                                (size_t)nnu +
                            (size_t)idx_nu;
                        const size_t idx3c =
                            (P_idx * (size_t)out_mu_dim + (size_t)mu_local) *
                                (size_t)out_nu_dim +
                            (size_t)nu_local;
                        out_eri3c[idx3c] = buf[cart_idx];

                        if (fill_symmetric && mu_sh != nu_sh)
                        {
                            const int nu_local_sym = nu_idx - mu_offset_base;
                            const int mu_local_sym = mu_idx - nu_offset_base;
                            const size_t idx3c_sym =
                                (P_idx * (size_t)out_mu_dim +
                                 (size_t)nu_local_sym) *
                                    (size_t)out_nu_dim +
                                (size_t)mu_local_sym;
                            out_eri3c[idx3c_sym] = buf[cart_idx];
                        }
                    }
                }
            }
        }
    }
}

static inline void QC_Launch_RI_3Center_Kernel(
    int threads, int n_tasks, const QC_RI_3C_TASK* tasks,
    const VECTOR* aux_centers, const int* aux_l_list, const float* aux_exps,
    const float* aux_coeffs, const int* aux_shell_offsets,
    const int* aux_shell_sizes, const int* aux_ao_offsets,
    const VECTOR* orb_centers, const int* orb_l_list, const float* orb_exps,
    const float* orb_coeffs, const int* orb_shell_offsets,
    const int* orb_shell_sizes, const int* orb_ao_offsets, int naux,
    int out_mu_dim, int out_nu_dim, int mu_offset_base, int nu_offset_base,
    bool fill_symmetric, const QC_RI_INTEGRAL_WORKSPACE& workspace,
    double* out_eri3c)
{
    if (n_tasks == 0) return;
    const int grid = (workspace.device.n_workers + threads - 1) / threads;
    Launch_Device_Kernel(
        QC_RI_3Center_Kernel, grid, threads, 0, 0, n_tasks, tasks, aux_centers,
        aux_l_list, aux_exps, aux_coeffs, aux_shell_offsets, aux_shell_sizes,
        aux_ao_offsets, orb_centers, orb_l_list, orb_exps, orb_coeffs,
        orb_shell_offsets, orb_shell_sizes, orb_ao_offsets, naux, out_mu_dim,
        out_nu_dim, mu_offset_base, nu_offset_base, fill_symmetric,
        workspace.device, out_eri3c);
}
