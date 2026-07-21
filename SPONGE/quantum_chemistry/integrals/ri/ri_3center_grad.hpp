#pragma once

// 三中心积分导数 d(P|μν)/dR 的内核
// 对每个 shell 三元组 (P_sh, mu_sh, nu_sh) 并行
// 计算 d/dA_P 和 d/dA_mu，d/dA_nu 由平移不变性得出
// 优化: primitive 循环在外，E/Boys/R 只算一次后对所有 Cartesian 分量查表收缩

#include "ri_2center_grad.hpp"

// 三中心积分导数内核。有界 worker 池按 stride 遍历辅助 shell，每个
// worker 独占动态 Cartesian/E/R/Boys workspace。
static __global__ void QC_RI_3Center_Grad_Kernel(
    const int naux_bas, const int norb_bas, const VECTOR* aux_centers,
    const int* aux_l_list, const float* aux_exps, const float* aux_coeffs,
    const int* aux_shell_offsets, const int* aux_shell_sizes,
    const int* aux_ao_offsets_cart, const int* aux_ao_offsets_sph,
    const VECTOR* orb_centers, const int* orb_l_list, const float* orb_exps,
    const float* orb_coeffs, const int* orb_shell_offsets,
    const int* orb_shell_sizes, const int* orb_ao_offsets_cart,
    const int* orb_ao_offsets_sph, int is_spherical, const float* aux_norms,
    const float* orb_norms, const float* U_aux, const float* U_orb,
    int naux_cart, int naux_sph, int nao_cart, int nao_sph,
    const double* D3_eff, const int* shell_atom_aux, const int* shell_atom_orb,
    QC_RI_DEVICE_WORKSPACE workspace, double* grad)
{
    const int nao = nao_sph;

    SIMPLE_DEVICE_FOR(worker_id, workspace.n_workers)
    {
        const int e_dim_a = workspace.e_dim_a;
        const int e_dim_b = workspace.e_dim_b;
        const int e_dim_n = workspace.e_dim_n;
        double* my_ws =
            workspace.cart + (size_t)worker_id * workspace.cart_stride;
        const size_t half_ws = workspace.cart_stride / 2;
        float* worker_e = workspace.e + (size_t)worker_id * workspace.e_stride;
        const size_t e_tensor_size =
            QC_RI_E_Tensor_Size(e_dim_a, e_dim_b, e_dim_n);
        float* E_Px = worker_e;
        float* E_Py = E_Px + e_tensor_size;
        float* E_Pz = E_Py + e_tensor_size;
        float* E_Kx = E_Pz + e_tensor_size;
        float* E_Ky = E_Kx + e_tensor_size;
        float* E_Kz = E_Ky + e_tensor_size;
        float* R_vals = workspace.r + (size_t)worker_id * workspace.r_stride;
        double* F_vals =
            workspace.boys + (size_t)worker_id * workspace.boys_stride;

        for (size_t P_sh_index = (size_t)worker_id;
             P_sh_index < (size_t)naux_bas;
             P_sh_index += (size_t)workspace.n_workers)
        {
            const int P_sh = (int)P_sh_index;
            const int lP = aux_l_list[P_sh];
            const int nP_cart = (int)QC_RI_Cartesian_Count(lP);
            const int offP_cart = aux_ao_offsets_cart[P_sh];
            const int offP_sph = aux_ao_offsets_sph[P_sh];
            const int nP_sph = 2 * lP + 1;
            const VECTOR A = aux_centers[P_sh];
            const int atom_P = shell_atom_aux[P_sh];

            for (int mu_sh = 0; mu_sh < norb_bas; mu_sh++)
            {
                const int lmu = orb_l_list[mu_sh];
                const int nmu_cart = (int)QC_RI_Cartesian_Count(lmu);
                const int offmu_cart = orb_ao_offsets_cart[mu_sh];
                const int offmu_sph = is_spherical ? orb_ao_offsets_sph[mu_sh]
                                                   : orb_ao_offsets_cart[mu_sh];
                const int nmu_sph = is_spherical ? (2 * lmu + 1) : nmu_cart;
                const VECTOR B = orb_centers[mu_sh];
                const int atom_mu = shell_atom_orb[mu_sh];

                for (int nu_sh = 0; nu_sh <= mu_sh; nu_sh++)
                {
                    const int lnu = orb_l_list[nu_sh];
                    const int nnu_cart = (int)QC_RI_Cartesian_Count(lnu);
                    const int offnu_cart = orb_ao_offsets_cart[nu_sh];
                    const int offnu_sph = is_spherical
                                              ? orb_ao_offsets_sph[nu_sh]
                                              : orb_ao_offsets_cart[nu_sh];
                    const int nnu_sph = is_spherical ? (2 * lnu + 1) : nnu_cart;
                    const VECTOR C = orb_centers[nu_sh];
                    const int atom_nu = shell_atom_orb[nu_sh];

                    const size_t n_cart_total =
                        (size_t)nP_cart * (size_t)nmu_cart * (size_t)nnu_cart;
                    double* d_cart_P = my_ws;
                    double* d_cart_mu = my_ws + half_ws;
                    for (size_t k = 0; k < n_cart_total * 3; k++)
                    {
                        d_cart_P[k] = 0.0;
                        d_cart_mu[k] = 0.0;
                    }

                    // primitive 循环 (外层)
                    for (int pP = 0; pP < aux_shell_sizes[P_sh]; pP++)
                    {
                        const float eP = aux_exps[aux_shell_offsets[P_sh] + pP];
                        const float cP =
                            aux_coeffs[aux_shell_offsets[P_sh] + pP];

                        // Bra E 系数: 算到 lP+1 以支持 d/dA_P
                        QC_RI_Compute_MD_Coeffs(E_Px, e_dim_a, e_dim_b, e_dim_n,
                                                lP + 1, 0, 0.0f, 0.0f,
                                                0.5f / eP);
                        QC_RI_Compute_MD_Coeffs(E_Py, e_dim_a, e_dim_b, e_dim_n,
                                                lP + 1, 0, 0.0f, 0.0f,
                                                0.5f / eP);
                        QC_RI_Compute_MD_Coeffs(E_Pz, e_dim_a, e_dim_b, e_dim_n,
                                                lP + 1, 0, 0.0f, 0.0f,
                                                0.5f / eP);

                        for (int p_mu = 0; p_mu < orb_shell_sizes[mu_sh];
                             p_mu++)
                        {
                            const float e_mu =
                                orb_exps[orb_shell_offsets[mu_sh] + p_mu];
                            const float c_mu =
                                orb_coeffs[orb_shell_offsets[mu_sh] + p_mu];

                            for (int p_nu = 0; p_nu < orb_shell_sizes[nu_sh];
                                 p_nu++)
                            {
                                const float e_nu =
                                    orb_exps[orb_shell_offsets[nu_sh] + p_nu];
                                const float c_nu =
                                    orb_coeffs[orb_shell_offsets[nu_sh] + p_nu];

                                const float g_ket = e_mu + e_nu;
                                const float BC2 = (B.x - C.x) * (B.x - C.x) +
                                                  (B.y - C.y) * (B.y - C.y) +
                                                  (B.z - C.z) * (B.z - C.z);
                                const float K_ket =
                                    expf(-e_mu * e_nu / g_ket * BC2);

                                const float Qx =
                                    (e_mu * B.x + e_nu * C.x) / g_ket;
                                const float Qy =
                                    (e_mu * B.y + e_nu * C.y) / g_ket;
                                const float Qz =
                                    (e_mu * B.z + e_nu * C.z) / g_ket;

                                // Ket E 系数: 算到 (lmu+1, lnu+1) 以支持
                                // d/dA_mu
                                QC_RI_Compute_MD_Coeffs(
                                    E_Kx, e_dim_a, e_dim_b, e_dim_n, lmu + 1,
                                    lnu + 1, Qx - B.x, Qx - C.x, 0.5f / g_ket);
                                QC_RI_Compute_MD_Coeffs(
                                    E_Ky, e_dim_a, e_dim_b, e_dim_n, lmu + 1,
                                    lnu + 1, Qy - B.y, Qy - C.y, 0.5f / g_ket);
                                QC_RI_Compute_MD_Coeffs(
                                    E_Kz, e_dim_a, e_dim_b, e_dim_n, lmu + 1,
                                    lnu + 1, Qz - B.z, Qz - C.z, 0.5f / g_ket);

                                // Boys + R 张量: L_tot+1 阶
                                const float alpha_pq =
                                    eP * g_ket / (eP + g_ket);
                                const float AQ2 = (A.x - Qx) * (A.x - Qx) +
                                                  (A.y - Qy) * (A.y - Qy) +
                                                  (A.z - Qz) * (A.z - Qz);
                                const float T_val = alpha_pq * AQ2;
                                const int L_tot = lP + lmu + lnu;
                                const int R_order = L_tot + 1;
                                const int R_n_stride = R_order + 1;
                                QC_RI_Compute_Boys_Double(F_vals, T_val,
                                                          R_order);
                                float AQ[3] = {A.x - Qx, A.y - Qy, A.z - Qz};
                                QC_RI_Compute_R_Tensor(R_vals, F_vals, alpha_pq,
                                                       AQ, R_order);

                                const double prefactor =
                                    (double)cP * (double)c_mu * (double)c_nu *
                                    (double)K_ket *
                                    (2.0 * CONSTANT_Pi * CONSTANT_Pi *
                                     sqrt(CONSTANT_Pi)) /
                                    ((double)eP * (double)g_ket *
                                     sqrt((double)(eP + g_ket)));

                                // 收缩函数
                                auto contract_3c =
                                    [&](int axP, int ayP, int azP, int ax_mu,
                                        int ay_mu, int az_mu, int ax_nu,
                                        int ay_nu, int az_nu) -> double
                                {
                                    if (axP < 0 || ayP < 0 || azP < 0 ||
                                        ax_mu < 0 || ay_mu < 0 || az_mu < 0 ||
                                        ax_nu < 0 || ay_nu < 0 || az_nu < 0)
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
                                                for (int tt = 0;
                                                     tt <= ax_mu + ax_nu; tt++)
                                                {
                                                    double eKx = (double)
                                                        E_Kx[QC_RI_E_Index(
                                                            ax_mu, ax_nu, tt,
                                                            e_dim_b, e_dim_n)];
                                                    if (eKx == 0.0) continue;
                                                    for (int uu = 0;
                                                         uu <= ay_mu + ay_nu;
                                                         uu++)
                                                    {
                                                        double eKy = (double)
                                                            E_Ky[QC_RI_E_Index(
                                                                ay_mu, ay_nu,
                                                                uu, e_dim_b,
                                                                e_dim_n)];
                                                        if (eKy == 0.0)
                                                            continue;
                                                        for (int vv = 0;
                                                             vv <=
                                                             az_mu + az_nu;
                                                             vv++)
                                                        {
                                                            double eKz = (double)
                                                                E_Kz[QC_RI_E_Index(
                                                                    az_mu,
                                                                    az_nu, vv,
                                                                    e_dim_b,
                                                                    e_dim_n)];
                                                            if (eKz == 0.0)
                                                                continue;
                                                            double sign =
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
                                                                        R_n_stride)];
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    return v_sum;
                                };

                                // Cartesian 分量循环 (内层)
                                for (int idxP = 0; idxP < nP_cart; idxP++)
                                {
                                    int lxP, lyP, lzP;
                                    QC_Get_Lxyz_Device(lP, idxP, lxP, lyP, lzP);

                                    for (int idx_mu = 0; idx_mu < nmu_cart;
                                         idx_mu++)
                                    {
                                        int lx_mu, ly_mu, lz_mu;
                                        QC_Get_Lxyz_Device(lmu, idx_mu, lx_mu,
                                                           ly_mu, lz_mu);

                                        for (int idx_nu = 0; idx_nu < nnu_cart;
                                             idx_nu++)
                                        {
                                            int lx_nu, ly_nu, lz_nu;
                                            QC_Get_Lxyz_Device(lnu, idx_nu,
                                                               lx_nu, ly_nu,
                                                               lz_nu);

                                            const size_t cidx =
                                                ((size_t)idxP * nmu_cart +
                                                 (size_t)idx_mu) *
                                                    nnu_cart +
                                                (size_t)idx_nu;

                                            // d(P|μν)/dA_P
                                            {
                                                double dx =
                                                    2.0 * (double)eP *
                                                    contract_3c(
                                                        lxP + 1, lyP, lzP,
                                                        lx_mu, ly_mu, lz_mu,
                                                        lx_nu, ly_nu, lz_nu);
                                                if (lxP > 0)
                                                    dx -= (double)lxP *
                                                          contract_3c(
                                                              lxP - 1, lyP, lzP,
                                                              lx_mu, ly_mu,
                                                              lz_mu, lx_nu,
                                                              ly_nu, lz_nu);
                                                double dy =
                                                    2.0 * (double)eP *
                                                    contract_3c(
                                                        lxP, lyP + 1, lzP,
                                                        lx_mu, ly_mu, lz_mu,
                                                        lx_nu, ly_nu, lz_nu);
                                                if (lyP > 0)
                                                    dy -= (double)lyP *
                                                          contract_3c(
                                                              lxP, lyP - 1, lzP,
                                                              lx_mu, ly_mu,
                                                              lz_mu, lx_nu,
                                                              ly_nu, lz_nu);
                                                double dz =
                                                    2.0 * (double)eP *
                                                    contract_3c(
                                                        lxP, lyP, lzP + 1,
                                                        lx_mu, ly_mu, lz_mu,
                                                        lx_nu, ly_nu, lz_nu);
                                                if (lzP > 0)
                                                    dz -= (double)lzP *
                                                          contract_3c(
                                                              lxP, lyP, lzP - 1,
                                                              lx_mu, ly_mu,
                                                              lz_mu, lx_nu,
                                                              ly_nu, lz_nu);
                                                d_cart_P[cidx * 3 + 0] +=
                                                    prefactor * dx;
                                                d_cart_P[cidx * 3 + 1] +=
                                                    prefactor * dy;
                                                d_cart_P[cidx * 3 + 2] +=
                                                    prefactor * dz;
                                            }

                                            // d(P|μν)/dA_mu
                                            {
                                                double dx =
                                                    2.0 * (double)e_mu *
                                                    contract_3c(
                                                        lxP, lyP, lzP,
                                                        lx_mu + 1, ly_mu, lz_mu,
                                                        lx_nu, ly_nu, lz_nu);
                                                if (lx_mu > 0)
                                                    dx -= (double)lx_mu *
                                                          contract_3c(
                                                              lxP, lyP, lzP,
                                                              lx_mu - 1, ly_mu,
                                                              lz_mu, lx_nu,
                                                              ly_nu, lz_nu);
                                                double dy =
                                                    2.0 * (double)e_mu *
                                                    contract_3c(
                                                        lxP, lyP, lzP, lx_mu,
                                                        ly_mu + 1, lz_mu, lx_nu,
                                                        ly_nu, lz_nu);
                                                if (ly_mu > 0)
                                                    dy -= (double)ly_mu *
                                                          contract_3c(
                                                              lxP, lyP, lzP,
                                                              lx_mu, ly_mu - 1,
                                                              lz_mu, lx_nu,
                                                              ly_nu, lz_nu);
                                                double dz =
                                                    2.0 * (double)e_mu *
                                                    contract_3c(
                                                        lxP, lyP, lzP, lx_mu,
                                                        ly_mu, lz_mu + 1, lx_nu,
                                                        ly_nu, lz_nu);
                                                if (lz_mu > 0)
                                                    dz -= (double)lz_mu *
                                                          contract_3c(
                                                              lxP, lyP, lzP,
                                                              lx_mu, ly_mu,
                                                              lz_mu - 1, lx_nu,
                                                              ly_nu, lz_nu);
                                                d_cart_mu[cidx * 3 + 0] +=
                                                    prefactor * dx;
                                                d_cart_mu[cidx * 3 + 1] +=
                                                    prefactor * dy;
                                                d_cart_mu[cidx * 3 + 2] +=
                                                    prefactor * dz;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Cart2Sph + 归一化 + D3_eff 收缩 → atomicAdd
                    for (int ps = 0; ps < nP_sph; ps++)
                    {
                        const int P_sph = offP_sph + ps;
                        const double normP_val = (double)aux_norms[P_sph];

                        for (int ms = 0; ms < nmu_sph; ms++)
                        {
                            const int mu_sph = offmu_sph + ms;
                            const double normMu = (double)orb_norms[mu_sph];

                            for (int ns = 0; ns < nnu_sph; ns++)
                            {
                                const int nu_sph = offnu_sph + ns;
                                const double normNu = (double)orb_norms[nu_sph];

                                double d_sph_P[3] = {0.0, 0.0, 0.0};
                                double d_sph_mu[3] = {0.0, 0.0, 0.0};

                                for (int pc = 0; pc < nP_cart; pc++)
                                {
                                    double u_p = (double)
                                        U_aux[(offP_cart + pc) * naux_sph +
                                              P_sph];
                                    if (u_p == 0.0) continue;

                                    for (int mc = 0; mc < nmu_cart; mc++)
                                    {
                                        double u_m =
                                            is_spherical
                                                ? (double)
                                                      U_orb[(offmu_cart + mc) *
                                                                nao_sph +
                                                            mu_sph]
                                                : (mc == ms ? 1.0 : 0.0);
                                        if (u_m == 0.0) continue;

                                        for (int nc = 0; nc < nnu_cart; nc++)
                                        {
                                            double u_n =
                                                is_spherical
                                                    ? (double)
                                                          U_orb[(offnu_cart +
                                                                 nc) *
                                                                    nao_sph +
                                                                nu_sph]
                                                    : (nc == ns ? 1.0 : 0.0);
                                            if (u_n == 0.0) continue;

                                            double w = u_p * u_m * u_n;
                                            const size_t cidx =
                                                ((size_t)pc * nmu_cart +
                                                 (size_t)mc) *
                                                    nnu_cart +
                                                (size_t)nc;
                                            d_sph_P[0] +=
                                                w * d_cart_P[cidx * 3 + 0];
                                            d_sph_P[1] +=
                                                w * d_cart_P[cidx * 3 + 1];
                                            d_sph_P[2] +=
                                                w * d_cart_P[cidx * 3 + 2];
                                            d_sph_mu[0] +=
                                                w * d_cart_mu[cidx * 3 + 0];
                                            d_sph_mu[1] +=
                                                w * d_cart_mu[cidx * 3 + 1];
                                            d_sph_mu[2] +=
                                                w * d_cart_mu[cidx * 3 + 2];
                                        }
                                    }
                                }

                                double norm_all = normP_val * normMu * normNu;
                                for (int d = 0; d < 3; d++)
                                {
                                    d_sph_P[d] *= norm_all;
                                    d_sph_mu[d] *= norm_all;
                                }

                                double dens =
                                    D3_eff[(long long)P_sph * nao * nao +
                                           (long long)mu_sph * nao + nu_sph];
                                if (mu_sh != nu_sh)
                                    dens +=
                                        D3_eff[(long long)P_sph * nao * nao +
                                               (long long)nu_sph * nao +
                                               mu_sph];

                                for (int d = 0; d < 3; d++)
                                {
                                    double g_P = dens * d_sph_P[d];
                                    double g_mu = dens * d_sph_mu[d];
                                    double g_nu = -(g_P + g_mu);

                                    atomicAdd(&grad[atom_P * 3 + d], g_P);
                                    atomicAdd(&grad[atom_mu * 3 + d], g_mu);
                                    atomicAdd(&grad[atom_nu * 3 + d], g_nu);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
