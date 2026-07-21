#pragma once

// RI (Density Fitting) 解析梯度驱动
//
// RI-J 梯度:
//   dE_J/dR = Σ_{P,μν} g_P D_μν d(P|μν)/dR
//           - (1/2) Σ_{PQ} g_P g_Q d(P|Q)/dR
//
// RI-K 梯度:
//   3c 部分: dE_K/dR|_{3c} = -exx Σ_{Q,μν} D3_K[Q,μν] d(Q|μν)/dR
//     D3_K[Q,μ,λ] = Σ_{P,ν,k} M^{-1/2}[QP] (B L)[P,ν,k] L[λ,k] D[μ,ν]
//   2c 部分: dE_K/dR|_{2c} = Σ_{PQ} D2_K[P,Q] d(P|Q)/dR
//     D2_K = U (F ⊙ (U^T Z_K U)) U^T   (Daleckii-Kreĭn 矩阵函数导数)
//     Z_K[P,Q] = -Σ_{m,n,l,k} (Q|ml) L[l,k] (B L)[P,n,k] D[m,n]
//     F[k,l] = { -1/2 λ_k^{-3/2}                    if k=l
//              { (λ_k^{-1/2} - λ_l^{-1/2})/(λ_k-λ_l) if k≠l

#include <cmath>
#include <cstring>
#include <vector>

#include "../integrals/ri/ri_2center_grad.hpp"
#include "../integrals/ri/ri_3center_grad.hpp"
#include "ri_metric_response.hpp"

// cblas.h 已通过 device_backend/cpu_api.h 引入，此处仅做特性检测
#if !defined(USE_GPU) && (defined(USE_MKL) || defined(USE_OPENBLAS))
#define RI_GRAD_HAS_BLAS 1
#else
#define RI_GRAD_HAS_BLAS 0
#endif

// 共用辅助函数

// 构建 RI-J 的三中心有效密度 D3_J[P,mu,nu] = g_P D_mu,nu。
static inline void QC_Init_D3_J(const int nao, const int naux,
                                const double* g_vec,
                                const float* coulomb_density,
                                std::vector<double>& D3_eff)
{
    const long long nao2 = (long long)nao * nao;
    D3_eff.resize((size_t)naux * nao2);
    std::vector<double> D_d(nao2);
    for (long long i = 0; i < nao2; i++) D_d[i] = (double)coulomb_density[i];
    for (int P = 0; P < naux; P++)
    {
        const double gP = g_vec[P];
        double* dst = D3_eff.data() + (long long)P * nao2;
        for (long long mn = 0; mn < nao2; mn++) dst[mn] = gP * D_d[mn];
    }
}

// 将一个自旋通道的 RI-K 贡献累加到 D3_eff。RHF 的 spin_density
// 已含占据因子 2；UHF 的 alpha/beta 密度分别传入，不得与总密度混用。
static inline void QC_Accumulate_D3_K_Channel(
    const int nao, const int naux, const float* spin_density,
    const double* metric_inv_sqrt,
    const float* B_factor,       // [M × rank] column-major, M=naux*nao
    const float* density_factor, // [nao × nao] row-major, P_sigma=L L^T
    const int factor_rank, const float exx_fraction,
    std::vector<double>& D3_eff)
{
    const long long nao2 = (long long)nao * nao;
    // D3_K[Q,μ,λ] = Σ_{P,ν,k} M^{-1/2}[QP] (B L)[P,ν,k] L[λ,k] D[μ,ν]
    // 分解:
    //   X_P[ν,λ] = Σ_i B_occ[P,ν,i] × C[λ,i]   (per-P matmul)
    //   Y_P[μ,λ] = Σ_ν D[μ,ν] × X_P[ν,λ]       (per-P matmul)
    //   D3_K[Q,ml] = -exx × Σ_P M^{-1/2}[Q,P] × Y[P,ml]  (matmul)
    if (exx_fraction != 0.0f && factor_rank > 0 && B_factor != nullptr &&
        density_factor != nullptr && spin_density != nullptr)
    {
        const int M_dim = naux * nao;
        const double neg_exx = -(double)exx_fraction;
        const double one_d = 1.0;

#if RI_GRAD_HAS_BLAS
        // BLAS 路径: 使用 cblas_dgemm 加速矩阵乘法

        // 准备 double 精度缓冲
        std::vector<double> D_d(nao2);
        for (long long i = 0; i < nao2; i++) D_d[i] = (double)spin_density[i];

        // Double representation of L[:, :rank].
        std::vector<double> C_d((size_t)nao * factor_rank);
        for (int lam = 0; lam < nao; lam++)
            for (int i = 0; i < factor_rank; i++)
                C_d[lam * factor_rank + i] =
                    (double)density_factor[lam * nao + i];

        // B_P_d: 行优先 [nao, nocc] 连续缓冲
        std::vector<double> B_P_d((size_t)nao * factor_rank);
        std::vector<double> X_P((size_t)nao * nao);
        std::vector<double> Y((size_t)naux * nao2, 0.0);

        for (int P = 0; P < naux; P++)
        {
            // 拷贝 B_P → 行优先 [nao, nocc]: B_P_d[ν*nocc+i]
            for (int i = 0; i < factor_rank; i++)
                for (int nu = 0; nu < nao; nu++)
                    B_P_d[(size_t)nu * factor_rank + i] = (double)
                        B_factor[(size_t)(P * nao + nu) + (size_t)M_dim * i];

            // X_P[ν,λ] = Σ_i B_P[ν,i] × C[λ,i]
            // 行优先: X_P(nao×nao) = B_P_d(nao×nocc) @ C_d^T(nocc×nao)
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nao, nao,
                        factor_rank, 1.0, B_P_d.data(), factor_rank,
                        C_d.data(), factor_rank, 0.0,
                        X_P.data(), nao);

            // Y_P[μ,λ] = Σ_ν D[μ,ν] × X_P[ν,λ]
            // 行优先: Y_P(nao×nao) = D(nao×nao) @ X_P(nao×nao)
            double* Y_P = Y.data() + (long long)P * nao2;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nao, nao,
                        nao, 1.0, D_d.data(), nao, X_P.data(), nao, 0.0, Y_P,
                        nao);
        }

        // Step 3: D3_eff += neg_exx × M^{-1/2} @ Y
        // 行优先: D3_eff[naux, nao²] += neg_exx * M[naux,naux] @ Y[naux,nao²]
        // 列优先视角: D3^T[nao²,naux] += neg_exx * Y^T[nao²,naux] *
        // M^T[naux,naux] M^{-1/2} 对称 → M^T = M
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, naux, (int)nao2,
                    naux, neg_exx, metric_inv_sqrt, naux, Y.data(), (int)nao2,
                    one_d, D3_eff.data(), (int)nao2);

#else
        // 回退: 标量循环 (无 BLAS 可用)
        std::vector<double> Y((size_t)naux * nao2, 0.0);
        std::vector<double> X_P((size_t)nao * nao, 0.0);
        for (int P = 0; P < naux; P++)
        {
            std::fill(X_P.begin(), X_P.end(), 0.0);
            for (int nu = 0; nu < nao; nu++)
                for (int i = 0; i < factor_rank; i++)
                {
                    double bval = (double)
                        B_factor[(size_t)(P * nao + nu) + (size_t)M_dim * i];
                    if (bval == 0.0) continue;
                    for (int lam = 0; lam < nao; lam++)
                        X_P[nu * nao + lam] +=
                            bval * (double)density_factor[lam * nao + i];
                }
            double* Y_P = Y.data() + (long long)P * nao2;
            for (int mu = 0; mu < nao; mu++)
                for (int nu = 0; nu < nao; nu++)
                {
                    double d = (double)spin_density[mu * nao + nu];
                    if (d == 0.0) continue;
                    for (int lam = 0; lam < nao; lam++)
                        Y_P[mu * nao + lam] += d * X_P[nu * nao + lam];
                }
        }
        for (int Q = 0; Q < naux; Q++)
            for (int P = 0; P < naux; P++)
            {
                double w = neg_exx * metric_inv_sqrt[(size_t)Q * naux + P];
                if (w == 0.0) continue;
                for (long long mn = 0; mn < nao2; mn++)
                    D3_eff[(long long)Q * nao2 + mn] +=
                        w * Y[(long long)P * nao2 + mn];
            }
#endif
    }
}

// D2_J 初始化: D2_eff[P,Q] = -0.5 * g_P * g_Q
static inline void QC_Init_D2_J(const int naux, const double* g_vec,
                                std::vector<double>& D2_eff)
{
    const size_t naux2 = (size_t)naux * naux;
    D2_eff.assign(naux2, 0.0);
    for (int P = 0; P < naux; P++)
        for (int Q = 0; Q < naux; Q++)
            D2_eff[(size_t)P * naux + Q] = -0.5 * g_vec[P] * g_vec[Q];
}

// D2_K 的 Daleckii-Kreĭn 矩阵函数导数: D2_K = U (F ⊙ (U^T Z_K U)) U^T
// 结果累加到 D2_eff。
static inline void QC_Accumulate_D2_K_DaleckiiKrein(
    const int naux, const int naux_eff, const double* Z_K, const double* eigval,
    const double* eigvec, const double exx_fraction,
    std::vector<double>& D2_eff)
{
    const size_t naux2 = (size_t)naux * naux;
    std::vector<double> UZU(naux2, 0.0);
    std::vector<double> tmp(naux2, 0.0);

#if RI_GRAD_HAS_BLAS
    // BLAS 路径: U^T Z_K U → D2_K = U (F ⊙ UZU) U^T
    // Z_K 行优先 [naux,naux], eigvec 列优先 [naux,naux]

    // tmp = Z_K @ U: Z_K 行优先, U 列优先
    // 行优先 Z_K × 列优先 U: 用 cblas_dgemm(RowMajor, N, N, ...)
    // Z_K[i,k] * U_col[k,j] = Z_K[i,k] * eigvec[k + j*naux]
    // cblas_dgemm(CblasRowMajor, N, N, naux, naux, naux, 1, Z_K, naux,
    // eigvec(列→行等价T), naux, 0, tmp, naux) eigvec 列优先 [naux,naux]
    // 用行优先读 = eigvec^T 所以 Z_K @ U = Z_K(行) @ eigvec^T(行)^T → RowMajor,
    // N, N 不对 正确: tmp[i,j] = Σ_k Z_K[i,k] * U[k,j] U[k,j] = eigvec[k +
    // j*naux] (col-major) = eigvec_rowmajor^T → tmp = Z_K @ U, 其中 U 列优先
    // CblasRowMajor: C = A × B, A 行优先 [M,K], B 行优先 [K,N]
    // 但 U 是列优先, 行优先读就是 U^T, 需要转置标记
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, naux, naux, naux, 1.0,
                Z_K, naux, eigvec, naux, 0.0, tmp.data(), naux);
    // 上面: RowMajor, N, T → tmp[i,j] = Σ_k Z_K[i,k] * eigvec_row[j,k]
    //   = Σ_k Z_K[i,k] * eigvec[j*naux+k] ← 不对，应该是 eigvec[k+j*naux]
    // eigvec 列优先: eigvec[k + j*naux]
    // CblasRowMajor 把指针读为行优先: B_row[j,k] = eigvec[j*naux + k]
    // CblasTrans → B^T[k,j] = B_row[j,k] = eigvec[j*naux + k]
    // tmp[i,j] = Σ_k Z_K[i,k] * B^T[k,j] = Σ_k Z_K[i,k] * eigvec[j*naux+k]
    // 但我们要 Σ_k Z_K[i,k] * eigvec[k + j*naux]
    // eigvec[j*naux+k] = eigvec[k + j*naux] ← 只有 naux*j+k vs k+j*naux →
    // 是同一个东西！ ✓ 正确

    // UZU = U^T @ tmp: UZU[i,j] = Σ_k U[k,i] * tmp[k,j]
    //   = Σ_k eigvec[k + i*naux] * tmp[k,j]
    // RowMajor: U_row(行优先读 eigvec) = eigvec[i*naux+k], 即 U^T[i,k]
    // UZU = U_row @ tmp (RowMajor, NoTrans, NoTrans)
    // UZU[i,j] = Σ_k U_row[i,k] * tmp[k,j] = Σ_k eigvec[i*naux+k] * tmp[k,j]
    //   = Σ_k eigvec[k + i*naux] * tmp[k,j] ✓
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, naux, naux, naux,
                1.0, eigvec, naux, tmp.data(), naux, 0.0, UZU.data(), naux);

#else
    // tmp = Z_K @ U
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            for (int k = 0; k < naux; k++)
                tmp[(size_t)i * naux + j] +=
                    Z_K[(size_t)i * naux + k] *
                    eigvec[(size_t)k + (size_t)j * naux];
    // UZU = U^T @ tmp
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            for (int k = 0; k < naux; k++)
                UZU[(size_t)i * naux + j] +=
                    eigvec[(size_t)k + (size_t)i * naux] *
                    tmp[(size_t)k * naux + j];
#endif

    // F ⊙ UZU (Hadamard 积).  Metric_InvSqrt uses a truncated spectral
    // function: f(lambda)=0 in the discarded subspace and lambda^(-1/2) in
    // the retained subspace.  Its Loewner matrix must use the same cutoff.
    //
    // For two retained eigenvalues, evaluating the divided difference as
    //     (x^(-1/2) - y^(-1/2)) / (x - y)
    // is both cancellation-prone and undefined for a degenerate pair.  The
    // algebraically equivalent expression below is finite at x == y and
    // naturally reduces to f'(x) = -1/(2*x^(3/2)).
    for (int k = 0; k < naux; k++)
    {
        for (int l = 0; l < naux; l++)
        {
            const double f =
                QC_RI_Truncated_InvSqrt_Loewner(eigval, naux, naux_eff, k, l);
            UZU[(size_t)k * naux + l] *= exx_fraction * f;
        }
    }

#if RI_GRAD_HAS_BLAS
    // D2_K = U @ (F⊙UZU) @ U^T
    // tmp = (F⊙UZU) @ U^T: tmp[i,j] = Σ_k UZU[i,k] * U[j,k]^T
    //   = Σ_k UZU[i,k] * eigvec[j + k*naux]
    // RowMajor 读 eigvec: eigvec_row[k,j] = eigvec[k*naux+j]
    // tmp = UZU @ eigvec_row^T → RowMajor, N, T
    // tmp[i,j] = Σ_k UZU[i,k] * eigvec_row[j,k] = Σ_k UZU[i,k] *
    // eigvec[j*naux+k]
    //   = Σ_k UZU[i,k] * eigvec[k*naux+j]? No.
    // eigvec_row[j,k] = eigvec[j*naux+k] (RowMajor 解释)
    // 我们要 eigvec[j + k*naux] ← 不等于 eigvec[j*naux+k]
    // 用 NoTrans: tmp = UZU @ eigvec_row(NoTrans)
    // tmp[i,j] = Σ_k UZU[i,k] * eigvec_row[k,j] = Σ_k UZU[i,k] *
    // eigvec[k*naux+j]
    //   = Σ_k UZU[i,k] * eigvec[j + k*naux] ✓ (因为 k*naux+j = j + k*naux)
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, naux, naux, naux,
                1.0, UZU.data(), naux, eigvec, naux, 0.0, tmp.data(), naux);

    // D2_eff += U @ tmp: D2_eff[i,j] += Σ_k U[i,k] * tmp[k,j]
    //   = Σ_k eigvec[i + k*naux] * tmp[k,j]
    //   = Σ_k eigvec_row[k,i]^T * tmp[k,j]
    // D2_eff += eigvec_row^T @ tmp → RowMajor, Trans, NoTrans
    // 但 eigvec_row^T[i,k] = eigvec_row[k,i] = eigvec[k*naux+i] = eigvec[i +
    // k*naux] ✓
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, naux, naux, naux, 1.0,
                eigvec, naux, tmp.data(), naux, 1.0, D2_eff.data(), naux);
#else
    // D2_K = U @ (F⊙UZU) @ U^T
    std::fill(tmp.begin(), tmp.end(), 0.0);
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            for (int k = 0; k < naux; k++)
                tmp[(size_t)i * naux + j] +=
                    UZU[(size_t)i * naux + k] *
                    eigvec[(size_t)j + (size_t)k * naux];
    for (int i = 0; i < naux; i++)
        for (int j = 0; j < naux; j++)
            for (int k = 0; k < naux; k++)
                D2_eff[(size_t)i * naux + j] +=
                    eigvec[(size_t)i + (size_t)k * naux] *
                    tmp[(size_t)k * naux + j];
#endif
}

// 上传 D2/D3 有效密度到 device，分配 workspace，启动 2c/3c 梯度内核。
// 基组指针均为 device 指针（已在 SCF 初始化时分配）。
// 仅 h_U_aux/h_U_orb 和 D2/D3 需要临时上传。
static inline bool QC_Launch_RI_Grad_Kernels(
    const QC_MOLECULE& mol, const QC_RI_WORKSPACE& ri, const float* d_orb_norms,
    const QC_GRAD_WORKSPACE& grad_ws, CONTROLLER* controller,
    const std::vector<double>& D2_eff, const std::vector<double>& D3_eff)
{
    // 局部别名: 保持函数体不变
    const int naux_bas = ri.naux_bas;
    const int norb_bas = mol.nbas;
    const bool is_spherical = mol.is_spherical;
    const int naux_cart = ri.naux_cart;
    const int naux = ri.naux;
    const int nao_cart = mol.nao_cart;
    const int nao = mol.nao;
    const VECTOR* d_aux_centers = ri.d_aux_centers;
    const int* d_aux_l_list = ri.d_aux_l_list;
    const float* d_aux_exps = ri.d_aux_exps;
    const float* d_aux_coeffs = ri.d_aux_coeffs;
    const int* d_aux_shell_offsets = ri.d_aux_shell_offsets;
    const int* d_aux_shell_sizes = ri.d_aux_shell_sizes;
    const int* d_aux_ao_offsets_cart = ri.d_aux_ao_offsets;
    const int* d_aux_ao_offsets_sph = ri.d_aux_ao_offsets_sph;
    const float* d_aux_norms = ri.d_aux_norms;
    const VECTOR* d_orb_centers = mol.d_centers;
    const int* d_orb_l_list = mol.d_l_list;
    const float* d_orb_exps = mol.d_exps;
    const float* d_orb_coeffs = mol.d_coeffs;
    const int* d_orb_shell_offsets = mol.d_shell_offsets;
    const int* d_orb_shell_sizes = mol.d_shell_sizes;
    const int* d_orb_ao_offsets_cart = mol.d_ao_offsets;
    const int* d_orb_ao_offsets_sph = mol.d_ao_offsets_sph;
    const float* h_U_aux = ri.h_U_aux.data();
    const float* h_U_orb = is_spherical ? ri.h_U_orb.data() : nullptr;
    const int* d_shell_atom_aux = grad_ws.d_shell_atom_aux;
    const int* d_shell_atom_orb = grad_ws.d_shell_atom;
    double* d_grad = grad_ws.d_grad;
    const int threads = 64;
    int max_aux_l = 0;
    for (int l : ri.h_aux_l_list) max_aux_l = std::max(max_aux_l, l);
    int max_orb_l = 0;
    for (int l : mol.h_l_list) max_orb_l = std::max(max_orb_l, l);

    QC_RI_INTEGRAL_WORKSPACE workspace_2c;
    QC_RI_INTEGRAL_WORKSPACE workspace_3c;
    if (naux_bas <= 0 || norb_bas <= 0 || naux_cart <= 0 || naux <= 0 ||
        nao_cart <= 0 || nao <= 0 ||
        !QC_RI_Build_2Center_Gradient_Workspace_Layout(naux_bas, max_aux_l,
                                                       &workspace_2c) ||
        !QC_RI_Build_3Center_Gradient_Workspace_Layout(
            naux_bas, max_aux_l, max_orb_l, &workspace_3c))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_RI_Gradient",
            "Reason:\n    invalid or overflowing RI gradient workspace "
            "dimensions: naux_bas=%d, norb_bas=%d, auxiliary/orbital angular "
            "momenta=%d/%d\n",
            naux_bas, norb_bas, max_aux_l, max_orb_l);
        return false;
    }

    size_t U_aux_count = 0;
    size_t U_orb_count = 0;
    size_t U_aux_bytes = 0;
    size_t U_orb_bytes = 0;
    size_t D2_bytes = 0;
    size_t D3_bytes = 0;
    size_t expected_D2_size = 0;
    size_t expected_D3_size = 0;
    if (!QC_RI_Checked_Mul_Size((size_t)naux_cart, (size_t)naux,
                                &U_aux_count) ||
        !QC_RI_Checked_Bytes(U_aux_count, sizeof(float), &U_aux_bytes) ||
        (h_U_orb != NULL &&
         (!QC_RI_Checked_Mul_Size((size_t)nao_cart, (size_t)nao,
                                  &U_orb_count) ||
          !QC_RI_Checked_Bytes(U_orb_count, sizeof(float), &U_orb_bytes))) ||
        !QC_RI_Checked_Mul_Size((size_t)naux, (size_t)naux,
                                &expected_D2_size) ||
        !QC_RI_Checked_Mul_Size((size_t)naux, (size_t)nao, &expected_D3_size) ||
        !QC_RI_Checked_Mul_Size(expected_D3_size, (size_t)nao,
                                &expected_D3_size) ||
        !QC_RI_Checked_Bytes(D2_eff.size(), sizeof(double), &D2_bytes) ||
        !QC_RI_Checked_Bytes(D3_eff.size(), sizeof(double), &D3_bytes))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_RI_Gradient",
            "Reason:\n    RI gradient upload buffer size overflows size_t\n");
        return false;
    }
    if (D2_eff.size() != expected_D2_size ||
        D3_eff.size() != expected_D3_size)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Build_RI_Gradient",
            "Reason:\n    RI effective-density dimensions are inconsistent: "
            "D2=%zu (expected %zu), D3=%zu (expected %zu)\n",
            D2_eff.size(), expected_D2_size, D3_eff.size(), expected_D3_size);
        return false;
    }

    float* d_U_aux = NULL;
    float* d_U_orb = NULL;
    double* d_D2 = NULL;
    double* d_D3 = NULL;
    auto free_uploads = [&]()
    {
        if (d_D3 != NULL) deviceFree(d_D3);
        if (d_D2 != NULL) deviceFree(d_D2);
        if (d_U_orb != NULL) deviceFree(d_U_orb);
        if (d_U_aux != NULL) deviceFree(d_U_aux);
        d_D3 = NULL;
        d_D2 = NULL;
        d_U_orb = NULL;
        d_U_aux = NULL;
    };

    if (h_U_aux == NULL || U_aux_bytes == 0 || D2_bytes == 0 || D3_bytes == 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Build_RI_Gradient",
            "Reason:\n    missing RI gradient upload data\n");
        return false;
    }
    if (!Device_Malloc_Safely((void**)&d_U_aux, U_aux_bytes)) return false;
    deviceMemcpy(d_U_aux, h_U_aux, U_aux_bytes, deviceMemcpyHostToDevice);
    if (h_U_orb != NULL && U_orb_bytes > 0)
    {
        if (!Device_Malloc_Safely((void**)&d_U_orb, U_orb_bytes))
        {
            free_uploads();
            return false;
        }
        deviceMemcpy(d_U_orb, h_U_orb, U_orb_bytes, deviceMemcpyHostToDevice);
    }

    if (!Device_Malloc_Safely((void**)&d_D2, D2_bytes) ||
        !QC_RI_Allocate_Integral_Workspace(&workspace_2c))
    {
        QC_RI_Free_Integral_Workspace(&workspace_2c);
        free_uploads();
        return false;
    }
    deviceMemcpy(d_D2, D2_eff.data(), D2_bytes, deviceMemcpyHostToDevice);

    const int grid_2c = (workspace_2c.device.n_workers + threads - 1) / threads;
    Launch_Device_Kernel(QC_RI_2Center_Grad_Kernel, grid_2c, threads, 0, 0,
                         naux_bas, d_aux_centers, d_aux_l_list, d_aux_exps,
                         d_aux_coeffs, d_aux_shell_offsets, d_aux_shell_sizes,
                         d_aux_ao_offsets_cart, d_aux_ao_offsets_sph,
                         d_aux_norms, d_U_aux, naux_cart, naux, d_D2,
                         d_shell_atom_aux, workspace_2c.device, d_grad);

    QC_RI_Free_Integral_Workspace(&workspace_2c);
    deviceFree(d_D2);
    d_D2 = NULL;

    if (!Device_Malloc_Safely((void**)&d_D3, D3_bytes) ||
        !QC_RI_Allocate_Integral_Workspace(&workspace_3c))
    {
        QC_RI_Free_Integral_Workspace(&workspace_3c);
        free_uploads();
        return false;
    }
    deviceMemcpy(d_D3, D3_eff.data(), D3_bytes, deviceMemcpyHostToDevice);

    const int grid_3c = (workspace_3c.device.n_workers + threads - 1) / threads;
    Launch_Device_Kernel(
        QC_RI_3Center_Grad_Kernel, grid_3c, threads, 0, 0, naux_bas, norb_bas,
        d_aux_centers, d_aux_l_list, d_aux_exps, d_aux_coeffs,
        d_aux_shell_offsets, d_aux_shell_sizes, d_aux_ao_offsets_cart,
        d_aux_ao_offsets_sph, d_orb_centers, d_orb_l_list, d_orb_exps,
        d_orb_coeffs, d_orb_shell_offsets, d_orb_shell_sizes,
        d_orb_ao_offsets_cart, d_orb_ao_offsets_sph, (int)is_spherical,
        d_aux_norms, d_orb_norms, d_U_aux, d_U_orb, naux_cart, naux, nao_cart,
        nao, d_D3, d_shell_atom_aux, d_shell_atom_orb, workspace_3c.device,
        d_grad);

    QC_RI_Free_Integral_Workspace(&workspace_3c);
    free_uploads();
    return true;
}

// 从 stored eri3c 累加一个自旋通道的 Z_K。各通道求和后只对
// metric inverse square root 做一次 Daleckii-Krein 导数。
static inline void QC_Accumulate_Z_K_Stored_Channel(
    const int nao, const int naux, const double* eri3c,
    const float* spin_density, const float* B_factor,
    const float* density_factor, const int factor_rank,
    std::vector<double>& Z_K)
{
    const long long nao2 = (long long)nao * nao;
    const size_t naux2 = (size_t)naux * naux;
    if (Z_K.size() != naux2) Z_K.assign(naux2, 0.0);

    if (factor_rank > 0 && spin_density != nullptr && B_factor != nullptr &&
        density_factor != nullptr)
    {
        const int M_dim = naux * nao;

        // Z_K[P,Q] = -Σ_{m,l} eri3c[Q,m,l] × V[P,m,l]
        //   V[P,m,l] = Σ_{n,i} B_occ[P,n,i] × C[l,i] × D[m,n]
        // 分解为矩阵乘法链:
        //   R[P,i,m] = Σ_n B_occ[P,n,i] × D[m,n]
        //   V[P,m,l] = Σ_i R[P,i,m] × C[l,i] = Σ_i D[m,n]·B_occ[P,n,i]·C[l,i]
        //   Z_K = -V_flat @ eri3c^T   (将 P,ml 展开)
        //
        // 实际用 V[naux, nao²] 和 eri3c[naux, nao²]:
        //   Z_K[P,Q] = -Σ_{ml} V[P, m*nao+l] × eri3c[Q, m*nao+l]

#if RI_GRAD_HAS_BLAS
        // BLAS 路径: V_P = D @ B_P @ C^T, Z_K = -V @ eri3c^T
        std::vector<double> D_d(nao2);
        for (long long i = 0; i < nao2; i++) D_d[i] = (double)spin_density[i];

        std::vector<double> C_d((size_t)nao * factor_rank);
        for (int lam = 0; lam < nao; lam++)
            for (int i = 0; i < factor_rank; i++)
                C_d[lam * factor_rank + i] =
                    (double)density_factor[lam * nao + i];

        std::vector<double> V(naux * nao2, 0.0);
        std::vector<double> B_P_d((size_t)nao * factor_rank);
        std::vector<double> X_P((size_t)nao * nao);
        for (int P = 0; P < naux; P++)
        {
            for (int i = 0; i < factor_rank; i++)
                for (int n = 0; n < nao; n++)
                    B_P_d[(size_t)n * factor_rank + i] = (double)
                        B_factor[(size_t)(P * nao + n) + (size_t)M_dim * i];

            // X_P[n,l] = Σ_i B_P[n,i] × C[l,i]
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, nao, nao,
                        factor_rank, 1.0, B_P_d.data(), factor_rank,
                        C_d.data(), factor_rank, 0.0,
                        X_P.data(), nao);

            // V_P[m,l] = Σ_n D[m,n] × X_P[n,l]
            double* V_P = V.data() + (long long)P * nao2;
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, nao, nao,
                        nao, 1.0, D_d.data(), nao, X_P.data(), nao, 0.0, V_P,
                        nao);
        }

        // Z_K = -V @ eri3c^T: [naux,nao²] @ [nao²,naux] → [naux,naux]
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, naux, naux,
                    (int)nao2, -1.0, V.data(), (int)nao2, eri3c, (int)nao2, 1.0,
                    Z_K.data(), naux);
#else
        // 标量回退
        std::vector<double> V(naux * nao2, 0.0);
        for (int P = 0; P < naux; P++)
        {
            for (int i = 0; i < factor_rank; i++)
            {
                for (int m = 0; m < nao; m++)
                {
                    double r = 0.0;
                    for (int n = 0; n < nao; n++)
                        r += (double)B_factor[(size_t)(P * nao + n) +
                                              (size_t)M_dim * i] *
                             (double)spin_density[m * nao + n];
                    for (int l = 0; l < nao; l++)
                        V[(long long)P * nao2 + m * nao + l] +=
                            r * (double)density_factor[l * nao + i];
                }
            }
        }

        for (int P = 0; P < naux; P++)
            for (int Q = 0; Q < naux; Q++)
            {
                double z = 0.0;
                for (long long ml = 0; ml < nao2; ml++)
                    z += V[(long long)P * nao2 + ml] *
                         eri3c[(long long)Q * nao2 + ml];
                Z_K[(size_t)P * naux + Q] -= z;
            }
#endif
    }
}
