#pragma once

#include "diis_coefficients.hpp"

// 批量 Tr(A_i · B_j): 在 device 上用 BLAS 计算，避免大量 D2H 拷贝
// h_out[i*mb+j] = Σ_k a_row[i][k] * b_row[j][k]
static void QC_Batched_Trace(BLAS_HANDLE blas_handle, int n, int ma,
                             const double* const* a_ptrs, int mb,
                             const double* const* b_ptrs, double* d_gather_a,
                             double* d_gather_b, double* d_dot_out,
                             std::vector<double>& h_out)
{
    h_out.assign(ma * mb, 0.0);
    if (ma == 0 || mb == 0) return;

    const size_t row_bytes = sizeof(double) * n;

    // D2D gather: 把分散的历史向量收集到连续 device 缓冲
    for (int i = 0; i < ma; i++)
        deviceMemcpy(d_gather_a + (size_t)i * n, a_ptrs[i], row_bytes,
                     deviceMemcpyDeviceToDevice);
    for (int j = 0; j < mb; j++)
        deviceMemcpy(d_gather_b + (size_t)j * n, b_ptrs[j], row_bytes,
                     deviceMemcpyDeviceToDevice);

    // Device BLAS: C[ma × mb] = A[ma × n] * B^T[n × mb]
    QC_Dgemm_NT(blas_handle, ma, mb, n, d_gather_a, n, d_gather_b, n, d_dot_out,
                mb);

    // 只拷回小结果矩阵 (ma × mb ≤ 36 doubles)
    deviceMemcpy(h_out.data(), d_dot_out, sizeof(double) * ma * mb,
                 deviceMemcpyDeviceToHost);
}

// Add alpha and (for UKS) beta traces.  A single coefficient vector must act
// on both spin channels because each Fock matrix depends on both densities.
static void QC_Batched_Trace_Spin_Summed(
    BLAS_HANDLE blas_handle, int n, int ma, const double* const* a_ptrs,
    int mb, const double* const* b_ptrs,
    const double* const* a_ptrs_beta, const double* const* b_ptrs_beta,
    double* d_gather_a, double* d_gather_b, double* d_dot_out,
    std::vector<double>& h_out)
{
    QC_Batched_Trace(blas_handle, n, ma, a_ptrs, mb, b_ptrs, d_gather_a,
                     d_gather_b, d_dot_out, h_out);
    if (a_ptrs_beta == nullptr || b_ptrs_beta == nullptr) return;

    std::vector<double> beta;
    QC_Batched_Trace(blas_handle, n, ma, a_ptrs_beta, mb, b_ptrs_beta,
                     d_gather_a, d_gather_b, d_dot_out, beta);
    for (size_t i = 0; i < h_out.size(); ++i) h_out[i] += beta[i];
}

// DIIS 误差构造
static void QC_Build_DIIS_Error_Double(BLAS_HANDLE blas_handle, int nao,
                                       const double* d_F, const float* d_P,
                                       const float* d_S, double* d_err,
                                       double* d_tmp1, double* d_tmp2,
                                       double* d_tmp3)
{
    const int nao2 = nao * nao;
    QC_Float_To_Double(nao2, d_P, d_tmp1);
    QC_Float_To_Double(nao2, d_S, d_tmp2);
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_F, nao, d_tmp1, nao, d_tmp3, nao);
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_tmp3, nao, d_tmp2, nao, d_err,
                nao);
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_tmp2, nao, d_tmp1, nao, d_tmp3,
                nao);
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_tmp3, nao, d_F, nao, d_tmp1, nao);
    QC_Double_Sub(nao2, d_err, d_tmp1, d_err);
}

static int QC_DIIS_History_Advance(int diis_space, int& hist_count,
                                   int& hist_head)
{
    if (diis_space <= 0) return -1;
    int write_idx = 0;
    if (hist_count < diis_space)
    {
        write_idx = (hist_head + hist_count) % diis_space;
        hist_count++;
    }
    else
    {
        write_idx = hist_head;
        hist_head = (hist_head + 1) % diis_space;
        hist_count = diis_space;
    }
    return write_idx;
}

static void QC_DIIS_History_Store_Channel(
    int nao2, int write_idx, double** d_f_hist, double** d_e_hist,
    double** d_d_hist, const double* d_f_new, const double* d_e_new,
    const float* d_P)
{
    const int bytes_d = sizeof(double) * nao2;
    deviceMemcpy(d_f_hist[write_idx], d_f_new, bytes_d,
                 deviceMemcpyDeviceToDevice);
    deviceMemcpy(d_e_hist[write_idx], d_e_new, bytes_d,
                 deviceMemcpyDeviceToDevice);
    QC_Float_To_Double(nao2, d_P, d_d_hist[write_idx]);
}

// EDIIS 外推
// E^EDIIS(c) = Σ c_i E_i - 0.5 Σ_ij c_i c_j Tr((F_i-F_j)(D_i-D_j))
static bool QC_EDIIS_Extrapolate(BLAS_HANDLE blas_handle, int nao2,
                                 int diis_space, int hist_count, int hist_head,
                                 double** d_f_hist, double** d_d_hist,
                                 double** d_f_hist_beta,
                                 double** d_d_hist_beta,
                                 double* energy_hist, double* d_gather_a,
                                 double* d_gather_b, double* d_dot_out,
                                 std::vector<double>& c_out)
{
    const int m = std::min(hist_count, diis_space);
    if (m < 2) return false;
    auto idx = [&](int i) { return (hist_head + i) % diis_space; };

    // Tr(F_i * D_j) via dgemm
    std::vector<const double*> fa(m), da(m), fb, db;
    if (d_f_hist_beta != nullptr && d_d_hist_beta != nullptr)
    {
        fb.resize(m);
        db.resize(m);
    }
    for (int i = 0; i < m; i++)
    {
        fa[i] = d_f_hist[idx(i)];
        da[i] = d_d_hist[idx(i)];
        if (!fb.empty())
        {
            fb[i] = d_f_hist_beta[idx(i)];
            db[i] = d_d_hist_beta[idx(i)];
        }
    }
    std::vector<double> FD;
    QC_Batched_Trace_Spin_Summed(
        blas_handle, nao2, m, fa.data(), m, da.data(),
        fb.empty() ? nullptr : fb.data(), db.empty() ? nullptr : db.data(),
        d_gather_a, d_gather_b, d_dot_out, FD);

    // H_ij = Tr((F_i-F_j)(D_i-D_j)) = FD[i,i] + FD[j,j] - FD[i,j] - FD[j,i]
    // 目标: min Σ c_i E_i - 0.5 Σ_ij c_i c_j H_ij = min g^T c + 0.5 c^T (-H) c
    std::vector<double> H(m * m), g(m);
    for (int i = 0; i < m; i++)
    {
        g[i] = energy_hist[idx(i)];
        for (int j = 0; j < m; j++)
        {
            double hij =
                FD[i * m + i] + FD[j * m + j] - FD[i * m + j] - FD[j * m + i];
            H[i * m + j] = -hij;  // QP 的 H 矩阵取负号
        }
    }

    return QC_SCF_Solve_Simplex_QP(m, H, g, c_out);
}

// ADIIS 外推
// E^ADIIS(c) = E_n + 2 Σ c_i Tr((D_i-D_n) F_n)
//            + Σ_ij c_i c_j Tr((D_i-D_n)(F_j-F_n))
static bool QC_ADIIS_Extrapolate(BLAS_HANDLE blas_handle, int nao2,
                                 int diis_space, int hist_count, int hist_head,
                                 double** d_f_hist, double** d_d_hist,
                                 double** d_f_hist_beta,
                                 double** d_d_hist_beta,
                                 double* d_gather_a, double* d_gather_b,
                                 double* d_dot_out,
                                 std::vector<double>& c_out)
{
    const int m = std::min(hist_count, diis_space);
    if (m < 2) return false;
    auto idx = [&](int i) { return (hist_head + i) % diis_space; };

    // Tr(D_i * F_j) via dgemm
    std::vector<const double*> da(m), fa(m), db, fb;
    if (d_f_hist_beta != nullptr && d_d_hist_beta != nullptr)
    {
        db.resize(m);
        fb.resize(m);
    }
    for (int i = 0; i < m; i++)
    {
        da[i] = d_d_hist[idx(i)];
        fa[i] = d_f_hist[idx(i)];
        if (!db.empty())
        {
            db[i] = d_d_hist_beta[idx(i)];
            fb[i] = d_f_hist_beta[idx(i)];
        }
    }
    std::vector<double> DF;
    QC_Batched_Trace_Spin_Summed(
        blas_handle, nao2, m, da.data(), m, fa.data(),
        db.empty() ? nullptr : db.data(), fb.empty() ? nullptr : fb.data(),
        d_gather_a, d_gather_b, d_dot_out, DF);

    std::vector<double> H, g;
    if (!QC_SCF_Build_ADIIS_QP(DF, m, H, g)) return false;
    return QC_SCF_Solve_Simplex_QP(m, H, g, c_out);
}

// CDIIS 外推
static bool QC_CDIIS_Extrapolate(BLAS_HANDLE blas_handle, int nao,
                                 int diis_space, int hist_count, int hist_head,
                                 double** d_f_hist, double** d_e_hist,
                                 double** d_f_hist_beta,
                                 double** d_e_hist_beta, double reg,
                                 double* d_f_out, double* d_f_out_beta,
                                 double* d_gather_a, double* d_gather_b,
                                 double* d_dot_out)
{
    if (hist_count < 2 || diis_space <= 0) return false;
    const int m = std::min(hist_count, diis_space);
    if (m < 2) return false;
    const int nao2 = nao * nao;
    auto hist_idx = [&](int logical_idx) -> int
    { return (hist_head + logical_idx) % diis_space; };

    // The UKS Pulay metric is the direct-sum inner product of both spin
    // commutators.  Solving this one system and applying its coefficients to
    // both Fock histories preserves the coupled alpha/beta SCF map.
    std::vector<const double*> ea(m), eb;
    if (d_e_hist_beta != nullptr) eb.resize(m);
    for (int i = 0; i < m; i++)
    {
        ea[i] = d_e_hist[hist_idx(i)];
        if (!eb.empty()) eb[i] = d_e_hist_beta[hist_idx(i)];
    }
    std::vector<double> EET;
    QC_Batched_Trace_Spin_Summed(
        blas_handle, nao2, m, ea.data(), m, ea.data(),
        eb.empty() ? nullptr : eb.data(), eb.empty() ? nullptr : eb.data(),
        d_gather_a, d_gather_b, d_dot_out, EET);
    std::vector<double> coefficients;
    if (!QC_SCF_Solve_CDIIS_Coefficients(EET, m, reg, coefficients))
        return false;

    deviceMemset(d_f_out, 0, sizeof(double) * nao2);
    if (d_f_out_beta != nullptr)
        deviceMemset(d_f_out_beta, 0, sizeof(double) * nao2);
    for (int i = 0; i < m; i++)
    {
        const double c = coefficients[i];
        QC_Double_Axpy(nao2, c, d_f_hist[hist_idx(i)], d_f_out);
        if (d_f_out_beta != nullptr && d_f_hist_beta != nullptr)
            QC_Double_Axpy(nao2, c, d_f_hist_beta[hist_idx(i)],
                           d_f_out_beta);
    }
    return true;
}

// MESA + CDIIS 外推
// MESA: 同时算 EDIIS 和 ADIIS，选密度变化更小的
// 切换: 误差范数大时用 MESA，小时用 CDIIS
static bool QC_MESA_Or_CDIIS_Extrapolate(
    BLAS_HANDLE blas_handle, int nao, int diis_space, int hist_count,
    int hist_head, double** d_f_hist, double** d_e_hist, double** d_d_hist,
    double** d_f_hist_beta, double** d_e_hist_beta,
    double** d_d_hist_beta, double* energy_hist, double reg, double enorm,
    double mesa_to_cdiis_threshold, double* d_f_out, double* d_f_out_beta,
    double* d_gather_a, double* d_gather_b, double* d_dot_out)
{
    const int m = std::min(hist_count, diis_space);
    if (m < 2) return false;
    const int nao2 = nao * nao;
    auto idx = [&](int i) { return (hist_head + i) % diis_space; };

    const QC_SCF_DIIS_Extrapolation_Plan plan =
        QC_SCF_Select_DIIS_Extrapolation_Plan(
            enorm, mesa_to_cdiis_threshold);
    if (plan.try_cdiis_first)
    {
        if (QC_CDIIS_Extrapolate(
                blas_handle, nao, diis_space, hist_count, hist_head,
                d_f_hist, d_e_hist, d_f_hist_beta, d_e_hist_beta, reg,
                d_f_out, d_f_out_beta, d_gather_a, d_gather_b, d_dot_out))
            return true;
    }
    if (!plan.try_bounded_simplex) return false;

    std::vector<double> c_ediis, c_adiis;
    bool ok_e = QC_EDIIS_Extrapolate(
        blas_handle, nao2, diis_space, hist_count, hist_head, d_f_hist,
        d_d_hist, d_f_hist_beta, d_d_hist_beta, energy_hist, d_gather_a,
        d_gather_b, d_dot_out, c_ediis);
    bool ok_a = QC_ADIIS_Extrapolate(
        blas_handle, nao2, diis_space, hist_count, hist_head, d_f_hist,
        d_d_hist, d_f_hist_beta, d_d_hist_beta, d_gather_a, d_gather_b,
        d_dot_out, c_adiis);
    if (!ok_e && !ok_a) return false;

    // Compare the actual spin-summed Frobenius density displacement, not a
    // proxy based only on the newest coefficient.
    std::vector<const double*> da(m), db;
    if (d_d_hist_beta != nullptr) db.resize(m);
    for (int i = 0; i < m; ++i)
    {
        da[i] = d_d_hist[idx(i)];
        if (!db.empty()) db[i] = d_d_hist_beta[idx(i)];
    }
    std::vector<double> density_gram;
    QC_Batched_Trace_Spin_Summed(
        blas_handle, nao2, m, da.data(), m, da.data(),
        db.empty() ? nullptr : db.data(), db.empty() ? nullptr : db.data(),
        d_gather_a, d_gather_b, d_dot_out, density_gram);
    auto density_displacement_squared = [&](const std::vector<double>& c)
    {
        double value = 0.0;
        for (int i = 0; i < m; ++i)
        {
            const double qi = c[i] - (i == m - 1 ? 1.0 : 0.0);
            for (int j = 0; j < m; ++j)
            {
                const double qj = c[j] - (j == m - 1 ? 1.0 : 0.0);
                value += qi * qj * density_gram[(size_t)i * m + j];
            }
        }
        if (!QC_SCF_DIIS_Double_Is_Finite(value))
            return std::numeric_limits<double>::infinity();
        return std::max(0.0, value);
    };

    const std::vector<double>& c_best = [&]() -> const std::vector<double>&
    {
        if (!ok_e) return c_adiis;
        if (!ok_a) return c_ediis;
        return density_displacement_squared(c_ediis) <=
                       density_displacement_squared(c_adiis)
                   ? c_ediis
                   : c_adiis;
    }();

    deviceMemset(d_f_out, 0, sizeof(double) * nao2);
    if (d_f_out_beta != nullptr)
        deviceMemset(d_f_out_beta, 0, sizeof(double) * nao2);
    for (int i = 0; i < m; i++)
    {
        QC_Double_Axpy(nao2, c_best[i], d_f_hist[idx(i)], d_f_out);
        if (d_f_out_beta != nullptr && d_f_hist_beta != nullptr)
            QC_Double_Axpy(nao2, c_best[i], d_f_hist_beta[idx(i)],
                           d_f_out_beta);
    }
    return true;
}

// SCF 中应用 DIIS
// MESA 算法: 远离收敛用 EDIIS/ADIIS，近收敛用 CDIIS
// 参考:
//   S. Lehtola, "OpenOrbitalOptimizer", arXiv:2503.23034 (2025).
//   X. Hu, W. Yang, J. Chem. Phys. 132, 054109 (2010). (ADIIS)
//   K. N. Kudin et al., J. Chem. Phys. 116, 8255 (2002). (EDIIS)
void QUANTUM_CHEMISTRY::Apply_DIIS(int iter)
{
    if (!scf_ws.runtime.use_diis || (iter + 1) < scf_ws.runtime.diis_start_iter)
        return;

    const int nao2 = (int)mol.nao2;
    const bool unrestricted = scf_ws.runtime.unrestricted;
    double* dF = scf_ws.alpha.d_F_double;
    double* dFb = unrestricted ? scf_ws.beta.d_F_double : nullptr;

    double current_energy = 0.0;
    deviceMemcpy(&current_energy, scf_ws.core.d_scf_energy, sizeof(double),
                 deviceMemcpyDeviceToHost);

    // Alpha and beta share one ring position.  This makes every history item
    // a complete UKS state rather than two independently indexed channels.
    const int write_idx = QC_DIIS_History_Advance(
        scf_ws.runtime.diis_space, scf_ws.diis.diis_hist_count,
        scf_ws.diis.diis_hist_head);
    if (write_idx < 0) return;

    deviceMemset(scf_ws.diis.d_diis_accum, 0, sizeof(double));
    QC_Build_DIIS_Error_Double(
        blas_handle, mol.nao, dF, scf_ws.alpha.d_P, scf_ws.core.d_S,
        scf_ws.diis.d_diis_err, scf_ws.ortho.d_dwork_nao2_2,
        scf_ws.ortho.d_dwork_nao2_3, scf_ws.ortho.d_dwork_nao2_4);
    QC_Double_Dot(nao2, scf_ws.diis.d_diis_err, scf_ws.diis.d_diis_err,
                  scf_ws.diis.d_diis_accum);
    QC_DIIS_History_Store_Channel(
        nao2, write_idx, scf_ws.diis.d_diis_f_hist.data(),
        scf_ws.diis.d_diis_e_hist.data(),
        scf_ws.diis.d_diis_d_hist.data(), dF, scf_ws.diis.d_diis_err,
        scf_ws.alpha.d_P);

    if (unrestricted)
    {
        QC_Build_DIIS_Error_Double(
            blas_handle, mol.nao, dFb, scf_ws.beta.d_P, scf_ws.core.d_S,
            scf_ws.diis.d_diis_err, scf_ws.ortho.d_dwork_nao2_2,
            scf_ws.ortho.d_dwork_nao2_3, scf_ws.ortho.d_dwork_nao2_4);
        QC_Double_Dot(nao2, scf_ws.diis.d_diis_err,
                      scf_ws.diis.d_diis_err, scf_ws.diis.d_diis_accum);
        QC_DIIS_History_Store_Channel(
            nao2, write_idx, scf_ws.diis.d_diis_f_hist_b.data(),
            scf_ws.diis.d_diis_e_hist_b.data(),
            scf_ws.diis.d_diis_d_hist_b.data(), dFb,
            scf_ws.diis.d_diis_err, scf_ws.beta.d_P);
    }
    scf_ws.diis.energy_hist[write_idx] = current_energy;

    double enorm_sq = 0.0;
    deviceMemcpy(&enorm_sq, scf_ws.diis.d_diis_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    const double enorm =
        std::sqrt(enorm_sq / (nao2 * (unrestricted ? 2.0 : 1.0)));
    if (!QC_SCF_DIIS_Double_Is_Finite(enorm))
    {
        scf_ws.diis.diis_hist_count = scf_ws.diis.diis_hist_head = 0;
        scf_ws.diis.last_enorm = std::numeric_limits<double>::max();
        return;
    }
    scf_ws.diis.last_enorm = enorm;

    bool extrapolated = false;
    while (scf_ws.diis.diis_hist_count >= 2)
    {
        extrapolated = QC_MESA_Or_CDIIS_Extrapolate(
            blas_handle, mol.nao, scf_ws.runtime.diis_space,
            scf_ws.diis.diis_hist_count, scf_ws.diis.diis_hist_head,
            scf_ws.diis.d_diis_f_hist.data(),
            scf_ws.diis.d_diis_e_hist.data(),
            scf_ws.diis.d_diis_d_hist.data(),
            unrestricted ? scf_ws.diis.d_diis_f_hist_b.data() : nullptr,
            unrestricted ? scf_ws.diis.d_diis_e_hist_b.data() : nullptr,
            unrestricted ? scf_ws.diis.d_diis_d_hist_b.data() : nullptr,
            scf_ws.diis.energy_hist.data(), scf_ws.runtime.diis_reg, enorm,
            scf_ws.diis.mesa_to_cdiis_threshold, dF, dFb,
            scf_ws.diis.d_gather_a, scf_ws.diis.d_gather_b,
            scf_ws.diis.d_dot_out);
        if (extrapolated) break;

        // A failed or unsafe solve is repaired by shrinking the subspace,
        // never by publishing non-finite or strongly cancelling Fock data.
        scf_ws.diis.diis_hist_head =
            (scf_ws.diis.diis_hist_head + 1) % scf_ws.runtime.diis_space;
        --scf_ws.diis.diis_hist_count;
    }

    if (!extrapolated) return;
    QC_Double_To_Float(nao2, dF, scf_ws.alpha.d_F);
    if (unrestricted) QC_Double_To_Float(nao2, dFb, scf_ws.beta.d_F);
}
