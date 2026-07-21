#pragma once

// Metric (P|Q) 分解（特征分解）：
// (P|Q)^{-1/2}  → 用于 RI-K B 张量构建
// (P|Q)^{-1}    → 用于 RI-J 拟合系数求解

#include "../../structure/matrix.h"

// Build both truncated metric inverses from one authoritative eigensystem.
// The inverse is formed as inv_sqrt * inv_sqrt, so RI-J, RI-K, and the metric
// response use exactly the same retained subspace.
static bool QC_RI_Build_Metric_Inverses(
    SOLVER_HANDLE solver_handle, BLAS_HANDLE blas_handle, int naux,
    const double* d_metric, double* d_inv_sqrt, double* d_inv,
    int* out_naux_eff, std::vector<double>* out_eigval,
    std::vector<double>* out_eigvec, std::vector<double>* out_inv_sqrt,
    std::vector<double>* out_inv, int* out_solver_api_status,
    int* out_solver_info, double lindep_thresh = 1e-10)
{
    if (out_naux_eff) *out_naux_eff = 0;
    if (out_solver_api_status) *out_solver_api_status = 0;
    if (out_solver_info) *out_solver_info = 0;
    if (naux <= 0 || naux > std::numeric_limits<int>::max() / naux ||
        d_metric == nullptr || d_inv_sqrt == nullptr || d_inv == nullptr ||
        out_naux_eff == nullptr || out_eigval == nullptr ||
        out_eigvec == nullptr || out_inv_sqrt == nullptr ||
        out_inv == nullptr || out_solver_api_status == nullptr ||
        out_solver_info == nullptr ||
        !Double_Memory_Is_Finite(&lindep_thresh) ||
        !(lindep_thresh > 0.0))
        return false;
    const size_t naux2 = (size_t)naux * naux;
    if (naux2 > std::numeric_limits<size_t>::max() / sizeof(double))
        return false;

    deviceMemcpy(d_inv_sqrt, d_metric, sizeof(double) * naux2,
                 deviceMemcpyDeviceToDevice);

    double* d_eigval = NULL;
    if (!Device_Malloc_Safely((void**)&d_eigval, sizeof(double) * naux) ||
        d_eigval == nullptr)
        return false;

    double* d_work = NULL;
    int lwork = 0;
    const int workspace_status = QC_Diagonalize_Double_Workspace_Size(
        solver_handle, naux, d_inv_sqrt, d_eigval, &d_work, &lwork);
    if (workspace_status != 0 || lwork <= 0 || d_work == nullptr)
    {
        *out_solver_api_status =
            workspace_status != 0 ? workspace_status
                                  : (lwork <= 0 ? -2 : -3);
        if (d_work) deviceFree(d_work);
        if (d_eigval) deviceFree(d_eigval);
        return false;
    }

    int info = 0;
    const int api_status = QC_Diagonalize_Double(
        solver_handle, naux, d_inv_sqrt, d_eigval, d_work, lwork, &info);
    *out_solver_api_status = api_status;
    *out_solver_info = info;
    if (api_status != 0 || info != 0)
    {
        if (d_work) deviceFree(d_work);
        if (d_eigval) deviceFree(d_eigval);
        return false;
    }

    std::vector<double> h_eigval(naux);
    deviceMemcpy(h_eigval.data(), d_eigval, sizeof(double) * naux,
                 deviceMemcpyDeviceToHost);

    std::vector<double> h_eigvec(naux2);
    deviceMemcpy(h_eigvec.data(), d_inv_sqrt, sizeof(double) * naux2,
                 deviceMemcpyDeviceToHost);
    for (int i = 0; i < naux; ++i)
    {
        if (!Double_Memory_Is_Finite(&h_eigval[i]) ||
            h_eigval[i] < -lindep_thresh ||
            (i > 0 && h_eigval[i] < h_eigval[i - 1]))
        {
            deviceFree(d_work);
            deviceFree(d_eigval);
            return false;
        }
    }
    for (size_t i = 0; i < h_eigvec.size(); ++i)
    {
        if (!Double_Memory_Is_Finite(&h_eigvec[i]))
        {
            deviceFree(d_work);
            deviceFree(d_eigval);
            return false;
        }
    }

    int n_skip = 0;
    for (int i = 0; i < naux; i++)
    {
        if (h_eigval[i] < lindep_thresh)
            n_skip++;
        else
            break;
    }
    const int naux_eff = naux - n_skip;
    if (naux_eff <= 0)
    {
        deviceFree(d_work);
        deviceFree(d_eigval);
        return false;
    }
    for (int i = n_skip; i < naux; ++i)
    {
        if (!(h_eigval[i] >= lindep_thresh) || !(h_eigval[i] > 0.0))
        {
            deviceFree(d_work);
            deviceFree(d_eigval);
            return false;
        }
    }

    double* d_V = NULL;
    const size_t retained_size = (size_t)naux * naux_eff;
    if (!Device_Malloc_Safely((void**)&d_V, sizeof(double) * retained_size) ||
        d_V == nullptr)
    {
        deviceFree(d_work);
        deviceFree(d_eigval);
        return false;
    }

    std::vector<double> h_V(retained_size);
    for (int k = 0; k < naux_eff; k++)
    {
        const double scale = pow(h_eigval[k + n_skip], -0.25);
        if (!Double_Memory_Is_Finite(&scale) || !(scale > 0.0))
        {
            deviceFree(d_V);
            deviceFree(d_work);
            deviceFree(d_eigval);
            return false;
        }
        for (int i = 0; i < naux; i++)
        {
            h_V[i + k * naux] = h_eigvec[i + (k + n_skip) * naux] * scale;
            if (!Double_Memory_Is_Finite(&h_V[i + k * naux]))
            {
                deviceFree(d_V);
                deviceFree(d_work);
                deviceFree(d_eigval);
                return false;
            }
        }
    }
    deviceMemcpy(d_V, h_V.data(), sizeof(double) * retained_size,
                 deviceMemcpyHostToDevice);

    const double one = 1.0, zero = 0.0;
    deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_T, naux, naux,
                    naux_eff, &one, d_V, naux, d_V, naux, &zero, d_inv_sqrt,
                    naux);
    deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, naux, naux,
                    naux, &one, d_inv_sqrt, naux, d_inv_sqrt, naux, &zero,
                    d_inv, naux);

    std::vector<double> h_inv_sqrt(naux2);
    std::vector<double> h_inv(naux2);
    deviceMemcpy(h_inv_sqrt.data(), d_inv_sqrt, sizeof(double) * naux2,
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(h_inv.data(), d_inv, sizeof(double) * naux2,
                 deviceMemcpyDeviceToHost);
    for (size_t i = 0; i < naux2; ++i)
    {
        if (!Double_Memory_Is_Finite(&h_inv_sqrt[i]) ||
            !Double_Memory_Is_Finite(&h_inv[i]))
        {
            deviceFree(d_V);
            deviceFree(d_work);
            deviceFree(d_eigval);
            return false;
        }
    }
    for (int i = 0; i < naux; ++i)
    {
        const double inv_sqrt_diagonal = h_inv_sqrt[(size_t)i * naux + i];
        const double inv_diagonal = h_inv[(size_t)i * naux + i];
        if (!(inv_sqrt_diagonal > 0.0) || !(inv_diagonal > 0.0))
        {
            deviceFree(d_V);
            deviceFree(d_work);
            deviceFree(d_eigval);
            return false;
        }
    }

    *out_naux_eff = naux_eff;
    *out_eigval = std::move(h_eigval);
    *out_eigvec = std::move(h_eigvec);
    *out_inv_sqrt = std::move(h_inv_sqrt);
    *out_inv = std::move(h_inv);

    if (d_V) deviceFree(d_V);
    if (d_work) deviceFree(d_work);
    if (d_eigval) deviceFree(d_eigval);

    return true;
}
