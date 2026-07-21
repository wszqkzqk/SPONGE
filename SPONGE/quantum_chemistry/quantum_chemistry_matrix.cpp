#include "quantum_chemistry.h"
#include "structure/cart2sph.hpp"
#include "structure/matrix.h"

#include <cfloat>

// 单精度 BLAS 包装

void QC_MatMul_RowRow_Blas(BLAS_HANDLE blas_handle, int m, int n, int kdim,
                           const float* A_row, const float* B_row, float* C_row)
{
    const float alpha = 1.0f;
    const float beta = 0.0f;
    deviceBlasSgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, n, m, kdim,
                    &alpha, B_row, n, A_row, kdim, &beta, C_row, n);
}

void QC_MatMul_RowCol_Blas(BLAS_HANDLE blas_handle, int m, int n, int kdim,
                           const float* A_row, const float* B_col, float* C_row)
{
    const float alpha = 1.0f;
    const float beta = 0.0f;
    deviceBlasSgemm(blas_handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_N, n, m, kdim,
                    &alpha, B_col, kdim, A_row, kdim, &beta, C_row, n);
}

void QC_Build_Density_Blas(BLAS_HANDLE blas_handle, int nao, int n_occ,
                           float density_factor, const float* C_row,
                           float* P_new_row)
{
    const int nao2 = (int)nao * (int)nao;
    deviceMemset(P_new_row, 0, sizeof(float) * nao2);
    if (n_occ <= 0 || density_factor == 0.0f) return;

    const float alpha = density_factor;
    const float beta = 0.0f;
    deviceBlasSgemm(blas_handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_N, nao, nao,
                    n_occ, &alpha, C_row, nao, C_row, nao, &beta, P_new_row,
                    nao);
}

void QC_Build_Density_Double_Blas(BLAS_HANDLE blas_handle, int nao, int n_occ,
                                  double density_factor, const double* C_row,
                                  int orbital_stride, double* P_new_row)
{
    const int nao2 = nao * nao;
    deviceMemset(P_new_row, 0, sizeof(double) * nao2);
    if (n_occ <= 0 || density_factor == 0.0) return;

    const double beta = 0.0;
    // C_row is nao x orbital_stride in row-major storage.  BLAS sees its
    // transpose, so op(T)*op(N) forms C_occ*C_occ^T in the same row-major AO
    // layout used throughout SCF.  Occupation is applied before the sole
    // double-to-float conversion of the resulting density.
    deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_N, nao, nao,
                    n_occ, &density_factor, C_row, orbital_stride, C_row,
                    orbital_stride, &beta, P_new_row, nao);
}

// 双精度 BLAS/Solver 包装

int QC_Diagonalize_Double_Workspace_Size(SOLVER_HANDLE solver_handle, int n,
                                         double* mat, double* w,
                                         double** work_ptr, int* lwork)
{
    if (work_ptr == NULL || lwork == NULL || mat == NULL || w == NULL || n <= 0)
    {
        if (lwork != NULL) *lwork = 0;
        return -1;
    }
#ifdef USE_GPU
    int stat = (int)deviceSolverDsyevdBufferSize(
        solver_handle, DEVICE_EIG_MODE_VECTOR, DEVICE_FILL_MODE_UPPER, n, mat,
        n, w, lwork);
    if (stat != SOLVER_SUCCESS || *lwork <= 0)
    {
        *lwork = 0;
        return (stat != SOLVER_SUCCESS) ? stat : -2;
    }
#elif defined(USE_MKL) || defined(USE_OPENBLAS)
    (void)solver_handle;
    double work_query = 0.0;
    lapack_int iwork_query = 0;
    const lapack_int query_info = LAPACKE_dsyevd_work(
        LAPACK_COL_MAJOR, 'V', 'U', (lapack_int)n, mat, (lapack_int)n, w,
        &work_query, (lapack_int)-1, &iwork_query, (lapack_int)-1);
    if (query_info != 0 || !Double_Memory_Is_Finite(&work_query) ||
        !(work_query >= 1.0) ||
        work_query > (double)std::numeric_limits<int>::max())
    {
        *lwork = 0;
        return (query_info != 0) ? (int)query_info : -2;
    }
    *lwork = (int)(work_query + 0.5);
#else
    (void)solver_handle;
    *lwork = 0;
    return -1;
#endif
    if (*work_ptr)
    {
        deviceFree(*work_ptr);
        *work_ptr = NULL;
    }
    Device_Malloc_Safely((void**)work_ptr, sizeof(double) * (*lwork));
    if (*work_ptr == NULL)
    {
        *lwork = 0;
        return -3;
    }
    deviceMemset(*work_ptr, 0, sizeof(double) * (*lwork));
    return 0;
}

int QC_Diagonalize_Double(SOLVER_HANDLE solver_handle, int n, double* mat,
                          double* w, double* work, int lwork, int* info)
{
#ifdef USE_GPU
    if (info == NULL || mat == NULL || w == NULL || work == NULL || n <= 0 ||
        lwork <= 0)
        return -1;
    // 1 int 的设备缓冲，直接 alloc/free 避免 function-local static 生命周期
    int* d_info = NULL;
    Device_Malloc_Safely((void**)&d_info, sizeof(int));
    if (d_info == NULL) return -2;
    const int pending_info = std::numeric_limits<int>::min();
    deviceMemcpy(d_info, &pending_info, sizeof(int), deviceMemcpyHostToDevice);
    const int api_status = (int)deviceSolverDsyevd(
        solver_handle, DEVICE_EIG_MODE_VECTOR, DEVICE_FILL_MODE_UPPER, n, mat,
        n, w, work, lwork, d_info);
    if (api_status != SOLVER_SUCCESS)
    {
        *info = pending_info;
        deviceFree(d_info);
        return api_status;
    }
    deviceMemcpy(info, d_info, sizeof(int), deviceMemcpyDeviceToHost);
    deviceFree(d_info);
    return api_status;
#elif defined(USE_MKL) || defined(USE_OPENBLAS)
    (void)solver_handle;
    if (info == NULL || mat == NULL || w == NULL || work == NULL || n <= 0 ||
        lwork <= 0)
        return -1;
    double work_query = 0.0;
    lapack_int iwork_query = 0;
    const lapack_int query_info = LAPACKE_dsyevd_work(
        LAPACK_COL_MAJOR, 'V', 'U', (lapack_int)n, mat, (lapack_int)n, w,
        &work_query, (lapack_int)-1, &iwork_query, (lapack_int)-1);
    if (query_info != 0 || iwork_query <= 0)
    {
        *info = (query_info != 0) ? (int)query_info : -1;
        return 0;
    }
    std::vector<lapack_int> iwork((size_t)iwork_query);
    *info = (int)LAPACKE_dsyevd_work(
        LAPACK_COL_MAJOR, 'V', 'U', (lapack_int)n, mat, (lapack_int)n, w, work,
        (lapack_int)lwork, iwork.data(), iwork_query);
    return 0;
#else
    (void)solver_handle;
    (void)n;
    (void)mat;
    (void)w;
    (void)work;
    (void)lwork;
    if (info != NULL) *info = -1;
    return -1;
#endif
}

void QC_Dgemm_NN(BLAS_HANDLE handle, int m, int n, int k, const double* A,
                 int lda, const double* B, int ldb, double* C, int ldc)
{
    const double one = 1.0, zero = 0.0;
    deviceBlasDgemm(handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, n, m, k, &one,
                    B, ldb, A, lda, &zero, C, ldc);
}

void QC_Dgemm_TN(BLAS_HANDLE handle, int m, int n, int k, const double* A,
                 int lda, const double* B, int ldb, double* C, int ldc)
{
    const double one = 1.0, zero = 0.0;
    deviceBlasDgemm(handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_T, n, m, k, &one,
                    B, ldb, A, lda, &zero, C, ldc);
}

void QC_Dgemm_NT(BLAS_HANDLE handle, int m, int n, int k, const double* A,
                 int lda, const double* B, int ldb, double* C, int ldc)
{
    const double one = 1.0, zero = 0.0;
    deviceBlasDgemm(handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_N, n, m, k, &one,
                    B, ldb, A, lda, &zero, C, ldc);
}

// RI BLAS 包装

void QC_Sgemm_NN(BLAS_HANDLE handle, int m, int n, int k, float alpha,
                 const float* A, int lda, const float* B, int ldb, float beta,
                 float* C, int ldc)
{
    // C[m×n] = alpha * A[m×k] * B[k×n] + beta * C (col-major)
    deviceBlasSgemm(handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, m, n, k, &alpha,
                    A, lda, B, ldb, &beta, C, ldc);
}

void QC_Sgemm_TN(BLAS_HANDLE handle, int m, int n, int k, float alpha,
                 const float* A, int lda, const float* B, int ldb, float beta,
                 float* C, int ldc)
{
    // C[m×n] = alpha * A^T[m×k] * B[k×n] + beta * C (col-major)
    deviceBlasSgemm(handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_N, m, n, k, &alpha,
                    A, lda, B, ldb, &beta, C, ldc);
}

void QC_Sgemm_RowMajor_NN(BLAS_HANDLE handle, int m, int n, int k, float alpha,
                          const float* A, int lda, const float* B, int ldb,
                          float beta, float* C, int ldc)
{
    // Row-major C[m×n] = alpha * A[m×k] * B[k×n] + beta * C[m×n].
    deviceBlasSgemm(handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, n, m, k, &alpha,
                    B, ldb, A, lda, &beta, C, ldc);
}

void QC_Sgemm_RowMajor_NT(BLAS_HANDLE handle, int m, int n, int k, float alpha,
                          const float* A, int lda, const float* B, int ldb,
                          float beta, float* C, int ldc)
{
    // Row-major C[m×n] = alpha * A[m×k] * B^T[k×n] + beta * C[m×n].
    deviceBlasSgemm(handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_N, n, m, k, &alpha,
                    B, ldb, A, lda, &beta, C, ldc);
}

// 常用通用矩阵函数包装

static __global__ void QC_Add_Matrix_Kernel(const int n, const float* A,
                                            const float* B, float* C)
{
    SIMPLE_DEVICE_FOR(idx, n) { C[idx] = A[idx] + B[idx]; }
}

void QC_Add_Matrix(int n, const float* A, const float* B, float* C)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Add_Matrix_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, A, B, C);
}

static __global__ void QC_Sub_Matrix_Kernel(const int n, const float* A,
                                            const float* B, float* C)
{
    SIMPLE_DEVICE_FOR(idx, n) { C[idx] = A[idx] - B[idx]; }
}

void QC_Sub_Matrix(int n, const float* A, const float* B, float* C)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Sub_Matrix_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, A, B, C);
}

static __global__ void QC_Scale_Matrix_By_Norms_Kernel(const int nao,
                                                       const float* norms,
                                                       float* M)
{
    const int total = nao * nao;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        int i = idx / nao;
        int j = idx - i * nao;
        M[idx] *= norms[i] * norms[j];
    }
}

void QC_Scale_Matrix_By_Norms(int nao, const float* norms, float* M)
{
    const int threads = 256;
    const int total = nao * nao;
    Launch_Device_Kernel(QC_Scale_Matrix_By_Norms_Kernel,
                         (total + threads - 1) / threads, threads, 0, 0, nao,
                         norms, M);
}

static __global__ void QC_Float_To_Double_Kernel(const int n, const float* src,
                                                 double* dst)
{
    SIMPLE_DEVICE_FOR(i, n) { dst[i] = (double)src[i]; }
}

void QC_Float_To_Double(int n, const float* src, double* dst)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Float_To_Double_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, src, dst);
}

static __global__ void QC_Double_To_Float_Kernel(const int n, const double* src,
                                                 float* dst)
{
    SIMPLE_DEVICE_FOR(i, n) { dst[i] = (float)src[i]; }
}

void QC_Double_To_Float(int n, const double* src, float* dst)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Double_To_Float_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, src, dst);
}

static __global__ void QC_Float_To_Double_Copy_Kernel(const int n,
                                                      const float* src,
                                                      double* dst)
{
    SIMPLE_DEVICE_FOR(i, n) { dst[i] = (double)src[i]; }
}

void QC_Float_To_Double_Copy(int n, const float* src, double* dst)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Float_To_Double_Copy_Kernel,
                         (n + threads - 1) / threads, threads, 0, 0, n, src,
                         dst);
}

static __global__ void QC_Symmetrize_Double_Matrix_Kernel(const int n,
                                                           double* matrix)
{
    SIMPLE_DEVICE_FOR(index, n * n)
    {
        const int row = index / n;
        const int column = index % n;
        if (row < column)
        {
            const double value =
                0.5 * (matrix[row * n + column] +
                       matrix[column * n + row]);
            matrix[row * n + column] = value;
            matrix[column * n + row] = value;
        }
    }
}

void QC_Symmetrize_Double_Matrix(int n, double* matrix)
{
    const int threads = 256;
    Launch_Device_Kernel(
        QC_Symmetrize_Double_Matrix_Kernel,
        Positive_Int_Ceil_Div(n * n, threads), threads, 0, 0, n, matrix);
}

static __global__ void
QC_Round_Symmetric_Double_Matrix_For_Nonincreasing_Linear_Objective_Kernel(
    const int n, const double* matrix, const double* coefficients,
    float* rounded_matrix)
{
    SIMPLE_DEVICE_FOR(index, n * n)
    {
        const int row = index / n;
        const int column = index % n;
        if (row <= column)
        {
            const int transpose = column * n + row;
            const double value =
                row == column
                    ? matrix[index]
                    : 0.5 * (matrix[index] + matrix[transpose]);
            const double coefficient =
                row == column
                    ? coefficients[index]
                    : coefficients[index] + coefficients[transpose];
            float rounded = static_cast<float>(value);
            if (coefficient * (static_cast<double>(rounded) - value) > 0.0)
                rounded = nextafterf(
                    rounded, coefficient > 0.0 ? -FLT_MAX : FLT_MAX);
            rounded_matrix[index] = rounded;
            rounded_matrix[transpose] = rounded;
        }
    }
}

void QC_Round_Symmetric_Double_Matrix_For_Nonincreasing_Linear_Objective(
    int n, const double* matrix, const double* coefficients,
    float* rounded_matrix)
{
    const int threads = 256;
    Launch_Device_Kernel(
        QC_Round_Symmetric_Double_Matrix_For_Nonincreasing_Linear_Objective_Kernel,
        Positive_Int_Ceil_Div(n * n, threads), threads, 0, 0, n, matrix,
        coefficients, rounded_matrix);
}

static __global__ void QC_Double_Dot_Kernel(const int n, const double* A,
                                            const double* B, double* out_sum)
{
    SIMPLE_DEVICE_FOR(i, n) { atomicAdd(out_sum, A[i] * B[i]); }
}

void QC_Double_Dot(int n, const double* A, const double* B, double* out_sum)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Double_Dot_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, A, B, out_sum);
}

static __global__ void QC_Double_Axpy_Kernel(const int n, const double coeff,
                                             const double* src, double* dst)
{
    SIMPLE_DEVICE_FOR(i, n) { dst[i] += coeff * src[i]; }
}

void QC_Double_Axpy(int n, double coeff, const double* src, double* dst)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Double_Axpy_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, coeff, src, dst);
}

static __global__ void QC_Double_Sub_Kernel(const int n, const double* A,
                                            const double* B, double* dst)
{
    SIMPLE_DEVICE_FOR(i, n) { dst[i] = A[i] - B[i]; }
}

void QC_Double_Sub(int n, const double* A, const double* B, double* dst)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Double_Sub_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, A, B, dst);
}

// 常用 SCF 矩阵函数包装

static __global__ void QC_Elec_Energy_Accumulate_Kernel(const int nao2,
                                                        const float* P,
                                                        const float* H_core,
                                                        const double* F,
                                                        double* out_sum)
{
    SIMPLE_DEVICE_FOR(idx, nao2)
    {
        atomicAdd(out_sum, 0.5 * (double)P[idx] *
                               ((double)H_core[idx] + (double)F[idx]));
    }
}

void QC_Elec_Energy_Accumulate(int nao2, const float* P, const float* H_core,
                               const double* F, double* out_sum)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Elec_Energy_Accumulate_Kernel,
                         (nao2 + threads - 1) / threads, threads, 0, 0, nao2, P,
                         H_core, F, out_sum);
}

static __global__ void QC_Mat_Dot_Accumulate_Kernel(const int nao2,
                                                    const float* A,
                                                    const float* B,
                                                    double* out_sum)
{
    SIMPLE_DEVICE_FOR(idx, nao2)
    {
        atomicAdd(out_sum, (double)A[idx] * (double)B[idx]);
    }
}

void QC_Mat_Dot_Accumulate(int nao2, const float* A, const float* B,
                           double* out_sum)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Mat_Dot_Accumulate_Kernel,
                         (nao2 + threads - 1) / threads, threads, 0, 0, nao2, A,
                         B, out_sum);
}

static __global__ void QC_Level_Shift_Kernel(const int n, const double ls,
                                             const double density_factor,
                                             const double* dS,
                                             const double* dSPS, double* dF)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        dF[i] += ls * (dS[i] - density_factor * dSPS[i]);
    }
}

void QC_Level_Shift(int n, double ls, double density_factor, const double* dS,
                    const double* dSPS, double* dF)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Level_Shift_Kernel, (n + threads - 1) / threads,
                         threads, 0, 0, n, ls, density_factor, dS, dSPS, dF);
}

static __global__ void QC_Build_X_Canonical_Kernel(
    const int nao, const int nao_eff, const double* eigvec_col,
    const double* eigval, const double lindep_thresh, double* X_row)
{
    SIMPLE_DEVICE_FOR(i, nao)
    {
        int col = 0;
        for (int k = 0; k < nao; k++)
        {
            if (eigval[k] < lindep_thresh) continue;
            X_row[i * nao + col] = eigvec_col[i + k * nao] / sqrt(eigval[k]);
            col++;
        }
    }
}

void QC_Build_X_Canonical(int nao, int nao_eff, const double* eigvec_col,
                          const double* eigval, double lindep_thresh,
                          double* X_row)
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Build_X_Canonical_Kernel,
                         (nao + threads - 1) / threads, threads, 0, 0, nao,
                         nao_eff, eigvec_col, eigval, lindep_thresh, X_row);
}

static __global__ void QC_Rect_Double_To_Padded_Float_Kernel(const int nao,
                                                             const int ne,
                                                             const double* src,
                                                             float* dst)
{
    SIMPLE_DEVICE_FOR(idx, nao * nao)
    {
        int i = idx / nao;
        int j = idx % nao;
        dst[idx] = (j < ne) ? (float)src[i * ne + j] : 0.0f;
    }
}

void QC_Rect_Double_To_Padded_Float(int nao, int ne, const double* src,
                                    float* dst)
{
    const int nao2 = nao * nao;
    const int threads = 256;
    Launch_Device_Kernel(QC_Rect_Double_To_Padded_Float_Kernel,
                         (nao2 + threads - 1) / threads, threads, 0, 0, nao, ne,
                         src, dst);
}

// 笛卡尔基转归一化实球谐基。所有角动量都由唯一的通用构造器生成。
std::vector<float> QC_Build_Cart2Sph_Mat_Host(const std::vector<int>& l_list,
                                              int nao_cart, int nao_sph)
{
    return qc_cart2sph::Build_Matrix(l_list, nao_cart, nao_sph);
}

static __global__ void QC_Cart2Sph_MatMul_UT_RowRow_Kernel(
    const int m, const int n, const int kdim, const float* U_row_k_m,
    const float* B_row_k_n, float* C_row_m_n)
{
    SIMPLE_DEVICE_FOR(idx, m * n)
    {
        int i = idx / n;
        int j = idx - i * n;
        double sum = 0.0;
        for (int k = 0; k < kdim; k++)
        {
            sum += (double)U_row_k_m[k * m + i] * (double)B_row_k_n[k * n + j];
        }
        C_row_m_n[i * n + j] = (float)sum;
    }
}

void QUANTUM_CHEMISTRY::Build_Cart2Sph_Matrix()
{
    int nao_c = mol.nao_cart;
    int nao_s = mol.nao_sph;
    std::vector<float> cart2sph_mat;
    try
    {
        cart2sph_mat = QC_Build_Cart2Sph_Mat_Host(mol.h_l_list, nao_c, nao_s);
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Initial",
            "Reason:\n    Failed to construct Cartesian-to-spherical basis "
            "matrix: %s\n",
            error.what());
    }
    if (cart2sph_mat.empty())
    {
        cart2sph.d_cart2sph_mat = nullptr;
    }
    else
    {
        Device_Malloc_Safely((void**)&cart2sph.d_cart2sph_mat,
                             sizeof(float) * cart2sph_mat.size());
        deviceMemcpy(cart2sph.d_cart2sph_mat, cart2sph_mat.data(),
                     sizeof(float) * cart2sph_mat.size(),
                     deviceMemcpyHostToDevice);
    }
}

void QUANTUM_CHEMISTRY::Cart2Sph_Single_Matrix(float* d_cart, float* d_sph)
{
    if (!mol.is_spherical) return;
    const int nao_c = mol.nao_cart;
    const int nao_s = mol.nao_sph;
    const int threads = 256;
    const int total = nao_s * nao_s;
    QC_MatMul_RowRow_Blas(blas_handle, nao_c, nao_s, nao_c, d_cart,
                          cart2sph.d_cart2sph_mat, cart2sph.d_cart2sph_1e_tmp);
    Launch_Device_Kernel(QC_Cart2Sph_MatMul_UT_RowRow_Kernel,
                         (total + threads - 1) / threads, threads, 0, 0, nao_s,
                         nao_s, nao_c, cart2sph.d_cart2sph_mat,
                         cart2sph.d_cart2sph_1e_tmp, d_sph);
}

void QUANTUM_CHEMISTRY::Cart2Sph_OneE_Integrals()
{
    if (!mol.is_spherical) return;
    Cart2Sph_Single_Matrix(cart2sph.d_S_cart, scf_ws.core.d_S);
    Cart2Sph_Single_Matrix(cart2sph.d_T_cart, scf_ws.core.d_T);
    Cart2Sph_Single_Matrix(cart2sph.d_V_cart, scf_ws.core.d_V);
}
