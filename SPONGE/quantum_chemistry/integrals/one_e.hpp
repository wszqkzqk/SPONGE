#pragma once

#include "../../common.h"
#include "../quantum_chemistry.h"
#include "../structure/cartesian_components.hpp"
#include "boys.hpp"

__host__ __device__ __forceinline__ static int HR_IDX_RUNTIME(int t, int u,
                                                              int v, int n,
                                                              int hr_base)
{
    return (((t * hr_base + u) * hr_base + v) * hr_base + n);
}

__device__ __forceinline__ static void QC_Get_Lxyz_Device(int l, int idx,
                                                          int& lx, int& ly,
                                                          int& lz)
{
    if (!qc_cartesian::Component(l, idx, lx, ly, lz)) lx = ly = lz = -1;
}

__forceinline__ static void QC_Get_Lxyz_Host(int l, int idx, int& lx, int& ly,
                                             int& lz)
{
    if (!qc_cartesian::Component(l, idx, lx, ly, lz)) lx = ly = lz = -1;
}

static __device__ void get_overlap1d_arr(int l1, int l2, float PA, float PB,
                                         float gamma, float res[7][7])
{
    res[0][0] = sqrtf(CONSTANT_Pi / gamma);
    for (int i = 0; i <= l1; i++)
    {
        for (int j = 0; j <= l2; j++)
        {
            if (i == 0 && j == 0) continue;
            if (j == 0)
            {
                float val = PA * res[i - 1][0];
                if (i > 1) val += (float)(i - 1) * 0.5f / gamma * res[i - 2][0];
                res[i][0] = val;
            }
            else
            {
                float val = PB * res[i][j - 1];
                if (i > 0) val += (float)i * 0.5f / gamma * res[i - 1][j - 1];
                if (j > 1) val += (float)(j - 1) * 0.5f / gamma * res[i][j - 2];
                res[i][j] = val;
            }
        }
    }
}

static __device__ float get_overlap1d_val(int l1, int l2, float PA, float PB,
                                          float gamma)
{
    float res[7][7];
    get_overlap1d_arr(l1, l2, PA, PB, gamma, res);
    return res[l1][l2];
}

static __device__ float get_kin1d(int l1, int l2, float PA, float PB,
                                  float gamma, float alpha, float beta,
                                  float res[7][7])
{
    get_overlap1d_arr(l1 + 1, l2 + 1, PA, PB, gamma, res);
    float t = 2.0f * alpha * beta * res[l1 + 1][l2 + 1];
    if (l2 > 0) t -= alpha * (float)l2 * res[l1 + 1][l2 - 1];
    if (l1 > 0) t -= beta * (float)l1 * res[l1 - 1][l2 + 1];
    if (l1 > 0 && l2 > 0)
        t += 0.5f * (float)l1 * (float)l2 * res[l1 - 1][l2 - 1];
    return t;
}

// Double-precision Boys function.  The shared implementation selects upward
// recursion by an explicit error bound and otherwise uses a tail-bounded
// positive series followed by stable downward recursion.
static __device__ void compute_boys_double(double* F, float t, int max_m)
{
    QC_Compute_Boys_Double(F, t, max_m);
}

static __device__ void compute_md_coeffs(float E[6][6][11], int la_max,
                                         int lb_max, float PA, float PB,
                                         float one_over_2p)
{
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            for (int n = 0; n < 11; n++) E[i][j][n] = 0.0f;
    E[0][0][0] = 1.0f;
    for (int la = 0; la <= la_max; la++)
    {
        for (int lb = 0; lb <= lb_max; lb++)
        {
            if (la == 0 && lb == 0) continue;
            if (la > 0)
            {
                int la_p = la - 1;
                for (int n = 0; n <= la + lb; n++)
                {
                    float val = PA * E[la_p][lb][n];
                    if (n > 0) val += one_over_2p * E[la_p][lb][n - 1];
                    if ((n + 1) <= la_p + lb)
                        val += (float)(n + 1) * E[la_p][lb][n + 1];
                    E[la][lb][n] = val;
                }
            }
            else
            {
                int lb_p = lb - 1;
                for (int n = 0; n <= la + lb; n++)
                {
                    float val = PB * E[la][lb_p][n];
                    if (n > 0) val += one_over_2p * E[la][lb_p][n - 1];
                    if ((n + 1) <= la + lb_p)
                        val += (float)(n + 1) * E[la][lb_p][n + 1];
                    E[la][lb][n] = val;
                }
            }
        }
    }
}

static __device__ void compute_r_tensor_1e(float* R, double* F, float alpha,
                                           float PC[3], int L_tot)
{
    int total_size = ONEE_MD_BASE * ONEE_MD_BASE * ONEE_MD_BASE * ONEE_MD_BASE;
    for (int i = 0; i < total_size; i++) R[i] = 0.0f;

    double m2a = -2.0 * (double)alpha;
    double fac = 1.0;
    for (int n = 0; n <= L_tot; n++)
    {
        R[ONEE_MD_IDX(0, 0, 0, n)] = (float)(fac * F[n]);
        fac *= m2a;
    }

    for (int N = 1; N <= L_tot; N++)
    {
        for (int t = 0; t <= N; t++)
        {
            for (int u = 0; u <= N - t; u++)
            {
                int v = N - t - u;
                int max_n = L_tot - N;
                for (int n = 0; n <= max_n; n++)
                {
                    double val = 0.0;
                    if (t > 0)
                    {
                        val =
                            (double)PC[0] * R[ONEE_MD_IDX(t - 1, u, v, n + 1)];
                        if (t > 1)
                            val += (double)(t - 1) *
                                   R[ONEE_MD_IDX(t - 2, u, v, n + 1)];
                    }
                    else if (u > 0)
                    {
                        val =
                            (double)PC[1] * R[ONEE_MD_IDX(t, u - 1, v, n + 1)];
                        if (u > 1)
                            val += (double)(u - 1) *
                                   R[ONEE_MD_IDX(t, u - 2, v, n + 1)];
                    }
                    else if (v > 0)
                    {
                        val =
                            (double)PC[2] * R[ONEE_MD_IDX(t, u, v - 1, n + 1)];
                        if (v > 1)
                            val += (double)(v - 1) *
                                   R[ONEE_MD_IDX(t, u, v - 2, n + 1)];
                    }
                    R[ONEE_MD_IDX(t, u, v, n)] = (float)val;
                }
            }
        }
    }
}

static __global__ void OneE_Kernel(
    const int n_tasks, const QC_ONE_E_TASK* tasks, const VECTOR* centers,
    const int* l_list, const float* exps, const float* coeffs,
    const int* shell_offsets, const int* shell_sizes, const int* ao_offsets,
    const int* atm, const float* env, int natm, float* out_S, float* out_T,
    float* out_V, int nao_total)
{
    SIMPLE_DEVICE_FOR(task_id, n_tasks)
    {
        QC_ONE_E_TASK sh_idx = tasks[task_id];
        int i_sh = sh_idx.x;
        int j_sh = sh_idx.y;

        int li = l_list[i_sh], lj = l_list[j_sh];
        int ni = (li + 1) * (li + 2) / 2, nj = (lj + 1) * (lj + 2) / 2;
        int off_i = ao_offsets[i_sh], off_j = ao_offsets[j_sh];
        const VECTOR A = centers[i_sh];
        const VECTOR B = centers[j_sh];
        float Ax = A.x, Ay = A.y, Az = A.z;
        float Bx = B.x, By = B.y, Bz = B.z;
        float dist_sq = (Ax - Bx) * (Ax - Bx) + (Ay - By) * (Ay - By) +
                        (Az - Bz) * (Az - Bz);
        float res_x[7][7], res_y[7][7], res_z[7][7];

        for (int idx_i = 0; idx_i < ni; idx_i++)
        {
            for (int idx_j = 0; idx_j < nj; idx_j++)
            {
                int lx_i, ly_i, lz_i, lx_j, ly_j, lz_j;
                QC_Get_Lxyz_Device(li, idx_i, lx_i, ly_i, lz_i);
                QC_Get_Lxyz_Device(lj, idx_j, lx_j, ly_j, lz_j);
                float total_S = 0.0f, total_T = 0.0f, total_V = 0.0f;
                for (int pi = 0; pi < shell_sizes[i_sh]; pi++)
                {
                    float ei = exps[shell_offsets[i_sh] + pi];
                    float ci = coeffs[shell_offsets[i_sh] + pi];
                    for (int pj = 0; pj < shell_sizes[j_sh]; pj++)
                    {
                        float ej = exps[shell_offsets[j_sh] + pj];
                        float cj = coeffs[shell_offsets[j_sh] + pj];
                        float g = ei + ej;
                        float Kab = expf(-ei * ej / g * dist_sq);
                        float Px = (ei * Ax + ej * Bx) / g;
                        float Py = (ei * Ay + ej * By) / g;
                        float Pz = (ei * Az + ej * Bz) / g;
                        float E_x[6][6][11], E_y[6][6][11], E_z[6][6][11];
                        compute_md_coeffs(E_x, li, lj, Px - Ax, Px - Bx,
                                          0.5f / g);
                        compute_md_coeffs(E_y, li, lj, Py - Ay, Py - By,
                                          0.5f / g);
                        compute_md_coeffs(E_z, li, lj, Pz - Az, Pz - Bz,
                                          0.5f / g);
                        // S
                        total_S +=
                            ci * cj * Kab *
                            get_overlap1d_val(lx_i, lx_j, Px - Ax, Px - Bx, g) *
                            get_overlap1d_val(ly_i, ly_j, Py - Ay, Py - By, g) *
                            get_overlap1d_val(lz_i, lz_j, Pz - Az, Pz - Bz, g);

                        // T
                        float tx = get_kin1d(lx_i, lx_j, Px - Ax, Px - Bx, g,
                                             ei, ej, res_x);
                        float sx = res_x[lx_i][lx_j];
                        float ty = get_kin1d(ly_i, ly_j, Py - Ay, Py - By, g,
                                             ei, ej, res_y);
                        float sy = res_y[ly_i][ly_j];
                        float tz = get_kin1d(lz_i, lz_j, Pz - Az, Pz - Bz, g,
                                             ei, ej, res_z);
                        float sz = res_z[lz_i][lz_j];
                        total_T += ci * cj * Kab *
                                   (tx * sy * sz + sx * ty * sz + sx * sy * tz);

                        // V
                        for (int iat = 0; iat < natm; iat++)
                        {
                            int ptr_coord = atm[iat * 6 + 1];
                            float Cx = env[ptr_coord], Cy = env[ptr_coord + 1],
                                  Cz = env[ptr_coord + 2];
                            float PC2 = (Px - Cx) * (Px - Cx) +
                                        (Py - Cy) * (Py - Cy) +
                                        (Pz - Cz) * (Pz - Cz);
                            float PC[3] = {Px - Cx, Py - Cy, Pz - Cz};
                            int L_tot = li + lj;
                            double F_vals[ONEE_MD_BASE];
                            float R_vals[ONEE_MD_BASE * ONEE_MD_BASE *
                                         ONEE_MD_BASE * ONEE_MD_BASE];
                            compute_boys_double(F_vals, g * PC2, L_tot);
                            compute_r_tensor_1e(R_vals, F_vals, g, PC, L_tot);

                            double v_sum = 0.0;
                            for (int tx = 0; tx <= lx_i + lx_j; tx++)
                            {
                                float ex = E_x[lx_i][lx_j][tx];
                                if (ex == 0.0f) continue;
                                for (int ty = 0; ty <= ly_i + ly_j; ty++)
                                {
                                    float ey = E_y[ly_i][ly_j][ty];
                                    if (ey == 0.0f) continue;
                                    for (int tz = 0; tz <= lz_i + lz_j; tz++)
                                    {
                                        float ez = E_z[lz_i][lz_j][tz];
                                        if (ez == 0.0f) continue;
                                        v_sum += (double)ex * (double)ey *
                                                 (double)ez *
                                                 (double)R_vals[ONEE_MD_IDX(
                                                     tx, ty, tz, 0)];
                                    }
                                }
                            }
                            total_V += ci * cj * Kab * -(float)atm[iat * 6] *
                                       (2.0f * CONSTANT_Pi / g) * (float)v_sum;
                        }
                    }
                }
                int idx = (int)(off_i + idx_i) * nao_total + (off_j + idx_j);
                out_S[idx] = total_S;
                out_T[idx] = total_T;
                out_V[idx] = total_V;
            }
        }
    }
}
