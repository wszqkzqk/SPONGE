#pragma once

#include "../one_e.hpp"

// RI two-/three-center energy and gradient kernels share this dynamic
// McMurchie-Davidson workspace implementation.  Keeping the tensor layout and
// recurrences here prevents the energy and gradient paths from drifting apart.

struct QC_RI_DEVICE_WORKSPACE
{
    int n_workers = 0;
    int e_dim_a = 0;
    int e_dim_b = 0;
    int e_dim_n = 0;
    size_t e_stride = 0;
    size_t r_stride = 0;
    size_t boys_stride = 0;
    size_t cart_stride = 0;
    float* e = NULL;
    float* r = NULL;
    double* boys = NULL;
    double* cart = NULL;
};

struct QC_RI_INTEGRAL_WORKSPACE
{
    QC_RI_DEVICE_WORKSPACE device;
    size_t e_bytes = 0;
    size_t r_bytes = 0;
    size_t boys_bytes = 0;
    size_t cart_bytes = 0;
};

static __host__ __device__ __forceinline__ size_t QC_RI_E_Tensor_Size(
    int dim_a, int dim_b, int dim_n)
{
    return (size_t)dim_a * (size_t)dim_b * (size_t)dim_n;
}

// Each worker simultaneously needs x/y/z tensors for the auxiliary bra and
// the auxiliary/ket partner.
static __host__ __device__ __forceinline__ size_t QC_RI_E_Workspace_Size(
    int dim_a, int dim_b, int dim_n)
{
    return 6 * QC_RI_E_Tensor_Size(dim_a, dim_b, dim_n);
}

static __host__ __device__ __forceinline__ size_t QC_RI_E_Index(
    int a, int b, int n, int dim_b, int dim_n)
{
    return ((size_t)a * (size_t)dim_b + (size_t)b) * (size_t)dim_n +
           (size_t)n;
}

static __host__ __device__ __forceinline__ size_t
QC_RI_R_Spatial_Size(int max_degree)
{
    size_t a = (size_t)max_degree + 1;
    size_t b = (size_t)max_degree + 2;
    size_t c = (size_t)max_degree + 3;
    if ((a & 1U) == 0)
        a /= 2;
    else
        b /= 2;
    if (a % 3 == 0)
        a /= 3;
    else if (b % 3 == 0)
        b /= 3;
    else
        c /= 3;
    return a * b * c;
}

// Spatial indices (t,u,v) are packed by total degree t+u+v.
static __host__ __device__ __forceinline__ size_t QC_RI_R_Spatial_Index(
    int t, int u, int v)
{
    const size_t degree = (size_t)(t + u + v);
    const size_t degree_offset =
        degree == 0 ? 0 : QC_RI_R_Spatial_Size((int)degree - 1);
    size_t a = (size_t)t;
    size_t b = 2 * degree + 3 - (size_t)t;
    if ((a & 1U) == 0)
        a /= 2;
    else
        b /= 2;
    return degree_offset + a * b + (size_t)u;
}

static __host__ __device__ __forceinline__ size_t QC_RI_R_Index(
    int t, int u, int v, int n, int n_stride)
{
    return QC_RI_R_Spatial_Index(t, u, v) * (size_t)n_stride + (size_t)n;
}

static __host__ __device__ __forceinline__ size_t QC_RI_R_Workspace_Size(
    int max_order)
{
    return QC_RI_R_Spatial_Size(max_order) * (size_t)(max_order + 1);
}

static __host__ __device__ __forceinline__ size_t
QC_RI_Boys_Workspace_Size(int max_m)
{
    return (size_t)max_m + 1;
}

static __host__ __device__ __forceinline__ size_t QC_RI_Cartesian_Count(int l)
{
    return ((size_t)l + 1) * ((size_t)l + 2) / 2;
}

static inline bool QC_RI_Checked_Add_Size(size_t a, size_t b, size_t* out)
{
    if (a > std::numeric_limits<size_t>::max() - b) return false;
    *out = a + b;
    return true;
}

static inline bool QC_RI_Checked_Mul_Size(size_t a, size_t b, size_t* out)
{
    if (a != 0 && b > std::numeric_limits<size_t>::max() / a) return false;
    *out = a * b;
    return true;
}

static inline bool QC_RI_Checked_Bytes(size_t count, size_t element_size,
                                       size_t* bytes)
{
    return QC_RI_Checked_Mul_Size(count, element_size, bytes);
}

static inline bool QC_RI_Checked_E_Size(int dim_a, int dim_b, int dim_n,
                                        size_t* size)
{
    if (dim_a <= 0 || dim_b <= 0 || dim_n <= 0) return false;
    size_t value = 0;
    if (!QC_RI_Checked_Mul_Size((size_t)dim_a, (size_t)dim_b, &value) ||
        !QC_RI_Checked_Mul_Size(value, (size_t)dim_n, &value) ||
        !QC_RI_Checked_Mul_Size(value, 6, &value))
        return false;
    *size = value;
    return true;
}

static inline bool QC_RI_Checked_R_Size(int max_order, size_t* size)
{
    if (max_order < 0) return false;
    size_t a = (size_t)max_order + 1;
    size_t b = (size_t)max_order + 2;
    size_t c = (size_t)max_order + 3;
    if ((a & 1U) == 0)
        a /= 2;
    else
        b /= 2;
    if (a % 3 == 0)
        a /= 3;
    else if (b % 3 == 0)
        b /= 3;
    else
        c /= 3;
    size_t spatial = 0;
    return QC_RI_Checked_Mul_Size(a, b, &spatial) &&
           QC_RI_Checked_Mul_Size(spatial, c, &spatial) &&
           QC_RI_Checked_Mul_Size(spatial, (size_t)max_order + 1, size);
}

static inline bool QC_RI_Checked_Cartesian_Count(int l, size_t* count)
{
    if (l < 0) return false;
    size_t value = 0;
    return QC_RI_Checked_Mul_Size((size_t)l + 1, (size_t)l + 2, &value) &&
           ((*count = value / 2) <= (size_t)std::numeric_limits<int>::max());
}

static inline bool QC_RI_Finalize_Workspace_Layout(
    int n_tasks, QC_RI_INTEGRAL_WORKSPACE* workspace)
{
    if (workspace == NULL || n_tasks < 0) return false;
    if (n_tasks == 0)
    {
        workspace->device.n_workers = 0;
        return true;
    }

    size_t e_worker_bytes = 0;
    size_t r_worker_bytes = 0;
    size_t boys_worker_bytes = 0;
    size_t cart_worker_bytes = 0;
    size_t per_worker_bytes = 0;
    if (!QC_RI_Checked_Bytes(workspace->device.e_stride, sizeof(float),
                             &e_worker_bytes) ||
        !QC_RI_Checked_Bytes(workspace->device.r_stride, sizeof(float),
                             &r_worker_bytes) ||
        !QC_RI_Checked_Bytes(workspace->device.boys_stride, sizeof(double),
                             &boys_worker_bytes) ||
        !QC_RI_Checked_Bytes(workspace->device.cart_stride, sizeof(double),
                             &cart_worker_bytes) ||
        !QC_RI_Checked_Add_Size(e_worker_bytes, r_worker_bytes,
                                &per_worker_bytes) ||
        !QC_RI_Checked_Add_Size(per_worker_bytes, boys_worker_bytes,
                                &per_worker_bytes) ||
        !QC_RI_Checked_Add_Size(per_worker_bytes, cart_worker_bytes,
                                &per_worker_bytes) ||
        per_worker_bytes == 0)
        return false;

    // Bound the temporary allocation while retaining enough workers for GPU
    // occupancy.  A worker processes tasks with a fixed stride when there are
    // more tasks than workspace slots.
    constexpr size_t workspace_budget = (size_t)64 * 1024 * 1024;
    constexpr size_t worker_cap = 4096;
    if (per_worker_bytes > workspace_budget) return false;
    size_t workers_by_budget = workspace_budget / per_worker_bytes;
    size_t n_workers = std::min((size_t)n_tasks,
                                std::min(worker_cap, workers_by_budget));
    if (n_workers == 0 ||
        n_workers > (size_t)std::numeric_limits<int>::max())
        return false;
    workspace->device.n_workers = (int)n_workers;

    return QC_RI_Checked_Mul_Size(e_worker_bytes, n_workers,
                                  &workspace->e_bytes) &&
           QC_RI_Checked_Mul_Size(r_worker_bytes, n_workers,
                                  &workspace->r_bytes) &&
           QC_RI_Checked_Mul_Size(boys_worker_bytes, n_workers,
                                  &workspace->boys_bytes) &&
           QC_RI_Checked_Mul_Size(cart_worker_bytes, n_workers,
                                  &workspace->cart_bytes);
}

static inline bool QC_RI_Build_2Center_Workspace_Layout(
    int n_tasks, int max_aux_l, QC_RI_INTEGRAL_WORKSPACE* workspace)
{
    if (workspace == NULL || max_aux_l < 0 ||
        max_aux_l > (std::numeric_limits<int>::max() - 1) / 2)
        return false;
    *workspace = QC_RI_INTEGRAL_WORKSPACE{};
    workspace->device.e_dim_a = max_aux_l + 1;
    workspace->device.e_dim_b = 1;
    workspace->device.e_dim_n = max_aux_l + 1;
    const int r_order = 2 * max_aux_l;
    if (!QC_RI_Checked_E_Size(workspace->device.e_dim_a,
                              workspace->device.e_dim_b,
                              workspace->device.e_dim_n,
                              &workspace->device.e_stride) ||
        !QC_RI_Checked_R_Size(r_order, &workspace->device.r_stride))
        return false;
    workspace->device.boys_stride = QC_RI_Boys_Workspace_Size(r_order);
    return QC_RI_Finalize_Workspace_Layout(n_tasks, workspace);
}

static inline bool QC_RI_Build_3Center_Workspace_Layout(
    int n_tasks, int max_aux_l, int max_orb_l,
    QC_RI_INTEGRAL_WORKSPACE* workspace)
{
    if (workspace == NULL || max_aux_l < 0 || max_orb_l < 0 ||
        max_aux_l > std::numeric_limits<int>::max() - 1 ||
        max_orb_l > (std::numeric_limits<int>::max() - 1) / 2 ||
        max_aux_l > std::numeric_limits<int>::max() - 2 * max_orb_l)
        return false;
    *workspace = QC_RI_INTEGRAL_WORKSPACE{};
    workspace->device.e_dim_a = std::max(max_aux_l, max_orb_l) + 1;
    workspace->device.e_dim_b = max_orb_l + 1;
    workspace->device.e_dim_n =
        std::max(max_aux_l + 1, 2 * max_orb_l + 1);
    const int r_order = max_aux_l + 2 * max_orb_l;
    size_t aux_cart = 0;
    size_t orb_cart = 0;
    size_t cart_stride = 0;
    if (!QC_RI_Checked_E_Size(workspace->device.e_dim_a,
                              workspace->device.e_dim_b,
                              workspace->device.e_dim_n,
                              &workspace->device.e_stride) ||
        !QC_RI_Checked_R_Size(r_order, &workspace->device.r_stride) ||
        !QC_RI_Checked_Cartesian_Count(max_aux_l, &aux_cart) ||
        !QC_RI_Checked_Cartesian_Count(max_orb_l, &orb_cart) ||
        !QC_RI_Checked_Mul_Size(aux_cart, orb_cart, &cart_stride) ||
        !QC_RI_Checked_Mul_Size(cart_stride, orb_cart, &cart_stride))
        return false;
    workspace->device.boys_stride = QC_RI_Boys_Workspace_Size(r_order);
    workspace->device.cart_stride = cart_stride;
    return QC_RI_Finalize_Workspace_Layout(n_tasks, workspace);
}

static inline bool QC_RI_Build_2Center_Gradient_Workspace_Layout(
    int n_tasks, int max_aux_l, QC_RI_INTEGRAL_WORKSPACE* workspace)
{
    if (workspace == NULL || max_aux_l < 0 ||
        max_aux_l > (std::numeric_limits<int>::max() - 2) / 2)
        return false;
    *workspace = QC_RI_INTEGRAL_WORKSPACE{};
    workspace->device.e_dim_a = max_aux_l + 2;
    workspace->device.e_dim_b = 1;
    workspace->device.e_dim_n = max_aux_l + 2;
    const int r_order = 2 * max_aux_l + 1;
    size_t aux_cart = 0;
    size_t cart_stride = 0;
    if (!QC_RI_Checked_E_Size(workspace->device.e_dim_a,
                              workspace->device.e_dim_b,
                              workspace->device.e_dim_n,
                              &workspace->device.e_stride) ||
        !QC_RI_Checked_R_Size(r_order, &workspace->device.r_stride) ||
        !QC_RI_Checked_Cartesian_Count(max_aux_l, &aux_cart) ||
        !QC_RI_Checked_Mul_Size(aux_cart, aux_cart, &cart_stride) ||
        !QC_RI_Checked_Mul_Size(cart_stride, 3, &cart_stride))
        return false;
    workspace->device.boys_stride = QC_RI_Boys_Workspace_Size(r_order);
    workspace->device.cart_stride = cart_stride;
    return QC_RI_Finalize_Workspace_Layout(n_tasks, workspace);
}

static inline bool QC_RI_Build_3Center_Gradient_Workspace_Layout(
    int n_tasks, int max_aux_l, int max_orb_l,
    QC_RI_INTEGRAL_WORKSPACE* workspace)
{
    if (workspace == NULL || max_aux_l < 0 || max_orb_l < 0 ||
        max_aux_l > std::numeric_limits<int>::max() - 2 ||
        max_orb_l > (std::numeric_limits<int>::max() - 3) / 2 ||
        max_aux_l >
            std::numeric_limits<int>::max() - 1 - 2 * max_orb_l)
        return false;
    *workspace = QC_RI_INTEGRAL_WORKSPACE{};
    workspace->device.e_dim_a = std::max(max_aux_l, max_orb_l) + 2;
    workspace->device.e_dim_b = max_orb_l + 2;
    workspace->device.e_dim_n =
        std::max(max_aux_l + 2, 2 * max_orb_l + 3);
    const int r_order = max_aux_l + 2 * max_orb_l + 1;
    size_t aux_cart = 0;
    size_t orb_cart = 0;
    size_t cart_stride = 0;
    if (!QC_RI_Checked_E_Size(workspace->device.e_dim_a,
                              workspace->device.e_dim_b,
                              workspace->device.e_dim_n,
                              &workspace->device.e_stride) ||
        !QC_RI_Checked_R_Size(r_order, &workspace->device.r_stride) ||
        !QC_RI_Checked_Cartesian_Count(max_aux_l, &aux_cart) ||
        !QC_RI_Checked_Cartesian_Count(max_orb_l, &orb_cart) ||
        !QC_RI_Checked_Mul_Size(aux_cart, orb_cart, &cart_stride) ||
        !QC_RI_Checked_Mul_Size(cart_stride, orb_cart, &cart_stride) ||
        !QC_RI_Checked_Mul_Size(cart_stride, 6, &cart_stride))
        return false;
    workspace->device.boys_stride = QC_RI_Boys_Workspace_Size(r_order);
    workspace->device.cart_stride = cart_stride;
    return QC_RI_Finalize_Workspace_Layout(n_tasks, workspace);
}

static inline void QC_RI_Free_Integral_Workspace(
    QC_RI_INTEGRAL_WORKSPACE* workspace)
{
    if (workspace == NULL) return;
    if (workspace->device.e != NULL) deviceFree(workspace->device.e);
    if (workspace->device.r != NULL) deviceFree(workspace->device.r);
    if (workspace->device.boys != NULL) deviceFree(workspace->device.boys);
    if (workspace->device.cart != NULL) deviceFree(workspace->device.cart);
    workspace->device.e = NULL;
    workspace->device.r = NULL;
    workspace->device.boys = NULL;
    workspace->device.cart = NULL;
}

static inline bool QC_RI_Allocate_Integral_Workspace(
    QC_RI_INTEGRAL_WORKSPACE* workspace)
{
    if (workspace == NULL) return false;
    if (workspace->device.n_workers == 0) return true;
    if (!Device_Malloc_Safely((void**)&workspace->device.e,
                              workspace->e_bytes) ||
        !Device_Malloc_Safely((void**)&workspace->device.r,
                              workspace->r_bytes) ||
        !Device_Malloc_Safely((void**)&workspace->device.boys,
                              workspace->boys_bytes) ||
        (workspace->cart_bytes != 0 &&
         !Device_Malloc_Safely((void**)&workspace->device.cart,
                               workspace->cart_bytes)))
    {
        QC_RI_Free_Integral_Workspace(workspace);
        return false;
    }
    return true;
}

static __device__ void QC_RI_Compute_MD_Coeffs(
    float* E, int dim_a, int dim_b, int dim_n, int la_max, int lb_max,
    float PA, float PB, float one_over_2p)
{
    const size_t tensor_size = QC_RI_E_Tensor_Size(dim_a, dim_b, dim_n);
    for (size_t i = 0; i < tensor_size; i++) E[i] = 0.0f;
    E[QC_RI_E_Index(0, 0, 0, dim_b, dim_n)] = 1.0f;
    for (int la = 0; la <= la_max; la++)
    {
        for (int lb = 0; lb <= lb_max; lb++)
        {
            if (la == 0 && lb == 0) continue;
            if (la > 0)
            {
                const int la_p = la - 1;
                for (int n = 0; n <= la + lb; n++)
                {
                    float value =
                        PA * E[QC_RI_E_Index(la_p, lb, n, dim_b, dim_n)];
                    if (n > 0)
                        value += one_over_2p *
                                 E[QC_RI_E_Index(la_p, lb, n - 1, dim_b,
                                                dim_n)];
                    if (n + 1 <= la_p + lb)
                        value += (float)(n + 1) *
                                 E[QC_RI_E_Index(la_p, lb, n + 1, dim_b,
                                                dim_n)];
                    E[QC_RI_E_Index(la, lb, n, dim_b, dim_n)] = value;
                }
            }
            else
            {
                const int lb_p = lb - 1;
                for (int n = 0; n <= la + lb; n++)
                {
                    float value =
                        PB * E[QC_RI_E_Index(la, lb_p, n, dim_b, dim_n)];
                    if (n > 0)
                        value += one_over_2p *
                                 E[QC_RI_E_Index(la, lb_p, n - 1, dim_b,
                                                dim_n)];
                    if (n + 1 <= la + lb_p)
                        value += (float)(n + 1) *
                                 E[QC_RI_E_Index(la, lb_p, n + 1, dim_b,
                                                dim_n)];
                    E[QC_RI_E_Index(la, lb, n, dim_b, dim_n)] = value;
                }
            }
        }
    }
}

static __device__ void QC_RI_Compute_Boys_Double(double* F, float t,
                                                  int max_m)
{
    QC_Compute_Boys_Double(F, t, max_m);
}

static __device__ void QC_RI_Compute_R_Tensor(float* R, const double* F,
                                               float alpha,
                                               const float PC[3], int max_order)
{
    const int n_stride = max_order + 1;
    const size_t total_size = QC_RI_R_Workspace_Size(max_order);
    for (size_t i = 0; i < total_size; i++) R[i] = 0.0f;

    const double minus_two_alpha = -2.0 * (double)alpha;
    double factor = 1.0;
    for (int n = 0; n <= max_order; n++)
    {
        R[QC_RI_R_Index(0, 0, 0, n, n_stride)] = (float)(factor * F[n]);
        factor *= minus_two_alpha;
    }

    for (int degree = 1; degree <= max_order; degree++)
    {
        for (int t = 0; t <= degree; t++)
        {
            for (int u = 0; u <= degree - t; u++)
            {
                const int v = degree - t - u;
                const int max_n = max_order - degree;
                for (int n = 0; n <= max_n; n++)
                {
                    double value = 0.0;
                    if (t > 0)
                    {
                        value = (double)PC[0] *
                                R[QC_RI_R_Index(t - 1, u, v, n + 1,
                                                n_stride)];
                        if (t > 1)
                            value += (double)(t - 1) *
                                     R[QC_RI_R_Index(t - 2, u, v, n + 1,
                                                     n_stride)];
                    }
                    else if (u > 0)
                    {
                        value = (double)PC[1] *
                                R[QC_RI_R_Index(t, u - 1, v, n + 1,
                                                n_stride)];
                        if (u > 1)
                            value += (double)(u - 1) *
                                     R[QC_RI_R_Index(t, u - 2, v, n + 1,
                                                     n_stride)];
                    }
                    else
                    {
                        value = (double)PC[2] *
                                R[QC_RI_R_Index(t, u, v - 1, n + 1,
                                                n_stride)];
                        if (v > 1)
                            value += (double)(v - 1) *
                                     R[QC_RI_R_Index(t, u, v - 2, n + 1,
                                                     n_stride)];
                    }
                    R[QC_RI_R_Index(t, u, v, n, n_stride)] = (float)value;
                }
            }
        }
    }
}
