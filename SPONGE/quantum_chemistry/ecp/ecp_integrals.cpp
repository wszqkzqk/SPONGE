// ECP 积分求值 (Type-1 local + Type-2 semi-local)
//
// Local n_k=2 terms retain the exact three-center Gaussian-overlap fast path.
// Other local radial powers and every semi-local term use the complete
// ECP-center angular series in structure/ecp.h. In particular, a projector-l
// integral contains q=l,l+2,... components; keeping only q=l is not an ECP
// integral for displaced or Cartesian higher-angular-momentum Gaussians.

// clang-format off
#include "../quantum_chemistry.h"
#include "../integrals/one_e.hpp"
#include "../structure/cart2sph.hpp"
#include "ecp_error_policy.hpp"
#include "ecp_integrals.h"
// clang-format on

// Numerical convergence estimates are checked after contraction, because ECP
// and basis coefficients can amplify a tiny primitive estimate by six or more
// orders of magnitude. Series values that would consume a material fraction
// of the final budget are reevaluated by the independent quadrature path. The
// tail/coarse-fine differences are practical estimators, not formal analytic
// upper bounds, so diagnostics deliberately call them estimates.
static __host__ __device__ __forceinline__ bool ecp_double_is_finite(
    double value)
{
    return qc_ecp_error_policy::Double_Is_Finite(value);
}

static __host__ __device__ bool ecp_contracted_error_is_acceptable(
    double magnitude, double error)
{
    return qc_ecp_error_policy::Contracted_Error_Is_Acceptable(magnitude,
                                                                error);
}

static __device__ bool ecp_primitive_result_is_usable(
    const QC_ECP_PRIMITIVE_RESULT& result)
{
    return result.converged && ecp_double_is_finite(result.value) &&
           ecp_double_is_finite(result.last_layer_bound) &&
           result.last_layer_bound >= 0.0;
}

static __device__ bool ecp_derivative_result_is_usable(
    const QC_ECP_PRIMITIVE_DERIVATIVE_RESULT& result)
{
    if (!result.converged || !ecp_double_is_finite(result.estimated_error) ||
        result.estimated_error < 0.0)
        return false;
    for (int direction = 0; direction < 3; ++direction)
    {
        if (!ecp_double_is_finite(result.derivative_i[direction]) ||
            !ecp_double_is_finite(result.derivative_j[direction]) ||
            !ecp_double_is_finite(result.derivative_i_error[direction]) ||
            !ecp_double_is_finite(result.derivative_j_error[direction]) ||
            result.derivative_i_error[direction] < 0.0 ||
            result.derivative_j_error[direction] < 0.0)
            return false;
    }
    return true;
}

static __device__ QC_ECP_EVALUATION_FAILURE_KIND ecp_primitive_failure_kind(
    bool converged)
{
    return converged ? QC_ECP_NONFINITE_VALUE_OR_ESTIMATE
                     : QC_ECP_PRIMITIVE_QUADRATURE_NOT_CONVERGED;
}

static __device__ bool ecp_series_error_is_small_enough(
    const QC_ECP_PRIMITIVE_RESULT& result, double contraction_scale)
{
    return qc_ecp_error_policy::Series_Estimate_Is_Small_Enough(
        ecp_primitive_result_is_usable(result), contraction_scale,
        result.last_layer_bound);
}

static __device__ QC_ECP_PRIMITIVE_RESULT ecp_semilocal_complete(
    float ei, float ej, float Ax, float Ay, float Az, float Bx, float By,
    float Bz, float Cx, float Cy, float Cz, float zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j, int ecp_l,
    double contraction_scale = 1.0)
{
    const double dx_i = static_cast<double>(Cx) - static_cast<double>(Ax);
    const double dy_i = static_cast<double>(Cy) - static_cast<double>(Ay);
    const double dz_i = static_cast<double>(Cz) - static_cast<double>(Az);
    const double dx_j = static_cast<double>(Cx) - static_cast<double>(Bx);
    const double dy_j = static_cast<double>(Cy) - static_cast<double>(By);
    const double dz_j = static_cast<double>(Cz) - static_cast<double>(Bz);
    const double eta =
        static_cast<double>(ei) + static_cast<double>(ej) + zeta;
    const double distance_i2 = dx_i * dx_i + dy_i * dy_i + dz_i * dz_i;
    const double distance_j2 = dx_j * dx_j + dy_j * dy_j + dz_j * dz_j;
    const double extent = qc_ecp_math::Semilocal_Series_Extent(
        static_cast<double>(ei), static_cast<double>(ej), distance_i2,
        distance_j2, eta, lx_i + ly_i + lz_i, lx_j + ly_j + lz_j, ecp_l);
    const int order = qc_ecp_math::Semilocal_Series_Order(
        static_cast<double>(ei), static_cast<double>(ej),
        distance_i2, distance_j2, eta,
        lx_i + ly_i + lz_i, lx_j + ly_j + lz_j, ecp_l);
    if (extent <= static_cast<double>(QC_ECP_MAX_SERIES_ORDER))
    {
        const QC_ECP_PRIMITIVE_RESULT series =
            qc_ecp_math::Semilocal_Primitive(
                static_cast<double>(ei), static_cast<double>(ej),
                static_cast<double>(Ax), static_cast<double>(Ay),
                static_cast<double>(Az), static_cast<double>(Bx),
                static_cast<double>(By), static_cast<double>(Bz),
                static_cast<double>(Cx), static_cast<double>(Cy),
                static_cast<double>(Cz), static_cast<double>(zeta), n_k, lx_i,
                ly_i, lz_i, lx_j, ly_j, lz_j, ecp_l, order);
        if (ecp_series_error_is_small_enough(series, contraction_scale))
            return series;
    }
    return qc_ecp_math::Semilocal_Primitive_Quadrature(
        static_cast<double>(ei), static_cast<double>(ej),
        static_cast<double>(Ax), static_cast<double>(Ay),
        static_cast<double>(Az), static_cast<double>(Bx),
        static_cast<double>(By), static_cast<double>(Bz),
        static_cast<double>(Cx), static_cast<double>(Cy),
        static_cast<double>(Cz), static_cast<double>(zeta), n_k, lx_i, ly_i,
        lz_i, lx_j, ly_j, lz_j, ecp_l);
}

static __device__ QC_ECP_PRIMITIVE_RESULT ecp_local_complete(
    float ei, float ej, float Ax, float Ay, float Az, float Bx, float By,
    float Bz, float Cx, float Cy, float Cz, float zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j,
    double contraction_scale = 1.0)
{
    if (n_k == 2)
    {
        return qc_ecp_math::Local_Primitive_N2(
            static_cast<double>(ei), static_cast<double>(ej),
            static_cast<double>(Ax), static_cast<double>(Ay),
            static_cast<double>(Az), static_cast<double>(Bx),
            static_cast<double>(By), static_cast<double>(Bz),
            static_cast<double>(Cx), static_cast<double>(Cy),
            static_cast<double>(Cz), static_cast<double>(zeta), lx_i, ly_i,
            lz_i, lx_j, ly_j, lz_j);
    }
    const double dx_i = static_cast<double>(Cx) - static_cast<double>(Ax);
    const double dy_i = static_cast<double>(Cy) - static_cast<double>(Ay);
    const double dz_i = static_cast<double>(Cz) - static_cast<double>(Az);
    const double dx_j = static_cast<double>(Cx) - static_cast<double>(Bx);
    const double dy_j = static_cast<double>(Cy) - static_cast<double>(By);
    const double dz_j = static_cast<double>(Cz) - static_cast<double>(Bz);
    const double eta =
        static_cast<double>(ei) + static_cast<double>(ej) + zeta;
    const int angular_degree =
        lx_i + ly_i + lz_i + lx_j + ly_j + lz_j;
    const double extent = qc_ecp_math::Local_Series_Extent(
        static_cast<double>(ei), static_cast<double>(ej), dx_i, dy_i, dz_i,
        dx_j, dy_j, dz_j, eta, angular_degree);
    if (extent <= static_cast<double>(QC_ECP_MAX_SERIES_ORDER))
    {
        const int order = qc_ecp_math::Local_Series_Order(
            static_cast<double>(ei), static_cast<double>(ej), dx_i, dy_i, dz_i,
            dx_j, dy_j, dz_j, eta, angular_degree);
        const QC_ECP_PRIMITIVE_RESULT series =
            qc_ecp_math::Local_Primitive_Series(
                static_cast<double>(ei), static_cast<double>(ej),
                static_cast<double>(Ax), static_cast<double>(Ay),
                static_cast<double>(Az), static_cast<double>(Bx),
                static_cast<double>(By), static_cast<double>(Bz),
                static_cast<double>(Cx), static_cast<double>(Cy),
                static_cast<double>(Cz), static_cast<double>(zeta), n_k, lx_i,
                ly_i, lz_i, lx_j, ly_j, lz_j, order);
        if (ecp_series_error_is_small_enough(series, contraction_scale))
            return series;
    }
    return qc_ecp_math::Local_Primitive_Quadrature(
        static_cast<double>(ei), static_cast<double>(ej),
        static_cast<double>(Ax), static_cast<double>(Ay),
        static_cast<double>(Az), static_cast<double>(Bx),
        static_cast<double>(By), static_cast<double>(Bz),
        static_cast<double>(Cx), static_cast<double>(Cy),
        static_cast<double>(Cz), static_cast<double>(zeta), n_k, lx_i, ly_i,
        lz_i, lx_j, ly_j, lz_j);
}

struct QC_ECP_DEVICE_FAILURE
{
    int kind;
    int task_id;
    int atom;
    int channel;
    int term;
    int primitive_i;
    int primitive_j;
    int cartesian_i;
    int cartesian_j;
    double value;
    double estimated_error;
};

static __device__ int claim_ecp_failure_kind(
    QC_ECP_DEVICE_FAILURE* failure, QC_ECP_EVALUATION_FAILURE_KIND kind)
{
#ifdef GPU_ARCH_NAME
    return atomicCAS(&failure->kind, QC_ECP_EVALUATION_OK,
                     static_cast<int>(kind));
#else
    int observed = QC_ECP_EVALUATION_OK;
#pragma omp critical(sponge_ecp_failure_claim)
    {
        observed = failure->kind;
        if (observed == QC_ECP_EVALUATION_OK)
            failure->kind = static_cast<int>(kind);
    }
    return observed;
#endif
}

static __device__ void record_ecp_failure(
    QC_ECP_DEVICE_FAILURE* failure, QC_ECP_EVALUATION_FAILURE_KIND kind,
    int task_id, int atom, int channel, int term, int primitive_i,
    int primitive_j, int cartesian_i, int cartesian_j, double value,
    double estimated_error)
{
    if (claim_ecp_failure_kind(failure, kind) != QC_ECP_EVALUATION_OK)
        return;
    failure->task_id = task_id;
    failure->atom = atom;
    failure->channel = channel;
    failure->term = term;
    failure->primitive_i = primitive_i;
    failure->primitive_j = primitive_j;
    failure->cartesian_i = cartesian_i;
    failure->cartesian_j = cartesian_j;
    failure->value = value;
    failure->estimated_error = estimated_error;
}

// 找到 local channel (l == l_max 或 l < 0)
static __device__ int find_local_channel(int ch_start, int ch_end, int l_max,
                                         const int* ecp_channel_l)
{
    for (int ich = ch_start; ich < ch_end; ich++)
    {
        int cl = ecp_channel_l[ich];
        if (cl < 0 || cl == l_max) return ich;
    }
    return -1;
}

// ECP 主 Kernel
static __global__ void ECP_Kernel(
    const int n_tasks, const int task_offset, const QC_ONE_E_TASK* tasks,
    const VECTOR* centers,
    const int* l_list, const float* exps_arr, const float* coeffs_arr,
    const int* shell_offsets, const int* shell_sizes, const int* ao_offsets,
    // ECP 参数
    const VECTOR* atom_coords, int natm, const int* ecp_l_max,
    const int* ecp_atom_channel_range, const int* ecp_channel_l,
    const int* ecp_channel_offsets, const int* ecp_channel_sizes,
    const float* ecp_d, const float* ecp_zeta, const int* ecp_n,
    // 输出
    float* out_V_ECP, int nao_total, QC_ECP_DEVICE_FAILURE* failure)
{
    SIMPLE_DEVICE_FOR(task_id, n_tasks)
    {
        const int global_task_id = task_offset + task_id;
        QC_ONE_E_TASK sh_idx = tasks[task_id];
        const int i_sh = sh_idx.x;
        const int j_sh = sh_idx.y;

        const int li = l_list[i_sh], lj = l_list[j_sh];
        const int ni = (li + 1) * (li + 2) / 2;
        const int nj = (lj + 1) * (lj + 2) / 2;
        const int off_i = ao_offsets[i_sh], off_j = ao_offsets[j_sh];
        const VECTOR A = centers[i_sh], B = centers[j_sh];
        const float Ax = A.x, Ay = A.y, Az = A.z;
        const float Bx = B.x, By = B.y, Bz = B.z;
        for (int idx_i = 0; idx_i < ni; idx_i++)
        {
            for (int idx_j = 0; idx_j < nj; idx_j++)
            {
                int lx_i, ly_i, lz_i, lx_j, ly_j, lz_j;
                QC_Get_Lxyz_Device(li, idx_i, lx_i, ly_i, lz_i);
                QC_Get_Lxyz_Device(lj, idx_j, lx_j, ly_j, lz_j);

                double total_ecp = 0.0;
                double total_error_bound = 0.0;
                double largest_weighted_error = -1.0;
                int error_atom = -1, error_channel = -1, error_term = -1;
                int error_primitive_i = -1, error_primitive_j = -1;

                // 遍历 ECP 原子
                for (int iat = 0; iat < natm; iat++)
                {
                    if (ecp_l_max[iat] < 0) continue;  // 该原子无 ECP

                    const float Cx = atom_coords[iat].x;
                    const float Cy = atom_coords[iat].y;
                    const float Cz = atom_coords[iat].z;
                    const int l_max = ecp_l_max[iat];
                    const int ch_start = ecp_atom_channel_range[iat];
                    const int ch_end = ecp_atom_channel_range[iat + 1];

                    int local_ch = find_local_channel(ch_start, ch_end, l_max,
                                                      ecp_channel_l);

                    for (int pi = 0; pi < shell_sizes[i_sh]; pi++)
                    {
                        const float ei = exps_arr[shell_offsets[i_sh] + pi];
                        const float ci = coeffs_arr[shell_offsets[i_sh] + pi];

                        for (int pj = 0; pj < shell_sizes[j_sh]; pj++)
                        {
                            const float ej = exps_arr[shell_offsets[j_sh] + pj];
                            const float cj =
                                coeffs_arr[shell_offsets[j_sh] + pj];
                            const double cc = static_cast<double>(ci) *
                                              static_cast<double>(cj);

                            // 1. Local 贡献: ⟨μ|U_L|ν⟩ (无投影)
                            if (local_ch >= 0)
                            {
                                const int t_off = ecp_channel_offsets[local_ch];
                                const int t_cnt = ecp_channel_sizes[local_ch];
                                for (int it = 0; it < t_cnt; it++)
                                {
                                    const float dk = ecp_d[t_off + it];
                                    const float zk = ecp_zeta[t_off + it];
                                    const double factor =
                                        cc * static_cast<double>(dk);
                                    if (!ecp_double_is_finite(factor))
                                    {
                                        record_ecp_failure(
                                            failure,
                                            QC_ECP_NONFINITE_VALUE_OR_ESTIMATE,
                                            global_task_id, iat, local_ch,
                                            t_off + it, pi, pj, idx_i, idx_j,
                                            factor, DBL_MAX);
                                        continue;
                                    }
                                    const QC_ECP_PRIMITIVE_RESULT evaluated =
                                        ecp_local_complete(
                                            ei, ej, Ax, Ay, Az, Bx, By, Bz, Cx,
                                            Cy, Cz, zk, ecp_n[t_off + it], lx_i,
                                            ly_i, lz_i, lx_j, ly_j, lz_j,
                                            qc_ecp_math::Abs(factor));
                                    if (!ecp_primitive_result_is_usable(
                                            evaluated))
                                    {
                                        record_ecp_failure(
                                            failure,
                                            ecp_primitive_failure_kind(
                                                evaluated.converged),
                                            global_task_id, iat, local_ch,
                                            t_off + it, pi, pj, idx_i, idx_j,
                                            evaluated.value,
                                            evaluated.last_layer_bound);
                                        continue;
                                    }
                                    total_ecp += factor * evaluated.value;
                                    const double weighted_error =
                                        qc_ecp_math::Abs(factor) *
                                        evaluated.last_layer_bound;
                                    total_error_bound += weighted_error;
                                    if (weighted_error > largest_weighted_error)
                                    {
                                        largest_weighted_error = weighted_error;
                                        error_atom = iat;
                                        error_channel = local_ch;
                                        error_term = t_off + it;
                                        error_primitive_i = pi;
                                        error_primitive_j = pj;
                                    }
                                }
                            }

                            // 2. Semi-local: ⟨μ|ΔU_l P_l|ν⟩
                            //    通道数据已存储 ΔU_l = U_l - U_L
                            for (int ich = ch_start; ich < ch_end; ich++)
                            {
                                const int ch_l = ecp_channel_l[ich];
                                if (ch_l < 0 || ch_l == l_max)
                                    continue;  // skip local

                                const int t_off = ecp_channel_offsets[ich];
                                const int t_cnt = ecp_channel_sizes[ich];
                                for (int it = 0; it < t_cnt; it++)
                                {
                                    const float dk = ecp_d[t_off + it];
                                    const float zk = ecp_zeta[t_off + it];
                                    const double factor =
                                        cc * static_cast<double>(dk);
                                    if (!ecp_double_is_finite(factor))
                                    {
                                        record_ecp_failure(
                                            failure,
                                            QC_ECP_NONFINITE_VALUE_OR_ESTIMATE,
                                            global_task_id, iat, ich,
                                            t_off + it, pi, pj, idx_i, idx_j,
                                            factor, DBL_MAX);
                                        continue;
                                    }
                                    const QC_ECP_PRIMITIVE_RESULT evaluated =
                                        ecp_semilocal_complete(
                                            ei, ej, Ax, Ay, Az, Bx, By, Bz, Cx,
                                            Cy, Cz, zk, ecp_n[t_off + it], lx_i,
                                            ly_i, lz_i, lx_j, ly_j, lz_j,
                                            ch_l, qc_ecp_math::Abs(factor));
                                    if (!ecp_primitive_result_is_usable(
                                            evaluated))
                                    {
                                        record_ecp_failure(
                                            failure,
                                            ecp_primitive_failure_kind(
                                                evaluated.converged),
                                            global_task_id, iat, ich,
                                            t_off + it, pi, pj, idx_i, idx_j,
                                            evaluated.value,
                                            evaluated.last_layer_bound);
                                        continue;
                                    }
                                    total_ecp += factor * evaluated.value;
                                    const double weighted_error =
                                        qc_ecp_math::Abs(factor) *
                                        evaluated.last_layer_bound;
                                    total_error_bound += weighted_error;
                                    if (weighted_error > largest_weighted_error)
                                    {
                                        largest_weighted_error = weighted_error;
                                        error_atom = iat;
                                        error_channel = ich;
                                        error_term = t_off + it;
                                        error_primitive_i = pi;
                                        error_primitive_j = pj;
                                    }
                                }
                            }
                        }
                    }
                }

                // 1e task 列表包含全部 (i,j) 和 (j,i)，无需手动对称化
                const int idx = (off_i + idx_i) * nao_total + (off_j + idx_j);
                const qc_ecp_error_policy::Matrix_Storage_Assessment storage =
                    qc_ecp_error_policy::Assess_Matrix_Storage(
                        total_ecp, total_error_bound);
                if (!storage.accepted)
                {
                    QC_ECP_EVALUATION_FAILURE_KIND kind =
                        QC_ECP_MATRIX_STORAGE_EXCEEDS_BUDGET;
                    if (!ecp_double_is_finite(total_ecp) ||
                        !ecp_double_is_finite(total_error_bound))
                        kind = QC_ECP_NONFINITE_VALUE_OR_ESTIMATE;
                    else if (!ecp_contracted_error_is_acceptable(
                                 qc_ecp_math::Abs(total_ecp),
                                 total_error_bound))
                        kind = QC_ECP_MATRIX_ESTIMATE_EXCEEDS_BUDGET;
                    record_ecp_failure(
                        failure, kind, global_task_id, error_atom,
                        error_channel, error_term, error_primitive_i,
                        error_primitive_j, idx_i, idx_j, storage.stored_value,
                        storage.total_error_estimate);
                    continue;
                }
                // Every ordered shell pair appears exactly once in the 1e task
                // list and Cartesian AO ranges do not overlap, so this element
                // has a single writer. Assignment avoids an unnecessary
                // second float rounding through atomic addition with zero.
                out_V_ECP[idx] = static_cast<float>(storage.stored_value);
            }
        }
    }
}

// ECP 积分驱动
static void QC_Copy_ECP_Failure_Context(
    const QC_MOLECULE& mol, const QC_INTEGRAL_TASKS& task_ctx,
    const QC_ECP_DEVICE_FAILURE& device_failure,
    QC_ECP_EVALUATION_FAILURE* failure)
{
    if (failure == nullptr) return;
    failure->kind = static_cast<QC_ECP_EVALUATION_FAILURE_KIND>(
        device_failure.kind);
    failure->task_id = device_failure.task_id;
    if (device_failure.task_id >= 0 &&
        device_failure.task_id < task_ctx.topo.n_1e_tasks)
    {
        const QC_ONE_E_TASK& task =
            task_ctx.topo.h_1e_tasks[device_failure.task_id];
        failure->shell_i = task.x;
        failure->shell_j = task.y;
    }
    failure->atom = device_failure.atom;
    failure->channel = device_failure.channel;
    failure->term = device_failure.term;
    failure->primitive_i = device_failure.primitive_i;
    failure->primitive_j = device_failure.primitive_j;
    failure->cartesian_i = device_failure.cartesian_i;
    failure->cartesian_j = device_failure.cartesian_j;
    failure->value = device_failure.value;
    failure->estimated_error = device_failure.estimated_error;
    if (device_failure.channel >= 0 &&
        device_failure.channel < static_cast<int>(mol.h_ecp_l.size()))
        failure->channel_l = mol.h_ecp_l[device_failure.channel];
    if (device_failure.term >= 0 &&
        device_failure.term < static_cast<int>(mol.h_ecp_n.size()))
        failure->n_k = mol.h_ecp_n[device_failure.term];
}

static bool QC_Allocate_ECP_Failure(QC_ECP_DEVICE_FAILURE** device_failure)
{
    QC_ECP_DEVICE_FAILURE initial = {};
    initial.kind = QC_ECP_EVALUATION_OK;
    initial.task_id = initial.atom = initial.channel = initial.term = -1;
    initial.primitive_i = initial.primitive_j = -1;
    initial.cartesian_i = initial.cartesian_j = -1;
    initial.value = initial.estimated_error = 0.0;
    if (!Device_Malloc_Safely(reinterpret_cast<void**>(device_failure),
                              sizeof(initial)))
        return false;
    deviceMemcpy(*device_failure, &initial, sizeof(initial),
                 deviceMemcpyHostToDevice);
    return true;
}

static bool QC_Report_ECP_Resource_Failure(
    QC_ECP_EVALUATION_FAILURE* failure)
{
    if (failure != nullptr)
        failure->kind = QC_ECP_RESOURCE_ALLOCATION_FAILED;
    return false;
}

static bool QC_Finalize_ECP_Failure(
    const QC_MOLECULE& mol, const QC_INTEGRAL_TASKS& task_ctx,
    QC_ECP_DEVICE_FAILURE* device_failure,
    QC_ECP_EVALUATION_FAILURE* failure)
{
    QC_ECP_DEVICE_FAILURE observed = {};
    deviceMemcpy(&observed, device_failure, sizeof(observed),
                 deviceMemcpyDeviceToHost);
    deviceFree(device_failure);
    if (observed.kind == QC_ECP_EVALUATION_OK) return true;
    QC_Copy_ECP_Failure_Context(mol, task_ctx, observed, failure);
    return false;
}

static __global__ void ECP_Commit_Gradient_Kernel(
    int component_count, const double* ecp_gradient, double* gradient)
{
    SIMPLE_DEVICE_FOR(component, component_count)
    {
        gradient[component] += ecp_gradient[component];
    }
}

static bool QC_Finalize_ECP_Gradient(
    const QC_MOLECULE& mol, const QC_INTEGRAL_TASKS& task_ctx,
    QC_ECP_DEVICE_FAILURE* device_failure, double* device_certificate,
    QC_ECP_DEVICE_FAILURE* device_context, std::size_t component_count,
    double* committed_gradient, QC_ECP_EVALUATION_FAILURE* failure)
{
    QC_ECP_DEVICE_FAILURE observed = {};
    std::vector<double> certificate(2 * component_count, 0.0);
    std::vector<QC_ECP_DEVICE_FAILURE> context(component_count);
    deviceMemcpy(&observed, device_failure, sizeof(observed),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(certificate.data(), device_certificate,
                 certificate.size() * sizeof(double),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(context.data(), device_context,
                 context.size() * sizeof(QC_ECP_DEVICE_FAILURE),
                 deviceMemcpyDeviceToHost);
    if (observed.kind != QC_ECP_EVALUATION_OK)
    {
        QC_Copy_ECP_Failure_Context(mol, task_ctx, observed, failure);
        deviceFree(device_failure);
        deviceFree(device_certificate);
        deviceFree(device_context);
        return false;
    }

    const double* error = certificate.data();
    const double* observable = certificate.data() + component_count;
    for (std::size_t component = 0; component < component_count; ++component)
    {
        const qc_ecp_error_policy::Gradient_Assessment assessment =
            qc_ecp_error_policy::Assess_Gradient(observable[component],
                                                 error[component]);
        if (assessment.accepted)
            continue;

        QC_ECP_DEVICE_FAILURE aggregate = context[component];
        aggregate.kind =
            ecp_double_is_finite(observable[component]) &&
                    ecp_double_is_finite(error[component])
                ? QC_ECP_GRADIENT_ESTIMATE_EXCEEDS_BUDGET
                : QC_ECP_NONFINITE_VALUE_OR_ESTIMATE;
        aggregate.value = observable[component];
        aggregate.estimated_error = error[component];
        QC_Copy_ECP_Failure_Context(mol, task_ctx, aggregate, failure);
        if (failure != nullptr)
        {
            failure->observable_atom = static_cast<int>(component / 3);
            failure->direction = static_cast<int>(component % 3);
        }
        deviceFree(device_failure);
        deviceFree(device_certificate);
        deviceFree(device_context);
        return false;
    }

    const int count = static_cast<int>(component_count);
    const double* device_observable = device_certificate + component_count;
    Launch_Device_Kernel(ECP_Commit_Gradient_Kernel, (count + 63) / 64, 64, 0,
                         0, count, device_observable, committed_gradient);
    deviceFree(device_failure);
    deviceFree(device_certificate);
    deviceFree(device_context);
    return true;
}

bool QC_Compute_V_ECP(const QC_MOLECULE& mol,
                      const QC_INTEGRAL_TASKS& task_ctx, float* d_V_ECP,
                      QC_ECP_EVALUATION_FAILURE* failure)
{
    if (!mol.has_ecp || mol.ecp_total_terms == 0) return true;

    QC_ECP_DEVICE_FAILURE* device_failure = nullptr;
    if (!QC_Allocate_ECP_Failure(&device_failure))
        return QC_Report_ECP_Resource_Failure(failure);

    const int n_total = task_ctx.topo.n_1e_tasks;
    const int chunk_size = ONE_E_BATCH_SIZE;

    for (int i = 0; i < n_total; i += chunk_size)
    {
        int current_chunk = std::min(chunk_size, n_total - i);
        const QC_ONE_E_TASK* task_ptr = task_ctx.buffers.d_1e_tasks + i;
        Launch_Device_Kernel(
            ECP_Kernel, (current_chunk + 63) / 64, 64, 0, 0, current_chunk, i,
            task_ptr, mol.d_centers, mol.d_l_list, mol.d_exps, mol.d_coeffs,
            mol.d_shell_offsets, mol.d_shell_sizes, mol.d_ao_offsets,
            mol.d_atom_coords, mol.natm, mol.d_ecp_l_max,
            mol.d_ecp_atom_channel_range, mol.d_ecp_l,
            mol.d_ecp_channel_offsets, mol.d_ecp_channel_sizes, mol.d_ecp_d,
            mol.d_ecp_zeta, mol.d_ecp_n, d_V_ECP, mol.nao_cart,
            device_failure);
    }
    return QC_Finalize_ECP_Failure(mol, task_ctx, device_failure, failure);
}

// ECP 梯度 Kernel
// 使用角动量平移求导: d/dA_x V = 2α V(lx_i+1) - lx_i V(lx_i-1)
// ECP 中心导数由平移不变性得: d/dC = -(d/dA + d/dB)
// 遵循 V_Grad_Kernel 约定:
//   bra 导数 × 2.0 × P → atom_i (因子 2 包含 ket 贡献)
//   ECP 中心导数 × 1.0 × P → ecp atom

static __device__ QC_ECP_PRIMITIVE_RESULT ecp_integral_for_term(
    float ei, float ej, float Ax, float Ay, float Az, float Bx, float By,
    float Bz, float Cx, float Cy, float Cz, float dist_sq_AB, float zeta,
    int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j, int lz_j,
    int ch_l, bool is_local, double contraction_scale)
{
    (void)dist_sq_AB;
    if (lx_i < 0 || ly_i < 0 || lz_i < 0 || lx_j < 0 || ly_j < 0 ||
        lz_j < 0)
    {
        QC_ECP_PRIMITIVE_RESULT zero;
        zero.converged = true;
        return zero;
    }
    if (is_local)
        return ecp_local_complete(ei, ej, Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz,
                                  zeta, n_k, lx_i, ly_i, lz_i, lx_j, ly_j,
                                  lz_j, contraction_scale);
    return ecp_semilocal_complete(ei, ej, Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz,
                                  zeta, n_k, lx_i, ly_i, lz_i, lx_j, ly_j,
                                  lz_j, ch_l, contraction_scale);
}

static __device__ QC_ECP_PRIMITIVE_RESULT ecp_integral_or_record_failure(
    float ei, float ej, float Ax, float Ay, float Az, float Bx, float By,
    float Bz, float Cx, float Cy, float Cz, float dist_sq_AB, float zeta,
    int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j, int lz_j,
    int ch_l, bool is_local, QC_ECP_DEVICE_FAILURE* failure, int task_id,
    int atom, int channel, int term, int primitive_i, int primitive_j,
    int cartesian_i, int cartesian_j, double contraction_scale)
{
    const QC_ECP_PRIMITIVE_RESULT evaluated = ecp_integral_for_term(
        ei, ej, Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, dist_sq_AB, zeta, n_k,
        lx_i, ly_i, lz_i, lx_j, ly_j, lz_j, ch_l, is_local,
        contraction_scale);
    if (!ecp_primitive_result_is_usable(evaluated))
        record_ecp_failure(
            failure, ecp_primitive_failure_kind(evaluated.converged), task_id,
            atom, channel, term, primitive_i, primitive_j, cartesian_i,
            cartesian_j, evaluated.value, evaluated.last_layer_bound);
    return evaluated;
}

static __device__ bool ecp_gradient_requires_direct_quadrature(
    float ei, float ej, float Ax, float Ay, float Az, float Bx, float By,
    float Bz, float Cx, float Cy, float Cz, float zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j, int ch_l,
    bool is_local)
{
    const double dx_i = static_cast<double>(Cx) - static_cast<double>(Ax);
    const double dy_i = static_cast<double>(Cy) - static_cast<double>(Ay);
    const double dz_i = static_cast<double>(Cz) - static_cast<double>(Az);
    const double dx_j = static_cast<double>(Cx) - static_cast<double>(Bx);
    const double dy_j = static_cast<double>(Cy) - static_cast<double>(By);
    const double dz_j = static_cast<double>(Cz) - static_cast<double>(Bz);
    const double eta =
        static_cast<double>(ei) + static_cast<double>(ej) + zeta;
    if (is_local)
    {
        // n_k=2 has an exact Gaussian-overlap path for every raised/lowered
        // integral and is always cheaper than numerical quadrature.
        if (n_k == 2) return false;
        const double extent = qc_ecp_math::Local_Series_Extent(
            static_cast<double>(ei), static_cast<double>(ej), dx_i, dy_i, dz_i,
            dx_j, dy_j, dz_j, eta,
            lx_i + ly_i + lz_i + lx_j + ly_j + lz_j + 1);
        return extent > static_cast<double>(QC_ECP_MAX_SERIES_ORDER);
    }
    const double distance_i2 = dx_i * dx_i + dy_i * dy_i + dz_i * dz_i;
    const double distance_j2 = dx_j * dx_j + dy_j * dy_j + dz_j * dz_j;
    const double extent = qc_ecp_math::Semilocal_Series_Extent(
        static_cast<double>(ei), static_cast<double>(ej), distance_i2,
        distance_j2, eta, lx_i + ly_i + lz_i + 1,
        lx_j + ly_j + lz_j + 1, ch_l);
    return extent > static_cast<double>(QC_ECP_MAX_SERIES_ORDER);
}

static __device__ QC_ECP_PRIMITIVE_DERIVATIVE_RESULT
ecp_direct_primitive_derivative(
    float ei, float ej, float Ax, float Ay, float Az, float Bx, float By,
    float Bz, float Cx, float Cy, float Cz, float zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j, int ch_l,
    bool is_local)
{
    if (is_local)
        return qc_ecp_math::Local_Primitive_Derivative_Quadrature(
            static_cast<double>(ei), static_cast<double>(ej),
            static_cast<double>(Ax), static_cast<double>(Ay),
            static_cast<double>(Az), static_cast<double>(Bx),
            static_cast<double>(By), static_cast<double>(Bz),
            static_cast<double>(Cx), static_cast<double>(Cy),
            static_cast<double>(Cz), static_cast<double>(zeta), n_k, lx_i,
            ly_i, lz_i, lx_j, ly_j, lz_j);
    return qc_ecp_math::Semilocal_Primitive_Derivative_Quadrature(
        static_cast<double>(ei), static_cast<double>(ej),
        static_cast<double>(Ax), static_cast<double>(Ay),
        static_cast<double>(Az), static_cast<double>(Bx),
        static_cast<double>(By), static_cast<double>(Bz),
        static_cast<double>(Cx), static_cast<double>(Cy),
        static_cast<double>(Cz), static_cast<double>(zeta), n_k, lx_i, ly_i,
        lz_i, lx_j, ly_j, lz_j, ch_l);
}

static __device__ void record_ecp_gradient_context(
    QC_ECP_DEVICE_FAILURE* context, int task_id, int ecp_atom, int channel,
    int term, int primitive_i, int primitive_j, int cartesian_i,
    int cartesian_j, double contribution, double error)
{
    if (ecp_double_is_finite(contribution) && ecp_double_is_finite(error) &&
        qc_ecp_math::Abs(contribution) == 0.0 && error == 0.0)
        return;
    record_ecp_failure(
        context, QC_ECP_GRADIENT_ESTIMATE_EXCEEDS_BUDGET, task_id, ecp_atom,
        channel, term, primitive_i, primitive_j, cartesian_i, cartesian_j,
        contribution, error);
}

static __global__ void ECP_Grad_Kernel(
    const int n_tasks, const int task_offset, const QC_ONE_E_TASK* tasks,
    const VECTOR* centers, const int* l_list, const float* exps_arr,
    const float* coeffs_arr, const int* shell_offsets, const int* shell_sizes,
    const int* ao_offsets,
    // ECP 参数
    const VECTOR* atom_coords, int natm, const int* ecp_l_max,
    const int* ecp_atom_channel_range, const int* ecp_channel_l,
    const int* ecp_channel_offsets, const int* ecp_channel_sizes,
    const float* ecp_d, const float* ecp_zeta, const int* ecp_n,
    // 密度 (已含 norm 权重)
    int nao_total, const int* shell_atom, const float* P_weighted,
    // 输出与最终 observable 误差证书
    double* gradient_error, double* gradient_observable,
    QC_ECP_DEVICE_FAILURE* gradient_context,
    QC_ECP_DEVICE_FAILURE* failure)
{
    SIMPLE_DEVICE_FOR(task_id, n_tasks)
    {
        const int global_task_id = task_offset + task_id;
        const QC_ONE_E_TASK sh_idx = tasks[task_id];
        const int i_sh = sh_idx.x;
        const int j_sh = sh_idx.y;

        const int li = l_list[i_sh], lj = l_list[j_sh];
        const int ni = (li + 1) * (li + 2) / 2;
        const int nj = (lj + 1) * (lj + 2) / 2;
        const int off_i = ao_offsets[i_sh], off_j = ao_offsets[j_sh];
        const int atom_i = shell_atom[i_sh];
        const VECTOR A = centers[i_sh], B = centers[j_sh];
        const float Ax = A.x, Ay = A.y, Az = A.z;
        const float Bx = B.x, By = B.y, Bz = B.z;
        const float dist_sq = (Ax - Bx) * (Ax - Bx) +
                              (Ay - By) * (Ay - By) +
                              (Az - Bz) * (Az - Bz);

        for (int idx_i = 0; idx_i < ni; ++idx_i)
        {
            for (int idx_j = 0; idx_j < nj; ++idx_j)
            {
                int lx_i, ly_i, lz_i, lx_j, ly_j, lz_j;
                QC_Get_Lxyz_Device(li, idx_i, lx_i, ly_i, lz_i);
                QC_Get_Lxyz_Device(lj, idx_j, lx_j, ly_j, lz_j);

                const int mu = off_i + idx_i;
                const int nu = off_j + idx_j;
                const double dp = static_cast<double>(
                    P_weighted[mu * nao_total + nu]);
                if (!ecp_double_is_finite(dp))
                {
                    record_ecp_failure(
                        failure, QC_ECP_NONFINITE_VALUE_OR_ESTIMATE,
                        global_task_id, -1, -1, -1, -1, -1, idx_i, idx_j, dp,
                        DBL_MAX);
                    continue;
                }
                const double abs_dp = qc_ecp_math::Abs(dp);

                // 遍历 ECP 原子
                for (int iat = 0; iat < natm; ++iat)
                {
                    if (ecp_l_max[iat] < 0) continue;

                    const float Cx = atom_coords[iat].x;
                    const float Cy = atom_coords[iat].y;
                    const float Cz = atom_coords[iat].z;
                    const int l_max = ecp_l_max[iat];
                    const int ch_start = ecp_atom_channel_range[iat];
                    const int ch_end = ecp_atom_channel_range[iat + 1];
                    const int local_ch = find_local_channel(
                        ch_start, ch_end, l_max, ecp_channel_l);

                    // 对每个 primitive pair 累积梯度。误差先按 ECP/basis
                    // coefficient 做绝对值收缩，随后再按密度和对称因子
                    // 累计到最终原子力 observable。
                    for (int pi = 0; pi < shell_sizes[i_sh]; ++pi)
                    {
                        const float ei = exps_arr[shell_offsets[i_sh] + pi];
                        const float ci = coeffs_arr[shell_offsets[i_sh] + pi];
                        for (int pj = 0; pj < shell_sizes[j_sh]; ++pj)
                        {
                            const float ej = exps_arr[shell_offsets[j_sh] + pj];
                            const float cj =
                                coeffs_arr[shell_offsets[j_sh] + pj];
                            const double cc = static_cast<double>(ci) *
                                              static_cast<double>(cj);
                            double dV_A[3] = {}, dV_B[3] = {};
                            double dV_A_error[3] = {}, dV_B_error[3] = {};
                            int l_i[3] = {lx_i, ly_i, lz_i};
                            int l_j[3] = {lx_j, ly_j, lz_j};

                            auto accumulate_grad = [&](float dk, float zk,
                                                       int n_k, int ch_l,
                                                       bool is_local,
                                                       int channel, int term)
                            {
                                const double cdk =
                                    cc * static_cast<double>(dk);
                                if (!ecp_double_is_finite(cdk))
                                {
                                    record_ecp_failure(
                                        failure,
                                        QC_ECP_NONFINITE_VALUE_OR_ESTIMATE,
                                        global_task_id, iat, channel, term, pi,
                                        pj, idx_i, idx_j, cdk, DBL_MAX);
                                    return;
                                }
                                const double abs_cdk = qc_ecp_math::Abs(cdk);
                                double term_A[3] = {}, term_B[3] = {};
                                double term_A_error[3] = {};
                                double term_B_error[3] = {};

                                if (ecp_gradient_requires_direct_quadrature(
                                        ei, ej, Ax, Ay, Az, Bx, By, Bz, Cx, Cy,
                                        Cz, zk, n_k, lx_i, ly_i, lz_i, lx_j,
                                        ly_j, lz_j, ch_l, is_local))
                                {
                                    const QC_ECP_PRIMITIVE_DERIVATIVE_RESULT
                                        evaluated =
                                            ecp_direct_primitive_derivative(
                                                ei, ej, Ax, Ay, Az, Bx, By, Bz,
                                                Cx, Cy, Cz, zk, n_k, lx_i, ly_i,
                                                lz_i, lx_j, ly_j, lz_j, ch_l,
                                                is_local);
                                    if (!ecp_derivative_result_is_usable(
                                            evaluated))
                                    {
                                        record_ecp_failure(
                                            failure,
                                            ecp_primitive_failure_kind(
                                                evaluated.converged),
                                            global_task_id, iat, channel, term,
                                            pi, pj, idx_i, idx_j,
                                            evaluated.value,
                                            evaluated.estimated_error);
                                        return;
                                    }
                                    for (int d = 0; d < 3; ++d)
                                    {
                                        term_A[d] =
                                            cdk * evaluated.derivative_i[d];
                                        term_B[d] =
                                            cdk * evaluated.derivative_j[d];
                                        term_A_error[d] =
                                            abs_cdk *
                                            evaluated.derivative_i_error[d];
                                        term_B_error[d] =
                                            abs_cdk *
                                            evaluated.derivative_j_error[d];
                                    }
                                }
                                else
                                {
                                    // d/dA_d = 2α V(l_i[d]+1)
                                    //           - l_i[d] V(l_i[d]-1), and
                                    // analogously for B. Each primitive error
                                    // is propagated with the absolute linear
                                    // coefficient; no cancellation is claimed.
                                    for (int d = 0; d < 3; ++d)
                                    {
                                        const int orig_i = l_i[d];
                                        const double plus_i_coefficient =
                                            2.0 * static_cast<double>(ei);
                                        const double minus_i_coefficient =
                                            static_cast<double>(orig_i);
                                        l_i[d] = orig_i + 1;
                                        const QC_ECP_PRIMITIVE_RESULT vp_i =
                                            ecp_integral_or_record_failure(
                                                ei, ej, Ax, Ay, Az, Bx, By, Bz,
                                                Cx, Cy, Cz, dist_sq, zk, n_k,
                                                l_i[0], l_i[1], l_i[2], l_j[0],
                                                l_j[1], l_j[2], ch_l, is_local,
                                                failure, global_task_id, iat,
                                                channel, term, pi, pj, idx_i,
                                                idx_j,
                                                2.0 * abs_dp * abs_cdk *
                                                    plus_i_coefficient);
                                        l_i[d] = orig_i - 1;
                                        const QC_ECP_PRIMITIVE_RESULT vm_i =
                                            ecp_integral_or_record_failure(
                                                ei, ej, Ax, Ay, Az, Bx, By, Bz,
                                                Cx, Cy, Cz, dist_sq, zk, n_k,
                                                l_i[0], l_i[1], l_i[2], l_j[0],
                                                l_j[1], l_j[2], ch_l, is_local,
                                                failure, global_task_id, iat,
                                                channel, term, pi, pj, idx_i,
                                                idx_j,
                                                2.0 * abs_dp * abs_cdk *
                                                    minus_i_coefficient);
                                        l_i[d] = orig_i;
                                        if (!ecp_primitive_result_is_usable(
                                                vp_i) ||
                                            !ecp_primitive_result_is_usable(
                                                vm_i))
                                            return;
                                        term_A[d] =
                                            cdk *
                                            (plus_i_coefficient * vp_i.value -
                                             minus_i_coefficient * vm_i.value);
                                        term_A_error[d] =
                                            qc_ecp_error_policy::
                                                Raising_Lowering_Error_Estimate(
                                                    abs_cdk,
                                                    static_cast<double>(ei),
                                                    orig_i,
                                                    vp_i.last_layer_bound,
                                                    vm_i.last_layer_bound);

                                        const int orig_j = l_j[d];
                                        const double plus_j_coefficient =
                                            2.0 * static_cast<double>(ej);
                                        const double minus_j_coefficient =
                                            static_cast<double>(orig_j);
                                        l_j[d] = orig_j + 1;
                                        const QC_ECP_PRIMITIVE_RESULT vp_j =
                                            ecp_integral_or_record_failure(
                                                ei, ej, Ax, Ay, Az, Bx, By, Bz,
                                                Cx, Cy, Cz, dist_sq, zk, n_k,
                                                l_i[0], l_i[1], l_i[2], l_j[0],
                                                l_j[1], l_j[2], ch_l, is_local,
                                                failure, global_task_id, iat,
                                                channel, term, pi, pj, idx_i,
                                                idx_j,
                                                abs_dp * abs_cdk *
                                                    plus_j_coefficient);
                                        l_j[d] = orig_j - 1;
                                        const QC_ECP_PRIMITIVE_RESULT vm_j =
                                            ecp_integral_or_record_failure(
                                                ei, ej, Ax, Ay, Az, Bx, By, Bz,
                                                Cx, Cy, Cz, dist_sq, zk, n_k,
                                                l_i[0], l_i[1], l_i[2], l_j[0],
                                                l_j[1], l_j[2], ch_l, is_local,
                                                failure, global_task_id, iat,
                                                channel, term, pi, pj, idx_i,
                                                idx_j,
                                                abs_dp * abs_cdk *
                                                    minus_j_coefficient);
                                        l_j[d] = orig_j;
                                        if (!ecp_primitive_result_is_usable(
                                                vp_j) ||
                                            !ecp_primitive_result_is_usable(
                                                vm_j))
                                            return;
                                        term_B[d] =
                                            cdk *
                                            (plus_j_coefficient * vp_j.value -
                                             minus_j_coefficient * vm_j.value);
                                        term_B_error[d] =
                                            qc_ecp_error_policy::
                                                Raising_Lowering_Error_Estimate(
                                                    abs_cdk,
                                                    static_cast<double>(ej),
                                                    orig_j,
                                                    vp_j.last_layer_bound,
                                                    vm_j.last_layer_bound);
                                    }
                                }

                                for (int d = 0; d < 3; ++d)
                                {
                                    dV_A[d] += term_A[d];
                                    dV_B[d] += term_B[d];
                                    dV_A_error[d] += term_A_error[d];
                                    dV_B_error[d] += term_B_error[d];

                                    const double bra_contribution =
                                        2.0 * dp * term_A[d];
                                    const double center_contribution =
                                        -dp * (term_A[d] + term_B[d]);
                                    const double bra_error =
                                        2.0 * abs_dp * term_A_error[d];
                                    const double center_error =
                                        abs_dp *
                                        (term_A_error[d] + term_B_error[d]);
                                    record_ecp_gradient_context(
                                        &gradient_context[atom_i * 3 + d],
                                        global_task_id, iat, channel, term, pi,
                                        pj, idx_i, idx_j, bra_contribution,
                                        bra_error);
                                    record_ecp_gradient_context(
                                        &gradient_context[iat * 3 + d],
                                        global_task_id, iat, channel, term, pi,
                                        pj, idx_i, idx_j, center_contribution,
                                        center_error);
                                }
                            };

                            if (local_ch >= 0)
                            {
                                const int t_off = ecp_channel_offsets[local_ch];
                                const int t_cnt = ecp_channel_sizes[local_ch];
                                for (int it = 0; it < t_cnt; ++it)
                                    accumulate_grad(
                                        ecp_d[t_off + it],
                                        ecp_zeta[t_off + it],
                                        ecp_n[t_off + it], -1, true, local_ch,
                                        t_off + it);
                            }

                            for (int ich = ch_start; ich < ch_end; ++ich)
                            {
                                const int ch_l = ecp_channel_l[ich];
                                if (ch_l < 0 || ch_l == l_max) continue;
                                const int t_off = ecp_channel_offsets[ich];
                                const int t_cnt = ecp_channel_sizes[ich];
                                for (int it = 0; it < t_cnt; ++it)
                                    accumulate_grad(
                                        ecp_d[t_off + it],
                                        ecp_zeta[t_off + it],
                                        ecp_n[t_off + it], ch_l, false, ich,
                                        t_off + it);
                            }

                            for (int d = 0; d < 3; ++d)
                            {
                                const int bra_component = atom_i * 3 + d;
                                const int center_component = iat * 3 + d;
                                const double bra_contribution =
                                    2.0 * dp * dV_A[d];
                                const double center_contribution =
                                    -dp * (dV_A[d] + dV_B[d]);
                                const double bra_error =
                                    2.0 * abs_dp * dV_A_error[d];
                                const double center_error =
                                    abs_dp *
                                    (dV_A_error[d] + dV_B_error[d]);

                                atomicAdd(&gradient_error[bra_component],
                                          bra_error);
                                atomicAdd(&gradient_error[center_component],
                                          center_error);
                                atomicAdd(&gradient_observable[bra_component],
                                          bra_contribution);
                                atomicAdd(
                                    &gradient_observable[center_component],
                                    center_contribution);
                            }
                        }
                    }
                }
            }
        }
    }
}

// 球谐→笛卡尔密度变换
void QC_Sph2Cart_Density_Host(int ns, int nc, const std::vector<float>& h_norms,
                              const std::vector<float>& h_C,
                              const std::vector<float>& h_M_sph,
                              std::vector<float>& h_M_cart)
{
    h_M_cart = qc_cart2sph::Transform_Weighted_Matrix_To_Cartesian(
        nc, ns, h_C, h_norms, h_M_sph);
}

// ECP 梯度驱动
bool QC_Compute_ECP_Gradient(const QC_MOLECULE& mol,
                             const QC_INTEGRAL_TASKS& task_ctx,
                             const int* d_shell_atom,
                             const float* d_P_cart_eff, double* d_grad,
                             QC_ECP_EVALUATION_FAILURE* failure)
{
    if (!mol.has_ecp || mol.ecp_total_terms == 0) return true;

    QC_ECP_DEVICE_FAILURE* device_failure = nullptr;
    if (!QC_Allocate_ECP_Failure(&device_failure))
        return QC_Report_ECP_Resource_Failure(failure);

    const std::size_t component_count =
        static_cast<std::size_t>(mol.natm) * 3;
    double* device_certificate = nullptr;
    QC_ECP_DEVICE_FAILURE* device_context = nullptr;
    if (!Device_Malloc_Safely(
            reinterpret_cast<void**>(&device_certificate),
            2 * component_count * sizeof(double)))
    {
        deviceFree(device_failure);
        return QC_Report_ECP_Resource_Failure(failure);
    }
    if (!Device_Malloc_Safely(reinterpret_cast<void**>(&device_context),
                              component_count *
                                  sizeof(QC_ECP_DEVICE_FAILURE)))
    {
        deviceFree(device_certificate);
        deviceFree(device_failure);
        return QC_Report_ECP_Resource_Failure(failure);
    }
    deviceMemset(device_certificate, 0,
                 2 * component_count * sizeof(double));
    std::vector<QC_ECP_DEVICE_FAILURE> initial_context(component_count);
    for (QC_ECP_DEVICE_FAILURE& item : initial_context)
    {
        item.kind = QC_ECP_EVALUATION_OK;
        item.task_id = item.atom = item.channel = item.term = -1;
        item.primitive_i = item.primitive_j = -1;
        item.cartesian_i = item.cartesian_j = -1;
        item.value = item.estimated_error = 0.0;
    }
    deviceMemcpy(device_context, initial_context.data(),
                 component_count * sizeof(QC_ECP_DEVICE_FAILURE),
                 deviceMemcpyHostToDevice);
    double* const device_error = device_certificate;
    double* const device_observable = device_certificate + component_count;

    const int nao_cart = mol.nao_cart;
    const int n_total = task_ctx.topo.n_1e_tasks;
    const int chunk_size = ONE_E_BATCH_SIZE;

    for (int i = 0; i < n_total; i += chunk_size)
    {
        int current_chunk = std::min(chunk_size, n_total - i);
        const QC_ONE_E_TASK* task_ptr = task_ctx.buffers.d_1e_tasks + i;
        Launch_Device_Kernel(
            ECP_Grad_Kernel, (current_chunk + 63) / 64, 64, 0, 0,
            current_chunk, i, task_ptr, mol.d_centers, mol.d_l_list,
            mol.d_exps, mol.d_coeffs, mol.d_shell_offsets, mol.d_shell_sizes,
            mol.d_ao_offsets, mol.d_atom_coords, mol.natm, mol.d_ecp_l_max,
            mol.d_ecp_atom_channel_range, mol.d_ecp_l,
            mol.d_ecp_channel_offsets, mol.d_ecp_channel_sizes, mol.d_ecp_d,
            mol.d_ecp_zeta, mol.d_ecp_n, nao_cart, d_shell_atom,
            d_P_cart_eff, device_error, device_observable, device_context,
            device_failure);
    }
    return QC_Finalize_ECP_Gradient(
        mol, task_ctx, device_failure, device_certificate, device_context,
        component_count, d_grad, failure);
}
