#define NO_GLOBAL_CONTROLLER

#include "quantum_chemistry/integrals/ri/ri_2center.hpp"
#include "quantum_chemistry/integrals/ri/ri_3center.hpp"

#include <climits>
#include <cmath>
#include <cstdio>
#include <vector>

static bool Check_Packed_R(int order)
{
    const size_t spatial = QC_RI_R_Spatial_Size(order);
    std::vector<int> seen(spatial, 0);
    for (int degree = 0; degree <= order; degree++)
        for (int t = 0; t <= degree; t++)
            for (int u = 0; u <= degree - t; u++)
            {
                const int v = degree - t - u;
                const size_t index = QC_RI_R_Spatial_Index(t, u, v);
                if (index >= spatial || seen[index] != 0) return false;
                seen[index] = 1;
            }
    for (int value : seen)
        if (value != 1) return false;
    return true;
}
static void Attach_Host_Storage(QC_RI_INTEGRAL_WORKSPACE& workspace,
                                std::vector<float>& e,
                                std::vector<float>& r,
                                std::vector<double>& boys,
                                std::vector<double>& cart)
{
    const size_t workers = (size_t)workspace.device.n_workers;
    e.resize(workers * workspace.device.e_stride);
    r.resize(workers * workspace.device.r_stride);
    boys.resize(workers * workspace.device.boys_stride);
    cart.resize(workers * workspace.device.cart_stride);
    workspace.device.e = e.data();
    workspace.device.r = r.data();
    workspace.device.boys = boys.data();
    workspace.device.cart = cart.empty() ? NULL : cart.data();
}

int main()
{
    QC_RI_INTEGRAL_WORKSPACE e2;
    QC_RI_INTEGRAL_WORKSPACE e3;
    QC_RI_INTEGRAL_WORKSPACE g2;
    QC_RI_INTEGRAL_WORKSPACE g3;
    if (!QC_RI_Build_2Center_Workspace_Layout(10000, 6, &e2) ||
        !QC_RI_Build_3Center_Workspace_Layout(10000, 6, 4, &e3) ||
        !QC_RI_Build_2Center_Gradient_Workspace_Layout(10000, 6, &g2) ||
        !QC_RI_Build_3Center_Gradient_Workspace_Layout(10000, 6, 4, &g3))
        return 1;

    if (e2.device.e_dim_a != 7 || e2.device.e_dim_b != 1 ||
        e2.device.e_dim_n != 7 ||
        e2.device.r_stride != QC_RI_R_Workspace_Size(12) ||
        e2.device.boys_stride != 13 || e2.device.cart_stride != 0)
        return 2;
    if (e3.device.e_dim_a != 7 || e3.device.e_dim_b != 5 ||
        e3.device.e_dim_n != 9 ||
        e3.device.r_stride != QC_RI_R_Workspace_Size(14) ||
        e3.device.boys_stride != 15 ||
        e3.device.cart_stride != (size_t)28 * 15 * 15)
        return 3;
    if (g2.device.e_dim_a != 8 || g2.device.e_dim_b != 1 ||
        g2.device.e_dim_n != 8 ||
        g2.device.r_stride != QC_RI_R_Workspace_Size(13) ||
        g2.device.boys_stride != 14 ||
        g2.device.cart_stride != (size_t)28 * 28 * 3)
        return 4;
    if (g3.device.e_dim_a != 8 || g3.device.e_dim_b != 6 ||
        g3.device.e_dim_n != 11 ||
        g3.device.r_stride != QC_RI_R_Workspace_Size(15) ||
        g3.device.boys_stride != 16 ||
        g3.device.cart_stride != (size_t)28 * 15 * 15 * 6)
        return 5;
    if (e2.device.n_workers <= 0 || e3.device.n_workers <= 0 ||
        g2.device.n_workers <= 0 || g3.device.n_workers <= 0 ||
        e2.e_bytes + e2.r_bytes + e2.boys_bytes + e2.cart_bytes >
            (size_t)64 * 1024 * 1024 ||
        e3.e_bytes + e3.r_bytes + e3.boys_bytes + e3.cart_bytes >
            (size_t)64 * 1024 * 1024 ||
        g2.e_bytes + g2.r_bytes + g2.boys_bytes + g2.cart_bytes >
            (size_t)64 * 1024 * 1024 ||
        g3.e_bytes + g3.r_bytes + g3.boys_bytes + g3.cart_bytes >
            (size_t)64 * 1024 * 1024)
        return 6;
    if (!Check_Packed_R(15)) return 7;

    QC_RI_INTEGRAL_WORKSPACE overflow;
    if (QC_RI_Build_2Center_Workspace_Layout(1, INT_MAX, &overflow) ||
        QC_RI_Build_3Center_Workspace_Layout(1, INT_MAX, INT_MAX, &overflow) ||
        QC_RI_Build_2Center_Gradient_Workspace_Layout(1, INT_MAX, &overflow) ||
        QC_RI_Build_3Center_Gradient_Workspace_Layout(1, INT_MAX, INT_MAX,
                                                     &overflow))
        return 8;

    std::vector<float> E(e3.device.e_stride / 6);
    QC_RI_Compute_MD_Coeffs(E.data(), e3.device.e_dim_a, e3.device.e_dim_b,
                            e3.device.e_dim_n, 6, 0, 0.0f, 0.0f, 0.25f);
    if (!std::isfinite(E[QC_RI_E_Index(6, 0, 6, e3.device.e_dim_b,
                                      e3.device.e_dim_n)]))
        return 9;

    std::vector<double> boys(e3.device.boys_stride);
    std::vector<float> R(e3.device.r_stride);
    const float displacement[3] = {0.3f, -0.2f, 0.1f};
    QC_RI_Compute_Boys_Double(boys.data(), 0.4f, 14);
    QC_RI_Compute_R_Tensor(R.data(), boys.data(), 0.7f, displacement, 14);
    for (float value : R)
        if (!std::isfinite(value)) return 10;

    const VECTOR aux_center[1] = {{0.0f, 0.0f, 0.0f}};
    const VECTOR orb_center[1] = {{0.2f, -0.1f, 0.3f}};
    const int aux_l[1] = {6};
    const int orb_l[1] = {4};
    const float aux_exp[1] = {0.8f};
    const float orb_exp[1] = {0.7f};
    const float coefficient[1] = {1.0f};
    const int shell_offset[1] = {0};
    const int shell_size[1] = {1};
    const int ao_offset[1] = {0};

    QC_RI_INTEGRAL_WORKSPACE kernel_e2;
    if (!QC_RI_Build_2Center_Workspace_Layout(1, 6, &kernel_e2)) return 11;
    std::vector<float> kernel_e2_e;
    std::vector<float> kernel_e2_r;
    std::vector<double> kernel_e2_boys;
    std::vector<double> kernel_e2_cart;
    Attach_Host_Storage(kernel_e2, kernel_e2_e, kernel_e2_r, kernel_e2_boys,
                        kernel_e2_cart);
    const QC_ONE_E_TASK task_2c[1] = {{0, 0}};
    std::vector<double> metric((size_t)28 * 28, NAN);
    QC_RI_2Center_Kernel(1, task_2c, aux_center, aux_l, aux_exp,
                         coefficient, shell_offset, shell_size, ao_offset, 28,
                         kernel_e2.device, metric.data());
    bool metric_nonzero = false;
    for (double value : metric)
    {
        if (!std::isfinite(value)) return 12;
        metric_nonzero = metric_nonzero || value != 0.0;
    }
    if (!metric_nonzero) return 13;

    QC_RI_INTEGRAL_WORKSPACE kernel_e3;
    if (!QC_RI_Build_3Center_Workspace_Layout(1, 6, 4, &kernel_e3)) return 14;
    std::vector<float> kernel_e3_e;
    std::vector<float> kernel_e3_r;
    std::vector<double> kernel_e3_boys;
    std::vector<double> kernel_e3_cart;
    Attach_Host_Storage(kernel_e3, kernel_e3_e, kernel_e3_r, kernel_e3_boys,
                        kernel_e3_cart);
    const QC_RI_3C_TASK task_3c[1] = {{0, 0, 0}};
    std::vector<double> eri3c((size_t)28 * 15 * 15, NAN);
    QC_RI_3Center_Kernel(
        1, task_3c, aux_center, aux_l, aux_exp, coefficient, shell_offset,
        shell_size, ao_offset, orb_center, orb_l, orb_exp, coefficient,
        shell_offset, shell_size, ao_offset, 28, 15, 15, 0, 0, false,
        kernel_e3.device, eri3c.data());
    bool eri3c_nonzero = false;
    for (double value : eri3c)
    {
        if (!std::isfinite(value)) return 15;
        eri3c_nonzero = eri3c_nonzero || value != 0.0;
    }
    if (!eri3c_nonzero) return 16;

    std::printf("energy_workers=%d/%d gradient_workers=%d/%d\n",
                e2.device.n_workers, e3.device.n_workers,
                g2.device.n_workers, g3.device.n_workers);
    return 0;
}
