#include "basis/basis.h"
#include "integrals/eri/common/direct_fock_kernels.hpp"
#include "integrals/ri/ri_2center.hpp"
#include "integrals/ri/ri_3center.hpp"
#include "integrals/ri/ri_metric.hpp"
#include "quantum_chemistry.h"
#include "scf/eigensolver_policy.hpp"
#include "structure/analytical_norms.hpp"

static __global__ void QC_Scale_RI_Metric_Kernel(const int naux,
                                                 const float* aux_norms,
                                                 double* metric)
{
    const int total = naux * naux;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        const int P = idx / naux;
        const int Q = idx - P * naux;
        metric[idx] *= (double)aux_norms[P] * (double)aux_norms[Q];
    }
}

static __global__ void QC_Scale_RI_3Center_Kernel(const int naux, const int nao,
                                                  const float* aux_norms,
                                                  const float* orb_norms,
                                                  double* eri3c)
{
    const long long total = (long long)naux * nao * nao;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        const int nu = (int)(idx % nao);
        const long long tmp = idx / nao;
        const int mu = (int)(tmp % nao);
        const int P = (int)(tmp / nao);
        eri3c[idx] *= (double)aux_norms[P] * (double)orb_norms[mu] *
                      (double)orb_norms[nu];
    }
}

static void QC_Build_RI_Aux_Norms(QUANTUM_CHEMISTRY* qc, CONTROLLER* controller)
{
    auto& ri = qc->scf_ws.ri;
    try
    {
        ri.h_aux_norms = qc_analytical_norms::Build(
            ri.h_aux_l_list, ri.h_aux_shell_sizes, ri.h_aux_shell_offsets,
            ri.h_aux_exps, ri.h_aux_coeffs, ri.h_aux_ao_offsets,
            ri.h_aux_ao_offsets_sph, true, ri.h_U_aux, ri.naux_cart, ri.naux);
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "QUANTUM_CHEMISTRY::RI_Precompute",
            "Reason:\n    Failed to build finite positive auxiliary-basis "
            "normalizations: %s\n",
            error.what());
    }

    if (ri.d_aux_norms) deviceFree(ri.d_aux_norms);
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_norms,
                                  (void*)ri.h_aux_norms.data(),
                                  sizeof(float) * ri.h_aux_norms.size());
}

void QUANTUM_CHEMISTRY::Initial_Auxiliary_Basis(CONTROLLER* controller)
{
    QC_BASIS_SET* aux_basis = QC_Get_JKFIT_Basis(orbital_basis_name.c_str());
    if (aux_basis == nullptr)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, "QUANTUM_CHEMISTRY::Initial",
            "Reason:\n    No JKFIT auxiliary basis available for orbital "
            "basis \"%s\". RI requires a matching auxiliary basis.\n",
            orbital_basis_name.c_str());
    }
    aux_basis->Initialize();

    auto& ri = scf_ws.ri;
    ri.naux_cart = 0;
    ri.naux = 0;
    ri.naux_bas = 0;

    // Auxiliary-basis identity follows the real element, never the ECP
    // effective charge.  The latter is still stored in atm for integral
    // kernels that inspect the nuclear charge field.
    for (int i = 0; i < mol.natm; ++i)
    {
        const int atomic_number = mol.Atomic_Number(i);
        const int effective_nuclear_charge =
            mol.Effective_Nuclear_Charge(i);
        const char* symbol =
            sponge_qc_elements::Symbol_From_Atomic_Number(atomic_number);
        if (symbol == nullptr)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "QUANTUM_CHEMISTRY::Initial",
                "Reason:\n    Unknown atomic number %d in auxiliary basis "
                "setup\n",
                atomic_number);
        }
        const std::string sym(symbol);

        auto it_basis = aux_basis->data.find(sym);
        if (it_basis == aux_basis->data.end())
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand, "QUANTUM_CHEMISTRY::Initial",
                "Reason:\n    JKFIT auxiliary basis not available for element "
                "%s (Z=%d)\n",
                sym.c_str(), atomic_number);
        }
        const auto& shells = it_basis->second;

        // 辅助基的 atm 条目
        int ptr_coord = ri.h_aux_env.size();
        ri.h_aux_env.push_back(0.0f);  // 坐标占位，后续 Update_Coordinates 填充
        ri.h_aux_env.push_back(0.0f);
        ri.h_aux_env.push_back(0.0f);
        ri.h_aux_atm.push_back(effective_nuclear_charge);
        ri.h_aux_atm.push_back(ptr_coord);
        ri.h_aux_atm.push_back(1);
        ri.h_aux_atm.push_back(0);
        ri.h_aux_atm.push_back(0);
        ri.h_aux_atm.push_back(0);

        for (const auto& shell : shells)
        {
            int cart_dimension = 0;
            int spherical_dimension = 0;
            try
            {
                const std::pair<int, int> dimensions =
                    qc_cart2sph::Int_Dimensions(shell.l);
                cart_dimension = dimensions.first;
                spherical_dimension = dimensions.second;
            }
            catch (const std::exception& error)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow, "QUANTUM_CHEMISTRY::Initial",
                    "Reason:\n    Invalid auxiliary-basis angular momentum "
                    "l=%d for element %s: %s\n",
                    shell.l, sym.c_str(), error.what());
            }
            if (shell.exps.empty() ||
                shell.exps.size() != shell.coeffs.size() ||
                shell.exps.size() > static_cast<std::size_t>(
                                        std::numeric_limits<int>::max() / 2) ||
                ri.h_aux_env.size() >
                    static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                        2U * shell.exps.size() ||
                ri.h_aux_exps.size() >
                    static_cast<std::size_t>(std::numeric_limits<int>::max()) -
                        shell.exps.size())
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow, "QUANTUM_CHEMISTRY::Initial",
                    "Reason:\n    Invalid or overflowing auxiliary primitive "
                    "data for l=%d shell of element %s\n",
                    shell.l, sym.c_str());
            }
            if (ri.naux_cart >
                    std::numeric_limits<int>::max() - cart_dimension ||
                ri.naux >
                    std::numeric_limits<int>::max() - spherical_dimension ||
                ri.naux_bas == std::numeric_limits<int>::max())
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow, "QUANTUM_CHEMISTRY::Initial",
                    "Reason:\n    Auxiliary-basis dimension overflows int "
                    "while adding l=%d shell for element %s\n",
                    shell.l, sym.c_str());
            }
            int ptr_exp = ri.h_aux_env.size();
            ri.h_aux_env.insert(ri.h_aux_env.end(), shell.exps.begin(),
                                shell.exps.end());
            int ptr_coeff = ri.h_aux_env.size();
            ri.h_aux_env.insert(ri.h_aux_env.end(), shell.coeffs.begin(),
                                shell.coeffs.end());

            ri.h_aux_bas.push_back(i);  // atom index
            ri.h_aux_bas.push_back(shell.l);
            ri.h_aux_bas.push_back(shell.exps.size());
            ri.h_aux_bas.push_back(1);  // nctr
            ri.h_aux_bas.push_back(0);
            ri.h_aux_bas.push_back(ptr_exp);
            ri.h_aux_bas.push_back(ptr_coeff);
            ri.h_aux_bas.push_back(0);

            ri.h_aux_l_list.push_back(shell.l);
            ri.h_aux_shell_sizes.push_back(shell.exps.size());
            ri.h_aux_shell_offsets.push_back(ri.h_aux_exps.size());

            ri.h_aux_exps.insert(ri.h_aux_exps.end(), shell.exps.begin(),
                                 shell.exps.end());
            ri.h_aux_coeffs.insert(ri.h_aux_coeffs.end(), shell.coeffs.begin(),
                                   shell.coeffs.end());
            ri.h_aux_centers.push_back(VECTOR(0.0f));

            ri.h_aux_ao_offsets.push_back(ri.naux_cart);
            ri.h_aux_ao_offsets_sph.push_back(ri.naux);

            ri.naux_cart += cart_dimension;
            ri.naux += spherical_dimension;
            ri.naux_bas++;
        }
    }

    try
    {
        ri.h_U_aux =
            QC_Build_Cart2Sph_Mat_Host(ri.h_aux_l_list, ri.naux_cart, ri.naux);
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Initial",
            "Reason:\n    Failed to construct auxiliary "
            "Cartesian-to-spherical matrix: %s\n",
            error.what());
    }

    // 拷贝到 device
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_l_list,
                                  (void*)ri.h_aux_l_list.data(),
                                  sizeof(int) * ri.h_aux_l_list.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_shell_offsets,
                                  (void*)ri.h_aux_shell_offsets.data(),
                                  sizeof(int) * ri.h_aux_shell_offsets.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_shell_sizes,
                                  (void*)ri.h_aux_shell_sizes.data(),
                                  sizeof(int) * ri.h_aux_shell_sizes.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_ao_offsets,
                                  (void*)ri.h_aux_ao_offsets.data(),
                                  sizeof(int) * ri.h_aux_ao_offsets.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_ao_offsets_sph,
                                  (void*)ri.h_aux_ao_offsets_sph.data(),
                                  sizeof(int) * ri.h_aux_ao_offsets_sph.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_exps,
                                  (void*)ri.h_aux_exps.data(),
                                  sizeof(float) * ri.h_aux_exps.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_coeffs,
                                  (void*)ri.h_aux_coeffs.data(),
                                  sizeof(float) * ri.h_aux_coeffs.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_centers,
                                  (void*)ri.h_aux_centers.data(),
                                  sizeof(VECTOR) * ri.h_aux_centers.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_atm,
                                  (void*)ri.h_aux_atm.data(),
                                  sizeof(int) * ri.h_aux_atm.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_bas,
                                  (void*)ri.h_aux_bas.data(),
                                  sizeof(int) * ri.h_aux_bas.size());
    Device_Malloc_And_Copy_Safely((void**)&ri.d_aux_env,
                                  (void*)ri.h_aux_env.data(),
                                  sizeof(float) * ri.h_aux_env.size());
}

void QUANTUM_CHEMISTRY::RI_Memory_Allocate()
{
    auto& ri = scf_ws.ri;
    const int naux = ri.naux;
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    const long long naux2 = (long long)naux * naux;
    const long long n3c = (long long)naux * nao2;

    // 根据用户设置或自动选择 stored/direct 模式
    const double mem_3c_mb = n3c * sizeof(double) / 1e6;
    if (ri.mode == QC_RI_WORKSPACE::DF_DIRECT)
    {
        ri.direct = true;
    }
    else if (ri.mode == QC_RI_WORKSPACE::DF_STORED)
    {
        ri.direct = false;
    }
    else  // DF_AUTO
    {
        ri.direct = (mem_3c_mb > 512.0);
    }

    // 二中心 metric（两种模式均需要）
    Device_Malloc_Safely((void**)&ri.d_metric, sizeof(double) * naux2);
    Device_Malloc_Safely((void**)&ri.d_metric_inv, sizeof(double) * naux2);
    Device_Malloc_Safely((void**)&ri.d_metric_inv_sqrt, sizeof(double) * naux2);
    deviceMemset(ri.d_metric, 0, sizeof(double) * naux2);

    // RI-J scratch
    Device_Malloc_Safely((void**)&ri.d_d_vec, sizeof(double) * naux);
    Device_Malloc_Safely((void**)&ri.d_g_vec, sizeof(double) * naux);

    // RI-J persistent scratch
    Device_Malloc_Safely((void**)&ri.d_P_double, sizeof(double) * nao2);
    Device_Malloc_Safely((void**)&ri.d_J_double, sizeof(double) * nao2);

    // RI-K scratch.  An ensemble density can have rank greater than the
    // integer occupied-orbital count, so dimension these buffers for the
    // mathematically complete nao-column density factor.
    if (dft.exx_fraction != 0.0f)
    {
        size_t n_bocc = 0;
        size_t n_bocc_bytes = 0;
        if (!QC_RI_Checked_Mul_Size((size_t)naux, (size_t)nao, &n_bocc) ||
            !QC_RI_Checked_Mul_Size(n_bocc, (size_t)nao, &n_bocc) ||
            n_bocc > (size_t)std::numeric_limits<int>::max() ||
            !QC_RI_Checked_Bytes(n_bocc, sizeof(float), &n_bocc_bytes))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::RI_Memory_Allocate",
                "Reason:\n    RI-K density-factor workspace exceeds the "
                "supported index range: naux=%d, nao=%d\n",
                naux, nao);
            return;
        }
        Device_Malloc_Safely((void**)&ri.d_density_factor_alpha,
                             sizeof(float) * nao2);
        if (scf_ws.runtime.unrestricted)
            Device_Malloc_Safely((void**)&ri.d_density_factor_beta,
                                 sizeof(float) * nao2);
        Device_Malloc_Safely((void**)&ri.d_density_factor_accum,
                             sizeof(double));
        Device_Malloc_Safely((void**)&ri.d_density_factor_failure,
                             sizeof(int));
        Device_Malloc_Safely((void**)&ri.d_B_occ, n_bocc_bytes);
        Device_Malloc_Safely((void**)&ri.d_B_flat, n_bocc_bytes);
        Device_Malloc_Safely((void**)&ri.d_K_scratch, sizeof(float) * nao2);
    }

    // cart2sph 轨道变换矩阵（stored 和 direct 模式都需要）
    if (mol.is_spherical)
    {
        try
        {
            ri.h_U_orb =
                QC_Build_Cart2Sph_Mat_Host(mol.h_l_list, mol.nao_cart, mol.nao);
        }
        catch (const std::exception& error)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::Initial",
                "Reason:\n    Failed to construct RI orbital "
                "Cartesian-to-spherical matrix: %s\n",
                error.what());
        }
    }
    else
    {
        ri.h_U_orb.clear();
    }

    if (!ri.direct)
    {
        // Stored 模式：预存 eri3c 和 B
        Device_Malloc_Safely((void**)&ri.d_eri3c, sizeof(double) * n3c);
        deviceMemset(ri.d_eri3c, 0, sizeof(double) * n3c);
        Device_Malloc_Safely((void**)&ri.d_B, sizeof(float) * n3c);
    }
    // Direct 模式：d_3c_buf 在 Build_Fock_RI_Direct 内临时分配/释放
}

void QUANTUM_CHEMISTRY::RI_Precompute()
{
    auto& ri = scf_ws.ri;
    const int naux = ri.naux;
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    const int threads = 256;

    QC_Build_RI_Aux_Norms(this, controller);

    const int Pc = ri.naux_cart, Ps = ri.naux;
    const int Mc = mol.nao_cart, Ms = mol.nao;
    int max_aux_l = 0;
    for (int l : ri.h_aux_l_list) max_aux_l = std::max(max_aux_l, l);
    int max_orb_l = 0;
    for (int l : mol.h_l_list) max_orb_l = std::max(max_orb_l, l);

    // 上传 U_aux/U_orb 到 device (float → double, 用于 DGEMM cart2sph)
    double* d_U_aux = NULL;
    double* d_U_orb = NULL;
    {
        Device_Malloc_Safely((void**)&d_U_aux, sizeof(double) * Pc * Ps);
        std::vector<double> h_Ud((size_t)Pc * Ps);
        for (int i = 0; i < Pc * Ps; i++) h_Ud[i] = (double)ri.h_U_aux[i];
        deviceMemcpy(d_U_aux, h_Ud.data(), sizeof(double) * Pc * Ps,
                     deviceMemcpyHostToDevice);
        if (mol.is_spherical && Mc > 0 && Ms > 0)
        {
            Device_Malloc_Safely((void**)&d_U_orb, sizeof(double) * Mc * Ms);
            h_Ud.resize((size_t)Mc * Ms);
            for (int i = 0; i < Mc * Ms; i++) h_Ud[i] = (double)ri.h_U_orb[i];
            deviceMemcpy(d_U_orb, h_Ud.data(), sizeof(double) * Mc * Ms,
                         deviceMemcpyHostToDevice);
        }
    }

    // 1. 计算二中心 metric (P|Q)
    {
        const long long n2c_cart = (long long)Pc * Pc;
        double* d_metric_cart = NULL;
        QC_ONE_E_TASK* d_2c_tasks = NULL;
        Device_Malloc_Safely((void**)&d_metric_cart, sizeof(double) * n2c_cart);
        deviceMemset(d_metric_cart, 0, sizeof(double) * n2c_cart);

        std::vector<QC_ONE_E_TASK> h_2c_tasks;
        for (int i = 0; i < ri.naux_bas; i++)
            for (int j = 0; j <= i; j++) h_2c_tasks.push_back({i, j});
        if (h_2c_tasks.size() > (size_t)std::numeric_limits<int>::max())
        {
            deviceFree(d_metric_cart);
            deviceFree(d_U_aux);
            deviceFree(d_U_orb);
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::RI_Precompute",
                "Reason:\n    RI two-center task count exceeds int range\n");
            return;
        }
        const int n_2c = (int)h_2c_tasks.size();

        QC_RI_INTEGRAL_WORKSPACE metric_workspace;
        if (!QC_RI_Build_2Center_Workspace_Layout(n_2c, max_aux_l,
                                                  &metric_workspace))
        {
            deviceFree(d_metric_cart);
            deviceFree(d_U_aux);
            deviceFree(d_U_orb);
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::RI_Precompute",
                "Reason:\n    RI two-center workspace dimensions overflow "
                "for max auxiliary angular momentum %d and %d tasks\n",
                max_aux_l, n_2c);
            return;
        }
        if (!QC_RI_Allocate_Integral_Workspace(&metric_workspace))
        {
            deviceFree(d_metric_cart);
            deviceFree(d_U_aux);
            deviceFree(d_U_orb);
            return;
        }

        Device_Malloc_Safely((void**)&d_2c_tasks, sizeof(QC_ONE_E_TASK) * n_2c);
        deviceMemcpy(d_2c_tasks, h_2c_tasks.data(),
                     sizeof(QC_ONE_E_TASK) * n_2c, deviceMemcpyHostToDevice);

        QC_Launch_RI_2Center_Kernel(
            threads, n_2c, d_2c_tasks, ri.d_aux_centers, ri.d_aux_l_list,
            ri.d_aux_exps, ri.d_aux_coeffs, ri.d_aux_shell_offsets,
            ri.d_aux_shell_sizes, ri.d_aux_ao_offsets, ri.naux_cart,
            metric_workspace, d_metric_cart);

        QC_RI_Free_Integral_Workspace(&metric_workspace);
        deviceFree(d_2c_tasks);

        // M_sph = U^T @ M_cart @ U  (device DGEMM)
        // d_U_aux 列优先 = U^T[Ps,Pc], lda=Ps; op=N → U^T, op=T → U
        // d_metric_cart 列优先 = M[Pc,Pc] (对称), lda=Pc
        const double one = 1.0, zero = 0.0;
        double* d_tmp = NULL;
        Device_Malloc_Safely((void**)&d_tmp, sizeof(double) * Ps * Pc);
        // Step A: tmp[Ps,Pc] = U^T @ M_cart
        deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, Ps, Pc,
                        Pc, &one, d_U_aux, Ps, d_metric_cart, Pc, &zero, d_tmp,
                        Ps);
        // Step B: M_sph[Ps,Ps] = tmp @ U = tmp @ (U^T)^T
        deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_T, Ps, Ps,
                        Pc, &one, d_tmp, Ps, d_U_aux, Ps, &zero, ri.d_metric,
                        Ps);
        deviceFree(d_tmp);
        deviceFree(d_metric_cart);
        Launch_Device_Kernel(QC_Scale_RI_Metric_Kernel,
                             ((long long)naux * naux + threads - 1) / threads,
                             threads, 0, 0, naux, ri.d_aux_norms, ri.d_metric);
        // 2c metric done
    }

    // The Coulomb metric must be finite everywhere and strictly positive on
    // its diagonal in both stored and direct modes.  Validate before any
    // screening or eigensystem construction can disguise invalid data.
    {
        std::vector<double> h_metric_validation((size_t)Ps * Ps);
        deviceMemcpy(h_metric_validation.data(), ri.d_metric,
                     sizeof(double) * h_metric_validation.size(),
                     deviceMemcpyDeviceToHost);
        for (size_t index = 0; index < h_metric_validation.size(); ++index)
        {
            const double value = h_metric_validation[index];
            if (!Double_Memory_Is_Finite(&value))
            {
                deviceFree(d_U_aux);
                deviceFree(d_U_orb);
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::RI_Precompute",
                    "Reason:\n    RI Coulomb metric contains a non-finite "
                    "value at flat index %zu: %.17g\n",
                    index, value);
                return;
            }
        }
        for (int p = 0; p < Ps; ++p)
        {
            const double diagonal = h_metric_validation[(size_t)p * Ps + p];
            if (!(diagonal > 0.0))
            {
                deviceFree(d_U_aux);
                deviceFree(d_U_orb);
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::RI_Precompute",
                    "Reason:\n    RI Coulomb metric diagonal %d/%d must be "
                    "strictly positive, got %.17g\n",
                    p, Ps, diagonal);
                return;
            }
        }
    }

    // ---- 2. 计算三中心积分 (P|μν)（仅 stored 模式）----
    if (!ri.direct)
    {
        size_t n3c_cart = 0;
        size_t n3c_cart_bytes = 0;
        if (!QC_RI_Checked_Mul_Size((size_t)ri.naux_cart, (size_t)mol.nao_cart,
                                    &n3c_cart) ||
            !QC_RI_Checked_Mul_Size(n3c_cart, (size_t)mol.nao_cart,
                                    &n3c_cart) ||
            !QC_RI_Checked_Bytes(n3c_cart, sizeof(double), &n3c_cart_bytes))
        {
            deviceFree(d_U_aux);
            deviceFree(d_U_orb);
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::RI_Precompute",
                "Reason:\n    RI stored three-center Cartesian tensor size "
                "overflows size_t\n");
            return;
        }
        double* d_eri3c_cart = NULL;
        QC_RI_3C_TASK* d_3c_tasks = NULL;
        Device_Malloc_Safely((void**)&d_eri3c_cart, n3c_cart_bytes);
        deviceMemset(d_eri3c_cart, 0, n3c_cart_bytes);

        // Schwarz 筛选: |(P|μν)| ≤ sqrt(P|P) × sqrt(μν|μν)
        // 辅助基 bound: 从 2c metric 对角元素取 sqrt（归一化后）
        std::vector<double> h_metric_diag(Ps);
        {
            std::vector<double> h_metric_full((size_t)Ps * Ps);
            deviceMemcpy(h_metric_full.data(), ri.d_metric,
                         sizeof(double) * Ps * Ps,
                         deviceMemcpyDeviceToHost);
            for (int p = 0; p < Ps; p++)
                h_metric_diag[p] =
                    sqrt(fabs(h_metric_full[(size_t)p * Ps + p]));
        }
        std::vector<float> aux_shell_bound(ri.naux_bas, 0.0f);
        for (int sh = 0; sh < ri.naux_bas; sh++)
        {
            const int off = mol.is_spherical
                                ? ri.h_aux_ao_offsets_sph[sh]
                                : ri.h_aux_ao_offsets[sh];
            const int dim = mol.is_spherical
                                ? (2 * ri.h_aux_l_list[sh] + 1)
                                : ((ri.h_aux_l_list[sh] + 1) *
                                   (ri.h_aux_l_list[sh] + 2) / 2);
            for (int i = 0; i < dim; i++)
                aux_shell_bound[sh] = std::max(
                    aux_shell_bound[sh], (float)h_metric_diag[off + i]);
        }
        std::vector<std::vector<float>> orb_pair_bound(
            mol.nbas, std::vector<float>(mol.nbas, 0.0f));
        for (int pid = 0; pid < task_ctx.topo.n_shell_pairs; pid++)
        {
            const auto& pair = task_ctx.topo.h_shell_pairs[pid];
            const float bound = task_ctx.topo.h_shell_pair_bounds[pid];
            orb_pair_bound[pair.x][pair.y] = bound;
            orb_pair_bound[pair.y][pair.x] = bound;
        }

        const float screen_tol = task_ctx.params.eri_shell_screen_tol;
        std::vector<QC_RI_3C_TASK> h_3c_tasks;
        for (int P = 0; P < ri.naux_bas; P++)
            for (int mu = 0; mu < mol.nbas; mu++)
                for (int nu = 0; nu <= mu; nu++)
                {
                    if (aux_shell_bound[P] * orb_pair_bound[mu][nu] <
                        screen_tol)
                        continue;
                    h_3c_tasks.push_back({P, mu, nu});
                }
        if (h_3c_tasks.size() > (size_t)std::numeric_limits<int>::max())
        {
            deviceFree(d_eri3c_cart);
            deviceFree(d_U_aux);
            deviceFree(d_U_orb);
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::RI_Precompute",
                "Reason:\n    RI three-center task count exceeds int range\n");
            return;
        }
        const int n_3c = (int)h_3c_tasks.size();

        QC_RI_INTEGRAL_WORKSPACE eri3c_workspace;
        if (!QC_RI_Build_3Center_Workspace_Layout(n_3c, max_aux_l, max_orb_l,
                                                  &eri3c_workspace))
        {
            deviceFree(d_eri3c_cart);
            deviceFree(d_U_aux);
            deviceFree(d_U_orb);
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::RI_Precompute",
                "Reason:\n    RI three-center workspace dimensions overflow "
                "for auxiliary/orbital angular momenta %d/%d and %d tasks\n",
                max_aux_l, max_orb_l, n_3c);
            return;
        }
        if (!QC_RI_Allocate_Integral_Workspace(&eri3c_workspace))
        {
            deviceFree(d_eri3c_cart);
            deviceFree(d_U_aux);
            deviceFree(d_U_orb);
            return;
        }

        if (n_3c > 0)
        {
            if (!Device_Malloc_Safely(
                    (void**)&d_3c_tasks,
                    sizeof(QC_RI_3C_TASK) * (size_t)n_3c))
            {
                QC_RI_Free_Integral_Workspace(&eri3c_workspace);
                deviceFree(d_eri3c_cart);
                deviceFree(d_U_aux);
                deviceFree(d_U_orb);
                return;
            }
            deviceMemcpy(d_3c_tasks, h_3c_tasks.data(),
                         sizeof(QC_RI_3C_TASK) * (size_t)n_3c,
                         deviceMemcpyHostToDevice);
        }

        QC_Launch_RI_3Center_Kernel(
            threads, n_3c, d_3c_tasks, ri.d_aux_centers, ri.d_aux_l_list,
            ri.d_aux_exps, ri.d_aux_coeffs, ri.d_aux_shell_offsets,
            ri.d_aux_shell_sizes, ri.d_aux_ao_offsets, mol.d_centers,
            mol.d_l_list, mol.d_exps, mol.d_coeffs, mol.d_shell_offsets,
            mol.d_shell_sizes, mol.d_ao_offsets, ri.naux_cart, mol.nao_cart,
            mol.nao_cart, 0, 0, true, eri3c_workspace, d_eri3c_cart);

        QC_RI_Free_Integral_Workspace(&eri3c_workspace);
        if (d_3c_tasks != NULL) deviceFree(d_3c_tasks);
        // 3c kernel done

        // Cart2Sph 变换: 全部在 device 上用 DGEMM
        // d_U_aux 列优先 = U_aux^T[Ps,Pc], lda=Ps; op=N → U^T, op=T → U
        // d_U_orb 列优先 = U_orb^T[Ms,Mc], lda=Ms; op=N → U_orb^T, op=T → U_orb
        // d_eri3c_cart 行优先 [Pc,Mc,Mc] = 列优先视为 T^T 各种切片
        const double one = 1.0, zero = 0.0;

        if (!mol.is_spherical)
        {
            // 仅 P 变换: T_sph[Ps,Mc²] = U_aux^T @ T_cart[Pc,Mc²]
            // 列优先: T_sph^T[Mc²,Ps] = T_cart^T[Mc²,Pc] @ U[Pc,Ps]
            deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_T,
                            Mc * Mc, Ps, Pc, &one, d_eri3c_cart, Mc * Mc,
                            d_U_aux, Ps, &zero, ri.d_eri3c, Mc * Mc);
        }
        else
        {
            // Step 1 (ν): T1[Pc*Mc, Ms] = T_cart[Pc*Mc, Mc] @ U_orb[Mc, Ms]
            // 列优先: T1^T[Ms, Pc*Mc] = U_orb^T[Ms,Mc] @ T_cart^T[Mc, Pc*Mc]
            const long long n_step1 = (long long)Pc * Mc * Ms;
            double* d_step1 = NULL;
            Device_Malloc_Safely((void**)&d_step1, sizeof(double) * n_step1);
            deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, Ms,
                            Pc * Mc, Mc, &one, d_U_orb, Ms, d_eri3c_cart, Mc,
                            &zero, d_step1, Ms);

            // Step 2 (μ): 对每个 P, T2_P[Ms,Ms] = U_orb^T @ T1_P (行优先乘法)
            // 列优先: T2_P^T[Ms,Ms] = T1_P^T[Ms,Mc] @ U_orb[Mc,Ms]
            const long long n_step2 = (long long)Pc * Ms * Ms;
            double* d_step2 = NULL;
            Device_Malloc_Safely((void**)&d_step2, sizeof(double) * n_step2);
            for (int P = 0; P < Pc; P++)
            {
                deviceBlasDgemm(
                    blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_T, Ms, Ms, Mc,
                    &one, d_step1 + (long long)P * Mc * Ms, Ms, d_U_orb, Ms,
                    &zero, d_step2 + (long long)P * Ms * Ms, Ms);
            }
            deviceFree(d_step1);

            // Step 3 (P): T_sph[Ps, Ms²] = U_aux^T @ T2[Pc, Ms²]
            // 列优先: T_sph^T[Ms²,Ps] = T2^T[Ms²,Pc] @ U[Pc,Ps]
            deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_T,
                            Ms * Ms, Ps, Pc, &one, d_step2, Ms * Ms, d_U_aux,
                            Ps, &zero, ri.d_eri3c, Ms * Ms);
            deviceFree(d_step2);
        }
        deviceFree(d_eri3c_cart);

        Launch_Device_Kernel(
            QC_Scale_RI_3Center_Kernel,
            ((long long)naux * nao * nao + threads - 1) / threads, threads, 0,
            0, naux, nao, ri.d_aux_norms, scf_ws.ortho.d_norms, ri.d_eri3c);
    }

    // 3. One validated eigensystem defines both truncated metric inverses.
    bool metric_inverses_ok = false;
    int metric_solver_api_status = 0;
    int metric_solver_info = 0;
    try
    {
        metric_inverses_ok = QC_RI_Build_Metric_Inverses(
            solver_handle, blas_handle, naux, ri.d_metric, ri.d_metric_inv_sqrt,
            ri.d_metric_inv, &ri.naux_eff, &ri.h_eigval, &ri.h_eigvec,
            &ri.h_metric_inv_sqrt, &ri.h_metric_inv,
            &metric_solver_api_status, &metric_solver_info);
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "QUANTUM_CHEMISTRY::RI_Precompute",
            "Reason:\n    RI Coulomb metric inverse construction threw an "
            "exception: %s\n",
            error.what());
        return;
    }
    if (metric_solver_api_status != 0 || metric_solver_info != 0)
    {
        QC_SCF_Require_Eigensolver_Success(
            QC_SCF_EIGENSOLVER_RI_LOEWNER,
            QC_SCF_EIGENSOLVER_CHANNEL_AUXILIARY, naux,
            metric_solver_api_status, metric_solver_info,
            [&](const QC_SCF_Eigensolver_Failure& failure)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::RI_Precompute",
                    "Reason:\n    eigensolver failed during %s for channel "
                    "%s: dimension=%d, api_status=%d, info=%d\n",
                    failure.stage_name, failure.channel_name,
                    failure.dimension, failure.api_status, failure.info);
            });
        return;
    }
    if (!metric_inverses_ok || ri.naux_eff <= 0 || ri.naux_eff > naux ||
        ri.h_eigval.size() != (size_t)naux ||
        ri.h_eigvec.size() != (size_t)naux * naux ||
        ri.h_metric_inv_sqrt.size() != (size_t)naux * naux ||
        ri.h_metric_inv.size() != (size_t)naux * naux)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "QUANTUM_CHEMISTRY::RI_Precompute",
            "Reason:\n    RI Coulomb metric eigendecomposition failed or "
            "retained an invalid dimension: naux=%d, naux_eff=%d\n",
            naux, ri.naux_eff);
        return;
    }
    const int n_skip = naux - ri.naux_eff;
    for (int i = 0; i < naux; i++)
    {
        const double value = ri.h_eigval[i];
        if (!Double_Memory_Is_Finite(&value) || (i >= n_skip && !(value > 0.0)))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::RI_Precompute",
                "Reason:\n    invalid RI Coulomb metric eigenvalue %d/%d: "
                "%.17g (retained=%d)\n",
                i, naux, value, i >= n_skip);
            return;
        }
    }
    for (size_t i = 0; i < ri.h_eigvec.size(); i++)
    {
        const double value = ri.h_eigvec[i];
        if (!Double_Memory_Is_Finite(&value))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::RI_Precompute",
                "Reason:\n    non-finite RI Coulomb metric eigenvector "
                "coefficient at flat index %zu: %.17g\n",
                i, value);
            return;
        }
    }

    // ---- 4. 构建 B 张量（仅 stored 模式）----
    if (!ri.direct)
    {
        double* d_B_double = NULL;
        const long long n3c = (long long)naux * nao2;
        Device_Malloc_Safely((void**)&d_B_double, sizeof(double) * n3c);

        // eri3c is stored with row-major indexing [P, mu*nao+nu], which is
        // equivalent to a col-major matrix [nao2 x naux]. Build B in the same
        // physical layout so downstream kernels can keep using row-major
        // indexing on [naux x nao2].
        const double one = 1.0, zero = 0.0;
        deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, nao2,
                        naux, naux, &one, ri.d_eri3c, nao2,
                        ri.d_metric_inv_sqrt, naux, &zero, d_B_double, nao2);
        QC_Double_To_Float((int)n3c, d_B_double, ri.d_B);
        deviceFree(d_B_double);
    }

    deviceFree(d_U_aux);
    deviceFree(d_U_orb);
}

// F[i] += scale * K[i]
static __global__ void QC_Scaled_Add_Kernel(const int n, const float scale,
                                            const float* src, float* dst)
{
    SIMPLE_DEVICE_FOR(idx, n) { dst[idx] += scale * src[idx]; }
}

#ifdef GPU_ARCH_NAME
static __device__ __forceinline__ bool QC_RI_Float_Is_Finite(float value)
{
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
}
#else
static __host__ __device__ __forceinline__ bool QC_RI_Float_Is_Finite(
    float value)
{
    unsigned int bits = 0;
    static_assert(
        sizeof(bits) == sizeof(value) && std::numeric_limits<float>::is_iec559,
        "SPONGE requires 32-bit IEEE-754 floats");
    memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return (bits & 0x7f800000U) != 0x7f800000U;
}
#endif

// Copy the current AO density into the symmetric double matrix consumed by
// the eigensolver.  The antisymmetric norm is accumulated independently: a
// silently symmetrized, corrupted density must never become an RI-K input.
static __global__ void QC_Prepare_RI_Density_Factor_Kernel(
    const int nao, const double density_scale, const float* density,
    double* symmetric_density, double* antisymmetric_squared_norm,
    int* nonfinite_index)
{
    const int total = nao * nao;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        const int mu = idx / nao;
        const int nu = idx - mu * nao;
        const float p_mn = density[idx];
        const float p_nm = density[nu * nao + mu];
        if (!QC_RI_Float_Is_Finite(p_mn)) atomicExch(nonfinite_index, idx);
        const double difference =
            density_scale * ((double)p_mn - (double)p_nm);
        atomicAdd(antisymmetric_squared_norm, difference * difference);
        symmetric_density[idx] =
            0.5 * density_scale * ((double)p_mn + (double)p_nm);
    }
}

// dsyevd returns column-major eigenvectors.  Pack all strictly positive
// eigenmodes into the row-major factor expected by the existing RI-K GEMMs.
static __global__ void QC_Build_RI_Density_Factor_Kernel(
    const int nao, const int first_positive, const int rank,
    const double* eigenvectors, const double* eigenvalues, float* factor)
{
    const int total = nao * nao;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        const int mu = idx / nao;
        const int column = idx - mu * nao;
        if (column < rank)
        {
            const int eigenvector_column = first_positive + column;
            factor[idx] =
                (float)(eigenvectors[mu + eigenvector_column * nao] *
                        sqrt(eigenvalues[eigenvector_column]));
        }
        else
        {
            factor[idx] = 0.0f;
        }
    }
}

bool QUANTUM_CHEMISTRY::Factor_RI_Spin_Density(
    const float* d_density, double density_scale, float* d_factor,
    int* factor_rank, QC_SCF_Eigensolver_Channel spin_channel)
{
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    auto& ri = scf_ws.ri;
    const char* channel_name = QC_SCF_Eigensolver_Channel_Name(spin_channel);
    if (d_density == nullptr || d_factor == nullptr || factor_rank == nullptr ||
        ri.d_density_factor_accum == nullptr ||
        ri.d_density_factor_failure == nullptr || nao <= 0 ||
        !Double_Memory_Is_Finite(&density_scale) || !(density_scale > 0.0))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
            "Reason:\n    incomplete or invalid RI-K density-factor input "
            "for channel %s: nao=%d, density_scale=%.17g\n",
            channel_name, nao, density_scale);
        return false;
    }

    deviceMemset(ri.d_density_factor_accum, 0, sizeof(double));
    deviceMemset(ri.d_density_factor_failure, -1, sizeof(int));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_Prepare_RI_Density_Factor_Kernel,
        Positive_Int_Ceil_Div(nao2, threads), threads, 0, 0, nao,
        density_scale, d_density, scf_ws.ortho.d_dwork_nao2_1,
        ri.d_density_factor_accum, ri.d_density_factor_failure);
    // Preserve the exact symmetrized input.  The eigensolver overwrites its
    // matrix with eigenvectors; this copy is used below to certify the actual
    // float factor consumed by RI-K, not merely the double eigenspectrum.
    deviceMemcpy(scf_ws.ortho.d_dwork_nao2_4,
                 scf_ws.ortho.d_dwork_nao2_1, sizeof(double) * nao2,
                 deviceMemcpyDeviceToDevice);

    int nonfinite_index = -1;
    double antisymmetric_squared_norm = 0.0;
    deviceMemcpy(&nonfinite_index, ri.d_density_factor_failure, sizeof(int),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&antisymmetric_squared_norm, ri.d_density_factor_accum,
                 sizeof(double), deviceMemcpyDeviceToHost);
    if (nonfinite_index >= 0)
    {
        float invalid_value = 0.0f;
        deviceMemcpy(&invalid_value, d_density + nonfinite_index,
                     sizeof(float), deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
            "Reason:\n    non-finite %s spin density in RI-K "
            "factorization at matrix index %d: %.9g\n",
            channel_name, nonfinite_index, (double)invalid_value);
        return false;
    }
    if (!Double_Memory_Is_Finite(&antisymmetric_squared_norm) ||
        antisymmetric_squared_norm < 0.0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
            "Reason:\n    non-finite %s spin-density symmetry residual "
            "during RI-K factorization: %.17g\n",
            channel_name, antisymmetric_squared_norm);
        return false;
    }

    const double antisymmetric_rms =
        sqrt(antisymmetric_squared_norm / (double)nao2);
    // Float AO densities necessarily carry rounding noise.  This bound is a
    // small fraction of the requested fixed-point tolerance and is used only
    // to decide whether the represented matrix is a valid PSD density; it is
    // not an SCF convergence relaxation.
    const double representation_rms_tolerance =
        std::max(1.0e-8, 0.25 * scf_ws.runtime.density_tol);
    if (!Double_Memory_Is_Finite(&antisymmetric_rms) ||
        antisymmetric_rms > representation_rms_tolerance)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
            "Reason:\n    %s spin density is not symmetric enough for "
            "RI-K: antisymmetric RMS %.9g exceeds %.9g\n",
            channel_name, antisymmetric_rms, representation_rms_tolerance);
        return false;
    }

    int info = 0;
    const int api_status = QC_Diagonalize_Double(
        solver_handle, nao, scf_ws.ortho.d_dwork_nao2_1,
        scf_ws.ortho.d_dW_double, scf_ws.ortho.d_solver_work_double,
        scf_ws.ortho.lwork_double, &info);
    const bool solver_ok = QC_SCF_Require_Eigensolver_Success(
        QC_SCF_EIGENSOLVER_DENSITY_FACTOR, spin_channel, nao, api_status, info,
        [&](const QC_SCF_Eigensolver_Failure& failure)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
                "Reason:\n    eigensolver failed during %s for channel %s: "
                "dimension=%d, api_status=%d, info=%d\n",
                failure.stage_name, failure.channel_name, failure.dimension,
                failure.api_status, failure.info);
        });
    if (!solver_ok) return false;

    std::vector<double> eigenvalues((size_t)nao);
    deviceMemcpy(eigenvalues.data(), scf_ws.ortho.d_dW_double,
                 sizeof(double) * (size_t)nao, deviceMemcpyDeviceToHost);
    double negative_squared_norm = 0.0;
    int first_positive = nao;
    for (int i = 0; i < nao; ++i)
    {
        const double value = eigenvalues[i];
        if (!Double_Memory_Is_Finite(&value) ||
            (i > 0 && value < eigenvalues[i - 1]))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
                "Reason:\n    invalid ordered eigenvalue during RI-K %s "
                "spin-density factorization: eigenvalue=%d, value=%.17g\n",
                channel_name, i, value);
            return false;
        }
        if (value < 0.0)
            negative_squared_norm += value * value;
        else if (value > 0.0 && first_positive == nao)
            first_positive = i;
    }
    const double negative_projection_rms =
        sqrt(negative_squared_norm / (double)nao2);
    const double total_projection_rms =
        sqrt(0.25 * antisymmetric_rms * antisymmetric_rms +
             negative_projection_rms * negative_projection_rms);
    if (!Double_Memory_Is_Finite(&total_projection_rms) ||
        total_projection_rms > representation_rms_tolerance)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
            "Reason:\n    %s spin density is not positive semidefinite "
            "within its numerical representation: projection RMS %.9g "
            "exceeds %.9g\n",
            channel_name, total_projection_rms,
            representation_rms_tolerance);
        return false;
    }

    *factor_rank = nao - first_positive;
    const int occupied_orbitals =
        spin_channel == QC_SCF_EIGENSOLVER_CHANNEL_BETA
            ? scf_ws.runtime.n_beta
            : scf_ws.runtime.n_alpha;
    if (occupied_orbitals > 0 && *factor_rank <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
            "Reason:\n    %s spin density has zero numerical rank for %d "
            "occupied orbital(s)\n",
            channel_name, occupied_orbitals);
        return false;
    }
    Launch_Device_Kernel(
        QC_Build_RI_Density_Factor_Kernel,
        Positive_Int_Ceil_Div(nao2, threads), threads, 0, 0, nao,
        first_positive, *factor_rank, scf_ws.ortho.d_dwork_nao2_1,
        scf_ws.ortho.d_dW_double, d_factor);

    QC_Float_To_Double(nao2, d_factor, scf_ws.ortho.d_dwork_nao2_2);
    QC_Dgemm_NT(blas_handle, nao, nao, *factor_rank,
                scf_ws.ortho.d_dwork_nao2_2, nao,
                scf_ws.ortho.d_dwork_nao2_2, nao,
                scf_ws.ortho.d_dwork_nao2_3, nao);
    QC_Double_Sub(nao2, scf_ws.ortho.d_dwork_nao2_3,
                  scf_ws.ortho.d_dwork_nao2_4,
                  scf_ws.ortho.d_dwork_nao2_2);
    deviceMemset(ri.d_density_factor_accum, 0, sizeof(double));
    QC_Double_Dot(nao2, scf_ws.ortho.d_dwork_nao2_2,
                  scf_ws.ortho.d_dwork_nao2_2,
                  ri.d_density_factor_accum);
    double reconstruction_squared_norm = 0.0;
    deviceMemcpy(&reconstruction_squared_norm, ri.d_density_factor_accum,
                 sizeof(double), deviceMemcpyDeviceToHost);
    const double reconstructed_density_rms =
        sqrt(reconstruction_squared_norm / (double)nao2 +
             0.25 * antisymmetric_rms * antisymmetric_rms);
    if (!Double_Memory_Is_Finite(&reconstructed_density_rms) ||
        reconstructed_density_rms > representation_rms_tolerance)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Factor_RI_Spin_Density",
            "Reason:\n    reconstructed RI-K factor for channel %s differs "
            "from the current spin density by RMS %.9g, exceeding %.9g\n",
            channel_name, reconstructed_density_rms,
            representation_rms_tolerance);
        return false;
    }
    return true;
}

// B_occ[P*nao + mu + i*M] → B_flat[mu + (P*nocc + i)*nao]
// 将 (naux*nao, nocc) 重排为 (nao, naux*nocc)，使 K = B_flat * B_flat^T
static __global__ void QC_RI_Permute_B_Occ_Kernel(const int total,
                                                  const int nao, const int naux,
                                                  const int nocc,
                                                  const float* __restrict__ src,
                                                  float* __restrict__ dst)
{
    const int M = naux * nao;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        const int mu = idx % nao;
        const int rem = idx / nao;
        const int P = rem / nocc;
        const int i = rem % nocc;
        dst[idx] = src[P * nao + mu + (long long)i * M];
    }
}

// Stored-mode Build_Fock_RI: 使用预存的 d_eri3c 和 d_B
static void Build_Fock_RI_Stored(QUANTUM_CHEMISTRY* qc);
// Direct-mode Build_Fock_RI: 每轮在线计算 3c 积分
static void Build_Fock_RI_Direct(QUANTUM_CHEMISTRY* qc, CONTROLLER* controller);

// RI-JK Fock 构建（入口）
void QUANTUM_CHEMISTRY::Build_Fock_RI()
{
    auto& ri = scf_ws.ri;
    const int nao2 = mol.nao2;
    const int threads = 256;

    // 0. F = H_core + Vxc
    if (dft.enable_dft) Build_DFT_VXC();

    Launch_Device_Kernel(QC_Init_Fock_Kernel, (nao2 + threads - 1) / threads,
                         threads, 0, 0, nao2, scf_ws.core.d_H_core, dft.d_Vxc,
                         dft.enable_dft, scf_ws.alpha.d_F);
    if (scf_ws.runtime.unrestricted)
    {
        Launch_Device_Kernel(QC_Init_Fock_Kernel,
                             (nao2 + threads - 1) / threads, threads, 0, 0,
                             nao2, scf_ws.core.d_H_core, dft.d_Vxc_beta,
                             dft.enable_dft, scf_ws.beta.d_F);
    }

    if (dft.exx_fraction != 0.0f)
    {
        if (scf_ws.runtime.n_alpha > 0 &&
            !Factor_RI_Spin_Density(
                scf_ws.alpha.d_P,
                scf_ws.runtime.unrestricted ? 1.0 : 0.5,
                ri.d_density_factor_alpha,
                &ri.density_factor_rank_alpha,
                QC_SCF_EIGENSOLVER_CHANNEL_ALPHA))
            return;
        if (scf_ws.runtime.unrestricted && scf_ws.runtime.n_beta > 0 &&
            !Factor_RI_Spin_Density(
                scf_ws.beta.d_P, 1.0, ri.d_density_factor_beta,
                &ri.density_factor_rank_beta,
                QC_SCF_EIGENSOLVER_CHANNEL_BETA))
            return;
    }

    if (ri.direct)
        Build_Fock_RI_Direct(this, controller);
    else
        Build_Fock_RI_Stored(this);

    // H + Vxc (+ RI-K) is accumulated through the existing float path.  Add
    // the double-precision RI-J contribution only after promoting that base,
    // so the Fock operator retains the same precision as the RI-J functional.
    auto assemble_fock_double = [&](QC_SCF_Spin_Channel& channel)
    {
        QC_Float_To_Double_Copy(nao2, channel.d_F, channel.d_F_double);
        QC_Double_Axpy(nao2, 1.0, ri.d_J_double, channel.d_F_double);
    };
    assemble_fock_double(scf_ws.alpha);
    if (scf_ws.runtime.unrestricted) assemble_fock_double(scf_ws.beta);
}

// Stored mode
// RI-K helper: B_occ → permute → single GEMM → K
// 将 naux 次小 SGEMM 替换为 1 次 permute kernel + 1 次大 SGEMM
static void RI_K_Build_Single_GEMM(BLAS_HANDLE blas_handle, QC_RI_WORKSPACE& ri,
                                   const int naux, const int nao,
                                   const int factor_rank,
                                   const float* d_density_factor, float* d_F,
                                   const float neg_exx, const int threads)
{
    const int nao2 = nao * nao;
    const int M = naux * nao;
    const float one_f = 1.0f, zero_f = 0.0f;

    // B_occ[M,rank] = B^T[M,nao] * L^T[nao,rank], P_sigma=L L^T.
    deviceBlasSgemm(blas_handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_T, M,
                    factor_rank, nao, &one_f, ri.d_B, nao, d_density_factor,
                    nao, &zero_f, ri.d_B_occ, M);

    // Permute: B_occ[(naux*nao), nocc] → B_flat[nao, (naux*nocc)]
    // B_occ[P*nao+mu, i] → B_flat[mu, P*nocc+i]
    const int total = nao * naux * factor_rank;
    Launch_Device_Kernel(QC_RI_Permute_B_Occ_Kernel,
                         (total + threads - 1) / threads, threads, 0, 0, total,
                         nao, naux, factor_rank, ri.d_B_occ, ri.d_B_flat);

    // K = B_flat × B_flat^T  (single large GEMM)
    // K[nao, nao] = B_flat[nao, naux*nocc] × B_flat^T[naux*nocc, nao]
    const int K_dim = naux * factor_rank;
    deviceBlasSgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_T, nao, nao,
                    K_dim, &one_f, ri.d_B_flat, nao, ri.d_B_flat, nao, &zero_f,
                    ri.d_K_scratch, nao);

    // F += neg_exx * K
    Launch_Device_Kernel(QC_Scaled_Add_Kernel, (nao2 + threads - 1) / threads,
                         threads, 0, 0, nao2, neg_exx, ri.d_K_scratch, d_F);
}

static void Build_Fock_RI_Stored(QUANTUM_CHEMISTRY* qc)
{
    auto& ri = qc->scf_ws.ri;
    auto& scf_ws = qc->scf_ws;
    auto& mol = qc->mol;
    auto& dft = qc->dft;
    auto blas_handle = qc->blas_handle;
    const int naux = ri.naux;
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    const int threads = 256;

    // ---- RI-J（使用持久化缓冲，无 malloc/free）----
    {
        const float* d_P_coul = scf_ws.runtime.unrestricted
                                    ? scf_ws.direct.d_Ptot
                                    : scf_ws.alpha.d_P;

        QC_Float_To_Double(nao2, d_P_coul, ri.d_P_double);

        const double one = 1.0, zero = 0.0;
        // d_vec[P] = eri3c^T * P
        deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_N, naux,
                        1, nao2, &one, ri.d_eri3c, nao2, ri.d_P_double, nao2,
                        &zero, ri.d_d_vec, naux);

        // g = (P|Q)^{-1} * d_vec
        deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, naux,
                        1, naux, &one, ri.d_metric_inv, naux, ri.d_d_vec, naux,
                        &zero, ri.d_g_vec, naux);

        // J = eri3c * g
        deviceBlasDgemm(blas_handle, DEVICE_BLAS_OP_N, DEVICE_BLAS_OP_N, nao2,
                        1, naux, &one, ri.d_eri3c, nao2, ri.d_g_vec, naux,
                        &zero, ri.d_J_double, nao2);

    }

    // ---- RI-K: permute + single GEMM（替代 naux 次小 SGEMM）----
    if (dft.exx_fraction != 0.0f)
    {
        const float neg_exx = -dft.exx_fraction;

        const int rank_a = ri.density_factor_rank_alpha;
        if (rank_a > 0)
            RI_K_Build_Single_GEMM(
                blas_handle, ri, naux, nao, rank_a,
                ri.d_density_factor_alpha, scf_ws.alpha.d_F, neg_exx, threads);

        if (scf_ws.runtime.unrestricted)
        {
            const int rank_b = ri.density_factor_rank_beta;
            if (rank_b > 0)
                RI_K_Build_Single_GEMM(
                    blas_handle, ri, naux, nao, rank_b,
                    ri.d_density_factor_beta, scf_ws.beta.d_F, neg_exx,
                    threads);
        }
    }
}

// Direct mode
// 在线计算 3c 积分，不存储 eri3c/B 张量。
// 策略：逐轨道 shell pair 启动 GPU 3c kernel，结果到临时缓冲，
//       host 端做 cart2sph + 收缩到 d_vec / B_occ。
// RI-J: pass1 累加 d_P, solve g, pass2 累加 J（需两趟 3c）
// RI-K: 一趟累加 B_occ，最后 K = B_occ^T B_occ
//
// 当前实现：先在 host 上完成全部 shell-pair 循环 + 收缩。
// 后续可优化为 GPU kernel batch + device-side 收缩。
static void Build_Fock_RI_Direct(QUANTUM_CHEMISTRY* qc, CONTROLLER* controller)
{
    auto& ri = qc->scf_ws.ri;
    auto& scf_ws = qc->scf_ws;
    auto& mol = qc->mol;
    auto& dft = qc->dft;
    auto blas_handle = qc->blas_handle;
    const int naux = ri.naux;
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    const int threads = 256;

    // 下载到 host
    std::vector<float> h_D_f(nao2);
    const float* d_P_coul =
        scf_ws.runtime.unrestricted ? scf_ws.direct.d_Ptot : scf_ws.alpha.d_P;
    deviceMemcpy(h_D_f.data(), d_P_coul, sizeof(float) * nao2,
                 deviceMemcpyDeviceToHost);
    std::vector<double> h_D(nao2);
    for (int i = 0; i < nao2; i++) h_D[i] = (double)h_D_f[i];

    // 复用 RI_Precompute 阶段缓存的 host metric (不随 SCF 迭代改变)
    const std::vector<double>& h_inv = ri.h_metric_inv;
    const std::vector<double>& h_inv_sqrt = ri.h_metric_inv_sqrt;
    std::vector<float> h_orb_norms(nao);
    deviceMemcpy(h_orb_norms.data(), scf_ws.ortho.d_norms, sizeof(float) * nao,
                 deviceMemcpyDeviceToHost);

    // 临时 3c 缓冲（direct 模式按 shell-pair 紧凑输出）
    int max_orb_l = 0;
    for (int sh = 0; sh < mol.nbas; sh++)
        max_orb_l = std::max(max_orb_l, mol.h_l_list[sh]);
    int max_aux_l = 0;
    for (int l : ri.h_aux_l_list) max_aux_l = std::max(max_aux_l, l);
    size_t max_cart_size = 0;
    size_t buf_3c_size = 0;
    size_t buf_3c_bytes = 0;
    if (!QC_RI_Checked_Cartesian_Count(max_orb_l, &max_cart_size) ||
        !QC_RI_Checked_Mul_Size((size_t)ri.naux_cart, max_cart_size,
                                &buf_3c_size) ||
        !QC_RI_Checked_Mul_Size(buf_3c_size, max_cart_size, &buf_3c_size) ||
        !QC_RI_Checked_Bytes(buf_3c_size, sizeof(double), &buf_3c_bytes))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_Fock_RI_Direct",
            "Reason:\n    RI direct three-center output buffer size "
            "overflows size_t\n");
        return;
    }
    const int max_cart = (int)max_cart_size;

    QC_RI_INTEGRAL_WORKSPACE eri3c_workspace;
    if (!QC_RI_Build_3Center_Workspace_Layout(ri.naux_bas, max_aux_l, max_orb_l,
                                              &eri3c_workspace))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_Fock_RI_Direct",
            "Reason:\n    RI direct three-center workspace dimensions "
            "overflow for auxiliary/orbital angular momenta %d/%d and %d "
            "tasks\n",
            max_aux_l, max_orb_l, ri.naux_bas);
        return;
    }
    if (!QC_RI_Allocate_Integral_Workspace(&eri3c_workspace)) return;

    double* d_3c_buf = NULL;
    if (!Device_Malloc_Safely((void**)&d_3c_buf, buf_3c_bytes))
    {
        QC_RI_Free_Integral_Workspace(&eri3c_workspace);
        return;
    }

    // 准备 RI-J: d_vec, J (host)
    std::vector<double> h_d_vec(naux, 0.0);
    std::vector<double> h_J(nao2, 0.0);

    // Prepare RI-K factors from the actual spin densities.  These factors
    // were validated and built by Build_Fock_RI before entering direct mode.
    const bool need_exx = (dft.exx_fraction != 0.0f);
    const int rank_a = need_exx ? ri.density_factor_rank_alpha : 0;
    const int rank_b =
        (need_exx && scf_ws.runtime.unrestricted)
            ? ri.density_factor_rank_beta
            : 0;

    std::vector<float> h_factor_a, h_factor_b;
    std::vector<double> h_B_occ_a, h_B_occ_b;
    if (rank_a > 0)
    {
        h_factor_a.resize((size_t)nao * nao);
        deviceMemcpy(h_factor_a.data(), ri.d_density_factor_alpha,
                     sizeof(float) * nao * nao, deviceMemcpyDeviceToHost);
        h_B_occ_a.assign((size_t)naux * nao * rank_a, 0.0);
    }
    if (rank_b > 0)
    {
        h_factor_b.resize((size_t)nao * nao);
        deviceMemcpy(h_factor_b.data(), ri.d_density_factor_beta,
                     sizeof(float) * nao * nao, deviceMemcpyDeviceToHost);
        h_B_occ_b.assign((size_t)naux * nao * rank_b, 0.0);
    }

    // pass 1: 逐 shell pair 计算 3c 积分
    // 用 GPU kernel 计算笛卡尔 3c block，下载，cart2sph，收缩
    // 对每个 shell pair (mu_sh, nu_sh):
    //   tasks = 所有辅助 shell P_sh
    //   launch kernel → d_3c_buf [naux_cart × dim_mu_c × dim_nu_c]
    //   download, cart2sph → block_sph [naux × dim_mu_s × dim_nu_s]
    //   d_vec[P] += Σ_{μν} block_sph[P,μ,ν] * D[off_mu+μ, off_nu+ν]
    //   B_occ[P,off_mu+μ,i] += Σ_ν block_B[P,μ,ν] * C[off_nu+ν,i]

    // 分配 GPU 3c 任务缓冲（一次性，最多所有 P_sh × 1 pair）
    QC_RI_3C_TASK* d_tasks = NULL;
    if (!Device_Malloc_Safely((void**)&d_tasks,
                              sizeof(QC_RI_3C_TASK) *
                                  (size_t)ri.naux_bas))
    {
        deviceFree(d_3c_buf);
        QC_RI_Free_Integral_Workspace(&eri3c_workspace);
        return;
    }

    for (int mu_sh = 0; mu_sh < mol.nbas; mu_sh++)
    {
        const int l_mu = mol.h_l_list[mu_sh];
        const int dmc = (int)QC_RI_Cartesian_Count(l_mu);
        const int dms = mol.is_spherical ? (2 * l_mu + 1) : dmc;
        const int off_mu_s = mol.is_spherical ? mol.h_ao_offsets_sph[mu_sh]
                                              : mol.h_ao_offsets[mu_sh];

        for (int nu_sh = 0; nu_sh <= mu_sh; nu_sh++)
        {
            const int l_nu = mol.h_l_list[nu_sh];
            const int dnc = (int)QC_RI_Cartesian_Count(l_nu);
            const int dns = mol.is_spherical ? (2 * l_nu + 1) : dnc;
            const int off_nu_s = mol.is_spherical ? mol.h_ao_offsets_sph[nu_sh]
                                                  : mol.h_ao_offsets[nu_sh];

            std::vector<QC_RI_3C_TASK> h_tasks;
            h_tasks.reserve(ri.naux_bas);
            for (int P = 0; P < ri.naux_bas; P++)
                h_tasks.push_back({P, mu_sh, nu_sh});
            const int task_count = (int)h_tasks.size();
            if (task_count > 0)
                deviceMemcpy(d_tasks, h_tasks.data(),
                             sizeof(QC_RI_3C_TASK) * (size_t)task_count,
                             deviceMemcpyHostToDevice);

            // 清零缓冲，启动 kernel
            const size_t buf_n =
                (size_t)ri.naux_cart * (size_t)dmc * (size_t)dnc;
            deviceMemset(d_3c_buf, 0, sizeof(double) * buf_n);
            QC_Launch_RI_3Center_Kernel(
                threads, task_count, d_tasks, ri.d_aux_centers,
                ri.d_aux_l_list, ri.d_aux_exps, ri.d_aux_coeffs,
                ri.d_aux_shell_offsets, ri.d_aux_shell_sizes,
                ri.d_aux_ao_offsets, mol.d_centers, mol.d_l_list, mol.d_exps,
                mol.d_coeffs, mol.d_shell_offsets, mol.d_shell_sizes,
                mol.d_ao_offsets, ri.naux_cart, dmc, dnc,
                mol.h_ao_offsets[mu_sh], mol.h_ao_offsets[nu_sh], false,
                eri3c_workspace, d_3c_buf);

            // 下载笛卡尔 block
            std::vector<double> h_block_cart(buf_n);
            deviceMemcpy(h_block_cart.data(), d_3c_buf, sizeof(double) * buf_n,
                         deviceMemcpyDeviceToHost);

            // Cart2sph: block_cart[Pc, dmc, dnc] → block_sph[Ps, dms, dns]
            std::vector<double> block_sph((size_t)naux * dms * dns, 0.0);
            if (!mol.is_spherical)
            {
                const int Pc = ri.naux_cart;
                for (int ps = 0; ps < naux; ps++)
                    for (int pc = 0; pc < Pc; pc++)
                    {
                        double u = (double)ri.h_U_aux[pc * naux + ps];
                        if (u == 0.0) continue;
                        for (int i = 0; i < dms; i++)
                            for (int j = 0; j < dns; j++)
                                block_sph[ps * dms * dns + i * dns + j] +=
                                    u * h_block_cart[(long long)pc * dmc * dnc +
                                                     (long long)i * dnc + j];
                    }
            }
            else
            {
                // 先提取笛卡尔 block: [Pc, dmc, dnc]
                const int Pc = ri.naux_cart;
                std::vector<double> bc((size_t)Pc * dmc * dnc, 0.0);
                for (int P = 0; P < Pc; P++)
                    for (int i = 0; i < dmc; i++)
                        for (int j = 0; j < dnc; j++)
                            bc[P * dmc * dnc + i * dnc + j] =
                                h_block_cart[(long long)P * dmc * dnc +
                                             (long long)i * dnc + j];

                // ν 变换: T1[Pc, dmc, dns]
                std::vector<double> t1((size_t)Pc * dmc * dns, 0.0);
                for (int P = 0; P < Pc; P++)
                    for (int i = 0; i < dmc; i++)
                        for (int js = 0; js < dns; js++)
                            for (int jc = 0; jc < dnc; jc++)
                                t1[P * dmc * dns + i * dns + js] +=
                                    bc[P * dmc * dnc + i * dnc + jc] *
                                    (double)ri.h_U_orb
                                        [(mol.h_ao_offsets[nu_sh] + jc) * nao +
                                         off_nu_s + js];

                // μ 变换: T2[Pc, dms, dns]
                std::vector<double> t2((size_t)Pc * dms * dns, 0.0);
                for (int P = 0; P < Pc; P++)
                    for (int is_ = 0; is_ < dms; is_++)
                        for (int js = 0; js < dns; js++)
                            for (int ic = 0; ic < dmc; ic++)
                                t2[P * dms * dns + is_ * dns + js] +=
                                    (double)ri.h_U_orb
                                        [(mol.h_ao_offsets[mu_sh] + ic) * nao +
                                         off_mu_s + is_] *
                                    t1[P * dmc * dns + ic * dns + js];

                // P 变换: block_sph[Ps, dms, dns]
                for (int ps = 0; ps < naux; ps++)
                    for (int pc = 0; pc < Pc; pc++)
                    {
                        double u = (double)ri.h_U_aux[pc * naux + ps];
                        if (u == 0.0) continue;
                        for (int mn = 0; mn < dms * dns; mn++)
                            block_sph[ps * dms * dns + mn] +=
                                u * t2[pc * dms * dns + mn];
                    }
            }

            for (int P = 0; P < naux; P++)
            {
                const double p_scale = (double)ri.h_aux_norms[P];
                for (int i = 0; i < dms; i++)
                {
                    const double mu_scale = (double)h_orb_norms[off_mu_s + i];
                    for (int j = 0; j < dns; j++)
                    {
                        block_sph[P * dms * dns + i * dns + j] *=
                            p_scale * mu_scale *
                            (double)h_orb_norms[off_nu_s + j];
                    }
                }
            }

            // inv_sqrt @ block → B_block [naux, dms, dns]
            std::vector<double> B_block((size_t)naux * dms * dns, 0.0);
            for (int ps = 0; ps < naux; ps++)
                for (int qs = 0; qs < naux; qs++)
                {
                    double w = h_inv_sqrt[ps * naux + qs];
                    if (w == 0.0) continue;
                    for (int mn = 0; mn < dms * dns; mn++)
                        B_block[ps * dms * dns + mn] +=
                            w * block_sph[qs * dms * dns + mn];
                }

            // RI-J: d_vec[P] += Σ_{μν} block_sph[P,μ,ν] * D[μ,ν]
            for (int P = 0; P < naux; P++)
                for (int i = 0; i < dms; i++)
                    for (int j = 0; j < dns; j++)
                    {
                        double val = block_sph[P * dms * dns + i * dns + j];
                        h_d_vec[P] +=
                            val * h_D[(off_mu_s + i) * nao + (off_nu_s + j)];
                        if (mu_sh != nu_sh)
                            h_d_vec[P] +=
                                val *
                                h_D[(off_nu_s + j) * nao + (off_mu_s + i)];
                    }

            // RI-K: B_occ[P, off_mu+i, occ] += Σ_j B_block[P,i,j] * C[off_nu+j,
            // occ]
            if (rank_a > 0)
            {
                for (int P = 0; P < naux; P++)
                    for (int i = 0; i < dms; i++)
                        for (int j = 0; j < dns; j++)
                        {
                            double b = B_block[P * dms * dns + i * dns + j];
                            if (b == 0.0) continue;
                            int mu_idx = off_mu_s + i;
                            int nu_idx = off_nu_s + j;
                            for (int oc = 0; oc < rank_a; oc++)
                            {
                                h_B_occ_a[(long long)P * nao * rank_a +
                                          mu_idx * rank_a + oc] +=
                                    b * (double)h_factor_a[nu_idx * nao + oc];
                                if (mu_sh != nu_sh)
                                    h_B_occ_a[(long long)P * nao * rank_a +
                                              nu_idx * rank_a + oc] +=
                                        b * (double)
                                                h_factor_a[mu_idx * nao + oc];
                            }
                        }
            }
            if (rank_b > 0)
            {
                for (int P = 0; P < naux; P++)
                    for (int i = 0; i < dms; i++)
                        for (int j = 0; j < dns; j++)
                        {
                            double b = B_block[P * dms * dns + i * dns + j];
                            if (b == 0.0) continue;
                            int mu_idx = off_mu_s + i;
                            int nu_idx = off_nu_s + j;
                            for (int oc = 0; oc < rank_b; oc++)
                            {
                                h_B_occ_b[(long long)P * nao * rank_b +
                                          mu_idx * rank_b + oc] +=
                                    b * (double)h_factor_b[nu_idx * nao + oc];
                                if (mu_sh != nu_sh)
                                    h_B_occ_b[(long long)P * nao * rank_b +
                                              nu_idx * rank_b + oc] +=
                                        b * (double)
                                                h_factor_b[mu_idx * nao + oc];
                            }
                        }
            }
        }
    }
    deviceFree(d_tasks);
    d_tasks = NULL;

    // RI-J: g = inv @ d, J[μν] = Σ_P g_P block_sph[P,μ,ν]
    // 但 pass2 需要再算一遍 3c 积分... 用 metric_inv 直接算:
    // J[μν] = Σ_P (Σ_Q inv[P,Q] d_Q) block_sph[P,μ,ν]
    // 我们已经在 pass1 中积累了 d_vec。现在需要 pass2。
    // 但 pass2 和 pass1 一样需要遍历所有 shell pair。
    // 优化：可以在 pass1 同时积累 J（如果 g 已知）。
    // 但 g 依赖 d_vec 的完整值，所以必须两趟。
    //
    // 方案 B：用 metric_inv_sqrt，一趟就够。
    // d' = inv_sqrt @ d → g' = d'（因为 inv = inv_sqrt @ inv_sqrt）
    // J = Σ_P g'_P * B_block[P] = 两趟都需要 B_block...
    //
    // 最终方案：pass1 积累 d_vec，算 g_vec，pass2 重新计算 3c 积累 J。
    // 为避免重复计算 3c，我们在 pass1 中同时积累 J 的贡献。
    //
    // 实际上更聪明的做法：
    // J[μν] = Σ_PQ (P|μν) (P|Q)^{-1} (Q|ρσ) D[ρσ]
    //       = Σ_P (P|μν) g_P  where g = inv @ d
    // 需要 g 才能算 J，但 g 需要完整的 d。
    // 所以确实需要两趟。
    //
    // 但注意 RI-K 的 B_occ 在 pass1 中已经一趟算完。
    // 对 RI-J，我们可以重新计算 3c 做 pass2。
    // 或者：在 pass1 中存储每个 block_sph[P, μ_sh, ν_sh] 到一个
    // 压缩格式中。但这等于变相存储了 3c 张量。
    //
    // 最简洁的做法：pass1 计算 d_vec。然后:
    //   g = inv @ d
    //   pass2: 重新计算每个 shell pair 的 3c block，
    //          J[off_mu+i, off_nu+j] += Σ_P g[P] * block[P,i,j]

    // 先求 g
    std::vector<double> h_g(naux, 0.0);
    for (int P = 0; P < naux; P++)
        for (int Q = 0; Q < naux; Q++)
            h_g[P] += h_inv[P * naux + Q] * h_d_vec[Q];

    // pass 2: RI-J 第二趟，重新计算 3c 积分以构建 J
    if (!Device_Malloc_Safely((void**)&d_tasks,
                              sizeof(QC_RI_3C_TASK) *
                                  (size_t)ri.naux_bas))
    {
        deviceFree(d_3c_buf);
        QC_RI_Free_Integral_Workspace(&eri3c_workspace);
        return;
    }
    for (int mu_sh = 0; mu_sh < mol.nbas; mu_sh++)
    {
        const int l_mu = mol.h_l_list[mu_sh];
        const int dmc = (int)QC_RI_Cartesian_Count(l_mu);
        const int dms = mol.is_spherical ? (2 * l_mu + 1) : dmc;
        const int off_mu_s = mol.is_spherical ? mol.h_ao_offsets_sph[mu_sh]
                                              : mol.h_ao_offsets[mu_sh];

        for (int nu_sh = 0; nu_sh <= mu_sh; nu_sh++)
        {
            const int l_nu = mol.h_l_list[nu_sh];
            const int dnc = (int)QC_RI_Cartesian_Count(l_nu);
            const int dns = mol.is_spherical ? (2 * l_nu + 1) : dnc;
            const int off_nu_s = mol.is_spherical ? mol.h_ao_offsets_sph[nu_sh]
                                                  : mol.h_ao_offsets[nu_sh];

            std::vector<QC_RI_3C_TASK> h_tasks;
            h_tasks.reserve(ri.naux_bas);
            for (int P = 0; P < ri.naux_bas; P++)
                h_tasks.push_back({P, mu_sh, nu_sh});
            const int task_count = (int)h_tasks.size();
            if (task_count > 0)
                deviceMemcpy(d_tasks, h_tasks.data(),
                             sizeof(QC_RI_3C_TASK) * (size_t)task_count,
                             deviceMemcpyHostToDevice);

            const size_t buf_n =
                (size_t)ri.naux_cart * (size_t)dmc * (size_t)dnc;
            deviceMemset(d_3c_buf, 0, sizeof(double) * buf_n);
            QC_Launch_RI_3Center_Kernel(
                threads, task_count, d_tasks, ri.d_aux_centers,
                ri.d_aux_l_list, ri.d_aux_exps, ri.d_aux_coeffs,
                ri.d_aux_shell_offsets, ri.d_aux_shell_sizes,
                ri.d_aux_ao_offsets, mol.d_centers, mol.d_l_list, mol.d_exps,
                mol.d_coeffs, mol.d_shell_offsets, mol.d_shell_sizes,
                mol.d_ao_offsets, ri.naux_cart, dmc, dnc,
                mol.h_ao_offsets[mu_sh], mol.h_ao_offsets[nu_sh], false,
                eri3c_workspace, d_3c_buf);

            std::vector<double> h_block_cart(buf_n);
            deviceMemcpy(h_block_cart.data(), d_3c_buf, sizeof(double) * buf_n,
                         deviceMemcpyDeviceToHost);

            // cart2sph (同 pass1 逻辑)
            std::vector<double> block_sph((size_t)naux * dms * dns, 0.0);
            if (!mol.is_spherical)
            {
                const int Pc = ri.naux_cart;
                for (int ps = 0; ps < naux; ps++)
                    for (int pc = 0; pc < Pc; pc++)
                    {
                        double u = (double)ri.h_U_aux[pc * naux + ps];
                        if (u == 0.0) continue;
                        for (int i = 0; i < dms; i++)
                            for (int j = 0; j < dns; j++)
                                block_sph[ps * dms * dns + i * dns + j] +=
                                    u * h_block_cart[(long long)pc * dmc * dnc +
                                                     (long long)i * dnc + j];
                    }
            }
            else
            {
                const int Pc = ri.naux_cart;
                std::vector<double> bc((size_t)Pc * dmc * dnc, 0.0);
                for (int P = 0; P < Pc; P++)
                    for (int i = 0; i < dmc; i++)
                        for (int j = 0; j < dnc; j++)
                            bc[P * dmc * dnc + i * dnc + j] =
                                h_block_cart[(long long)P * dmc * dnc +
                                             (long long)i * dnc + j];
                std::vector<double> t1((size_t)Pc * dmc * dns, 0.0);
                for (int P = 0; P < Pc; P++)
                    for (int i = 0; i < dmc; i++)
                        for (int js = 0; js < dns; js++)
                            for (int jc = 0; jc < dnc; jc++)
                                t1[P * dmc * dns + i * dns + js] +=
                                    bc[P * dmc * dnc + i * dnc + jc] *
                                    (double)ri.h_U_orb
                                        [(mol.h_ao_offsets[nu_sh] + jc) * nao +
                                         off_nu_s + js];
                std::vector<double> t2((size_t)Pc * dms * dns, 0.0);
                for (int P = 0; P < Pc; P++)
                    for (int is_ = 0; is_ < dms; is_++)
                        for (int js = 0; js < dns; js++)
                            for (int ic = 0; ic < dmc; ic++)
                                t2[P * dms * dns + is_ * dns + js] +=
                                    (double)ri.h_U_orb
                                        [(mol.h_ao_offsets[mu_sh] + ic) * nao +
                                         off_mu_s + is_] *
                                    t1[P * dmc * dns + ic * dns + js];
                for (int ps = 0; ps < naux; ps++)
                    for (int pc = 0; pc < Pc; pc++)
                    {
                        double u = (double)ri.h_U_aux[pc * naux + ps];
                        if (u == 0.0) continue;
                        for (int mn = 0; mn < dms * dns; mn++)
                            block_sph[ps * dms * dns + mn] +=
                                u * t2[pc * dms * dns + mn];
                    }
            }

            for (int P = 0; P < naux; P++)
            {
                const double p_scale = (double)ri.h_aux_norms[P];
                for (int i = 0; i < dms; i++)
                {
                    const double mu_scale = (double)h_orb_norms[off_mu_s + i];
                    for (int j = 0; j < dns; j++)
                    {
                        block_sph[P * dms * dns + i * dns + j] *=
                            p_scale * mu_scale *
                            (double)h_orb_norms[off_nu_s + j];
                    }
                }
            }

            // J[off_mu+i, off_nu+j] += Σ_P g[P] * block_sph[P,i,j]
            for (int i = 0; i < dms; i++)
                for (int j = 0; j < dns; j++)
                {
                    double val = 0.0;
                    for (int P = 0; P < naux; P++)
                        val += h_g[P] * block_sph[P * dms * dns + i * dns + j];
                    h_J[(off_mu_s + i) * nao + (off_nu_s + j)] += val;
                    if (mu_sh != nu_sh)
                        h_J[(off_nu_s + j) * nao + (off_mu_s + i)] += val;
                }
        }
    }
    deviceFree(d_tasks);
    deviceFree(d_3c_buf);
    QC_RI_Free_Integral_Workspace(&eri3c_workspace);

    // Keep direct-mode RI-J in the same double buffer used by stored mode.
    // Build_Fock_RI combines it with each spin channel after RI-K is complete.
    deviceMemcpy(ri.d_J_double, h_J.data(), sizeof(double) * nao2,
                 deviceMemcpyHostToDevice);

    // RI-K
    if (need_exx)
    {
        const float neg_exx = -dft.exx_fraction;

        auto build_K_from_B_occ =
            [&](const std::vector<double>& h_B_occ, int nocc, float* d_F)
        {
            if (nocc <= 0) return;
            // K[μ,ν] = Σ_{P,i} B_occ[P,μ,i] * B_occ[P,ν,i]
            std::vector<float> h_K(nao2, 0.0f);
            for (int P = 0; P < naux; P++)
                for (int mu = 0; mu < nao; mu++)
                    for (int nu = 0; nu <= mu; nu++)
                    {
                        double sum = 0.0;
                        for (int oc = 0; oc < nocc; oc++)
                            sum += h_B_occ[(long long)P * nao * nocc +
                                           mu * nocc + oc] *
                                   h_B_occ[(long long)P * nao * nocc +
                                           nu * nocc + oc];
                        h_K[mu * nao + nu] += (float)sum;
                        if (mu != nu) h_K[nu * nao + mu] += (float)sum;
                    }
            // 上传并加到 F
            float* d_K = NULL;
            Device_Malloc_Safely((void**)&d_K, sizeof(float) * nao2);
            deviceMemcpy(d_K, h_K.data(), sizeof(float) * nao2,
                         deviceMemcpyHostToDevice);
            Launch_Device_Kernel(QC_Scaled_Add_Kernel,
                                 (nao2 + threads - 1) / threads, threads, 0, 0,
                                 nao2, neg_exx, d_K, d_F);
            deviceFree(d_K);
        };

        build_K_from_B_occ(h_B_occ_a, rank_a, scf_ws.alpha.d_F);
        if (scf_ws.runtime.unrestricted)
            build_K_from_B_occ(h_B_occ_b, rank_b, scf_ws.beta.d_F);
    }
}
