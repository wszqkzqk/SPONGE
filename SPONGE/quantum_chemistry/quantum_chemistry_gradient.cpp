// clang-format off
#include "quantum_chemistry.h"
#include "gradient/grad_one_e.hpp"
#include "gradient/grad_workspace.h"
#include "gradient/gradient.hpp"
#include "integrals/eri/common/direct_fock_kernels.hpp"
#include "integrals/eri/common/eri_rys.hpp"
#include "integrals/eri/eri_backend.hpp"
#include "integrals/ri/ri_3center.hpp"
#include "gradient/grad_eri.hpp"
#include "gradient/grad_ri.hpp"
#include "ecp/ecp_integrals.h"
// clang-format on

std::vector<float> QC_Build_Cart2Sph_Mat_Host(const std::vector<int>& l_list,
                                              int nao_cart, int nao_sph);

static __global__ void QC_Weight_By_Norms_Kernel(int nao, int matrix_size,
                                                 const float* P,
                                                 const float* norms, float* out)
{
    SIMPLE_DEVICE_FOR(idx, matrix_size)
    {
        int i = idx / nao, j = idx % nao;
        out[idx] = P[idx] * norms[i] * norms[j];
    }
}

// Prepare Cartesian contractions explicitly instead of relying on a side
// effect of one particular one-electron-gradient backend.  U is shared by P
// and W, so copy it and the AO norms only once per gradient evaluation.
static void QC_Prepare_Spherical_Gradient_Densities(
    int nao_sph, int nao_cart, const float* d_norms,
    const float* d_transform_cart_sph, const float* d_P_sph, float* d_P_cart,
    const float* d_W_sph, float* d_W_cart)
{
    std::vector<float> h_norms((size_t)nao_sph);
    std::vector<float> h_transform((size_t)nao_cart * nao_sph);
    std::vector<float> h_P((size_t)nao_sph * nao_sph);
    deviceMemcpy(h_norms.data(), d_norms, sizeof(float) * nao_sph,
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(h_transform.data(), d_transform_cart_sph,
                 sizeof(float) * (size_t)nao_cart * nao_sph,
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(h_P.data(), d_P_sph,
                 sizeof(float) * (size_t)nao_sph * nao_sph,
                 deviceMemcpyDeviceToHost);

    std::vector<float> h_P_cart;
    QC_Sph2Cart_Density_Host(nao_sph, nao_cart, h_norms, h_transform, h_P,
                             h_P_cart);
    deviceMemcpy(d_P_cart, h_P_cart.data(),
                 sizeof(float) * (size_t)nao_cart * nao_cart,
                 deviceMemcpyHostToDevice);

    if (d_W_sph != nullptr)
    {
        std::vector<float> h_W((size_t)nao_sph * nao_sph);
        deviceMemcpy(h_W.data(), d_W_sph,
                     sizeof(float) * (size_t)nao_sph * nao_sph,
                     deviceMemcpyDeviceToHost);
        std::vector<float> h_W_cart;
        QC_Sph2Cart_Density_Host(nao_sph, nao_cart, h_norms, h_transform, h_W,
                                 h_W_cart);
        deviceMemcpy(d_W_cart, h_W_cart.data(),
                     sizeof(float) * (size_t)nao_cart * nao_cart,
                     deviceMemcpyHostToDevice);
    }
}

void QUANTUM_CHEMISTRY::Compute_Gradient(
    VECTOR* local_frc, const VECTOR* global_crd, const int* global_to_local,
    int owned_atom_numbers, const VECTOR box_length, int need_virial,
    LTMatrix3* local_atom_virial)
{
    if (!is_initialized) return;
    (void)global_crd;
    const int natm = mol.natm;
    const int nao = mol.nao;
    const int nao2 = mol.nao2;

    deviceMemset(grad_ws.d_grad, 0, sizeof(double) * natm * 3);

    // 1. Build W from the accepted P and the physical F[P] that produced the
    // accepted energy.  This remains valid for fractional ensemble densities
    // and avoids replacing them with an unrelated integer Aufbau projector.
    {
        const double* alpha_fock = scf_ws.alpha.d_F_for_grad;
        const double* beta_fock = scf_ws.beta.d_F_for_grad;
        if (alpha_fock == nullptr ||
            (scf_ws.runtime.unrestricted && beta_fock == nullptr))
        {
            Build_Fock(1, true);
            alpha_fock = scf_ws.alpha.d_F_double;
            beta_fock = scf_ws.beta.d_F_double;
        }

        const int nao_effective =
            scf_ws.ortho.nao_eff > 0 ? scf_ws.ortho.nao_eff : nao;
        double* d_overlap_pseudoinverse =
            scf_ws.ortho.d_dwork_nao2_1;
        double* d_density_double = scf_ws.ortho.d_dwork_nao2_2;
        double* d_fock_density = scf_ws.ortho.d_dwork_nao2_3;
        double* d_weighted_density_double =
            scf_ws.ortho.d_dwork_nao2_4;
        QC_Build_Overlap_Pseudoinverse(
            blas_handle, nao, nao_effective, scf_ws.ortho.d_X,
            d_overlap_pseudoinverse);
        QC_Build_Energy_Weighted_Density_From_PF(
            blas_handle, nao, d_overlap_pseudoinverse, alpha_fock,
            scf_ws.alpha.d_P, grad_ws.d_W_density, d_density_double,
            d_fock_density, d_weighted_density_double, false);

        if (scf_ws.runtime.unrestricted && grad_ws.d_W_density_beta)
        {
            QC_Build_Energy_Weighted_Density_From_PF(
                blas_handle, nao, d_overlap_pseudoinverse, beta_fock,
                scf_ws.beta.d_P, grad_ws.d_W_density, d_density_double,
                d_fock_density, d_weighted_density_double, true);
        }
    }

    // 2. 核排斥梯度
    {
        const int threads = 256;
        const VECTOR box_bohr(box_length.x * CONSTANT_ANGSTROM_TO_BOHR,
                              box_length.y * CONSTANT_ANGSTROM_TO_BOHR,
                              box_length.z * CONSTANT_ANGSTROM_TO_BOHR);
        Launch_Device_Kernel(QC_Nuclear_Gradient_Kernel,
                             Positive_Int_Ceil_Div(natm, threads), threads, 0,
                             0,
                             natm, mol.d_Z, mol.d_atm, mol.d_env, box_bohr,
                             periodic_boundary, grad_ws.d_grad);
    }

    // 3. 单电子积分导数: Tr[P·dH/dR] - Tr[W·dS/dR]
    {
#ifdef USE_GPU
        constexpr bool cartesian_one_e_backend = true;
#else
        constexpr bool cartesian_one_e_backend = false;
#endif
        // The GPU one-electron backend consumes Cartesian P/W.  ECP
        // derivatives consume Cartesian P on both CPU and GPU.  Prepare those
        // dependencies before selecting the one-electron backend so the CPU
        // spherical path cannot leave ECP with stale workspace contents.
        if (mol.is_spherical && (cartesian_one_e_backend || mol.has_ecp))
        {
            QC_Prepare_Spherical_Gradient_Densities(
                mol.nao, mol.nao_cart, scf_ws.ortho.d_norms,
                cart2sph.d_cart2sph_mat, scf_ws.direct.d_P_coul,
                grad_ws.d_P_cart,
                cartesian_one_e_backend ? grad_ws.d_W_density : nullptr,
                cartesian_one_e_backend ? grad_ws.d_W_cart : nullptr);
        }

#ifndef USE_GPU
        if (mol.is_spherical)
        {
            std::vector<int> h_shell_atom(mol.nbas);
            for (int ish = 0; ish < mol.nbas; ish++)
                h_shell_atom[ish] = mol.h_bas[ish * 8 + 0];
            std::vector<float> h_cart2sph;
            try
            {
                h_cart2sph = QC_Build_Cart2Sph_Mat_Host(
                    mol.h_l_list, mol.nao_cart, mol.nao_sph);
            }
            catch (const std::exception& error)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow,
                    "QUANTUM_CHEMISTRY::Build_Full_Gradient",
                    "Reason:\n    Failed to construct one-electron gradient "
                    "Cartesian-to-spherical matrix: %s\n",
                    error.what());
                return;
            }

            QC_Build_OneE_Gradient_Spherical_CPU(
                task_ctx.topo.h_1e_tasks, mol.h_centers, mol.h_l_list,
                mol.h_exps, mol.h_coeffs, mol.h_shell_offsets,
                mol.h_shell_sizes, mol.h_ao_offsets, mol.h_ao_offsets_sph,
                mol.h_atm, mol.h_env, h_shell_atom, scf_ws.direct.d_P_coul,
                grad_ws.d_W_density, scf_ws.ortho.d_norms, h_cart2sph.data(),
                mol.natm, mol.nao_sph, grad_ws.d_grad);
        }
        else
#endif
        {
            const int n_tasks = task_ctx.topo.n_1e_tasks;
            const float* d_P_use;
            const float* d_W_use;
            const float* d_norms_use;
            int nao_1e;

            if (mol.is_spherical)
            {
                d_P_use = grad_ws.d_P_cart;
                d_W_use = grad_ws.d_W_cart;
                d_norms_use = grad_ws.d_norms_ones;
                nao_1e = mol.nao_cart;
            }
            else
            {
                d_P_use = scf_ws.direct.d_P_coul;
                d_W_use = grad_ws.d_W_density;
                d_norms_use = scf_ws.ortho.d_norms;
                nao_1e = mol.nao;
            }

            Launch_Device_Kernel(OneE_ST_Grad_Kernel, (n_tasks + 63) / 64, 64,
                                 0, 0, n_tasks, task_ctx.buffers.d_1e_tasks,
                                 mol.d_centers, mol.d_l_list, mol.d_exps,
                                 mol.d_coeffs, mol.d_shell_offsets,
                                 mol.d_shell_sizes, mol.d_ao_offsets, nao_1e,
                                 grad_ws.d_shell_atom, d_P_use, d_W_use,
                                 d_norms_use, grad_ws.d_grad);

            const int v_total = n_tasks * mol.natm;
            Launch_Device_Kernel(
                OneE_V_Grad_Kernel, (v_total + 63) / 64, 64, 0, 0, n_tasks,
                task_ctx.buffers.d_1e_tasks, mol.d_centers, mol.d_l_list,
                mol.d_exps, mol.d_coeffs, mol.d_shell_offsets,
                mol.d_shell_sizes, mol.d_ao_offsets, mol.d_atm, mol.d_env,
                mol.natm, nao_1e, grad_ws.d_shell_atom, d_P_use, d_norms_use,
                grad_ws.d_grad);
        }
    }

    // 3b. ECP 梯度
    if (mol.has_ecp)
    {
        if (!mol.is_spherical)
        {
            // 非球谐基: 需要计算 P_cart = P .* (norms * norms')
            Launch_Device_Kernel(QC_Weight_By_Norms_Kernel,
                                 Positive_Int_Ceil_Div(mol.nao_cart2, 256), 256,
                                 0, 0, mol.nao_cart, mol.nao_cart2,
                                 scf_ws.direct.d_P_coul, scf_ws.ortho.d_norms,
                                 grad_ws.d_P_cart);
        }
        // 球谐与笛卡尔两条路径都在上面显式准备 d_P_cart。
        QC_ECP_EVALUATION_FAILURE failure;
        if (!QC_Compute_ECP_Gradient(mol, task_ctx, grad_ws.d_shell_atom,
                                     grad_ws.d_P_cart, grad_ws.d_grad,
                                     &failure))
        {
            if (failure.kind == QC_ECP_RESOURCE_ALLOCATION_FAILED)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::Compute_Gradient",
                    "Reason:\n    %s\n",
                    QC_ECP_Evaluation_Failure_Kind_Name(failure.kind));
                return;
            }
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Compute_Gradient",
                "Reason:\n    ECP gradient evaluation failed: %s. Maximum "
                "angular series order q=%d; context: ECP atom %d, channel %d "
                "(l=%d), term %d (n_k=%d), one-electron task %d, primitives "
                "(%d,%d), shells (%d,%d), Cartesian functions (%d,%d), "
                "observable atom/direction (%d,%d); accumulated signed ECP "
                "observable %.17g, contracted absolute error estimate %.9g\n",
                QC_ECP_Evaluation_Failure_Kind_Name(failure.kind),
                QC_ECP_MAX_SERIES_ORDER, failure.atom, failure.channel,
                failure.channel_l, failure.term, failure.n_k,
                failure.task_id, failure.primitive_i, failure.primitive_j,
                failure.shell_i, failure.shell_j, failure.cartesian_i,
                failure.cartesian_j, failure.observable_atom,
                failure.direction, failure.value, failure.estimated_error);
            return;
        }
    }

    // 4. 双电子积分导数: Tr[Γ·dERI/dR]
    // grad_eri 内部在 Cartesian shell buffer 上计算导数积分，
    // is_spherical 时内部做 cart2sph，因此始终传入 SCF AO 基的密度和 norms。
    const float grad_shell_screen_tol =
        fmaxf(task_ctx.params.eri_shell_screen_tol, 1.0e-7f);
    const float grad_prim_screen_tol =
        fmaxf(task_ctx.params.direct_eri_prim_screen_tol, 1.0e-7f);
    if (scf_ws.ri.enabled)
    {
        Build_RI_Gradient();
    }
    else
    {
#ifndef USE_GPU
        QC_Build_ERI_Gradient_CPU(task_ctx, mol, cart2sph, scf_ws,
                                  dft.exx_fraction, grad_shell_screen_tol,
                                  grad_prim_screen_tol, grad_ws.d_shell_atom,
                                  grad_ws.d_grad);
#else
        // GPU ERI gradient: reuse screening infrastructure, launch gradient
        // kernel
        {
            // 0. Refresh pair density for screening (Fock build may have
            //    left stale incremental values)
            {
                const int threads_pd = 256;
                const bool need_exx_pd = (dft.exx_fraction != 0.0f);
                Launch_Device_Kernel(
                    QC_Build_Shell_Pair_Density_Kernel,
                    (task_ctx.topo.n_shell_pairs + threads_pd - 1) / threads_pd,
                    threads_pd, 0, 0, task_ctx.topo.n_shell_pairs,
                    task_ctx.buffers.d_shell_pairs, mol.d_ao_offsets,
                    mol.d_ao_offsets_sph, mol.d_l_list, mol.is_spherical, nao,
                    scf_ws.direct.d_P_coul, scf_ws.direct.d_pair_density_coul,
                    need_exx_pd ? scf_ws.alpha.d_P : (const float*)nullptr,
                    scf_ws.direct.d_pair_density_exx,
                    (need_exx_pd && scf_ws.runtime.unrestricted)
                        ? scf_ws.beta.d_P
                        : (const float*)nullptr,
                    scf_ws.direct.d_pair_density_exx_b);
            }

            // 1. Run screening (same as Fock build)
            deviceMemset(task_ctx.buffers.d_screen_counts, 0,
                         sizeof(int) * task_ctx.topo.n_combos);

            const int cp_needed = task_ctx.topo.n_combos + 1;
            deviceMemcpy(grad_ws.d_combo_prefix_grad,
                         (void*)task_ctx.topo.combo_prefix,
                         sizeof(int) * cp_needed, deviceMemcpyHostToDevice);

            const float exx_a = scf_ws.runtime.unrestricted
                                    ? dft.exx_fraction
                                    : (0.5f * dft.exx_fraction);
            const float exx_b =
                scf_ws.runtime.unrestricted ? dft.exx_fraction : 0.0f;

            QC_Launch_Screen(
                task_ctx.topo.total_quartets, task_ctx.buffers.d_combos,
                grad_ws.d_combo_prefix_grad, task_ctx.topo.n_combos,
                task_ctx.buffers.d_sorted_pair_ids,
                task_ctx.buffers.d_shell_pairs,
                task_ctx.buffers.d_shell_pair_bounds,
                scf_ws.direct.d_pair_density_coul,
                scf_ws.direct.d_pair_density_exx,
                scf_ws.runtime.unrestricted ? scf_ws.direct.d_pair_density_exx_b
                                            : (const float*)nullptr,
                grad_shell_screen_tol, exx_a, exx_b,
                task_ctx.buffers.d_screened_tasks,
                task_ctx.buffers.d_screen_counts);

            int h_counts[QC_INTEGRAL_TASKS::MAX_COMBOS] = {};
            deviceMemcpy(h_counts, task_ctx.buffers.d_screen_counts,
                         sizeof(int) * task_ctx.topo.n_combos,
                         deviceMemcpyDeviceToHost);

            // Count total screened tasks
            int total_screened = 0;
            for (int ci = 0; ci < task_ctx.topo.n_combos; ci++)
            {
                total_screened += h_counts[ci];
            }

            if (total_screened > 0)
            {
                // 2. 清零多副本梯度缓冲 (在 Memory_Allocate 中已预分配)
                const int grad_size = natm * 3;
                const size_t copies_needed =
                    (size_t)QC_GRAD_N_COPIES * (size_t)grad_size;
                deviceMemset(grad_ws.d_grad_copies, 0,
                             sizeof(double) * copies_needed);

                const int gamma_buf_size = grad_ws.grad_gamma_buf_size;

                // 3. Launch gradient kernel per combo
                const int threads = QC_GRAD_ERI_THREADS;
                for (int ci = 0; ci < task_ctx.topo.n_combos; ci++)
                {
                    const int n = h_counts[ci];
                    if (n == 0) continue;
                    const QC_ERI_TASK* d_tasks =
                        task_ctx.buffers.d_screened_tasks +
                        task_ctx.topo.h_combos[ci].output_offset;
                    const int max_blocks_pool =
                        grad_ws.grad_gamma_pool_slots / threads;
                    const int blocks = std::min({QC_GRAD_GAMMA_POOL_BLOCKS,
                                                 Positive_Int_Ceil_Div(n,
                                                                       threads),
                                                 max_blocks_pool});

                    Launch_Device_Kernel(
                        QC_ERI_Grad_Kernel, blocks, threads, 0, 0, n, d_tasks,
                        mol.d_atm, mol.d_bas, mol.d_env, mol.d_ao_offsets,
                        mol.d_ao_offsets_sph, scf_ws.ortho.d_norms,
                        task_ctx.buffers.d_shell_pair_bounds,
                        scf_ws.direct.d_pair_density_coul,
                        scf_ws.direct.d_pair_density_exx,
                        scf_ws.runtime.unrestricted
                            ? scf_ws.direct.d_pair_density_exx_b
                            : (const float*)nullptr,
                        grad_shell_screen_tol, scf_ws.direct.d_P_coul,
                        scf_ws.alpha.d_P,
                        scf_ws.runtime.unrestricted ? scf_ws.beta.d_P
                                                    : (const float*)nullptr,
                        exx_a, exx_b, nao, mol.nao_sph, mol.is_spherical,
                        cart2sph.d_cart2sph_mat, grad_ws.d_grad_gamma_pool,
                        gamma_buf_size, grad_ws.d_shell_atom,
                        grad_ws.d_grad_copies, QC_GRAD_N_COPIES, natm,
                        grad_prim_screen_tol);
                }

                // 4. Reduce gradient copies
                Launch_Device_Kernel(QC_Reduce_Grad_Copies_Kernel,
                                     (grad_size + 255) / 256, 256, 0, 0,
                                     grad_size, QC_GRAD_N_COPIES,
                                     grad_ws.d_grad_copies, grad_ws.d_grad);
            }
        }
#endif
    }

    // 5. DFT XC 网格梯度
    if (dft.enable_dft) Build_DFT_XC_Gradient();

    std::vector<double> h_grad((size_t)natm * 3);
    deviceMemcpy(h_grad.data(), grad_ws.d_grad, sizeof(double) * h_grad.size(),
                 deviceMemcpyDeviceToHost);
    for (int i = 0; i < natm; i++)
    {
        for (int component = 0; component < 3; component++)
        {
            const double value = h_grad[(size_t)i * 3 + component];
            const double force_kcal_mol_angstrom =
                -value * (double)CONSTANT_HARTREE_TO_KCAL_MOL *
                (double)CONSTANT_ANGSTROM_TO_BOHR;
            const float force_to_accumulate = (float)force_kcal_mol_angstrom;
            if (!Double_Memory_Is_Finite(&value) ||
                !Double_Memory_Is_Finite(&force_kcal_mol_angstrom) ||
                !Float_Memory_Is_Finite(&force_to_accumulate))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::Compute_Gradient",
                    "Reason:\n    non-finite QC gradient for global atom %d "
                    "component %d or unrepresentable MD force: %.17g "
                    "Hartree/Bohr, %.17g kcal/mol/Angstrom\n",
                    atom_local[i], component, value, force_kcal_mol_angstrom);
                return;
            }
        }
    }

    // Every QC atom must have exactly one owning PP rank.  Non-owning ranks
    // still skip the local write in the kernel, but the global invariant is
    // validated here so no force can disappear silently.
    if (global_to_local == nullptr || owned_atom_numbers < 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Compute_Gradient",
            "Reason:\n    invalid domain-decomposition mapping for QC force "
            "writeback\n");
        return;
    }
    std::vector<int> h_global_to_local((size_t)atom_numbers);
    deviceMemcpy(h_global_to_local.data(), global_to_local,
                 sizeof(int) * h_global_to_local.size(),
                 deviceMemcpyDeviceToHost);
    std::vector<int> local_owner(natm, 0);
    std::vector<int> owner_count(natm, 0);
    for (int i = 0; i < natm; i++)
    {
        const int global_atom = atom_local[i];
        if (global_atom < 0 || global_atom >= atom_numbers)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Compute_Gradient",
                "Reason:\n    QC atom %d has invalid global atom index %d "
                "for a %d-atom MD system\n",
                i, global_atom, atom_numbers);
            return;
        }
        const int local_atom = h_global_to_local[global_atom];
        local_owner[i] =
            local_atom >= 0 && local_atom < owned_atom_numbers ? 1 : 0;
    }
#ifdef USE_MPI
    if (CONTROLLER::PP_MPI_size > 1)
    {
        MPI_Allreduce(local_owner.data(), owner_count.data(), natm, MPI_INT,
                      MPI_SUM, CONTROLLER::pp_comm);
    }
    else
#endif
    {
        owner_count = local_owner;
    }
    for (int i = 0; i < natm; i++)
    {
        if (owner_count[i] != 1)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Compute_Gradient",
                "Reason:\n    QC global atom %d has %d owning PP ranks; "
                "exactly one is required for force writeback\n",
                atom_local[i], owner_count[i]);
            return;
        }
    }

    // 6. 将梯度写入 MD 力数组
    {
        const int threads = 256;
        Launch_Device_Kernel(
            QC_Writeback_Gradient_Kernel,
            Positive_Int_Ceil_Div(natm, threads),
            threads, 0, 0, natm, d_atom_local, global_to_local,
            owned_atom_numbers, grad_ws.d_grad, mol.d_atom_coords, local_frc,
            need_virial, local_atom_virial);
    }
}

// RI (Density Fitting) 解析梯度
// 支持 stored 模式（直接下载预存的 eri3c）和 direct 模式（逐 shell pair
// 重新计算 3c 积分）。 两种模式最终都构建相同的 D3_eff / D2_eff
// 有效密度，调用相同的梯度内核。
void QUANTUM_CHEMISTRY::Build_RI_Gradient()
{
    auto& ri = scf_ws.ri;
    const int natm = mol.natm;
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    const int naux = ri.naux;
    const bool need_exx = (dft.exx_fraction != 0.0f);

    if (naux <= 0 || ri.naux_eff <= 0 || ri.naux_eff > naux ||
        ri.h_metric_inv_sqrt.size() != (size_t)naux * naux ||
        ri.h_eigval.size() != (size_t)naux ||
        ri.h_eigvec.size() != (size_t)naux * naux ||
        (!ri.direct && need_exx &&
         (ri.d_B == nullptr || ri.d_B_occ == nullptr)))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Build_RI_Gradient",
            "Reason:\n    incomplete RI metric/exchange workspace: "
            "naux=%d, naux_eff=%d, direct=%d, exact_exchange=%d\n",
            naux, ri.naux_eff, ri.direct, need_exx);
        return;
    }

    if (need_exx)
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

    // 复用 RI_Precompute 阶段缓存的 host metric
    const std::vector<double>& h_metric_inv_sqrt = ri.h_metric_inv_sqrt;
    std::vector<float> h_P(nao2);
    std::vector<float> h_orb_norms(nao);

    const float* d_P_coul =
        scf_ws.runtime.unrestricted ? scf_ws.direct.d_Ptot : scf_ws.alpha.d_P;
    deviceMemcpy(h_P.data(), d_P_coul, sizeof(float) * nao2,
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(h_orb_norms.data(), scf_ws.ortho.d_norms, sizeof(float) * nao,
                 deviceMemcpyDeviceToHost);

    size_t M_size = 0;
    if (!QC_RI_Checked_Mul_Size((size_t)naux, (size_t)nao, &M_size) ||
        M_size > (size_t)std::numeric_limits<int>::max())
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_RI_Gradient",
            "Reason:\n    RI flattened auxiliary/orbital dimension "
            "overflows int: naux=%d, nao=%d\n",
            naux, nao);
        return;
    }
    const int M = (int)M_size;

    struct RI_GRAD_SPIN_CHANNEL
    {
        int factor_rank = 0;
        const float* d_factor = nullptr;
        std::vector<float> density;
        std::vector<float> factor;
        std::vector<double> B_occ_double;
        std::vector<float> B_occ;
    };
    std::vector<RI_GRAD_SPIN_CHANNEL> spin_channels;

    auto add_spin_channel =
        [&](int factor_rank, const float* d_density, const float* d_factor)
    {
        if (!need_exx || factor_rank <= 0) return;
        RI_GRAD_SPIN_CHANNEL channel;
        channel.factor_rank = factor_rank;
        channel.d_factor = d_factor;
        channel.density.resize(nao2);
        channel.factor.resize((size_t)nao * nao);
        deviceMemcpy(channel.density.data(), d_density, sizeof(float) * nao2,
                     deviceMemcpyDeviceToHost);
        deviceMemcpy(channel.factor.data(), d_factor,
                     sizeof(float) * (size_t)nao * nao,
                     deviceMemcpyDeviceToHost);
        spin_channels.push_back(std::move(channel));
    };

    add_spin_channel(ri.density_factor_rank_alpha, scf_ws.alpha.d_P,
                     ri.d_density_factor_alpha);
    if (scf_ws.runtime.unrestricted)
        add_spin_channel(ri.density_factor_rank_beta, scf_ws.beta.d_P,
                         ri.d_density_factor_beta);

    // 计算 max shell cart sizes（workspace 分配用，从 host 数据）
    int max_aux_cart = 0;
    int max_aux_l = 0;
    for (int i = 0; i < ri.naux_bas; i++)
    {
        max_aux_l = std::max(max_aux_l, ri.h_aux_l_list[i]);
        int nc = (int)QC_RI_Cartesian_Count(ri.h_aux_l_list[i]);
        if (nc > max_aux_cart) max_aux_cart = nc;
    }
    int max_orb_cart = 0;
    int max_orb_l = 0;
    for (int i = 0; i < mol.nbas; i++)
    {
        max_orb_l = std::max(max_orb_l, mol.h_l_list[i]);
        int nc = (int)QC_RI_Cartesian_Count(mol.h_l_list[i]);
        if (nc > max_orb_cart) max_orb_cart = nc;
    }

    // 启动梯度内核的公共 lambda（两种模式共用）
    auto launch_grad_kernels = [&](const std::vector<double>& D2_eff,
                                   const std::vector<double>& D3_eff) -> bool
    {
        return QC_Launch_RI_Grad_Kernels(mol, ri, scf_ws.ortho.d_norms, grad_ws,
                                         controller, D2_eff, D3_eff);
    };

    if (ri.direct)
    {
        // Direct 模式：增量累积，不存储完整 eri3c 张量
        //   Pass 1: 逐 shell pair 累积 d_vec + B_occ
        //   Pass 2 (仅 EXX): 逐 shell pair 累积 Z_K
        //   内存: O(naux·nao·nocc) + O(naux²)，而非 O(naux·nao²)

        size_t n_shell_pairs_twice = 0;
        if (!QC_RI_Checked_Mul_Size((size_t)mol.nbas, (size_t)mol.nbas + 1,
                                    &n_shell_pairs_twice))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_RI_Gradient",
                "Reason:\n    RI orbital shell-pair count overflows size_t\n");
            return;
        }
        const size_t n_shell_pairs = n_shell_pairs_twice / 2;

        // 复用 RI_Precompute 阶段缓存
        const std::vector<double>& h_metric_inv = ri.h_metric_inv;

        // 转换密度为 double
        std::vector<double> h_D(nao2);
        for (int i = 0; i < nao2; i++) h_D[i] = (double)h_P[i];

        // GPU 3c 缓冲
        const size_t max_cart = (size_t)max_orb_cart;
        size_t buf_3c_size = 0;
        size_t buf_3c_bytes = 0;
        if (!QC_RI_Checked_Mul_Size((size_t)ri.naux_cart, max_cart,
                                    &buf_3c_size) ||
            !QC_RI_Checked_Mul_Size(buf_3c_size, max_cart, &buf_3c_size) ||
            !QC_RI_Checked_Bytes(buf_3c_size, sizeof(double), &buf_3c_bytes))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_RI_Gradient",
                "Reason:\n    RI direct-gradient three-center output buffer "
                "size overflows size_t\n");
            return;
        }

        QC_RI_INTEGRAL_WORKSPACE eri3c_workspace;
        if (!QC_RI_Build_3Center_Workspace_Layout(ri.naux_bas, max_aux_l,
                                                  max_orb_l, &eri3c_workspace))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "QUANTUM_CHEMISTRY::Build_RI_Gradient",
                "Reason:\n    RI direct-gradient three-center workspace "
                "dimensions overflow for auxiliary/orbital angular momenta "
                "%d/%d and %d tasks\n",
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
        std::vector<QC_RI_3C_TASK> h_tasks;
        h_tasks.reserve(ri.naux_bas);
        QC_RI_3C_TASK* d_tasks = NULL;
        if (!Device_Malloc_Safely((void**)&d_tasks,
                                  sizeof(QC_RI_3C_TASK) * (size_t)ri.naux_bas))
        {
            deviceFree(d_3c_buf);
            QC_RI_Free_Integral_Workspace(&eri3c_workspace);
            return;
        }
        const int threads = 256;

        // 辅助 lambda: 计算一个 shell pair 的 block_sph
        // 复用 cart2sph + 归一化逻辑
        auto compute_block_sph = [&](int mu_sh, int nu_sh, int dmc, int dnc,
                                     int dms, int dns, int off_mu_s,
                                     int off_nu_s, std::vector<double>& out)
        {
            h_tasks.clear();
            for (int auxiliary_shell = 0; auxiliary_shell < ri.naux_bas;
                 ++auxiliary_shell)
                h_tasks.push_back({auxiliary_shell, mu_sh, nu_sh});
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
            std::vector<double> h_bc(buf_n);
            deviceMemcpy(h_bc.data(), d_3c_buf, sizeof(double) * buf_n,
                         deviceMemcpyDeviceToHost);

            out.assign((size_t)naux * dms * dns, 0.0);
            const int Pc = ri.naux_cart;
            if (!mol.is_spherical)
            {
                for (int ps = 0; ps < naux; ps++)
                    for (int pc = 0; pc < Pc; pc++)
                    {
                        double u = (double)ri.h_U_aux[pc * naux + ps];
                        if (u == 0.0) continue;
                        for (int ij = 0; ij < dms * dns; ij++)
                            out[ps * dms * dns + ij] +=
                                u * h_bc[(long long)pc * dmc * dnc + ij];
                    }
            }
            else
            {
                std::vector<double> t1((size_t)Pc * dmc * dns, 0.0);
                for (int P = 0; P < Pc; P++)
                    for (int i = 0; i < dmc; i++)
                        for (int js = 0; js < dns; js++)
                            for (int jc = 0; jc < dnc; jc++)
                                t1[P * dmc * dns + i * dns + js] +=
                                    h_bc[(long long)P * dmc * dnc + i * dnc +
                                         jc] *
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
                            out[ps * dms * dns + mn] +=
                                u * t2[pc * dms * dns + mn];
                    }
            }
            for (int P = 0; P < naux; P++)
            {
                double ps = (double)ri.h_aux_norms[P];
                for (int i = 0; i < dms; i++)
                {
                    double ms = ps * (double)h_orb_norms[off_mu_s + i];
                    for (int j = 0; j < dns; j++)
                        out[P * dms * dns + i * dns + j] *=
                            ms * (double)h_orb_norms[off_nu_s + j];
                }
            }
        };

        // Pass 1: 累积 d_vec 和 B_occ, 同时缓存 3c 积分块
        std::vector<double> h_d_vec(naux, 0.0);
        for (auto& channel : spin_channels)
            channel.B_occ_double.assign((size_t)M * channel.factor_rank, 0.0);

        // 缓存 3c 积分块以复用于 Pass 2 (避免重复 GPU kernel 调用)
        std::vector<std::vector<double>> blk_cache(need_exx ? n_shell_pairs
                                                            : 0);

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
                const int off_nu_s = mol.is_spherical
                                         ? mol.h_ao_offsets_sph[nu_sh]
                                         : mol.h_ao_offsets[nu_sh];

                std::vector<double> blk;
                compute_block_sph(mu_sh, nu_sh, dmc, dnc, dms, dns, off_mu_s,
                                  off_nu_s, blk);

                // d_vec[P] += block·D
                for (int P = 0; P < naux; P++)
                    for (int i = 0; i < dms; i++)
                        for (int j = 0; j < dns; j++)
                        {
                            double v = blk[P * dms * dns + i * dns + j];
                            h_d_vec[P] +=
                                v * h_D[(off_mu_s + i) * nao + (off_nu_s + j)];
                            if (mu_sh != nu_sh)
                                h_d_vec[P] +=
                                    v *
                                    h_D[(off_nu_s + j) * nao + (off_mu_s + i)];
                        }

                // Accumulate B times the validated factor of each actual
                // spin density.  For RKS the factor represents P/2; for UKS
                // it represents the channel density itself.
                if (need_exx)
                {
                    // B_block = metric_inv_sqrt @ block
                    std::vector<double> B_blk((size_t)naux * dms * dns, 0.0);
                    for (int P = 0; P < naux; P++)
                        for (int Q = 0; Q < naux; Q++)
                        {
                            double w = h_metric_inv_sqrt[(size_t)P * naux + Q];
                            if (w == 0.0) continue;
                            for (int mn = 0; mn < dms * dns; mn++)
                                B_blk[P * dms * dns + mn] +=
                                    w * blk[Q * dms * dns + mn];
                        }
                    // B_occ[(P*nao+μ)+M*k] += Σ_j B_blk[P,i,j] L[ν,k]
                    for (auto& channel : spin_channels)
                        for (int P = 0; P < naux; P++)
                            for (int i = 0; i < dms; i++)
                                for (int j = 0; j < dns; j++)
                                {
                                    const double b =
                                        B_blk[P * dms * dns + i * dns + j];
                                    if (b == 0.0) continue;
                                    const int mu_idx = off_mu_s + i;
                                    const int nu_idx = off_nu_s + j;
                                    for (int oc = 0;
                                         oc < channel.factor_rank; oc++)
                                    {
                                        channel.B_occ_double[(size_t)(P * nao +
                                                                      mu_idx) +
                                                             (size_t)M * oc] +=
                                            b *
                                            (double)
                                                channel.factor[nu_idx * nao +
                                                               oc];
                                        if (mu_sh != nu_sh)
                                            channel
                                                .B_occ_double[(size_t)(P * nao +
                                                                       nu_idx) +
                                                              (size_t)M * oc] +=
                                                b * (double)channel.factor
                                                        [mu_idx * nao + oc];
                                    }
                                }

                    // 缓存 blk 用于 Pass 2
                    const size_t pair_idx =
                        (size_t)mu_sh * ((size_t)mu_sh + 1) / 2 + (size_t)nu_sh;
                    blk_cache[pair_idx] = std::move(blk);
                }
            }
        }

        // g_vec = metric_inv @ d_vec
        std::vector<double> h_g(naux, 0.0);
        for (int P = 0; P < naux; P++)
            for (int Q = 0; Q < naux; Q++)
                h_g[P] += h_metric_inv[(size_t)P * naux + Q] * h_d_vec[Q];

        // Pass 2 (仅 EXX): 按自旋通道累积 Z_K，复用同一份 3c block。
        std::vector<double> h_Z_K((size_t)naux * naux, 0.0);
        auto accumulate_direct_z_k = [&](const RI_GRAD_SPIN_CHANNEL& channel)
        {
            const int nocc = channel.factor_rank;
            // R[P', oc, m]
            std::vector<double> R((size_t)naux * nocc * nao, 0.0);
            for (int Pp = 0; Pp < naux; Pp++)
                for (int oc = 0; oc < nocc; oc++)
                    for (int m = 0; m < nao; m++)
                    {
                        double sum = 0.0;
                        for (int n = 0; n < nao; n++)
                            sum += channel.B_occ_double[(size_t)(Pp * nao + n) +
                                                        (size_t)M * oc] *
                                   (double)channel.density[m * nao + n];
                        R[(long long)Pp * nocc * nao + oc * nao + m] = sum;
                    }

            const int max_sh_sph =
                mol.is_spherical ? (2 * max_orb_l + 1) : max_orb_cart;
            std::vector<double> T((size_t)naux * max_sh_sph * nocc);
            std::vector<double> Tt((size_t)naux * max_sh_sph * nocc);

            for (int mu_sh = 0; mu_sh < mol.nbas; mu_sh++)
            {
                const int l_mu = mol.h_l_list[mu_sh];
                const int dms = mol.is_spherical
                                    ? (2 * l_mu + 1)
                                    : ((int)QC_RI_Cartesian_Count(l_mu));
                const int off_mu_s = mol.is_spherical
                                         ? mol.h_ao_offsets_sph[mu_sh]
                                         : mol.h_ao_offsets[mu_sh];
                for (int nu_sh = 0; nu_sh <= mu_sh; nu_sh++)
                {
                    const int l_nu = mol.h_l_list[nu_sh];
                    const int dns = mol.is_spherical
                                        ? (2 * l_nu + 1)
                                        : ((int)QC_RI_Cartesian_Count(l_nu));
                    const int off_nu_s = mol.is_spherical
                                             ? mol.h_ao_offsets_sph[nu_sh]
                                             : mol.h_ao_offsets[nu_sh];

                    // 从缓存取出 blk (避免重复计算 3c 积分)
                    const size_t pair_idx =
                        (size_t)mu_sh * ((size_t)mu_sh + 1) / 2 + (size_t)nu_sh;
                    const std::vector<double>& blk = blk_cache[pair_idx];

                    // T[Q,i,k] = Σ_j blk[Q,i,j] * L[off_nu+j,k]
                    const size_t T_size = (size_t)naux * dms * nocc;
                    std::fill_n(T.begin(), T_size, 0.0);
                    for (int Q = 0; Q < naux; Q++)
                        for (int i = 0; i < dms; i++)
                            for (int j = 0; j < dns; j++)
                            {
                                double v = blk[Q * dms * dns + i * dns + j];
                                if (v == 0.0) continue;
                                for (int oc = 0; oc < nocc; oc++)
                                    T[(long long)Q * dms * nocc + i * nocc +
                                      oc] +=
                                        v * (double)channel.factor
                                                [(off_nu_s + j) * nao + oc];
                            }

                    // Z_K[P',Q'] += -Σ_{i,oc} R[P',oc,off_mu+i] * T[Q',i,oc]
                    for (int Pp = 0; Pp < naux; Pp++)
                        for (int Q = 0; Q < naux; Q++)
                        {
                            double z = 0.0;
                            for (int i = 0; i < dms; i++)
                                for (int oc = 0; oc < nocc; oc++)
                                    z += R[(long long)Pp * nocc * nao +
                                           oc * nao + (off_mu_s + i)] *
                                         T[(long long)Q * dms * nocc +
                                           i * nocc + oc];
                            h_Z_K[(size_t)Pp * naux + Q] -= z;
                        }

                    if (mu_sh != nu_sh)
                    {
                        const size_t Tt_size = (size_t)naux * dns * nocc;
                        std::fill_n(Tt.begin(), Tt_size, 0.0);
                        for (int Q = 0; Q < naux; Q++)
                            for (int i = 0; i < dms; i++)
                                for (int j = 0; j < dns; j++)
                                {
                                    double v = blk[Q * dms * dns + i * dns + j];
                                    if (v == 0.0) continue;
                                    for (int oc = 0; oc < nocc; oc++)
                                        Tt[(long long)Q * dns * nocc +
                                           j * nocc + oc] +=
                                            v *
                                            (double)channel.factor
                                                [(off_mu_s + i) * nao + oc];
                                }
                        for (int Pp = 0; Pp < naux; Pp++)
                            for (int Q = 0; Q < naux; Q++)
                            {
                                double z = 0.0;
                                for (int j = 0; j < dns; j++)
                                    for (int oc = 0; oc < nocc; oc++)
                                        z += R[(long long)Pp * nocc * nao +
                                               oc * nao + (off_nu_s + j)] *
                                             Tt[(long long)Q * dns * nocc +
                                                j * nocc + oc];
                                h_Z_K[(size_t)Pp * naux + Q] -= z;
                            }
                    }
                }
            }
        };
        for (const auto& channel : spin_channels)
            accumulate_direct_z_k(channel);
        blk_cache.clear();  // 释放缓存

        // B_occ: double -> float (Pass 2 使用 double 完成后再转换)
        for (auto& channel : spin_channels)
        {
            const size_t count = (size_t)M * channel.factor_rank;
            channel.B_occ.resize(count);
            for (size_t idx = 0; idx < count; idx++)
                channel.B_occ[idx] = (float)channel.B_occ_double[idx];
        }

        deviceFree(d_tasks);
        deviceFree(d_3c_buf);
        QC_RI_Free_Integral_Workspace(&eri3c_workspace);

        // 构建 D3_eff 和 D2_eff
        std::vector<double> D3_eff;
        QC_Init_D3_J(nao, naux, h_g.data(), h_P.data(), D3_eff);
        for (const auto& channel : spin_channels)
            QC_Accumulate_D3_K_Channel(nao, naux, channel.density.data(),
                                       h_metric_inv_sqrt.data(),
                                       channel.B_occ.data(),
                                       channel.factor.data(),
                                       channel.factor_rank, dft.exx_fraction,
                                       D3_eff);

        std::vector<double> D2_eff;
        QC_Init_D2_J(naux, h_g.data(), D2_eff);
        if (!spin_channels.empty())
            QC_Accumulate_D2_K_DaleckiiKrein(
                naux, ri.naux_eff, h_Z_K.data(), ri.h_eigval.data(),
                ri.h_eigvec.data(), dft.exx_fraction, D2_eff);

        if (!launch_grad_kernels(D2_eff, D3_eff)) return;
    }
    else
    {
        // Stored 模式：下载预存的 eri3c，调用原始函数
        std::vector<double> h_eri3c((size_t)naux * nao2);
        deviceMemcpy(h_eri3c.data(), ri.d_eri3c,
                     sizeof(double) * (size_t)naux * nao2,
                     deviceMemcpyDeviceToHost);

        // h_g_vec 仅 stored 模式需要 (direct 模式自行计算 h_g)
        std::vector<double> h_g_vec(naux);
        deviceMemcpy(h_g_vec.data(), ri.d_g_vec, sizeof(double) * naux,
                     deviceMemcpyDeviceToHost);

        // 使用共享 device scratch 依次构建并下载各自旋通道 B_occ。
        if (need_exx && ri.d_B != nullptr && ri.d_B_occ != nullptr)
        {
            const float one_f = 1.0f, zero_f = 0.0f;
            for (auto& channel : spin_channels)
            {
                deviceBlasSgemm(blas_handle, DEVICE_BLAS_OP_T, DEVICE_BLAS_OP_T,
                                M, channel.factor_rank, nao, &one_f, ri.d_B,
                                nao, channel.d_factor, nao, &zero_f,
                                ri.d_B_occ, M);
                const size_t count = (size_t)M * channel.factor_rank;
                channel.B_occ.resize(count);
                deviceMemcpy(channel.B_occ.data(), ri.d_B_occ,
                             sizeof(float) * count, deviceMemcpyDeviceToHost);
            }
        }

        // 构建 D3_eff 和 D2_eff
        std::vector<double> D3_eff;
        QC_Init_D3_J(nao, naux, h_g_vec.data(), h_P.data(), D3_eff);
        for (const auto& channel : spin_channels)
            QC_Accumulate_D3_K_Channel(nao, naux, channel.density.data(),
                                       h_metric_inv_sqrt.data(),
                                       channel.B_occ.data(),
                                       channel.factor.data(),
                                       channel.factor_rank, dft.exx_fraction,
                                       D3_eff);

        std::vector<double> h_Z_K((size_t)naux * naux, 0.0);
        for (const auto& channel : spin_channels)
            QC_Accumulate_Z_K_Stored_Channel(
                nao, naux, h_eri3c.data(), channel.density.data(),
                channel.B_occ.data(), channel.factor.data(),
                channel.factor_rank, h_Z_K);

        std::vector<double> D2_eff;
        QC_Init_D2_J(naux, h_g_vec.data(), D2_eff);
        if (!spin_channels.empty())
            QC_Accumulate_D2_K_DaleckiiKrein(
                naux, ri.naux_eff, h_Z_K.data(), ri.h_eigval.data(),
                ri.h_eigvec.data(), dft.exx_fraction, D2_eff);

        if (!launch_grad_kernels(D2_eff, D3_eff)) return;
    }
}
