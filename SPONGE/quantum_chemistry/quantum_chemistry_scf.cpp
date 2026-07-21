// clang-format off
#include "integrals/one_e.hpp"
#include "integrals/eri/common/eri_rys.hpp"
#include "quantum_chemistry.h"
#include "scf/accumulate_energy.hpp"
#include "scf/apply_diis.hpp"
#include "scf/build_fock.hpp"
#include "scf/diag_density.hpp"
#include "scf/energy_validation_policy.hpp"
#include "scf/ensemble_converge.hpp"
#include "scf/mix_converge.hpp"
#include "scf/pre_scf.hpp"
#include "scf/workspace.hpp"
#include "structure/matrix.h"
// clang-format on

#ifdef GPU_ARCH_NAME
static __device__ __forceinline__ bool QC_SCF_Float_Is_Finite(float value)
{
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
}
#else
static __host__ __device__ __forceinline__ bool QC_SCF_Float_Is_Finite(
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

static __global__ void QC_Validate_Density_Kernel(int nao2,
                                                  const float* alpha_density,
                                                  const float* beta_density,
                                                  int* failure)
{
    SIMPLE_DEVICE_FOR(i, nao2)
    {
        if (!QC_SCF_Float_Is_Finite(alpha_density[i])) atomicExch(failure, i);
        if (beta_density != NULL && !QC_SCF_Float_Is_Finite(beta_density[i]))
            atomicExch(failure, nao2 + i);
    }
}

void QUANTUM_CHEMISTRY::Solve_SCF(const VECTOR* crd, const VECTOR box_length,
                                  bool need_energy, int md_step,
                                  bool commit_sampling_state,
                                  std::uint64_t coordinate_generation)
{
    if (!is_initialized) return;

    const size_t density_bytes = sizeof(float) * mol.nao2;
    deviceMemcpy(scf_ws.alpha.d_P, d_accepted_alpha_density, density_bytes,
                 deviceMemcpyDeviceToDevice);
    if (scf_ws.runtime.unrestricted)
    {
        deviceMemcpy(scf_ws.beta.d_P, d_accepted_beta_density, density_bytes,
                     deviceMemcpyDeviceToDevice);
        QC_Add_Matrix((int)mol.nao2, scf_ws.alpha.d_P, scf_ws.beta.d_P,
                      scf_ws.direct.d_Ptot);
    }
    const size_t env_bytes = sizeof(float) * mol.h_env.size();
    deviceMemcpy(mol.d_env, d_accepted_env, env_bytes,
                 deviceMemcpyDeviceToDevice);
    if (scf_ws.ri.enabled)
    {
        const size_t aux_env_bytes = sizeof(float) * scf_ws.ri.h_aux_env.size();
        deviceMemcpy(scf_ws.ri.d_aux_env, d_accepted_aux_env, aux_env_bytes,
                     deviceMemcpyDeviceToDevice);
    }
    need_initial_guess = accepted_need_initial_guess;
    const bool restored_density_is_ensemble =
        accepted_density_is_ensemble;

    if (coordinate_generation == accepted_coordinate_generation)
        Refresh_Coordinate_Derived_State_From_Env();
    else if (!Update_Coordinates_From_MD(crd, box_length))
        return;
    if (!Validate_Nuclear_Geometry(box_length)) return;
    if (dft.enable_dft) Update_DFT_Grid();

    Reset_SCF_State();
    if (restored_density_is_ensemble)
    {
        // A fractional ensemble density is not an idempotent occupied-space
        // projector, so S-SPS/occ is not its virtual projector and must never
        // be used for level shifting.  Start with one unshifted raw F[P] map;
        // it either proves an ordinary fixed point or immediately resumes the
        // ensemble active-set solve from the accepted density.
        scf_ws.runtime.convergence.level_shift_stage =
            QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT;
        scf_ws.runtime.convergence.verifying_physical_fixed_point = true;
        scf_ws.runtime.level_shift = 0.0;
    }

    // 解析计算 norms（不依赖 1e 积分的 S 矩阵）
    Compute_Analytical_Norms();
    Compute_OneE_Integrals();
    Compute_ECP_Matrix();
    if (need_energy) Compute_Nuclear_Repulsion(box_length);
    Prepare_Integrals();
    Build_Shell_Pair_Bounds();
    if (scf_ws.ri.enabled) RI_Precompute();
    if (!Build_Overlap_X()) return;

    if (need_initial_guess)
    {
        if (!Build_Initial_Guess()) return;
        need_initial_guess = false;
    }

    // DFT gets a short shifted warmup for the SAP-to-KS transition.  Candidate
    // generation keeps the configured stabilization until an unshifted,
    // non-DIIS physical map verifies the fixed point.
    const int dft_warmup = dft.enable_dft ? 3 : 0;

    bool converged = false;
    for (int iter = 0; iter < scf_ws.runtime.max_scf_iter; ++iter)
    {
        const bool ensemble_iteration =
            scf_ws.ensemble.phase != QC_SCF_ENSEMBLE_INACTIVE;
        const bool physical_iteration =
            !ensemble_iteration &&
            scf_ws.runtime.convergence.verifying_physical_fixed_point;
        Build_Fock(iter, physical_iteration || ensemble_iteration);
        Accumulate_SCF_Energy(iter);
        double iteration_energy = 0.0;
        double delta_energy = 0.0;
        deviceMemcpy(&iteration_energy, scf_ws.core.d_scf_energy,
                     sizeof(double), deviceMemcpyDeviceToHost);
        deviceMemcpy(&delta_energy, scf_ws.runtime.d_delta_e, sizeof(double),
                     deviceMemcpyDeviceToHost);
        const bool finite_energy_observation =
            QC_SCF_Require_Finite_Energy_Observation(
                md_step, iter + 1, iter > 0, iteration_energy, delta_energy,
                [&](const QC_SCF_Energy_Validation_Failure& failure)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorSimulationBreakDown,
                        "QUANTUM_CHEMISTRY::Solve_SCF",
                        "Reason:\n    non-finite SCF %s at MD step %d, "
                        "iteration %d: %.17g Hartree\n",
                        failure.quantity_name, failure.md_step,
                        failure.iteration, failure.value);
                });
        if (!finite_energy_observation)
            return;

        // 缓存 DIIS 前的 Fock 供梯度使用（避免梯度中重建 Fock）
        if (need_gradient && scf_ws.alpha.d_F_for_grad)
        {
            deviceMemcpy(scf_ws.alpha.d_F_for_grad, scf_ws.alpha.d_F_double,
                         sizeof(double) * mol.nao2, deviceMemcpyDeviceToDevice);
            if (scf_ws.runtime.unrestricted && scf_ws.beta.d_F_for_grad)
                deviceMemcpy(scf_ws.beta.d_F_for_grad, scf_ws.beta.d_F_double,
                             sizeof(double) * mol.nao2,
                             deviceMemcpyDeviceToDevice);
        }

        // Select the current homotopy stage before extrapolation.  While a
        // positive shift is active, Apply_DIIS uses only bounded EDIIS/ADIIS;
        // unconstrained CDIIS is enabled after continuation reaches zero.
        scf_ws.runtime.level_shift = QC_SCF_Level_Shift_For_Map(
            scf_ws.runtime.configured_level_shift,
            scf_ws.runtime.convergence,
            physical_iteration || ensemble_iteration);
        if (!physical_iteration && !ensemble_iteration &&
            !(dft.enable_dft && iter < dft_warmup))
            Apply_DIIS(iter);

        if (!Diagonalize_And_Build_Density()) return;
        if (ensemble_iteration)
        {
            const int ensemble_result = Advance_Ensemble_Search(
                iter, md_step, iteration_energy, delta_energy);
            if (ensemble_result == 2)
            {
                converged = true;
                break;
            }
            if (ensemble_result == 1) continue;
        }
        if (Check_Convergence(iter, md_step, iteration_energy, delta_energy,
                              physical_iteration))
        {
            converged = true;
            break;
        }
    }

    double final_energy = 0.0;
    deviceMemcpy(&final_energy, scf_ws.core.d_scf_energy, sizeof(double),
                 deviceMemcpyDeviceToHost);
    const bool finite_final_energy = QC_SCF_Require_Finite_Energy_Observation(
        md_step, scf_ws.runtime.max_scf_iter, false, final_energy, 0.0,
        [&](const QC_SCF_Energy_Validation_Failure& failure)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Solve_SCF",
                "Reason:\n    non-finite final SCF %s at MD step %d after "
                "at most %d iteration(s): %.17g Hartree\n",
                failure.quantity_name, failure.md_step,
                scf_ws.runtime.max_scf_iter, failure.value);
        });
    if (!finite_final_energy)
        return;

    deviceMemset(d_scf_validation_failure, -1, sizeof(int));
    const int validation_threads = 256;
    Launch_Device_Kernel(
        QC_Validate_Density_Kernel,
        Positive_Int_Ceil_Div(static_cast<int>(mol.nao2), validation_threads),
        validation_threads, 0, 0, (int)mol.nao2, scf_ws.alpha.d_P,
        scf_ws.runtime.unrestricted ? scf_ws.beta.d_P : (const float*)NULL,
        d_scf_validation_failure);
    int validation_failure = -1;
    deviceMemcpy(&validation_failure, d_scf_validation_failure, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (validation_failure >= 0)
    {
        const bool beta_channel = validation_failure >= (int)mol.nao2;
        const int matrix_index =
            validation_failure - (beta_channel ? (int)mol.nao2 : 0);
        float invalid_density = 0.0f;
        const float* density =
            beta_channel ? scf_ws.beta.d_P : scf_ws.alpha.d_P;
        deviceMemcpy(&invalid_density, density + matrix_index, sizeof(float),
                     deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "QUANTUM_CHEMISTRY::Solve_SCF",
            "Reason:\n    non-finite %s SCF density at MD step %d, matrix "
            "index %d: %.9g\n",
            beta_channel ? "beta" : "alpha", md_step, matrix_index,
            (double)invalid_density);
        return;
    }
    if (!converged)
    {
        double final_delta_energy = 0.0;
        double final_density_residual = 0.0;
        deviceMemcpy(&final_delta_energy, scf_ws.runtime.d_delta_e,
                     sizeof(double), deviceMemcpyDeviceToHost);
        deviceMemcpy(&final_density_residual,
                     scf_ws.runtime.d_density_residual, sizeof(double),
                     deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "QUANTUM_CHEMISTRY::Solve_SCF",
            "Reason:\n    SCF failed to converge at MD step %d within %d "
            "iteration(s); last finite energy is %.17g Hartree, dE is "
            "%+.9g Hartree, and density RMS residual is %.9g\n",
            md_step, scf_ws.runtime.max_scf_iter, final_energy,
            final_delta_energy, final_density_residual);
        return;
    }

    if (commit_sampling_state)
    {
        deviceMemcpy(d_accepted_alpha_density, scf_ws.alpha.d_P, density_bytes,
                     deviceMemcpyDeviceToDevice);
        if (scf_ws.runtime.unrestricted)
        {
            deviceMemcpy(d_accepted_beta_density, scf_ws.beta.d_P,
                         density_bytes, deviceMemcpyDeviceToDevice);
        }
        deviceMemcpy(d_accepted_env, mol.d_env, env_bytes,
                     deviceMemcpyDeviceToDevice);
        if (scf_ws.ri.enabled)
        {
            const size_t aux_env_bytes =
                sizeof(float) * scf_ws.ri.h_aux_env.size();
            deviceMemcpy(d_accepted_aux_env, scf_ws.ri.d_aux_env, aux_env_bytes,
                         deviceMemcpyDeviceToDevice);
        }
        accepted_need_initial_guess = need_initial_guess;
        accepted_density_is_ensemble = scf_ws.ensemble.confirmed;
        accepted_coordinate_generation = coordinate_generation;
    }
}

static __global__ void QC_Accumulate_Total_Energy_Kernel(
    float* local_atom_energy, int local_atom, float energy)
{
#ifdef USE_GPU
    if (blockIdx.x != 0 || threadIdx.x != 0) return;
#endif
    atomicAdd(local_atom_energy + local_atom, energy);
}

void QUANTUM_CHEMISTRY::Accumulate_Energy(float* local_atom_energy,
                                          const int* global_to_local,
                                          int owned_atom_numbers)
{
    if (!is_initialized) return;
    if (atom_local.empty())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Accumulate_Energy",
            "Reason:\n    an initialized QC calculation has no QC atoms\n");
        return;
    }

    const int anchor_global_atom = atom_local.front();
    int anchor_local_atom = -1;
    deviceMemcpy(&anchor_local_atom, global_to_local + anchor_global_atom,
                 sizeof(int), deviceMemcpyDeviceToHost);
    const int owns_anchor =
        anchor_local_atom >= 0 && anchor_local_atom < owned_atom_numbers;
    int owner_count = owns_anchor;
#ifdef USE_MPI
    if (CONTROLLER::PP_MPI_size > 1)
    {
        MPI_Allreduce(&owns_anchor, &owner_count, 1, MPI_INT, MPI_SUM,
                      CONTROLLER::pp_comm);
    }
#endif
    if (owner_count != 1)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Accumulate_Energy",
            "Reason:\n    QC energy anchor global atom %d has %d owning PP "
            "ranks; exactly one is required\n",
            anchor_global_atom, owner_count);
        return;
    }

    double energy_hartree = 0.0;
    deviceMemcpy(&energy_hartree, scf_ws.core.d_scf_energy, sizeof(double),
                 deviceMemcpyDeviceToHost);
    const double energy_kcal_mol =
        energy_hartree * (double)CONSTANT_HARTREE_TO_KCAL_MOL;
    const float energy_to_accumulate = (float)energy_kcal_mol;
    if (!Double_Memory_Is_Finite(&energy_hartree) ||
        !Double_Memory_Is_Finite(&energy_kcal_mol) ||
        !Float_Memory_Is_Finite(&energy_to_accumulate))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Accumulate_Energy",
            "Reason:\n    QC energy %.17g Hartree cannot be represented as "
            "a finite SPONGE potential in kcal/mol\n",
            energy_hartree);
        return;
    }

    if (owns_anchor)
    {
        Launch_Device_Kernel(QC_Accumulate_Total_Energy_Kernel, 1, 1, 0, 0,
                             local_atom_energy, anchor_local_atom,
                             energy_to_accumulate);
    }
}

void QUANTUM_CHEMISTRY::Compute_Spin_Square()
{
    const int nao = mol.nao;
    const int nao2 = mol.nao2;

    // <S²> = s(s+1) + N_beta - Tr(P_alpha · S · P_beta · S)
    // 使用 ortho 的 double workspace 作为临时缓冲，避免污染 Fock 矩阵
    double* d_tmp1 = scf_ws.ortho.d_dwork_nao2_1;
    double* d_tmp2 = scf_ws.ortho.d_dwork_nao2_2;
    double* d_tmp3 = scf_ws.ortho.d_dwork_nao2_3;

    // 提升到 double: dPa = P_alpha, dS = S
    QC_Float_To_Double(nao2, scf_ws.alpha.d_P, d_tmp1);
    QC_Float_To_Double(nao2, scf_ws.core.d_S, d_tmp2);

    // d_tmp3 = P_alpha * S
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_tmp1, nao, d_tmp2, nao, d_tmp3,
                nao);

    // d_tmp1 = P_beta (提升)
    QC_Float_To_Double(nao2, scf_ws.beta.d_P, d_tmp1);

    // d_tmp1 = (P_alpha * S) * P_beta -> 复用: d_tmp4 借用 d_dwork_nao2_4
    double* d_tmp4 = scf_ws.ortho.d_dwork_nao2_4;
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_tmp3, nao, d_tmp1, nao, d_tmp4,
                nao);

    // Tr(P_alpha * S * P_beta * S) = Σ_ij (P_alpha·S·P_beta)_ij * S_ij
    double trace = 0.0;
    double* d_accum = scf_ws.diis.d_diis_accum;
    if (d_accum == NULL)
    {
        Device_Malloc_Safely((void**)&d_accum, sizeof(double));
    }
    deviceMemset(d_accum, 0, sizeof(double));
    QC_Double_Dot(nao2, d_tmp4, d_tmp2, d_accum);
    deviceMemcpy(&trace, d_accum, sizeof(double), deviceMemcpyDeviceToHost);
    if (scf_ws.diis.d_diis_accum == NULL)
    {
        deviceFree(d_accum);
    }

    double s = 0.5 * (scf_ws.runtime.n_alpha - scf_ws.runtime.n_beta);
    scf_ws.runtime.spin_square_exact = s * (s + 1.0);
    scf_ws.runtime.spin_square =
        scf_ws.runtime.spin_square_exact + scf_ws.runtime.n_beta - trace;
}
