#pragma once

static __global__ void QC_Accumulate_Density_Residual_Kernel(
    const int n, const float* density_a, const float* density_a_new,
    const float* density_b, const float* density_b_new, double* squared_norm)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        const double delta_a =
            (double)density_a_new[i] - (double)density_a[i];
        double value = delta_a * delta_a;
        if (density_b != nullptr)
        {
            const double delta_b =
                (double)density_b_new[i] - (double)density_b[i];
            value += delta_b * delta_b;
        }
        atomicAdd(squared_norm, value);
    }
}

static __global__ void QC_Finalize_Density_Residual_Kernel(
    const int matrix_element_count, const int spin_channel_count,
    double* squared_norm)
{
    const double sample_count =
        (double)matrix_element_count * (double)spin_channel_count;
    squared_norm[0] = sqrt(squared_norm[0] / sample_count);
}

bool QUANTUM_CHEMISTRY::Check_Convergence(int iter, int md_step,
                                          double h_energy, double h_delta_e,
                                          bool physical_iteration)
{
    const int nao2 = (int)mol.nao2;

    // Compare P_{n+1}, produced by the just-diagonalized Fock matrix, with the
    // P_n that generated the reported energy and physical (pre-DIIS) Fock.
    // Energy alone can cross a tiny delta accidentally while the density is
    // still oscillating, especially for open-shell DFT.
    deviceMemset(scf_ws.runtime.d_density_residual, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_Accumulate_Density_Residual_Kernel,
        Positive_Int_Ceil_Div(nao2, threads), threads, 0, 0, nao2,
        scf_ws.alpha.d_P, scf_ws.alpha.d_P_new,
        scf_ws.runtime.unrestricted ? scf_ws.beta.d_P : (const float*)nullptr,
        scf_ws.runtime.unrestricted ? scf_ws.beta.d_P_new
                                    : (const float*)nullptr,
        scf_ws.runtime.d_density_residual);
    Launch_Device_Kernel(QC_Finalize_Density_Residual_Kernel, 1, 1, 0, 0,
                         nao2, scf_ws.runtime.unrestricted ? 2 : 1,
                         scf_ws.runtime.d_density_residual);
    double h_density_residual = 0.0;
    deviceMemcpy(&h_density_residual, scf_ws.runtime.d_density_residual,
                 sizeof(double), deviceMemcpyDeviceToHost);
    if (!Double_Memory_Is_Finite(&h_density_residual))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Check_Convergence",
            "Reason:\n    non-finite SCF density residual at MD step %d, "
            "iteration %d\n",
            md_step, iter + 1);
        return false;
    }

    const QC_SCF_Convergence_Decision decision = QC_SCF_Observe_Iteration(
        scf_ws.runtime.convergence, physical_iteration, iter > 0, h_delta_e,
        h_density_residual, scf_ws.runtime.energy_tol,
        scf_ws.runtime.density_tol,
        scf_ws.runtime.configured_level_shift);
    if (decision.reset_diis_history)
    {
        scf_ws.diis.diis_hist_count = scf_ws.diis.diis_hist_head = 0;
        scf_ws.diis.last_enorm = 1.0e10;
    }

    const int required_consecutive_stable_iterations = 1;
    scf_ws.runtime.convergence_streak =
        physical_iteration
            ? scf_ws.runtime.convergence.physical_streak
            : scf_ws.runtime.convergence.accelerated_streak;
    const bool converged = decision.converged;
    if (physical_iteration && !converged && decision.reset_diis_history)
        Start_Ensemble_Probe(h_energy, h_density_residual);
    const int h_converged = converged ? 1 : 0;
    deviceMemcpy(scf_ws.runtime.d_converged, &h_converged, sizeof(int),
                 deviceMemcpyHostToDevice);

    // Only publish P_{n+1} when another iteration will rebuild E and F from
    // it.  On termination, retain P_n so the final E, physical F_for_grad, and
    // density all belong to the same SCF iteration.
    if (!converged)
    {
        deviceMemcpy(scf_ws.alpha.d_P, scf_ws.alpha.d_P_new,
                     sizeof(float) * nao2,
                     deviceMemcpyDeviceToDevice);
        if (scf_ws.runtime.unrestricted)
        {
            deviceMemcpy(scf_ws.beta.d_P, scf_ws.beta.d_P_new,
                         sizeof(float) * nao2, deviceMemcpyDeviceToDevice);
            QC_Add_Matrix(nao2, scf_ws.alpha.d_P, scf_ws.beta.d_P,
                          scf_ws.direct.d_Ptot);
        }
    }

    if (scf_ws.runtime.print_iter && CONTROLLER::MPI_rank == 0)
    {
        const char* fock_build = "full";
        if (scf_ws.runtime.fock_build_mode == QC_SCF_FOCK_BUILD_INCREMENTAL)
            fock_build = "incremental";
        else if (scf_ws.runtime.fock_build_mode == QC_SCF_FOCK_BUILD_RI)
            fock_build = "ri";
        FILE* out = (scf_output_file != NULL) ? scf_output_file : stdout;
        fprintf(out,
                "Step %6d | SCF Iter %3d | E(Ha)=%.12f | dE(Ha)=%+.6e "
                "| dP(rms)=%.6e | map=%s | fock=%s | shift=%.6e "
                "| stable=%d/%d",
                md_step, iter + 1, h_energy, h_delta_e, h_density_residual,
                physical_iteration ? "physical" : "accelerated",
                fock_build, scf_ws.runtime.level_shift,
                scf_ws.runtime.convergence_streak,
                required_consecutive_stable_iterations);
        fprintf(out, "\n");
        fflush(out);
    }

    return converged;
}
