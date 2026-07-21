#pragma once

static __global__ void QC_SCF_Store_Ensemble_Line_Kernel(
    const int n, const float* density, const float* density_new, float* origin,
    float* direction)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        origin[i] = density[i];
        direction[i] = density_new[i] - density[i];
    }
}

static __global__ void QC_SCF_Build_Spectral_Orbital_Direction_Kernel(
    const int n, const float* origin, const float* corrected,
    float* direction)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        direction[i] = corrected[i] - origin[i];
    }
}

static __global__ void QC_SCF_Set_Stored_Ensemble_Line_Trial_Kernel(
    const int n, const float* origin, const float* direction,
    const double fraction, float* density)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        density[i] = static_cast<float>(
            static_cast<double>(origin[i]) +
            fraction * static_cast<double>(direction[i]));
    }
}

static __global__ void QC_SCF_Build_Double_Ensemble_Line_Trial_Kernel(
    const int n, const float* origin, const float* direction,
    const double fraction, double* density)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        density[i] = static_cast<double>(origin[i]) +
                     fraction * static_cast<double>(direction[i]);
    }
}

static __global__ void QC_SCF_Accumulate_Fock_Trial_Change_Kernel(
    const int n, const double* origin_fock, const float* origin_density,
    const float* trial_density, double* result)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        atomicAdd(result,
                  origin_fock[i] *
                      (static_cast<double>(trial_density[i]) -
                       static_cast<double>(origin_density[i])));
    }
}

static __global__ void QC_SCF_Accumulate_Particle_Count_Kernel(
    const int n, const float* overlap, const float* density, double* result)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        atomicAdd(result, static_cast<double>(overlap[i]) *
                              static_cast<double>(density[i]));
    }
}

static __global__ void QC_SCF_Scale_Orbitals_By_Occupation_Kernel(
    const int matrix_element_count, const int orbital_count,
    const double* orbitals, const double* occupations,
    double* scaled_orbitals)
{
    SIMPLE_DEVICE_FOR(i, matrix_element_count)
    {
        scaled_orbitals[i] =
            orbitals[i] * occupations[i % orbital_count];
    }
}

static __global__ void QC_SCF_Extract_Diagonal_Kernel(
    const int dimension, const double* matrix, double* diagonal)
{
    SIMPLE_DEVICE_FOR(i, dimension)
    {
        diagonal[i] = matrix[(long long)i * dimension + i];
    }
}

static __global__ void QC_SCF_Accumulate_Fock_Direction_Kernel(
    const int n, const double* fock, const float* direction, double* result)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        atomicAdd(result, fock[i] * static_cast<double>(direction[i]));
    }
}

static __global__ void QC_SCF_Accumulate_Ensemble_Direction_Norm_Kernel(
    const int n, const float* direction_alpha,
    const float* direction_beta, double* result)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        const double alpha = static_cast<double>(direction_alpha[i]);
        double value = alpha * alpha;
        if (direction_beta != NULL)
        {
            const double beta = static_cast<double>(direction_beta[i]);
            value += beta * beta;
        }
        atomicAdd(result, value);
    }
}

static __global__ void QC_SCF_Accumulate_Ensemble_Difference_Norm_Kernel(
    const int n, const float* density_alpha, const float* atom_alpha,
    const float* density_beta, const float* atom_beta, double* result)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        const double alpha = static_cast<double>(density_alpha[i]) -
                             static_cast<double>(atom_alpha[i]);
        double value = alpha * alpha;
        if (density_beta != NULL)
        {
            const double beta = static_cast<double>(density_beta[i]) -
                                static_cast<double>(atom_beta[i]);
            value += beta * beta;
        }
        atomicAdd(result, value);
    }
}

static __global__ void QC_SCF_Finalize_Ensemble_Norm_Kernel(
    int matrix_element_count, int spin_channel_count, double* result)
{
    result[0] = sqrt(result[0] /
                     (static_cast<double>(matrix_element_count) *
                      static_cast<double>(spin_channel_count)));
}

static void QC_SCF_Store_Ensemble_Line(QC_SCF_WORKSPACE& workspace,
                                       int matrix_element_count)
{
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Store_Ensemble_Line_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, workspace.alpha.d_P, workspace.alpha.d_P_new,
        workspace.ensemble.d_origin_alpha,
        workspace.ensemble.d_direction_alpha);
    if (workspace.runtime.unrestricted)
        Launch_Device_Kernel(
            QC_SCF_Store_Ensemble_Line_Kernel,
            Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0,
            0, matrix_element_count, workspace.beta.d_P,
            workspace.beta.d_P_new, workspace.ensemble.d_origin_beta,
            workspace.ensemble.d_direction_beta);
}

static void QC_SCF_Set_Stored_Ensemble_Line_Trial(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count, double fraction)
{
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Set_Stored_Ensemble_Line_Trial_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, workspace.ensemble.d_origin_alpha,
        workspace.ensemble.d_direction_alpha, fraction,
        workspace.alpha.d_P);
    if (workspace.runtime.unrestricted)
    {
        Launch_Device_Kernel(
            QC_SCF_Set_Stored_Ensemble_Line_Trial_Kernel,
            Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0,
            0, matrix_element_count, workspace.ensemble.d_origin_beta,
            workspace.ensemble.d_direction_beta, fraction,
            workspace.beta.d_P);
        QC_Add_Matrix(matrix_element_count, workspace.alpha.d_P,
                      workspace.beta.d_P, workspace.direct.d_Ptot);
    }
}

static bool QC_SCF_Store_Ensemble_Corrective_Line(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count)
{
    QC_SCF_Ensemble_Workspace& ensemble = workspace.ensemble;
    if (ensemble.line_weight_direction.size() !=
        ensemble.active_atoms.size())
        return false;

    double* accumulator = workspace.ortho.d_dwork_nao2_1;
    double* atom_double = workspace.ortho.d_dwork_nao2_2;
    auto store_channel = [&](bool beta)
    {
        const float* density =
            beta ? workspace.beta.d_P : workspace.alpha.d_P;
        float* origin = beta ? ensemble.d_origin_beta
                             : ensemble.d_origin_alpha;
        float* direction = beta ? ensemble.d_direction_beta
                                : ensemble.d_direction_alpha;
        deviceMemcpy(origin, density, sizeof(float) * matrix_element_count,
                     deviceMemcpyDeviceToDevice);
        deviceMemset(accumulator, 0,
                     sizeof(double) * matrix_element_count);
        for (size_t i = 0; i < ensemble.active_atoms.size(); ++i)
        {
            const double coefficient = ensemble.line_weight_direction[i];
            if (!QC_SCF_Ensemble_Double_Is_Finite(coefficient)) return false;
            if (coefficient == 0.0) continue;
            const QC_SCF_Ensemble_Atom& atom = ensemble.active_atoms[i];
            QC_Float_To_Double(matrix_element_count,
                               beta ? atom.d_beta : atom.d_alpha,
                               atom_double);
            QC_Double_Axpy(matrix_element_count, coefficient, atom_double,
                           accumulator);
        }
        QC_Double_To_Float(matrix_element_count, accumulator, direction);
        return true;
    };

    if (!store_channel(false)) return false;
    return !workspace.runtime.unrestricted || store_channel(true);
}

static bool QC_SCF_Set_Ensemble_Trial(QC_SCF_WORKSPACE& workspace,
                                      int matrix_element_count,
    double fraction)
{
    QC_SCF_Ensemble_Workspace& ensemble = workspace.ensemble;
    if (ensemble.line_weight_direction.size() !=
            ensemble.active_atoms.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(fraction) || fraction < 0.0 ||
        fraction > ensemble.maximum_fraction + 1.0e-12)
        return false;
    std::vector<double> origin_weights;
    origin_weights.reserve(ensemble.active_atoms.size());
    for (const QC_SCF_Ensemble_Atom& atom : ensemble.active_atoms)
        origin_weights.push_back(atom.weight);
    std::vector<double> validated_trial_weights;
    if (!QC_SCF_Build_Ensemble_Line_Weights(
            origin_weights, ensemble.line_weight_direction, fraction,
            validated_trial_weights))
        return false;

    // The signed derivative is F:stored_direction, so every raw sample must
    // lie on that very same represented line.  Reconstructing the sample from
    // active weights performs an independent rounded sum; near degeneracy its
    // few-ULP chord error is amplified by small fractions and cancellation in
    // the energy components.  Use the exact committed float P as the origin,
    // apply origin+x*stored_direction in double, then cast only once.  The
    // validated weights remain the transaction metadata committed after the
    // raw sample is accepted; the sampled P itself is retained verbatim.
    QC_SCF_Set_Stored_Ensemble_Line_Trial(workspace, matrix_element_count,
                                          fraction);
    return true;
}

static double QC_SCF_Ensemble_Directional_Derivative(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count)
{
    deviceMemset(workspace.ensemble.d_accum, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Accumulate_Fock_Direction_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, workspace.alpha.d_F_double,
        workspace.ensemble.d_direction_alpha, workspace.ensemble.d_accum);
    if (workspace.runtime.unrestricted)
        Launch_Device_Kernel(
            QC_SCF_Accumulate_Fock_Direction_Kernel,
            Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0,
            0, matrix_element_count, workspace.beta.d_F_double,
            workspace.ensemble.d_direction_beta,
            workspace.ensemble.d_accum);
    double derivative = 0.0;
    deviceMemcpy(&derivative, workspace.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return derivative;
}

static double QC_SCF_Ensemble_Direction_RMS(QC_SCF_WORKSPACE& workspace,
                                            int matrix_element_count)
{
    deviceMemset(workspace.ensemble.d_accum, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Accumulate_Ensemble_Direction_Norm_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, workspace.ensemble.d_direction_alpha,
        workspace.runtime.unrestricted
            ? workspace.ensemble.d_direction_beta
            : static_cast<const float*>(NULL),
        workspace.ensemble.d_accum);
    Launch_Device_Kernel(QC_SCF_Finalize_Ensemble_Norm_Kernel, 1, 1, 0, 0,
                         matrix_element_count,
                         workspace.runtime.unrestricted ? 2 : 1,
                         workspace.ensemble.d_accum);
    double result = 0.0;
    deviceMemcpy(&result, workspace.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return result;
}

static double QC_SCF_Ensemble_Particle_Count(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count,
    const float* density)
{
    deviceMemset(workspace.ensemble.d_accum, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Accumulate_Particle_Count_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, workspace.core.d_S, density,
        workspace.ensemble.d_accum);
    double result = 0.0;
    deviceMemcpy(&result, workspace.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return result;
}

static double QC_SCF_Spectral_Origin_Fock_Trial_Change(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count,
    const double* origin_fock, const float* origin_density,
    const float* trial_density)
{
    deviceMemset(workspace.ensemble.d_accum, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Accumulate_Fock_Trial_Change_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, origin_fock, origin_density, trial_density,
        workspace.ensemble.d_accum);
    double result = 0.0;
    deviceMemcpy(&result, workspace.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return result;
}

static double QC_SCF_Spectral_Origin_Fock_Direction(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count,
    const double* origin_fock, const float* direction)
{
    deviceMemset(workspace.ensemble.d_accum, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Accumulate_Fock_Direction_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, origin_fock, direction,
        workspace.ensemble.d_accum);
    double result = 0.0;
    deviceMemcpy(&result, workspace.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return result;
}

static bool QC_SCF_Set_Spectral_Orbital_Trial(
    QC_SCF_WORKSPACE& workspace, int matrix_dimension,
    int matrix_element_count, double fraction)
{
    QC_SCF_Ensemble_Workspace& ensemble = workspace.ensemble;
    if (matrix_dimension <= 0 ||
        matrix_element_count != matrix_dimension * matrix_dimension ||
        !QC_SCF_Ensemble_Double_Is_Finite(fraction) || !(fraction > 0.0) ||
        fraction > 1.0 || ensemble.d_origin_alpha == NULL ||
        ensemble.d_direction_alpha == NULL ||
        ensemble.d_spectral_origin_fock_alpha == NULL)
        return false;

    const double channel_rounding_budget =
        workspace.runtime.energy_tol /
        (workspace.runtime.unrestricted ? 2.0 : 1.0);
    if (!QC_SCF_Ensemble_Double_Is_Finite(channel_rounding_budget) ||
        !(channel_rounding_budget > 0.0))
        return false;

    double* double_target = workspace.ortho.d_dwork_nao2_1;
    const int threads = 256;
    auto set_channel = [&](const float* origin, const float* direction,
                           const double* origin_fock,
                           float* trial_density) -> bool
    {
        if (origin == NULL || direction == NULL || origin_fock == NULL ||
            trial_density == NULL)
            return false;
        Launch_Device_Kernel(
            QC_SCF_Build_Double_Ensemble_Line_Trial_Kernel,
            Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0,
            0, matrix_element_count, origin, direction, fraction,
            double_target);
        QC_Double_To_Float(matrix_element_count, double_target,
                           trial_density);
        double fixed_fock_change =
            QC_SCF_Spectral_Origin_Fock_Trial_Change(
                workspace, matrix_element_count, origin_fock, origin,
                trial_density);
        const double ideal_fixed_fock_change =
            fraction * QC_SCF_Spectral_Origin_Fock_Direction(
                           workspace, matrix_element_count, origin_fock,
                           direction);
        if (!QC_SCF_Ensemble_Double_Is_Finite(fixed_fock_change) ||
            !QC_SCF_Ensemble_Double_Is_Finite(ideal_fixed_fock_change))
            return false;

        // Nearest rounding is unbiased and remains the normal path.  Only a
        // represented fraction whose origin-F objective falls behind its
        // exact affine target by more than the per-channel energy budget
        // receives the symmetric one-ULP correction.  Comparing with zero
        // would miss a several-microhartree rounding error whenever the
        // represented trial still happens to be a net descent step.
        double rounding_error =
            fixed_fock_change - ideal_fixed_fock_change;
        if (rounding_error > channel_rounding_budget)
        {
            QC_Round_Symmetric_Double_Matrix_For_Nonincreasing_Linear_Objective(
                matrix_dimension, double_target, origin_fock,
                trial_density);
            fixed_fock_change =
                QC_SCF_Spectral_Origin_Fock_Trial_Change(
                workspace, matrix_element_count, origin_fock, origin,
                trial_density);
            rounding_error =
                fixed_fock_change - ideal_fixed_fock_change;
            if (!QC_SCF_Ensemble_Double_Is_Finite(fixed_fock_change) ||
                !QC_SCF_Ensemble_Double_Is_Finite(rounding_error) ||
                rounding_error > channel_rounding_budget)
                return false;
        }

        const double origin_particles = QC_SCF_Ensemble_Particle_Count(
            workspace, matrix_element_count, origin);
        const double trial_particles = QC_SCF_Ensemble_Particle_Count(
            workspace, matrix_element_count, trial_density);
        return QC_SCF_Ensemble_Double_Is_Finite(origin_particles) &&
               QC_SCF_Ensemble_Double_Is_Finite(trial_particles) &&
               origin_particles >= 0.0 && trial_particles >= 0.0 &&
               std::fabs(trial_particles - origin_particles) <=
                   workspace.runtime.density_tol *
                       std::max(1.0, origin_particles);
    };

    if (!set_channel(ensemble.d_origin_alpha,
                     ensemble.d_direction_alpha,
                     ensemble.d_spectral_origin_fock_alpha,
                     workspace.alpha.d_P))
        return false;
    if (workspace.runtime.unrestricted)
    {
        if (!set_channel(ensemble.d_origin_beta,
                         ensemble.d_direction_beta,
                         ensemble.d_spectral_origin_fock_beta,
                         workspace.beta.d_P))
            return false;
        QC_Add_Matrix(matrix_element_count, workspace.alpha.d_P,
                      workspace.beta.d_P, workspace.direct.d_Ptot);
    }
    return true;
}

static double QC_SCF_Ensemble_Fock_Density_Inner_Product(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count,
    const float* density_alpha, const float* density_beta)
{
    deviceMemset(workspace.ensemble.d_accum, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Accumulate_Fock_Direction_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, workspace.alpha.d_F_double, density_alpha,
        workspace.ensemble.d_accum);
    if (workspace.runtime.unrestricted)
        Launch_Device_Kernel(
            QC_SCF_Accumulate_Fock_Direction_Kernel,
            Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0,
            0, matrix_element_count, workspace.beta.d_F_double, density_beta,
            workspace.ensemble.d_accum);
    double result = 0.0;
    deviceMemcpy(&result, workspace.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return result;
}

static double QC_SCF_Ensemble_Density_Difference_RMS(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count,
    const float* density_alpha, const float* atom_alpha,
    const float* density_beta, const float* atom_beta)
{
    deviceMemset(workspace.ensemble.d_accum, 0, sizeof(double));
    const int threads = 256;
    Launch_Device_Kernel(
        QC_SCF_Accumulate_Ensemble_Difference_Norm_Kernel,
        Positive_Int_Ceil_Div(matrix_element_count, threads), threads, 0, 0,
        matrix_element_count, density_alpha, atom_alpha,
        workspace.runtime.unrestricted ? density_beta : (const float*)NULL,
        workspace.runtime.unrestricted ? atom_beta : (const float*)NULL,
        workspace.ensemble.d_accum);
    Launch_Device_Kernel(QC_SCF_Finalize_Ensemble_Norm_Kernel, 1, 1, 0, 0,
                         matrix_element_count,
                         workspace.runtime.unrestricted ? 2 : 1,
                         workspace.ensemble.d_accum);
    double result = 0.0;
    deviceMemcpy(&result, workspace.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return result;
}

static bool QC_SCF_Evaluate_Ensemble_Active_Set_FW_Gap(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count, double& gap)
{
    const QC_SCF_Ensemble_Workspace& ensemble = workspace.ensemble;
    std::vector<double> linear_values;
    std::vector<double> weights;
    linear_values.reserve(ensemble.active_atoms.size());
    weights.reserve(ensemble.active_atoms.size());
    for (const QC_SCF_Ensemble_Atom& atom : ensemble.active_atoms)
    {
        linear_values.push_back(QC_SCF_Ensemble_Fock_Density_Inner_Product(
            workspace, matrix_element_count, atom.d_alpha, atom.d_beta));
        weights.push_back(atom.weight);
    }
    return QC_SCF_Ensemble_Active_Set_FW_Gap(linear_values, weights, gap);
}

static void QC_SCF_Clear_Ensemble_Active_Set(QC_SCF_WORKSPACE& workspace)
{
    for (QC_SCF_Ensemble_Atom& atom : workspace.ensemble.active_atoms)
    {
        if (atom.d_alpha != NULL) deviceFree(atom.d_alpha);
        if (atom.d_beta != NULL) deviceFree(atom.d_beta);
    }
    workspace.ensemble.active_atoms.clear();
    workspace.ensemble.line_weight_direction.clear();
    workspace.ensemble.corrective_inverse_hessian.clear();
    workspace.ensemble.corrective_reference_weights.clear();
    workspace.ensemble.corrective_reference_linear_values.clear();
    workspace.ensemble.corrective_line_used_quasi_newton = false;
    workspace.ensemble.corrective_line_used_pairwise = false;
    workspace.ensemble.corrective_force_pairwise = false;
}

static int QC_SCF_Add_Ensemble_Atom(QC_SCF_WORKSPACE& workspace,
                                    int matrix_element_count,
                                    const float* density_alpha,
                                    const float* density_beta, double weight)
{
    if (!QC_SCF_Ensemble_Double_Is_Finite(weight) || !(weight >= 0.0))
        return -1;
    QC_SCF_Ensemble_Atom atom;
    Device_Malloc_Safely((void**)&atom.d_alpha,
                         sizeof(float) * matrix_element_count);
    deviceMemcpy(atom.d_alpha, density_alpha,
                 sizeof(float) * matrix_element_count,
                 deviceMemcpyDeviceToDevice);
    if (workspace.runtime.unrestricted)
    {
        Device_Malloc_Safely((void**)&atom.d_beta,
                             sizeof(float) * matrix_element_count);
        deviceMemcpy(atom.d_beta, density_beta,
                     sizeof(float) * matrix_element_count,
                     deviceMemcpyDeviceToDevice);
    }
    atom.weight = weight;
    workspace.ensemble.active_atoms.push_back(atom);
    return static_cast<int>(workspace.ensemble.active_atoms.size()) - 1;
}

static int QC_SCF_Find_Ensemble_Atom(QC_SCF_WORKSPACE& workspace,
                                     int matrix_element_count,
                                     const float* density_alpha,
                                     const float* density_beta)
{
    for (size_t i = 0; i < workspace.ensemble.active_atoms.size(); ++i)
    {
        const QC_SCF_Ensemble_Atom& atom =
            workspace.ensemble.active_atoms[i];
        const double difference = QC_SCF_Ensemble_Density_Difference_RMS(
            workspace, matrix_element_count, density_alpha, atom.d_alpha,
            density_beta, atom.d_beta);
        if (QC_SCF_Ensemble_Density_Is_Same_Vertex(difference))
            return static_cast<int>(i);
    }
    return -1;
}

static void QC_SCF_Reconstruct_Ensemble_Density(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count)
{
    double* accumulator = workspace.ortho.d_dwork_nao2_1;
    double* atom_double = workspace.ortho.d_dwork_nao2_2;
    auto reconstruct_channel = [&](bool beta)
    {
        deviceMemset(accumulator, 0,
                     sizeof(double) * matrix_element_count);
        for (const QC_SCF_Ensemble_Atom& atom :
             workspace.ensemble.active_atoms)
        {
            QC_Float_To_Double(matrix_element_count,
                               beta ? atom.d_beta : atom.d_alpha,
                               atom_double);
            QC_Double_Axpy(matrix_element_count, atom.weight, atom_double,
                           accumulator);
        }
        QC_Double_To_Float(matrix_element_count, accumulator,
                           beta ? workspace.beta.d_P
                                : workspace.alpha.d_P);
    };
    reconstruct_channel(false);
    if (workspace.runtime.unrestricted)
    {
        reconstruct_channel(true);
        QC_Add_Matrix(matrix_element_count, workspace.alpha.d_P,
                      workspace.beta.d_P, workspace.direct.d_Ptot);
    }
}

static bool QC_SCF_Prune_Ensemble_Zero_Weight_Atoms(
    QC_SCF_WORKSPACE& workspace)
{
    QC_SCF_Ensemble_Workspace& ensemble = workspace.ensemble;
    size_t positive_count = 0;
    for (const QC_SCF_Ensemble_Atom& atom : ensemble.active_atoms)
        if (atom.weight > 0.0) ++positive_count;
    if (positive_count == 0) return false;
    if (positive_count == ensemble.active_atoms.size()) return true;

    std::vector<QC_SCF_Ensemble_Atom> positive_atoms;
    positive_atoms.reserve(positive_count);
    for (QC_SCF_Ensemble_Atom& atom : ensemble.active_atoms)
    {
        if (atom.weight > 0.0)
            positive_atoms.push_back(atom);
        else
        {
            if (atom.d_alpha != NULL) deviceFree(atom.d_alpha);
            if (atom.d_beta != NULL) deviceFree(atom.d_beta);
        }
    }
    ensemble.active_atoms.swap(positive_atoms);
    QC_SCF_Initialize_Corrective_Inverse_Hessian(
        ensemble.active_atoms.size(), ensemble.corrective_inverse_hessian);
    ensemble.corrective_reference_weights.clear();
    ensemble.corrective_reference_linear_values.clear();
    return true;
}

static bool QC_SCF_Commit_Ensemble_Corrective_Weights(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count, double fraction)
{
    QC_SCF_Ensemble_Workspace& ensemble = workspace.ensemble;
    if (ensemble.line_weight_direction.size() !=
            ensemble.active_atoms.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(fraction) || fraction < 0.0 ||
        fraction > ensemble.maximum_fraction + 1.0e-12)
        return false;

    std::vector<double> origin_weights;
    origin_weights.reserve(ensemble.active_atoms.size());
    for (const QC_SCF_Ensemble_Atom& atom : ensemble.active_atoms)
        origin_weights.push_back(atom.weight);
    std::vector<double> committed_weights;
    if (!QC_SCF_Build_Ensemble_Line_Weights(
            origin_weights, ensemble.line_weight_direction, fraction,
            committed_weights))
        return false;

    for (size_t i = 0; i < ensemble.active_atoms.size(); ++i)
        ensemble.active_atoms[i].weight = committed_weights[i];
    if (!QC_SCF_Prune_Ensemble_Zero_Weight_Atoms(workspace)) return false;
    // P is already the exact float density whose raw E/F and derivative were
    // sampled at this fraction.  Commit the matching weights but retain this
    // same matrix verbatim; the next line stores P itself as its actual origin.
    // Reconstruct only on an explicit transactional rollback.
    ensemble.line_weight_direction.clear();
    return true;
}

static bool QC_SCF_Prepare_Ensemble_Corrective_Line(
    QC_SCF_WORKSPACE& workspace, int matrix_element_count, double energy,
    bool include_aufbau_vertex)
{
    QC_SCF_Ensemble_Workspace& ensemble = workspace.ensemble;
    if (ensemble.active_atoms.empty()) return false;

    // Make the current global Aufbau minimizer part of the active hull before
    // projecting the energy gradient.  A newly discovered vertex starts on
    // the simplex boundary with zero weight; the tangent-cone projection may
    // then increase it while jointly correcting every existing weight.
    if (include_aufbau_vertex)
    {
        // Only an exactly reproduced float density is the same simplex
        // vertex.  A small density RMS is not an energy-equivalence test and
        // can still carry a global FW gap larger than the strict tolerance.
        int aufbau_atom = QC_SCF_Find_Ensemble_Atom(
            workspace, matrix_element_count, workspace.alpha.d_P_new,
            workspace.runtime.unrestricted ? workspace.beta.d_P_new
                                           : (const float*)NULL);
        if (aufbau_atom < 0)
        {
            aufbau_atom = QC_SCF_Add_Ensemble_Atom(
                workspace, matrix_element_count, workspace.alpha.d_P_new,
                workspace.runtime.unrestricted ? workspace.beta.d_P_new
                                               : (const float*)NULL,
                0.0);
            if (aufbau_atom < 0) return false;
            QC_SCF_Initialize_Corrective_Inverse_Hessian(
                ensemble.active_atoms.size(),
                ensemble.corrective_inverse_hessian);
            ensemble.corrective_reference_weights.clear();
            ensemble.corrective_reference_linear_values.clear();
        }
    }

    std::vector<double> linear_values;
    std::vector<double> weights;
    linear_values.reserve(ensemble.active_atoms.size());
    weights.reserve(ensemble.active_atoms.size());
    for (const QC_SCF_Ensemble_Atom& atom : ensemble.active_atoms)
    {
        linear_values.push_back(QC_SCF_Ensemble_Fock_Density_Inner_Product(
            workspace, matrix_element_count, atom.d_alpha, atom.d_beta));
        weights.push_back(atom.weight);
    }
    if (!QC_SCF_Ensemble_Active_Set_FW_Gap(
            linear_values, weights, ensemble.active_fw_gap))
        return false;

    bool direction_ready = false;
    bool used_quasi_newton = false;
    bool used_pairwise = false;
    const bool force_pairwise =
        !include_aufbau_vertex && ensemble.corrective_force_pairwise;
    if (force_pairwise)
    {
        direction_ready = QC_SCF_Build_Pairwise_Corrective_Direction(
            linear_values, weights, ensemble.line_weight_direction);
        used_pairwise = direction_ready;
    }
    ensemble.corrective_force_pairwise = false;
    if (!direction_ready && !force_pairwise && !include_aufbau_vertex &&
        ensemble.corrective_reference_weights.size() == weights.size() &&
        ensemble.corrective_reference_linear_values.size() == weights.size())
    {
        std::vector<double> updated_inverse_hessian =
            ensemble.corrective_inverse_hessian;
        direction_ready = QC_SCF_Build_Quasi_Newton_Corrective_Direction(
            linear_values, weights,
            ensemble.corrective_reference_linear_values,
            ensemble.corrective_reference_weights, updated_inverse_hessian,
            ensemble.line_weight_direction);
        if (direction_ready)
        {
            used_quasi_newton = true;
            ensemble.corrective_inverse_hessian.swap(
                updated_inverse_hessian);
        }
    }
    if (!direction_ready && !force_pairwise)
        direction_ready = QC_SCF_Build_Fully_Corrective_Direction(
            linear_values, weights, ensemble.line_weight_direction);
    ensemble.corrective_line_used_quasi_newton = used_quasi_newton;
    ensemble.corrective_line_used_pairwise = used_pairwise;
    if (!direction_ready ||
        !QC_SCF_Store_Ensemble_Corrective_Line(workspace,
                                               matrix_element_count))
        return false;

    double host_derivative = 0.0;
    for (size_t i = 0; i < linear_values.size(); ++i)
        host_derivative +=
            linear_values[i] * ensemble.line_weight_direction[i];
    const double derivative =
        QC_SCF_Ensemble_Directional_Derivative(workspace,
                                                matrix_element_count);
    const double direction_rms =
        QC_SCF_Ensemble_Direction_RMS(workspace, matrix_element_count);
    if (!QC_SCF_Ensemble_Double_Is_Finite(host_derivative) ||
        !QC_SCF_Ensemble_Double_Is_Finite(derivative) ||
        !QC_SCF_Ensemble_Double_Is_Finite(direction_rms) ||
        !(host_derivative < -1.0e-12) ||
        !(derivative < -1.0e-12) ||
        !(direction_rms > 0.0))
        return false;

    ensemble.line_origin_derivative = derivative;
    ensemble.line_origin_energy = energy;
    ensemble.direction_density_rms = direction_rms;
    ensemble.line_derivative_tolerance =
        std::min(workspace.runtime.energy_tol,
                 std::max(1.0e-10, -0.1 * derivative));
    ensemble.line_energy_guard_multiplier = 2.0;
    ensemble.maximum_fraction = 1.0;
    ensemble.current_fraction = 1.0;
    ensemble.probe_evaluations = 0;
    ensemble.corrective_reference_weights = weights;
    ensemble.corrective_reference_linear_values = linear_values;
    if (!QC_SCF_Set_Ensemble_Trial(workspace, matrix_element_count,
                                   ensemble.current_fraction))
        return false;
    ensemble.phase = QC_SCF_ENSEMBLE_PROBE_UPPER;
    return true;
}

double QUANTUM_CHEMISTRY::Ensemble_Commutator_RMS()
{
    const int nao2 = static_cast<int>(mol.nao2);
    deviceMemset(scf_ws.ensemble.d_accum, 0, sizeof(double));
    QC_Build_DIIS_Error_Double(
        blas_handle, mol.nao, scf_ws.alpha.d_F_double, scf_ws.alpha.d_P,
        scf_ws.core.d_S, scf_ws.ensemble.d_commutator,
        scf_ws.ortho.d_dwork_nao2_2, scf_ws.ortho.d_dwork_nao2_3,
        scf_ws.ortho.d_dwork_nao2_4);
    QC_Double_Dot(nao2, scf_ws.ensemble.d_commutator,
                  scf_ws.ensemble.d_commutator, scf_ws.ensemble.d_accum);
    if (scf_ws.runtime.unrestricted)
    {
        QC_Build_DIIS_Error_Double(
            blas_handle, mol.nao, scf_ws.beta.d_F_double,
            scf_ws.beta.d_P, scf_ws.core.d_S,
            scf_ws.ensemble.d_commutator,
            scf_ws.ortho.d_dwork_nao2_2, scf_ws.ortho.d_dwork_nao2_3,
            scf_ws.ortho.d_dwork_nao2_4);
        QC_Double_Dot(nao2, scf_ws.ensemble.d_commutator,
                      scf_ws.ensemble.d_commutator,
                      scf_ws.ensemble.d_accum);
    }
    double squared_norm = 0.0;
    deviceMemcpy(&squared_norm, scf_ws.ensemble.d_accum, sizeof(double),
                 deviceMemcpyDeviceToHost);
    return std::sqrt(std::max(
        0.0, squared_norm /
                 (static_cast<double>(nao2) *
                  (scf_ws.runtime.unrestricted ? 2.0 : 1.0))));
}

bool QUANTUM_CHEMISTRY::Start_Ensemble_Probe(double energy,
                                              double density_residual)
{
    if (scf_ws.ensemble.phase != QC_SCF_ENSEMBLE_INACTIVE ||
        scf_ws.runtime.level_shift != 0.0 ||
        !(density_residual > scf_ws.runtime.density_tol))
        return false;

    const int nao2 = static_cast<int>(mol.nao2);
    QC_SCF_Store_Ensemble_Line(scf_ws, nao2);
    const double derivative =
        QC_SCF_Ensemble_Directional_Derivative(scf_ws, nao2);
    if (!QC_SCF_Ensemble_Double_Is_Finite(derivative) ||
        !(derivative < -scf_ws.runtime.energy_tol))
        return false;

    QC_SCF_Ensemble_Workspace& ensemble = scf_ws.ensemble;
    QC_SCF_Clear_Ensemble_Active_Set(scf_ws);
    if (QC_SCF_Add_Ensemble_Atom(
            scf_ws, nao2, scf_ws.alpha.d_P,
            scf_ws.runtime.unrestricted ? scf_ws.beta.d_P
                                        : (const float*)NULL,
            1.0) < 0)
        return false;
    ensemble.committed_energy = energy;
    if (QC_SCF_Prepare_Ensemble_Corrective_Line(scf_ws, nao2, energy, true))
        return true;
    QC_SCF_Clear_Ensemble_Active_Set(scf_ws);
    return false;
}

int QUANTUM_CHEMISTRY::Advance_Ensemble_Search(int iter, int md_step,
                                                double energy,
                                                double delta_energy)
{
    (void)delta_energy;
    const int nao2 = static_cast<int>(mol.nao2);
    QC_SCF_Ensemble_Workspace& ensemble = scf_ws.ensemble;

    auto print_iteration = [&](const char* phase, double fraction,
                               double line_metric, double commutator,
                               bool reports_global_fw_gap = false)
    {
        if (!scf_ws.runtime.print_iter || CONTROLLER::MPI_rank != 0) return;
        FILE* out = (scf_output_file != NULL) ? scf_output_file : stdout;
        fprintf(out,
                "Step %6d | SCF Ensemble %3d | E(Ha)=%.12f | phase=%s "
                "| x=%.12e | %s(Ha)=%.12e | comm=%.12e | active-gap=%.12e "
                "| active=%zu\n",
                md_step, iter + 1, energy, phase, fraction,
                reports_global_fw_gap ? "global-gap" : "g", line_metric,
                commutator, ensemble.active_fw_gap,
                ensemble.active_atoms.size());
        fflush(out);
    };

    auto abort_search = [&]()
    {
        // Never expose P_new from a failed trial to ordinary SCF.  Restore the
        // last committed active-set density while its atoms are still alive,
        // then keep one explicit ensemble recovery phase so the next loop
        // rebuilds raw E/F from that P and skips Check_Convergence entirely.
        QC_SCF_Reconstruct_Ensemble_Density(scf_ws, nao2);
        ensemble.interior_minimum_confirmed = false;
        QC_SCF_Clear_Ensemble_Active_Set(scf_ws);
        ensemble.current_fraction = 0.0;
        ensemble.phase = QC_SCF_ENSEMBLE_RECOVER_COMMITTED;
        return 1;
    };

    auto retry_corrective_line = [&]()
    {
        // A quasi-Newton direction is only a proposal.  If its raw line map
        // cannot establish an energy-consistent minimum, restore the exact
        // committed density and restart the same active face with the
        // tangent-cone projected direction.  A projected-line failure remains
        // fatal; this prevents an endless alternation between two policies.
        const bool retry_projected =
            ensemble.corrective_line_used_quasi_newton;
        const bool retry_pairwise =
            !retry_projected &&
            !ensemble.corrective_line_used_pairwise;
        if (!retry_projected && !retry_pairwise) return abort_search();
        QC_SCF_Reconstruct_Ensemble_Density(scf_ws, nao2);
        ensemble.line_weight_direction.clear();
        QC_SCF_Initialize_Corrective_Inverse_Hessian(
            ensemble.active_atoms.size(),
            ensemble.corrective_inverse_hessian);
        ensemble.corrective_reference_weights.clear();
        ensemble.corrective_reference_linear_values.clear();
        ensemble.corrective_line_used_quasi_newton = false;
        ensemble.corrective_line_used_pairwise = false;
        ensemble.corrective_force_pairwise = retry_pairwise;
        ensemble.current_fraction = 0.0;
        ensemble.phase = QC_SCF_ENSEMBLE_VERIFY_COMMITTED;
        return 1;
    };

    auto prepare_spectral_orbital_correction =
        [&](double origin_commutator) -> bool
    {
        QC_SCF_Ensemble_Workspace& correction = scf_ws.ensemble;
        const int nao = mol.nao;
        const int ne =
            scf_ws.ortho.nao_eff > 0 ? scf_ws.ortho.nao_eff : nao;
        if (correction.active_atoms.empty() || ne <= 0 || ne > nao ||
            !QC_SCF_Ensemble_Double_Is_Finite(origin_commutator) ||
            !(origin_commutator > scf_ws.runtime.density_tol) ||
            correction.d_spectral_origin_fock_alpha == NULL ||
            (scf_ws.runtime.unrestricted &&
             correction.d_spectral_origin_fock_beta == NULL))
            return false;

        // Raw trial evaluations overwrite each channel's physical Fock
        // matrix.  Preserve the symmetric operator that defined this
        // spectral line so every later representational rounding decision is
        // measured against the same fixed-F objective.
        deviceMemcpy(correction.d_spectral_origin_fock_alpha,
                     scf_ws.alpha.d_F_double, sizeof(double) * nao2,
                     deviceMemcpyDeviceToDevice);
        if (scf_ws.runtime.unrestricted)
            deviceMemcpy(correction.d_spectral_origin_fock_beta,
                         scf_ws.beta.d_F_double, sizeof(double) * nao2,
                         deviceMemcpyDeviceToDevice);

        deviceMemcpy(correction.d_origin_alpha, scf_ws.alpha.d_P,
                     sizeof(float) * nao2, deviceMemcpyDeviceToDevice);
        if (scf_ws.runtime.unrestricted)
            deviceMemcpy(correction.d_origin_beta, scf_ws.beta.d_P,
                         sizeof(float) * nao2,
                         deviceMemcpyDeviceToDevice);

        const std::vector<double> original_particle_counts =
            scf_ws.runtime.unrestricted
                ? std::vector<double>{
                      QC_SCF_Ensemble_Particle_Count(
                          scf_ws, nao2, correction.d_origin_alpha),
                      QC_SCF_Ensemble_Particle_Count(
                          scf_ws, nao2, correction.d_origin_beta)}
                : std::vector<double>{QC_SCF_Ensemble_Particle_Count(
                      scf_ws, nao2, correction.d_origin_alpha)};
        const std::vector<double> nominal_particle_counts =
            scf_ws.runtime.unrestricted
                ? std::vector<double>{
                      static_cast<double>(scf_ws.runtime.n_alpha),
                      static_cast<double>(scf_ws.runtime.n_beta)}
                : std::vector<double>{
                      2.0 * static_cast<double>(scf_ws.runtime.n_alpha)};
        for (size_t channel = 0; channel < original_particle_counts.size();
             ++channel)
        {
            const double represented = original_particle_counts[channel];
            const double nominal = nominal_particle_counts[channel];
            if (!QC_SCF_Ensemble_Double_Is_Finite(represented) ||
                std::fabs(represented - nominal) >
                    scf_ws.runtime.density_tol * std::max(1.0, nominal))
                return false;
        }

        auto restore_origin = [&]()
        {
            deviceMemcpy(scf_ws.alpha.d_P, correction.d_origin_alpha,
                         sizeof(float) * nao2,
                         deviceMemcpyDeviceToDevice);
            if (scf_ws.runtime.unrestricted)
            {
                deviceMemcpy(scf_ws.beta.d_P, correction.d_origin_beta,
                             sizeof(float) * nao2,
                             deviceMemcpyDeviceToDevice);
                QC_Add_Matrix(nao2, scf_ws.alpha.d_P, scf_ws.beta.d_P,
                              scf_ws.direct.d_Ptot);
            }
        };

        auto build_spectral_orbital_channel =
            [&](QC_SCF_Spin_Channel& channel, double particle_count,
                double maximum_occupation,
                QC_SCF_Eigensolver_Channel spin_channel,
                double& spectrum_preserving_linear_change) -> bool
        {
            spectrum_preserving_linear_change = 0.0;
            double* dC = scf_ws.ortho.d_dwork_nao2_1;
            double* dS_or_N = scf_ws.ortho.d_dwork_nao2_2;
            double* dFp_or_P = scf_ws.ortho.d_dwork_nao2_3;
            double* dTmp = scf_ws.ortho.d_dwork_nao2_4;

            auto fixed_fock_value = [&](const double* density) -> double
            {
                deviceMemset(scf_ws.ensemble.d_accum, 0, sizeof(double));
                QC_Double_Dot(nao2, channel.d_F_double, density,
                              scf_ws.ensemble.d_accum);
                double value = 0.0;
                deviceMemcpy(&value, scf_ws.ensemble.d_accum,
                             sizeof(double), deviceMemcpyDeviceToHost);
                return value;
            };

            // Reproduce the current raw generalized Fock eigenbasis without
            // modifying the cached physical F.  All matrices remain double
            // through N=C^T S P S C.  Diagonalizing N gives the invariant
            // natural-occupation spectrum; assigning it in descending order
            // to the ascending Fock orbitals is the fixed-F orbital-rotation
            // minimizer and removes all F/P coherences.
            deviceMemcpy(dC, channel.d_F_double, sizeof(double) * nao2,
                         deviceMemcpyDeviceToDevice);
            QC_Dgemm_NN(blas_handle, nao, ne, nao, dC, nao,
                        scf_ws.ortho.d_X, nao, dS_or_N, ne);
            QC_Dgemm_TN(blas_handle, ne, ne, nao, scf_ws.ortho.d_X,
                        nao, dS_or_N, ne, dFp_or_P, ne);
            QC_Symmetrize_Double_Matrix(ne, dFp_or_P);
            const int threads = 256;
            int info = 0;
            const int api_status = QC_Diagonalize_Double(
                solver_handle, ne, dFp_or_P, scf_ws.ortho.d_dW_double,
                scf_ws.ortho.d_solver_work_double,
                scf_ws.ortho.lwork_double, &info);
            const bool solver_ok = QC_SCF_Require_Eigensolver_Success(
                QC_SCF_EIGENSOLVER_FOCK, spin_channel, ne, api_status, info,
                [&](const QC_SCF_Eigensolver_Failure& failure)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorSimulationBreakDown,
                        "QUANTUM_CHEMISTRY::Advance_Ensemble_Search",
                        "Reason:\n    eigensolver failed while building the "
                        "spectral orbital correction for channel %s: "
                        "dimension=%d, "
                        "api_status=%d, info=%d\n",
                        failure.channel_name, failure.dimension,
                        failure.api_status, failure.info);
                });
            if (!solver_ok) return false;

            std::vector<double> fock_eigenvalues(ne);
            deviceMemcpy(fock_eigenvalues.data(),
                         scf_ws.ortho.d_dW_double,
                         sizeof(double) * ne, deviceMemcpyDeviceToHost);

            QC_Dgemm_NT(blas_handle, nao, ne, ne, scf_ws.ortho.d_X,
                        nao, dFp_or_P, ne, dC, ne);
            QC_Float_To_Double(nao2, scf_ws.core.d_S, dS_or_N);
            QC_Float_To_Double(nao2, channel.d_P, dFp_or_P);
            QC_Dgemm_NN(blas_handle, nao, nao, nao, dS_or_N, nao,
                        dFp_or_P, nao, dTmp, nao);
            QC_Dgemm_NN(blas_handle, nao, nao, nao, dTmp, nao, dS_or_N,
                        nao, dFp_or_P, nao);
            QC_Dgemm_NN(blas_handle, nao, ne, nao, dFp_or_P, nao, dC,
                        ne, dTmp, ne);
            QC_Dgemm_TN(blas_handle, ne, ne, nao, dC, ne, dTmp, ne,
                        dS_or_N, ne);

            QC_Symmetrize_Double_Matrix(ne, dS_or_N);
            Launch_Device_Kernel(
                QC_SCF_Extract_Diagonal_Kernel,
                Positive_Int_Ceil_Div(ne, threads), threads, 0, 0, ne,
                dS_or_N, dTmp);
            std::vector<double> fock_basis_density_diagonal(ne);
            deviceMemcpy(fock_basis_density_diagonal.data(), dTmp,
                         sizeof(double) * ne, deviceMemcpyDeviceToHost);
            info = 0;
            const int occupation_api_status = QC_Diagonalize_Double(
                solver_handle, ne, dS_or_N,
                scf_ws.ortho.d_dW_double,
                scf_ws.ortho.d_solver_work_double,
                scf_ws.ortho.lwork_double, &info);
            const bool occupation_solver_ok =
                QC_SCF_Require_Eigensolver_Success(
                    QC_SCF_EIGENSOLVER_ENSEMBLE_OCCUPATION, spin_channel,
                    ne, occupation_api_status, info,
                    [&](const QC_SCF_Eigensolver_Failure& failure)
                    {
                        controller->Throw_Formatted_SPONGE_Error(
                            spongeErrorSimulationBreakDown,
                            "QUANTUM_CHEMISTRY::Advance_Ensemble_Search",
                            "Reason:\n    eigensolver failed during %s for "
                            "channel %s: dimension=%d, api_status=%d, "
                            "info=%d\n",
                            failure.stage_name, failure.channel_name,
                            failure.dimension, failure.api_status,
                            failure.info);
                    });
            if (!occupation_solver_ok) return false;

            std::vector<double> natural_occupation_eigenvalues(ne);
            deviceMemcpy(natural_occupation_eigenvalues.data(),
                         scf_ws.ortho.d_dW_double,
                         sizeof(double) * ne, deviceMemcpyDeviceToHost);
            if (!QC_SCF_Spectrum_Preserving_Orbital_Rotation_Change(
                    fock_eigenvalues, fock_basis_density_diagonal,
                    natural_occupation_eigenvalues,
                    spectrum_preserving_linear_change))
                return false;
            std::vector<double> occupations;
            if (!QC_SCF_Prepare_Spectral_Orbital_Occupations(
                    natural_occupation_eigenvalues, particle_count,
                    maximum_occupation, scf_ws.runtime.density_tol,
                    occupations))
                return false;
            deviceMemcpy(scf_ws.ortho.d_dW_double, occupations.data(),
                         sizeof(double) * ne, deviceMemcpyHostToDevice);
            Launch_Device_Kernel(
                QC_SCF_Scale_Orbitals_By_Occupation_Kernel,
                Positive_Int_Ceil_Div(nao * ne, threads), threads, 0, 0,
                nao * ne, ne, dC, scf_ws.ortho.d_dW_double, dFp_or_P);
            QC_Dgemm_NT(blas_handle, nao, nao, ne, dFp_or_P, ne, dC,
                        ne, dTmp, nao);
            QC_Symmetrize_Double_Matrix(nao, dTmp);
            const double projected_fixed_fock = fixed_fock_value(dTmp);
            // Nearest rounding is unbiased and normally preserves the
            // projected endpoint.  Only use one-sided, objective-aware
            // rounding when the float representation raises F:P relative to
            // that same double endpoint.  Comparing with the off-manifold
            // origin would misclassify the separate nominal-particle repair
            // as rounding error and accumulate a trace bias.
            QC_Double_To_Float(nao2, dTmp, channel.d_P);
            QC_Float_To_Double(nao2, channel.d_P, dS_or_N);
            double float_fixed_fock = fixed_fock_value(dS_or_N);
            const double channel_rounding_budget =
                scf_ws.runtime.energy_tol /
                (scf_ws.runtime.unrestricted ? 2.0 : 1.0);
            double rounding_error =
                float_fixed_fock - projected_fixed_fock;
            if (rounding_error >
                channel_rounding_budget)
            {
                QC_Round_Symmetric_Double_Matrix_For_Nonincreasing_Linear_Objective(
                    nao, dTmp, channel.d_F_double, channel.d_P);
                QC_Float_To_Double(nao2, channel.d_P, dS_or_N);
                float_fixed_fock = fixed_fock_value(dS_or_N);
                rounding_error =
                    float_fixed_fock - projected_fixed_fock;
            }
            if (!QC_SCF_Ensemble_Double_Is_Finite(projected_fixed_fock) ||
                !QC_SCF_Ensemble_Double_Is_Finite(float_fixed_fock) ||
                !QC_SCF_Ensemble_Double_Is_Finite(rounding_error) ||
                rounding_error > channel_rounding_budget)
                return false;
            return true;
        };

        double spectrum_preserving_linear_change = 0.0;
        double alpha_spectrum_preserving_linear_change = 0.0;
        bool built = build_spectral_orbital_channel(
            scf_ws.alpha, nominal_particle_counts[0],
            scf_ws.runtime.unrestricted ? 1.0 : 2.0,
            QC_SCF_EIGENSOLVER_CHANNEL_ALPHA,
            alpha_spectrum_preserving_linear_change);
        spectrum_preserving_linear_change +=
            alpha_spectrum_preserving_linear_change;
        if (built && scf_ws.runtime.unrestricted)
        {
            double beta_spectrum_preserving_linear_change = 0.0;
            built = build_spectral_orbital_channel(
                scf_ws.beta, nominal_particle_counts[1],
                1.0, QC_SCF_EIGENSOLVER_CHANNEL_BETA,
                beta_spectrum_preserving_linear_change);
            spectrum_preserving_linear_change +=
                beta_spectrum_preserving_linear_change;
        }
        if (!built)
        {
            restore_origin();
            return false;
        }

        const int threads = 256;
        Launch_Device_Kernel(
            QC_SCF_Build_Spectral_Orbital_Direction_Kernel,
            Positive_Int_Ceil_Div(nao2, threads), threads, 0, 0, nao2,
            correction.d_origin_alpha, scf_ws.alpha.d_P,
            correction.d_direction_alpha);
        if (scf_ws.runtime.unrestricted)
        {
            Launch_Device_Kernel(
                QC_SCF_Build_Spectral_Orbital_Direction_Kernel,
                Positive_Int_Ceil_Div(nao2, threads), threads, 0, 0, nao2,
                correction.d_origin_beta, scf_ws.beta.d_P,
                correction.d_direction_beta);
            QC_Add_Matrix(nao2, scf_ws.alpha.d_P, scf_ws.beta.d_P,
                          scf_ws.direct.d_Ptot);
        }

        double total_linear_change =
            QC_SCF_Spectral_Origin_Fock_Direction(
                scf_ws, nao2,
                correction.d_spectral_origin_fock_alpha,
                correction.d_direction_alpha);
        if (scf_ws.runtime.unrestricted)
            total_linear_change +=
                QC_SCF_Spectral_Origin_Fock_Direction(
                    scf_ws, nao2,
                    correction.d_spectral_origin_fock_beta,
                    correction.d_direction_beta);
        const double repair_linear_component =
            total_linear_change - spectrum_preserving_linear_change;

        const std::vector<double> corrected_particle_counts =
            scf_ws.runtime.unrestricted
                ? std::vector<double>{QC_SCF_Ensemble_Particle_Count(
                                          scf_ws, nao2, scf_ws.alpha.d_P),
                                      QC_SCF_Ensemble_Particle_Count(
                                          scf_ws, nao2, scf_ws.beta.d_P)}
                : std::vector<double>{QC_SCF_Ensemble_Particle_Count(
                      scf_ws, nao2, scf_ws.alpha.d_P)};
        const double corrected_commutator = Ensemble_Commutator_RMS();
        const double direction_rms =
            QC_SCF_Ensemble_Direction_RMS(scf_ws, nao2);
        const bool properties_valid =
            QC_SCF_Spectral_Orbital_Properties_Are_Valid(
                original_particle_counts, corrected_particle_counts,
                spectrum_preserving_linear_change, corrected_commutator,
                scf_ws.runtime.energy_tol, scf_ws.runtime.density_tol);
        const bool direction_valid =
            QC_SCF_Ensemble_Double_Is_Finite(direction_rms) &&
            direction_rms > 0.0 &&
            QC_SCF_Ensemble_Double_Is_Finite(total_linear_change) &&
            QC_SCF_Ensemble_Double_Is_Finite(repair_linear_component);
        if (!properties_valid || !direction_valid)
        {
            restore_origin();
            return false;
        }

        correction.line_origin_energy = energy;
        correction.direction_density_rms = direction_rms;
        correction.line_energy_guard_multiplier = 2.0;
        correction.current_fraction = 1.0;
        correction.maximum_fraction = 1.0;
        correction.probe_evaluations = 0;
        correction.spectral_orbital_origin_commutator = origin_commutator;
        correction.spectral_orbital_repair_linear_component =
            repair_linear_component;
        correction.phase = QC_SCF_ENSEMBLE_PROBE_SPECTRAL_ORBITAL;
        return true;
    };

    auto verify_global_kkt_and_prepare_next =
        [&](bool repeated_verification) -> int
    {
        if (!repeated_verification)
        {
            // This is the first raw F[P]/E[P] evaluation after rebuilding P
            // from committed active-set weights.  Advance the committed
            // energy here, never from the preceding line-search sample.
            if (!QC_SCF_Ensemble_Energy_Within_Line_Guard(
                    ensemble.line_origin_energy, energy,
                    scf_ws.runtime.energy_tol,
                    ensemble.line_energy_guard_multiplier))
                return abort_search();
            ensemble.committed_energy = energy;
        }

        // Re-form the standard Frank-Wolfe direction P_Aufbau-P.  Corrective
        // active-hull directions accelerate optimization, but only this global
        // linear-oracle gap is a KKT certificate.
        QC_SCF_Store_Ensemble_Line(scf_ws, nao2);
        const double global_direction_rms =
            QC_SCF_Ensemble_Direction_RMS(scf_ws, nao2);
        const double global_directional_derivative =
            QC_SCF_Ensemble_Directional_Derivative(scf_ws, nao2);
        double global_fw_gap = 0.0;
        const bool global_gap_valid = QC_SCF_Ensemble_Global_FW_Gap(
            global_directional_derivative, global_fw_gap);
        const double commutator = Ensemble_Commutator_RMS();
        double active_fw_gap = 0.0;
        const bool active_gap_valid =
            QC_SCF_Evaluate_Ensemble_Active_Set_FW_Gap(
                scf_ws, nao2, active_fw_gap);
        ensemble.active_fw_gap = active_fw_gap;
        print_iteration(repeated_verification ? "verify-kkt-repeat"
                                              : "verify-kkt",
                        ensemble.current_fraction, global_fw_gap,
                        commutator, true);

        if (!global_gap_valid ||
            !QC_SCF_Ensemble_Double_Is_Finite(global_direction_rms) ||
            !QC_SCF_Ensemble_Double_Is_Finite(commutator) ||
            !active_gap_valid)
            return abort_search();

        const bool kkt_stationary =
            global_fw_gap <= scf_ws.runtime.energy_tol &&
            active_fw_gap <= scf_ws.runtime.energy_tol &&
            commutator <= scf_ws.runtime.density_tol;
        if (kkt_stationary && !repeated_verification)
        {
            // Keep exactly the same P and rebuild raw F[P] once more.  This
            // makes final E, F, P and the KKT certificate independently
            // reproducible rather than accepting a one-off diagonalization.
            ensemble.verification_global_fw_gap = global_fw_gap;
            ensemble.verification_commutator = commutator;
            ensemble.verification_active_fw_gap = active_fw_gap;
            ensemble.phase = QC_SCF_ENSEMBLE_VERIFY_KKT_REPEAT;
            return 1;
        }
        if (kkt_stationary)
        {
            // The second raw build is a reproducibility check, not another
            // optimization step.  P is intentionally unchanged; require E,
            // the global FW gap, and the generalized commutator to reproduce.
            const bool energy_repeated =
                std::fabs(energy - ensemble.committed_energy) <=
                scf_ws.runtime.energy_tol;
            const bool global_gap_repeated =
                std::fabs(global_fw_gap -
                          ensemble.verification_global_fw_gap) <=
                scf_ws.runtime.energy_tol;
            const bool commutator_repeated =
                std::fabs(commutator - ensemble.verification_commutator) <=
                scf_ws.runtime.density_tol;
            const bool active_gap_repeated =
                std::fabs(active_fw_gap -
                          ensemble.verification_active_fw_gap) <=
                scf_ws.runtime.energy_tol;
            if (!energy_repeated || !global_gap_repeated ||
                !commutator_repeated || !active_gap_repeated)
                return abort_search();
            ensemble.confirmed = true;
            ensemble.phase = QC_SCF_ENSEMBLE_INACTIVE;
            const int converged = 1;
            deviceMemcpy(scf_ws.runtime.d_converged, &converged, sizeof(int),
                         deviceMemcpyHostToDevice);
            return 2;
        }

        // A requested repeat must reproduce the stationary certificate.  If
        // it does not, the previous build was not an independently verified
        // solution and cannot be turned into another optimization line.
        if (repeated_verification) return abort_search();

        // The occupation and orbital KKT blocks both need progress.  Select
        // the block with the larger residual normalized by its unchanged
        // final tolerance (Gauss-Southwell scheduling).  Thus a spectral step
        // cannot starve a larger linear gap, while an unresolved commutator
        // cannot be postponed until the active hull is solved to completion.
        // This changes only scheduling; the final joint KKT test above remains
        // unchanged.
        if (QC_SCF_Should_Start_Spectral_Orbital_Correction(
                global_fw_gap, active_fw_gap, commutator,
                scf_ws.runtime.energy_tol, scf_ws.runtime.density_tol,
                repeated_verification))
        {
            if (!prepare_spectral_orbital_correction(commutator))
                return abort_search();
            return 1;
        }

        // Correct the current finite hull only to a forcing tolerance
        // proportional to the current global gap.  Solving every intermediate
        // hull to the final tolerance wastes raw Fock builds because the next
        // global vertex changes that hull.  This target only schedules the next
        // oracle call; the final certificate above still tests the global gap,
        // active gap, and commutator separately at their original tolerances.
        double active_correction_target = 0.0;
        if (!QC_SCF_Ensemble_Active_Correction_Target(
                global_fw_gap,
                scf_ws.runtime.energy_tol, active_correction_target))
            return abort_search();
        if (active_fw_gap > active_correction_target)
        {
            if (!QC_SCF_Prepare_Ensemble_Corrective_Line(
                    scf_ws, nao2, energy, false))
                return abort_search();
            return 1;
        }

        if (!(global_directional_derivative < -1.0e-12) ||
            !(global_direction_rms > scf_ws.runtime.density_tol))
            return abort_search();
        if (!QC_SCF_Prepare_Ensemble_Corrective_Line(
                scf_ws, nao2, energy, true))
            return abort_search();
        return 1;
    };

    if (!QC_SCF_Ensemble_Double_Is_Finite(energy)) return abort_search();

    if (ensemble.phase == QC_SCF_ENSEMBLE_RECOVER_COMMITTED)
    {
        // This iteration's raw E/F now belongs to the restored committed P.
        // Do not publish its diagonalized P_new in the same iteration; hand
        // the unchanged P back to ordinary SCF on the following map.
        print_iteration("recover-committed", 0.0, 0.0, 0.0);
        ensemble.phase = QC_SCF_ENSEMBLE_INACTIVE;
        scf_ws.diis.diis_hist_count = scf_ws.diis.diis_hist_head = 0;
        scf_ws.diis.last_enorm = 1.0e10;
        scf_ws.runtime.convergence.level_shift_stage =
            QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT;
        scf_ws.runtime.convergence.verifying_physical_fixed_point = true;
        scf_ws.runtime.level_shift = 0.0;
        return 1;
    }

    if (ensemble.phase == QC_SCF_ENSEMBLE_PROBE_SPECTRAL_ORBITAL)
    {
        ++ensemble.probe_evaluations;
        const double correction_derivative =
            QC_SCF_Ensemble_Directional_Derivative(scf_ws, nao2);
        const double corrected_commutator = Ensemble_Commutator_RMS();
        print_iteration("spectral-orbital-trial", ensemble.current_fraction,
                        correction_derivative, corrected_commutator);
        if (QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
                ensemble.line_origin_energy, energy,
                ensemble.current_fraction,
                ensemble.spectral_orbital_repair_linear_component,
                ensemble.spectral_orbital_origin_commutator,
                corrected_commutator, scf_ws.runtime.energy_tol,
                scf_ws.runtime.density_tol,
                ensemble.line_energy_guard_multiplier))
        {
            // Only the raw observation commits the trial.  The previous hull
            // remains the rollback source until this point; now restart from
            // the accepted legal ensemble density as one weight-one atom.
            QC_SCF_Clear_Ensemble_Active_Set(scf_ws);
            if (QC_SCF_Add_Ensemble_Atom(
                    scf_ws, nao2, scf_ws.alpha.d_P,
                    scf_ws.runtime.unrestricted ? scf_ws.beta.d_P
                                                : (const float*)NULL,
                    1.0) < 0)
                return abort_search();
            ensemble.committed_energy = energy;
            ensemble.line_origin_energy = energy;
            ensemble.current_fraction = 0.0;
            return verify_global_kkt_and_prepare_next(false);
        }

        double next_fraction = 0.0;
        if (!QC_SCF_Next_Spectral_Orbital_Fraction(
                ensemble.current_fraction, ensemble.probe_evaluations,
                next_fraction))
            return abort_search();
        ensemble.current_fraction = next_fraction;
        if (!QC_SCF_Set_Spectral_Orbital_Trial(
                scf_ws, mol.nao, nao2, next_fraction))
            return abort_search();
        return 1;
    }

    if (ensemble.phase == QC_SCF_ENSEMBLE_VERIFY_KKT_REPEAT)
        return verify_global_kkt_and_prepare_next(true);
    if (ensemble.phase == QC_SCF_ENSEMBLE_VERIFY_COMMITTED)
        return verify_global_kkt_and_prepare_next(false);

    const double derivative =
        QC_SCF_Ensemble_Directional_Derivative(scf_ws, nao2);
    if (!QC_SCF_Ensemble_Double_Is_Finite(derivative)) return abort_search();

    if (ensemble.phase == QC_SCF_ENSEMBLE_PROBE_UPPER)
    {
        ++ensemble.probe_evaluations;
        print_iteration("probe-upper", ensemble.current_fraction, derivative,
                        0.0);
        if (derivative > ensemble.line_derivative_tolerance)
        {
            const QC_SCF_Ensemble_Root_Decision decision =
                QC_SCF_Start_Ensemble_Root(
                    ensemble.line_origin_derivative,
                    ensemble.line_origin_energy,
                    derivative, energy, ensemble.direction_density_rms,
                    ensemble.line_derivative_tolerance, ensemble.bracket, false,
                    ensemble.current_fraction);
            if (decision.status != QC_SCF_ENSEMBLE_ROOT_EVALUATE)
                return retry_corrective_line();
            ensemble.current_fraction = decision.next_fraction;
            if (!QC_SCF_Set_Ensemble_Trial(scf_ws, nao2,
                                           ensemble.current_fraction))
                return abort_search();
            ensemble.phase = QC_SCF_ENSEMBLE_EVALUATE_ROOT;
            return 1;
        }

        // A non-positive derivative at gamma_max is a valid drop/swap step.
        // If the independently accumulated energy disagrees, backtrack on
        // this same fixed line; never turn a numerical mismatch into a silent
        // energy-increasing sequence of different directions.
        if (!QC_SCF_Ensemble_Energy_Within_Line_Guard(
                ensemble.line_origin_energy, energy,
                scf_ws.runtime.energy_tol))
        {
            if (ensemble.probe_evaluations >=
                    QC_SCF_ENSEMBLE_MAX_ROOT_EVALUATIONS ||
                ensemble.current_fraction *
                        ensemble.direction_density_rms <=
                    scf_ws.runtime.density_tol)
                return retry_corrective_line();
            ensemble.current_fraction *= 0.5;
            if (!QC_SCF_Set_Ensemble_Trial(scf_ws, nao2,
                                           ensemble.current_fraction))
                return abort_search();
            return 1;
        }

        if (!QC_SCF_Commit_Ensemble_Corrective_Weights(
                scf_ws, nao2, ensemble.current_fraction))
            return abort_search();
        return verify_global_kkt_and_prepare_next(false);
    }

    if (ensemble.phase == QC_SCF_ENSEMBLE_VERIFY_LINE_BEST)
    {
        print_iteration("verify-line-best", ensemble.current_fraction,
                        derivative, 0.0);
        if (energy > ensemble.bracket.best_energy +
                         0.25 * scf_ws.runtime.energy_tol ||
            !QC_SCF_Ensemble_Energy_Within_Line_Guard(
                ensemble.line_origin_energy, energy,
                scf_ws.runtime.energy_tol))
            return retry_corrective_line();
        if (!QC_SCF_Commit_Ensemble_Corrective_Weights(
                scf_ws, nao2, ensemble.current_fraction))
            return retry_corrective_line();
        return verify_global_kkt_and_prepare_next(false);
    }

    if (ensemble.phase == QC_SCF_ENSEMBLE_VERIFY_LINE_STATIONARITY)
    {
        print_iteration("verify-line-stationarity", ensemble.current_fraction,
                        derivative, 0.0);
        if (std::fabs(derivative -
                      ensemble.bracket.best_stationarity_derivative) >
                ensemble.line_derivative_tolerance ||
            energy > ensemble.bracket.best_stationarity_energy +
                         0.25 * scf_ws.runtime.energy_tol ||
            !QC_SCF_Ensemble_Energy_Within_Line_Guard(
                ensemble.line_origin_energy, energy,
                scf_ws.runtime.energy_tol, 4.0))
            return retry_corrective_line();
        if (!QC_SCF_Commit_Ensemble_Corrective_Weights(
                scf_ws, nao2, ensemble.current_fraction))
            return retry_corrective_line();
        ensemble.line_energy_guard_multiplier = 4.0;
        return verify_global_kkt_and_prepare_next(false);
    }

    if (ensemble.phase != QC_SCF_ENSEMBLE_EVALUATE_ROOT)
        return abort_search();

    print_iteration("root", ensemble.current_fraction, derivative, 0.0);
    const QC_SCF_Ensemble_Root_Decision decision =
        QC_SCF_Observe_Ensemble_Trial(
            ensemble.bracket, ensemble.current_fraction, derivative, energy,
            ensemble.line_derivative_tolerance, scf_ws.runtime.energy_tol,
            scf_ws.runtime.density_tol);
    if (decision.status == QC_SCF_ENSEMBLE_ROOT_EVALUATE)
    {
        ensemble.current_fraction = decision.next_fraction;
        if (!QC_SCF_Set_Ensemble_Trial(scf_ws, nao2,
                                       ensemble.current_fraction))
            return abort_search();
        return 1;
    }
    if (decision.status == QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY)
    {
        ensemble.current_fraction = decision.next_fraction;
        if (!QC_SCF_Set_Ensemble_Trial(scf_ws, nao2,
                                       ensemble.current_fraction))
            return abort_search();
        ensemble.phase = QC_SCF_ENSEMBLE_VERIFY_LINE_BEST;
        return 1;
    }
    if (decision.status == QC_SCF_ENSEMBLE_ROOT_USE_BEST_STATIONARITY)
    {
        ensemble.current_fraction = decision.next_fraction;
        if (!QC_SCF_Set_Ensemble_Trial(scf_ws, nao2,
                                       ensemble.current_fraction))
            return abort_search();
        ensemble.phase = QC_SCF_ENSEMBLE_VERIFY_LINE_STATIONARITY;
        return 1;
    }
    if (decision.status != QC_SCF_ENSEMBLE_ROOT_STATIONARY)
    {
        print_iteration("root-policy-failed", ensemble.current_fraction,
                        derivative, 0.0);
        return retry_corrective_line();
    }
    if (!QC_SCF_Ensemble_Energy_Within_Line_Guard(
            ensemble.line_origin_energy, energy,
            scf_ws.runtime.energy_tol))
        return retry_corrective_line();
    if (!QC_SCF_Commit_Ensemble_Corrective_Weights(
            scf_ws, nao2, ensemble.current_fraction))
    {
        print_iteration("active-commit-failed", ensemble.current_fraction,
                        derivative, 0.0);
        return retry_corrective_line();
    }
    ensemble.interior_minimum_confirmed = true;
    return verify_global_kkt_and_prepare_next(false);
}
