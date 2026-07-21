#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "quantum_chemistry/scf/convergence_policy.hpp"
#include "quantum_chemistry/scf/diis_coefficients.hpp"
#include "quantum_chemistry/scf/eigensolver_policy.hpp"
#include "quantum_chemistry/scf/ensemble_policy.hpp"

namespace
{

bool Close(double actual, double expected, double tolerance = 1.0e-12)
{
    return QC_SCF_DIIS_Double_Is_Finite(actual) &&
           std::fabs(actual - expected) <= tolerance;
}

double Double_From_Bits(std::uint64_t bit_pattern)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bit_pattern));
#endif
    double value = 0.0;
    std::memcpy(&value, &bit_pattern, sizeof(value));
    return value;
}

bool Eigensolver_Failure_Matches(QC_SCF_Eigensolver_Stage stage,
                                 QC_SCF_Eigensolver_Channel channel,
                                 int dimension, int api_status, int info,
                                 const char* stage_name,
                                 const char* channel_name)
{
    QC_SCF_Eigensolver_Failure observed = {};
    int calls = 0;
    const bool accepted = QC_SCF_Require_Eigensolver_Success(
        stage, channel, dimension, api_status, info,
        [&](const QC_SCF_Eigensolver_Failure& failure)
        {
            observed = failure;
            ++calls;
        });
    return !accepted && calls == 1 && observed.stage == stage &&
           observed.channel == channel && observed.dimension == dimension &&
           observed.api_status == api_status && observed.info == info &&
           std::strcmp(observed.stage_name, stage_name) == 0 &&
           std::strcmp(observed.channel_name, channel_name) == 0;
}

bool Workspace_Failure_Matches(int dimension, int api_status,
                               int workspace_size, bool workspace_available)
{
    QC_SCF_Eigensolver_Workspace_Failure observed = {};
    int calls = 0;
    const bool accepted = QC_SCF_Require_Eigensolver_Workspace(
        dimension, api_status, workspace_size, workspace_available,
        [&](const QC_SCF_Eigensolver_Workspace_Failure& failure)
        {
            observed = failure;
            ++calls;
        });
    return !accepted && calls == 1 &&
           observed.stage == QC_SCF_EIGENSOLVER_WORKSPACE &&
           observed.channel == QC_SCF_EIGENSOLVER_CHANNEL_SHARED &&
           observed.dimension == dimension &&
           observed.api_status == api_status &&
           observed.workspace_size == workspace_size &&
           observed.workspace_available == workspace_available &&
           std::strcmp(observed.stage_name,
                       "shared double-precision eigensolver workspace query") ==
               0 &&
           std::strcmp(observed.channel_name, "shared") == 0;
}

std::vector<double> Matrix_Multiply(const std::vector<double>& a, int rows,
                                    int inner,
                                    const std::vector<double>& b, int cols)
{
    std::vector<double> result(rows * cols, 0.0);
    for (int i = 0; i < rows; ++i)
        for (int k = 0; k < inner; ++k)
            for (int j = 0; j < cols; ++j)
                result[i * cols + j] +=
                    a[i * inner + k] * b[k * cols + j];
    return result;
}

std::vector<double> Matrix_Transpose(const std::vector<double>& matrix,
                                     int rows, int cols)
{
    std::vector<double> result(cols * rows, 0.0);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result[j * rows + i] = matrix[i * cols + j];
    return result;
}

double Trace_Product(const std::vector<double>& a,
                     const std::vector<double>& b, int n)
{
    double result = 0.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            result += a[i * n + j] * b[j * n + i];
    return result;
}

double Frobenius_Norm(const std::vector<double>& matrix)
{
    double squared = 0.0;
    for (double value : matrix) squared += value * value;
    return std::sqrt(squared);
}

}  // namespace

int main()
{
    if (QC_SCF_DIIS_Double_Is_Finite(
            Double_From_Bits(UINT64_C(0x7ff8000000000000))) ||
        QC_SCF_DIIS_Double_Is_Finite(
            Double_From_Bits(UINT64_C(0x7ff0000000000000))))
    {
        std::fprintf(stderr,
                     "non-finite DIIS values survived fast-math checking\n");
        return EXIT_FAILURE;
    }

    int unexpected_failure_calls = 0;
    if (!QC_SCF_Require_Eigensolver_Success(
            QC_SCF_EIGENSOLVER_FOCK, QC_SCF_EIGENSOLVER_CHANNEL_ALPHA, 7,
            0, 0,
            [&](const QC_SCF_Eigensolver_Failure&)
            { ++unexpected_failure_calls; }) ||
        unexpected_failure_calls != 0 ||
        !QC_SCF_Require_Eigensolver_Workspace(
            7, 0, 128, true,
            [&](const QC_SCF_Eigensolver_Workspace_Failure&)
            { ++unexpected_failure_calls; }) ||
        unexpected_failure_calls != 0)
    {
        std::fprintf(stderr, "successful eigensolver state was rejected\n");
        return EXIT_FAILURE;
    }

    if (!Eigensolver_Failure_Matches(
            QC_SCF_EIGENSOLVER_OVERLAP,
            QC_SCF_EIGENSOLVER_CHANNEL_OVERLAP, 19, 0, 4,
            "overlap orthogonalization", "overlap") ||
        !Eigensolver_Failure_Matches(
            QC_SCF_EIGENSOLVER_FOCK, QC_SCF_EIGENSOLVER_CHANNEL_ALPHA, 23,
            -17, 0, "SCF Fock diagonalization", "alpha") ||
        !Eigensolver_Failure_Matches(
            QC_SCF_EIGENSOLVER_FOCK, QC_SCF_EIGENSOLVER_CHANNEL_BETA, 23, 0,
            8, "SCF Fock diagonalization", "beta") ||
        !Eigensolver_Failure_Matches(
            QC_SCF_EIGENSOLVER_RI_LOEWNER,
            QC_SCF_EIGENSOLVER_CHANNEL_AUXILIARY, 31, -29, -5,
            "RI Coulomb-metric Loewner factorization", "auxiliary") ||
        !Eigensolver_Failure_Matches(
            QC_SCF_EIGENSOLVER_ENSEMBLE_OCCUPATION,
            QC_SCF_EIGENSOLVER_CHANNEL_BETA, 17, 0, 3,
            "ensemble natural-occupation diagonalization", "beta"))
    {
        std::fprintf(stderr,
                     "eigensolver failure metadata was not propagated\n");
        return EXIT_FAILURE;
    }

    if (!Workspace_Failure_Matches(37, -41, 256, true) ||
        !Workspace_Failure_Matches(37, 0, 0, true) ||
        !Workspace_Failure_Matches(37, 0, -9, true) ||
        !Workspace_Failure_Matches(37, 0, 256, false))
    {
        std::fprintf(stderr,
                     "eigensolver workspace failure was not propagated\n");
        return EXIT_FAILURE;
    }

    // Alpha alone would give (1/2, 1/2).  Adding the beta commutators changes
    // the shared UKS optimum to (2/103, 101/103).
    const std::vector<double> alpha_gram = {1.0, 0.0, 0.0, 1.0};
    const std::vector<double> beta_gram = {100.0, 0.0, 0.0, 1.0};
    std::vector<double> joint_gram(4);
    for (int i = 0; i < 4; ++i)
        joint_gram[i] = alpha_gram[i] + beta_gram[i];
    std::vector<double> coefficients;
    if (!QC_SCF_Solve_CDIIS_Coefficients(joint_gram, 2, 0.0,
                                          coefficients) ||
        !Close(coefficients[0], 2.0 / 103.0) ||
        !Close(coefficients[1], 101.0 / 103.0))
    {
        std::fprintf(stderr, "coupled CDIIS coefficients are incorrect\n");
        return EXIT_FAILURE;
    }

    if (QC_SCF_DIIS_Coefficients_Are_Usable({100.0, -99.0}) ||
        QC_SCF_Solve_CDIIS_Coefficients({1.0, 1.0, 1.0, 1.0}, 2, 0.0,
                                         coefficients))
    {
        std::fprintf(stderr, "unsafe/singular CDIIS solve was accepted\n");
        return EXIT_FAILURE;
    }

    const QC_SCF_DIIS_Extrapolation_Plan small_residual_plan =
        QC_SCF_Select_DIIS_Extrapolation_Plan(1.0e-4, 1.0e-1);
    const QC_SCF_DIIS_Extrapolation_Plan large_residual_plan =
        QC_SCF_Select_DIIS_Extrapolation_Plan(1.0, 1.0e-1);
    if (!small_residual_plan.try_cdiis_first ||
        !small_residual_plan.try_bounded_simplex ||
        large_residual_plan.try_cdiis_first ||
        !large_residual_plan.try_bounded_simplex)
    {
        std::fprintf(stderr,
                     "CDIIS failure cannot reach the bounded simplex "
                     "fallback\n");
        return EXIT_FAILURE;
    }

    // The old clip-and-renormalize iteration returned (0.68, 0.32, 0) for
    // this convex problem even though the exact boundary solution is
    // (0.8, 0.2, 0).  The active-face solver must satisfy the full KKT system.
    const std::vector<double> identity_qp = {
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0,
    };
    if (!QC_SCF_Solve_Simplex_QP(
            3, identity_qp, {-0.8, -0.2, 1.0}, coefficients) ||
        !Close(coefficients[0], 0.8) || !Close(coefficients[1], 0.2) ||
        !Close(coefficients[2], 0.0))
    {
        std::fprintf(stderr, "simplex QP boundary/KKT solution is incorrect\n");
        return EXIT_FAILURE;
    }

    // A completely flat (and hence singular on every multi-point face)
    // objective is valid.  Deterministic tie-breaking must retain the newest
    // history point instead of jumping to an arbitrary old density.
    if (!QC_SCF_Solve_Simplex_QP(
            3, std::vector<double>(9, 0.0),
            std::vector<double>(3, 0.0), coefficients) ||
        !Close(coefficients[0], 0.0) || !Close(coefficients[1], 0.0) ||
        !Close(coefficients[2], 1.0))
    {
        std::fprintf(stderr, "degenerate simplex QP handling is incorrect\n");
        return EXIT_FAILURE;
    }

    // For DF={{1,-1/2},{-1/2,0}}, ADIIS is -c0+2*c0^2 and therefore has
    // c0=1/4.  Passing B instead of B+B^T to a 0.5*c^T*H*c solver gives the
    // old, incorrect c0=1/2 result.
    std::vector<double> adiis_hessian;
    std::vector<double> adiis_linear;
    if (!QC_SCF_Build_ADIIS_QP({1.0, -0.5, -0.5, 0.0}, 2,
                               adiis_hessian, adiis_linear) ||
        !Close(adiis_hessian[0], 4.0) ||
        !Close(adiis_hessian[1], adiis_hessian[2]) ||
        !QC_SCF_Solve_Simplex_QP(2, adiis_hessian, adiis_linear,
                                 coefficients) ||
        !Close(coefficients[0], 0.25) || !Close(coefficients[1], 0.75))
    {
        std::fprintf(stderr, "ADIIS symmetric Hessian/factor is incorrect\n");
        return EXIT_FAILURE;
    }

    // A nonlinear UKS response need not give an exactly symmetric B.  The
    // scalar quadratic depends on its symmetric part, so the QP Hessian must
    // contain B+B^T, not just 2*B or one triangle.
    if (!QC_SCF_Build_ADIIS_QP(
            {0.0, 1.0, 0.0,
             0.0, 0.0, 0.0,
             0.0, 0.0, 0.0},
            3, adiis_hessian, adiis_linear) ||
        !Close(adiis_hessian[1], 1.0) ||
        !Close(adiis_hessian[3], 1.0))
    {
        std::fprintf(stderr,
                     "non-symmetric ADIIS response was not symmetrized\n");
        return EXIT_FAILURE;
    }

    // Exhaustive face comparison must also handle an indefinite Hessian; the
    // global solution of this problem is the first vertex, not its interior
    // stationary point.
    if (!QC_SCF_Solve_Simplex_QP(2, {-2.0, 0.0, 0.0, 2.0},
                                 {0.0, 0.0}, coefficients) ||
        !Close(coefficients[0], 1.0) || !Close(coefficients[1], 0.0))
    {
        std::fprintf(stderr,
                     "indefinite simplex QP global minimum is incorrect\n");
        return EXIT_FAILURE;
    }

    coefficients = {1.0};
    if (QC_SCF_Solve_Simplex_QP(
            QC_SCF_SIMPLEX_QP_MAX_DIMENSION + 1,
            std::vector<double>(
                (QC_SCF_SIMPLEX_QP_MAX_DIMENSION + 1) *
                    (QC_SCF_SIMPLEX_QP_MAX_DIMENSION + 1),
                0.0),
            std::vector<double>(QC_SCF_SIMPLEX_QP_MAX_DIMENSION + 1, 0.0),
            coefficients) ||
        !coefficients.empty())
    {
        std::fprintf(stderr, "oversized simplex QP did not request fallback\n");
        return EXIT_FAILURE;
    }

    const double quiet_nan =
        Double_From_Bits(UINT64_C(0x7ff8000000000000));
    const double positive_infinity =
        Double_From_Bits(UINT64_C(0x7ff0000000000000));
    if (QC_SCF_Solve_Simplex_QP(
            2, {1.0, quiet_nan, quiet_nan, 1.0}, {0.0, 0.0},
            coefficients) ||
        QC_SCF_Solve_Simplex_QP(2, {1.0, 0.0, 0.0, 1.0},
                                {0.0, positive_infinity}, coefficients) ||
        QC_SCF_Build_ADIIS_QP({1.0, 0.0, quiet_nan, 0.0}, 2,
                              adiis_hessian, adiis_linear))
    {
        std::fprintf(stderr,
                     "non-finite simplex/ADIIS data survived validation\n");
        return EXIT_FAILURE;
    }

    constexpr double configured_shift = 1.5;
    QC_SCF_Convergence_State state;
    double expected_shift = configured_shift;
    QC_SCF_Convergence_Decision decision;
    for (int stage = 0;
         stage < QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT; ++stage)
    {
        if (!Close(QC_SCF_Level_Shift_For_Stage(configured_shift, stage),
                   expected_shift) ||
            !Close(QC_SCF_Level_Shift_For_Map(configured_shift, state, false),
                   expected_shift))
        {
            std::fprintf(stderr,
                         "level-shift continuation stage is incorrect\n");
            return EXIT_FAILURE;
        }
        decision = QC_SCF_Observe_Iteration(
            state, false, true, 1.0e-8, 1.0e-8, 1.0e-6, 1.0e-6,
            configured_shift);
        if (decision.converged || !decision.reset_diis_history ||
            !decision.advanced_level_shift ||
            state.level_shift_stage != stage + 1 ||
            state.verifying_physical_fixed_point !=
                (stage + 1 == QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT))
        {
            std::fprintf(stderr,
                         "stable shifted stage did not advance safely\n");
            return EXIT_FAILURE;
        }

        // A lower-stage residual rebound must wait at that stage.  It must
        // neither restore the previous shift nor permanently latch shifting
        // on after the finite continuation sequence reaches zero.
        if (stage == 0)
        {
            decision = QC_SCF_Observe_Iteration(
                state, false, true, 1.0e-2, 1.0e-2, 1.0e-6, 1.0e-6,
                configured_shift);
            if (decision.reset_diis_history ||
                decision.advanced_level_shift ||
                state.level_shift_stage != 1 ||
                !Close(QC_SCF_Level_Shift_For_Map(configured_shift, state,
                                                  false),
                       configured_shift * QC_SCF_LEVEL_SHIFT_STAGE_RATIO))
            {
                std::fprintf(stderr,
                             "continuation rebound changed shift stage\n");
                return EXIT_FAILURE;
            }
        }
        expected_shift *= QC_SCF_LEVEL_SHIFT_STAGE_RATIO;
    }

    if (!Close(QC_SCF_Level_Shift_For_Map(configured_shift, state, false),
               0.0) ||
        !state.verifying_physical_fixed_point)
    {
        std::fprintf(stderr,
                     "finite continuation did not request a raw zero-shift "
                     "map\n");
        return EXIT_FAILURE;
    }

    // With shifting disabled there is no shifted candidate, so a stable
    // unshifted accelerated map still has to request the same raw physical
    // verification explicitly.
    QC_SCF_Convergence_State unshifted_state;
    decision = QC_SCF_Observe_Iteration(unshifted_state, false, true, 1.0e-8,
                                        1.0e-8, 1.0e-6, 1.0e-6, 0.0);
    if (decision.converged || decision.advanced_level_shift ||
        !unshifted_state.verifying_physical_fixed_point)
    {
        std::fprintf(stderr,
                     "unshifted accelerated candidate handling is incorrect\n");
        return EXIT_FAILURE;
    }

    // The first physical dE compares against an accelerated state and is not
    // a physical-to-physical quantity.  The unshifted, non-extrapolated
    // density map itself is the fixed-point oracle and must be sufficient.
    decision = QC_SCF_Observe_Iteration(state, true, true, 1.0, 1.0e-8,
                                        1.0e-6, 1.0e-6,
                                        configured_shift);
    if (!decision.converged || state.physical_streak != 1)
    {
        std::fprintf(stderr, "physical fixed point was not confirmed\n");
        return EXIT_FAILURE;
    }

    QC_SCF_Convergence_State false_candidate;
    false_candidate.level_shift_stage =
        QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT;
    false_candidate.verifying_physical_fixed_point = true;
    decision = QC_SCF_Observe_Iteration(false_candidate, true, true, 0.0,
                                        1.0e-2, 1.0e-6, 1.0e-6,
                                        configured_shift);
    if (!decision.reset_diis_history ||
        false_candidate.verifying_physical_fixed_point ||
        false_candidate.level_shift_stage !=
            QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT ||
        !Close(QC_SCF_Level_Shift_For_Map(configured_shift, false_candidate,
                                          false),
               0.0))
    {
        std::fprintf(stderr, "false physical candidate reused DIIS history\n");
        return EXIT_FAILURE;
    }

    if (!Close(QC_SCF_Level_Shift_Density_Factor(false), 0.5) ||
        !Close(QC_SCF_Level_Shift_Density_Factor(true), 1.0) ||
        !Close(QC_SCF_Level_Shift_For_Map(configured_shift, state, true),
               0.0) ||
        !Close(QC_SCF_Level_Shift_For_Map(0.0, state, false), 0.0) ||
        !Close(QC_SCF_Level_Shift_For_Map(
                   Double_From_Bits(UINT64_C(0x7ff8000000000000)), state,
                   false),
               0.0) ||
        !Close(QC_SCF_Level_Shift_For_Stage(configured_shift, -1), 0.0))
    {
        std::fprintf(stderr, "level-shift policy is incorrect\n");
        return EXIT_FAILURE;
    }

    if (QC_SCF_Should_Use_Incremental_Fock(0, false) ||
        QC_SCF_Should_Use_Incremental_Fock(1, false) ||
        !QC_SCF_Should_Use_Incremental_Fock(2, false) ||
        QC_SCF_Should_Use_Incremental_Fock(200, true))
    {
        std::fprintf(stderr,
                     "physical Fock confirmation did not force a full "
                     "rebuild\n");
        return EXIT_FAILURE;
    }

    QC_SCF_Ensemble_Bracket ensemble_bracket;
    QC_SCF_Ensemble_Root_Decision ensemble_decision =
        QC_SCF_Start_Ensemble_Root(-0.6, 0.09, 1.4, 0.49, 0.25, 1.0e-8,
                                   ensemble_bracket);
    if (ensemble_decision.status != QC_SCF_ENSEMBLE_ROOT_EVALUATE ||
        !Close(ensemble_decision.next_fraction, 0.3))
    {
        std::fprintf(stderr, "ensemble root bracket was not constructed\n");
        return EXIT_FAILURE;
    }

    QC_SCF_Ensemble_Bracket endpoint_root_bracket;
    endpoint_root_bracket.lower_fraction = 0.0;
    endpoint_root_bracket.upper_fraction = 1.0;
    endpoint_root_bracket.lower_derivative = -0.003;
    endpoint_root_bracket.upper_derivative = 0.997;
    if (!Close(QC_SCF_Ensemble_Safeguarded_Secant(endpoint_root_bracket),
               0.003))
    {
        std::fprintf(stderr,
                     "endpoint-near secant was clipped by a fixed margin\n");
        return EXIT_FAILURE;
    }
    endpoint_root_bracket.upper_derivative = quiet_nan;
    const double midpoint_fallback =
        QC_SCF_Ensemble_Safeguarded_Secant(endpoint_root_bracket);
    if (!Close(midpoint_fallback, 0.5) ||
        !(midpoint_fallback > endpoint_root_bracket.lower_fraction) ||
        !(midpoint_fallback < endpoint_root_bracket.upper_fraction))
    {
        std::fprintf(stderr,
                     "invalid secant did not use a strict-interior midpoint\n");
        return EXIT_FAILURE;
    }

    // An internal exact line root can commit immediately.  Commit retains the
    // already sampled P verbatim while updating only its active weights, so
    // the sample's raw F/E/P_new can be used for the first global KKT
    // observation.  The extra same-P raw build remains reserved for the final
    // global certificate.
    ensemble_decision = QC_SCF_Observe_Ensemble_Trial(
        ensemble_bracket, 0.3, 0.0, 0.0, 1.0e-8, 1.0e-8, 1.0e-6);
    if (ensemble_decision.status != QC_SCF_ENSEMBLE_ROOT_STATIONARY)
    {
        std::fprintf(stderr, "ensemble root was not confirmed\n");
        return EXIT_FAILURE;
    }

    QC_SCF_Ensemble_Bracket rejected_bracket;
    if (QC_SCF_Start_Ensemble_Root(-0.6, 0.09, -0.1, 0.01, 0.25,
                                   1.0e-8, rejected_bracket)
            .status != QC_SCF_ENSEMBLE_ROOT_REJECTED ||
        QC_SCF_Start_Ensemble_Root(
            quiet_nan, 0.09, 1.4, 0.49, 0.25, 1.0e-8,
            rejected_bracket)
                .status != QC_SCF_ENSEMBLE_ROOT_REJECTED)
    {
        std::fprintf(stderr,
                     "invalid/non-bracketing ensemble candidate survived\n");
        return EXIT_FAILURE;
    }

    ensemble_decision = QC_SCF_Start_Ensemble_Root(
        -0.6, 0.09, 1.4, 0.49, 0.25, 1.0e-8, rejected_bracket);
    ensemble_decision = QC_SCF_Observe_Ensemble_Trial(
        rejected_bracket, ensemble_decision.next_fraction, 0.0, 0.2,
        1.0e-8, 1.0e-8, 1.0e-6);
    if (ensemble_decision.status == QC_SCF_ENSEMBLE_ROOT_STATIONARY)
    {
        std::fprintf(stderr,
                     "ensemble root without an interior energy minimum was "
                     "accepted\n");
        return EXIT_FAILURE;
    }

    // Pairwise Frank-Wolfe lines end at the away atom's active weight, not
    // necessarily at one.  The generalized bracket must preserve that bound.
    QC_SCF_Ensemble_Bracket bounded_bracket;
    ensemble_decision = QC_SCF_Start_Ensemble_Root(
        -1.0, 0.0, 1.0, 0.1, 0.2, 1.0e-8, bounded_bracket, false, 0.25);
    if (ensemble_decision.status != QC_SCF_ENSEMBLE_ROOT_EVALUATE ||
        !Close(ensemble_decision.next_fraction, 0.125) ||
        !Close(bounded_bracket.upper_fraction, 0.25))
    {
        std::fprintf(stderr, "bounded pairwise root was not constructed\n");
        return EXIT_FAILURE;
    }

    // Even a bracket narrower than the density tolerance accepts a genuine
    // stationary root; its raw repeat occurs after commit.
    QC_SCF_Ensemble_Bracket narrow_bracket;
    ensemble_decision = QC_SCF_Start_Ensemble_Root(
        -1.0, 0.0, 1.0, 0.1, 0.2, 1.0e-8, narrow_bracket, false, 1.0e-7);
    ensemble_decision = QC_SCF_Observe_Ensemble_Trial(
        narrow_bracket, ensemble_decision.next_fraction, 0.0, -0.01,
        1.0e-8, 1.0e-8, 1.0e-6);
    if (ensemble_decision.status != QC_SCF_ENSEMBLE_ROOT_STATIONARY)
    {
        std::fprintf(stderr, "narrow root repeat was not accepted\n");
        return EXIT_FAILURE;
    }

    // At float-density resolution, an energy-inconsistent bracket contains
    // no further information.  It must reject the current line immediately
    // instead of resampling indistinguishable mixtures until the hard cap.
    QC_SCF_Ensemble_Bracket exhausted_bracket;
    ensemble_decision = QC_SCF_Start_Ensemble_Root(
        -1.0, 0.0, 1.0, 0.1, 1.0e-7, 1.0e-8,
        exhausted_bracket);
    ensemble_decision = QC_SCF_Observe_Ensemble_Trial(
        exhausted_bracket, ensemble_decision.next_fraction, 1.0e-4, 0.1,
        1.0e-8, 1.0e-8, 1.0e-6);
    if (ensemble_decision.status != QC_SCF_ENSEMBLE_ROOT_FAILED)
    {
        std::fprintf(stderr,
                     "exhausted inconsistent root kept resampling\n");
        return EXIT_FAILURE;
    }

    // If F-based stationarity and the independently sampled energy disagree,
    // retain the actual descending sample instead of accepting a rising root
    // or looping forever at float density resolution.
    QC_SCF_Ensemble_Bracket inexact_bracket;
    ensemble_decision = QC_SCF_Start_Ensemble_Root(
        -1.0, 0.0, 1.0, 0.1, 0.2, 1.0e-8, inexact_bracket);
    (void)QC_SCF_Observe_Ensemble_Trial(
        inexact_bracket, 0.4, -0.1, -0.01, 1.0e-8, 1.0e-8, 1.0e-6);
    ensemble_decision = QC_SCF_Observe_Ensemble_Trial(
        inexact_bracket, 0.5, 0.0, 0.01, 1.0e-8, 1.0e-8, 1.0e-6);
    if (ensemble_decision.status !=
            QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY ||
        !Close(ensemble_decision.next_fraction, 0.4))
    {
        std::fprintf(stderr, "inexact line did not retain its best sample\n");
        return EXIT_FAILURE;
    }

    std::vector<double> corrective_direction;
    if (!QC_SCF_Build_Fully_Corrective_Direction(
            {3.0, 2.0, 0.0}, {0.6, 0.4, 0.0},
            corrective_direction) ||
        corrective_direction.size() != 3 ||
        !Close(corrective_direction[0], -0.6) ||
        !Close(corrective_direction[1], -0.15) ||
        !Close(corrective_direction[2], 0.75) ||
        !Close(corrective_direction[0] + corrective_direction[1] +
                   corrective_direction[2],
               0.0) ||
        !Close(0.6 + corrective_direction[0], 0.0) ||
        !Close(0.4 + corrective_direction[1], 0.25))
    {
        std::fprintf(stderr,
                     "fully-corrective simplex direction was not feasible\n");
        return EXIT_FAILURE;
    }
    if (!QC_SCF_Build_Pairwise_Corrective_Direction(
            {10.0, -100.0, 9.0}, {1.0, 0.0, 0.0},
            corrective_direction) ||
        corrective_direction.size() != 3 ||
        !Close(corrective_direction[0], -1.0) ||
        !Close(corrective_direction[1], 1.0) ||
        !Close(corrective_direction[2], 0.0))
    {
        std::fprintf(stderr,
                     "pairwise face fallback chose an invalid vertex pair\n");
        return EXIT_FAILURE;
    }
    if (!QC_SCF_Build_Pairwise_Corrective_Direction(
            {1.0, 100.0, -1.0}, {0.5, 0.0, 0.5},
            corrective_direction) ||
        !Close(corrective_direction[0], -0.5) ||
        !Close(corrective_direction[1], 0.0) ||
        !Close(corrective_direction[2], 0.5))
    {
        std::fprintf(stderr,
                     "pairwise fallback selected a zero-weight away vertex\n");
        return EXIT_FAILURE;
    }
    if (!QC_SCF_Build_Fully_Corrective_Direction(
            {10.0, -100.0, 9.0}, {1.0, 0.0, 0.0},
            corrective_direction) ||
        corrective_direction.size() != 3 ||
        !Close(corrective_direction[0], -1.0) ||
        !Close(corrective_direction[1], 1.0) ||
        !Close(corrective_direction[2], 0.0))
    {
        std::fprintf(stderr,
                     "fully-corrective tangent-cone projection admitted an "
                     "infeasible zero-weight direction\n");
        return EXIT_FAILURE;
    }

    double active_fw_gap = 0.0;
    if (!QC_SCF_Ensemble_Active_Set_FW_Gap(
            {10.0, -100.0, 9.0}, {1.0, 0.0, 0.0}, active_fw_gap) ||
        !Close(active_fw_gap, 110.0) ||
        QC_SCF_Ensemble_Active_Set_FW_Gap(
            {10.0, quiet_nan, 9.0}, {1.0, 0.0, 0.0}, active_fw_gap))
    {
        std::fprintf(stderr, "active-hull Frank-Wolfe gap is incorrect\n");
        return EXIT_FAILURE;
    }

    // Inexact active-hull forcing only decides when to ask the global oracle
    // for another vertex.  It must track a material global gap, but collapse
    // exactly to the original final tolerance near convergence.
    double active_correction_target = 0.0;
    if (!QC_SCF_Ensemble_Active_Correction_Target(
            8.0e-4, 1.0e-6, active_correction_target) ||
        !Close(active_correction_target, 2.0e-4) ||
        // This active gap is good enough to request a new global vertex, but
        // is deliberately not a valid final active-gap certificate.
        !(1.0e-4 <= active_correction_target) || !(1.0e-4 > 1.0e-6) ||
        !QC_SCF_Ensemble_Active_Correction_Target(
            2.0e-6, 1.0e-6, active_correction_target) ||
        !Close(active_correction_target, 1.0e-6) ||
        QC_SCF_Ensemble_Active_Correction_Target(
            -1.0, 1.0e-6, active_correction_target) ||
        QC_SCF_Ensemble_Active_Correction_Target(
            1.0, 1.0e-6, active_correction_target, 1.0))
    {
        std::fprintf(stderr,
                     "active-hull forcing target changed final tolerance\n");
        return EXIT_FAILURE;
    }

    // A tiny positive weight on a very unfavorable vertex makes the first
    // projected step tiny because that step must reach the adjacent simplex
    // face without crossing it.  It is still essential: after the bad vertex
    // is dropped, the next fixed-hull direction is material.  Rejecting the
    // first direction merely because its density norm is below the final SCF
    // tolerance permanently blocks the real correction.
    const std::vector<double> tiny_face_weights =
        {1.0 - 2.0e-7, 2.0e-7, 0.0};
    if (!QC_SCF_Build_Fully_Corrective_Direction(
            {0.0, 100.0, -1.0}, tiny_face_weights,
            corrective_direction))
    {
        std::fprintf(stderr, "tiny simplex face could not be reached\n");
        return EXIT_FAILURE;
    }
    std::vector<double> pairwise_tiny_direction;
    std::vector<double> pairwise_after_tiny_drop;
    if (!QC_SCF_Build_Pairwise_Corrective_Direction(
            {0.0, 100.0, -1.0}, tiny_face_weights,
            pairwise_tiny_direction) ||
        !Close(pairwise_tiny_direction[0], 0.0) ||
        !Close(pairwise_tiny_direction[1], -2.0e-7) ||
        !Close(pairwise_tiny_direction[2], 2.0e-7) ||
        !QC_SCF_Build_Ensemble_Line_Weights(
            tiny_face_weights, pairwise_tiny_direction, 1.0,
            pairwise_after_tiny_drop) ||
        pairwise_after_tiny_drop[1] != 0.0)
    {
        std::fprintf(stderr,
                     "pairwise fallback did not exactly drop tiny face\n");
        return EXIT_FAILURE;
    }
    double tiny_direction_norm = 0.0;
    for (double value : corrective_direction)
        tiny_direction_norm += value * value;
    tiny_direction_norm = std::sqrt(tiny_direction_norm);
    std::vector<double> after_tiny_drop;
    if (!(tiny_direction_norm > 0.0) || !(tiny_direction_norm < 1.0e-6) ||
        !QC_SCF_Build_Ensemble_Line_Weights(
            tiny_face_weights, corrective_direction, 1.0,
            after_tiny_drop) ||
        after_tiny_drop[1] != 0.0 ||
        !QC_SCF_Build_Fully_Corrective_Direction(
            {0.0, 100.0, -1.0}, after_tiny_drop,
            corrective_direction))
    {
        std::fprintf(stderr, "tiny simplex face did not expose correction\n");
        return EXIT_FAILURE;
    }
    double exposed_direction_norm = 0.0;
    for (double value : corrective_direction)
        exposed_direction_norm += value * value;
    if (!(std::sqrt(exposed_direction_norm) > 0.5))
    {
        std::fprintf(stderr,
                     "face drop concealed the material corrective step\n");
        return EXIT_FAILURE;
    }

    // Trial validation and weight commit share this exact
    // clamp-and-normalize operation.  Repeating it must produce
    // bitwise-identical weights, and a clamped micro-weight must not leave the
    // simplex sum below one.
    const std::vector<double> micro_weights =
        {1.0 - 2.0e-13, 2.0e-13, 0.0};
    const std::vector<double> micro_direction = {-0.5, 0.0, 0.5};
    std::vector<double> trial_weights;
    std::vector<double> committed_weights;
    if (!QC_SCF_Build_Ensemble_Line_Weights(
            micro_weights, micro_direction, 1.0, trial_weights) ||
        !QC_SCF_Build_Ensemble_Line_Weights(
            micro_weights, micro_direction, 1.0, committed_weights) ||
        trial_weights != committed_weights || trial_weights[1] != 0.0 ||
        !Close(trial_weights[0] + trial_weights[1] + trial_weights[2],
               1.0) ||
        QC_SCF_Build_Ensemble_Line_Weights(
            micro_weights, {0.0, 0.0, 1.0e-3}, 1.0,
            committed_weights))
    {
        std::fprintf(stderr,
                     "ensemble trial/commit weights are not reproducible\n");
        return EXIT_FAILURE;
    }

    if (!QC_SCF_Ensemble_Density_Is_Same_Vertex(0.0) ||
        QC_SCF_Ensemble_Density_Is_Same_Vertex(1.0e-7) ||
        QC_SCF_Ensemble_Density_Is_Same_Vertex(quiet_nan))
    {
        std::fprintf(stderr,
                     "near-duplicate density was merged as one vertex\n");
        return EXIT_FAILURE;
    }

    // One reduced-coordinate inverse-BFGS update on a convex quadratic must
    // produce a tangent, feasible descent direction.  A failed secant update
    // is transactional and leaves the accepted inverse Hessian untouched so
    // the caller can safely use the projected-gradient fallback.
    std::vector<double> inverse_hessian;
    QC_SCF_Initialize_Corrective_Inverse_Hessian(3, inverse_hessian);
    std::vector<double> quasi_newton_direction;
    const std::vector<double> qn_weights = {0.4, 0.35, 0.25};
    const std::vector<double> qn_linear = {0.0, -0.1, 0.2};
    if (!QC_SCF_Build_Quasi_Newton_Corrective_Direction(
            qn_linear, qn_weights, {0.0, -0.3, 0.2},
            {0.5, 0.25, 0.25}, inverse_hessian,
            quasi_newton_direction))
    {
        std::fprintf(stderr, "valid inverse-BFGS update was rejected\n");
        return EXIT_FAILURE;
    }
    double qn_sum = 0.0;
    double qn_derivative = 0.0;
    for (size_t i = 0; i < qn_weights.size(); ++i)
    {
        qn_sum += quasi_newton_direction[i];
        qn_derivative += qn_linear[i] * quasi_newton_direction[i];
        if (qn_weights[i] + quasi_newton_direction[i] < -1.0e-12)
        {
            std::fprintf(stderr,
                         "inverse-BFGS direction crossed simplex boundary\n");
            return EXIT_FAILURE;
        }
    }
    if (!Close(qn_sum, 0.0) || !(qn_derivative < 0.0))
    {
        std::fprintf(stderr,
                     "inverse-BFGS direction is not tangent/descent\n");
        return EXIT_FAILURE;
    }
    const std::vector<double> accepted_inverse_hessian = inverse_hessian;
    if (QC_SCF_Build_Quasi_Newton_Corrective_Direction(
            qn_linear, qn_weights, qn_linear, qn_weights,
            inverse_hessian, quasi_newton_direction) ||
        inverse_hessian != accepted_inverse_hessian ||
        !QC_SCF_Build_Fully_Corrective_Direction(
            qn_linear, qn_weights, quasi_newton_direction))
    {
        std::fprintf(stderr,
                     "unsafe inverse-BFGS update polluted fallback state\n");
        return EXIT_FAILURE;
    }
    const std::vector<double> pairwise_failure_weights = {0.4, 0.6};
    const std::vector<double> pairwise_failure_hessian = inverse_hessian;
    if (QC_SCF_Build_Pairwise_Corrective_Direction(
            {1.0, 1.0}, pairwise_failure_weights,
            quasi_newton_direction) ||
        pairwise_failure_weights != std::vector<double>({0.4, 0.6}) ||
        inverse_hessian != pairwise_failure_hessian)
    {
        std::fprintf(stderr,
                     "failed pairwise proposal mutated committed state\n");
        return EXIT_FAILURE;
    }

    double global_fw_gap = 0.0;
    if (!QC_SCF_Ensemble_Global_FW_Gap(-5.0e-7, global_fw_gap) ||
        !Close(global_fw_gap, 5.0e-7) ||
        !QC_SCF_Ensemble_Global_FW_Gap(1.742e-6, global_fw_gap) ||
        !Close(global_fw_gap, 0.0) ||
        QC_SCF_Ensemble_Global_FW_Gap(
            Double_From_Bits(UINT64_C(0x7ff8000000000000)), global_fw_gap))
    {
        std::fprintf(stderr, "global FW gap used the wrong sign convention\n");
        return EXIT_FAILURE;
    }

    // Schedule the KKT block with the larger residual normalized by its final
    // tolerance.  The old spectral-only state must return to FW, while a
    // commutator-dominated state must not wait for exact hull convergence.
    // A final-certificate repeat never starts a new optimization step.
    if (!QC_SCF_Should_Start_Spectral_Orbital_Correction(
            5.0e-7, 5.0e-7, 2.0e-5, 1.0e-6, 1.0e-6, false) ||
        !QC_SCF_Should_Start_Spectral_Orbital_Correction(
            4.0e-6, 5.0e-6, 5.0e-5, 1.0e-6, 1.0e-6, false) ||
        QC_SCF_Should_Start_Spectral_Orbital_Correction(
            2.5e-5, 2.5e-5, 3.0e-6, 1.0e-6, 1.0e-6, false) ||
        QC_SCF_Should_Start_Spectral_Orbital_Correction(
            5.0e-7, 2.5e-5, 2.0e-5, 1.0e-6, 1.0e-6, false) ||
        QC_SCF_Should_Start_Spectral_Orbital_Correction(
            5.0e-7, 5.0e-7, 5.0e-7, 1.0e-6, 1.0e-6, false) ||
        QC_SCF_Should_Start_Spectral_Orbital_Correction(
            5.0e-7, 5.0e-7, 2.0e-5, 1.0e-6, 1.0e-6, true) ||
        QC_SCF_Should_Start_Spectral_Orbital_Correction(
            quiet_nan, 5.0e-7, 2.0e-5, 1.0e-6, 1.0e-6, false))
    {
        std::fprintf(stderr,
                     "spectral orbital oracle used the wrong trigger\n");
        return EXIT_FAILURE;
    }

    std::vector<double> corrected_occupations;
    if (!QC_SCF_Prepare_Spectral_Orbital_Occupations(
            {0.3, 0.7}, 1.0, 1.0, 1.0e-10,
            corrected_occupations) ||
        corrected_occupations.size() != 2 ||
        !Close(corrected_occupations[0], 0.7) ||
        !Close(corrected_occupations[1], 0.3) ||
        !QC_SCF_Prepare_Spectral_Orbital_Occupations(
            {-5.0e-7, 1.0000005}, 1.0, 1.0, 1.0e-6,
            corrected_occupations) ||
        !Close(corrected_occupations[0], 1.0) ||
        !Close(corrected_occupations[1], 0.0) ||
        !QC_SCF_Prepare_Spectral_Orbital_Occupations(
            {0.3, 0.7000001}, 1.0, 1.0, 1.0e-6,
            corrected_occupations) ||
        !Close(corrected_occupations[0], 0.70000005) ||
        !Close(corrected_occupations[1], 0.29999995) ||
        !QC_SCF_Prepare_Spectral_Orbital_Occupations(
            {-0.09, 0.0, 0.9}, 1.0, 1.0, 0.2,
            corrected_occupations) ||
        corrected_occupations.size() != 3 ||
        !Close(corrected_occupations[0], 0.95) ||
        !Close(corrected_occupations[1], 0.05) ||
        !Close(corrected_occupations[2], 0.0) ||
        QC_SCF_Prepare_Spectral_Orbital_Occupations(
            {-2.0e-6, 1.000002}, 1.0, 1.0, 1.0e-6,
            corrected_occupations) ||
        QC_SCF_Prepare_Spectral_Orbital_Occupations(
            {0.6, 0.3}, 1.0, 1.0, 1.0e-6,
            corrected_occupations) ||
        QC_SCF_Prepare_Spectral_Orbital_Occupations(
            {quiet_nan, 0.0}, 1.0, 1.0, 1.0e-6,
            corrected_occupations))
    {
        std::fprintf(stderr,
                     "spectral occupations lost order, bounds, or particles\n");
        return EXIT_FAILURE;
    }

    // A float AO density can be slightly off the nominal particle manifold.
    // Separating the invariant spectrum-preserving rotation from the capped-
    // simplex particle repair prevents the latter's positive chemical-
    // potential contribution from being mistaken for an uphill rotation.
    const std::vector<double> off_manifold_fock_eigenvalues = {-10.0, -1.0};
    const std::vector<double> off_manifold_fock_basis_diagonal = {
        0.7999999, 0.2000006};
    const std::vector<double> off_manifold_natural_occupations = {
        0.2000005, 0.8};
    double spectrum_preserving_change = 0.0;
    if (!QC_SCF_Spectrum_Preserving_Orbital_Rotation_Change(
            off_manifold_fock_eigenvalues,
            off_manifold_fock_basis_diagonal,
            off_manifold_natural_occupations,
            spectrum_preserving_change) ||
        !Close(spectrum_preserving_change, -9.0e-7) ||
        !QC_SCF_Prepare_Spectral_Orbital_Occupations(
            off_manifold_natural_occupations, 1.0, 1.0, 1.0e-6,
            corrected_occupations))
    {
        std::fprintf(stderr,
                     "spectrum-preserving orbital change is invalid\n");
        return EXIT_FAILURE;
    }
    double projected_endpoint_change = 0.0;
    for (size_t i = 0; i < corrected_occupations.size(); ++i)
        projected_endpoint_change +=
            off_manifold_fock_eigenvalues[i] *
            (corrected_occupations[i] -
             off_manifold_fock_basis_diagonal[i]);
    if (!(projected_endpoint_change > 1.0e-6) ||
        !QC_SCF_Spectral_Orbital_Properties_Are_Valid(
            {1.0000005}, {1.0}, spectrum_preserving_change, 0.0,
            1.0e-6, 1.0e-6) ||
        QC_SCF_Spectral_Orbital_Properties_Are_Valid(
            {1.0000005}, {1.0}, projected_endpoint_change, 0.0,
            1.0e-6, 1.0e-6) ||
        QC_SCF_Spectrum_Preserving_Orbital_Rotation_Change(
            {-1.0, -2.0}, off_manifold_fock_basis_diagonal,
            off_manifold_natural_occupations,
            spectrum_preserving_change))
    {
        std::fprintf(stderr,
                     "particle repair was conflated with orbital rotation\n");
        return EXIT_FAILURE;
    }

    const double float_ulp_at_one = std::ldexp(1.0, -23);
    const double rounds_up = 1.0 + 0.75 * float_ulp_at_one;
    const double rounds_down = 1.0 + 0.25 * float_ulp_at_one;
    const float positive_from_up =
        QC_SCF_Round_Float_For_Nonincreasing_Linear_Objective(rounds_up,
                                                               2.0);
    const float negative_from_up =
        QC_SCF_Round_Float_For_Nonincreasing_Linear_Objective(rounds_up,
                                                               -2.0);
    const float positive_from_down =
        QC_SCF_Round_Float_For_Nonincreasing_Linear_Objective(rounds_down,
                                                               2.0);
    const float negative_from_down =
        QC_SCF_Round_Float_For_Nonincreasing_Linear_Objective(rounds_down,
                                                               -2.0);
    if (2.0 * (static_cast<double>(positive_from_up) - rounds_up) > 0.0 ||
        -2.0 * (static_cast<double>(negative_from_up) - rounds_up) > 0.0 ||
        2.0 * (static_cast<double>(positive_from_down) - rounds_down) > 0.0 ||
        -2.0 * (static_cast<double>(negative_from_down) - rounds_down) >
            0.0 ||
        !Close(positive_from_up, 1.0) ||
        !Close(negative_from_up, 1.0 + float_ulp_at_one) ||
        !Close(positive_from_down, 1.0) ||
        !Close(negative_from_down, 1.0 + float_ulp_at_one))
    {
        std::fprintf(stderr,
                     "fixed-F-aware float rounding increased its objective\n");
        return EXIT_FAILURE;
    }

    // Non-orthogonal AO example with one discarded AO direction (ne=2<3).
    // C^T S C=I, while Q contains a legal off-diagonal ensemble coherence.
    // Preserve Q's natural-occupation spectrum, pair its larger eigenvalue
    // with the lower Fock orbital, and remove the generalized commutator.  The
    // fixed-F value must not increase, and every backtracking point remains
    // inside the convex ensemble domain.
    const std::vector<double> overlap = {
        4.0, 0.0, 0.0,
        0.0, 9.0, 0.0,
        0.0, 0.0, 16.0,
    };
    const std::vector<double> orbitals = {
        0.5, 0.0,
        0.0, 1.0 / 3.0,
        0.0, 0.0,
    };
    const std::vector<double> occupation_matrix = {
        0.7, 0.2,
        0.2, 0.3,
    };
    const double spectral_discriminant = std::sqrt(0.32);
    const double high_occupation = 0.5 * (1.0 + spectral_discriminant);
    const double low_occupation = 0.5 * (1.0 - spectral_discriminant);
    const std::vector<double> spectral_occupation = {
        high_occupation, 0.0,
        0.0, low_occupation,
    };
    const std::vector<double> orbitals_t =
        Matrix_Transpose(orbitals, 3, 2);
    const std::vector<double> density = Matrix_Multiply(
        Matrix_Multiply(orbitals, 3, 2, occupation_matrix, 2), 3, 2,
        orbitals_t, 3);
    const std::vector<double> corrected_density = Matrix_Multiply(
        Matrix_Multiply(orbitals, 3, 2, spectral_occupation, 2), 3, 2,
        orbitals_t, 3);
    const std::vector<double> recovered_occupation = Matrix_Multiply(
        Matrix_Multiply(
            Matrix_Multiply(
                Matrix_Multiply(orbitals_t, 2, 3, overlap, 3), 2, 3,
                density, 3),
            2, 3, overlap, 3),
        2, 3, orbitals, 2);
    const std::vector<double> fock = {
        -3.2, 0.0, 0.0,
        0.0, -1.8, 0.0,
        0.0, 0.0, 16.0,
    };
    auto generalized_commutator =
        [&](const std::vector<double>& p)
    {
        const std::vector<double> fps = Matrix_Multiply(
            Matrix_Multiply(fock, 3, 3, p, 3), 3, 3, overlap, 3);
        const std::vector<double> spf = Matrix_Multiply(
            Matrix_Multiply(overlap, 3, 3, p, 3), 3, 3, fock, 3);
        std::vector<double> result(9, 0.0);
        for (int i = 0; i < 9; ++i) result[i] = fps[i] - spf[i];
        return result;
    };
    const std::vector<double> original_commutator =
        generalized_commutator(density);
    const std::vector<double> corrected_commutator =
        generalized_commutator(corrected_density);
    constexpr double line_fraction = 0.4;
    std::vector<double> line_density(9, 0.0);
    for (int i = 0; i < 9; ++i)
        line_density[i] = density[i] +
                          line_fraction *
                              (corrected_density[i] - density[i]);
    const std::vector<double> line_commutator =
        generalized_commutator(line_density);
    bool line_commutator_scaled = true;
    for (int i = 0; i < 9; ++i)
        line_commutator_scaled =
            line_commutator_scaled &&
            Close(line_commutator[i],
                  (1.0 - line_fraction) * original_commutator[i]);
    std::vector<double> line_occupation(4, 0.0);
    for (int i = 0; i < 4; ++i)
        line_occupation[i] =
            occupation_matrix[i] +
            line_fraction *
                (spectral_occupation[i] - occupation_matrix[i]);
    const double line_determinant =
        line_occupation[0] * line_occupation[3] -
        line_occupation[1] * line_occupation[2];
    const double complement_determinant =
        (1.0 - line_occupation[0]) *
            (1.0 - line_occupation[3]) -
        line_occupation[1] * line_occupation[2];
    const double original_spectral_determinant =
        occupation_matrix[0] * occupation_matrix[3] -
        occupation_matrix[1] * occupation_matrix[2];
    const double corrected_spectral_determinant =
        spectral_occupation[0] * spectral_occupation[3];
    const double fixed_fock_change =
        Trace_Product(fock, corrected_density, 3) -
        Trace_Product(fock, density, 3);
    if (!Close(recovered_occupation[0], occupation_matrix[0]) ||
        !Close(recovered_occupation[1], occupation_matrix[1]) ||
        !Close(recovered_occupation[2], occupation_matrix[2]) ||
        !Close(recovered_occupation[3], occupation_matrix[3]) ||
        !Close(Trace_Product(overlap, density, 3), 1.0) ||
        !Close(Trace_Product(overlap, corrected_density, 3), 1.0) ||
        !Close(Trace_Product(overlap, line_density, 3), 1.0) ||
        !Close(high_occupation + low_occupation, 1.0) ||
        !Close(original_spectral_determinant,
               corrected_spectral_determinant) ||
        !(fixed_fock_change < 0.0) ||
        !(Frobenius_Norm(original_commutator) > 0.5) ||
        !Close(Frobenius_Norm(corrected_commutator), 0.0) ||
        !line_commutator_scaled || !(line_determinant > 0.0) ||
        !(complement_determinant > 0.0) ||
        !QC_SCF_Spectral_Orbital_Properties_Are_Valid(
            {Trace_Product(overlap, density, 3)},
            {Trace_Product(overlap, corrected_density, 3)},
            fixed_fock_change, Frobenius_Norm(corrected_commutator),
            1.0e-6, 1.0e-6) ||
        QC_SCF_Spectral_Orbital_Properties_Are_Valid(
            {1.0}, {1.0}, 2.0e-6, 0.0, 1.0e-6, 1.0e-6))
    {
        std::fprintf(stderr,
                     "spectral orbital correction violated its invariants\n");
        return EXIT_FAILURE;
    }

    double next_correction_fraction = 0.0;
    if (!QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
            -10.0, -10.000001, 1.0, 0.0, 2.0e-5, 8.0e-6, 1.0e-6,
            1.0e-6) ||
        !QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
            -10.0, -9.999999, 1.0, 0.0, 2.0e-5, 8.0e-6, 1.0e-6,
            1.0e-6) ||
        QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
            -10.0, -9.999999, 1.0, 0.0, 2.0e-5, 2.1e-5, 1.0e-6,
            1.0e-6) ||
        QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
            -10.0, -9.999997, 1.0, 0.0, 2.0e-5, 8.0e-6, 1.0e-6,
            1.0e-6) ||
        !QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
            -10.0, -10.000001, 1.0, 0.0, 2.0e-5, 2.1e-5, 1.0e-6,
            1.0e-6) ||
        !QC_SCF_Next_Spectral_Orbital_Fraction(
            1.0, 1, next_correction_fraction) ||
        !Close(next_correction_fraction, 0.5) ||
        QC_SCF_Next_Spectral_Orbital_Fraction(
            1.0, QC_SCF_ENSEMBLE_MAX_SPECTRAL_ORBITAL_BACKTRACKS,
            next_correction_fraction))
    {
        std::fprintf(stderr,
                     "spectral orbital transaction policy is unsafe\n");
        return EXIT_FAILURE;
    }

    // At the observed RI-stored failure, the mandatory nominal-particle
    // repair contributes +6.151012069e-6 Ha to the fixed-origin-F line.  The
    // quarter step is within the unchanged 2 microhartree guard after
    // removing x times that normal component; the half step remains outside.
    const double spectral_guard_origin = -128.790142030047;
    const double spectral_repair_component = 6.151012069160523e-6;
    if (!QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
            spectral_guard_origin, -128.790139775499, 0.25,
            spectral_repair_component, 4.0e-6, 2.0e-6, 1.0e-6,
            1.0e-6) ||
        QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
            spectral_guard_origin, -128.790136873772, 0.5,
            spectral_repair_component, 4.0e-6, 2.0e-6, 1.0e-6,
            1.0e-6))
    {
        std::fprintf(stderr,
                     "spectral repair component changed the raw energy "
                     "guard\n");
        return EXIT_FAILURE;
    }

    // A lower sample from an older line is not the acceptance baseline for a
    // later corrective step.  The later sample descends from its own committed
    // origin and must survive even when float active-set reconstruction made
    // that origin a few microhartree higher than the historical sample.
    const double historical_trial = -128.7901460;
    const double committed_origin = -128.7901410;
    const double current_trial = -128.7901415;
    if (!QC_SCF_Ensemble_Energy_Within_Line_Guard(
            committed_origin, current_trial, 1.0e-6) ||
        current_trial <= historical_trial + 2.0e-6 ||
        QC_SCF_Ensemble_Energy_Within_Line_Guard(
            committed_origin, -128.7901380, 1.0e-6))
    {
        std::fprintf(stderr,
                     "ensemble line guard reused a stale trial-energy "
                     "baseline\n");
        return EXIT_FAILURE;
    }

    // S=I: both a doubly occupied RKS density with factor 1/2 and a singly
    // occupied UKS density with factor 1 give the same virtual projector.
    const double restricted_projector = 1.0 - 0.5 * 2.0;
    const double unrestricted_projector = 1.0 - 1.0 * 1.0;
    if (!Close(restricted_projector, 0.0) ||
        !Close(unrestricted_projector, 0.0))
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
