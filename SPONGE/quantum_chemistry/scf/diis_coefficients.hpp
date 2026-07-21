#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <utility>
#include <vector>

static inline bool QC_SCF_DIIS_Double_Is_Finite(double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires a 64-bit double representation");
    std::memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return ((bits >> 52U) & 0x7ffU) != 0x7ffU;
}

static inline bool QC_SCF_DIIS_Coefficients_Are_Usable(
    const std::vector<double>& coefficients, bool require_simplex = false)
{
    if (coefficients.empty()) return false;
    double sum = 0.0;
    double l1_norm = 0.0;
    for (double coefficient : coefficients)
    {
        if (!QC_SCF_DIIS_Double_Is_Finite(coefficient)) return false;
        if (require_simplex && coefficient < -1.0e-10) return false;
        sum += coefficient;
        l1_norm += std::fabs(coefficient);
    }
    if (!QC_SCF_DIIS_Double_Is_Finite(sum) ||
        !QC_SCF_DIIS_Double_Is_Finite(l1_norm))
        return false;
    if (std::fabs(sum - 1.0) > 1.0e-8 * std::max(1.0, l1_norm)) return false;

    // Very large cancelling coefficients turn roundoff in otherwise harmless
    // history vectors into a new SCF displacement.  Reject the extrapolation;
    // the caller will remove the oldest history point and retry.
    return l1_norm <= 50.0;
}

// EDIIS/ADIIS spaces are intentionally small (normally six vectors).  An
// exhaustive face search is deterministic and, unlike clipping followed by
// renormalization, actually solves the simplex-constrained quadratic problem.
// Larger histories are reduced by the caller before retrying extrapolation.
static constexpr int QC_SCF_SIMPLEX_QP_MAX_DIMENSION = 12;

struct QC_SCF_DIIS_Extrapolation_Plan
{
    bool try_cdiis_first = false;
    bool try_bounded_simplex = true;
};

// Near self-consistency CDIIS is the preferred extrapolator, but failure of
// its unconstrained linear solve must not discard the complete history until
// the code falls back to a raw Fock map.  The EDIIS/ADIIS simplex solve is a
// bounded, finite alternative using the same history and is safe at both
// large and small residuals.
static inline QC_SCF_DIIS_Extrapolation_Plan
QC_SCF_Select_DIIS_Extrapolation_Plan(double commutator_rms,
                                      double cdiis_threshold)
{
    return {commutator_rms <= cdiis_threshold, true};
}

static inline bool QC_SCF_Solve_Dense_Linear_System(
    std::vector<double> matrix, std::vector<double> rhs, int dimension,
    std::vector<double>& solution)
{
    solution.clear();
    if (dimension <= 0 || matrix.size() != (size_t)dimension * dimension ||
        rhs.size() != (size_t)dimension)
        return false;

    const std::vector<double> original_matrix = matrix;
    const std::vector<double> original_rhs = rhs;
    double matrix_norm = 0.0;
    for (double value : matrix)
    {
        if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;
        matrix_norm = std::max(matrix_norm, std::fabs(value));
    }
    for (double value : rhs)
        if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;

    const double pivot_tolerance =
        512.0 * std::numeric_limits<double>::epsilon() * dimension *
        std::max(1.0, matrix_norm);
    for (int column = 0; column < dimension; ++column)
    {
        int pivot = column;
        double pivot_abs =
            std::fabs(matrix[(size_t)column * dimension + column]);
        for (int row = column + 1; row < dimension; ++row)
        {
            const double candidate =
                std::fabs(matrix[(size_t)row * dimension + column]);
            if (candidate > pivot_abs)
            {
                pivot = row;
                pivot_abs = candidate;
            }
        }
        if (!(pivot_abs > pivot_tolerance)) return false;
        if (pivot != column)
        {
            for (int j = 0; j < dimension; ++j)
                std::swap(matrix[(size_t)column * dimension + j],
                          matrix[(size_t)pivot * dimension + j]);
            std::swap(rhs[column], rhs[pivot]);
        }

        const double diagonal =
            matrix[(size_t)column * dimension + column];
        for (int row = column + 1; row < dimension; ++row)
        {
            const double factor =
                matrix[(size_t)row * dimension + column] / diagonal;
            if (!QC_SCF_DIIS_Double_Is_Finite(factor)) return false;
            matrix[(size_t)row * dimension + column] = 0.0;
            for (int j = column + 1; j < dimension; ++j)
            {
                double& value = matrix[(size_t)row * dimension + j];
                value -= factor * matrix[(size_t)column * dimension + j];
                if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;
            }
            rhs[row] -= factor * rhs[column];
            if (!QC_SCF_DIIS_Double_Is_Finite(rhs[row])) return false;
        }
    }

    solution.assign(dimension, 0.0);
    for (int row = dimension - 1; row >= 0; --row)
    {
        double value = rhs[row];
        for (int j = row + 1; j < dimension; ++j)
            value -= matrix[(size_t)row * dimension + j] * solution[j];
        const double diagonal = matrix[(size_t)row * dimension + row];
        if (!(std::fabs(diagonal) > pivot_tolerance)) return false;
        solution[row] = value / diagonal;
        if (!QC_SCF_DIIS_Double_Is_Finite(solution[row])) return false;
    }

    double solution_norm = 0.0;
    double residual_norm = 0.0;
    for (double value : solution)
        solution_norm = std::max(solution_norm, std::fabs(value));
    for (int row = 0; row < dimension; ++row)
    {
        double residual = -original_rhs[row];
        for (int j = 0; j < dimension; ++j)
            residual += original_matrix[(size_t)row * dimension + j] *
                        solution[j];
        if (!QC_SCF_DIIS_Double_Is_Finite(residual)) return false;
        residual_norm = std::max(residual_norm, std::fabs(residual));
    }
    const double residual_tolerance =
        1.0e-10 *
        (1.0 + matrix_norm * solution_norm * (double)dimension);
    return QC_SCF_DIIS_Double_Is_Finite(residual_tolerance) &&
           residual_norm <= residual_tolerance;
}

// Solve min 0.5*c^T*H*c + g^T*c subject to c_i >= 0 and sum(c)=1.
//
// Every non-empty simplex face is considered.  The equality-constrained KKT
// system on a nonsingular face gives its only stationary point; singular faces
// can be skipped because any flat family of stationary points reaches a proper
// sub-face, which is enumerated independently.  Comparing all primal/dual KKT
// candidates (including every vertex) therefore finds the global quadratic
// minimum, even when the symmetric Hessian is indefinite.  All arithmetic is
// scaled and revalidated so a failed/ill-conditioned face is never published.
static inline bool QC_SCF_Solve_Simplex_QP(
    int dimension, const std::vector<double>& hessian,
    const std::vector<double>& linear, std::vector<double>& coefficients)
{
    coefficients.clear();
    const int m = dimension;
    if (m <= 0 || m > QC_SCF_SIMPLEX_QP_MAX_DIMENSION ||
        hessian.size() != (size_t)m * m || linear.size() != (size_t)m)
        return false;

    std::vector<double> symmetric_hessian((size_t)m * m, 0.0);
    double problem_scale = 0.0;
    for (int i = 0; i < m; ++i)
    {
        if (!QC_SCF_DIIS_Double_Is_Finite(linear[i])) return false;
        problem_scale = std::max(problem_scale, std::fabs(linear[i]));
        for (int j = 0; j < m; ++j)
        {
            const double hij = hessian[(size_t)i * m + j];
            const double hji = hessian[(size_t)j * m + i];
            if (!QC_SCF_DIIS_Double_Is_Finite(hij) ||
                !QC_SCF_DIIS_Double_Is_Finite(hji))
                return false;
            const double value = 0.5 * hij + 0.5 * hji;
            if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;
            symmetric_hessian[(size_t)i * m + j] = value;
            problem_scale = std::max(problem_scale, std::fabs(value));
        }
    }

    std::vector<double> scaled_hessian = symmetric_hessian;
    std::vector<double> scaled_linear = linear;
    if (problem_scale > 0.0)
    {
        for (double& value : scaled_hessian) value /= problem_scale;
        for (double& value : scaled_linear) value /= problem_scale;
    }

    const double primal_tolerance = 1.0e-10 * (m + 1.0);
    const double kkt_tolerance = 1.0e-10 * (m + 1.0);
    const std::uint64_t face_end = UINT64_C(1) << m;
    bool have_best = false;
    double best_objective = 0.0;
    double best_newest_distance = 0.0;
    std::vector<double> best;

    for (std::uint64_t face = 1; face < face_end; ++face)
    {
        std::vector<int> active;
        for (int i = 0; i < m; ++i)
            if ((face & (UINT64_C(1) << i)) != 0) active.push_back(i);

        const int active_count = (int)active.size();
        const int system_dimension = active_count + 1;
        std::vector<double> matrix(
            (size_t)system_dimension * system_dimension, 0.0);
        std::vector<double> rhs(system_dimension, 0.0);
        for (int row = 0; row < active_count; ++row)
        {
            const int i = active[row];
            rhs[row] = -scaled_linear[i];
            for (int column = 0; column < active_count; ++column)
            {
                const int j = active[column];
                matrix[(size_t)row * system_dimension + column] =
                    scaled_hessian[(size_t)i * m + j];
            }
            matrix[(size_t)row * system_dimension + active_count] = 1.0;
            matrix[(size_t)active_count * system_dimension + row] = 1.0;
        }
        rhs[active_count] = 1.0;

        std::vector<double> solution;
        if (!QC_SCF_Solve_Dense_Linear_System(
                std::move(matrix), std::move(rhs), system_dimension, solution))
            continue;

        std::vector<double> candidate(m, 0.0);
        bool primal_ok = true;
        for (int row = 0; row < active_count; ++row)
        {
            double value = solution[row];
            if (!QC_SCF_DIIS_Double_Is_Finite(value) ||
                value < -primal_tolerance)
            {
                primal_ok = false;
                break;
            }
            if (std::fabs(value) <= primal_tolerance) value = 0.0;
            candidate[active[row]] = value;
        }
        if (!primal_ok) continue;

        double coefficient_sum = 0.0;
        for (double value : candidate) coefficient_sum += value;
        if (!QC_SCF_DIIS_Double_Is_Finite(coefficient_sum) ||
            !(coefficient_sum > 0.0) ||
            std::fabs(coefficient_sum - 1.0) > 4.0 * primal_tolerance)
            continue;
        for (double& value : candidate) value /= coefficient_sum;

        std::vector<double> gradient(m, 0.0);
        int support_count = 0;
        double active_gradient = 0.0;
        bool finite_gradient = true;
        for (int i = 0; i < m; ++i)
        {
            double value = scaled_linear[i];
            for (int j = 0; j < m; ++j)
                value += scaled_hessian[(size_t)i * m + j] * candidate[j];
            if (!QC_SCF_DIIS_Double_Is_Finite(value))
            {
                finite_gradient = false;
                break;
            }
            gradient[i] = value;
            if (candidate[i] > primal_tolerance)
            {
                active_gradient += value;
                ++support_count;
            }
        }
        if (!finite_gradient || support_count == 0) continue;
        active_gradient /= support_count;

        bool kkt_ok = true;
        for (int i = 0; i < m; ++i)
        {
            if (candidate[i] > primal_tolerance)
            {
                if (std::fabs(gradient[i] - active_gradient) > kkt_tolerance)
                {
                    kkt_ok = false;
                    break;
                }
            }
            else if (gradient[i] < active_gradient - kkt_tolerance)
            {
                kkt_ok = false;
                break;
            }
        }
        if (!kkt_ok ||
            !QC_SCF_DIIS_Coefficients_Are_Usable(candidate, true))
            continue;

        double objective = 0.0;
        for (int i = 0; i < m; ++i)
        {
            objective += scaled_linear[i] * candidate[i];
            for (int j = 0; j < m; ++j)
                objective += 0.5 * candidate[i] *
                             scaled_hessian[(size_t)i * m + j] * candidate[j];
        }
        if (!QC_SCF_DIIS_Double_Is_Finite(objective)) continue;

        double newest_distance = 0.0;
        for (int i = 0; i < m; ++i)
        {
            const double target = (i == m - 1) ? 1.0 : 0.0;
            const double delta = candidate[i] - target;
            newest_distance += delta * delta;
        }

        bool replace = !have_best;
        if (have_best)
        {
            const double objective_tolerance =
                1.0e-12 *
                (1.0 + std::fabs(objective) + std::fabs(best_objective));
            if (objective < best_objective - objective_tolerance)
            {
                replace = true;
            }
            else if (std::fabs(objective - best_objective) <=
                     objective_tolerance)
            {
                const double distance_tolerance = 1.0e-12;
                if (newest_distance <
                    best_newest_distance - distance_tolerance)
                {
                    replace = true;
                }
                else if (std::fabs(newest_distance - best_newest_distance) <=
                         distance_tolerance)
                {
                    for (int i = m - 1; i >= 0; --i)
                    {
                        if (candidate[i] > best[i] + 1.0e-14)
                        {
                            replace = true;
                            break;
                        }
                        if (candidate[i] < best[i] - 1.0e-14) break;
                    }
                }
            }
        }
        if (replace)
        {
            have_best = true;
            best_objective = objective;
            best_newest_distance = newest_distance;
            best = candidate;
        }
    }

    if (!have_best) return false;
    coefficients = std::move(best);
    return QC_SCF_DIIS_Coefficients_Are_Usable(coefficients, true);
}

// Build the QP for
//   E_ADIIS = E_n + 2*A^T*c + c^T*B*c.
// Since QC_SCF_Solve_Simplex_QP uses 0.5*c^T*H*c, H must be B+B^T (2B only
// when B is exactly symmetric).  Normalizing all traces by one positive scale
// preserves the minimizer and prevents overflow in the finite checks below.
static inline bool QC_SCF_Build_ADIIS_QP(
    const std::vector<double>& density_fock_trace, int history_count,
    std::vector<double>& hessian, std::vector<double>& linear)
{
    hessian.clear();
    linear.clear();
    const int m = history_count;
    if (m < 2 || density_fock_trace.size() != (size_t)m * m) return false;

    double trace_scale = 0.0;
    for (double value : density_fock_trace)
    {
        if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;
        trace_scale = std::max(trace_scale, std::fabs(value));
    }
    std::vector<double> trace = density_fock_trace;
    if (trace_scale > 0.0)
        for (double& value : trace) value /= trace_scale;

    const int newest = m - 1;
    const double newest_trace = trace[(size_t)newest * m + newest];
    std::vector<double> response((size_t)m * m, 0.0);
    hessian.assign((size_t)m * m, 0.0);
    linear.assign(m, 0.0);
    for (int i = 0; i < m; ++i)
    {
        linear[i] =
            2.0 * (trace[(size_t)i * m + newest] - newest_trace);
        if (!QC_SCF_DIIS_Double_Is_Finite(linear[i])) return false;
        for (int j = 0; j < m; ++j)
        {
            const double value =
                trace[(size_t)i * m + j] -
                trace[(size_t)i * m + newest] -
                trace[(size_t)newest * m + j] + newest_trace;
            if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;
            response[(size_t)i * m + j] = value;
        }
    }
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < m; ++j)
        {
            const double value = response[(size_t)i * m + j] +
                                 response[(size_t)j * m + i];
            if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;
            hessian[(size_t)i * m + j] = value;
        }
    return true;
}

// Solve min c^T B c subject to sum(c)=1 for a precomputed (possibly
// spin-summed) Pulay error Gram matrix.  B is normalized before the user's
// regularizer is applied, so qc_diis_reg has the same meaning far from and
// close to convergence.
static inline bool QC_SCF_Solve_CDIIS_Coefficients(
    const std::vector<double>& error_gram, int history_count,
    double regularization, std::vector<double>& coefficients)
{
    coefficients.clear();
    const int m = history_count;
    if (m < 2 || error_gram.size() != (size_t)m * m ||
        !QC_SCF_DIIS_Double_Is_Finite(regularization) ||
        regularization < 0.0)
        return false;

    double gram_scale = 0.0;
    for (int i = 0; i < m; ++i)
    {
        const double diagonal = error_gram[(size_t)i * m + i];
        if (!QC_SCF_DIIS_Double_Is_Finite(diagonal)) return false;
        gram_scale = std::max(gram_scale, std::fabs(diagonal));
    }
    if (!(gram_scale > std::numeric_limits<double>::min())) return false;

    const int n = m + 1;
    std::vector<double> matrix((size_t)n * n, 0.0);
    std::vector<double> rhs(n, 0.0);
    rhs[m] = -1.0;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            const double a = error_gram[(size_t)i * m + j];
            const double b = error_gram[(size_t)j * m + i];
            if (!QC_SCF_DIIS_Double_Is_Finite(a) ||
                !QC_SCF_DIIS_Double_Is_Finite(b))
                return false;
            matrix[(size_t)i * n + j] = 0.5 * (a + b) / gram_scale;
        }
        matrix[(size_t)i * n + i] += regularization;
        matrix[(size_t)i * n + m] = -1.0;
        matrix[(size_t)m * n + i] = -1.0;
    }

    const std::vector<double> original_matrix = matrix;
    const std::vector<double> original_rhs = rhs;
    double matrix_norm = 0.0;
    for (double value : matrix)
    {
        if (!QC_SCF_DIIS_Double_Is_Finite(value)) return false;
        matrix_norm = std::max(matrix_norm, std::fabs(value));
    }
    const double pivot_tolerance =
        128.0 * std::numeric_limits<double>::epsilon() * n *
        std::max(1.0, matrix_norm);

    for (int column = 0; column < n; ++column)
    {
        int pivot = column;
        double pivot_abs =
            std::fabs(matrix[(size_t)column * n + column]);
        for (int row = column + 1; row < n; ++row)
        {
            const double candidate =
                std::fabs(matrix[(size_t)row * n + column]);
            if (candidate > pivot_abs)
            {
                pivot = row;
                pivot_abs = candidate;
            }
        }
        if (!(pivot_abs > pivot_tolerance)) return false;
        if (pivot != column)
        {
            for (int j = 0; j < n; ++j)
                std::swap(matrix[(size_t)column * n + j],
                          matrix[(size_t)pivot * n + j]);
            std::swap(rhs[column], rhs[pivot]);
        }

        const double diagonal = matrix[(size_t)column * n + column];
        for (int row = column + 1; row < n; ++row)
        {
            const double factor =
                matrix[(size_t)row * n + column] / diagonal;
            matrix[(size_t)row * n + column] = 0.0;
            for (int j = column + 1; j < n; ++j)
                matrix[(size_t)row * n + j] -=
                    factor * matrix[(size_t)column * n + j];
            rhs[row] -= factor * rhs[column];
        }
    }

    std::vector<double> solution(n, 0.0);
    for (int row = n - 1; row >= 0; --row)
    {
        double value = rhs[row];
        for (int j = row + 1; j < n; ++j)
            value -= matrix[(size_t)row * n + j] * solution[j];
        const double diagonal = matrix[(size_t)row * n + row];
        if (!(std::fabs(diagonal) > pivot_tolerance)) return false;
        solution[row] = value / diagonal;
        if (!QC_SCF_DIIS_Double_Is_Finite(solution[row])) return false;
    }

    double solution_norm = 0.0;
    double residual_norm = 0.0;
    for (double value : solution)
        solution_norm = std::max(solution_norm, std::fabs(value));
    for (int row = 0; row < n; ++row)
    {
        double residual = -original_rhs[row];
        for (int j = 0; j < n; ++j)
            residual += original_matrix[(size_t)row * n + j] * solution[j];
        residual_norm = std::max(residual_norm, std::fabs(residual));
    }
    const double residual_tolerance =
        1.0e-9 * (1.0 + matrix_norm * solution_norm * n);
    if (!QC_SCF_DIIS_Double_Is_Finite(residual_norm) ||
        residual_norm > residual_tolerance)
        return false;

    coefficients.assign(solution.begin(), solution.begin() + m);
    return QC_SCF_DIIS_Coefficients_Are_Usable(coefficients);
}
