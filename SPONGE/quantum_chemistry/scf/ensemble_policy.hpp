#pragma once

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

static inline bool QC_SCF_Ensemble_Double_Is_Finite(double value)
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

// Compare a line-search sample only with the raw energy of that line's
// committed origin.  An older sample from another line must never become a
// permanent threshold: reconstructing an active-set density in float can
// change its independently accumulated energy by a few microhartree even when
// the represented ensemble is unchanged.
static inline bool QC_SCF_Ensemble_Energy_Within_Line_Guard(
    double origin_energy, double sample_energy, double energy_tolerance,
    double tolerance_multiplier = 2.0)
{
    return QC_SCF_Ensemble_Double_Is_Finite(origin_energy) &&
           QC_SCF_Ensemble_Double_Is_Finite(sample_energy) &&
           QC_SCF_Ensemble_Double_Is_Finite(energy_tolerance) &&
           QC_SCF_Ensemble_Double_Is_Finite(tolerance_multiplier) &&
           energy_tolerance > 0.0 && tolerance_multiplier >= 2.0 &&
           sample_energy <=
               origin_energy + tolerance_multiplier * energy_tolerance;
}

// The natural occupations are recovered from a float AO density.  First reject
// any material bound or particle-count violation.  Then remove only the
// admitted representation noise by projecting onto the nearest capped simplex
// with the nominal particle count.  This prevents tolerance-sized trace and
// bound errors from becoming the input to the next orbital correction and
// accumulating across accepted float trials.  The projected occupations are
// assigned in descending order to the ascending Fock orbitals.
static inline bool QC_SCF_Prepare_Spectral_Orbital_Occupations(
    const std::vector<double>& raw_eigenvalues,
    double expected_particle_count, double maximum_occupation,
    double tolerance, std::vector<double>& occupations)
{
    occupations.clear();
    if (raw_eigenvalues.empty() ||
        !QC_SCF_Ensemble_Double_Is_Finite(expected_particle_count) ||
        !QC_SCF_Ensemble_Double_Is_Finite(maximum_occupation) ||
        !QC_SCF_Ensemble_Double_Is_Finite(tolerance) ||
        expected_particle_count < 0.0 || !(maximum_occupation > 0.0) ||
        !(tolerance > 0.0))
        return false;

    const double capacity =
        maximum_occupation * static_cast<double>(raw_eigenvalues.size());
    const double particle_tolerance =
        tolerance * std::max(1.0, expected_particle_count);
    if (!QC_SCF_Ensemble_Double_Is_Finite(capacity) ||
        expected_particle_count > capacity)
        return false;

    double raw_sum = 0.0;
    for (size_t i = 0; i < raw_eigenvalues.size(); ++i)
    {
        const double value = raw_eigenvalues[i];
        if (!QC_SCF_Ensemble_Double_Is_Finite(value) ||
            value < -tolerance || value > maximum_occupation + tolerance)
        {
            occupations.clear();
            return false;
        }
        raw_sum += value;
    }
    if (!QC_SCF_Ensemble_Double_Is_Finite(raw_sum) ||
        std::fabs(raw_sum - expected_particle_count) > particle_tolerance)
    {
        occupations.clear();
        return false;
    }

    occupations.assign(raw_eigenvalues.size(), 0.0);
    if (expected_particle_count == 0.0) return true;
    if (expected_particle_count == capacity)
    {
        std::fill(occupations.begin(), occupations.end(),
                  maximum_occupation);
        return true;
    }

    // The Euclidean projection onto {0 <= n_i <= maximum_occupation,
    // sum(n_i) = expected_particle_count} has the unique water-filling form
    // n_i = clamp(raw_i - shift, 0, maximum_occupation).  Its summed
    // occupation is continuous and monotone in shift, so a bracketed binary
    // search finds the common KKT shift without order-dependent clipping.
    const auto projected_sum = [&](double shift)
    {
        double sum = 0.0;
        for (double value : raw_eigenvalues)
            sum += std::max(
                0.0, std::min(maximum_occupation, value - shift));
        return sum;
    };
    const auto raw_bounds = std::minmax_element(raw_eigenvalues.begin(),
                                                 raw_eigenvalues.end());
    double lower_shift = *raw_bounds.first - maximum_occupation;
    double upper_shift = *raw_bounds.second;
    for (int iteration = 0; iteration < 128; ++iteration)
    {
        const double shift =
            lower_shift + 0.5 * (upper_shift - lower_shift);
        if (shift == lower_shift || shift == upper_shift) break;
        if (projected_sum(shift) > expected_particle_count)
            lower_shift = shift;
        else
            upper_shift = shift;
    }
    const double shift =
        lower_shift + 0.5 * (upper_shift - lower_shift);
    double projected_particle_count = 0.0;
    for (size_t i = 0; i < raw_eigenvalues.size(); ++i)
    {
        occupations[i] = std::max(
            0.0,
            std::min(maximum_occupation, raw_eigenvalues[i] - shift));
        projected_particle_count += occupations[i];
    }

    // Binary search is already at double resolution.  Repair its final few
    // arithmetic ULPs without changing any active bound; this is numerical
    // closure of the KKT solution, not a second projection heuristic.
    double residual = expected_particle_count - projected_particle_count;
    if (residual > 0.0)
    {
        for (double& value : occupations)
        {
            const double adjustment =
                std::min(residual, maximum_occupation - value);
            value += adjustment;
            residual -= adjustment;
            if (!(residual > 0.0)) break;
        }
    }
    else if (residual < 0.0)
    {
        for (double& value : occupations)
        {
            const double adjustment = std::min(-residual, value);
            value -= adjustment;
            residual += adjustment;
            if (!(residual < 0.0)) break;
        }
    }

    projected_particle_count = 0.0;
    for (double value : occupations)
    {
        if (!QC_SCF_Ensemble_Double_Is_Finite(value) || value < 0.0 ||
            value > maximum_occupation)
        {
            occupations.clear();
            return false;
        }
        projected_particle_count += value;
    }
    const double projection_roundoff =
        64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, capacity);
    if (!QC_SCF_Ensemble_Double_Is_Finite(projected_particle_count) ||
        std::fabs(projected_particle_count - expected_particle_count) >
            projection_roundoff)
    {
        occupations.clear();
        return false;
    }

    std::sort(occupations.begin(), occupations.end(),
              [](double left, double right) { return left > right; });
    return true;
}

// In the Fock eigenbasis, a spectrum-preserving orbital rotation replaces
// the diagonal of the current occupation matrix N by the same eigenvalues in
// descending order against ascending orbital energies.  Evaluating the
// change term by term avoids subtracting two large AO-space F:P values.
static inline bool QC_SCF_Spectrum_Preserving_Orbital_Rotation_Change(
    const std::vector<double>& ascending_fock_eigenvalues,
    const std::vector<double>& fock_basis_density_diagonal,
    const std::vector<double>& raw_natural_occupations,
    double& linear_change)
{
    linear_change = 0.0;
    const size_t count = ascending_fock_eigenvalues.size();
    if (count == 0 || fock_basis_density_diagonal.size() != count ||
        raw_natural_occupations.size() != count)
        return false;

    std::vector<double> descending_occupations = raw_natural_occupations;
    for (size_t i = 0; i < count; ++i)
    {
        if (!QC_SCF_Ensemble_Double_Is_Finite(
                ascending_fock_eigenvalues[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(
                fock_basis_density_diagonal[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(descending_occupations[i]) ||
            (i > 0 && ascending_fock_eigenvalues[i] <
                          ascending_fock_eigenvalues[i - 1]))
            return false;
    }
    std::sort(descending_occupations.begin(),
              descending_occupations.end(),
              [](double left, double right) { return left > right; });

    double compensation = 0.0;
    for (size_t i = 0; i < count; ++i)
    {
        const double term =
            ascending_fock_eigenvalues[i] *
            (descending_occupations[i] -
             fock_basis_density_diagonal[i]);
        const double corrected_term = term - compensation;
        const double updated_change = linear_change + corrected_term;
        compensation =
            (updated_change - linear_change) - corrected_term;
        linear_change = updated_change;
    }
    return QC_SCF_Ensemble_Double_Is_Finite(linear_change);
}

// Round to the nearest float unless that rounding raises coefficient*value.
// In that one case, move by exactly one float ULP toward the non-increasing
// side.  Applying this once per diagonal or symmetric off-diagonal AO pair
// keeps the float density symmetric and prevents representation rounding from
// reversing a valid double-precision fixed-F step.
#ifdef GPU_ARCH_NAME
#define QC_SCF_ENSEMBLE_HOST_DEVICE __host__ __device__
#else
#define QC_SCF_ENSEMBLE_HOST_DEVICE
#endif
static QC_SCF_ENSEMBLE_HOST_DEVICE inline float
QC_SCF_Round_Float_For_Nonincreasing_Linear_Objective(
    double value, double coefficient)
{
    float rounded = static_cast<float>(value);
    if (coefficient * (static_cast<double>(rounded) - value) > 0.0)
        rounded = nextafterf(rounded,
                             coefficient > 0.0 ? -FLT_MAX : FLT_MAX);
    return rounded;
}
#undef QC_SCF_ENSEMBLE_HOST_DEVICE

// Validate the spectrum-preserving orbital rotation independently from the
// capped-simplex repair that returns a tolerance-close float density to the
// nominal particle manifold.  The repair is certified by particle counts
// here and by the subsequent raw energy/commutator transaction.
static inline bool QC_SCF_Spectral_Orbital_Properties_Are_Valid(
    const std::vector<double>& original_particle_counts,
    const std::vector<double>& corrected_particle_counts,
    double spectrum_preserving_linear_change, double corrected_commutator,
    double energy_tolerance, double density_tolerance)
{
    if (original_particle_counts.empty() ||
        corrected_particle_counts.size() != original_particle_counts.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(
            spectrum_preserving_linear_change) ||
        !QC_SCF_Ensemble_Double_Is_Finite(corrected_commutator) ||
        !QC_SCF_Ensemble_Double_Is_Finite(energy_tolerance) ||
        !QC_SCF_Ensemble_Double_Is_Finite(density_tolerance) ||
        !(energy_tolerance > 0.0) || !(density_tolerance > 0.0) ||
        corrected_commutator < 0.0 ||
        spectrum_preserving_linear_change > energy_tolerance ||
        corrected_commutator > density_tolerance)
        return false;

    for (size_t i = 0; i < original_particle_counts.size(); ++i)
    {
        const double original = original_particle_counts[i];
        const double corrected = corrected_particle_counts[i];
        if (!QC_SCF_Ensemble_Double_Is_Finite(original) ||
            !QC_SCF_Ensemble_Double_Is_Finite(corrected) || original < 0.0 ||
            corrected < 0.0 ||
            std::fabs(corrected - original) >
                density_tolerance * std::max(1.0, original))
            return false;
    }
    return true;
}

static inline bool QC_SCF_Should_Start_Spectral_Orbital_Correction(
    double global_fw_gap, double active_fw_gap, double commutator,
    double energy_tolerance, double density_tolerance,
    bool repeated_verification)
{
    if (repeated_verification ||
        !QC_SCF_Ensemble_Double_Is_Finite(global_fw_gap) ||
        !QC_SCF_Ensemble_Double_Is_Finite(active_fw_gap) ||
        !QC_SCF_Ensemble_Double_Is_Finite(commutator) ||
        !QC_SCF_Ensemble_Double_Is_Finite(energy_tolerance) ||
        !QC_SCF_Ensemble_Double_Is_Finite(density_tolerance) ||
        global_fw_gap < 0.0 || active_fw_gap < 0.0 || commutator < 0.0 ||
        !(energy_tolerance > 0.0) || !(density_tolerance > 0.0) ||
        !(commutator > density_tolerance))
        return false;
    const double normalized_linear_gap =
        std::max(global_fw_gap, active_fw_gap) / energy_tolerance;
    const double normalized_commutator = commutator / density_tolerance;
    return QC_SCF_Ensemble_Double_Is_Finite(normalized_linear_gap) &&
           QC_SCF_Ensemble_Double_Is_Finite(normalized_commutator) &&
           normalized_commutator > normalized_linear_gap;
}

static inline bool QC_SCF_Ensemble_Global_FW_Gap(
    double directional_derivative, double& global_fw_gap)
{
    global_fw_gap = 0.0;
    if (!QC_SCF_Ensemble_Double_Is_Finite(directional_derivative))
        return false;
    // The Aufbau density minimizes the current linearized energy.  Therefore
    // only a negative directional derivative is a KKT violation; a positive
    // value is finite-precision oracle overshoot, not a gap of that size.
    global_fw_gap = std::max(0.0, -directional_derivative);
    return true;
}

static inline bool QC_SCF_Spectral_Orbital_Observation_Is_Acceptable(
    double origin_energy, double sample_energy, double sample_fraction,
    double repair_linear_component, double origin_commutator,
    double sample_commutator, double energy_tolerance,
    double density_tolerance, double energy_guard_multiplier = 2.0)
{
    // Float interpolation makes the represented spectral line piecewise
    // constant at small fractions, so its independently accumulated energy
    // can move by a few microhartree even while the generalized commutator
    // makes real progress.  Stay inside the existing energy guard and commit
    // only when at least one of the two physical residuals improves.  The
    // independent final KKT certificate still applies the original energy,
    // global-gap, active-gap, and commutator tolerances.
    // The origin may differ from the nominal particle manifold only by
    // admitted float representation noise.  Returning to that manifold is a
    // mandatory normal correction, not an orbital-rotation ascent.  Remove
    // its fixed-origin-F first-order contribution from the guard comparison;
    // the unchanged guard still measures the spectrum-preserving step,
    // float-line error, and nonlinear remainder.  The caller continues to
    // commit sample_energy itself, never this comparison-only value.
    const double guarded_sample_energy =
        sample_energy - sample_fraction * repair_linear_component;
    return QC_SCF_Ensemble_Double_Is_Finite(sample_fraction) &&
           sample_fraction > 0.0 && sample_fraction <= 1.0 &&
           QC_SCF_Ensemble_Double_Is_Finite(repair_linear_component) &&
           QC_SCF_Ensemble_Double_Is_Finite(guarded_sample_energy) &&
           QC_SCF_Ensemble_Energy_Within_Line_Guard(
               origin_energy, guarded_sample_energy, energy_tolerance,
               energy_guard_multiplier) &&
           QC_SCF_Ensemble_Double_Is_Finite(origin_commutator) &&
           QC_SCF_Ensemble_Double_Is_Finite(sample_commutator) &&
           QC_SCF_Ensemble_Double_Is_Finite(density_tolerance) &&
           density_tolerance > 0.0 && origin_commutator > density_tolerance &&
           sample_commutator >= 0.0 &&
           (sample_energy <= origin_energy ||
            sample_commutator < origin_commutator);
}

static constexpr int QC_SCF_ENSEMBLE_MAX_SPECTRAL_ORBITAL_BACKTRACKS = 20;

static inline bool QC_SCF_Next_Spectral_Orbital_Fraction(
    double current_fraction, int completed_evaluations, double& next_fraction)
{
    next_fraction = 0.0;
    if (!QC_SCF_Ensemble_Double_Is_Finite(current_fraction) ||
        !(current_fraction > 0.0) || current_fraction > 1.0 ||
        completed_evaluations < 0 ||
        completed_evaluations >=
            QC_SCF_ENSEMBLE_MAX_SPECTRAL_ORBITAL_BACKTRACKS)
        return false;

    // The KKT commutator tolerance is a final residual criterion, not a
    // minimum AO-density line step.  In particular, the RI response can turn
    // a sub-tolerance density displacement into a material commutator change.
    // Keep the transaction finite via the existing backtrack cap and let the
    // unchanged raw-energy/residual acceptance policy judge every candidate.
    const double candidate = 0.5 * current_fraction;
    if (!QC_SCF_Ensemble_Double_Is_Finite(candidate) || !(candidate > 0.0))
        return false;
    next_fraction = candidate;
    return true;
}

// Build the simplex weights represented by one corrective-line sample.  Both
// the trial and commit paths must use this exact operation: independently
// clamping tiny weights in only one path changes the reconstructed float
// density and invalidates a line-search sample after it is committed.
static inline bool QC_SCF_Build_Ensemble_Line_Weights(
    const std::vector<double>& origin_weights,
    const std::vector<double>& direction, double fraction,
    std::vector<double>& normalized_weights,
    double weight_tolerance = 1.0e-12,
    double normalization_tolerance = 1.0e-10)
{
    normalized_weights.assign(origin_weights.size(), 0.0);
    if (origin_weights.empty() || direction.size() != origin_weights.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(fraction) || fraction < 0.0 ||
        !QC_SCF_Ensemble_Double_Is_Finite(weight_tolerance) ||
        !(weight_tolerance > 0.0) ||
        !QC_SCF_Ensemble_Double_Is_Finite(normalization_tolerance) ||
        !(normalization_tolerance > 0.0))
        return false;

    double weight_sum = 0.0;
    for (size_t i = 0; i < origin_weights.size(); ++i)
    {
        if (!QC_SCF_Ensemble_Double_Is_Finite(origin_weights[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(direction[i]))
            return false;
        const double weight = origin_weights[i] + fraction * direction[i];
        if (!QC_SCF_Ensemble_Double_Is_Finite(weight) ||
            weight < -weight_tolerance)
            return false;
        normalized_weights[i] =
            weight <= weight_tolerance ? 0.0 : weight;
        weight_sum += normalized_weights[i];
    }
    if (!QC_SCF_Ensemble_Double_Is_Finite(weight_sum) ||
        !(weight_sum > 0.0) ||
        std::fabs(weight_sum - 1.0) > normalization_tolerance)
        return false;

    for (double& weight : normalized_weights) weight /= weight_sum;
    return true;
}

// A density-space proximity threshold cannot identify an active-set vertex:
// two distinct float matrices may be close in RMS while their Fock linear
// values differ by more than the requested KKT tolerance.  Only bytewise-
// equivalent float matrices produce an exactly zero accumulated RMS here.
static inline bool QC_SCF_Ensemble_Density_Is_Same_Vertex(
    double density_difference_rms)
{
    return QC_SCF_Ensemble_Double_Is_Finite(density_difference_rms) &&
           density_difference_rms == 0.0;
}

// A failed raw Aufbau map can expose a genuine zero-temperature ensemble
// minimum.  Along D(x) = D0 + x (D1-D0), Janak's theorem gives
// dE/dx = sum_sigma Tr[F_sigma(D(x)) (D1_sigma-D0_sigma)].  An interior
// ensemble is considered only when this derivative changes from negative at
// x=0 to positive at x=1; no electronic temperature or global smearing is
// introduced.
struct QC_SCF_Ensemble_Bracket
{
    double lower_fraction = 0.0;
    double upper_fraction = 1.0;
    double lower_derivative = 0.0;
    double upper_derivative = 0.0;
    double lower_energy = 0.0;
    double upper_energy = 0.0;
    double line_origin_energy = 0.0;
    double endpoint_energy_minimum = 0.0;
    double best_energy = 0.0;
    double best_fraction = 0.0;
    double best_derivative = 0.0;
    double best_stationarity_energy = 0.0;
    double best_stationarity_fraction = 0.0;
    double best_stationarity_derivative = 0.0;
    double direction_density_rms = 0.0;
    double previous_fraction = 0.0;
    int evaluations = 0;
    bool prior_interior_minimum = false;
};

enum QC_SCF_Ensemble_Root_Status
{
    QC_SCF_ENSEMBLE_ROOT_REJECTED = 0,
    QC_SCF_ENSEMBLE_ROOT_EVALUATE = 1,
    QC_SCF_ENSEMBLE_ROOT_STATIONARY = 2,
    QC_SCF_ENSEMBLE_ROOT_FAILED = 3,
    QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY = 4,
    QC_SCF_ENSEMBLE_ROOT_USE_BEST_STATIONARITY = 5,
};

struct QC_SCF_Ensemble_Root_Decision
{
    QC_SCF_Ensemble_Root_Status status = QC_SCF_ENSEMBLE_ROOT_REJECTED;
    double next_fraction = 0.0;
};

static constexpr int QC_SCF_ENSEMBLE_MAX_ROOT_EVALUATIONS = 64;

static inline double QC_SCF_Ensemble_Safeguarded_Secant(
    const QC_SCF_Ensemble_Bracket& bracket)
{
    const double lower = bracket.lower_fraction;
    const double upper = bracket.upper_fraction;
    const double width =
        upper - lower;
    const double denominator =
        bracket.upper_derivative - bracket.lower_derivative;
    double trial = lower + 0.5 * width;
    if (QC_SCF_Ensemble_Double_Is_Finite(denominator) && denominator > 0.0)
    {
        const double secant =
            lower - bracket.lower_derivative * width / denominator;
        if (QC_SCF_Ensemble_Double_Is_Finite(secant) && secant > lower &&
            secant < upper)
            trial = secant;
    }

    // Do not impose a fixed fractional margin: roots near an endpoint are
    // common for a well-preconditioned corrective line, and clipping 0.003 to
    // 0.1 then 0.01 wastes two raw Fock builds.  Fall back only when the secant
    // is not a strict interior double.  nextafter handles a rounded endpoint;
    // if no representable interior point exists, Observe_Ensemble_Trial's
    // density-resolution check terminates the line.
    if (!(trial > lower) || !(trial < upper)) trial = lower + 0.5 * width;
    if (!(trial > lower)) trial = std::nextafter(lower, upper);
    if (!(trial < upper)) trial = std::nextafter(upper, lower);
    return trial;
}

static inline QC_SCF_Ensemble_Root_Decision QC_SCF_Start_Ensemble_Root(
    double lower_derivative, double lower_energy, double upper_derivative,
    double upper_energy, double direction_density_rms,
    double derivative_tolerance, QC_SCF_Ensemble_Bracket& bracket,
    bool prior_interior_minimum = false, double upper_fraction = 1.0)
{
    bracket = {};
    if (!QC_SCF_Ensemble_Double_Is_Finite(lower_derivative) ||
        !QC_SCF_Ensemble_Double_Is_Finite(upper_derivative) ||
        !QC_SCF_Ensemble_Double_Is_Finite(lower_energy) ||
        !QC_SCF_Ensemble_Double_Is_Finite(upper_energy) ||
        !QC_SCF_Ensemble_Double_Is_Finite(direction_density_rms) ||
        !QC_SCF_Ensemble_Double_Is_Finite(derivative_tolerance) ||
        !QC_SCF_Ensemble_Double_Is_Finite(upper_fraction) ||
        !(direction_density_rms > 0.0) ||
        !(derivative_tolerance > 0.0) ||
        !(upper_fraction > 0.0) ||
        !(lower_derivative < -derivative_tolerance) ||
        !(upper_derivative > derivative_tolerance))
        return {};

    bracket.upper_fraction = upper_fraction;
    bracket.lower_derivative = lower_derivative;
    bracket.upper_derivative = upper_derivative;
    bracket.lower_energy = lower_energy;
    bracket.upper_energy = upper_energy;
    bracket.line_origin_energy = lower_energy;
    bracket.endpoint_energy_minimum =
        std::min(lower_energy, upper_energy);
    if (lower_energy <= upper_energy)
    {
        bracket.best_energy = lower_energy;
        bracket.best_fraction = 0.0;
        bracket.best_derivative = lower_derivative;
    }
    else
    {
        bracket.best_energy = upper_energy;
        bracket.best_fraction = upper_fraction;
        bracket.best_derivative = upper_derivative;
    }
    bracket.direction_density_rms = direction_density_rms;
    bracket.previous_fraction = upper_fraction;
    bracket.best_stationarity_energy = upper_energy;
    bracket.best_stationarity_fraction = upper_fraction;
    bracket.best_stationarity_derivative = upper_derivative;
    bracket.prior_interior_minimum = prior_interior_minimum;
    return {QC_SCF_ENSEMBLE_ROOT_EVALUATE,
            QC_SCF_Ensemble_Safeguarded_Secant(bracket)};
}

static inline QC_SCF_Ensemble_Root_Decision QC_SCF_Observe_Ensemble_Trial(
    QC_SCF_Ensemble_Bracket& bracket, double fraction, double derivative,
    double energy, double derivative_tolerance, double energy_tolerance,
    double density_tolerance)
{
    if (!QC_SCF_Ensemble_Double_Is_Finite(fraction) ||
        !QC_SCF_Ensemble_Double_Is_Finite(derivative) ||
        !QC_SCF_Ensemble_Double_Is_Finite(energy) ||
        !QC_SCF_Ensemble_Double_Is_Finite(derivative_tolerance) ||
        !QC_SCF_Ensemble_Double_Is_Finite(energy_tolerance) ||
        !QC_SCF_Ensemble_Double_Is_Finite(density_tolerance) ||
        !(derivative_tolerance > 0.0) || !(energy_tolerance > 0.0) ||
        !(density_tolerance > 0.0))
        return {QC_SCF_ENSEMBLE_ROOT_FAILED, 0.0};

    if (!(fraction > bracket.lower_fraction) ||
        !(fraction < bracket.upper_fraction))
    {
        const double scale =
            std::max(1.0, std::max(std::fabs(fraction),
                                   std::max(std::fabs(bracket.lower_fraction),
                                            std::fabs(bracket.upper_fraction))));
        const double resolution =
            16.0 * std::numeric_limits<double>::epsilon() * scale;
        const bool exhausted_lower =
            std::fabs(fraction - bracket.lower_fraction) <= resolution;
        const bool exhausted_upper =
            std::fabs(fraction - bracket.upper_fraction) <= resolution;
        if ((exhausted_lower || exhausted_upper) &&
            bracket.best_fraction > 0.0)
            return {QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY,
                    bracket.best_fraction};
        return {QC_SCF_ENSEMBLE_ROOT_FAILED, 0.0};
    }

    ++bracket.evaluations;
    if (bracket.evaluations > QC_SCF_ENSEMBLE_MAX_ROOT_EVALUATIONS)
    {
        if (bracket.best_fraction > 0.0 &&
            bracket.best_energy < bracket.line_origin_energy)
            return {QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY,
                    bracket.best_fraction};
        return {QC_SCF_ENSEMBLE_ROOT_FAILED, 0.0};
    }

    // A derivative bracket g(0)<0<g(gamma_max) is the primary proof of an
    // interior line minimum.  Requiring that minimum to beat an endpoint by
    // a full SCF tolerance is mathematically wrong near convergence, where
    // the remaining decrease must tend to zero.  Energy is retained as an
    // independent consistency guard: two independently accumulated endpoint
    // energies can each carry one tolerance of numerical uncertainty, so a
    // root may be at most their combined uncertainty above the better end.
    if (energy < bracket.best_energy)
    {
        bracket.best_energy = energy;
        bracket.best_fraction = fraction;
        bracket.best_derivative = derivative;
    }
    if (fraction > 0.0 &&
        std::fabs(derivative) <
            std::fabs(bracket.best_stationarity_derivative))
    {
        bracket.best_stationarity_energy = energy;
        bracket.best_stationarity_fraction = fraction;
        bracket.best_stationarity_derivative = derivative;
    }
    const bool has_interior_energy_minimum =
        energy <= bracket.endpoint_energy_minimum + 2.0 * energy_tolerance;
    const double density_step =
        std::fabs(fraction - bracket.previous_fraction) *
        bracket.direction_density_rms;
    const double bracket_density_width =
        (bracket.upper_fraction - bracket.lower_fraction) *
        bracket.direction_density_rms;
    const bool derivative_stationary =
        std::fabs(derivative) <= derivative_tolerance;
    const bool density_stationary =
        density_step <= density_tolerance ||
        bracket_density_width <= density_tolerance;

    if (derivative_stationary && has_interior_energy_minimum)
    {
        // The smooth directional derivative locates the variational root,
        // while raw energies are observations of the actually represented
        // float densities.  If another sampled density on this same line has
        // a strictly lower raw energy, verify and commit that exact sample.
        if (bracket.best_fraction > 0.0 && bracket.best_energy < energy)
            return {QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY,
                    bracket.best_fraction};
        return {QC_SCF_ENSEMBLE_ROOT_STATIONARY, fraction};
    }
    // The caller commits this sampled density and immediately performs a raw
    // F[P]/E[P] build before the global KKT check.  Repeating the trial once
    // here as well would perform two identical internal confirmations; only
    // the final global certificate requires its own additional same-P repeat.
    // At float density resolution, two adjacent representable active-set
    // mixtures can straddle the derivative root while their independently
    // accumulated energies differ from the line origin by a few SCF energy
    // tolerances.  The sign-changing derivative bracket remains the primary
    // variational certificate.  Select the actually sampled point with the
    // smallest |dE/dgamma|, then require the caller to repeat that exact P
    // before committing it.  The four-tolerance cap is confined to this
    // resolution-exhausted branch and still rejects a macroscopic rise.
    const bool resolution_energy_consistent =
        bracket.best_stationarity_energy <=
        bracket.line_origin_energy + 4.0 * energy_tolerance;
    if (density_stationary &&
        bracket.best_stationarity_fraction > 0.0 &&
        resolution_energy_consistent)
        return {QC_SCF_ENSEMBLE_ROOT_USE_BEST_STATIONARITY,
                bracket.best_stationarity_fraction};
    if (derivative_stationary && !has_interior_energy_minimum &&
        bracket.best_fraction > 0.0 &&
        bracket.best_fraction <= bracket.upper_fraction &&
        bracket.best_energy < energy)
        return {QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY,
                bracket.best_fraction};
    if (!derivative_stationary && density_stationary &&
        bracket.best_fraction > 0.0 &&
        bracket.best_fraction <= bracket.upper_fraction &&
        bracket.best_energy <=
            bracket.endpoint_energy_minimum + 2.0 * energy_tolerance)
        return {QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY,
                bracket.best_fraction};
    if (density_stationary && bracket.best_fraction > 0.0 &&
        bracket.best_energy < bracket.line_origin_energy)
        return {QC_SCF_ENSEMBLE_ROOT_USE_BEST_ENERGY,
                bracket.best_fraction};
    // Once the bracket is narrower than the representable density resolution,
    // another fraction cannot provide new SCF information.  If neither a
    // stationary, energy-consistent sample nor an actually descending sample
    // exists, fail this line immediately so the caller can reject/restart a
    // quasi-Newton proposal.  Iterating to the numerical fraction limit only
    // repeats the same float density and can consume the entire SCF budget.
    if (density_stationary)
        return {QC_SCF_ENSEMBLE_ROOT_FAILED, 0.0};
    if (derivative < 0.0)
    {
        bracket.lower_fraction = fraction;
        bracket.lower_derivative = derivative;
        bracket.lower_energy = energy;
    }
    else
    {
        bracket.upper_fraction = fraction;
        bracket.upper_derivative = derivative;
        bracket.upper_energy = energy;
    }
    bracket.previous_fraction = fraction;
    bracket.endpoint_energy_minimum =
        std::min(bracket.endpoint_energy_minimum, energy);

    return {QC_SCF_ENSEMBLE_ROOT_EVALUATE,
            QC_SCF_Ensemble_Safeguarded_Secant(bracket)};
}

static inline bool QC_SCF_Ensemble_Active_Set_FW_Gap(
    const std::vector<double>& linear_values,
    const std::vector<double>& weights, double& gap,
    double weight_tolerance = 1.0e-12)
{
    gap = std::numeric_limits<double>::max();
    if (linear_values.empty() || linear_values.size() != weights.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(weight_tolerance) ||
        !(weight_tolerance > 0.0))
        return false;

    double weight_sum = 0.0;
    double weighted_linear_value = 0.0;
    double minimum_linear_value = std::numeric_limits<double>::max();
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (!QC_SCF_Ensemble_Double_Is_Finite(linear_values[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(weights[i]) ||
            weights[i] < -weight_tolerance)
            return false;
        const double weight = std::max(0.0, weights[i]);
        weight_sum += weight;
        weighted_linear_value += weight * linear_values[i];
        minimum_linear_value =
            std::min(minimum_linear_value, linear_values[i]);
    }
    if (!QC_SCF_Ensemble_Double_Is_Finite(weight_sum) ||
        !(weight_sum > 0.0))
        return false;

    gap = weighted_linear_value / weight_sum - minimum_linear_value;
    if (gap < 0.0 && gap >= -64.0 * weight_tolerance) gap = 0.0;
    return QC_SCF_Ensemble_Double_Is_Finite(gap) && gap >= 0.0;
}

// Fully correcting every intermediate active hull to the final SCF tolerance
// defeats the purpose of a global Frank-Wolfe oracle: most of that work is
// invalidated as soon as the next Aufbau vertex is admitted.  Use the standard
// inexact-oracle forcing rule for that scheduling decision only.  The final
// KKT test remains independent and still requires both gaps to satisfy the
// original tolerance.
static inline bool QC_SCF_Ensemble_Active_Correction_Target(
    double global_fw_gap, double final_tolerance, double& target,
    double forcing_fraction = 0.25)
{
    target = std::numeric_limits<double>::max();
    if (!QC_SCF_Ensemble_Double_Is_Finite(global_fw_gap) ||
        !QC_SCF_Ensemble_Double_Is_Finite(final_tolerance) ||
        !QC_SCF_Ensemble_Double_Is_Finite(forcing_fraction) ||
        global_fw_gap < 0.0 || !(final_tolerance > 0.0) ||
        !(forcing_fraction > 0.0) || !(forcing_fraction < 1.0))
        return false;
    target = std::max(final_tolerance, forcing_fraction * global_fw_gap);
    return QC_SCF_Ensemble_Double_Is_Finite(target);
}

static inline void QC_SCF_Initialize_Corrective_Inverse_Hessian(
    size_t active_atom_count, std::vector<double>& inverse_hessian)
{
    const size_t dimension = active_atom_count > 0 ? active_atom_count - 1 : 0;
    inverse_hessian.assign(dimension * dimension, 0.0);
    for (size_t i = 0; i < dimension; ++i)
        inverse_hessian[i * dimension + i] = 1.0;
}

// Apply one inverse-BFGS secant update in independent simplex coordinates and
// form a feasible quasi-Newton direction.  The caller keeps the active hull
// fixed between the reference and current raw Fock builds.  If curvature or a
// boundary constraint makes the update unsafe, returning false deliberately
// falls back to the tangent-cone projected gradient below.
static inline bool QC_SCF_Build_Quasi_Newton_Corrective_Direction(
    const std::vector<double>& linear_values,
    const std::vector<double>& weights,
    const std::vector<double>& reference_linear_values,
    const std::vector<double>& reference_weights,
    std::vector<double>& inverse_hessian, std::vector<double>& direction,
    double weight_tolerance = 1.0e-12)
{
    direction.assign(weights.size(), 0.0);
    if (weights.size() < 2 || linear_values.size() != weights.size() ||
        reference_linear_values.size() != weights.size() ||
        reference_weights.size() != weights.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(weight_tolerance) ||
        !(weight_tolerance > 0.0))
        return false;

    const size_t dimension = weights.size() - 1;
    if (inverse_hessian.size() != dimension * dimension) return false;
    for (size_t i = 0; i < weights.size(); ++i)
        if (!QC_SCF_Ensemble_Double_Is_Finite(linear_values[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(weights[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(reference_linear_values[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(reference_weights[i]))
            return false;
    std::vector<double> updated_inverse_hessian = inverse_hessian;
    std::vector<double> step(dimension, 0.0);
    std::vector<double> gradient_change(dimension, 0.0);
    std::vector<double> reduced_gradient(dimension, 0.0);
    double step_norm_squared = 0.0;
    double change_norm_squared = 0.0;
    double curvature = 0.0;
    for (size_t i = 0; i < dimension; ++i)
    {
        const size_t atom = i + 1;
        step[i] = weights[atom] - reference_weights[atom];
        gradient_change[i] =
            (linear_values[atom] - linear_values[0]) -
            (reference_linear_values[atom] -
             reference_linear_values[0]);
        reduced_gradient[i] = linear_values[atom] - linear_values[0];
        step_norm_squared += step[i] * step[i];
        change_norm_squared += gradient_change[i] * gradient_change[i];
        curvature += step[i] * gradient_change[i];
    }
    const double curvature_floor =
        1.0e-10 * std::sqrt(step_norm_squared * change_norm_squared);
    if (!QC_SCF_Ensemble_Double_Is_Finite(curvature) ||
        !(step_norm_squared > 0.0) || !(change_norm_squared > 0.0) ||
        !(curvature > std::max(1.0e-16, curvature_floor)))
        return false;

    std::vector<double> hessian_times_change(dimension, 0.0);
    for (size_t i = 0; i < dimension; ++i)
        for (size_t j = 0; j < dimension; ++j)
            hessian_times_change[i] +=
                updated_inverse_hessian[i * dimension + j] *
                gradient_change[j];
    double change_hessian_change = 0.0;
    for (size_t i = 0; i < dimension; ++i)
        change_hessian_change +=
            gradient_change[i] * hessian_times_change[i];
    if (!QC_SCF_Ensemble_Double_Is_Finite(change_hessian_change) ||
        !(change_hessian_change >= 0.0))
        return false;

    const double rank_one_scale =
        (curvature + change_hessian_change) / (curvature * curvature);
    for (size_t i = 0; i < dimension; ++i)
        for (size_t j = 0; j < dimension; ++j)
        {
            updated_inverse_hessian[i * dimension + j] +=
                rank_one_scale * step[i] * step[j] -
                (hessian_times_change[i] * step[j] +
                 step[i] * hessian_times_change[j]) /
                    curvature;
            if (!QC_SCF_Ensemble_Double_Is_Finite(
                    updated_inverse_hessian[i * dimension + j]))
                return false;
        }

    double independent_sum = 0.0;
    for (size_t i = 0; i < dimension; ++i)
    {
        double value = 0.0;
        for (size_t j = 0; j < dimension; ++j)
            value -= updated_inverse_hessian[i * dimension + j] *
                     reduced_gradient[j];
        direction[i + 1] = value;
        independent_sum += value;
    }
    direction[0] = -independent_sum;

    double maximum_step = std::numeric_limits<double>::max();
    double derivative = 0.0;
    double direction_sum = 0.0;
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (!QC_SCF_Ensemble_Double_Is_Finite(direction[i]) ||
            weights[i] < -weight_tolerance ||
            (weights[i] <= weight_tolerance && direction[i] < 0.0))
            return false;
        derivative += linear_values[i] * direction[i];
        direction_sum += direction[i];
        if (direction[i] < 0.0)
            maximum_step =
                std::min(maximum_step, weights[i] / -direction[i]);
    }
    if (!QC_SCF_Ensemble_Double_Is_Finite(derivative) ||
        !(derivative < -1.0e-14) ||
        std::fabs(direction_sum) > 1.0e-10 ||
        !QC_SCF_Ensemble_Double_Is_Finite(maximum_step) ||
        !(maximum_step > 0.0))
        return false;

    for (double& value : direction) value *= maximum_step;
    inverse_hessian.swap(updated_inverse_hessian);
    return true;
}

// Project -grad(E) onto the tangent cone of the current active simplex.  All
// positive-weight atoms participate; zero-weight atoms enter exactly when
// their linear value lies below the active mean.  The result is scaled so a
// unit line step reaches (but never crosses) the first simplex face.  Thus one
// corrective line simultaneously reoptimizes every active weight and can drop
// obsolete atoms, instead of transferring a single pair of weights per raw
// Fock build.
static inline bool QC_SCF_Build_Fully_Corrective_Direction(
    const std::vector<double>& linear_values,
    const std::vector<double>& weights, std::vector<double>& direction,
    double weight_tolerance = 1.0e-12)
{
    direction.assign(weights.size(), 0.0);
    if (linear_values.empty() || linear_values.size() != weights.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(weight_tolerance) ||
        !(weight_tolerance > 0.0))
        return false;

    std::vector<unsigned char> eligible(weights.size(), 0);
    int eligible_count = 0;
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (!QC_SCF_Ensemble_Double_Is_Finite(linear_values[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(weights[i]) ||
            weights[i] < -weight_tolerance)
            return false;
        if (weights[i] > weight_tolerance)
        {
            eligible[i] = 1;
            ++eligible_count;
        }
    }
    if (eligible_count == 0) return false;

    double eligible_sum = 0.0;
    std::vector<size_t> zero_indices;
    zero_indices.reserve(weights.size());
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (eligible[i])
            eligible_sum += linear_values[i];
        else
            zero_indices.push_back(i);
    }
    std::sort(zero_indices.begin(), zero_indices.end(),
              [&](size_t lhs, size_t rhs)
              { return linear_values[lhs] < linear_values[rhs]; });

    double mean = eligible_sum / static_cast<double>(eligible_count);
    for (const size_t index : zero_indices)
    {
        // Add boundary atoms one at a time in increasing-gradient order.  A
        // bulk add against the old mean can make a later zero-weight atom's
        // projected direction negative after the mean moves, which is outside
        // the simplex tangent cone.
        if (linear_values[index] < mean)
        {
            eligible[index] = 1;
            eligible_sum += linear_values[index];
            ++eligible_count;
            mean = eligible_sum / static_cast<double>(eligible_count);
        }
        else
            break;
    }

    double direction_sum = 0.0;
    int correction_index = -1;
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (!eligible[i]) continue;
        direction[i] = mean - linear_values[i];
        direction_sum += direction[i];
        if (correction_index < 0 && weights[i] > weight_tolerance)
            correction_index = static_cast<int>(i);
    }
    if (correction_index < 0) return false;
    direction[correction_index] -= direction_sum;

    double maximum_step = std::numeric_limits<double>::max();
    double squared_norm = 0.0;
    for (size_t i = 0; i < weights.size(); ++i)
    {
        squared_norm += direction[i] * direction[i];
        if (direction[i] < 0.0)
            maximum_step =
                std::min(maximum_step, weights[i] / -direction[i]);
    }
    if (!QC_SCF_Ensemble_Double_Is_Finite(squared_norm) ||
        !(squared_norm > 0.0) ||
        !QC_SCF_Ensemble_Double_Is_Finite(maximum_step) ||
        !(maximum_step > 0.0))
        return false;

    for (double& value : direction) value *= maximum_step;
    double scaled_sum = 0.0;
    for (size_t i = 0; i < weights.size(); ++i)
    {
        scaled_sum += direction[i];
        if (weights[i] + direction[i] < -weight_tolerance) return false;
    }
    if (std::fabs(scaled_sum) > 1.0e-10) return false;
    return true;
}

// A strict pairwise Frank-Wolfe face step used only to globalize a rejected
// multi-vertex corrective line.  Move weight from the largest-linear-value
// positive vertex to the global minimum-linear-value vertex.  Scaling by the
// away weight makes a unit step exactly feasible and exposes a lower-
// dimensional face when the endpoint is accepted.
static inline bool QC_SCF_Build_Pairwise_Corrective_Direction(
    const std::vector<double>& linear_values,
    const std::vector<double>& weights, std::vector<double>& direction,
    double weight_tolerance = 1.0e-12)
{
    direction.assign(weights.size(), 0.0);
    if (weights.size() < 2 || linear_values.size() != weights.size() ||
        !QC_SCF_Ensemble_Double_Is_Finite(weight_tolerance) ||
        !(weight_tolerance > 0.0))
        return false;

    size_t toward = 0;
    size_t away = weights.size();
    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (!QC_SCF_Ensemble_Double_Is_Finite(linear_values[i]) ||
            !QC_SCF_Ensemble_Double_Is_Finite(weights[i]) ||
            weights[i] < -weight_tolerance)
            return false;
        if (linear_values[i] < linear_values[toward]) toward = i;
        if (weights[i] > weight_tolerance &&
            (away == weights.size() ||
             linear_values[i] > linear_values[away]))
            away = i;
    }
    if (away == weights.size() || toward == away ||
        !(linear_values[toward] < linear_values[away]))
        return false;

    const double transfer = weights[away];
    if (!QC_SCF_Ensemble_Double_Is_Finite(transfer) || !(transfer > 0.0))
        return false;
    direction[away] = -transfer;
    direction[toward] = transfer;
    return true;
}
