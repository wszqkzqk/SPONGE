#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

static inline bool QC_SCF_Convergence_Double_Is_Finite(double value)
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

// Level shifting is a homotopy, not a binary stabilization switch.  Each
// positive stage must converge before moving to the next one, and the last
// stage is exactly zero rather than the result of repeated floating-point
// multiplication.  The two reduced positive stages remain large enough to
// stabilize the orbital response.  Once the last positive stage is stable,
// its density is handed directly to one unshifted, non-extrapolated physical
// map.  Requiring an intervening zero-shift accelerated solve can destroy that
// stable density when the ordinary Aufbau map is unstable near a degeneracy.
static constexpr int QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT = 3;
static constexpr double QC_SCF_LEVEL_SHIFT_STAGE_RATIO = 0.25;

static inline double QC_SCF_Level_Shift_For_Stage(double configured_shift,
                                                  int stage)
{
    if (!QC_SCF_Convergence_Double_Is_Finite(configured_shift) ||
        !(configured_shift > 0.0) || stage < 0 ||
        stage >= QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT)
        return 0.0;

    double shift = configured_shift;
    for (int i = 0; i < stage; ++i)
        shift *= QC_SCF_LEVEL_SHIFT_STAGE_RATIO;
    return shift;
}

// State for the continuation plus the two-stage SCF stopping criterion.
// Accelerated iterations (DIIS and/or level shifting) may advance the
// homotopy or identify an unshifted candidate, but only one unshifted,
// non-extrapolated F[P] -> P' map is allowed to confirm a physical fixed
// point.
struct QC_SCF_Convergence_State
{
    int accelerated_streak = 0;
    int physical_streak = 0;
    int level_shift_stage = 0;
    bool verifying_physical_fixed_point = false;
};

struct QC_SCF_Convergence_Decision
{
    bool converged = false;
    bool reset_diis_history = false;
    bool advanced_level_shift = false;
};

static inline QC_SCF_Convergence_Decision QC_SCF_Observe_Iteration(
    QC_SCF_Convergence_State& state, bool physical_iteration,
    bool has_previous_energy, double delta_energy, double density_residual,
    double energy_tolerance, double density_tolerance,
    double configured_level_shift)
{
    QC_SCF_Convergence_Decision decision;
    const bool density_stable = density_residual < density_tolerance;
    const bool energy_stable =
        has_previous_energy && std::fabs(delta_energy) < energy_tolerance;

    if (!physical_iteration)
    {
        state.physical_streak = 0;
        if (!(density_stable && energy_stable))
        {
            state.accelerated_streak = 0;
            return decision;
        }

        const double active_shift = QC_SCF_Level_Shift_For_Stage(
            configured_level_shift, state.level_shift_stage);
        if (active_shift > 0.0)
        {
            // Preserve P_{n+1} as the initial density for the next homotopy
            // stage, but never mix Fock/density history across two different
            // level-shift maps.  A residual rebound at the lower stage merely
            // delays its next transition; it cannot latch the old shift back
            // on permanently.
            ++state.level_shift_stage;
            state.accelerated_streak = 0;
            if (state.level_shift_stage >=
                QC_SCF_LEVEL_SHIFT_POSITIVE_STAGE_COUNT)
            {
                // The converged final shifted density is already the best
                // candidate supplied by the homotopy.  Verify it with a raw
                // zero-shift F[P] map before allowing zero-shift DIIS to move
                // it away from the stable branch.
                state.verifying_physical_fixed_point = true;
            }
            decision.reset_diis_history = true;
            decision.advanced_level_shift = true;
            return decision;
        }

        ++state.accelerated_streak;

        // An accelerated map is only a candidate generator.  Its Fock matrix
        // can be extrapolated, so a small step here is not proof of a physical
        // SCF fixed point.  This path supplies the candidate when level
        // shifting is disabled; a completed positive-shift continuation goes
        // directly to the same raw physical verification above.
        if (state.accelerated_streak >= 1)
        {
            state.verifying_physical_fixed_point = true;
        }
        return decision;
    }

    state.accelerated_streak = 0;
    // This residual is ||P[F(P)] - P|| for an unshifted, non-extrapolated
    // Fock matrix built from P.  It is therefore the physical fixed-point
    // oracle itself; a second raw SCF map is neither a stronger test nor
    // necessary for energy alignment.  In the float32 density path it can
    // amplify roundoff along a weakly contractive direction and reject an
    // already verified fixed point.  The caller retains P (rather than P') on
    // termination, so the published energy, Fock matrix, and density remain
    // from the same physical iteration.
    if (density_stable)
    {
        state.physical_streak = 1;
        decision.converged = true;
    }
    else
    {
        // Do not immediately reuse a history which produced a false fixed
        // point candidate.  The next accelerated phase must build a fresh
        // coupled DIIS subspace around the physical-map result.
        state.physical_streak = 0;
        state.verifying_physical_fixed_point = false;
        decision.reset_diis_history = true;
    }
    return decision;
}

static inline double QC_SCF_Level_Shift_For_Map(
    double configured_shift, const QC_SCF_Convergence_State& state,
    bool physical_iteration)
{
    if (physical_iteration) return 0.0;
    return QC_SCF_Level_Shift_For_Stage(configured_shift,
                                        state.level_shift_stage);
}

static inline double QC_SCF_Level_Shift_Density_Factor(bool unrestricted)
{
    // P has occupation two in a restricted calculation and occupation one in
    // each unrestricted spin channel.  S - S P S / occupation is the AO
    // representation of the virtual-space projector.
    return unrestricted ? 1.0 : 0.5;
}

static inline bool QC_SCF_Should_Use_Incremental_Fock(
    int iter, bool force_full_rebuild)
{
    // The first two SCF maps must be built from the full density because the
    // SAP-to-SCF transition is too large for useful delta-density screening.
    // A physical fixed-point check also needs a fresh F[P], independent of
    // any rounding or screening error accumulated by incremental updates.
    return iter >= 2 && !force_full_rebuild;
}
