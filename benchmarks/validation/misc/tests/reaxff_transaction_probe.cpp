#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include "manybody/reaxff/reaxff_evaluation.h"

namespace
{

bool Equal_Vector(const VECTOR& lhs, const VECTOR& rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

bool Equal_Matrix(const LTMatrix3& lhs, const LTMatrix3& rhs)
{
    return lhs.a11 == rhs.a11 && lhs.a21 == rhs.a21 &&
           lhs.a22 == rhs.a22 && lhs.a31 == rhs.a31 &&
           lhs.a32 == rhs.a32 && lhs.a33 == rhs.a33;
}

template <typename T, int N, typename Equal>
bool Equal_Array(const T (&lhs)[N], const T (&rhs)[N], Equal equal)
{
    for (int i = 0; i < N; i++)
        if (!equal(lhs[i], rhs[i])) return false;
    return true;
}

template <int N>
bool Equal_Floats(const float (&lhs)[N], const float (&rhs)[N])
{
    for (int i = 0; i < N; i++)
        if (lhs[i] != rhs[i]) return false;
    return true;
}

}  // namespace

int main()
{
    constexpr int atoms = 3;
    constexpr int history_capacity = 5;
    constexpr int history_count = 2;
    const int local_to_global[atoms] = {2, 0, 1};

    int checked_coordinate_count = -1;
    if (ReaxFFEvaluation::Checked_Coordinate_Count(
            atoms, 2, sizeof(LTMatrix3), &checked_coordinate_count) !=
            ReaxFFEvaluation::COORDINATE_COUNT_OK ||
        checked_coordinate_count != atoms + 2)
    {
        std::fprintf(stderr,
                     "owned+ghost coordinate extent was not checked exactly\n");
        return EXIT_FAILURE;
    }
    if (ReaxFFEvaluation::Checked_Coordinate_Count(
            std::numeric_limits<int>::max(), 1, sizeof(LTMatrix3),
            &checked_coordinate_count) !=
        ReaxFFEvaluation::COORDINATE_COUNT_INT_OVERFLOW)
    {
        std::fprintf(stderr, "coordinate int overflow was not rejected\n");
        return EXIT_FAILURE;
    }
    const std::size_t overflowing_record_size =
        std::numeric_limits<std::size_t>::max() / 2U + 1U;
    if (ReaxFFEvaluation::Checked_Coordinate_Count(
            1, 1, overflowing_record_size, &checked_coordinate_count) !=
        ReaxFFEvaluation::COORDINATE_COUNT_SIZE_OVERFLOW)
    {
        std::fprintf(stderr, "coordinate size_t overflow was not rejected\n");
        return EXIT_FAILURE;
    }

    const VECTOR angle_u = {1.0f, 0.0f, 0.0f};
    const VECTOR angle_v = {1.0f, 1.0e-4f, 0.0f};
    REAXFF_ANGLE_GEOMETRY tiny_angle;
    if (!Compute_ReaxFF_Angle_Geometry(angle_u, angle_v, &tiny_angle))
    {
        std::fprintf(stderr, "tiny nonzero ReaxFF angle was rejected\n");
        return EXIT_FAILURE;
    }
    const float tiny_hb_energy_factor =
        tiny_angle.half_sine_squared * tiny_angle.half_sine_squared;
    const float tiny_hb_angle_derivative =
        tiny_angle.half_sine_squared * tiny_angle.sine;
    const float tiny_hb_force_component = static_cast<float>(
        -static_cast<double>(tiny_hb_angle_derivative) *
        tiny_angle.dtheta_dv[1]);
    const double expected_half_sine_squared =
        0.5 * (1.0 - 1.0 / std::sqrt(1.0 + 1.0e-8));
    if (!(tiny_angle.half_sine_squared > 0.0f) ||
        !(tiny_hb_energy_factor > 0.0f) ||
        !(tiny_hb_angle_derivative > 0.0f) ||
        tiny_hb_force_component == 0.0f ||
        !std::isfinite(tiny_hb_force_component) ||
        std::fabs(static_cast<double>(tiny_angle.half_sine_squared) -
                  expected_half_sine_squared) >
            1.0e-6 * expected_half_sine_squared)
    {
        std::fprintf(stderr,
                     "tiny ReaxFF hydrogen-bond angle lost its half-angle "
                     "energy/force factor: q=%g E-factor=%g d-factor=%g "
                     "force=%g expected-q=%.17g\n",
                     tiny_angle.half_sine_squared, tiny_hb_energy_factor,
                     tiny_hb_angle_derivative, tiny_hb_force_component,
                     expected_half_sine_squared);
        return EXIT_FAILURE;
    }

    VECTOR committed_frc[atoms] = {
        {10.0f, 11.0f, 12.0f},
        {20.0f, 21.0f, 22.0f},
        {30.0f, 31.0f, 32.0f}};
    float committed_energy[atoms] = {100.0f, 200.0f, 300.0f};
    LTMatrix3 committed_virial[atoms] = {
        {1, 2, 3, 4, 5, 6},
        {7, 8, 9, 10, 11, 12},
        {13, 14, 15, 16, 17, 18}};
    float committed_local_charge[atoms] = {41.0f, 42.0f, 43.0f};
    float committed_global_charge[atoms] = {51.0f, 52.0f, 53.0f};
    float committed_s_history[history_capacity * atoms] = {};
    float committed_t_history[history_capacity * atoms] = {};
    for (int frame = 0; frame < history_capacity; frame++)
    {
        for (int atom = 0; atom < atoms; atom++)
        {
            committed_s_history[frame * atoms + atom] =
                1000.0f + 10.0f * frame + atom;
            committed_t_history[frame * atoms + atom] =
                2000.0f + 10.0f * frame + atom;
        }
    }

    const VECTOR original_frc[atoms] = {committed_frc[0], committed_frc[1],
                                        committed_frc[2]};
    const float original_energy[atoms] = {committed_energy[0],
                                          committed_energy[1],
                                          committed_energy[2]};
    const LTMatrix3 original_virial[atoms] = {
        committed_virial[0], committed_virial[1], committed_virial[2]};
    const float original_local_charge[atoms] = {
        committed_local_charge[0], committed_local_charge[1],
        committed_local_charge[2]};
    const float original_global_charge[atoms] = {
        committed_global_charge[0], committed_global_charge[1],
        committed_global_charge[2]};
    float original_s_history[history_capacity * atoms] = {};
    float original_t_history[history_capacity * atoms] = {};
    for (int i = 0; i < history_capacity * atoms; i++)
    {
        original_s_history[i] = committed_s_history[i];
        original_t_history[i] = committed_t_history[i];
    }

    VECTOR staged_frc[atoms] = {
        {0.1f, 0.2f, 0.3f},
        {0.4f, 0.5f, 0.6f},
        {0.7f, 0.8f, 0.9f}};
    float staged_energy[atoms] = {1.0f,
                                  std::numeric_limits<float>::quiet_NaN(), 3.0f};
    LTMatrix3 staged_virial[atoms] = {
        {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f},
        {0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f},
        {1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f}};
    float staged_charge[atoms] = {-1.0f, 2.0f, -3.0f};
    float candidate_s[atoms] = {5.0f, 6.0f, 7.0f};
    float candidate_t[atoms] = {8.0f, 9.0f, 10.0f};
    int error[REAXFF_GEOMETRY_ERROR_SIZE] = {};

    ReaxFFEvaluation::Validate_Staging_Kernel(
        atoms, staged_frc, staged_energy, staged_virial, staged_charge,
        committed_frc, committed_energy, committed_virial, true, true,
        candidate_s, candidate_t, committed_s_history, committed_t_history,
        history_count, true, error);
    if (error[0] != ReaxFFEvaluation::NONFINITE_STAGED_RESULT)
    {
        std::fprintf(stderr, "non-finite staging was not rejected: %d\n",
                     error[0]);
        return EXIT_FAILURE;
    }

    const bool unchanged =
        Equal_Array(committed_frc, original_frc, Equal_Vector) &&
        Equal_Floats(committed_energy, original_energy) &&
        Equal_Array(committed_virial, original_virial, Equal_Matrix) &&
        Equal_Floats(committed_local_charge, original_local_charge) &&
        Equal_Floats(committed_global_charge, original_global_charge) &&
        Equal_Floats(committed_s_history, original_s_history) &&
        Equal_Floats(committed_t_history, original_t_history);
    if (!unchanged)
    {
        std::fprintf(stderr,
                     "failed staging validation modified committed state\n");
        return EXIT_FAILURE;
    }

    staged_energy[1] = 2.0f;
    for (int& value : error) value = 0;
    ReaxFFEvaluation::Validate_Staging_Kernel(
        atoms, staged_frc, staged_energy, staged_virial, staged_charge,
        committed_frc, committed_energy, committed_virial, true, true,
        candidate_s, candidate_t, committed_s_history, committed_t_history,
        history_count, true, error);
    if (error[0] != ReaxFFEvaluation::EVALUATION_OK)
    {
        std::fprintf(stderr, "valid staging was rejected: %d\n", error[0]);
        return EXIT_FAILURE;
    }

    ReaxFFEvaluation::Commit_Kernel(
        atoms, local_to_global, staged_frc, staged_energy, staged_virial,
        staged_charge, committed_frc, committed_energy, committed_virial,
        committed_local_charge, committed_global_charge, true, true,
        candidate_s, candidate_t, committed_s_history, committed_t_history,
        history_count, history_capacity, true);

    for (int local = 0; local < atoms; local++)
    {
        const VECTOR expected_force = original_frc[local] + staged_frc[local];
        const LTMatrix3 expected_virial =
            original_virial[local] + staged_virial[local];
        if (!Equal_Vector(committed_frc[local], expected_force) ||
            committed_energy[local] !=
                original_energy[local] + staged_energy[local] ||
            !Equal_Matrix(committed_virial[local], expected_virial) ||
            committed_local_charge[local] != staged_charge[local] ||
            committed_global_charge[local_to_global[local]] !=
                staged_charge[local] ||
            committed_s_history[history_count * atoms +
                                local_to_global[local]] != candidate_s[local] ||
            committed_t_history[history_count * atoms +
                                local_to_global[local]] != candidate_t[local])
        {
            std::fprintf(stderr, "valid transaction commit failed at atom %d\n",
                         local);
            return EXIT_FAILURE;
        }
    }
    for (int i = 0; i < history_count * atoms; i++)
    {
        if (committed_s_history[i] != original_s_history[i] ||
            committed_t_history[i] != original_t_history[i])
        {
            std::fprintf(stderr, "history prefix changed during append\n");
            return EXIT_FAILURE;
        }
    }

    float s_history_before_no_commit[history_capacity * atoms] = {};
    float t_history_before_no_commit[history_capacity * atoms] = {};
    for (int i = 0; i < history_capacity * atoms; i++)
    {
        s_history_before_no_commit[i] = committed_s_history[i];
        t_history_before_no_commit[i] = committed_t_history[i];
    }
    ReaxFFEvaluation::Commit_Kernel(
        atoms, local_to_global, staged_frc, staged_energy, staged_virial,
        staged_charge, committed_frc, committed_energy, committed_virial,
        committed_local_charge, committed_global_charge, true, true,
        candidate_s, candidate_t, committed_s_history, committed_t_history,
        history_count + 1, history_capacity, false);
    if (!Equal_Floats(committed_s_history, s_history_before_no_commit) ||
        !Equal_Floats(committed_t_history, t_history_before_no_commit))
    {
        std::fprintf(stderr,
                     "commit_sampling_state=false modified EEQ history\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
