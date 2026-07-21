#pragma once

#include "reaxff_geometry.h"

namespace ReaxFFEvaluation
{

#ifdef GPU_ARCH_NAME
#define REAXFF_EVALUATION_SKIP_ATOM() return
#else
#define REAXFF_EVALUATION_SKIP_ATOM() continue
#endif

enum ERROR_CODE
{
    EVALUATION_OK = 0,
    INVALID_LOCAL_TO_GLOBAL = 100,
    DUPLICATE_GLOBAL_ATOM = 101,
    NONFINITE_COORDINATE = 102,
    INVALID_OR_INCONSISTENT_TYPE = 103,
    INVALID_HYDROGEN_FLAG = 104,
    INVALID_NEIGHBOR_ROW = 105,
    INVALID_NEIGHBOR_ATOM = 106,
    NONFINITE_STAGED_RESULT = 107,
    NONFINITE_COMMITTED_RESULT = 108,
    NONFINITE_HISTORY = 109
};

enum COORDINATE_COUNT_STATUS
{
    COORDINATE_COUNT_OK = 0,
    COORDINATE_COUNT_NEGATIVE = 1,
    COORDINATE_COUNT_INT_OVERFLOW = 2,
    COORDINATE_COUNT_SIZE_OVERFLOW = 3
};

// Validate DD's owned+ghost coordinate extent without allocating or touching
// device state.  bytes_per_coordinate is explicit so every caller can prove
// the allocation it is about to make; ReaxFF staging passes the size of its
// largest per-coordinate record.
static inline COORDINATE_COUNT_STATUS Checked_Coordinate_Count(
    int owned_count, int ghost_count, std::size_t bytes_per_coordinate,
    int* coordinate_count)
{
    if (coordinate_count == NULL || owned_count < 0 || ghost_count < 0)
        return COORDINATE_COUNT_NEGATIVE;
    *coordinate_count = 0;
    if (owned_count > INT_MAX - ghost_count)
        return COORDINATE_COUNT_INT_OVERFLOW;
    const int checked_count = owned_count + ghost_count;
    if (checked_count > 0 &&
        bytes_per_coordinate >
            std::numeric_limits<std::size_t>::max() /
                static_cast<std::size_t>(checked_count))
    {
        return COORDINATE_COUNT_SIZE_OVERFLOW;
    }
    *coordinate_count = checked_count;
    return COORDINATE_COUNT_OK;
}

// Validate every DD-local identity and every index that ReaxFF kernels will
// dereference before the first Hamiltonian kernel is launched.  ReaxFF is
// currently restricted to one PP rank with no ghosts, so local_atom_numbers is
// also the global atom count and a range-checked, duplicate-free mapping is a
// complete bijection.
static __global__ void Preflight_Kernel(
    int atom_numbers, const VECTOR* crd, const int* atom_local,
    const int* reference_type_local, const int* reference_type_global,
    int atom_type_numbers, const int* bond_order_type, const int* bond_type,
    const int* vdw_type, const int* ovun_type, const int* angle_type,
    const int* torsion_type, const int* hb_type,
    const int* is_hydrogen_local, const int* is_hydrogen_global,
    const ATOM_GROUP* half_nl, int half_capacity,
    const ATOM_GROUP* full_nl, int full_capacity, int* seen_global,
    int* evaluation_error)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        const VECTOR coordinate = crd[i];
        if (!ReaxFF_Vector_Is_Finite(coordinate))
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, NONFINITE_COORDINATE, i);
            REAXFF_EVALUATION_SKIP_ATOM();
        }

        const int global_id = atom_local[i];
        if (global_id < 0 || global_id >= atom_numbers)
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, INVALID_LOCAL_TO_GLOBAL, i, global_id);
            REAXFF_EVALUATION_SKIP_ATOM();
        }
        if (atomicAdd(seen_global + global_id, 1) != 0)
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, DUPLICATE_GLOBAL_ATOM, i, global_id);
            REAXFF_EVALUATION_SKIP_ATOM();
        }

        const int type = reference_type_local[i];
        if (type < 0 || type >= atom_type_numbers ||
            type != reference_type_global[global_id] ||
            type != bond_order_type[i] || type != bond_type[i] ||
            type != vdw_type[i] || type != ovun_type[i] ||
            type != angle_type[i] || type != torsion_type[i] ||
            type != hb_type[i])
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, INVALID_OR_INCONSISTENT_TYPE, i, global_id,
                type);
            REAXFF_EVALUATION_SKIP_ATOM();
        }

        const int hydrogen = is_hydrogen_local[i];
        if ((hydrogen != 0 && hydrogen != 1) ||
            hydrogen != is_hydrogen_global[global_id])
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, INVALID_HYDROGEN_FLAG, i, global_id,
                hydrogen);
            REAXFF_EVALUATION_SKIP_ATOM();
        }

        const ATOM_GROUP half = half_nl[i];
        if (half.atom_numbers < 0 || half.atom_numbers > half_capacity ||
            half.ghost_numbers != 0 ||
            (half.atom_numbers > 0 && half.atom_serial == NULL))
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, INVALID_NEIGHBOR_ROW, i, 0,
                half.atom_numbers, half.ghost_numbers);
            REAXFF_EVALUATION_SKIP_ATOM();
        }
        bool half_valid = true;
        for (int entry = 0; entry < half.atom_numbers; entry++)
        {
            const int neighbor = half.atom_serial[entry];
            if (neighbor < 0 || neighbor >= atom_numbers || neighbor == i)
            {
                Record_ReaxFF_Geometry_Error(
                    evaluation_error, INVALID_NEIGHBOR_ATOM, i, neighbor, 0,
                    entry);
                half_valid = false;
                break;
            }
        }
        if (!half_valid) REAXFF_EVALUATION_SKIP_ATOM();

        const ATOM_GROUP full = full_nl[i];
        if (full.atom_numbers < 0 || full.atom_numbers > full_capacity ||
            full.ghost_numbers != 0 ||
            (full.atom_numbers > 0 && full.atom_serial == NULL))
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, INVALID_NEIGHBOR_ROW, i, 1,
                full.atom_numbers, full.ghost_numbers);
            REAXFF_EVALUATION_SKIP_ATOM();
        }
        for (int entry = 0; entry < full.atom_numbers; entry++)
        {
            const int neighbor = full.atom_serial[entry];
            if (neighbor < 0 || neighbor >= atom_numbers || neighbor == i)
            {
                Record_ReaxFF_Geometry_Error(
                    evaluation_error, INVALID_NEIGHBOR_ATOM, i, neighbor, 1,
                    entry);
                break;
            }
        }
    }
}

// Check both the private ReaxFF result and the exact value that would be
// visible after addition to the already accumulated DD buffers.  No public
// state is touched by this kernel.
static __global__ void Validate_Staging_Kernel(
    int atom_numbers, const VECTOR* staged_frc, const float* staged_energy,
    const LTMatrix3* staged_virial, const float* staged_charge,
    const VECTOR* committed_frc, const float* committed_energy,
    const LTMatrix3* committed_virial, bool commit_energy, bool commit_virial,
    const float* candidate_s, const float* candidate_t,
    const float* committed_s_history, const float* committed_t_history,
    int history_count, bool commit_history, int* evaluation_error)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        const VECTOR staged_force = staged_frc[i];
        const float staged_atom_energy = staged_energy[i];
        const LTMatrix3 staged_atom_virial = staged_virial[i];
        if (!ReaxFF_Vector_Is_Finite(staged_force) ||
            !ReaxFF_Float_Is_Finite(staged_atom_energy) ||
            !ReaxFF_Matrix_Is_Finite(staged_atom_virial) ||
            !ReaxFF_Float_Is_Finite(staged_charge[i]))
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, NONFINITE_STAGED_RESULT, i);
            REAXFF_EVALUATION_SKIP_ATOM();
        }

        const VECTOR final_force = committed_frc[i] + staged_force;
        if (!ReaxFF_Vector_Is_Finite(committed_frc[i]) ||
            !ReaxFF_Vector_Is_Finite(final_force))
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, NONFINITE_COMMITTED_RESULT, i, 0);
            REAXFF_EVALUATION_SKIP_ATOM();
        }
        if (commit_energy &&
            (!ReaxFF_Float_Is_Finite(committed_energy[i]) ||
             !ReaxFF_Float_Is_Finite(committed_energy[i] +
                                     staged_atom_energy)))
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, NONFINITE_COMMITTED_RESULT, i, 1);
            REAXFF_EVALUATION_SKIP_ATOM();
        }
        if (commit_virial &&
            (!ReaxFF_Matrix_Is_Finite(committed_virial[i]) ||
             !ReaxFF_Matrix_Is_Finite(committed_virial[i] +
                                      staged_atom_virial)))
        {
            Record_ReaxFF_Geometry_Error(
                evaluation_error, NONFINITE_COMMITTED_RESULT, i, 2);
            REAXFF_EVALUATION_SKIP_ATOM();
        }

        if (commit_history)
        {
            if (!ReaxFF_Float_Is_Finite(candidate_s[i]) ||
                !ReaxFF_Float_Is_Finite(candidate_t[i]))
            {
                Record_ReaxFF_Geometry_Error(
                    evaluation_error, NONFINITE_HISTORY, i, history_count);
                REAXFF_EVALUATION_SKIP_ATOM();
            }
            for (int frame = 0; frame < history_count; frame++)
            {
                const int history_index = frame * atom_numbers + i;
                if (!ReaxFF_Float_Is_Finite(
                        committed_s_history[history_index]) ||
                    !ReaxFF_Float_Is_Finite(
                        committed_t_history[history_index]))
                {
                    Record_ReaxFF_Geometry_Error(
                        evaluation_error, NONFINITE_HISTORY, i, frame);
                    break;
                }
            }
        }
    }
}

// This is the sole publication point for a ReaxFF force evaluation.  The
// preflight and staging validator prove every indexed write and arithmetic
// result before this kernel is launched.  Each local/global/history element has
// exactly one writer.
static __global__ void Commit_Kernel(
    int atom_numbers, const int* atom_local, const VECTOR* staged_frc,
    const float* staged_energy, const LTMatrix3* staged_virial,
    const float* staged_charge, VECTOR* committed_frc,
    float* committed_energy, LTMatrix3* committed_virial,
    float* committed_local_charge, float* committed_global_charge,
    bool commit_energy, bool commit_virial, const float* candidate_s,
    const float* candidate_t, float* committed_s_history,
    float* committed_t_history, int history_count, int history_capacity,
    bool commit_history)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        committed_frc[i] = committed_frc[i] + staged_frc[i];
        if (commit_energy)
            committed_energy[i] += staged_energy[i];
        if (commit_virial)
            committed_virial[i] = committed_virial[i] + staged_virial[i];

        const int global_id = atom_local[i];
        const float charge = staged_charge[i];
        committed_local_charge[i] = charge;
        committed_global_charge[global_id] = charge;

        if (commit_history)
        {
            if (history_count < history_capacity)
            {
                const int target = history_count * atom_numbers + global_id;
                committed_s_history[target] = candidate_s[i];
                committed_t_history[target] = candidate_t[i];
            }
            else
            {
                for (int frame = 0; frame < history_capacity - 1; frame++)
                {
                    const int target = frame * atom_numbers + global_id;
                    const int source = target + atom_numbers;
                    committed_s_history[target] = committed_s_history[source];
                    committed_t_history[target] = committed_t_history[source];
                }
                const int target =
                    (history_capacity - 1) * atom_numbers + global_id;
                committed_s_history[target] = candidate_s[i];
                committed_t_history[target] = candidate_t[i];
            }
        }
    }
}

}  // namespace ReaxFFEvaluation

#undef REAXFF_EVALUATION_SKIP_ATOM
