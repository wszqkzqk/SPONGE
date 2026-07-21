#pragma once

#include <cstdint>
#include <cstring>

#include "eigensolver_policy.hpp"

#include "../ecp/ecp_integrals.h"
#include "../structure/analytical_norms.hpp"

// 坐标同步
// 从 MD 坐标更新 QC 的原子环境与壳层中心（含周期边界修正）
#ifdef GPU_ARCH_NAME
static __device__ __forceinline__ bool QC_Geometry_Float_Is_Finite(float value)
{
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
}

static __device__ __forceinline__ unsigned int QC_Geometry_Float_Bits(
    float value)
{
    return __float_as_uint(value);
}
#else
static __host__ __device__ __forceinline__ bool QC_Geometry_Float_Is_Finite(
    float value)
{
    unsigned int bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires a 32-bit float representation");
    memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return (bits & 0x7f800000U) != 0x7f800000U;
}

static __host__ __device__ __forceinline__ unsigned int QC_Geometry_Float_Bits(
    float value)
{
    unsigned int bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}
#endif

enum QC_GEOMETRY_FAILURE_STAGE
{
    QC_GEOMETRY_FAILURE_NONE = -1,
    QC_GEOMETRY_FAILURE_RAW_COORDINATE = 0,
    QC_GEOMETRY_FAILURE_PREVIOUS_COORDINATE = 1,
    QC_GEOMETRY_FAILURE_ANCHOR_DISPLACEMENT = 2,
    QC_GEOMETRY_FAILURE_IMAGE_DISPLACEMENT = 3,
    QC_GEOMETRY_FAILURE_UNWRAPPED_COORDINATE = 4,
    QC_GEOMETRY_FAILURE_ANGSTROM_TO_BOHR = 5,
    QC_GEOMETRY_FAILURE_PAIR_DISPLACEMENT = 6,
};

static __host__ __device__ __forceinline__ float QC_Geometry_Component(
    const VECTOR& value, int axis)
{
    return axis == 0 ? value.x : (axis == 1 ? value.y : value.z);
}

static __device__ __forceinline__ bool QC_Claim_Geometry_Failure(int* failure,
                                                                 int stage)
{
#ifdef GPU_ARCH_NAME
    return atomicCAS(failure, QC_GEOMETRY_FAILURE_NONE, stage) ==
           QC_GEOMETRY_FAILURE_NONE;
#else
    bool claimed = false;
#pragma omp critical(sponge_qc_geometry_failure)
    {
        if (failure[0] == QC_GEOMETRY_FAILURE_NONE)
        {
            failure[0] = stage;
            claimed = true;
        }
    }
    return claimed;
#endif
}

static __device__ __forceinline__ void QC_Record_Geometry_Failure(
    int* failure, int stage, int local_atom, int axis, float value)
{
    if (!QC_Claim_Geometry_Failure(failure, stage)) return;
    failure[1] = local_atom;
    failure[2] = axis;
    failure[3] = (int)QC_Geometry_Float_Bits(value);
}

static __device__ __forceinline__ void QC_Record_Geometry_Pair_Failure(
    int* failure, int local_atom_i, int local_atom_j, int axis, float value)
{
    if (!QC_Claim_Geometry_Failure(
            failure, QC_GEOMETRY_FAILURE_PAIR_DISPLACEMENT))
        return;
    failure[1] = local_atom_i;
    failure[2] = local_atom_j;
    failure[3] = axis;
    failure[4] = (int)QC_Geometry_Float_Bits(value);
}

static __device__ __forceinline__ bool QC_Validate_Geometry_Vector(
    const VECTOR& value, int* failure, int stage, int local_atom)
{
    for (int axis = 0; axis < 3; ++axis)
    {
        const float component = QC_Geometry_Component(value, axis);
        if (!QC_Geometry_Float_Is_Finite(component))
        {
            QC_Record_Geometry_Failure(failure, stage, local_atom, axis,
                                       component);
            return false;
        }
    }
    return true;
}

static __global__ void QC_Initialize_Env_From_Crd_Kernel(
    const int natm, const int* atom_local, const VECTOR* crd, const int* atm,
    float* env, const float to_bohr, const VECTOR box_length,
    const int periodic_boundary, int* geometry_failure)
{
    SIMPLE_DEVICE_FOR(i, natm)
    {
        const VECTOR anchor = crd[atom_local[0]];
        const VECTOR r = crd[atom_local[i]];
        if (!QC_Validate_Geometry_Vector(
                anchor, geometry_failure,
                QC_GEOMETRY_FAILURE_RAW_COORDINATE, 0) ||
            !QC_Validate_Geometry_Vector(
                r, geometry_failure, QC_GEOMETRY_FAILURE_RAW_COORDINATE, i))
            continue;
        // Quantum-chemistry integrals are translationally invariant, while
        // their coordinates are stored in single precision.  Keep the first
        // QC atom at the origin so an arbitrary box-sized translation cannot
        // consume the significant bits needed by center differences.
        const VECTOR unwrapped =
            periodic_boundary
                ? Get_Periodic_Displacement(r, anchor, box_length)
                : r - anchor;
        if (!QC_Validate_Geometry_Vector(
                unwrapped, geometry_failure,
                QC_GEOMETRY_FAILURE_ANCHOR_DISPLACEMENT, i))
            continue;
        const int ptr_coord = atm[i * 6 + 1];
        const VECTOR converted(unwrapped.x * to_bohr,
                               unwrapped.y * to_bohr,
                               unwrapped.z * to_bohr);
        if (!QC_Validate_Geometry_Vector(
                converted, geometry_failure,
                QC_GEOMETRY_FAILURE_ANGSTROM_TO_BOHR, i))
            continue;
        env[ptr_coord + 0] = converted.x;
        env[ptr_coord + 1] = converted.y;
        env[ptr_coord + 2] = converted.z;
    }
}

static __global__ void QC_Update_Env_From_Crd_Kernel(
    const int natm, const int* atom_local, const VECTOR* crd, const int* atm,
    float* env, const float to_bohr, const VECTOR box_length,
    const int periodic_boundary, int* geometry_failure)
{
    SIMPLE_DEVICE_FOR(i, natm)
    {
        const VECTOR anchor = crd[atom_local[0]];
        const int md_idx = atom_local[i];
        const VECTOR r = crd[md_idx];
        if (!QC_Validate_Geometry_Vector(
                anchor, geometry_failure,
                QC_GEOMETRY_FAILURE_RAW_COORDINATE, 0) ||
            !QC_Validate_Geometry_Vector(
                r, geometry_failure, QC_GEOMETRY_FAILURE_RAW_COORDINATE, i))
            continue;
        const int ptr_coord = atm[i * 6 + 1];
        const VECTOR previous_bohr(env[ptr_coord + 0], env[ptr_coord + 1],
                                   env[ptr_coord + 2]);
        if (!QC_Validate_Geometry_Vector(
                previous_bohr, geometry_failure,
                QC_GEOMETRY_FAILURE_PREVIOUS_COORDINATE, i))
            continue;
        const VECTOR prev(previous_bohr.x / to_bohr,
                          previous_bohr.y / to_bohr,
                          previous_bohr.z / to_bohr);
        if (!QC_Validate_Geometry_Vector(
                prev, geometry_failure,
                QC_GEOMETRY_FAILURE_PREVIOUS_COORDINATE, i))
            continue;
        const VECTOR relative =
            periodic_boundary
                ? Get_Periodic_Displacement(r, anchor, box_length)
                : r - anchor;
        if (!QC_Validate_Geometry_Vector(
                relative, geometry_failure,
                QC_GEOMETRY_FAILURE_ANCHOR_DISPLACEMENT, i))
            continue;
        // Select the periodic image nearest the previously accepted relative
        // coordinate.  This preserves trajectory continuity and MC rollback
        // semantics without reintroducing an absolute translation.
        const VECTOR dr =
            periodic_boundary
                ? Get_Periodic_Displacement(relative, prev, box_length)
                : relative - prev;
        if (!QC_Validate_Geometry_Vector(
                dr, geometry_failure,
                QC_GEOMETRY_FAILURE_IMAGE_DISPLACEMENT, i))
            continue;
        const VECTOR unwrapped(prev.x + dr.x, prev.y + dr.y, prev.z + dr.z);
        if (!QC_Validate_Geometry_Vector(
                unwrapped, geometry_failure,
                QC_GEOMETRY_FAILURE_UNWRAPPED_COORDINATE, i))
            continue;
        const VECTOR converted(unwrapped.x * to_bohr,
                               unwrapped.y * to_bohr,
                               unwrapped.z * to_bohr);
        if (!QC_Validate_Geometry_Vector(
                converted, geometry_failure,
                QC_GEOMETRY_FAILURE_ANGSTROM_TO_BOHR, i))
            continue;
        env[ptr_coord + 0] = converted.x;
        env[ptr_coord + 1] = converted.y;
        env[ptr_coord + 2] = converted.z;
    }
}

static __global__ void QC_Update_Atom_Coords_From_Env_Kernel(
    const int natm, const int* atm, const float* env, VECTOR* atom_coords)
{
    SIMPLE_DEVICE_FOR(iat, natm)
    {
        const int ptr_coord = atm[iat * 6 + 1];
        atom_coords[iat] = {env[ptr_coord + 0], env[ptr_coord + 1],
                            env[ptr_coord + 2]};
    }
}

static __global__ void QC_Update_Centers_From_Env_Kernel(const int nbas,
                                                         const int* bas,
                                                         const int* atm,
                                                         const float* env,
                                                         VECTOR* centers)
{
    SIMPLE_DEVICE_FOR(ish, nbas)
    {
        const int iatm = bas[ish * 8 + 0];
        const int ptr_coord = atm[iatm * 6 + 1];
        centers[ish] = {env[ptr_coord + 0], env[ptr_coord + 1],
                        env[ptr_coord + 2]};
    }
}

static inline bool QC_Geometry_Host_Float_Is_Finite(float value)
{
    unsigned int bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return (bits & 0x7f800000U) != 0x7f800000U;
}

static bool QC_Validate_Periodic_Box_For_Geometry(
    CONTROLLER* controller, const VECTOR& box_length, int periodic_boundary,
    const char* operation)
{
    if (!periodic_boundary) return true;
    static const char axis_names[] = {'x', 'y', 'z'};
    for (int axis = 0; axis < 3; ++axis)
    {
        const float value = QC_Geometry_Component(box_length, axis);
        if (!QC_Geometry_Host_Float_Is_Finite(value) || !(value > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, operation,
                "Reason:\n    invalid periodic QC box at raw input: axis %c, "
                "value %.9g Angstrom; each component must be finite and > 0\n",
                axis_names[axis], (double)value);
            return false;
        }
        const float converted = value * CONSTANT_ANGSTROM_TO_BOHR;
        if (!QC_Geometry_Host_Float_Is_Finite(converted))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, operation,
                "Reason:\n    invalid periodic QC box after Angstrom-to-Bohr "
                "conversion: axis %c, raw value %.9g Angstrom, converted "
                "value %.9g Bohr\n",
                axis_names[axis], (double)value, (double)converted);
            return false;
        }
    }
    return true;
}

static bool QC_Report_Geometry_Transform_Failure(
    CONTROLLER* controller, const int* d_geometry_failure,
    const std::vector<int>& atom_local, const char* operation)
{
    int failure[5] = {QC_GEOMETRY_FAILURE_NONE, -1, -1, 0, 0};
    deviceMemcpy(failure, d_geometry_failure, sizeof(failure),
                 deviceMemcpyDeviceToHost);
    if (failure[0] == QC_GEOMETRY_FAILURE_NONE) return true;

    static const char* stage_names[] = {
        "raw-coordinate input",
        "previous accepted coordinate",
        "anchor subtraction/periodic displacement",
        "continuous-image displacement",
        "unwrapped coordinate accumulation",
        "Angstrom-to-Bohr conversion",
        "pair displacement",
    };
    const int stage = failure[0];
    if (stage == QC_GEOMETRY_FAILURE_PAIR_DISPLACEMENT)
    {
        const int local_atom_i = failure[1];
        const int local_atom_j = failure[2];
        const int global_atom_i =
            local_atom_i >= 0 &&
                    local_atom_i < static_cast<int>(atom_local.size())
                ? atom_local[local_atom_i]
                : -1;
        const int global_atom_j =
            local_atom_j >= 0 &&
                    local_atom_j < static_cast<int>(atom_local.size())
                ? atom_local[local_atom_j]
                : -1;
        const int axis = failure[3];
        const char axis_name = axis == 0 ? 'x' : (axis == 1 ? 'y' : 'z');
        const unsigned int value_bits =
            static_cast<unsigned int>(failure[4]);
        float value = 0.0f;
        std::memcpy(&value, &value_bits, sizeof(value));
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, operation,
            "Reason:\n    non-finite QC geometry during pair displacement: "
            "local atoms %d/%d (global atoms %d/%d), axis %c, value %.9g "
            "Bohr\n",
            local_atom_i, local_atom_j, global_atom_i, global_atom_j,
            axis_name, (double)value);
        return false;
    }
    const int local_atom = failure[1];
    const int axis = failure[2];
    const char* stage_name =
        stage >= QC_GEOMETRY_FAILURE_RAW_COORDINATE &&
                stage <= QC_GEOMETRY_FAILURE_PAIR_DISPLACEMENT
            ? stage_names[stage]
            : "unknown coordinate transformation";
    const int global_atom =
        local_atom >= 0 && local_atom < static_cast<int>(atom_local.size())
            ? atom_local[local_atom]
            : -1;
    const char axis_name = axis == 0 ? 'x' : (axis == 1 ? 'y' : 'z');
    const unsigned int value_bits = static_cast<unsigned int>(failure[3]);
    float value = 0.0f;
    std::memcpy(&value, &value_bits, sizeof(value));
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown, operation,
        "Reason:\n    non-finite QC geometry during %s: local atom %d "
        "(global atom %d), axis %c, value %.9g\n",
        stage_name, local_atom, global_atom, axis_name, (double)value);
    return false;
}

void QUANTUM_CHEMISTRY::Refresh_Coordinate_Derived_State_From_Env()
{
    const int threads = 256;
    Launch_Device_Kernel(QC_Update_Centers_From_Env_Kernel,
                         Positive_Int_Ceil_Div(mol.nbas, threads), threads, 0,
                         0,
                         mol.nbas, mol.d_bas, mol.d_atm, mol.d_env,
                         mol.d_centers);
    Launch_Device_Kernel(QC_Update_Atom_Coords_From_Env_Kernel,
                         Positive_Int_Ceil_Div(mol.natm, threads), threads, 0,
                         0,
                         mol.natm, mol.d_atm, mol.d_env, mol.d_atom_coords);
    if (scf_ws.ri.enabled)
    {
        auto& ri = scf_ws.ri;
        Launch_Device_Kernel(QC_Update_Centers_From_Env_Kernel,
                             Positive_Int_Ceil_Div(ri.naux_bas, threads),
                             threads, 0, 0, ri.naux_bas, ri.d_aux_bas,
                             ri.d_aux_atm,
                             ri.d_aux_env, ri.d_aux_centers);
    }
}

bool QUANTUM_CHEMISTRY::Initialize_Coordinates_From_MD(
    const VECTOR* crd, const VECTOR box_length)
{
    if (!QC_Validate_Periodic_Box_For_Geometry(
            controller, box_length, periodic_boundary,
            "QUANTUM_CHEMISTRY::Initialize_Coordinates_From_MD"))
        return false;
    const int threads = 256;
    deviceMemset(d_nuclear_geometry_failure, -1, 5 * sizeof(int));
    Launch_Device_Kernel(QC_Initialize_Env_From_Crd_Kernel,
                         Positive_Int_Ceil_Div(mol.natm, threads), threads, 0,
                         0,
                         mol.natm, d_atom_local, crd, mol.d_atm, mol.d_env,
                         CONSTANT_ANGSTROM_TO_BOHR, box_length,
                         periodic_boundary, d_nuclear_geometry_failure);
    if (!QC_Report_Geometry_Transform_Failure(
            controller, d_nuclear_geometry_failure, atom_local,
            "QUANTUM_CHEMISTRY::Initialize_Coordinates_From_MD"))
        return false;
    Launch_Device_Kernel(QC_Update_Centers_From_Env_Kernel,
                         Positive_Int_Ceil_Div(mol.nbas, threads), threads, 0,
                         0,
                         mol.nbas, mol.d_bas, mol.d_atm, mol.d_env,
                         mol.d_centers);
    Launch_Device_Kernel(QC_Update_Atom_Coords_From_Env_Kernel,
                         Positive_Int_Ceil_Div(mol.natm, threads), threads, 0,
                         0,
                         mol.natm, mol.d_atm, mol.d_env, mol.d_atom_coords);

    if (scf_ws.ri.enabled)
    {
        auto& ri = scf_ws.ri;
        deviceMemset(d_nuclear_geometry_failure, -1, 5 * sizeof(int));
        Launch_Device_Kernel(QC_Initialize_Env_From_Crd_Kernel,
                             Positive_Int_Ceil_Div(mol.natm, threads), threads,
                             0, 0, mol.natm, d_atom_local, crd, ri.d_aux_atm,
                             ri.d_aux_env, CONSTANT_ANGSTROM_TO_BOHR,
                             box_length, periodic_boundary,
                             d_nuclear_geometry_failure);
        if (!QC_Report_Geometry_Transform_Failure(
                controller, d_nuclear_geometry_failure, atom_local,
                "QUANTUM_CHEMISTRY::Initialize_Coordinates_From_MD"))
            return false;
        Launch_Device_Kernel(QC_Update_Centers_From_Env_Kernel,
                             Positive_Int_Ceil_Div(ri.naux_bas, threads),
                             threads, 0, 0, ri.naux_bas, ri.d_aux_bas,
                             ri.d_aux_atm,
                             ri.d_aux_env, ri.d_aux_centers);
    }
    return true;
}

bool QUANTUM_CHEMISTRY::Update_Coordinates_From_MD(
    const VECTOR* crd, const VECTOR box_length)
{
    if (!QC_Validate_Periodic_Box_For_Geometry(
            controller, box_length, periodic_boundary,
            "QUANTUM_CHEMISTRY::Update_Coordinates_From_MD"))
        return false;
    const int threads = 256;
    deviceMemset(d_nuclear_geometry_failure, -1, 5 * sizeof(int));
    Launch_Device_Kernel(QC_Update_Env_From_Crd_Kernel,
                         Positive_Int_Ceil_Div(mol.natm, threads), threads, 0,
                         0,
                         mol.natm, d_atom_local, crd, mol.d_atm, mol.d_env,
                         CONSTANT_ANGSTROM_TO_BOHR, box_length,
                         periodic_boundary, d_nuclear_geometry_failure);
    if (!QC_Report_Geometry_Transform_Failure(
            controller, d_nuclear_geometry_failure, atom_local,
            "QUANTUM_CHEMISTRY::Update_Coordinates_From_MD"))
        return false;
    Launch_Device_Kernel(QC_Update_Centers_From_Env_Kernel,
                         Positive_Int_Ceil_Div(mol.nbas, threads), threads, 0,
                         0,
                         mol.nbas, mol.d_bas, mol.d_atm, mol.d_env,
                         mol.d_centers);
    // 同步原子坐标 (按原子索引, ECP 内核使用)
    Launch_Device_Kernel(QC_Update_Atom_Coords_From_Env_Kernel,
                         Positive_Int_Ceil_Div(mol.natm, threads), threads, 0,
                         0,
                         mol.natm, mol.d_atm, mol.d_env, mol.d_atom_coords);

    // RI: 同步辅助基坐标（原子相同，atm 格式相同）
    if (scf_ws.ri.enabled)
    {
        auto& ri = scf_ws.ri;
        deviceMemset(d_nuclear_geometry_failure, -1, 5 * sizeof(int));
        Launch_Device_Kernel(
            QC_Update_Env_From_Crd_Kernel,
            Positive_Int_Ceil_Div(mol.natm, threads), threads, 0, 0, mol.natm,
            d_atom_local, crd, ri.d_aux_atm,
            ri.d_aux_env, CONSTANT_ANGSTROM_TO_BOHR, box_length,
            periodic_boundary, d_nuclear_geometry_failure);
        if (!QC_Report_Geometry_Transform_Failure(
                controller, d_nuclear_geometry_failure, atom_local,
                "QUANTUM_CHEMISTRY::Update_Coordinates_From_MD"))
            return false;
        Launch_Device_Kernel(QC_Update_Centers_From_Env_Kernel,
                             Positive_Int_Ceil_Div(ri.naux_bas, threads),
                             threads, 0, 0, ri.naux_bas, ri.d_aux_bas,
                             ri.d_aux_atm,
                             ri.d_aux_env, ri.d_aux_centers);
    }
    return true;
}

static __device__ __forceinline__ bool QC_Claim_Nuclear_Overlap(
    int* overlap_pair, int atom_i)
{
#ifdef GPU_ARCH_NAME
    return atomicCAS(overlap_pair, -1, atom_i) == -1;
#else
    bool claimed = false;
#pragma omp critical(sponge_qc_nuclear_overlap)
    {
        if (overlap_pair[0] == -1)
        {
            overlap_pair[0] = atom_i;
            claimed = true;
        }
    }
    return claimed;
#endif
}

static __global__ void QC_Validate_Nuclear_Geometry_Kernel(
    const int natm, const int* atm, const float* env, const VECTOR box_length,
    const int periodic_boundary, int* overlap_pair, int* geometry_failure)
{
    SIMPLE_DEVICE_FOR(i, natm)
    {
        const int ptr_i = atm[i * 6 + 1];
        const VECTOR ri(env[ptr_i + 0], env[ptr_i + 1], env[ptr_i + 2]);
        for (int j = i + 1; j < natm; j++)
        {
            const int ptr_j = atm[j * 6 + 1];
            const VECTOR rj(env[ptr_j + 0], env[ptr_j + 1], env[ptr_j + 2]);
            const VECTOR dr =
                periodic_boundary
                    ? Get_Periodic_Displacement(ri, rj, box_length)
                    : ri - rj;
            bool finite_displacement = true;
            for (int axis = 0; axis < 3; ++axis)
            {
                const float component = QC_Geometry_Component(dr, axis);
                if (!QC_Geometry_Float_Is_Finite(component))
                {
                    QC_Record_Geometry_Pair_Failure(
                        geometry_failure, i, j, axis, component);
                    finite_displacement = false;
                    break;
                }
            }
            if (!finite_displacement) continue;
            const double r2 = (double)dr.x * dr.x + (double)dr.y * dr.y +
                              (double)dr.z * dr.z;
            if (r2 == 0.0 && QC_Claim_Nuclear_Overlap(overlap_pair, i))
                overlap_pair[1] = j;
        }
    }
}

bool QUANTUM_CHEMISTRY::Validate_Nuclear_Geometry(const VECTOR box_length)
{
    if (!QC_Validate_Periodic_Box_For_Geometry(
            controller, box_length, periodic_boundary,
            "QUANTUM_CHEMISTRY::Validate_Nuclear_Geometry"))
        return false;
    deviceMemset(d_nuclear_overlap_pair, -1, 2 * sizeof(int));
    deviceMemset(d_nuclear_geometry_failure, -1, 5 * sizeof(int));
    const int threads = 256;
    const VECTOR box_bohr(box_length.x * CONSTANT_ANGSTROM_TO_BOHR,
                          box_length.y * CONSTANT_ANGSTROM_TO_BOHR,
                          box_length.z * CONSTANT_ANGSTROM_TO_BOHR);
    Launch_Device_Kernel(
        QC_Validate_Nuclear_Geometry_Kernel,
        Positive_Int_Ceil_Div(mol.natm, threads), threads, 0, 0, mol.natm,
        mol.d_atm, mol.d_env, box_bohr, periodic_boundary,
        d_nuclear_overlap_pair, d_nuclear_geometry_failure);
    if (!QC_Report_Geometry_Transform_Failure(
            controller, d_nuclear_geometry_failure, atom_local,
            "QUANTUM_CHEMISTRY::Validate_Nuclear_Geometry"))
        return false;
    int overlap_pair[2] = {-1, -1};
    deviceMemcpy(overlap_pair, d_nuclear_overlap_pair, 2 * sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (overlap_pair[0] >= 0)
    {
        const int atom_i = overlap_pair[0];
        const int atom_j = overlap_pair[1];
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Validate_Nuclear_Geometry",
            "Reason:\n    QC nuclei %d (global atom %d) and %d (global atom "
            "%d) overlap exactly under %s boundary geometry\n",
            atom_i, atom_local[atom_i], atom_j, atom_local[atom_j],
            periodic_boundary ? "periodic" : "non-periodic");
        return false;
    }
    return true;
}

// SCF 状态重置
// 清零收敛标志、能量缓存并重置 DIIS 历史
// 保留密度矩阵 P（复用上一步 MD 的收敛密度作为初猜）
void QUANTUM_CHEMISTRY::Reset_SCF_State()
{
    scf_ws.diis.diis_hist_count = scf_ws.diis.diis_hist_head = 0;
    scf_ws.diis.last_enorm = 1.0e10;
    scf_ws.runtime.convergence_streak = 0;
    scf_ws.runtime.convergence = QC_SCF_Convergence_State{};
    QC_SCF_Clear_Ensemble_Active_Set(scf_ws);
    scf_ws.ensemble.phase = QC_SCF_ENSEMBLE_INACTIVE;
    scf_ws.ensemble.bracket = QC_SCF_Ensemble_Bracket{};
    scf_ws.ensemble.line_origin_derivative = 0.0;
    scf_ws.ensemble.line_origin_energy = 0.0;
    scf_ws.ensemble.direction_density_rms = 0.0;
    scf_ws.ensemble.line_derivative_tolerance =
        scf_ws.runtime.energy_tol;
    scf_ws.ensemble.line_energy_guard_multiplier = 2.0;
    scf_ws.ensemble.current_fraction = 0.0;
    scf_ws.ensemble.maximum_fraction = 1.0;
    scf_ws.ensemble.committed_energy = 0.0;
    scf_ws.ensemble.verification_global_fw_gap = 0.0;
    scf_ws.ensemble.verification_commutator = 0.0;
    scf_ws.ensemble.active_fw_gap = 0.0;
    scf_ws.ensemble.spectral_orbital_origin_commutator = 0.0;
    scf_ws.ensemble.spectral_orbital_repair_linear_component = 0.0;
    scf_ws.ensemble.probe_evaluations = 0;
    scf_ws.ensemble.confirmed = false;
    scf_ws.ensemble.interior_minimum_confirmed = false;

    deviceMemset(scf_ws.core.d_scf_energy, 0, sizeof(double));
    deviceMemset(scf_ws.runtime.d_prev_energy, 0, sizeof(double));
    deviceMemset(scf_ws.runtime.d_delta_e, 0, sizeof(double));
    deviceMemset(scf_ws.runtime.d_density_residual, 0, sizeof(double));
    deviceMemset(scf_ws.runtime.d_e, 0, sizeof(double));
    if (scf_ws.runtime.unrestricted)
        deviceMemset(scf_ws.runtime.d_e_b, 0, sizeof(double));
    deviceMemset(scf_ws.runtime.d_pvxc, 0, sizeof(double));
    deviceMemset(scf_ws.runtime.d_converged, 0, sizeof(int));

    // Reset incremental Fock state
    const size_t fb = sizeof(float) * mol.nao2;
    const size_t db = sizeof(double) * mol.nao2;
    if (scf_ws.direct.d_P_coul_prev)
        deviceMemset(scf_ws.direct.d_P_coul_prev, 0, fb);
    if (scf_ws.direct.d_F_eri_accum)
        deviceMemset(scf_ws.direct.d_F_eri_accum, 0, db);
    if (scf_ws.direct.d_F_eri_accum_f)
        deviceMemset(scf_ws.direct.d_F_eri_accum_f, 0, fb);
    if (scf_ws.direct.d_P_exx_prev)
        deviceMemset(scf_ws.direct.d_P_exx_prev, 0, fb);
    if (scf_ws.direct.d_P_exx_b_prev)
        deviceMemset(scf_ws.direct.d_P_exx_b_prev, 0, fb);
    if (scf_ws.direct.d_F_eri_b_accum)
        deviceMemset(scf_ws.direct.d_F_eri_b_accum, 0, db);
    if (scf_ws.direct.d_F_eri_b_accum_f)
        deviceMemset(scf_ws.direct.d_F_eri_b_accum_f, 0, fb);
}

// 单电子积分
// 计算 S/T/V 单电子积分，并在球谐基下执行笛卡尔到球谐变换
void QUANTUM_CHEMISTRY::Compute_OneE_Integrals()
{
    const int nao_c = mol.nao_cart;
    float* p_S = mol.is_spherical ? cart2sph.d_S_cart : scf_ws.core.d_S;
    float* p_T = mol.is_spherical ? cart2sph.d_T_cart : scf_ws.core.d_T;
    float* p_V = mol.is_spherical ? cart2sph.d_V_cart : scf_ws.core.d_V;

    deviceMemset(p_S, 0, sizeof(float) * nao_c * nao_c);
    deviceMemset(p_T, 0, sizeof(float) * nao_c * nao_c);
    deviceMemset(p_V, 0, sizeof(float) * nao_c * nao_c);

    const int chunk_size = ONE_E_BATCH_SIZE;
    for (int i = 0; i < task_ctx.topo.n_1e_tasks; i += chunk_size)
    {
        int current_chunk = std::min(chunk_size, task_ctx.topo.n_1e_tasks - i);
        QC_ONE_E_TASK* task_ptr = task_ctx.buffers.d_1e_tasks + i;
        Launch_Device_Kernel(
            OneE_Kernel, (current_chunk + 63) / 64, 64, 0, 0, current_chunk,
            task_ptr, mol.d_centers, mol.d_l_list, mol.d_exps, mol.d_coeffs,
            mol.d_shell_offsets, mol.d_shell_sizes, mol.d_ao_offsets, mol.d_atm,
            mol.d_env, mol.natm, p_S, p_T, p_V, nao_c);
    }
    Cart2Sph_OneE_Integrals();
}

// ECP 矩阵
// 计算 V_ECP 并变换到有效基下, 归一化后加入 H_core
void QUANTUM_CHEMISTRY::Compute_ECP_Matrix()
{
    if (!mol.has_ecp) return;

    const int nao_c = mol.nao_cart;
    const int nao = mol.nao;

    // 在 Cartesian 基下计算 V_ECP
    float* d_V_ECP_cart = scf_ws.core.d_V_ECP;
    if (mol.is_spherical)
    {
        // 复用 cart2sph 临时缓冲 (d_V_cart)
        d_V_ECP_cart = cart2sph.d_V_cart;
    }
    deviceMemset(d_V_ECP_cart, 0, sizeof(float) * nao_c * nao_c);
    QC_ECP_EVALUATION_FAILURE failure;
    if (!QC_Compute_V_ECP(mol, task_ctx, d_V_ECP_cart, &failure))
    {
        if (failure.kind == QC_ECP_RESOURCE_ALLOCATION_FAILED)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Compute_ECP_Matrix",
                "Reason:\n    %s\n",
                QC_ECP_Evaluation_Failure_Kind_Name(failure.kind));
            return;
        }
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Compute_ECP_Matrix",
            "Reason:\n    ECP matrix evaluation failed: %s. Maximum angular "
            "series order q=%d; context: ECP atom %d, channel %d (l=%d), "
            "term %d (n_k=%d), "
            "one-electron task %d, shells (%d,%d), primitives (%d,%d), "
            "Cartesian functions (%d,%d); reported value %.17g, contracted "
            "error estimate %.9g\n",
            QC_ECP_Evaluation_Failure_Kind_Name(failure.kind),
            QC_ECP_MAX_SERIES_ORDER, failure.atom, failure.channel,
            failure.channel_l, failure.term, failure.n_k, failure.task_id,
            failure.shell_i, failure.shell_j, failure.primitive_i,
            failure.primitive_j, failure.cartesian_i, failure.cartesian_j,
            failure.value, failure.estimated_error);
        return;
    }

    // 球谐变换
    if (mol.is_spherical)
        Cart2Sph_Single_Matrix(d_V_ECP_cart, scf_ws.core.d_V_ECP);
}

// 核排斥能
// 累加核间库仑排斥能，结果写入设备侧 d_nuc_energy_dev
static __global__ void QC_Accumulate_Nuclear_Repulsion_Kernel(
    const int natm, const int* z_nuc, const int* atm, const float* env,
    double* e_nuc, const VECTOR box_length, const int periodic_boundary)
{
    SIMPLE_DEVICE_FOR(i, natm)
    {
        const int ptr_i = atm[i * 6 + 1];
        const double zi = (double)z_nuc[i];
        const VECTOR ri(env[ptr_i + 0], env[ptr_i + 1], env[ptr_i + 2]);
        double local = 0.0;
        for (int j = i + 1; j < natm; j++)
        {
            const int ptr_j = atm[j * 6 + 1];
            const double zj = (double)z_nuc[j];
            const VECTOR rj(env[ptr_j + 0], env[ptr_j + 1], env[ptr_j + 2]);
            const VECTOR dr =
                periodic_boundary
                    ? Get_Periodic_Displacement(ri, rj, box_length)
                    : ri - rj;
            const double r = sqrt((double)dr.x * dr.x + (double)dr.y * dr.y +
                                  (double)dr.z * dr.z);
            local += zi * zj / r;
        }
        atomicAdd(e_nuc, local);
    }
}

void QUANTUM_CHEMISTRY::Compute_Nuclear_Repulsion(const VECTOR box_length)
{
    deviceMemset(scf_ws.core.d_nuc_energy_dev, 0, sizeof(double));
    const int threads = 256;
    const VECTOR box_bohr(box_length.x * CONSTANT_ANGSTROM_TO_BOHR,
                          box_length.y * CONSTANT_ANGSTROM_TO_BOHR,
                          box_length.z * CONSTANT_ANGSTROM_TO_BOHR);
    Launch_Device_Kernel(QC_Accumulate_Nuclear_Repulsion_Kernel,
                         Positive_Int_Ceil_Div(mol.natm, threads), threads, 0,
                         0,
                         mol.natm, mol.d_Z, mol.d_atm, mol.d_env,
                         scf_ws.core.d_nuc_energy_dev, box_bohr,
                         periodic_boundary);
}

// 积分预处理
// 归一化单电子积分并构建 Hcore；双电子积分在 Build_Fock 中 direct 计算

void QUANTUM_CHEMISTRY::Compute_Analytical_Norms()
{
    try
    {
        std::vector<float> cart2sph_host;
        if (mol.is_spherical)
            cart2sph_host =
                QC_Build_Cart2Sph_Mat_Host(mol.h_l_list, mol.nao_cart, mol.nao);
        const std::vector<float> norms = qc_analytical_norms::Build(
            mol.h_l_list, mol.h_shell_sizes, mol.h_shell_offsets, mol.h_exps,
            mol.h_coeffs, mol.h_ao_offsets, mol.h_ao_offsets_sph,
            mol.is_spherical != 0, cart2sph_host, mol.nao_cart, mol.nao);
        deviceMemcpy(scf_ws.ortho.d_norms, norms.data(),
                     sizeof(float) * norms.size(), deviceMemcpyHostToDevice);
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "QUANTUM_CHEMISTRY::Compute_Analytical_Norms",
            "Reason:\n    Failed to build finite positive AO "
            "normalizations: %s\n",
            error.what());
    }
}

void QUANTUM_CHEMISTRY::Build_Shell_Pair_Bounds()
{
    if (task_ctx.topo.n_shell_pairs <= 0) return;
    // 固定 grid 大小，scratch 池按池槽数分配 (见 QC_BOUNDS_POOL_SLOTS)
    const int threads = 64;
    const int blocks = QC_BOUNDS_POOL_SLOTS / threads;
    const int n = task_ctx.topo.n_shell_pairs;
    Launch_Device_Kernel(
        QC_Build_Shell_Pair_Bounds_Kernel, blocks, threads, 0, 0, n,
        task_ctx.buffers.d_shell_pairs, mol.d_atm, mol.d_bas, mol.d_env,
        mol.d_ao_offsets, mol.d_ao_offsets_sph, scf_ws.ortho.d_norms,
        mol.is_spherical, cart2sph.d_cart2sph_mat, mol.nao_sph,
        task_ctx.buffers.d_shell_pair_bounds, scf_ws.direct.d_hr_pool,
        task_ctx.params.eri_hr_base, task_ctx.params.eri_hr_size,
        task_ctx.params.eri_shell_buf_size,
        task_ctx.params.eri_prim_screen_tol);
    task_ctx.topo.h_shell_pair_bounds.resize((size_t)n);
    deviceMemcpy(task_ctx.topo.h_shell_pair_bounds.data(),
                 task_ctx.buffers.d_shell_pair_bounds, sizeof(float) * n,
                 deviceMemcpyDeviceToHost);
}

static __global__ void QC_Scale_OneE_And_Build_Hcore_Kernel(
    const int nao, const float* norms, float* S, float* T, float* V,
    float* V_ECP, float* H_core)
{
    const int total = nao * nao;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        int i = idx / nao;
        int j = idx - i * nao;
        float scale = norms[i] * norms[j];
        S[idx] *= scale;
        T[idx] *= scale;
        V[idx] *= scale;
        H_core[idx] = T[idx] + V[idx];
        if (V_ECP)
        {
            V_ECP[idx] *= scale;
            H_core[idx] += V_ECP[idx];
        }
    }
}

void QUANTUM_CHEMISTRY::Prepare_Integrals()
{
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    const int threads = 256;

    // norms 已在 Compute_Analytical_Norms 中计算完毕
    // 归一化 S/T/V(/V_ECP) 并构建 Hcore = T + V (+ V_ECP)
    Launch_Device_Kernel(QC_Scale_OneE_And_Build_Hcore_Kernel,
                         Positive_Int_Ceil_Div(nao2, threads), threads, 0, 0,
                         nao,
                         scf_ws.ortho.d_norms, scf_ws.core.d_S, scf_ws.core.d_T,
                         scf_ws.core.d_V, scf_ws.core.d_V_ECP,
                         scf_ws.core.d_H_core);
}

// 重叠正交化矩阵
// 对重叠矩阵 S 做 double 精度本征分解，并构建正交化变换矩阵 X
bool QUANTUM_CHEMISTRY::Build_Overlap_X()
{
    const int nao = mol.nao;
    const int nao2 = mol.nao2;

    QC_Float_To_Double(nao2, scf_ws.core.d_S, scf_ws.ortho.d_dwork_nao2_1);
    // Integral kernels evaluate the two AO triangles independently.  Work
    // with their symmetric part everywhere: the generalized eigensolver,
    // particle counts, commutators, and ensemble occupation transforms must
    // all use one and the same overlap metric.
    QC_Symmetrize_Double_Matrix(nao,
                                scf_ws.ortho.d_dwork_nao2_1);
    QC_Double_To_Float(nao2, scf_ws.ortho.d_dwork_nao2_1,
                       scf_ws.core.d_S);

    int info = 0;
    const int api_status = QC_Diagonalize_Double(
        solver_handle, nao, scf_ws.ortho.d_dwork_nao2_1,
        scf_ws.ortho.d_dW_double, scf_ws.ortho.d_solver_work_double,
        scf_ws.ortho.lwork_double, &info);
    const bool solver_ok = QC_SCF_Require_Eigensolver_Success(
        QC_SCF_EIGENSOLVER_OVERLAP, QC_SCF_EIGENSOLVER_CHANNEL_OVERLAP,
        nao, api_status, info,
        [&](const QC_SCF_Eigensolver_Failure& failure)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Build_Overlap_X",
                "Reason:\n    eigensolver failed during %s for channel %s: "
                "dimension=%d, api_status=%d, info=%d\n",
                failure.stage_name, failure.channel_name, failure.dimension,
                failure.api_status, failure.info);
        });
    if (!solver_ok) return false;

    std::vector<double> h_W(nao);
    deviceMemcpy(h_W.data(), scf_ws.ortho.d_dW_double, sizeof(double) * nao,
                 deviceMemcpyDeviceToHost);
    const double lindep_thresh = scf_ws.ortho.lindep_threshold;
    if (!Double_Memory_Is_Finite(&lindep_thresh) ||
        !(lindep_thresh > 0.0))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Build_Overlap_X",
            "Reason:\n    invalid overlap linear-dependence threshold: "
            "%.17g\n",
            lindep_thresh);
        return false;
    }
    int nao_eff = 0;
    for (int k = 0; k < nao; k++)
    {
        if (!Double_Memory_Is_Finite(&h_W[k]))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Build_Overlap_X",
                "Reason:\n    non-finite overlap eigenvalue: channel=overlap, "
                "dimension=%d, eigenvalue=%d, value=%.17g\n",
                nao, k, h_W[k]);
            return false;
        }
        if (h_W[k] < -lindep_thresh ||
            (k > 0 && h_W[k] < h_W[k - 1]))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Build_Overlap_X",
                "Reason:\n    invalid overlap eigenspectrum: "
                "channel=overlap, dimension=%d, eigenvalue=%d, "
                "value=%.17g, previous=%.17g, threshold=%.17g\n",
                nao, k, h_W[k], k > 0 ? h_W[k - 1] : h_W[k],
                lindep_thresh);
            return false;
        }
        if (h_W[k] >= lindep_thresh) nao_eff++;
    }
    const int required_orbitals =
        std::max(scf_ws.runtime.n_alpha, scf_ws.runtime.n_beta);
    if (nao_eff <= 0 || nao_eff < required_orbitals)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Build_Overlap_X",
            "Reason:\n    overlap retained space cannot hold the occupied "
            "orbitals: channel=overlap, dimension=%d, retained=%d, "
            "required=%d, threshold=%.17g\n",
            nao, nao_eff, required_orbitals, lindep_thresh);
        return false;
    }
    scf_ws.ortho.nao_eff = nao_eff;

    QC_Double_To_Float(nao, scf_ws.ortho.d_dW_double, scf_ws.ortho.d_W);

    deviceMemset(scf_ws.ortho.d_X, 0, sizeof(double) * nao2);
    QC_Build_X_Canonical(nao, nao_eff, scf_ws.ortho.d_dwork_nao2_1,
                         scf_ws.ortho.d_dW_double, lindep_thresh,
                         scf_ws.ortho.d_X);
    return true;
}
