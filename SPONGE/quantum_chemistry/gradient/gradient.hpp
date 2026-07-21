#pragma once

#include "grad_nuclear.hpp"
#include "grad_wdensity.hpp"

// 力/维里回写 kernel
// 将 QC 原子梯度（Hartree/Bohr）转换为 MD 力（kcal/mol/Å）
// 并累加到 MD 力数组；可选地计算维里张量
static __global__ void QC_Writeback_Gradient_Kernel(
    const int natm, const int* qc_to_global, const int* global_to_local,
    const int owned_atom_numbers, const double* grad,
    const VECTOR* qc_coords_bohr, VECTOR* local_frc, const int need_virial,
    LTMatrix3* local_atom_virial)
{
    SIMPLE_DEVICE_FOR(i, natm)
    {
        const int global_atom = qc_to_global[i];
        const int local_atom = global_to_local[global_atom];
        if (local_atom < 0 || local_atom >= owned_atom_numbers) continue;
        const double ha_bohr_to_kcal_ang =
            (double)CONSTANT_HARTREE_TO_KCAL_MOL *
            (double)CONSTANT_ANGSTROM_TO_BOHR;
        const float fx = (float)(-grad[i * 3 + 0] * ha_bohr_to_kcal_ang);
        const float fy = (float)(-grad[i * 3 + 1] * ha_bohr_to_kcal_ang);
        const float fz = (float)(-grad[i * 3 + 2] * ha_bohr_to_kcal_ang);
        atomicAdd(&local_frc[local_atom].x, fx);
        atomicAdd(&local_frc[local_atom].y, fy);
        atomicAdd(&local_frc[local_atom].z, fz);

        if (need_virial && local_atom_virial != NULL)
        {
            const VECTOR f_vec = {fx, fy, fz};
            const VECTOR r_bohr = qc_coords_bohr[i];
            const float bohr_to_angstrom =
                1.0f / (float)CONSTANT_ANGSTROM_TO_BOHR;
            const VECTOR r_vec = {r_bohr.x * bohr_to_angstrom,
                                  r_bohr.y * bohr_to_angstrom,
                                  r_bohr.z * bohr_to_angstrom};
            atomicAdd(local_atom_virial + local_atom,
                      Get_Virial_From_Force_Dis(f_vec, r_vec));
        }
    }
}
