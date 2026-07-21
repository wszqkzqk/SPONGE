#ifndef REAXFF_EEQ_H
#define REAXFF_EEQ_H

#include "../../common.h"
#include "../../control.h"

namespace ReaxFFEEQ
{

constexpr float COULOMB_CONSTANT = 332.05221729f;

struct Pair_Force_Result
{
    VECTOR force = {0.0f, 0.0f, 0.0f};
    bool active = false;
};

// Canonical shielded EEQ pair force used by the production kernel and its
// numerical regression probe.  Activity is based on exact input zero, not a
// magnitude threshold.  At exact overlap the r^3 shielding makes dH/dr
// O(r^2), so the Cartesian force has the unique zero limit.
static __device__ __host__ __forceinline__ Pair_Force_Result
Evaluate_Pair_Force(const VECTOR& drij, float qi, float qj, float shield,
                    float cutoff)
{
    Pair_Force_Result result;
    if (qi == 0.0f || qj == 0.0f) return result;
    result.active = true;

    const float r2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
    if (r2 == 0.0f) return result;
    const float r = sqrtf(r2);
    if (r >= cutoff) return result;

    const float inv_cutoff = 1.0f / cutoff;
    const float x = r * inv_cutoff;
    const float x2 = x * x;
    const float x3 = x2 * x;
    const float x4 = x2 * x2;
    const float x5 = x4 * x;
    const float x6 = x5 * x;
    const float x7 = x6 * x;
    const float taper =
        20.0f * x7 - 70.0f * x6 + 84.0f * x5 - 35.0f * x4 + 1.0f;
    const float dtaper_dr =
        inv_cutoff * (140.0f * x6 - 420.0f * x5 + 420.0f * x4 - 140.0f * x3);

    const float u = r2 * r + shield;
    const float inv_u_cbrt = 1.0f / cbrtf(u);
    const float dH_dr = COULOMB_CONSTANT *
                        (dtaper_dr * inv_u_cbrt - taper * r2 * inv_u_cbrt / u);
    const float force_mag = -qi * qj * dH_dr / r;
    result.force = {force_mag * drij.x, force_mag * drij.y, force_mag * drij.z};
    return result;
}

}  // namespace ReaxFFEEQ

struct REAXFF_EEQ
{
    int is_initialized = 0;
    CONTROLLER* controller = NULL;

    int atom_numbers = 0;
    int atom_type_numbers = 0;

    float* h_chi = NULL;
    float* h_eta = NULL;
    float* h_gamma = NULL;
    float* h_shield = NULL;

    float* d_chi = NULL;
    float* d_eta = NULL;
    float* d_gamma = NULL;
    float* d_shield = NULL;

    int* h_atom_type = NULL;
    // Immutable input-order table. d_atom_type is the current DD-local view.
    int* d_atom_type_global = NULL;
    int* d_atom_type = NULL;

    // Solver temporary arrays (pre-allocated)
    float* d_b = NULL;
    float* d_r = NULL;
    float* d_p = NULL;
    float* d_Ap = NULL;
    float* d_q = NULL;

    // Additional temporary arrays for CG
    float* d_s = NULL;
    float* d_t = NULL;
    float* d_z = NULL;  // Jacobi preconditioner temporary
    float* d_temp_sum = NULL;
    int* h_h_numnbrs = NULL;
    int* h_h_firstnbrs = NULL;
    int* d_h_numnbrs = NULL;
    int* d_h_firstnbrs = NULL;
    int* d_h_jlist = NULL;
    float* d_h_val = NULL;
    int h_matrix_capacity = 0;

    // Device-side CG scalar buffers (GPU optimization: avoid host sync)
    float* d_rr_old = NULL;
    float* d_rr_new = NULL;
    float* d_pAp_buf = NULL;
    float* d_cg_alpha = NULL;
    float* d_cg_beta = NULL;
    int* d_solver_error = NULL;

    // Convergence parameters
    float tolerance = 1e-4f;
    int max_iter = 1000;
    // Charge history for polynomial extrapolation (LAMMPS-style)
    enum
    {
        HIST_SIZE = 5
    };
    int nprev = 0;
    float* d_s_hist = NULL;
    float* d_t_hist = NULL;

    float h_energy = 0.0f;
    // Host-side result of the current private evaluation.  REAXFF publishes it
    // to h_energy only after every ReaxFF module has succeeded.
    float pending_energy = 0.0f;

    void Initial(CONTROLLER* controller, int atom_numbers,
                 const char* parameter_in_file, const char* type_in_file);
    void Calculate_Charges(int atom_numbers, const int* atom_local,
                           float* d_charge, const VECTOR* d_crd,
                           const LTMatrix3 cell,
                           const LTMatrix3 rcell, const ATOM_GROUP* fnl_d_nl,
                           float cutoff, float* d_energy = NULL,
                           VECTOR* frc = NULL, int need_virial = 0,
                           LTMatrix3* atom_virial = NULL);
    void Step_Print(CONTROLLER* controller);
    void Print_Charges(const float* d_charge);
};

#endif
