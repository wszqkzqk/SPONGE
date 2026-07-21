#include "eeq.h"

#include "atom_identity.h"
#include "reaxff_geometry.h"
#include "reaxff_input.h"

#ifdef USE_CPU
#define EEQ_SIMPLE_DEVICE_FOR(i, N)                           \
    PRAGMA(omp parallel for schedule(static) if ((N) >= 512)) \
    for (int i = 0; i < N; i++)
#else
#define EEQ_SIMPLE_DEVICE_FOR(i, N) SIMPLE_DEVICE_FOR(i, N)
#endif

#ifndef USE_CPU
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>
#endif

void REAXFF_EEQ::Initial(CONTROLLER* controller, int atom_numbers,
                         const char* parameter_in_file,
                         const char* type_in_file)
{
    if (parameter_in_file == NULL || type_in_file == NULL)
    {
        controller->printf(
            "REAXFF_EEQ IS NOT INITIALIZED (missing input files)\n\n");
        return;
    }

    controller->printf("START INITIALIZING REAXFF_EEQ\n");
    REAXFF_INPUT_ERROR input_error;
    REAXFF_FORCE_FIELD_IR force_field;
    if (!ReaxFF_Parse_Force_Field_File(parameter_in_file, &force_field,
                                       &input_error))
    {
        const std::string reason = input_error.Describe();
        controller->Throw_Formatted_SPONGE_Error(
            input_error.kind == REAXFF_INPUT_OPEN_ERROR
                ? spongeErrorOpenFileFailed
                : spongeErrorBadFileFormat,
            "REAXFF_EEQ::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }
    std::vector<int> atom_type;
    if (!ReaxFF_Parse_Type_File_Path(type_in_file, atom_numbers, force_field,
                                     &atom_type, NULL, &input_error))
    {
        const std::string reason = input_error.Describe();
        controller->Throw_Formatted_SPONGE_Error(
            input_error.kind == REAXFF_INPUT_OPEN_ERROR
                ? spongeErrorOpenFileFailed
                : spongeErrorBadFileFormat,
            "REAXFF_EEQ::Initial", "Reason:\n\t%s", reason.c_str());
        return;
    }

    const int n_atom_types = static_cast<int>(force_field.atom_types.size());
    int pair_parameter_count = 0;
    if (!ReaxFF_Checked_Dense_Table_Count(
            n_atom_types, 2, sizeof(float), &pair_parameter_count))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "REAXFF_EEQ::Initial",
            "Reason:\n\tatom type count %d exceeds the supported EEQ pair "
            "table extent in file %s",
            n_atom_types, parameter_in_file);
        return;
    }

    std::vector<float> chi(n_atom_types);
    std::vector<float> eta(n_atom_types);
    std::vector<float> gamma(n_atom_types);
    std::vector<float> shield(pair_parameter_count);
    for (int i = 0; i < n_atom_types; i++)
    {
        const REAXFF_ATOM_TYPE_IR& atom = force_field.atom_types[i];
        gamma[i] = atom.values[0][5];
        chi[i] = atom.values[1][5] * CONSTANT_EV_TO_KCAL_MOL;
        eta[i] = atom.values[1][6] * CONSTANT_EV_TO_KCAL_MOL * 2.0f;
        if (!Float_Memory_Is_Finite(&chi[i]) ||
            !Float_Memory_Is_Finite(&eta[i]) || !(eta[i] > 0.0f) ||
            !Float_Memory_Is_Finite(&gamma[i]) || !(gamma[i] > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "REAXFF_EEQ::Initial",
                "Reason:\n\tEEQ parameters for atom type %d must remain "
                "finite and eta/gamma must be positive in file %s",
                i + 1, parameter_in_file);
            return;
        }
    }
    for (int i = 0; i < n_atom_types; i++)
    {
        for (int j = 0; j < n_atom_types; j++)
        {
            const double gamma_product =
                static_cast<double>(gamma[i]) * static_cast<double>(gamma[j]);
            const double shield_double = pow(gamma_product, -1.5);
            const float shield_value = static_cast<float>(shield_double);
            if (!(gamma_product > 0.0) ||
                !Double_Memory_Is_Finite(&gamma_product) ||
                !(shield_double > 0.0) ||
                !Double_Memory_Is_Finite(&shield_double) ||
                !(shield_value > 0.0f) ||
                !Float_Memory_Is_Finite(&shield_value))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat, "REAXFF_EEQ::Initial",
                    "Reason:\n\tEEQ shielding parameters produce a "
                    "non-positive or unrepresentable value in file %s",
                    parameter_in_file);
                return;
            }
            shield[i * n_atom_types + j] = shield_value;
        }
    }

    float* staged_chi = NULL;
    float* staged_eta = NULL;
    float* staged_gamma = NULL;
    float* staged_shield = NULL;
    int* staged_atom_type = NULL;
    Malloc_Safely((void**)&staged_chi, sizeof(float) * n_atom_types);
    Malloc_Safely((void**)&staged_eta, sizeof(float) * n_atom_types);
    Malloc_Safely((void**)&staged_gamma, sizeof(float) * n_atom_types);
    Malloc_Safely((void**)&staged_shield,
                  sizeof(float) * pair_parameter_count);
    Malloc_Safely((void**)&staged_atom_type, sizeof(int) * atom_numbers);
    memcpy(staged_chi, chi.data(), sizeof(float) * n_atom_types);
    memcpy(staged_eta, eta.data(), sizeof(float) * n_atom_types);
    memcpy(staged_gamma, gamma.data(), sizeof(float) * n_atom_types);
    memcpy(staged_shield, shield.data(),
           sizeof(float) * pair_parameter_count);
    memcpy(staged_atom_type, atom_type.data(), sizeof(int) * atom_numbers);

    this->controller = controller;
    this->atom_numbers = atom_numbers;
    this->atom_type_numbers = n_atom_types;
    h_chi = staged_chi;
    h_eta = staged_eta;
    h_gamma = staged_gamma;
    h_shield = staged_shield;
    h_atom_type = staged_atom_type;

    Device_Malloc_And_Copy_Safely((void**)&d_chi, h_chi,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_eta, h_eta,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_gamma, h_gamma,
                                  sizeof(float) * n_atom_types);
    Device_Malloc_And_Copy_Safely((void**)&d_shield, h_shield,
                                  sizeof(float) * pair_parameter_count);
    Device_Malloc_And_Copy_Safely((void**)&d_atom_type_global, h_atom_type,
                                  sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_atom_type, sizeof(int) * atom_numbers);

    Device_Malloc_Safely((void**)&d_b, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_r, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_p, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_Ap, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_q, sizeof(float) * atom_numbers);

    Device_Malloc_Safely((void**)&d_s, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_t, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_z, sizeof(float) * atom_numbers);
    Device_Malloc_Safely((void**)&d_temp_sum, sizeof(float));
    Malloc_Safely((void**)&h_h_numnbrs, sizeof(int) * atom_numbers);
    Malloc_Safely((void**)&h_h_firstnbrs, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_h_numnbrs, sizeof(int) * atom_numbers);
    Device_Malloc_Safely((void**)&d_h_firstnbrs, sizeof(int) * atom_numbers);
    deviceMemset(d_q, 0, sizeof(float) * atom_numbers);
    deviceMemset(d_s, 0, sizeof(float) * atom_numbers);
    deviceMemset(d_t, 0, sizeof(float) * atom_numbers);

    // Device-side CG scalar buffers
    Device_Malloc_Safely((void**)&d_rr_old, sizeof(float));
    Device_Malloc_Safely((void**)&d_rr_new, sizeof(float));
    Device_Malloc_Safely((void**)&d_pAp_buf, sizeof(float));
    Device_Malloc_Safely((void**)&d_cg_alpha, sizeof(float));
    Device_Malloc_Safely((void**)&d_cg_beta, sizeof(float));
    Device_Malloc_Safely((void**)&d_solver_error, sizeof(int));

    // Charge history for extrapolation
    Device_Malloc_Safely((void**)&d_s_hist,
                         sizeof(float) * HIST_SIZE * atom_numbers);
    Device_Malloc_Safely((void**)&d_t_hist,
                         sizeof(float) * HIST_SIZE * atom_numbers);
    nprev = 0;

    is_initialized = 1;
    controller->Step_Print_Initial("REAXFF_EEQ", "%14.7e");
    controller->printf("END INITIALIZING REAXFF_EEQ\n\n");
}

// =====================================================================
// Shared kernels (CPU + GPU)
// =====================================================================

static __global__ void EEQ_Matrix_Vector_Multiply(
    int atom_numbers, const int* __restrict__ firstnbrs,
    const int* __restrict__ numnbrs, const int* __restrict__ jlist,
    const float* __restrict__ h_val, const int* __restrict__ atom_types,
    const float* __restrict__ eta, const float* __restrict__ p,
    float* __restrict__ Ap)
{
    EEQ_SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type_i = atom_types[i];
        float sum = eta[type_i] * p[i];
        int begin = firstnbrs[i];
        int end = begin + numnbrs[i];
        for (int idx = begin; idx < end; idx++)
        {
            sum += h_val[idx] * p[jlist[idx]];
        }
        Ap[i] = sum;
    }
}

static __global__ void EEQ_Count_H_Matrix_Entries(
    int atom_numbers, const VECTOR* crd, const int* atom_types,
    const float* shield, int atom_type_numbers, const ATOM_GROUP* nl,
    const LTMatrix3 cell, const LTMatrix3 rcell, float cutoff, int* numnbrs)
{
    EEQ_SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int count = 0;
        ATOM_GROUP nl_i = nl[i];
        VECTOR ri = crd[i];
        int type_i = atom_types[i];
        for (int j_idx = 0; j_idx < nl_i.atom_numbers; j_idx++)
        {
            int atom_j = nl_i.atom_serial[j_idx];
            int type_j = atom_types[atom_j];
            VECTOR rj = crd[atom_j];
            VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
            float r2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
            float r = sqrtf(r2);
            if (r < cutoff)
            {
                float shield_ij = shield[type_i * atom_type_numbers + type_j];
                if (shield_ij >= 0.0f) count++;
            }
        }
        numnbrs[i] = count;
    }
}

static __global__ void EEQ_Fill_H_Matrix(
    int atom_numbers, const VECTOR* crd, const int* atom_types,
    const float* shield, int atom_type_numbers, const ATOM_GROUP* nl,
    const LTMatrix3 cell, const LTMatrix3 rcell, float cutoff,
    const int* firstnbrs, int* jlist, float* h_val)
{
    EEQ_SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        ATOM_GROUP nl_i = nl[i];
        VECTOR ri = crd[i];
        int type_i = atom_types[i];
        int write_idx = firstnbrs[i];
        for (int j_idx = 0; j_idx < nl_i.atom_numbers; j_idx++)
        {
            int atom_j = nl_i.atom_serial[j_idx];
            int type_j = atom_types[atom_j];
            VECTOR rj = crd[atom_j];
            VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
            float r2 = drij.x * drij.x + drij.y * drij.y + drij.z * drij.z;
            float r = sqrtf(r2);
            if (r < cutoff)
            {
                float x = r / cutoff;
                float x2 = x * x;
                float x4 = x2 * x2;
                float x5 = x4 * x;
                float x6 = x5 * x;
                float x7 = x6 * x;
                float taper =
                    20.0f * x7 - 70.0f * x6 + 84.0f * x5 - 35.0f * x4 + 1.0f;
                float shield_ij = shield[type_i * atom_type_numbers + type_j];
                jlist[write_idx] = atom_j;
                h_val[write_idx] = taper * (ReaxFFEEQ::COULOMB_CONSTANT /
                                            cbrtf(r2 * r + shield_ij));
                write_idx++;
            }
        }
    }
}

static __global__ void Vector_Update_P(int n, float* p, const float* r,
                                       float beta)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { p[i] = r[i] + beta * p[i]; }
}

static __global__ void Vector_Update_X_R(int n, float* x, float* r,
                                         const float* p, const float* Ap,
                                         float alpha)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n)
    {
        x[i] += alpha * p[i];
        r[i] -= alpha * Ap[i];
    }
}

static __global__ void Vector_Subtract(int n, float* out, const float* a,
                                       const float* b)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { out[i] = a[i] - b[i]; }
}

static __global__ void Vector_Copy(int n, float* dst, const float* src)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { dst[i] = src[i]; }
}

static __global__ void Setup_B_Chi(int n, float* b, const int* atom_types,
                                   const float* chi)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { b[i] = -chi[atom_types[i]]; }
}

static __global__ void Setup_B_One(int n, float* b)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { b[i] = 1.0f; }
}

static __global__ void Vector_Scale_Add(int n, float* q, const float* t,
                                        const float* s, float mu)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { q[i] = t[i] + mu * s[i]; }
}

static __global__ void EEQ_Convert_Charge_Unit(int n, float* q_out,
                                               const float* q_in, float scale)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { q_out[i] = q_in[i] * scale; }
}

static __global__ void Elementwise_Multiply(int n, float* out, const float* a,
                                            const float* b)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n) { out[i] = a[i] * b[i]; }
}

#ifndef USE_CPU

static __device__ __forceinline__ float EEQ_Warp_Reduce_Sum(float value)
{
    for (int offset = warpSize >> 1; offset > 0; offset >>= 1)
    {
        value += deviceShflDown(FULL_MASK, value, offset);
    }
    return value;
}

static __device__ __forceinline__ float EEQ_Block_Reduce_Sum(float value)
{
    __shared__ float warp_sums[32];
    int lane = threadIdx.x & (warpSize - 1);
    int warp_id = threadIdx.x >> 5;
    int warp_count = (blockDim.x + warpSize - 1) / warpSize;

    value = EEQ_Warp_Reduce_Sum(value);
    if (lane == 0) warp_sums[warp_id] = value;
    __syncthreads();

    float block_sum = (threadIdx.x < warp_count) ? warp_sums[lane] : 0.0f;
    if (warp_id == 0)
    {
        block_sum = EEQ_Warp_Reduce_Sum(block_sum);
    }
    return block_sum;
}

static __global__ void Dot_Product_Reduce_Kernel(int n, const float* a,
                                                 const float* __restrict__ b,
                                                 float* out)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    float value = 0.0f;
    if (i < n) value = a[i] * b[i];
    float block_sum = EEQ_Block_Reduce_Sum(value);
    if (threadIdx.x == 0) atomicAdd(out, block_sum);
}

static __global__ void Initialize_Preconditioned_CG_State(
    int n, const float* __restrict__ r, float* __restrict__ z,
    float* __restrict__ p, const float* __restrict__ eta,
    const int* __restrict__ atom_types, float* rz_out)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    float value = 0.0f;
    if (i < n)
    {
        float diag = eta[atom_types[i]];
        float zi = (diag != 0.0f) ? r[i] / diag : r[i];
        z[i] = zi;
        p[i] = zi;
        value = r[i] * zi;
    }
    float block_sum = EEQ_Block_Reduce_Sum(value);
    if (threadIdx.x == 0) atomicAdd(rz_out, block_sum);
}

static __global__ void Update_X_R_Precondition_Dot_Kernel(
    int n, float* __restrict__ x, float* __restrict__ r,
    const float* __restrict__ p, const float* __restrict__ Ap,
    float* __restrict__ z, const float* __restrict__ eta,
    const int* __restrict__ atom_types, const float* __restrict__ d_alpha,
    float* rz_out)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    float value = 0.0f;
    float alpha = *d_alpha;
    if (i < n)
    {
        float ri = r[i] - alpha * Ap[i];
        x[i] += alpha * p[i];
        r[i] = ri;
        float diag = eta[atom_types[i]];
        float zi = (diag != 0.0f) ? ri / diag : ri;
        z[i] = zi;
        value = ri * zi;
    }
    float block_sum = EEQ_Block_Reduce_Sum(value);
    if (threadIdx.x == 0) atomicAdd(rz_out, block_sum);
}

#endif

// Polynomial extrapolation coefficients (oldest-to-newest order)
// nprev=k uses row k-1: Newton forward difference formula
static const float EXTRAP_COEFFS[5][5] = {
    {1.0f, 0.0f, 0.0f, 0.0f, 0.0f},     {-1.0f, 2.0f, 0.0f, 0.0f, 0.0f},
    {1.0f, -3.0f, 3.0f, 0.0f, 0.0f},    {-1.0f, 4.0f, -6.0f, 4.0f, 0.0f},
    {1.0f, -5.0f, 10.0f, -10.0f, 5.0f},
};

// Jacobi preconditioner: z[i] = r[i] / eta[type_i]
static __global__ void Jacobi_Precondition(int n, float* z, const float* r,
                                           const float* eta,
                                           const int* atom_types)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n)
    {
        float diag = eta[atom_types[i]];
        z[i] = (diag != 0.0f) ? r[i] / diag : r[i];
    }
}

static __global__ void EEQ_Distribute_Energy_Kernel(
    int n, float* d_energy, const float* d_charge, const int* atom_types,
    const float* d_chi, const float* d_eta, const float* d_Aq)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n)
    {
        int type_i = atom_types[i];
        float qi = d_charge[i];
        float e_pol_i = d_chi[type_i] * qi + 0.5f * d_eta[type_i] * qi * qi;
        float e_ele_i = 0.5f * qi * (d_Aq[i] - d_eta[type_i] * qi);
        float en_i = e_pol_i + e_ele_i;
        atomicAdd(&d_energy[i], en_i);
    }
}

static __global__ void EEQ_Calculate_Force_Kernel(
    int atom_numbers, const VECTOR* crd, const int* atom_types,
    const float* shield, int atom_type_numbers, const float* d_charge,
    VECTOR* frc, const ATOM_GROUP* nl, const LTMatrix3 cell,
    const LTMatrix3 rcell, float cutoff, LTMatrix3* atom_virial)
{
    EEQ_SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        int type_i = atom_types[i];
        float qi = d_charge[i];
        // An exactly neutral site has no pair force.  Tiny nonzero EEQ
        // charges are still part of the Hamiltonian and must not be silently
        // truncated.
        if (qi != 0.0f)
        {
            ATOM_GROUP nl_i = nl[i];
            VECTOR ri = crd[i];

            for (int j_idx = 0; j_idx < nl_i.atom_numbers; j_idx++)
            {
                int atom_j = nl_i.atom_serial[j_idx];
                if (atom_j <= i) continue;

                float qj = d_charge[atom_j];
                if (qj == 0.0f) continue;

                int type_j = atom_types[atom_j];

                VECTOR rj = crd[atom_j];
                VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
                const float shield_ij =
                    shield[type_i * atom_type_numbers + type_j];
                const ReaxFFEEQ::Pair_Force_Result pair =
                    ReaxFFEEQ::Evaluate_Pair_Force(drij, qi, qj, shield_ij,
                                                   cutoff);
                const VECTOR fij = pair.force;

                atomicAdd(&frc[i].x, fij.x);
                atomicAdd(&frc[i].y, fij.y);
                atomicAdd(&frc[i].z, fij.z);
                atomicAdd(&frc[atom_j].x, -fij.x);
                atomicAdd(&frc[atom_j].y, -fij.y);
                atomicAdd(&frc[atom_j].z, -fij.z);
                if (atom_virial)
                {
                    atomicAdd(atom_virial + i,
                              Get_Virial_From_Force_Dis(fij, drij));
                }
            }
        }
    }
}

static __global__ void EEQ_Calculate_Epol_Kernel(int n, float* out,
                                                 const int* types,
                                                 const float* chi,
                                                 const float* eta,
                                                 const float* q)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n)
    {
        int t = types[i];
        out[i] = chi[t] * q[i] + 0.5f * eta[t] * q[i] * q[i];
    }
}

static __global__ void EEQ_Calculate_Eele_Kernel(int n, float* out,
                                                 const int* types,
                                                 const float* eta,
                                                 const float* q,
                                                 const float* Aq)
{
    EEQ_SIMPLE_DEVICE_FOR(i, n)
    {
        out[i] = 0.5f * q[i] * (Aq[i] - eta[types[i]] * q[i]);
    }
}

// =====================================================================
// GPU-only kernels: device-side CG scalar operations
// =====================================================================
#ifndef USE_CPU

enum EEQ_SOLVER_ERROR
{
    EEQ_SOLVER_OK = 0,
    EEQ_SOLVER_NONFINITE = 1,
    EEQ_SOLVER_BREAKDOWN = 2
};

static __global__ void CG_Compute_Alpha_Kernel(const float* rr_old,
                                               const float* pAp, float* alpha,
                                               float tolerance_squared,
                                               int* solver_error)
{
    if (threadIdx.x == 0 && blockIdx.x == 0)
    {
        const float numerator = *rr_old;
        const float denominator = *pAp;
        if (!ReaxFF_Float_Is_Finite(numerator) ||
            !ReaxFF_Float_Is_Finite(denominator))
        {
            *solver_error = EEQ_SOLVER_NONFINITE;
            *alpha = 0.0f;
        }
        else if (fabsf(numerator) <= tolerance_squared)
        {
            *alpha = 0.0f;
        }
        else if (!(denominator > 0.0f))
        {
            *solver_error = EEQ_SOLVER_BREAKDOWN;
            *alpha = 0.0f;
        }
        else
        {
            const float value = numerator / denominator;
            if (!ReaxFF_Float_Is_Finite(value))
            {
                *solver_error = EEQ_SOLVER_NONFINITE;
                *alpha = 0.0f;
            }
            else
            {
                *alpha = value;
            }
        }
    }
}

static __global__ void CG_Compute_Beta_Kernel(float* rr_old,
                                              const float* rr_new, float* beta,
                                              float tolerance_squared,
                                              int* solver_error)
{
    if (threadIdx.x == 0 && blockIdx.x == 0)
    {
        const float old_val = *rr_old;
        const float new_val = *rr_new;
        if (!ReaxFF_Float_Is_Finite(old_val) ||
            !ReaxFF_Float_Is_Finite(new_val))
        {
            *solver_error = EEQ_SOLVER_NONFINITE;
            *beta = 0.0f;
        }
        else if (fabsf(new_val) <= tolerance_squared)
        {
            *beta = 0.0f;
        }
        else if (!(old_val > 0.0f) || !(new_val >= 0.0f))
        {
            *solver_error = EEQ_SOLVER_BREAKDOWN;
            *beta = 0.0f;
        }
        else
        {
            const float value = new_val / old_val;
            if (!ReaxFF_Float_Is_Finite(value))
            {
                *solver_error = EEQ_SOLVER_NONFINITE;
                *beta = 0.0f;
            }
            else
            {
                *beta = value;
            }
        }
        *rr_old = new_val;
    }
}

static __global__ void CG_Update_X_R_Kernel(int n, float* x, float* r,
                                            const float* p, const float* Ap,
                                            const float* d_alpha)
{
    float alpha = *d_alpha;
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n)
    {
        x[i] += alpha * p[i];
        r[i] -= alpha * Ap[i];
    }
}

static __global__ void CG_Update_P_Kernel(int n, float* p, const float* r,
                                          const float* d_beta)
{
    float beta = *d_beta;
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n)
    {
        p[i] = r[i] + beta * p[i];
    }
}

#endif  // !USE_CPU

// =====================================================================
// Calculate_Charges implementation
// =====================================================================

void REAXFF_EEQ::Calculate_Charges(int atom_numbers, const int* atom_local,
                                   float* d_charge, const VECTOR* d_crd,
                                   const LTMatrix3 cell,
                                   const LTMatrix3 rcell,
                                   const ATOM_GROUP* fnl_d_nl, float cutoff,
                                   float* d_energy, VECTOR* frc,
                                   int need_virial, LTMatrix3* atom_virial)
{
    if (!is_initialized || fnl_d_nl == NULL) return;
    if (atom_numbers != this->atom_numbers || atom_numbers <= 0 ||
        atom_local == NULL || d_charge == NULL || d_crd == NULL ||
        !(cutoff > 0.0f) || !Float_Memory_Is_Finite(&cutoff))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF_EEQ::Calculate_Charges",
            "Invalid staged EEQ evaluation: local/global atom counts are "
            "%d/%d, local-to-global map=%p, staged charge=%p, coordinates=%p, "
            "and cutoff=%g.",
            atom_numbers, this->atom_numbers,
            static_cast<const void*>(atom_local), static_cast<void*>(d_charge),
            static_cast<const void*>(d_crd), cutoff);
        return;
    }

    dim3 blockSize = {std::min(160u, CONTROLLER::device_max_thread)};
    dim3 gridSize = {(atom_numbers + blockSize.x - 1) / blockSize.x};
    const float tolerance_squared = tolerance * tolerance;
    if (!(tolerance > 0.0f) || !Float_Memory_Is_Finite(&tolerance) ||
        !(tolerance_squared > 0.0f) ||
        !Float_Memory_Is_Finite(&tolerance_squared) || max_iter <= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF_EEQ::Calculate_Charges",
            "Invalid EEQ solver controls: tolerance=%g and max_iter=%d.",
            tolerance, max_iter);
    }

    // ---- Build H matrix CSR ----
    Launch_Device_Kernel(EEQ_Count_H_Matrix_Entries, gridSize, blockSize, 0,
                         NULL, atom_numbers, d_crd, d_atom_type, d_shield,
                         atom_type_numbers, fnl_d_nl, cell, rcell, cutoff,
                         d_h_numnbrs);

    int total_nnz = 0;
#ifndef USE_CPU
    {
        thrust::device_ptr<int> d_numnbrs_ptr(d_h_numnbrs);
        thrust::device_ptr<int> d_firstnbrs_ptr(d_h_firstnbrs);
        const long long total_nnz_wide = thrust::reduce(
            d_numnbrs_ptr, d_numnbrs_ptr + atom_numbers, 0LL,
            thrust::plus<long long>());
        if (total_nnz_wide < 0 || total_nnz_wide > INT_MAX)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "REAXFF_EEQ::Calculate_Charges",
                "EEQ sparse matrix requires %lld entries, exceeding the "
                "signed-int CSR limit.",
                total_nnz_wide);
        }
        total_nnz = static_cast<int>(total_nnz_wide);
        thrust::exclusive_scan(d_numnbrs_ptr, d_numnbrs_ptr + atom_numbers,
                               d_firstnbrs_ptr);
    }
#else
    deviceMemcpy(h_h_numnbrs, d_h_numnbrs, sizeof(int) * atom_numbers,
                 deviceMemcpyDeviceToHost);
    long long total_nnz_wide = 0;
    for (int i = 0; i < atom_numbers; i++)
    {
        if (h_h_numnbrs[i] < 0 ||
            total_nnz_wide > INT_MAX - h_h_numnbrs[i])
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorOverflow, "REAXFF_EEQ::Calculate_Charges",
                "EEQ sparse-matrix row %d makes the signed-int CSR count "
                "unrepresentable (prefix=%lld, row=%d).",
                i, total_nnz_wide, h_h_numnbrs[i]);
        }
        h_h_firstnbrs[i] = static_cast<int>(total_nnz_wide);
        total_nnz_wide += h_h_numnbrs[i];
    }
    total_nnz = static_cast<int>(total_nnz_wide);
    deviceMemcpy(d_h_firstnbrs, h_h_firstnbrs, sizeof(int) * atom_numbers,
                 deviceMemcpyHostToDevice);
#endif

    if (total_nnz > h_matrix_capacity)
    {
        if (d_h_jlist != NULL) deviceFree(d_h_jlist);
        if (d_h_val != NULL) deviceFree(d_h_val);
        h_matrix_capacity = total_nnz;
        if (h_matrix_capacity > 0)
        {
            Device_Malloc_Safely((void**)&d_h_jlist,
                                 sizeof(int) * h_matrix_capacity);
            Device_Malloc_Safely((void**)&d_h_val,
                                 sizeof(float) * h_matrix_capacity);
        }
    }
    if (total_nnz > 0)
    {
        Launch_Device_Kernel(EEQ_Fill_H_Matrix, gridSize, blockSize, 0, NULL,
                             atom_numbers, d_crd, d_atom_type, d_shield,
                             atom_type_numbers, fnl_d_nl, cell, rcell, cutoff,
                             d_h_firstnbrs, d_h_jlist, d_h_val);
    }

    // ---- CG solver ----
#ifndef USE_CPU
    // GPU path: Jacobi-preconditioned CG, device-side scalars
    auto solve = [&](float* x, float* b_in, bool warm, const char* system_name)
    {
        if (!warm)
        {
            deviceMemset(x, 0, sizeof(float) * atom_numbers);
            Launch_Device_Kernel(Vector_Copy, gridSize, blockSize, 0, NULL,
                                 atom_numbers, d_r, b_in);
        }
        else
        {
            Launch_Device_Kernel(EEQ_Matrix_Vector_Multiply, gridSize,
                                 blockSize, 0, NULL, atom_numbers,
                                 d_h_firstnbrs, d_h_numnbrs, d_h_jlist, d_h_val,
                                 d_atom_type, d_eta, x, d_Ap);
            Launch_Device_Kernel(Vector_Subtract, gridSize, blockSize, 0, NULL,
                                 atom_numbers, d_r, b_in, d_Ap);
        }

        deviceMemset(d_rr_old, 0, sizeof(float));
        deviceMemset(d_solver_error, 0, sizeof(int));
        Initialize_Preconditioned_CG_State<<<gridSize, blockSize>>>(
            atom_numbers, d_r, d_z, d_p, d_eta, d_atom_type, d_rr_old);

        const int check_interval = 5;
        float h_rz = 0;
        bool converged = false;

        for (int iter = 0; iter < max_iter; iter++)
        {
            // Check convergence every check_interval iterations
            if (iter % check_interval == 0)
            {
                deviceMemcpy(&h_rz, d_rr_old, sizeof(float),
                             deviceMemcpyDeviceToHost);
                if (!Float_Memory_Is_Finite(&h_rz) || !(h_rz >= 0.0f)) break;
                if (h_rz <= tolerance_squared)
                {
                    converged = true;
                    break;
                }
            }

            Launch_Device_Kernel(EEQ_Matrix_Vector_Multiply, gridSize,
                                 blockSize, 0, NULL, atom_numbers,
                                 d_h_firstnbrs, d_h_numnbrs, d_h_jlist, d_h_val,
                                 d_atom_type, d_eta, d_p, d_Ap);

            deviceMemset(d_pAp_buf, 0, sizeof(float));
            Dot_Product_Reduce_Kernel<<<gridSize, blockSize>>>(
                atom_numbers, d_p, d_Ap, d_pAp_buf);

            // alpha = rz_old / pAp (on device)
            CG_Compute_Alpha_Kernel<<<1, 1>>>(
                d_rr_old, d_pAp_buf, d_cg_alpha, tolerance_squared,
                d_solver_error);

            deviceMemset(d_rr_new, 0, sizeof(float));
            Update_X_R_Precondition_Dot_Kernel<<<gridSize, blockSize>>>(
                atom_numbers, x, d_r, d_p, d_Ap, d_z, d_eta, d_atom_type,
                d_cg_alpha, d_rr_new);

            // beta = rz_new/rz_old, rz_old = rz_new (on device)
            CG_Compute_Beta_Kernel<<<1, 1>>>(
                d_rr_old, d_rr_new, d_cg_beta, tolerance_squared,
                d_solver_error);

            // p = z + beta*p
            CG_Update_P_Kernel<<<gridSize, blockSize>>>(atom_numbers, d_p, d_z,
                                                        d_cg_beta);
        }

        int solver_error = EEQ_SOLVER_OK;
        deviceMemcpy(&h_rz, d_rr_old, sizeof(float), deviceMemcpyDeviceToHost);
        deviceMemcpy(&solver_error, d_solver_error, sizeof(int),
                     deviceMemcpyDeviceToHost);
        converged = converged ||
                    (solver_error == EEQ_SOLVER_OK &&
                     Float_Memory_Is_Finite(&h_rz) && h_rz >= 0.0f &&
                     h_rz <= tolerance_squared);
        if (!converged || solver_error != EEQ_SOLVER_OK)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "REAXFF_EEQ::Calculate_Charges",
                "EEQ %s solve failed to converge to a finite residual: "
                "rM^-1r=%g, tolerance^2=%g, solver_error=%d, max_iter=%d. "
                "No staged result was published.",
                system_name, h_rz, tolerance_squared, solver_error, max_iter);
        }
    };
#else
    // CPU path: Jacobi-preconditioned CG with host-side scalars
    auto solve = [&](float* x, float* b_in, bool warm, const char* system_name)
    {
        if (!warm)
        {
            deviceMemset(x, 0, sizeof(float) * atom_numbers);
            deviceMemcpy(d_r, b_in, sizeof(float) * atom_numbers,
                         deviceMemcpyDeviceToDevice);
        }
        else
        {
            Launch_Device_Kernel(EEQ_Matrix_Vector_Multiply, gridSize,
                                 blockSize, 0, NULL, atom_numbers,
                                 d_h_firstnbrs, d_h_numnbrs, d_h_jlist, d_h_val,
                                 d_atom_type, d_eta, x, d_Ap);
            Launch_Device_Kernel(Vector_Subtract, gridSize, blockSize, 0, NULL,
                                 atom_numbers, d_r, b_in, d_Ap);
        }

        // z = M^{-1} * r
        Launch_Device_Kernel(Jacobi_Precondition, gridSize, blockSize, 0, NULL,
                             atom_numbers, d_z, d_r, d_eta, d_atom_type);

        deviceMemcpy(d_p, d_z, sizeof(float) * atom_numbers,
                     deviceMemcpyDeviceToDevice);

        float rz_old = 0, rz_new = 0;
        Launch_Device_Kernel(Elementwise_Multiply, gridSize, blockSize, 0, NULL,
                             atom_numbers, d_q, d_r, d_z);
        Sum_Of_List(d_q, d_temp_sum, atom_numbers);
        deviceMemcpy(&rz_old, d_temp_sum, sizeof(float),
                     deviceMemcpyDeviceToHost);
        bool converged = Float_Memory_Is_Finite(&rz_old) && rz_old >= 0.0f &&
                         rz_old <= tolerance_squared;

        for (int iter = 0; iter < max_iter && !converged; iter++)
        {
            if (!Float_Memory_Is_Finite(&rz_old) || !(rz_old >= 0.0f)) break;

            Launch_Device_Kernel(EEQ_Matrix_Vector_Multiply, gridSize,
                                 blockSize, 0, NULL, atom_numbers,
                                 d_h_firstnbrs, d_h_numnbrs, d_h_jlist, d_h_val,
                                 d_atom_type, d_eta, d_p, d_Ap);

            float p_dot_Ap = 0;
            Launch_Device_Kernel(Elementwise_Multiply, gridSize, blockSize, 0,
                                 NULL, atom_numbers, d_q, d_p, d_Ap);
            Sum_Of_List(d_q, d_temp_sum, atom_numbers);
            deviceMemcpy(&p_dot_Ap, d_temp_sum, sizeof(float),
                         deviceMemcpyDeviceToHost);

            if (!Float_Memory_Is_Finite(&p_dot_Ap) || !(p_dot_Ap > 0.0f))
                break;
            const float alpha = rz_old / p_dot_Ap;
            if (!Float_Memory_Is_Finite(&alpha)) break;
            Launch_Device_Kernel(Vector_Update_X_R, gridSize, blockSize, 0,
                                 NULL, atom_numbers, x, d_r, d_p, d_Ap, alpha);

            // z = M^{-1} * r
            Launch_Device_Kernel(Jacobi_Precondition, gridSize, blockSize, 0,
                                 NULL, atom_numbers, d_z, d_r, d_eta,
                                 d_atom_type);

            Launch_Device_Kernel(Elementwise_Multiply, gridSize, blockSize, 0,
                                 NULL, atom_numbers, d_q, d_r, d_z);
            Sum_Of_List(d_q, d_temp_sum, atom_numbers);
            deviceMemcpy(&rz_new, d_temp_sum, sizeof(float),
                         deviceMemcpyDeviceToHost);

            if (!Float_Memory_Is_Finite(&rz_new) || !(rz_new >= 0.0f)) break;
            if (rz_new <= tolerance_squared)
            {
                rz_old = rz_new;
                converged = true;
                break;
            }

            const float beta = rz_new / rz_old;
            if (!Float_Memory_Is_Finite(&beta)) break;
            Launch_Device_Kernel(Vector_Update_P, gridSize, blockSize, 0, NULL,
                                 atom_numbers, d_p, d_z, beta);
            rz_old = rz_new;
        }
        if (!converged)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "REAXFF_EEQ::Calculate_Charges",
                "EEQ %s solve failed to converge to a finite residual: "
                "rM^-1r=%g, tolerance^2=%g, max_iter=%d. No staged result "
                "was published.",
                system_name, rz_old, tolerance_squared, max_iter);
        }
    };
#endif

    Launch_Device_Kernel(Setup_B_Chi, gridSize, blockSize, 0, NULL,
                         atom_numbers, d_b, d_atom_type, d_chi);
    bool warm = nprev > 0;
    if (warm)
    {
        const float* c = EXTRAP_COEFFS[nprev - 1];
        ReaxFFAtomIdentity::Gather_Float_History_By_Global_Id(
            atom_numbers, atom_local, d_t, d_t_hist, this->atom_numbers, nprev,
            c);
    }
    solve(d_t, d_b, warm, "electronegativity");

    Launch_Device_Kernel(Setup_B_One, gridSize, blockSize, 0, NULL,
                         atom_numbers, d_b);
    if (warm)
    {
        const float* c = EXTRAP_COEFFS[nprev - 1];
        ReaxFFAtomIdentity::Gather_Float_History_By_Global_Id(
            atom_numbers, atom_local, d_s, d_s_hist, this->atom_numbers, nprev,
            c);
    }
    solve(d_s, d_b, warm, "charge-constraint");

    // d_s/d_t remain private candidate frames.  The REAXFF owner publishes
    // them in the same final kernel as force/energy/virial/charge, and only
    // when commit_sampling_state is true.

    float sum_t = 0, sum_s = 0;
    Sum_Of_List(d_t, d_temp_sum, atom_numbers);
    deviceMemcpy(&sum_t, d_temp_sum, sizeof(float), deviceMemcpyDeviceToHost);
    Sum_Of_List(d_s, d_temp_sum, atom_numbers);
    deviceMemcpy(&sum_s, d_temp_sum, sizeof(float), deviceMemcpyDeviceToHost);

    if (!Float_Memory_Is_Finite(&sum_t) ||
        !Float_Memory_Is_Finite(&sum_s) || !(sum_s > 0.0f))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF_EEQ::Calculate_Charges",
            "EEQ charge constraint is singular or non-finite "
            "(sum_t=%g, sum_s=%g). No staged result was published.",
            sum_t, sum_s);
    }

    const float Qtot = 0.0f;
    const float mu = (Qtot - sum_t) / sum_s;
    if (!Float_Memory_Is_Finite(&mu))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF_EEQ::Calculate_Charges",
            "EEQ charge-constraint multiplier is non-finite (mu=%g). No "
            "staged result was published.",
            mu);
    }

    Launch_Device_Kernel(Vector_Scale_Add, gridSize, blockSize, 0, NULL,
                         atom_numbers, d_q, d_t, d_s, mu);

    Launch_Device_Kernel(EEQ_Matrix_Vector_Multiply, gridSize, blockSize, 0,
                         NULL, atom_numbers, d_h_firstnbrs, d_h_numnbrs,
                         d_h_jlist, d_h_val, d_atom_type, d_eta, d_q, d_Ap);

    Launch_Device_Kernel(EEQ_Calculate_Epol_Kernel, gridSize, blockSize, 0,
                         NULL, atom_numbers, d_r, d_atom_type, d_chi, d_eta,
                         d_q);

    float sum_epol = 0;
    Sum_Of_List(d_r, d_temp_sum, atom_numbers);
    deviceMemcpy(&sum_epol, d_temp_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);

    Launch_Device_Kernel(EEQ_Calculate_Eele_Kernel, gridSize, blockSize, 0,
                         NULL, atom_numbers, d_r, d_atom_type, d_eta, d_q,
                         d_Ap);

    float sum_eele = 0;
    Sum_Of_List(d_r, d_temp_sum, atom_numbers);
    deviceMemcpy(&sum_eele, d_temp_sum, sizeof(float),
                 deviceMemcpyDeviceToHost);

    pending_energy = sum_epol + sum_eele;
    if (!Float_Memory_Is_Finite(&sum_epol) ||
        !Float_Memory_Is_Finite(&sum_eele) ||
        !Float_Memory_Is_Finite(&pending_energy))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "REAXFF_EEQ::Calculate_Charges",
            "EEQ produced non-finite staged energies (polarization=%g, "
            "electrostatic=%g, total=%g).",
            sum_epol, sum_eele, pending_energy);
    }
    if (d_energy != NULL)
    {
        Launch_Device_Kernel(EEQ_Distribute_Energy_Kernel, gridSize, blockSize,
                             0, NULL, atom_numbers, d_energy, d_q, d_atom_type,
                             d_chi, d_eta, d_Ap);
    }

    if (frc != NULL)
    {
        Launch_Device_Kernel(EEQ_Calculate_Force_Kernel, gridSize, blockSize, 0,
                             NULL, atom_numbers, d_crd, d_atom_type, d_shield,
                             atom_type_numbers, d_q, frc, fnl_d_nl, cell, rcell,
                             cutoff, need_virial ? atom_virial : NULL);
    }

    Launch_Device_Kernel(EEQ_Convert_Charge_Unit, gridSize, blockSize, 0, NULL,
                         atom_numbers, d_charge, d_q,
                         CONSTANT_SPONGE_CHARGE_SCALE);
}

void REAXFF_EEQ::Step_Print(CONTROLLER* controller)

{
    if (!is_initialized) return;
    controller->Step_Print("REAXFF_EEQ", h_energy, true);
}

void REAXFF_EEQ::Print_Charges(const float* d_charge)
{
    if (!is_initialized) return;
    float* h_q = NULL;
    Malloc_Safely((void**)&h_q, sizeof(float) * atom_numbers);
    deviceMemcpy(h_q, d_charge, sizeof(float) * atom_numbers,
                 deviceMemcpyDeviceToHost);

    FILE* fp = fopen("eeq_charges.txt", "w");
    if (fp)
    {
        for (int i = 0; i < atom_numbers; i++)
        {
            float q_elementary = h_q[i] / CONSTANT_SPONGE_CHARGE_SCALE;
            fprintf(fp, "%d %.6f\n", i + 1, q_elementary);
        }
        fclose(fp);
    }
    free(h_q);
}
