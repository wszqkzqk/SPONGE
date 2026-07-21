#pragma once

// Energy-weighted AO density used by the overlap (Pulay) response.
//
// For a stationary, possibly fractional zero-temperature density,
//
//     F C = S C epsilon,       P = C f C^T,
//     W = C f epsilon C^T = S^+ F P.
//
// This identity is valid for integer and ensemble occupations alike.  It also
// keeps W on the same P/F generation as the accepted SCF energy, whereas
// rebuilding W from the most recently diagonalized integer Aufbau orbitals
// silently changes an ensemble density back into a projector.  S^+ is the
// canonical-overlap pseudoinverse X X^T on the retained AO subspace.

static __global__ void QC_Symmetrize_Energy_Weighted_Density_Kernel(
    const int nao, const double* weighted_density, const int accumulate,
    float* output)
{
    const int total = nao * nao;
    SIMPLE_DEVICE_FOR(idx, total)
    {
        const int mu = idx / nao;
        const int nu = idx - mu * nao;
        const float value =
            (float)(0.5 * (weighted_density[idx] +
                          weighted_density[nu * nao + mu]));
        if (accumulate)
            output[idx] += value;
        else
            output[idx] = value;
    }
}

static inline void QC_Build_Overlap_Pseudoinverse(
    BLAS_HANDLE blas_handle, int nao, int nao_effective, const double* d_X,
    double* d_overlap_pseudoinverse)
{
    // X is row-major [nao, nao] with only the first nao_effective columns
    // populated.  Respect its padded row stride while forming X X^T.
    QC_Dgemm_NT(blas_handle, nao, nao, nao_effective, d_X, nao, d_X, nao,
                d_overlap_pseudoinverse, nao);
}

static inline void QC_Build_Energy_Weighted_Density_From_PF(
    BLAS_HANDLE blas_handle, int nao, const double* d_overlap_pseudoinverse,
    const double* d_fock, const float* d_density, float* d_weighted_density,
    double* d_density_double, double* d_fock_density,
    double* d_weighted_density_double, bool accumulate)
{
    const int nao2 = nao * nao;
    QC_Float_To_Double(nao2, d_density, d_density_double);
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_fock, nao, d_density_double, nao,
                d_fock_density, nao);
    QC_Dgemm_NN(blas_handle, nao, nao, nao, d_overlap_pseudoinverse, nao,
                d_fock_density, nao, d_weighted_density_double, nao);

    const int threads = 256;
    Launch_Device_Kernel(
        QC_Symmetrize_Energy_Weighted_Density_Kernel,
        Positive_Int_Ceil_Div(nao2, threads), threads, 0, 0, nao,
        d_weighted_density_double, accumulate ? 1 : 0, d_weighted_density);
}
