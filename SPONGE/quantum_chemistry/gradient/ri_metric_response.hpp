#pragma once

#include <cmath>

// Loewner coefficient for the same truncated inverse-square-root spectral
// function used by RI metric construction. Eigenvalues are ordered ascending;
// the first naux-naux_eff modes belong to the fixed discarded subspace.
static inline double QC_RI_Truncated_InvSqrt_Loewner(
    const double* eigval, const int naux, const int naux_eff, const int k,
    const int l)
{
    const int n_skip = naux - naux_eff;
    const bool k_active = (k >= n_skip);
    const bool l_active = (l >= n_skip);

    if (k_active && l_active)
    {
        const double sqrt_k = std::sqrt(eigval[k]);
        const double sqrt_l = std::sqrt(eigval[l]);
        // Stable for both distinct and exactly degenerate retained eigenvalues.
        return -1.0 / (sqrt_k * sqrt_l * (sqrt_k + sqrt_l));
    }
    if (k_active != l_active)
    {
        const double value_k = k_active ? 1.0 / std::sqrt(eigval[k]) : 0.0;
        const double value_l = l_active ? 1.0 / std::sqrt(eigval[l]) : 0.0;
        return (value_k - value_l) / (eigval[k] - eigval[l]);
    }
    return 0.0;
}
