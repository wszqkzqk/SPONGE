#pragma once

struct PME_Excluded_Radial_Kernels
{
    // erf(beta * r) / r, with its finite r -> 0 limit.
    float energy;
    // (erf(beta * r) - 2 beta r exp(-(beta r)^2) / sqrt(pi)) / r^3,
    // the scalar multiplying q_i q_j * displacement.
    float force;
};

static constexpr float PME_TWO_DIVIDED_BY_SQRT_PI = 1.1283791670218446f;

__host__ __device__ __forceinline__ PME_Excluded_Radial_Kernels
Get_PME_Excluded_Radial_Kernels(const float distance_squared, const float beta)
{
    const float beta_squared = beta * beta;
    const float beta_cubed = beta_squared * beta;
    if (distance_squared == 0.0f)
    {
        return {PME_TWO_DIVIDED_BY_SQRT_PI * beta,
                PME_TWO_DIVIDED_BY_SQRT_PI * beta_cubed * (2.0f / 3.0f)};
    }

    const float distance = sqrtf(distance_squared);
    const float beta_distance = beta * distance;
    if (fabsf(beta_distance) <= 0.25f)
    {
        // Direct evaluation of the force numerator subtracts two O(x)
        // quantities to obtain O(x^3).  These Taylor forms retain the finite
        // limit and avoid that cancellation for small x = beta * r.
        const float x2 = beta_distance * beta_distance;
        const float energy_series =
            1.0f + x2 * (-1.0f / 3.0f +
                         x2 * (1.0f / 10.0f +
                               x2 * (-1.0f / 42.0f +
                                     x2 * (1.0f / 216.0f - x2 / 1320.0f))));
        const float force_series =
            2.0f / 3.0f +
            x2 * (-2.0f / 5.0f +
                  x2 * (1.0f / 7.0f + x2 * (-1.0f / 27.0f + x2 / 132.0f)));
        return {PME_TWO_DIVIDED_BY_SQRT_PI * beta * energy_series,
                PME_TWO_DIVIDED_BY_SQRT_PI * beta_cubed * force_series};
    }

    const float erf_beta_distance = erff(beta_distance);
    const float force_numerator =
        erf_beta_distance - beta_distance * PME_TWO_DIVIDED_BY_SQRT_PI *
                                expf(-beta_distance * beta_distance);
    return {erf_beta_distance / distance,
            force_numerator / (distance_squared * distance)};
}
