#pragma once

// Scalar pieces of the Becke atom-partition construction. Keep this header
// independent of the molecular-grid implementation so the exact production
// formulas can be exercised by a small host-side regression probe.
#if defined(__CUDACC__) || defined(__HIPCC__)
#define QC_BECKE_HOST_DEVICE __host__ __device__
#else
#define QC_BECKE_HOST_DEVICE
#endif

static inline QC_BECKE_HOST_DEVICE double QC_Becke_Shape(double mu)
{
    // Becke's three applications of f(x)=x(3-x^2)/2, followed by
    // s=(1-f)/2. Evolve delta=1-f directly: forming 1-f at the end loses the
    // small positive factor next to mu=+1.
    if (mu <= -1.0) return 1.0;
    if (mu >= 1.0) return 0.0;
    if (mu < 0.0)
    {
        // Evaluate the small complementary factor at -mu.  Besides keeping
        // s(mu)+s(-mu)=1, this avoids losing 1+mu while subtracting a value
        // close to -1 from one.
        double complement = 1.0 + mu;
        for (int iteration = 0; iteration < 3; ++iteration)
            complement =
                0.5 * complement * complement * (3.0 - complement);
        return 1.0 - 0.5 * complement;
    }
    double delta = 1.0 - mu;
    for (int iteration = 0; iteration < 3; ++iteration)
        delta = 0.5 * delta * delta * (3.0 - delta);
    return 0.5 * delta;
}

static inline QC_BECKE_HOST_DEVICE void QC_Becke_Shape_And_Derivative(
    double mu, double& shape, double& derivative)
{
    if (mu <= -1.0)
    {
        shape = 1.0;
        derivative = 0.0;
        return;
    }
    if (mu >= 1.0)
    {
        shape = 0.0;
        derivative = 0.0;
        return;
    }
    const bool use_complement = mu < 0.0;
    double delta = use_complement ? 1.0 + mu : 1.0 - mu;
    double ddelta_dmu = use_complement ? 1.0 : -1.0;
    for (int iteration = 0; iteration < 3; ++iteration)
    {
        const double old_delta = delta;
        ddelta_dmu *= 1.5 * old_delta * (2.0 - old_delta);
        delta = 0.5 * old_delta * old_delta * (3.0 - old_delta);
    }
    if (use_complement)
    {
        shape = 1.0 - 0.5 * delta;
        derivative = -0.5 * ddelta_dmu;
    }
    else
    {
        shape = 0.5 * delta;
        derivative = 0.5 * ddelta_dmu;
    }
}

static inline QC_BECKE_HOST_DEVICE double
QC_Becke_Size_Adjustment_Coefficient(double radius_a, double radius_b)
{
    // Becke's heteronuclear correction is
    //
    //   a = (r_b/r_a - r_a/r_b)/4,
    //
    // limited to [-1/2,1/2] as prescribed by the partition definition. Test
    // the limiting boundary before forming either radius ratio. Thus every
    // division below has a quotient in [1,1+sqrt(2)) and cannot overflow,
    // even for the widest representable finite radius pair.
    constexpr double ratio_limit = 2.4142135623730950488;  // 1 + sqrt(2)
    if (radius_a == radius_b) return 0.0;
    if (radius_a > radius_b)
    {
        if (radius_a / ratio_limit >= radius_b) return -0.5;
        const double ratio = radius_a / radius_b;
        return 0.25 * (1.0 / ratio - ratio);
    }
    if (radius_b / ratio_limit >= radius_a) return 0.5;
    const double ratio = radius_b / radius_a;
    return 0.25 * (ratio - 1.0 / ratio);
}

static inline QC_BECKE_HOST_DEVICE double QC_Becke_Size_Adjusted_Mu(
    double mu, double radius_a, double radius_b)
{
    const double adjustment =
        QC_Becke_Size_Adjustment_Coefficient(radius_a, radius_b);
    return mu + adjustment * (1.0 - mu * mu);
}

// Construct one directed pair factor's complete nuclear response. The grid
// point follows grid_owner, so it has a third coordinate response in addition
// to the two explicit nuclear responses. Build that response as the exact
// floating-point negative sum of the first two contributions; this preserves
// pairwise translational invariance instead of merely relying on three
// algebraically equivalent, independently rounded expressions.
static inline QC_BECKE_HOST_DEVICE void QC_Becke_Pair_Response(
    double scale, double mu, double inverse_rab, double uax, double uay,
    double uaz, double ubx, double uby, double ubz, double hx, double hy,
    double hz, double response_a[3], double response_b[3],
    double response_grid[3])
{
    const double ua[3] = {uax, uay, uaz};
    const double ub[3] = {ubx, uby, ubz};
    const double h[3] = {hx, hy, hz};
    for (int axis = 0; axis < 3; ++axis)
    {
        response_a[axis] =
            scale * (-ua[axis] - mu * h[axis]) * inverse_rab;
        response_b[axis] =
            scale * (ub[axis] + mu * h[axis]) * inverse_rab;
        response_grid[axis] = -(response_a[axis] + response_b[axis]);
    }
}

#undef QC_BECKE_HOST_DEVICE
