#pragma once

#include <cmath>

// First-order forward automatic differentiation.  This is analytic chain
// differentiation, not numerical finite differencing.  Keeping the value and
// derivative paths together is especially useful for the five independent UKS
// variables and makes endpoint identities auditable.
template <int N>
struct QC_XC_Dual
{
    double value;
    double deriv[N];

    static inline __host__ __device__ QC_XC_Dual Constant(double v)
    {
        QC_XC_Dual result;
        result.value = v;
        for (int i = 0; i < N; ++i) result.deriv[i] = 0.0;
        return result;
    }

    static inline __host__ __device__ QC_XC_Dual Variable(double v,
                                                           int component)
    {
        QC_XC_Dual result = Constant(v);
        result.deriv[component] = 1.0;
        return result;
    }
};

static inline __host__ __device__ double QC_XC_Product(double a, double b)
{
    // Tail algebra can contain a representational zero multiplied by an
    // unbounded chain factor.  Its representable contribution is zero; doing
    // the hardware 0*Inf operation would instead poison every derivative with
    // NaN under fast-math.
    return (a == 0.0 || b == 0.0) ? 0.0 : a * b;
}

static inline __host__ __device__ double QC_XC_Coefficient_Times_Power(
    double coefficient, double base, double power)
{
    if (coefficient == 0.0) return 0.0;
    const double magnitude =
        exp(log(fabs(coefficient)) + power * log(base));
    return coefficient < 0.0 ? -magnitude : magnitude;
}

static inline __host__ __device__ double QC_XC_Exp_Difference(
    double log_positive, double log_negative)
{
    if (log_positive == log_negative) return 0.0;
    if (log_positive > log_negative)
    {
        const double log_magnitude =
            log_positive + log(-expm1(log_negative - log_positive));
        return exp(log_magnitude);
    }
    const double log_magnitude =
        log_negative + log(-expm1(log_positive - log_negative));
    return -exp(log_magnitude);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator+(
    const QC_XC_Dual<N>& a, const QC_XC_Dual<N>& b)
{
    QC_XC_Dual<N> result;
    result.value = a.value + b.value;
    for (int i = 0; i < N; ++i) result.deriv[i] = a.deriv[i] + b.deriv[i];
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator-(
    const QC_XC_Dual<N>& a, const QC_XC_Dual<N>& b)
{
    QC_XC_Dual<N> result;
    result.value = a.value - b.value;
    for (int i = 0; i < N; ++i) result.deriv[i] = a.deriv[i] - b.deriv[i];
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator-(
    const QC_XC_Dual<N>& a)
{
    QC_XC_Dual<N> result;
    result.value = -a.value;
    for (int i = 0; i < N; ++i) result.deriv[i] = -a.deriv[i];
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator*(
    const QC_XC_Dual<N>& a, const QC_XC_Dual<N>& b)
{
    QC_XC_Dual<N> result;
    result.value = a.value * b.value;
    for (int i = 0; i < N; ++i)
        result.deriv[i] = QC_XC_Product(a.deriv[i], b.value) +
                          QC_XC_Product(a.value, b.deriv[i]);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator/(
    const QC_XC_Dual<N>& a, const QC_XC_Dual<N>& b)
{
    QC_XC_Dual<N> result;
    result.value = a.value / b.value;
    for (int i = 0; i < N; ++i)
        result.deriv[i] =
            (a.deriv[i] - QC_XC_Product(result.value, b.deriv[i])) /
            b.value;
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator+(
    const QC_XC_Dual<N>& a, double b)
{
    QC_XC_Dual<N> result = a;
    result.value += b;
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator+(double a,
                                                           const QC_XC_Dual<N>& b)
{
    return b + a;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator-(
    const QC_XC_Dual<N>& a, double b)
{
    QC_XC_Dual<N> result = a;
    result.value -= b;
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator-(double a,
                                                           const QC_XC_Dual<N>& b)
{
    return QC_XC_Dual<N>::Constant(a) - b;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator*(
    const QC_XC_Dual<N>& a, double b)
{
    QC_XC_Dual<N> result;
    result.value = a.value * b;
    for (int i = 0; i < N; ++i)
        result.deriv[i] = QC_XC_Product(a.deriv[i], b);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator*(double a,
                                                           const QC_XC_Dual<N>& b)
{
    return b * a;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator/(
    const QC_XC_Dual<N>& a, double b)
{
    return a * (1.0 / b);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> operator/(double a,
                                                           const QC_XC_Dual<N>& b)
{
    return QC_XC_Dual<N>::Constant(a) / b;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Log(
    const QC_XC_Dual<N>& x)
{
    QC_XC_Dual<N> result;
    result.value = log(x.value);
    for (int i = 0; i < N; ++i) result.deriv[i] = x.deriv[i] / x.value;
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Log1p(
    const QC_XC_Dual<N>& x)
{
    QC_XC_Dual<N> result;
    result.value = log1p(x.value);
    const double scale = 1.0 / (1.0 + x.value);
    for (int i = 0; i < N; ++i) result.deriv[i] = x.deriv[i] * scale;
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Exp(
    const QC_XC_Dual<N>& x)
{
    QC_XC_Dual<N> result;
    result.value = exp(x.value);
    // Once exp underflows exactly, every representable first derivative is
    // also zero.  Avoid the indeterminate floating operation 0 * infinity.
    if (result.value == 0.0)
    {
        for (int i = 0; i < N; ++i) result.deriv[i] = 0.0;
    }
    else
    {
        for (int i = 0; i < N; ++i)
            result.deriv[i] = QC_XC_Product(result.value, x.deriv[i]);
    }
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Expm1(
    const QC_XC_Dual<N>& x)
{
    QC_XC_Dual<N> result;
    result.value = expm1(x.value);
    const double scale = exp(x.value);
    for (int i = 0; i < N; ++i)
        result.deriv[i] = QC_XC_Product(scale, x.deriv[i]);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Sqrt(
    const QC_XC_Dual<N>& x)
{
    QC_XC_Dual<N> result;
    result.value = sqrt(x.value);
    const double scale = 0.5 / result.value;
    for (int i = 0; i < N; ++i)
        result.deriv[i] = QC_XC_Product(scale, x.deriv[i]);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Cbrt(
    const QC_XC_Dual<N>& x)
{
    QC_XC_Dual<N> result;
    result.value = cbrt(x.value);
    const double scale = 1.0 / (3.0 * result.value * result.value);
    for (int i = 0; i < N; ++i)
        result.deriv[i] = QC_XC_Product(scale, x.deriv[i]);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Atan(
    const QC_XC_Dual<N>& x)
{
    QC_XC_Dual<N> result;
    result.value = atan(x.value);
    const double scale = 1.0 / (1.0 + x.value * x.value);
    for (int i = 0; i < N; ++i) result.deriv[i] = scale * x.deriv[i];
    return result;
}

// x^(4/3) is continuously differentiable at x=0 even though a naive cbrt
// chain is not.  This endpoint definition is used by spin interpolation.
template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PowFourThirds(
    const QC_XC_Dual<N>& x)
{
    if (x.value == 0.0) return QC_XC_Dual<N>::Constant(0.0);
    const QC_XC_Dual<N> root = QC_XC_Cbrt(x);
    return x * root;
}

// Used only to evaluate phi at a fully polarized endpoint.  Its missing-spin
// derivative is overwritten with the exact +infinity below; setting the local
// derivative of x^(2/3) to zero here lets all other finite partials be obtained
// without producing an intermediate NaN.
template <int N>
static inline __host__ __device__ QC_XC_Dual<N>
QC_XC_PowTwoThirds_EndpointAware(const QC_XC_Dual<N>& x)
{
    if (x.value == 0.0) return QC_XC_Dual<N>::Constant(0.0);
    const QC_XC_Dual<N> root = QC_XC_Cbrt(x);
    return root * root;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_LogisticMinus(
    const QC_XC_Dual<N>& log_ratio)
{
    QC_XC_Dual<N> result;
    if (log_ratio.value >= 0.0)
    {
        const double t = exp(-log_ratio.value);
        result.value = t / (1.0 + t);
    }
    else
    {
        const double t = exp(log_ratio.value);
        result.value = 1.0 / (1.0 + t);
    }
    const double factor = -result.value * (1.0 - result.value);
    if (factor == 0.0)
    {
        for (int i = 0; i < N; ++i) result.deriv[i] = 0.0;
    }
    else
    {
        for (int i = 0; i < N; ++i)
            result.deriv[i] = factor * log_ratio.deriv[i];
    }
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_LogAddExp(
    const QC_XC_Dual<N>& a, const QC_XC_Dual<N>& b)
{
    if (a.value >= b.value) return a + QC_XC_Log1p(QC_XC_Exp(b - a));
    return b + QC_XC_Log1p(QC_XC_Exp(a - b));
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Pow43_Density(
    const QC_XC_Dual<N>& rho)
{
    return rho * QC_XC_Cbrt(rho);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Slater_Impl(
    const QC_XC_Dual<N>& rho)
{
    const double cx = 0.75 * cbrt(3.0 / CONSTANT_Pi);
    return -cx * QC_XC_Pow43_Density(rho);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PW92_Impl(
    const QC_XC_Dual<N>& rho, const double parameter[6])
{
    const double rs_factor = cbrt(3.0 / (4.0 * CONSTANT_Pi));
    const QC_XC_Dual<N> x = QC_XC_Cbrt(rho);
    const QC_XC_Dual<N> inverse_sqrt_rs = QC_XC_Sqrt(x / rs_factor);
    if (inverse_sqrt_rs.value <= 1.0)
    {
        // Exact low-density rearrangement.  Forming rs first makes drs/drho
        // overflow around rho=1e-235 although epsilon and its final density
        // derivative remain representable.  In y=1/sqrt(rs), Q is a regular
        // polynomial and the large (1+p1*rs) factor is combined with log1p.
        const QC_XC_Dual<N> y = inverse_sqrt_rs;
        const QC_XC_Dual<N> y2 = y * y;
        const QC_XC_Dual<N> y3 = y2 * y;
        const QC_XC_Dual<N> y4 = y2 * y2;
        const QC_XC_Dual<N> denominator =
            parameter[5] + parameter[4] * y + parameter[3] * y2 +
            parameter[2] * y3;
        const QC_XC_Dual<N> u =
            0.5 * y4 / (parameter[0] * denominator);
        const QC_XC_Dual<N> logarithm = QC_XC_Log1p(u);
        return -2.0 * parameter[0] *
               (logarithm + parameter[1] * logarithm / y2);
    }

    const QC_XC_Dual<N> sqrt_rs = 1.0 / inverse_sqrt_rs;
    const QC_XC_Dual<N> rs = sqrt_rs * sqrt_rs;
    const QC_XC_Dual<N> polynomial =
        sqrt_rs *
        (parameter[2] +
         sqrt_rs *
             (parameter[3] +
              sqrt_rs * (parameter[4] + parameter[5] * sqrt_rs)));
    const QC_XC_Dual<N> logarithm =
        QC_XC_Log1p(0.5 / (parameter[0] * polynomial));
    return -2.0 * parameter[0] * (1.0 + parameter[1] * rs) * logarithm;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PW92_Unpolarized_Impl(
    const QC_XC_Dual<N>& rho)
{
    const double p[6] = {0.03109070, 0.21370, 7.59570,
                         3.5876,     1.63820, 0.49294};
    return QC_XC_PW92_Impl(rho, p);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PW92_Polarized_Impl(
    const QC_XC_Dual<N>& rho)
{
    const double p[6] = {0.01554535, 0.20548, 14.11890,
                         6.1977,     3.36620, 0.62517};
    return QC_XC_PW92_Impl(rho, p);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PW92_Alpha_Impl(
    const QC_XC_Dual<N>& rho)
{
    const double p[6] = {0.01688690, 0.11125, 10.35700,
                         3.6231,     0.88026, 0.49671};
    return QC_XC_PW92_Impl(rho, p);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PW92_Spin_Impl(
    const QC_XC_Dual<N>& rho, const QC_XC_Dual<N>& zeta)
{
    const QC_XC_Dual<N> plus = 1.0 + zeta;
    const QC_XC_Dual<N> minus = 1.0 - zeta;
    const double denominator = cbrt(2.0) * 2.0 - 2.0;
    const QC_XC_Dual<N> fz =
        (QC_XC_PowFourThirds(plus) + QC_XC_PowFourThirds(minus) - 2.0) /
        denominator;
    const QC_XC_Dual<N> z2 = zeta * zeta;
    const QC_XC_Dual<N> z4 = z2 * z2;
    const QC_XC_Dual<N> ec0 = QC_XC_PW92_Unpolarized_Impl(rho);
    const QC_XC_Dual<N> ec1 = QC_XC_PW92_Polarized_Impl(rho);
    const QC_XC_Dual<N> ec2 =
        QC_XC_PW92_Alpha_Impl(rho) / 1.70992093416136561756;
    return ec0 + fz * (z4 * (ec1 - ec0) - (1.0 - z4) * ec2);
}

// Canonical VWN form.  In the density tail its three O(1/sqrt(rs))
// contributions cancel analytically, leaving O(1/rs).  Evaluating those terms
// separately eventually loses every significant bit even if each logarithm
// uses log1p.  For y=1/sqrt(rs) <= 1/64 we therefore use the Taylor series of
// the *combined* canonical expression.  Sixteen terms put the first omitted
// contribution below double roundoff at the switch (the closest logarithmic
// singularity is more than sixteen times farther away); this is a numerical
// representation branch, never a density cutoff.
template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_VWN5_Canonical_Impl(
    const QC_XC_Dual<N>& rho, double amplitude, double x0, double b, double c,
    const double tail_coefficient[15])
{
    const double rs_factor = cbrt(3.0 / (4.0 * CONSTANT_Pi));
    const QC_XC_Dual<N> density_root = QC_XC_Cbrt(rho);
    // Form y directly.  At rho~1e-300, drs/drho is larger than DBL_MAX even
    // though dy/drho and the final VWN potential are representable.
    const QC_XC_Dual<N> y = QC_XC_Sqrt(density_root / rs_factor);
    if (y.value <= 1.0 / 64.0)
    {
        QC_XC_Dual<N> polynomial =
            QC_XC_Dual<N>::Constant(tail_coefficient[14]);
        for (int degree = 15; degree >= 2; --degree)
            polynomial = tail_coefficient[degree - 2] + y * polynomial;
        return y * y * polynomial;
    }

    const QC_XC_Dual<N> s = 1.0 / y;
    const QC_XC_Dual<N> rs = s * s;
    const double x0_polynomial = x0 * x0 + b * x0 + c;
    const double q = sqrt(4.0 * c - b * b);
    const QC_XC_Dual<N> denominator = rs + b * s + c;
    const QC_XC_Dual<N> log_rs_ratio =
        -QC_XC_Log1p((b * s + c) / rs);
    const QC_XC_Dual<N> difference =
        -(2.0 * x0 + b) * s + x0 * x0 - c;
    const QC_XC_Dual<N> log_shifted_ratio =
        QC_XC_Log1p(difference / denominator);
    const QC_XC_Dual<N> angle = QC_XC_Atan(q / (2.0 * s + b));
    return amplitude *
           (log_rs_ratio + 2.0 * b * angle / q -
            (b * x0 / x0_polynomial) *
                (log_shifted_ratio + 2.0 * (b + 2.0 * x0) * angle / q));
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_VWN5_Impl(
    const QC_XC_Dual<N>& rho)
{
    const double tail[15] = {
        -4.14330420340463840e-1,  1.03044597895496882e+0,
        -2.01037782769109637e-1, -7.39792424753635025e+0,
        2.47130744235946514e+1,  -1.06044091189988243e+1,
        -2.05165033831426540e+2, 7.86457099943648298e+2,
        -5.15243890859578198e+2, -6.57740109338758327e+3,
        2.80277812514475752e+4,  -2.44448069630531737e+4,
        -2.26144593290141294e+5, 1.06078304738386467e+6,
        -1.14730748927092429e+6};
    return QC_XC_VWN5_Canonical_Impl(rho, 0.0310907, -0.10498, 3.72744,
                                     12.9352, tail);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_VWN5_Polarized_Impl(
    const QC_XC_Dual<N>& rho)
{
    const double tail[15] = {
        -3.16385748745275000e-1, 1.49693954639905259e+0,
        -5.07203501815328620e+0, 1.24301868281404169e+1,
        -1.20755530929833771e+1, -8.72511457763403641e+1,
        7.02569449302436744e+2,  -3.18383743569203545e+3,
        1.00818196747590545e+4,  -1.76709923510527734e+4,
        -3.73453272718465206e+4, 5.13398157643526203e+5,
        -2.78785661805350475e+6, 1.03364802974056079e+7,
        -2.43689113747697158e+7};
    return QC_XC_VWN5_Canonical_Impl(rho, 0.01554535, -0.32500, 7.06042,
                                     18.0578, tail);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_VWN5_Alpha_Impl(
    const QC_XC_Dual<N>& rho)
{
    const double tail[15] = {
        -1.31817657838303637e+0, 9.93968384905211222e-1,
        7.72792778478138764e+0,  -1.47482989405025382e+1,
        -5.30974258492646314e+1, 1.88473242941541827e+2,
        3.31349728810417901e+2,  -2.23947110188419538e+3,
        -1.16753131736472255e+3, 2.50285830883085764e+4,
        -1.32973553369941383e+4, -2.61526416929140774e+5,
        4.22897599064846116e+5,  2.50111258221168583e+6,
        -7.46425041899847891e+6};
    const double inverse_pi2 = 1.0 / (CONSTANT_Pi * CONSTANT_Pi);
    return QC_XC_VWN5_Canonical_Impl(rho, inverse_pi2, -0.0047584, 1.13107,
                                     13.0045, tail);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_VWN5_Spin_Impl(
    const QC_XC_Dual<N>& rho_a, const QC_XC_Dual<N>& rho_b)
{
    const QC_XC_Dual<N> rho = rho_a + rho_b;
    const QC_XC_Dual<N> zeta = (rho_a - rho_b) / rho;
    const QC_XC_Dual<N> plus = 1.0 + zeta;
    const QC_XC_Dual<N> minus = 1.0 - zeta;
    const QC_XC_Dual<N> z2 = zeta * zeta;
    const QC_XC_Dual<N> z4 = z2 * z2;
    const QC_XC_Dual<N> fz_numerator =
        QC_XC_PowFourThirds(plus) + QC_XC_PowFourThirds(minus) - 2.0;
    const QC_XC_Dual<N> eps0 = QC_XC_VWN5_Impl(rho);
    const QC_XC_Dual<N> eps1 = QC_XC_VWN5_Polarized_Impl(rho);
    const QC_XC_Dual<N> eps_alpha = QC_XC_VWN5_Alpha_Impl(rho);
    const double cbrt2 = cbrt(2.0);
    const QC_XC_Dual<N> eps =
        eps0 - eps_alpha * fz_numerator * (1.0 - z4) * (3.0 / 16.0) +
        (eps1 - eps0) * fz_numerator * z4 / (2.0 * (cbrt2 - 1.0));
    return rho * eps;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PBE_X_Impl(
    const QC_XC_Dual<N>& rho, const QC_XC_Dual<N>& sigma)
{
    const double cx = 0.75 * cbrt(3.0 / CONSTANT_Pi);
    const double kappa = 0.804;
    const double mu = 0.2195149727645171;
    const double kf_factor = cbrt(3.0 * CONSTANT_Pi * CONSTANT_Pi);
    const double reduced_factor = mu / (4.0 * kappa * kf_factor * kf_factor);
    const double rho13 = cbrt(rho.value);
    const double rho43 = rho.value * rho13;

    double inverse_enhancement = 1.0;
    if (sigma.value > 0.0)
    {
        const double log_ratio =
            log(reduced_factor) + log(sigma.value) -
            (8.0 / 3.0) * log(rho.value);
        if (log_ratio >= 0.0)
        {
            const double inverse_ratio = exp(-log_ratio);
            inverse_enhancement = inverse_ratio / (1.0 + inverse_ratio);
        }
        else
        {
            const double ratio = exp(log_ratio);
            inverse_enhancement = 1.0 / (1.0 + ratio);
        }
    }

    const double enhancement =
        1.0 + kappa - kappa * inverse_enhancement;
    const double vrho =
        -cx * rho13 *
        ((4.0 / 3.0) * enhancement -
         (8.0 / 3.0) * kappa * inverse_enhancement *
             (1.0 - inverse_enhancement));
    double vsigma = 0.0;
    if (inverse_enhancement != 0.0)
    {
        const double log_magnitude =
            log(cx * kappa * reduced_factor) -
            (4.0 / 3.0) * log(rho.value) +
            2.0 * log(inverse_enhancement);
        vsigma = -exp(log_magnitude);
    }

    QC_XC_Dual<N> result = QC_XC_Dual<N>::Constant(-cx * rho43 * enhancement);
    for (int i = 0; i < N; ++i)
        result.deriv[i] = QC_XC_Product(vrho, rho.deriv[i]) +
                          QC_XC_Product(vsigma, sigma.deriv[i]);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_AsinhExp(
    const QC_XC_Dual<N>& log_value)
{
    QC_XC_Dual<N> result;
    double scale;
    if (log_value.value >= 0.0)
    {
        const double inverse_square = exp(-2.0 * log_value.value);
        result.value =
            log_value.value + log1p(sqrt(1.0 + inverse_square));
        scale = 1.0 / sqrt(1.0 + inverse_square);
    }
    else
    {
        const double value = exp(log_value.value);
        result.value = asinh(value);
        scale = value / sqrt(1.0 + value * value);
    }
    for (int i = 0; i < N; ++i)
        result.deriv[i] = scale * log_value.deriv[i];
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_B88_Impl(
    const QC_XC_Dual<N>& rho, const QC_XC_Dual<N>& sigma)
{
    const double beta = 0.0042;
    const double cbrt2 = cbrt(2.0);
    const double cx = 0.75 * cbrt(3.0 / CONSTANT_Pi);
    const QC_XC_Dual<N> rho43 = QC_XC_Pow43_Density(rho);
    QC_XC_Dual<N> result = -cx * rho43;
    if (sigma.value == 0.0)
    {
        const double slope =
            -exp(log(beta * cbrt2) - (4.0 / 3.0) * log(rho.value));
        for (int i = 0; i < N; ++i)
            result.deriv[i] += slope * sigma.deriv[i];
        return result;
    }

    const QC_XC_Dual<N> log_q =
        log(cbrt2) + 0.5 * QC_XC_Log(sigma) -
        (4.0 / 3.0) * QC_XC_Log(rho);
    const QC_XC_Dual<N> asinh_q = QC_XC_AsinhExp(log_q);
    const QC_XC_Dual<N> log_denominator =
        QC_XC_LogAddExp(-log_q, log(6.0 * beta) + QC_XC_Log(asinh_q));
    const QC_XC_Dual<N> log_correction =
        log(beta) + 0.5 * QC_XC_Log(sigma) - log_denominator;
    return result - QC_XC_Exp(log_correction);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PBE_C_Impl(
    const QC_XC_Dual<N>& rho, const QC_XC_Dual<N>& sigma)
{
    const double gamma = (1.0 - log(2.0)) / (CONSTANT_Pi * CONSTANT_Pi);
    const double beta = 0.06672455060314922;
    const double beta_gamma = beta / gamma;
    const double d2 =
        (1.0 / 144.0) * pow(3.0, 5.0 / 3.0) *
        pow(CONSTANT_Pi, 1.0 / 3.0);
    const QC_XC_Dual<N> eps = QC_XC_PW92_Unpolarized_Impl(rho);

    QC_XC_Dual<N> y;
    if (sigma.value == 0.0)
    {
        QC_XC_Dual<N> result = rho * eps;
        const double vsigma =
            exp(log(beta * d2) - (4.0 / 3.0) * log(rho.value));
        for (int i = 0; i < N; ++i)
            result.deriv[i] += QC_XC_Product(vsigma, sigma.deriv[i]);
        return result;
    }
    else
    {
        // y = numerator/(tail_metric+numerator), evaluated as a logistic so
        // two positive, representable inputs cannot both disappear through
        // intermediate product underflow.
        const QC_XC_Dual<N> log_tail_to_numerator =
            QC_XC_Log(QC_XC_Expm1(-eps / gamma)) +
            (7.0 / 3.0) * QC_XC_Log(rho) -
            log(beta_gamma * d2) - QC_XC_Log(sigma);
        y = QC_XC_LogisticMinus(log_tail_to_numerator);
    }
    const QC_XC_Dual<N> r = y / (1.0 - y + y * y);
    const QC_XC_Dual<N> eps_plus_h =
        gamma * QC_XC_Log1p((1.0 - r) * QC_XC_Expm1(eps / gamma));
    QC_XC_Dual<N> result = rho * eps_plus_h;

    // Recover dE/dsigma directly.  In the small-gradient limit y may round to
    // zero although y/sigma remains finite and physically important.
    const double y_value = y.value;
    const double denominator = 1.0 - y_value + y_value * y_value;
    double dy_dsigma;
    if (y_value == 0.0)
    {
        const double log_tail =
            log(expm1(-eps.value / gamma)) +
            (7.0 / 3.0) * log(rho.value);
        dy_dsigma = exp(log(beta_gamma * d2) - log_tail);
    }
    else if (y_value == 1.0)
    {
        dy_dsigma = 0.0;
    }
    else
    {
        dy_dsigma =
            exp(log(y_value) + log1p(-y_value) - log(sigma.value));
    }
    const double dr_dy =
        (1.0 - y_value * y_value) / (denominator * denominator);
    const double m = expm1(eps.value / gamma);
    const double vsigma =
        rho.value * gamma * m / (1.0 + (1.0 - r.value) * m) *
        (-dr_dy * dy_dsigma);
    for (int i = 0; i < N; ++i)
        if (sigma.deriv[i] != 0.0 && rho.deriv[i] == 0.0)
            result.deriv[i] = QC_XC_Product(vsigma, sigma.deriv[i]);
    return result;
}

// LYP in a tail-stable algebra.  x=cbrt(rho) is used to combine every
// exp(-c/x)*rho^(-n/3) product before evaluation, so no infinity-times-zero
// intermediate is created at positive low density.
template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_LYP_Impl(
    const QC_XC_Dual<N>& rho, const QC_XC_Dual<N>& sigma)
{
    const double a = 0.04918;
    const double b = 0.132;
    const double c = 0.2533;
    const double d = 0.349;
    const double constant_term =
        0.3 * pow(3.0, 2.0 / 3.0) * pow(CONSTANT_Pi, 4.0 / 3.0);
    const QC_XC_Dual<N> x = QC_XC_Cbrt(rho);
    const QC_XC_Dual<N> denominator = x + d;
    const QC_XC_Dual<N> x2 = x * x;
    const QC_XC_Dual<N> x4 = x2 * x2;
    const QC_XC_Dual<N> exponential = QC_XC_Exp(-c / x);
    const QC_XC_Dual<N> local = -a * x4 / denominator;
    const QC_XC_Dual<N> constant =
        -a * b * constant_term * exponential * x4 / denominator;

    const QC_XC_Dual<N> numerator =
        3.0 * x * denominator + 7.0 * d * x + 7.0 * c * denominator;
    const QC_XC_Dual<N> scaled_exponential =
        QC_XC_Exp(-c / x - 5.0 * QC_XC_Log(x));
    const QC_XC_Dual<N> denominator_squared = denominator * denominator;
    const QC_XC_Dual<N> gradient =
        (a * b / 72.0) * sigma * scaled_exponential * numerator /
        denominator_squared;
    return local + constant + gradient;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_PBE_C_Spin_Impl(
    const QC_XC_Dual<N>& rho_a, const QC_XC_Dual<N>& rho_b,
    const QC_XC_Dual<N>& sigma_total)
{
    const double gamma = (1.0 - log(2.0)) / (CONSTANT_Pi * CONSTANT_Pi);
    const double beta = 0.06672455060314922;
    const double beta_gamma = beta / gamma;
    const QC_XC_Dual<N> rho = rho_a + rho_b;
    const QC_XC_Dual<N> zeta = (rho_a - rho_b) / rho;
    const QC_XC_Dual<N> eps = QC_XC_PW92_Spin_Impl(rho, zeta);

    const QC_XC_Dual<N> phi =
        0.5 * (QC_XC_PowTwoThirds_EndpointAware(1.0 + zeta) +
               QC_XC_PowTwoThirds_EndpointAware(1.0 - zeta));
    const QC_XC_Dual<N> phi2 = phi * phi;
    const QC_XC_Dual<N> phi3 = phi2 * phi;
    const QC_XC_Dual<N> kf =
        cbrt(3.0 * CONSTANT_Pi * CONSTANT_Pi) * QC_XC_Cbrt(rho);
    const QC_XC_Dual<N> ks2 = (4.0 / CONSTANT_Pi) * kf;
    if (sigma_total.value == 0.0)
    {
        QC_XC_Dual<N> result = rho * eps;
        const double vsigma =
            beta * CONSTANT_Pi * phi.value /
            (16.0 * kf.value * rho.value);
        for (int i = 0; i < N; ++i)
            result.deriv[i] += QC_XC_Product(vsigma, sigma_total.deriv[i]);
        return result;
    }
    const QC_XC_Dual<N> em1 = QC_XC_Expm1(-eps / (gamma * phi3));
    const QC_XC_Dual<N> log_tail_to_numerator =
        QC_XC_Log(em1) + log(4.0) + 2.0 * QC_XC_Log(phi) +
        QC_XC_Log(ks2) + 2.0 * QC_XC_Log(rho) - log(beta_gamma) -
        QC_XC_Log(sigma_total);
    const QC_XC_Dual<N> y =
        QC_XC_LogisticMinus(log_tail_to_numerator);
    const QC_XC_Dual<N> r = y / (1.0 - y + y * y);
    const QC_XC_Dual<N> eps_plus_h =
        gamma * phi3 *
        QC_XC_Log1p((1.0 - r) *
                    QC_XC_Expm1(eps / (gamma * phi3)));
    QC_XC_Dual<N> result = rho * eps_plus_h;

    const double y_value = y.value;
    const double denominator = 1.0 - y_value + y_value * y_value;
    double dy_dsigma;
    if (y_value == 0.0)
    {
        const double log_tail = log(em1.value) + log(4.0) +
                                2.0 * log(phi.value) + log(ks2.value) +
                                2.0 * log(rho.value);
        dy_dsigma = exp(log(beta_gamma) - log_tail);
    }
    else if (y_value == 1.0)
    {
        dy_dsigma = 0.0;
    }
    else
    {
        dy_dsigma = exp(log(y_value) + log1p(-y_value) -
                        log(sigma_total.value));
    }
    const double dr_dy =
        (1.0 - y_value * y_value) / (denominator * denominator);
    const double m = expm1(eps.value / (gamma * phi3.value));
    const double vsigma =
        rho.value * gamma * phi3.value * m /
        (1.0 + (1.0 - r.value) * m) * (-dr_dy * dy_dsigma);
    for (int i = 0; i < N; ++i)
        if (sigma_total.deriv[i] != 0.0 && rho.deriv[i] == 0.0)
            result.deriv[i] =
                QC_XC_Product(vsigma, sigma_total.deriv[i]);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_LYP_Spin_Impl(
    const QC_XC_Dual<N>& rho_a, const QC_XC_Dual<N>& rho_b,
    const QC_XC_Dual<N>& sigma_aa, const QC_XC_Dual<N>& sigma_ab,
    const QC_XC_Dual<N>& sigma_bb)
{
    const double a_lyp = 0.04918;
    const double b_lyp = 0.132;
    const double c = 0.2533;
    const double d = 0.349;
    const double cbrt2 = cbrt(2.0);
    const double scale83 = 4.0 * cbrt2 * cbrt2;
    const double gfac = pow(3.0 * CONSTANT_Pi * CONSTANT_Pi, 2.0 / 3.0);

    const QC_XC_Dual<N> rho = rho_a + rho_b;
    const QC_XC_Dual<N> plus = 2.0 * rho_a / rho;
    const QC_XC_Dual<N> minus = 2.0 * rho_b / rho;
    const QC_XC_Dual<N> one_minus_z2 = 4.0 * rho_a * rho_b / (rho * rho);
    const QC_XC_Dual<N> plus2 = plus * plus;
    const QC_XC_Dual<N> minus2 = minus * minus;
    const QC_XC_Dual<N> spin_power_sum =
        QC_XC_PowFourThirds(plus * plus) +
        QC_XC_PowFourThirds(minus * minus);
    const QC_XC_Dual<N> sigma_total =
        sigma_aa + 2.0 * sigma_ab + sigma_bb;
    const QC_XC_Dual<N> sigma_sum = sigma_aa + sigma_bb;
    const QC_XC_Dual<N> sigma_weighted =
        sigma_aa * plus + sigma_bb * minus;

    const QC_XC_Dual<N> t92_constant =
        -(3.0 / 20.0) * gfac * one_minus_z2 * spin_power_sum;

    QC_XC_Dual<N> gradient0 =
        -sigma_total * (47.0 * one_minus_z2 / 72.0 - 2.0 / 3.0);
    QC_XC_Dual<N> gradient1 =
        sigma_total * one_minus_z2 * (7.0 / 216.0);

    gradient0 =
        gradient0 + cbrt2 * one_minus_z2 * scale83 * 2.5 * sigma_sum / 32.0;
    gradient1 =
        gradient1 - cbrt2 * one_minus_z2 * scale83 * sigma_sum /
                        (32.0 * 54.0);
    gradient0 =
        gradient0 - cbrt2 * one_minus_z2 * scale83 * 11.0 *
                        sigma_weighted / 576.0;
    gradient1 =
        gradient1 + cbrt2 * one_minus_z2 * scale83 * sigma_weighted /
                        (576.0 * 3.0);
    gradient0 =
        gradient0 -
        (cbrt2 * scale83 / 8.0) *
            ((2.0 / 3.0) * sigma_sum - 0.25 * plus2 * sigma_bb -
             0.25 * minus2 * sigma_aa);

    const QC_XC_Dual<N> x = QC_XC_Cbrt(rho);
    const QC_XC_Dual<N> denominator = x + d;
    const QC_XC_Dual<N> x2 = x * x;
    const QC_XC_Dual<N> x4 = x2 * x2;
    const QC_XC_Dual<N> exponential = QC_XC_Exp(-c / x);
    const QC_XC_Dual<N> exp_x4 =
        QC_XC_Exp(-c / x - 4.0 * QC_XC_Log(x));
    const QC_XC_Dual<N> exp_x5 =
        QC_XC_Exp(-c / x - 5.0 * QC_XC_Log(x));

    const QC_XC_Dual<N> local =
        -a_lyp * one_minus_z2 * x4 / denominator;
    const QC_XC_Dual<N> constant =
        a_lyp * b_lyp * exponential * x4 * t92_constant / denominator;
    const QC_XC_Dual<N> gradient =
        a_lyp * b_lyp *
        (exp_x4 * (gradient0 + 3.0 * d * gradient1 / denominator) /
             denominator +
         exp_x5 * (3.0 * c * gradient1) / denominator);
    return local + constant + gradient;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_RKS_Energy_Impl(
    QC_METHOD method, const QC_XC_Dual<N>& rho, const QC_XC_Dual<N>& sigma)
{
    switch (method)
    {
        case QC_METHOD::LDA:
            return QC_XC_Slater_Impl(rho) + rho * QC_XC_VWN5_Impl(rho);
        case QC_METHOD::PBE:
            return QC_XC_PBE_X_Impl(rho, sigma) + QC_XC_PBE_C_Impl(rho, sigma);
        case QC_METHOD::BLYP:
            return QC_XC_B88_Impl(rho, sigma) + QC_XC_LYP_Impl(rho, sigma);
        case QC_METHOD::PBE0:
            return 0.75 * QC_XC_PBE_X_Impl(rho, sigma) +
                   QC_XC_PBE_C_Impl(rho, sigma);
        case QC_METHOD::B3LYP:
            return 0.08 * QC_XC_Slater_Impl(rho) +
                   0.72 * QC_XC_B88_Impl(rho, sigma) +
                   0.81 * QC_XC_LYP_Impl(rho, sigma) +
                   0.19 * rho * QC_XC_VWN5_Impl(rho);
        default:
            return QC_XC_Dual<N>::Constant(QC_XC_Quiet_NaN());
    }
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_GGA_X_Spin_Channel(
    QC_METHOD method, const QC_XC_Dual<N>& rho_spin,
    const QC_XC_Dual<N>& sigma_spin)
{
    if (rho_spin.value == 0.0) return QC_XC_Dual<N>::Constant(0.0);
    const QC_XC_Dual<N> doubled_rho = 2.0 * rho_spin;
    const QC_XC_Dual<N> quadrupled_sigma = 4.0 * sigma_spin;
    if (method == QC_METHOD::PBE || method == QC_METHOD::PBE0)
        return 0.5 * QC_XC_PBE_X_Impl(doubled_rho, quadrupled_sigma);
    return 0.5 * QC_XC_B88_Impl(doubled_rho, quadrupled_sigma);
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_Slater_Spin_Impl(
    const QC_XC_Dual<N>& rho_a, const QC_XC_Dual<N>& rho_b)
{
    QC_XC_Dual<N> result = QC_XC_Dual<N>::Constant(0.0);
    if (rho_a.value != 0.0) result = result + 0.5 * QC_XC_Slater_Impl(2.0 * rho_a);
    if (rho_b.value != 0.0) result = result + 0.5 * QC_XC_Slater_Impl(2.0 * rho_b);
    return result;
}

template <int N>
static inline __host__ __device__ QC_XC_Dual<N> QC_XC_UKS_Energy_Impl(
    QC_METHOD method, const QC_XC_Dual<N>& rho_a,
    const QC_XC_Dual<N>& rho_b, const QC_XC_Dual<N>& sigma_aa,
    const QC_XC_Dual<N>& sigma_ab, const QC_XC_Dual<N>& sigma_bb)
{
    const QC_XC_Dual<N> rho = rho_a + rho_b;
    const QC_XC_Dual<N> sigma_total =
        sigma_aa + 2.0 * sigma_ab + sigma_bb;
    const QC_XC_Dual<N> slater = QC_XC_Slater_Spin_Impl(rho_a, rho_b);
    switch (method)
    {
        case QC_METHOD::LDA:
            return slater + QC_XC_VWN5_Spin_Impl(rho_a, rho_b);
        case QC_METHOD::PBE:
            return QC_XC_GGA_X_Spin_Channel(method, rho_a, sigma_aa) +
                   QC_XC_GGA_X_Spin_Channel(method, rho_b, sigma_bb) +
                   QC_XC_PBE_C_Spin_Impl(rho_a, rho_b, sigma_total);
        case QC_METHOD::BLYP:
            return QC_XC_GGA_X_Spin_Channel(method, rho_a, sigma_aa) +
                   QC_XC_GGA_X_Spin_Channel(method, rho_b, sigma_bb) +
                   QC_XC_LYP_Spin_Impl(rho_a, rho_b, sigma_aa, sigma_ab,
                                       sigma_bb);
        case QC_METHOD::PBE0:
            return 0.75 *
                       (QC_XC_GGA_X_Spin_Channel(method, rho_a, sigma_aa) +
                        QC_XC_GGA_X_Spin_Channel(method, rho_b, sigma_bb)) +
                   QC_XC_PBE_C_Spin_Impl(rho_a, rho_b, sigma_total);
        case QC_METHOD::B3LYP:
            return 0.08 * slater +
                   0.72 *
                       (QC_XC_GGA_X_Spin_Channel(method, rho_a, sigma_aa) +
                        QC_XC_GGA_X_Spin_Channel(method, rho_b, sigma_bb)) +
                   0.81 * QC_XC_LYP_Spin_Impl(rho_a, rho_b, sigma_aa,
                                               sigma_ab, sigma_bb) +
                   0.19 * QC_XC_VWN5_Spin_Impl(rho_a, rho_b);
        default:
            return QC_XC_Dual<N>::Constant(QC_XC_Quiet_NaN());
    }
}

static inline __host__ __device__ void QC_VXC_Analytical_RKS(
    QC_METHOD method, double rho, double sigma, double& energy, double& vrho,
    double& vsigma)
{
    energy = vrho = vsigma = 0.0;
    if (!QC_XC_Method_Is_Supported(method) ||
        !QC_XC_Inputs_Valid_RKS(rho, sigma))
    {
        energy = vrho = vsigma = QC_XC_Quiet_NaN();
        return;
    }
    if (rho == 0.0) return;

    const QC_XC_Dual<2> rho_dual = QC_XC_Dual<2>::Variable(rho, 0);
    const QC_XC_Dual<2> sigma_dual = QC_XC_Dual<2>::Variable(sigma, 1);
    const QC_XC_Dual<2> result =
        QC_XC_RKS_Energy_Impl(method, rho_dual, sigma_dual);
    energy = result.value;
    vrho = result.deriv[0];
    vsigma = method == QC_METHOD::LDA ? 0.0 : result.deriv[1];

    // At sigma=0 the PBE exchange and correlation slopes both scale as
    // rho^(-4/3).  Each can overflow long before their analytically combined
    // coefficient does (notably near rho~1e-235).  Combine constants first,
    // then apply the density power once.
    if (sigma == 0.0 && (method == QC_METHOD::PBE || method == QC_METHOD::PBE0))
    {
        const double cx = 0.75 * cbrt(3.0 / CONSTANT_Pi);
        const double mu = 0.2195149727645171;
        const double beta = 0.06672455060314922;
        const double kf_factor = cbrt(3.0 * CONSTANT_Pi * CONSTANT_Pi);
        const double d2 =
            (1.0 / 144.0) * pow(3.0, 5.0 / 3.0) *
            pow(CONSTANT_Pi, 1.0 / 3.0);
        const double exchange_weight = method == QC_METHOD::PBE ? 1.0 : 0.75;
        const double combined =
            beta * d2 - exchange_weight * cx * mu /
                            (4.0 * kf_factor * kf_factor);
        vsigma = QC_XC_Coefficient_Times_Power(combined, rho, -4.0 / 3.0);
    }
}

static inline __host__ __device__ void QC_VXC_Analytical_UKS(
    QC_METHOD method, double rho_a, double rho_b, double sigma_aa,
    double sigma_ab, double sigma_bb, double& energy, double& vrho_a,
    double& vrho_b, double& vsigma_aa, double& vsigma_ab, double& vsigma_bb)
{
    energy = vrho_a = vrho_b = vsigma_aa = vsigma_ab = vsigma_bb = 0.0;
    if (!QC_XC_Method_Is_Supported(method) ||
        !QC_XC_Inputs_Valid_UKS(rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb))
    {
        energy = vrho_a = vrho_b = vsigma_aa = vsigma_ab = vsigma_bb =
            QC_XC_Quiet_NaN();
        return;
    }
    if (rho_a + rho_b == 0.0) return;

    const QC_XC_Dual<5> ra = QC_XC_Dual<5>::Variable(rho_a, 0);
    const QC_XC_Dual<5> rb = QC_XC_Dual<5>::Variable(rho_b, 1);
    const QC_XC_Dual<5> saa = QC_XC_Dual<5>::Variable(sigma_aa, 2);
    const QC_XC_Dual<5> sab = QC_XC_Dual<5>::Variable(sigma_ab, 3);
    const QC_XC_Dual<5> sbb = QC_XC_Dual<5>::Variable(sigma_bb, 4);
    QC_XC_Dual<5> result = QC_XC_UKS_Energy_Impl(method, ra, rb, saa, sab, sbb);

    energy = result.value;
    vrho_a = result.deriv[0];
    vrho_b = result.deriv[1];
    vsigma_aa = method == QC_METHOD::LDA ? 0.0 : result.deriv[2];
    vsigma_ab = method == QC_METHOD::LDA ? 0.0 : result.deriv[3];
    vsigma_bb = method == QC_METHOD::LDA ? 0.0 : result.deriv[4];

    const unsigned infinite_mask = QC_XC_UKS_Expected_Infinite_Output_Mask(
        method, rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb);
    if ((infinite_mask & QC_XC_UKS_VRHO_A_BIT) != 0u)
        vrho_a = QC_XC_Positive_Infinity();
    if ((infinite_mask & QC_XC_UKS_VRHO_B_BIT) != 0u)
        vrho_b = QC_XC_Positive_Infinity();

    // The same zero-gradient overflow cancellation must be done separately
    // for each spin-scaled exchange channel in UKS.
    if ((method == QC_METHOD::PBE || method == QC_METHOD::PBE0) &&
        sigma_aa == 0.0 && sigma_ab == 0.0 && sigma_bb == 0.0)
    {
        const double rho = rho_a + rho_b;
        const double plus = 2.0 * rho_a / rho;
        const double minus = 2.0 * rho_b / rho;
        const double phi =
            0.5 * (pow(plus, 2.0 / 3.0) + pow(minus, 2.0 / 3.0));
        const double cx = 0.75 * cbrt(3.0 / CONSTANT_Pi);
        const double mu = 0.2195149727645171;
        const double beta = 0.06672455060314922;
        const double kf_factor = cbrt(3.0 * CONSTANT_Pi * CONSTANT_Pi);
        const double correlation_constant =
            beta * CONSTANT_Pi * phi / (16.0 * kf_factor);
        const double exchange_constant =
            (method == QC_METHOD::PBE ? 1.0 : 0.75) * cx * mu /
            (4.0 * kf_factor * kf_factor * cbrt(2.0));
        const double log_correlation =
            log(correlation_constant) - (4.0 / 3.0) * log(rho);
        const double correlation = exp(log_correlation);
        vsigma_ab = 2.0 * correlation;
        if (rho_a == 0.0)
            vsigma_aa = correlation;
        else
            vsigma_aa = QC_XC_Exp_Difference(
                log_correlation,
                log(exchange_constant) - (4.0 / 3.0) * log(rho_a));
        if (rho_b == 0.0)
            vsigma_bb = correlation;
        else
            vsigma_bb = QC_XC_Exp_Difference(
                log_correlation,
                log(exchange_constant) - (4.0 / 3.0) * log(rho_b));
    }
}
