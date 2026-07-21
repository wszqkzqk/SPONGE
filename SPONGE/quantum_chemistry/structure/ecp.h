#pragma once

#include <cfloat>
#include <cmath>
#include <initializer_list>
#include <map>
#include <string>
#include <vector>

#if defined(__CUDACC__) || defined(__HIPCC__)
#define QC_ECP_HOST_DEVICE __host__ __device__
#else
#define QC_ECP_HOST_DEVICE
#endif

// ECP 用户选项 ("auto" 按基组自动选、"none" 禁用、显式指定集合)
enum class QC_ECP_TYPE
{
    AUTO,
    NONE,
    DEF2_ECP,
    LANL2DZ,
};

// ECP 单项: d_k · r^(n_k-2) · exp(-ζ_k · r²)
struct QC_ECP_TERM
{
    float d_k;     // 收缩系数
    int n_k;       // r 幂次 (实际幂为 n_k - 2)
    float zeta_k;  // 高斯指数
};

// ECP 角动量通道
struct QC_ECP_CHANNEL
{
    int l;  // 角动量 (-1 表示 local / ul)
    std::vector<QC_ECP_TERM> terms;
};

// 单个原子的 ECP 参数
struct QC_ECP_ATOM_DATA
{
    int n_core = 0;  // 被替代的核心电子数
    int l_max = -1;  // 最大角动量 (local channel 的 l 值)
    std::vector<QC_ECP_CHANNEL> channels;
    // channels[l_max] = local (U_L)
    // channels[0..l_max-1] = semi-local (U_l)
};

// The current semi-local integral backend implements projectors through f
// (l=3).  A g-labelled local channel (l_max=4) is still representable because
// the local radial operator is evaluated without an angular projector.
inline constexpr int QC_ECP_MAX_SEMILOCAL_L = 3;

inline constexpr bool QC_ECP_L_Max_Is_Supported(int l_max) noexcept
{
    return l_max >= 0 && l_max <= QC_ECP_MAX_SEMILOCAL_L + 1;
}

// A primitive Gaussian displaced from an ECP center has an infinite angular
// expansion. Terms in a projector-l channel occur at q=l,l+2,... . The final
// four layers are checked, so exhausting this bound is a reported numerical
// failure rather than a silently truncated integral.
inline constexpr int QC_ECP_MAX_SERIES_ORDER = 64;
inline constexpr int QC_ECP_MAX_PROJECTOR_COMPONENTS =
    2 * QC_ECP_MAX_SEMILOCAL_L + 1;
inline constexpr int QC_ECP_MAX_SERIES_LAYERS =
    QC_ECP_MAX_SERIES_ORDER / 2 + 1;
inline constexpr int QC_ECP_MAX_QUADRATURE_RADIAL = 64;
inline constexpr int QC_ECP_SEMILOCAL_QUADRATURE_RADIAL = 64;
inline constexpr int QC_ECP_SEMILOCAL_QUADRATURE_RADIAL_COARSE = 52;
inline constexpr int QC_ECP_LOCAL_QUADRATURE_RADIAL_COARSE = 52;
inline constexpr int QC_ECP_MAX_QUADRATURE_POLAR = 64;
inline constexpr int QC_ECP_MAX_QUADRATURE_AZIMUTHAL = 128;

struct QC_ECP_PRIMITIVE_RESULT
{
    double value = 0.0;
    double last_layer_bound = 0.0;
    bool converged = false;
};

struct QC_ECP_PRIMITIVE_DERIVATIVE_RESULT
{
    double value = 0.0;
    double derivative_i[3] = {};
    double derivative_j[3] = {};
    double derivative_i_error[3] = {};
    double derivative_j_error[3] = {};
    double estimated_error = 0.0;
    bool converged = false;
};

namespace qc_ecp_math
{
inline constexpr double PI = 3.141592653589793238462643383279502884;
inline constexpr double SQRT_FOUR_PI = 3.544907701811032054596334966682290365;
inline constexpr double INV_SQRT_FOUR_PI =
    0.2820947917738781434740397257803862929;

QC_ECP_HOST_DEVICE inline double Abs(double value)
{
    return value < 0.0 ? -value : value;
}

QC_ECP_HOST_DEVICE inline int Select_Series_Order(double leading_order,
                                                   int parity)
{
    const int candidates[5] = {16, 24, 32, 48,
                               QC_ECP_MAX_SERIES_ORDER};
    for (int i = 0; i < 5; ++i)
    {
        int order = candidates[i];
        if ((order & 1) != parity) --order;
        if (static_cast<double>(order) >= leading_order) return order;
    }
    return QC_ECP_MAX_SERIES_ORDER -
           ((QC_ECP_MAX_SERIES_ORDER & 1) != parity ? 1 : 0);
}

// The integrated Taylor layers have a Poisson-like peak near
// 2 alpha^2 |C-A|^2/eta.  Requiring twelve standard deviations beyond that
// peak prevents a tiny pre-peak layer from being mistaken for convergence.
QC_ECP_HOST_DEVICE inline double Semilocal_Series_Extent(
    double alpha_i, double alpha_j, double distance_i2, double distance_j2,
    double eta, int angular_i, int angular_j, int projector_l)
{
    const double peak_i = 2.0 * alpha_i * alpha_i * distance_i2 / eta;
    const double peak_j = 2.0 * alpha_j * alpha_j * distance_j2 / eta;
    const double extent_i = peak_i + 12.0 * sqrt(peak_i + 1.0) + angular_i +
                            projector_l;
    const double extent_j = peak_j + 12.0 * sqrt(peak_j + 1.0) + angular_j +
                            projector_l;
    return extent_i > extent_j ? extent_i : extent_j;
}

QC_ECP_HOST_DEVICE inline int Semilocal_Series_Order(
    double alpha_i, double alpha_j, double distance_i2, double distance_j2,
    double eta, int angular_i, int angular_j, int projector_l)
{
    return Select_Series_Order(
        Semilocal_Series_Extent(alpha_i, alpha_j, distance_i2, distance_j2,
                                eta, angular_i, angular_j, projector_l),
        projector_l & 1);
}

QC_ECP_HOST_DEVICE inline double Local_Series_Extent(
    double alpha_i, double alpha_j, double dx_i, double dy_i, double dz_i,
    double dx_j, double dy_j, double dz_j, double eta, int angular_degree)
{
    const double wx = alpha_i * dx_i + alpha_j * dx_j;
    const double wy = alpha_i * dy_i + alpha_j * dy_j;
    const double wz = alpha_i * dz_i + alpha_j * dz_j;
    const double peak = 2.0 * (wx * wx + wy * wy + wz * wz) / eta;
    return peak + 12.0 * sqrt(peak + 1.0) + angular_degree;
}

QC_ECP_HOST_DEVICE inline int Local_Series_Order(
    double alpha_i, double alpha_j, double dx_i, double dy_i, double dz_i,
    double dx_j, double dy_j, double dz_j, double eta, int angular_degree)
{
    return Select_Series_Order(
        Local_Series_Extent(alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j,
                            dz_j, eta, angular_degree),
        0);
}

QC_ECP_HOST_DEVICE inline double Integer_Power(double base, int exponent)
{
    double value = 1.0;
    for (int i = 0; i < exponent; ++i) value *= base;
    return value;
}

QC_ECP_HOST_DEVICE inline double Binomial(int n, int k)
{
    if (k < 0 || k > n) return 0.0;
    if (k > n - k) k = n - k;
    double value = 1.0;
    for (int i = 1; i <= k; ++i)
        value *= static_cast<double>(n - k + i) / static_cast<double>(i);
    return value;
}

// Integral over the unit sphere of x^a y^b z^c. Odd exponents vanish; for
// a=2A,b=2B,c=2C the even moment is
//   4 pi (2A-1)!!(2B-1)!!(2C-1)!!/(2A+2B+2C+1)!! .
QC_ECP_HOST_DEVICE inline double Sphere_Monomial_Moment(int a, int b, int c)
{
    if (a < 0 || b < 0 || c < 0 || (a & 1) || (b & 1) || (c & 1))
        return 0.0;
    double value = 4.0 * PI;
    for (int i = 1; i <= a / 2; ++i) value *= static_cast<double>(2 * i - 1);
    for (int i = 1; i <= b / 2; ++i) value *= static_cast<double>(2 * i - 1);
    for (int i = 1; i <= c / 2; ++i) value *= static_cast<double>(2 * i - 1);
    for (int i = 1; i <= (a + b + c) / 2; ++i)
        value /= static_cast<double>(2 * i + 1);
    return value;
}

// Integral Y_lm(real, orthonormal) x^a y^b z^c dOmega. The real-harmonic
// order is m=-l,...,+l. Writing the harmonics as normalized solid-harmonic
// polynomials reduces every value to exact sphere monomial moments.
QC_ECP_HOST_DEVICE inline double Real_Spherical_Monomial_Moment(
    int l, int m_index, int a, int b, int c)
{
    if (l == 0)
        return INV_SQRT_FOUR_PI * Sphere_Monomial_Moment(a, b, c);
    if (l == 1)
    {
        const double k = 0.48860251190291992159;  // sqrt(3/(4 pi))
        if (m_index == 0)
            return k * Sphere_Monomial_Moment(a, b + 1, c);
        if (m_index == 1)
            return k * Sphere_Monomial_Moment(a, b, c + 1);
        return k * Sphere_Monomial_Moment(a + 1, b, c);
    }
    if (l == 2)
    {
        if (m_index == 0)
            return 1.0925484305920790705 *
                   Sphere_Monomial_Moment(a + 1, b + 1, c);
        if (m_index == 1)
            return 1.0925484305920790705 *
                   Sphere_Monomial_Moment(a, b + 1, c + 1);
        if (m_index == 2)
            return 0.31539156525252000603 *
                   (2.0 * Sphere_Monomial_Moment(a, b, c + 2) -
                    Sphere_Monomial_Moment(a + 2, b, c) -
                    Sphere_Monomial_Moment(a, b + 2, c));
        if (m_index == 3)
            return 1.0925484305920790705 *
                   Sphere_Monomial_Moment(a + 1, b, c + 1);
        return 0.54627421529603953527 *
               (Sphere_Monomial_Moment(a + 2, b, c) -
                Sphere_Monomial_Moment(a, b + 2, c));
    }
    if (l == 3)
    {
        if (m_index == 0)
            return 0.59004358992664351035 *
                   (3.0 * Sphere_Monomial_Moment(a + 2, b + 1, c) -
                    Sphere_Monomial_Moment(a, b + 3, c));
        if (m_index == 1)
            return 2.8906114426405538302 *
                   Sphere_Monomial_Moment(a + 1, b + 1, c + 1);
        if (m_index == 2)
            return 0.45704579946446573616 *
                   (4.0 * Sphere_Monomial_Moment(a, b + 1, c + 2) -
                    Sphere_Monomial_Moment(a + 2, b + 1, c) -
                    Sphere_Monomial_Moment(a, b + 3, c));
        if (m_index == 3)
            return 0.37317633259011539141 *
                   (2.0 * Sphere_Monomial_Moment(a, b, c + 3) -
                    3.0 * Sphere_Monomial_Moment(a + 2, b, c + 1) -
                    3.0 * Sphere_Monomial_Moment(a, b + 2, c + 1));
        if (m_index == 4)
            return 0.45704579946446573616 *
                   (4.0 * Sphere_Monomial_Moment(a + 1, b, c + 2) -
                    Sphere_Monomial_Moment(a + 3, b, c) -
                    Sphere_Monomial_Moment(a + 1, b + 2, c));
        if (m_index == 5)
            return 1.4453057213202769151 *
                   (Sphere_Monomial_Moment(a + 2, b, c + 1) -
                    Sphere_Monomial_Moment(a, b + 2, c + 1));
        return 0.59004358992664351035 *
               (Sphere_Monomial_Moment(a + 3, b, c) -
                3.0 * Sphere_Monomial_Moment(a + 1, b + 2, c));
    }
    return 0.0;
}

// Coefficients of (d+x)^L exp(-2 alpha d x) = sum_p f_p x^p.
QC_ECP_HOST_DEVICE inline void One_Dimensional_Series(
    int angular_momentum, double alpha, double displacement, int max_order,
    double* output)
{
    double exponential[QC_ECP_MAX_SERIES_ORDER + 1];
    exponential[0] = 1.0;
    const double ratio = -2.0 * alpha * displacement;
    for (int p = 1; p <= max_order; ++p)
        exponential[p] = exponential[p - 1] * ratio / static_cast<double>(p);
    for (int p = 0; p <= max_order; ++p)
    {
        double value = 0.0;
        const int k_max = angular_momentum < p ? angular_momentum : p;
        for (int k = 0; k <= k_max; ++k)
            value += Binomial(angular_momentum, k) *
                     Integer_Power(displacement, angular_momentum - k) *
                     exponential[p - k];
        output[p] = value;
    }
}

// Coefficients after multiplying the bra and ket factors in one Cartesian
// direction. Local angular integration acts once on chi_i*chi_j and is not a
// product of two P_0 projections.
QC_ECP_HOST_DEVICE inline bool Combined_One_Dimensional_Series(
    int angular_i, double alpha_i, double displacement_i, int angular_j,
    double alpha_j, double displacement_j, int max_order, double* output)
{
    const int polynomial_degree = angular_i + angular_j;
    if (polynomial_degree > max_order) return false;
    double polynomial[QC_ECP_MAX_SERIES_ORDER + 1] = {};
    for (int i = 0; i <= angular_i; ++i)
        for (int j = 0; j <= angular_j; ++j)
            polynomial[i + j] +=
                Binomial(angular_i, i) *
                Integer_Power(displacement_i, angular_i - i) *
                Binomial(angular_j, j) *
                Integer_Power(displacement_j, angular_j - j);
    double exponential[QC_ECP_MAX_SERIES_ORDER + 1];
    exponential[0] = 1.0;
    const double ratio =
        -2.0 * (alpha_i * displacement_i + alpha_j * displacement_j);
    for (int p = 1; p <= max_order; ++p)
        exponential[p] = exponential[p - 1] * ratio / static_cast<double>(p);
    for (int p = 0; p <= max_order; ++p)
    {
        double value = 0.0;
        const int k_max = polynomial_degree < p ? polynomial_degree : p;
        for (int k = 0; k <= k_max; ++k)
            value += polynomial[k] * exponential[p - k];
        output[p] = value;
    }
    return true;
}

QC_ECP_HOST_DEVICE inline void Build_Projection_Series(
    int projector_l, int lx, int ly, int lz, double alpha, double dx,
    double dy, double dz, int max_order,
    double output[QC_ECP_MAX_PROJECTOR_COMPONENTS]
                 [QC_ECP_MAX_SERIES_LAYERS])
{
    double fx[QC_ECP_MAX_SERIES_ORDER + 1];
    double fy[QC_ECP_MAX_SERIES_ORDER + 1];
    double fz[QC_ECP_MAX_SERIES_ORDER + 1];
    One_Dimensional_Series(lx, alpha, dx, max_order, fx);
    One_Dimensional_Series(ly, alpha, dy, max_order, fy);
    One_Dimensional_Series(lz, alpha, dz, max_order, fz);
    const int component_count = 2 * projector_l + 1;
    const int layer_count = (max_order - projector_l) / 2 + 1;
    for (int m = 0; m < component_count; ++m)
        for (int layer = 0; layer < layer_count; ++layer)
            output[m][layer] = 0.0;
    for (int layer = 0; layer < layer_count; ++layer)
    {
        const int q = projector_l + 2 * layer;
        for (int a = 0; a <= q; ++a)
            for (int b = 0; b <= q - a; ++b)
            {
                const int c = q - a - b;
                const double coefficient = fx[a] * fy[b] * fz[c];
                for (int m = 0; m < component_count; ++m)
                    output[m][layer] +=
                        coefficient * Real_Spherical_Monomial_Moment(
                                          projector_l, m, a, b, c);
            }
    }
}

QC_ECP_HOST_DEVICE inline bool Build_Local_Angular_Series(
    int lx_i, int ly_i, int lz_i, double alpha_i, double dx_i, double dy_i,
    double dz_i, int lx_j, int ly_j, int lz_j, double alpha_j, double dx_j,
    double dy_j, double dz_j, int max_order, double* output)
{
    double fx[QC_ECP_MAX_SERIES_ORDER + 1];
    double fy[QC_ECP_MAX_SERIES_ORDER + 1];
    double fz[QC_ECP_MAX_SERIES_ORDER + 1];
    if (!Combined_One_Dimensional_Series(lx_i, alpha_i, dx_i, lx_j,
                                         alpha_j, dx_j, max_order, fx) ||
        !Combined_One_Dimensional_Series(ly_i, alpha_i, dy_i, ly_j,
                                         alpha_j, dy_j, max_order, fy) ||
        !Combined_One_Dimensional_Series(lz_i, alpha_i, dz_i, lz_j,
                                         alpha_j, dz_j, max_order, fz))
        return false;
    for (int layer = 0; layer <= max_order / 2; ++layer)
    {
        const int q = 2 * layer;
        double angular = 0.0;
        for (int a = 0; a <= q; ++a)
            for (int b = 0; b <= q - a; ++b)
            {
                const int c = q - a - b;
                angular += fx[a] * fy[b] * fz[c] *
                           Sphere_Monomial_Moment(a, b, c);
            }
        output[layer] = angular;
    }
    return true;
}

// Integral int_0^infinity r^power exp(-eta r^2) dr, evaluated without a
// gamma-library dependency so the same policy is available on CPU/GPU.
QC_ECP_HOST_DEVICE inline double Gaussian_Radial_Moment(int power, double eta)
{
    if ((power & 1) == 0)
    {
        double value = sqrt(PI) / (2.0 * sqrt(eta));
        for (int i = 1; i <= power / 2; ++i)
            value *= static_cast<double>(2 * i - 1) / (2.0 * eta);
        return value;
    }
    double value = 1.0 / (2.0 * eta);
    for (int i = 1; i <= (power - 1) / 2; ++i)
        value *= static_cast<double>(i) / eta;
    return value;
}

// Full-line Gaussian moment used by the exact local n_k=2 overlap path.
QC_ECP_HOST_DEVICE inline double Gaussian_Line_Moment(int power, double eta)
{
    if (power < 0 || (power & 1)) return 0.0;
    double value = sqrt(PI / eta);
    for (int i = 1; i <= power / 2; ++i)
        value *= static_cast<double>(2 * i - 1) / (2.0 * eta);
    return value;
}

QC_ECP_HOST_DEVICE inline double Gaussian_Polynomial_Overlap_1D(
    int angular_i, int angular_j, double center_minus_i,
    double center_minus_j, double eta)
{
    double value = 0.0;
    for (int i = 0; i <= angular_i; ++i)
        for (int j = 0; j <= angular_j; ++j)
            value += Binomial(angular_i, i) *
                     Integer_Power(center_minus_i, angular_i - i) *
                     Binomial(angular_j, j) *
                     Integer_Power(center_minus_j, angular_j - j) *
                     Gaussian_Line_Moment(i + j, eta);
    return value;
}

// Exact integral for a local n_k=2 primitive.  The ECP Gaussian is combined
// with the two basis Gaussians before integrating the Cartesian polynomials.
// The exponent is written in pair-distance form, which is non-negative and
// avoids cancellation for remote or very diffuse primitives.  There is no
// magnitude cutoff: a small Gaussian prefactor can be amplified by the exact
// volume/polynomial factor into a material integral.
QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_RESULT Local_Primitive_N2(
    double alpha_i, double alpha_j, double Ax, double Ay, double Az,
    double Bx, double By, double Bz, double Cx, double Cy, double Cz,
    double zeta, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j, int lz_j)
{
    QC_ECP_PRIMITIVE_RESULT result;
    if (!(alpha_i > 0.0) || !(alpha_j > 0.0) || zeta < 0.0 || lx_i < 0 ||
        ly_i < 0 || lz_i < 0 || lx_j < 0 || ly_j < 0 || lz_j < 0)
        return result;
    const double eta = alpha_i + alpha_j + zeta;
    const double Qx = (alpha_i * Ax + alpha_j * Bx + zeta * Cx) / eta;
    const double Qy = (alpha_i * Ay + alpha_j * By + zeta * Cy) / eta;
    const double Qz = (alpha_i * Az + alpha_j * Bz + zeta * Cz) / eta;
    const double ABx = Ax - Bx, ABy = Ay - By, ABz = Az - Bz;
    const double ACx = Ax - Cx, ACy = Ay - Cy, ACz = Az - Cz;
    const double BCx = Bx - Cx, BCy = By - Cy, BCz = Bz - Cz;
    const double exponent =
        (alpha_i * alpha_j *
             (ABx * ABx + ABy * ABy + ABz * ABz) +
         alpha_i * zeta * (ACx * ACx + ACy * ACy + ACz * ACz) +
         alpha_j * zeta * (BCx * BCx + BCy * BCy + BCz * BCz)) /
        eta;
    const double overlap_x = Gaussian_Polynomial_Overlap_1D(
        lx_i, lx_j, Qx - Ax, Qx - Bx, eta);
    const double overlap_y = Gaussian_Polynomial_Overlap_1D(
        ly_i, ly_j, Qy - Ay, Qy - By, eta);
    const double overlap_z = Gaussian_Polynomial_Overlap_1D(
        lz_i, lz_j, Qz - Az, Qz - Bz, eta);
    result.value = exp(-exponent) * overlap_x * overlap_y * overlap_z;
    result.converged = true;
    return result;
}

QC_ECP_HOST_DEVICE inline void Kahan_Add(double term, double& sum,
                                          double& correction)
{
    const double adjusted = term - correction;
    const double updated = sum + adjusted;
    correction = (updated - sum) - adjusted;
    sum = updated;
}

QC_ECP_HOST_DEVICE inline bool Tail_Is_Converged(const double* layers,
                                                  int layer_count,
                                                  double value,
                                                  double& bound)
{
    bound = 0.0;
    const int first = layer_count > 4 ? layer_count - 4 : 0;
    for (int layer = first; layer < layer_count; ++layer)
        bound += Abs(layers[layer]);
    const double tolerance = 2.0e-11 * (1.0 + Abs(value));
    return layer_count >= 4 && bound <= tolerance;
}

QC_ECP_HOST_DEVICE inline double Real_Spherical_Harmonic(
    int l, int m_index, double x, double y, double z)
{
    if (l == 0) return INV_SQRT_FOUR_PI;
    if (l == 1)
    {
        if (m_index == 0) return 0.48860251190291992159 * y;
        if (m_index == 1) return 0.48860251190291992159 * z;
        return 0.48860251190291992159 * x;
    }
    if (l == 2)
    {
        if (m_index == 0) return 1.0925484305920790705 * x * y;
        if (m_index == 1) return 1.0925484305920790705 * y * z;
        if (m_index == 2)
            return 0.31539156525252000603 *
                   (2.0 * z * z - x * x - y * y);
        if (m_index == 3) return 1.0925484305920790705 * x * z;
        return 0.54627421529603953527 * (x * x - y * y);
    }
    if (m_index == 0)
        return 0.59004358992664351035 * y * (3.0 * x * x - y * y);
    if (m_index == 1) return 2.8906114426405538302 * x * y * z;
    if (m_index == 2)
        return 0.45704579946446573616 * y *
               (4.0 * z * z - x * x - y * y);
    if (m_index == 3)
        return 0.37317633259011539141 * z *
               (2.0 * z * z - 3.0 * x * x - 3.0 * y * y);
    if (m_index == 4)
        return 0.45704579946446573616 * x *
               (4.0 * z * z - x * x - y * y);
    if (m_index == 5)
        return 1.4453057213202769151 * z * (x * x - y * y);
    return 0.59004358992664351035 * x * (x * x - 3.0 * y * y);
}

QC_ECP_HOST_DEVICE inline void Gauss_Legendre_Rule(int count, double lower,
                                                    double upper,
                                                    double* nodes,
                                                    double* weights)
{
    const int half = (count + 1) / 2;
    const double midpoint = 0.5 * (lower + upper);
    const double half_width = 0.5 * (upper - lower);
    for (int root = 0; root < half; ++root)
    {
        double z = cos(PI * (static_cast<double>(root) + 0.75) /
                       (static_cast<double>(count) + 0.5));
        double derivative = 0.0;
        for (int iteration = 0; iteration < 32; ++iteration)
        {
            double p0 = 1.0, p1 = z;
            if (count == 0) p1 = p0;
            for (int degree = 2; degree <= count; ++degree)
            {
                const double next =
                    ((2.0 * degree - 1.0) * z * p1 -
                     (degree - 1.0) * p0) /
                    static_cast<double>(degree);
                p0 = p1;
                p1 = next;
            }
            derivative = static_cast<double>(count) * (z * p1 - p0) /
                         (z * z - 1.0);
            const double next_z = z - p1 / derivative;
            if (Abs(next_z - z) <= 4.0 * DBL_EPSILON)
            {
                z = next_z;
                break;
            }
            z = next_z;
        }
        const double weight =
            half_width * 2.0 / ((1.0 - z * z) * derivative * derivative);
        nodes[root] = midpoint - half_width * z;
        weights[root] = weight;
        const int mirror = count - 1 - root;
        if (mirror != root)
        {
            nodes[mirror] = midpoint + half_width * z;
            weights[mirror] = weight;
        }
    }
}

// Nodes and weights for integral_0^infinity exp(-t) f(t) dt.  The Newton
// iteration follows the Laguerre recurrence and the root weight identity
// w_i=1/(x_i L'_n(x_i)^2), avoiding tables and host-only special functions.
QC_ECP_HOST_DEVICE inline void Gauss_Laguerre_Rule(int count, double* nodes,
                                                    double* weights)
{
    double root = 0.0;
    for (int index = 0; index < count; ++index)
    {
        if (index == 0)
            root = 3.0 / (1.0 + 2.4 * static_cast<double>(count));
        else if (index == 1)
            root += 15.0 / (1.0 + 2.5 * static_cast<double>(count));
        else
        {
            const double previous_index = static_cast<double>(index - 1);
            root += ((1.0 + 2.55 * previous_index) /
                     (1.9 * previous_index)) *
                    (root - nodes[index - 2]);
        }
        double derivative = 0.0;
        for (int iteration = 0; iteration < 32; ++iteration)
        {
            double previous = 1.0;
            double polynomial = 1.0 - root;
            if (count == 0) polynomial = previous;
            for (int degree = 2; degree <= count; ++degree)
            {
                const double next =
                    ((2.0 * degree - 1.0 - root) * polynomial -
                     (degree - 1.0) * previous) /
                    static_cast<double>(degree);
                previous = polynomial;
                polynomial = next;
            }
            derivative = static_cast<double>(count) *
                         (polynomial - previous) / root;
            const double next_root = root - polynomial / derivative;
            if (Abs(next_root - root) <=
                8.0 * DBL_EPSILON * (1.0 + Abs(next_root)))
            {
                root = next_root;
                break;
            }
            root = next_root;
        }
        // Recompute the derivative at the accepted root.
        double previous = 1.0;
        double polynomial = 1.0 - root;
        for (int degree = 2; degree <= count; ++degree)
        {
            const double next =
                ((2.0 * degree - 1.0 - root) * polynomial -
                 (degree - 1.0) * previous) /
                static_cast<double>(degree);
            previous = polynomial;
            polynomial = next;
        }
        derivative = static_cast<double>(count) *
                     (polynomial - previous) / root;
        nodes[index] = root;
        weights[index] = 1.0 / (root * derivative * derivative);
    }
}

struct Angular_Frame
{
    double axis_x;
    double axis_y;
    double axis_z;
    double frame_x_x;
    double frame_x_y;
    double frame_x_z;
    double frame_y_x;
    double frame_y_y;
    double frame_y_z;
};

QC_ECP_HOST_DEVICE inline Angular_Frame Build_Angular_Frame(
    double vx, double vy, double vz, double norm)
{
    Angular_Frame frame = {0.0, 0.0, 1.0, 1.0, 0.0,
                           0.0, 0.0, 1.0, 0.0};
    if (!(norm > 0.0)) return frame;
    frame.axis_x = vx / norm;
    frame.axis_y = vy / norm;
    frame.axis_z = vz / norm;
    const double norm_xy =
        sqrt(frame.axis_x * frame.axis_x + frame.axis_y * frame.axis_y);
    if (norm_xy > 0.1)
    {
        frame.frame_x_x = -frame.axis_y / norm_xy;
        frame.frame_x_y = frame.axis_x / norm_xy;
        frame.frame_x_z = 0.0;
        frame.frame_y_x = -frame.axis_z * frame.axis_x / norm_xy;
        frame.frame_y_y = -frame.axis_z * frame.axis_y / norm_xy;
        frame.frame_y_z = norm_xy;
    }
    else
    {
        const double norm_yz =
            sqrt(frame.axis_y * frame.axis_y + frame.axis_z * frame.axis_z);
        frame.frame_x_x = 0.0;
        frame.frame_x_y = frame.axis_z / norm_yz;
        frame.frame_x_z = -frame.axis_y / norm_yz;
        frame.frame_y_x = -norm_yz;
        frame.frame_y_y = frame.axis_x * frame.axis_y / norm_yz;
        frame.frame_y_z = frame.axis_x * frame.axis_z / norm_yz;
    }
    return frame;
}

QC_ECP_HOST_DEVICE inline void Rotate_Angular_Direction(
    const Angular_Frame& frame, double polar, double cosine_phi,
    double sine_phi, double& omega_x, double& omega_y, double& omega_z)
{
    const double transverse =
        sqrt((1.0 - polar * polar) > 0.0 ? (1.0 - polar * polar) : 0.0);
    const double frame_x = transverse * cosine_phi;
    const double frame_y = transverse * sine_phi;
    omega_x = frame_x * frame.frame_x_x + frame_y * frame.frame_y_x +
              polar * frame.axis_x;
    omega_y = frame_x * frame.frame_x_y + frame_y * frame.frame_y_y +
              polar * frame.axis_y;
    omega_z = frame_x * frame.frame_x_z + frame_y * frame.frame_y_z +
              polar * frame.axis_z;
}

// For a sharply peaked scaled exponential exp[-k(1+u)], t=k(1+u)
// converts the polar integral to a Laguerre-weighted integral.  For k>20 the
// omitted t>2k tail is below exp(-40); moderate k remains on finite-interval
// Gauss-Legendre quadrature.  This prevents two capped grids from jointly
// missing the endpoint layer.
QC_ECP_HOST_DEVICE inline bool Scaled_Polar_Quadrature_Point(
    int index, double sharpness, const double* legendre_nodes,
    const double* legendre_weights, const double* laguerre_nodes,
    const double* laguerre_weights, double& polar, double& weight,
    double& exponential)
{
    if (sharpness > 20.0)
    {
        const double t = laguerre_nodes[index];
        if (t >= 2.0 * sharpness) return false;
        polar = -1.0 + t / sharpness;
        weight = laguerre_weights[index] / sharpness;
        exponential = 1.0;
        return true;
    }
    polar = legendre_nodes[index];
    weight = legendre_weights[index];
    exponential = exp(-sharpness * (polar + 1.0));
    return true;
}

QC_ECP_HOST_DEVICE inline void Angular_Projection_Quadrature(
    double radius, int projector_l, int lx, int ly, int lz, double alpha,
    double dx, double dy, double dz, int polar_count,
    const double* polar_nodes, const double* polar_weights,
    const double* laguerre_nodes, const double* laguerre_weights,
    int azimuthal_count, const double* cosine_phi, const double* sine_phi,
    double* projection)
{
    const int component_count = 2 * projector_l + 1;
    double correction[QC_ECP_MAX_PROJECTOR_COMPONENTS] = {};
    for (int m = 0; m < component_count; ++m) projection[m] = 0.0;
    const double distance = sqrt(dx * dx + dy * dy + dz * dz);
    const Angular_Frame frame = Build_Angular_Frame(dx, dy, dz, distance);
    const double sharpness = 2.0 * alpha * radius * distance;
    const double phi_weight = 2.0 * PI / static_cast<double>(azimuthal_count);
    for (int iu = 0; iu < polar_count; ++iu)
    {
        double polar = 0.0, polar_weight = 0.0, exponential = 0.0;
        if (!Scaled_Polar_Quadrature_Point(
                iu, sharpness, polar_nodes, polar_weights, laguerre_nodes,
                laguerre_weights, polar, polar_weight, exponential))
            continue;
        for (int iphi = 0; iphi < azimuthal_count; ++iphi)
        {
            double x = 0.0, y = 0.0, z = 0.0;
            Rotate_Angular_Direction(frame, polar, cosine_phi[iphi],
                                     sine_phi[iphi], x, y, z);
            const double basis =
                Integer_Power(dx + radius * x, lx) *
                Integer_Power(dy + radius * y, ly) *
                Integer_Power(dz + radius * z, lz) * exponential;
            const double weight = polar_weight * phi_weight;
            for (int m = 0; m < component_count; ++m)
                Kahan_Add(weight * basis *
                              Real_Spherical_Harmonic(projector_l, m, x, y, z),
                          projection[m], correction[m]);
        }
    }
}

QC_ECP_HOST_DEVICE inline void Scaled_Primitive_Basis_And_Derivative(
    double radius, double omega_x, double omega_y, double omega_z, int lx,
    int ly, int lz, double alpha, double dx, double dy, double dz,
    double exponential, double& basis, double derivative[3])
{
    const double yx = dx + radius * omega_x;
    const double yy = dy + radius * omega_y;
    const double yz = dz + radius * omega_z;
    const double px = Integer_Power(yx, lx);
    const double py = Integer_Power(yy, ly);
    const double pz = Integer_Power(yz, lz);
    const double polynomial = px * py * pz;
    basis = polynomial * exponential;
    derivative[0] =
        (2.0 * alpha * yx * polynomial -
         static_cast<double>(lx) *
             (lx > 0 ? Integer_Power(yx, lx - 1) * py * pz : 0.0)) *
        exponential;
    derivative[1] =
        (2.0 * alpha * yy * polynomial -
         static_cast<double>(ly) *
             (ly > 0 ? px * Integer_Power(yy, ly - 1) * pz : 0.0)) *
        exponential;
    derivative[2] =
        (2.0 * alpha * yz * polynomial -
         static_cast<double>(lz) *
             (lz > 0 ? px * py * Integer_Power(yz, lz - 1) : 0.0)) *
        exponential;
}

// Evaluate a primitive's angular projection and its three basis-center
// derivatives at the same quadrature points.  With y=C-A+r*Omega,
//
//   d chi / d A_k = (2 alpha y_k P - l_k P_{l_k-1}) exp(-alpha y^2).
//
// The exponential below is peak-scaled in exactly the same way as the base
// projection.  The complementary scale is supplied by the radial quadrature,
// so these are derivatives of the physical Gaussian rather than derivatives
// of the scaled representation.
QC_ECP_HOST_DEVICE inline void Angular_Projection_Derivative_Quadrature(
    double radius, int projector_l, int lx, int ly, int lz, double alpha,
    double dx, double dy, double dz, int polar_count,
    const double* polar_nodes, const double* polar_weights,
    const double* laguerre_nodes, const double* laguerre_weights,
    int azimuthal_count, const double* cosine_phi, const double* sine_phi,
    double* projection,
    double derivative[3][QC_ECP_MAX_PROJECTOR_COMPONENTS])
{
    const int component_count = 2 * projector_l + 1;
    double correction[QC_ECP_MAX_PROJECTOR_COMPONENTS] = {};
    double derivative_correction[3][QC_ECP_MAX_PROJECTOR_COMPONENTS] = {};
    for (int m = 0; m < component_count; ++m)
    {
        projection[m] = 0.0;
        for (int direction = 0; direction < 3; ++direction)
            derivative[direction][m] = 0.0;
    }
    const double distance = sqrt(dx * dx + dy * dy + dz * dz);
    const Angular_Frame frame = Build_Angular_Frame(dx, dy, dz, distance);
    const double sharpness = 2.0 * alpha * radius * distance;
    const double phi_weight = 2.0 * PI / static_cast<double>(azimuthal_count);
    for (int iu = 0; iu < polar_count; ++iu)
    {
        double polar = 0.0, polar_weight = 0.0, exponential = 0.0;
        if (!Scaled_Polar_Quadrature_Point(
                iu, sharpness, polar_nodes, polar_weights, laguerre_nodes,
                laguerre_weights, polar, polar_weight, exponential))
            continue;
        for (int iphi = 0; iphi < azimuthal_count; ++iphi)
        {
            double omega_x = 0.0, omega_y = 0.0, omega_z = 0.0;
            Rotate_Angular_Direction(frame, polar, cosine_phi[iphi],
                                     sine_phi[iphi], omega_x, omega_y,
                                     omega_z);
            double basis = 0.0;
            double derivative_basis[3];
            Scaled_Primitive_Basis_And_Derivative(
                radius, omega_x, omega_y, omega_z, lx, ly, lz, alpha, dx, dy,
                dz, exponential, basis, derivative_basis);
            const double weight = polar_weight * phi_weight;
            for (int m = 0; m < component_count; ++m)
            {
                const double harmonic = Real_Spherical_Harmonic(
                    projector_l, m, omega_x, omega_y, omega_z);
                Kahan_Add(weight * basis * harmonic, projection[m],
                          correction[m]);
                for (int direction = 0; direction < 3; ++direction)
                    Kahan_Add(weight * derivative_basis[direction] * harmonic,
                              derivative[direction][m],
                              derivative_correction[direction][m]);
            }
        }
    }
}

// Bra and ket use the same angular grid.  Accumulating both projections in a
// single traversal avoids reevaluating nodes, weights, trigonometry and real
// spherical harmonics in the expensive quadrature fallback.
QC_ECP_HOST_DEVICE inline void
Angular_Projection_Pair_Derivative_Quadrature(
    double radius, int projector_l, int lx_i, int ly_i, int lz_i,
    double alpha_i, double dx_i, double dy_i, double dz_i, int lx_j, int ly_j,
    int lz_j, double alpha_j, double dx_j, double dy_j, double dz_j,
    int polar_count, const double* polar_nodes, const double* polar_weights,
    const double* laguerre_nodes, const double* laguerre_weights,
    int azimuthal_count, const double* cosine_phi, const double* sine_phi,
    double* projection_i,
    double derivative_i[3][QC_ECP_MAX_PROJECTOR_COMPONENTS],
    double* projection_j,
    double derivative_j[3][QC_ECP_MAX_PROJECTOR_COMPONENTS])
{
    const int component_count = 2 * projector_l + 1;
    double correction_i[QC_ECP_MAX_PROJECTOR_COMPONENTS] = {};
    double correction_j[QC_ECP_MAX_PROJECTOR_COMPONENTS] = {};
    double derivative_correction_i[3][QC_ECP_MAX_PROJECTOR_COMPONENTS] = {};
    double derivative_correction_j[3][QC_ECP_MAX_PROJECTOR_COMPONENTS] = {};
    for (int m = 0; m < component_count; ++m)
    {
        projection_i[m] = 0.0;
        projection_j[m] = 0.0;
        for (int direction = 0; direction < 3; ++direction)
        {
            derivative_i[direction][m] = 0.0;
            derivative_j[direction][m] = 0.0;
        }
    }
    const double distance_i =
        sqrt(dx_i * dx_i + dy_i * dy_i + dz_i * dz_i);
    const double distance_j =
        sqrt(dx_j * dx_j + dy_j * dy_j + dz_j * dz_j);
    const Angular_Frame frame_i =
        Build_Angular_Frame(dx_i, dy_i, dz_i, distance_i);
    const Angular_Frame frame_j =
        Build_Angular_Frame(dx_j, dy_j, dz_j, distance_j);
    const double sharpness_i = 2.0 * alpha_i * radius * distance_i;
    const double sharpness_j = 2.0 * alpha_j * radius * distance_j;
    const double phi_weight = 2.0 * PI / static_cast<double>(azimuthal_count);
    for (int iu = 0; iu < polar_count; ++iu)
    {
        double polar_i = 0.0, polar_j = 0.0;
        double polar_weight_i = 0.0, polar_weight_j = 0.0;
        double exponential_i = 0.0, exponential_j = 0.0;
        const bool valid_i = Scaled_Polar_Quadrature_Point(
            iu, sharpness_i, polar_nodes, polar_weights, laguerre_nodes,
            laguerre_weights, polar_i, polar_weight_i, exponential_i);
        const bool valid_j = Scaled_Polar_Quadrature_Point(
            iu, sharpness_j, polar_nodes, polar_weights, laguerre_nodes,
            laguerre_weights, polar_j, polar_weight_j, exponential_j);
        if (!valid_i && !valid_j) continue;
        for (int iphi = 0; iphi < azimuthal_count; ++iphi)
        {
            double basis_i = 0.0, basis_j = 0.0;
            double basis_derivative_i[3] = {};
            double basis_derivative_j[3] = {};
            double omega_i_x = 0.0, omega_i_y = 0.0, omega_i_z = 0.0;
            double omega_j_x = 0.0, omega_j_y = 0.0, omega_j_z = 0.0;
            if (valid_i)
            {
                Rotate_Angular_Direction(
                    frame_i, polar_i, cosine_phi[iphi], sine_phi[iphi],
                    omega_i_x, omega_i_y, omega_i_z);
                Scaled_Primitive_Basis_And_Derivative(
                    radius, omega_i_x, omega_i_y, omega_i_z, lx_i, ly_i, lz_i,
                    alpha_i, dx_i, dy_i, dz_i, exponential_i, basis_i,
                    basis_derivative_i);
            }
            if (valid_j)
            {
                Rotate_Angular_Direction(
                    frame_j, polar_j, cosine_phi[iphi], sine_phi[iphi],
                    omega_j_x, omega_j_y, omega_j_z);
                Scaled_Primitive_Basis_And_Derivative(
                    radius, omega_j_x, omega_j_y, omega_j_z, lx_j, ly_j, lz_j,
                    alpha_j, dx_j, dy_j, dz_j, exponential_j, basis_j,
                    basis_derivative_j);
            }
            for (int m = 0; m < component_count; ++m)
            {
                if (valid_i)
                {
                    const double harmonic_i = Real_Spherical_Harmonic(
                        projector_l, m, omega_i_x, omega_i_y, omega_i_z);
                    const double weight_i = polar_weight_i * phi_weight;
                    Kahan_Add(weight_i * basis_i * harmonic_i,
                              projection_i[m], correction_i[m]);
                    for (int direction = 0; direction < 3; ++direction)
                        Kahan_Add(
                            weight_i * basis_derivative_i[direction] *
                                harmonic_i,
                            derivative_i[direction][m],
                            derivative_correction_i[direction][m]);
                }
                if (valid_j)
                {
                    const double harmonic_j = Real_Spherical_Harmonic(
                        projector_l, m, omega_j_x, omega_j_y, omega_j_z);
                    const double weight_j = polar_weight_j * phi_weight;
                    Kahan_Add(weight_j * basis_j * harmonic_j,
                              projection_j[m], correction_j[m]);
                    for (int direction = 0; direction < 3; ++direction)
                        Kahan_Add(
                            weight_j * basis_derivative_j[direction] *
                                harmonic_j,
                            derivative_j[direction][m],
                            derivative_correction_j[direction][m]);
                }
            }
        }
    }
}

QC_ECP_HOST_DEVICE inline double Semilocal_Primitive_Quadrature_Order(
    double alpha_i, double alpha_j, double dx_i, double dy_i, double dz_i,
    double dx_j, double dy_j, double dz_j, double zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j, int projector_l,
    int radial_count, int polar_count, int azimuthal_count)
{
    double radial_nodes[QC_ECP_MAX_QUADRATURE_RADIAL];
    double radial_weights[QC_ECP_MAX_QUADRATURE_RADIAL];
    double polar_nodes[QC_ECP_MAX_QUADRATURE_POLAR];
    double polar_weights[QC_ECP_MAX_QUADRATURE_POLAR];
    double laguerre_nodes[QC_ECP_MAX_QUADRATURE_POLAR];
    double laguerre_weights[QC_ECP_MAX_QUADRATURE_POLAR];
    double cosine_phi[QC_ECP_MAX_QUADRATURE_AZIMUTHAL];
    double sine_phi[QC_ECP_MAX_QUADRATURE_AZIMUTHAL];
    const double distance_i = sqrt(dx_i * dx_i + dy_i * dy_i + dz_i * dz_i);
    const double distance_j = sqrt(dx_j * dx_j + dy_j * dy_j + dz_j * dz_j);
    const double eta = alpha_i + alpha_j + zeta;
    const double radial_center =
        (alpha_i * distance_i + alpha_j * distance_j) / eta;
    const double radial_half_width = 12.0 / sqrt(eta);
    const double lower = radial_center > radial_half_width
                             ? radial_center - radial_half_width
                             : 0.0;
    const double upper = radial_center + radial_half_width;
    Gauss_Legendre_Rule(radial_count, lower, upper, radial_nodes,
                        radial_weights);
    Gauss_Legendre_Rule(polar_count, -1.0, 1.0, polar_nodes, polar_weights);
    Gauss_Laguerre_Rule(polar_count, laguerre_nodes, laguerre_weights);
    for (int i = 0; i < azimuthal_count; ++i)
    {
        const double phi = 2.0 * PI * (static_cast<double>(i) + 0.5) /
                           static_cast<double>(azimuthal_count);
        cosine_phi[i] = cos(phi);
        sine_phi[i] = sin(phi);
    }
    double minimum_exponent =
        alpha_i * distance_i * distance_i +
        alpha_j * distance_j * distance_j -
        (alpha_i * distance_i + alpha_j * distance_j) *
            (alpha_i * distance_i + alpha_j * distance_j) / eta;
    if (minimum_exponent < 0.0) minimum_exponent = 0.0;
    const double scale = exp(-minimum_exponent);
    double value = 0.0, correction = 0.0;
    for (int ir = 0; ir < radial_count; ++ir)
    {
        const double radius = radial_nodes[ir];
        double bra[QC_ECP_MAX_PROJECTOR_COMPONENTS];
        double ket[QC_ECP_MAX_PROJECTOR_COMPONENTS];
        Angular_Projection_Quadrature(
            radius, projector_l, lx_i, ly_i, lz_i, alpha_i, dx_i, dy_i, dz_i,
            polar_count, polar_nodes, polar_weights, laguerre_nodes,
            laguerre_weights, azimuthal_count, cosine_phi, sine_phi, bra);
        Angular_Projection_Quadrature(
            radius, projector_l, lx_j, ly_j, lz_j, alpha_j, dx_j, dy_j, dz_j,
            polar_count, polar_nodes, polar_weights, laguerre_nodes,
            laguerre_weights, azimuthal_count, cosine_phi, sine_phi, ket);
        double angular = 0.0;
        for (int m = 0; m < 2 * projector_l + 1; ++m)
            angular += bra[m] * ket[m];
        const double displacement = radius - radial_center;
        const double radial =
            Integer_Power(radius, n_k) * exp(-eta * displacement * displacement);
        Kahan_Add(radial_weights[ir] * scale * radial * angular, value,
                  correction);
    }
    return value;
}

QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_DERIVATIVE_RESULT
Semilocal_Primitive_Derivative_Quadrature_Order(
    double alpha_i, double alpha_j, double dx_i, double dy_i, double dz_i,
    double dx_j, double dy_j, double dz_j, double zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j, int projector_l,
    int radial_count, int polar_count, int azimuthal_count)
{
    QC_ECP_PRIMITIVE_DERIVATIVE_RESULT result;
    double radial_nodes[QC_ECP_MAX_QUADRATURE_RADIAL];
    double radial_weights[QC_ECP_MAX_QUADRATURE_RADIAL];
    double polar_nodes[QC_ECP_MAX_QUADRATURE_POLAR];
    double polar_weights[QC_ECP_MAX_QUADRATURE_POLAR];
    double laguerre_nodes[QC_ECP_MAX_QUADRATURE_POLAR];
    double laguerre_weights[QC_ECP_MAX_QUADRATURE_POLAR];
    double cosine_phi[QC_ECP_MAX_QUADRATURE_AZIMUTHAL];
    double sine_phi[QC_ECP_MAX_QUADRATURE_AZIMUTHAL];
    const double distance_i = sqrt(dx_i * dx_i + dy_i * dy_i + dz_i * dz_i);
    const double distance_j = sqrt(dx_j * dx_j + dy_j * dy_j + dz_j * dz_j);
    const double eta = alpha_i + alpha_j + zeta;
    const double radial_center =
        (alpha_i * distance_i + alpha_j * distance_j) / eta;
    const double radial_half_width = 12.0 / sqrt(eta);
    const double lower = radial_center > radial_half_width
                             ? radial_center - radial_half_width
                             : 0.0;
    const double upper = radial_center + radial_half_width;
    Gauss_Legendre_Rule(radial_count, lower, upper, radial_nodes,
                        radial_weights);
    Gauss_Legendre_Rule(polar_count, -1.0, 1.0, polar_nodes, polar_weights);
    Gauss_Laguerre_Rule(polar_count, laguerre_nodes, laguerre_weights);
    for (int i = 0; i < azimuthal_count; ++i)
    {
        const double phi = 2.0 * PI * (static_cast<double>(i) + 0.5) /
                           static_cast<double>(azimuthal_count);
        cosine_phi[i] = cos(phi);
        sine_phi[i] = sin(phi);
    }
    double minimum_exponent =
        alpha_i * distance_i * distance_i +
        alpha_j * distance_j * distance_j -
        (alpha_i * distance_i + alpha_j * distance_j) *
            (alpha_i * distance_i + alpha_j * distance_j) / eta;
    if (minimum_exponent < 0.0) minimum_exponent = 0.0;
    const double scale = exp(-minimum_exponent);
    double value_correction = 0.0;
    double derivative_i_correction[3] = {};
    double derivative_j_correction[3] = {};
    for (int ir = 0; ir < radial_count; ++ir)
    {
        const double radius = radial_nodes[ir];
        double bra[QC_ECP_MAX_PROJECTOR_COMPONENTS];
        double ket[QC_ECP_MAX_PROJECTOR_COMPONENTS];
        double derivative_bra[3][QC_ECP_MAX_PROJECTOR_COMPONENTS];
        double derivative_ket[3][QC_ECP_MAX_PROJECTOR_COMPONENTS];
        Angular_Projection_Pair_Derivative_Quadrature(
            radius, projector_l, lx_i, ly_i, lz_i, alpha_i, dx_i, dy_i, dz_i,
            lx_j, ly_j, lz_j, alpha_j, dx_j, dy_j, dz_j, polar_count,
            polar_nodes, polar_weights, laguerre_nodes, laguerre_weights,
            azimuthal_count, cosine_phi, sine_phi, bra, derivative_bra, ket,
            derivative_ket);
        double angular = 0.0, angular_correction = 0.0;
        double angular_i[3] = {}, angular_j[3] = {};
        double angular_i_correction[3] = {}, angular_j_correction[3] = {};
        for (int m = 0; m < 2 * projector_l + 1; ++m)
        {
            Kahan_Add(bra[m] * ket[m], angular, angular_correction);
            for (int direction = 0; direction < 3; ++direction)
            {
                Kahan_Add(derivative_bra[direction][m] * ket[m],
                          angular_i[direction],
                          angular_i_correction[direction]);
                Kahan_Add(bra[m] * derivative_ket[direction][m],
                          angular_j[direction],
                          angular_j_correction[direction]);
            }
        }
        const double displacement = radius - radial_center;
        const double radial = radial_weights[ir] * scale *
                              Integer_Power(radius, n_k) *
                              exp(-eta * displacement * displacement);
        Kahan_Add(radial * angular, result.value, value_correction);
        for (int direction = 0; direction < 3; ++direction)
        {
            Kahan_Add(radial * angular_i[direction],
                      result.derivative_i[direction],
                      derivative_i_correction[direction]);
            Kahan_Add(radial * angular_j[direction],
                      result.derivative_j[direction],
                      derivative_j_correction[direction]);
        }
    }
    result.converged = true;
    return result;
}

QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_DERIVATIVE_RESULT
Semilocal_Primitive_Derivative_Quadrature(
    double alpha_i, double alpha_j, double Ax, double Ay, double Az,
    double Bx, double By, double Bz, double Cx, double Cy, double Cz,
    double zeta, int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j,
    int lz_j, int projector_l)
{
    QC_ECP_PRIMITIVE_DERIVATIVE_RESULT result;
    if (n_k < 0 || projector_l < 0 ||
        projector_l > QC_ECP_MAX_SEMILOCAL_L)
        return result;
    const double dx_i = Cx - Ax, dy_i = Cy - Ay, dz_i = Cz - Az;
    const double dx_j = Cx - Bx, dy_j = Cy - By, dz_j = Cz - Bz;
    const int angular_i = lx_i + ly_i + lz_i + projector_l + 1;
    const int angular_j = lx_j + ly_j + lz_j + projector_l + 1;
    const int angular_degree = angular_i > angular_j ? angular_i : angular_j;
    int fine_polar = 28 + 2 * angular_degree;
    if (fine_polar & 1) ++fine_polar;
    if (fine_polar > QC_ECP_MAX_QUADRATURE_POLAR)
    {
        result.estimated_error = DBL_MAX;
        return result;
    }
    const int coarse_polar = fine_polar - 8;
    const int fine_azimuthal = 2 * fine_polar;
    const int coarse_azimuthal = 2 * coarse_polar;
    const QC_ECP_PRIMITIVE_DERIVATIVE_RESULT coarse =
        Semilocal_Primitive_Derivative_Quadrature_Order(
            alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
            lx_i, ly_i, lz_i, lx_j, ly_j, lz_j, projector_l,
            QC_ECP_SEMILOCAL_QUADRATURE_RADIAL_COARSE,
            coarse_polar, coarse_azimuthal);
    const QC_ECP_PRIMITIVE_DERIVATIVE_RESULT fine =
        Semilocal_Primitive_Derivative_Quadrature_Order(
            alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
            lx_i, ly_i, lz_i, lx_j, ly_j, lz_j, projector_l,
            QC_ECP_SEMILOCAL_QUADRATURE_RADIAL, fine_polar, fine_azimuthal);
    result = fine;
    // This entry point certifies derivatives.  The matrix/value-only path has
    // its own coarse/fine comparison below, so a value component which is not
    // consumed by the gradient must not invalidate six otherwise converged
    // derivative components.
    result.converged = true;
    result.estimated_error = 0.0;
    for (int direction = 0; direction < 3; ++direction)
    {
        double difference =
            Abs(fine.derivative_i[direction] - coarse.derivative_i[direction]);
        result.derivative_i_error[direction] = difference;
        if (difference > result.estimated_error)
            result.estimated_error = difference;
        if (difference >
            2.0e-9 * (1.0 + Abs(fine.derivative_i[direction])))
            result.converged = false;
        difference =
            Abs(fine.derivative_j[direction] - coarse.derivative_j[direction]);
        result.derivative_j_error[direction] = difference;
        if (difference > result.estimated_error)
            result.estimated_error = difference;
        if (difference >
            2.0e-9 * (1.0 + Abs(fine.derivative_j[direction])))
            result.converged = false;
    }
    return result;
}

// Local ECP terms contain one angular integral over chi_i*chi_j.  Keeping the
// product intact is essential: two independent l=0 projections would discard
// angular components which couple to a scalar local operator.
QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_DERIVATIVE_RESULT
Local_Primitive_Derivative_Quadrature_Order(
    double alpha_i, double alpha_j, double dx_i, double dy_i, double dz_i,
    double dx_j, double dy_j, double dz_j, double zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j, int radial_count,
    int polar_count, int azimuthal_count)
{
    QC_ECP_PRIMITIVE_DERIVATIVE_RESULT result;
    double radial_nodes[QC_ECP_MAX_QUADRATURE_RADIAL];
    double radial_weights[QC_ECP_MAX_QUADRATURE_RADIAL];
    double polar_nodes[QC_ECP_MAX_QUADRATURE_POLAR];
    double polar_weights[QC_ECP_MAX_QUADRATURE_POLAR];
    double laguerre_nodes[QC_ECP_MAX_QUADRATURE_POLAR];
    double laguerre_weights[QC_ECP_MAX_QUADRATURE_POLAR];
    double cosine_phi[QC_ECP_MAX_QUADRATURE_AZIMUTHAL];
    double sine_phi[QC_ECP_MAX_QUADRATURE_AZIMUTHAL];
    const double wx = alpha_i * dx_i + alpha_j * dx_j;
    const double wy = alpha_i * dy_i + alpha_j * dy_j;
    const double wz = alpha_i * dz_i + alpha_j * dz_j;
    const double w_norm = sqrt(wx * wx + wy * wy + wz * wz);
    const Angular_Frame frame = Build_Angular_Frame(wx, wy, wz, w_norm);
    const double eta = alpha_i + alpha_j + zeta;
    const double radial_center = w_norm / eta;
    const int radial_power =
        n_k + lx_i + ly_i + lz_i + lx_j + ly_j + lz_j + 1;
    const double radial_peak =
        0.5 * (radial_center +
               sqrt(radial_center * radial_center +
                    2.0 * static_cast<double>(radial_power) / eta));
    const double radial_half_width = 12.0 / sqrt(eta);
    const double lower = radial_center > radial_half_width
                             ? radial_center - radial_half_width
                             : 0.0;
    const double upper = radial_peak + radial_half_width;
    Gauss_Legendre_Rule(radial_count, lower, upper, radial_nodes,
                        radial_weights);
    Gauss_Legendre_Rule(polar_count, -1.0, 1.0, polar_nodes, polar_weights);
    Gauss_Laguerre_Rule(polar_count, laguerre_nodes, laguerre_weights);
    for (int i = 0; i < azimuthal_count; ++i)
    {
        const double phi = 2.0 * PI * (static_cast<double>(i) + 0.5) /
                           static_cast<double>(azimuthal_count);
        cosine_phi[i] = cos(phi);
        sine_phi[i] = sin(phi);
    }
    double minimum_exponent =
        alpha_i * (dx_i * dx_i + dy_i * dy_i + dz_i * dz_i) +
        alpha_j * (dx_j * dx_j + dy_j * dy_j + dz_j * dz_j) -
        w_norm * w_norm / eta;
    if (minimum_exponent < 0.0) minimum_exponent = 0.0;
    const double scale = exp(-minimum_exponent);
    const double phi_weight = 2.0 * PI / static_cast<double>(azimuthal_count);
    double value_correction = 0.0;
    double derivative_i_correction[3] = {};
    double derivative_j_correction[3] = {};
    for (int ir = 0; ir < radial_count; ++ir)
    {
        const double radius = radial_nodes[ir];
        const double angular_sharpness = 2.0 * radius * w_norm;
        double angular = 0.0, angular_correction = 0.0;
        double angular_i[3] = {}, angular_j[3] = {};
        double angular_i_correction[3] = {}, angular_j_correction[3] = {};
        for (int iu = 0; iu < polar_count; ++iu)
        {
            double polar = 0.0, polar_weight = 0.0;
            double angular_exponential = 0.0;
            if (!Scaled_Polar_Quadrature_Point(
                    iu, angular_sharpness, polar_nodes, polar_weights,
                    laguerre_nodes, laguerre_weights, polar, polar_weight,
                    angular_exponential))
                continue;
            for (int iphi = 0; iphi < azimuthal_count; ++iphi)
            {
                double omega_x = 0.0, omega_y = 0.0, omega_z = 0.0;
                Rotate_Angular_Direction(frame, polar, cosine_phi[iphi],
                                         sine_phi[iphi], omega_x, omega_y,
                                         omega_z);
                const double yix = dx_i + radius * omega_x;
                const double yiy = dy_i + radius * omega_y;
                const double yiz = dz_i + radius * omega_z;
                const double yjx = dx_j + radius * omega_x;
                const double yjy = dy_j + radius * omega_y;
                const double yjz = dz_j + radius * omega_z;
                const double pix = Integer_Power(yix, lx_i);
                const double piy = Integer_Power(yiy, ly_i);
                const double piz = Integer_Power(yiz, lz_i);
                const double pjx = Integer_Power(yjx, lx_j);
                const double pjy = Integer_Power(yjy, ly_j);
                const double pjz = Integer_Power(yjz, lz_j);
                const double polynomial_i = pix * piy * piz;
                const double polynomial_j = pjx * pjy * pjz;
                const double exponential = angular_exponential;
                const double basis = polynomial_i * polynomial_j * exponential;
                const double derivative_basis_i[3] = {
                    (2.0 * alpha_i * yix * polynomial_i -
                     static_cast<double>(lx_i) *
                         (lx_i > 0
                              ? Integer_Power(yix, lx_i - 1) * piy * piz
                              : 0.0)) *
                        polynomial_j * exponential,
                    (2.0 * alpha_i * yiy * polynomial_i -
                     static_cast<double>(ly_i) *
                         (ly_i > 0
                              ? pix * Integer_Power(yiy, ly_i - 1) * piz
                              : 0.0)) *
                        polynomial_j * exponential,
                    (2.0 * alpha_i * yiz * polynomial_i -
                     static_cast<double>(lz_i) *
                         (lz_i > 0
                              ? pix * piy * Integer_Power(yiz, lz_i - 1)
                              : 0.0)) *
                        polynomial_j * exponential};
                const double derivative_basis_j[3] = {
                    polynomial_i *
                        (2.0 * alpha_j * yjx * polynomial_j -
                         static_cast<double>(lx_j) *
                             (lx_j > 0
                                  ? Integer_Power(yjx, lx_j - 1) * pjy * pjz
                                  : 0.0)) *
                        exponential,
                    polynomial_i *
                        (2.0 * alpha_j * yjy * polynomial_j -
                         static_cast<double>(ly_j) *
                             (ly_j > 0
                                  ? pjx * Integer_Power(yjy, ly_j - 1) * pjz
                                  : 0.0)) *
                        exponential,
                    polynomial_i *
                        (2.0 * alpha_j * yjz * polynomial_j -
                         static_cast<double>(lz_j) *
                             (lz_j > 0
                                  ? pjx * pjy * Integer_Power(yjz, lz_j - 1)
                                  : 0.0)) *
                        exponential};
                const double weight = polar_weight * phi_weight;
                Kahan_Add(weight * basis, angular, angular_correction);
                for (int direction = 0; direction < 3; ++direction)
                {
                    Kahan_Add(weight * derivative_basis_i[direction],
                              angular_i[direction],
                              angular_i_correction[direction]);
                    Kahan_Add(weight * derivative_basis_j[direction],
                              angular_j[direction],
                              angular_j_correction[direction]);
                }
            }
        }
        const double displacement = radius - radial_center;
        const double radial = radial_weights[ir] * scale *
                              Integer_Power(radius, n_k) *
                              exp(-eta * displacement * displacement);
        Kahan_Add(radial * angular, result.value, value_correction);
        for (int direction = 0; direction < 3; ++direction)
        {
            Kahan_Add(radial * angular_i[direction],
                      result.derivative_i[direction],
                      derivative_i_correction[direction]);
            Kahan_Add(radial * angular_j[direction],
                      result.derivative_j[direction],
                      derivative_j_correction[direction]);
        }
    }
    result.converged = true;
    return result;
}

QC_ECP_HOST_DEVICE inline double Local_Primitive_Quadrature_Order(
    double alpha_i, double alpha_j, double dx_i, double dy_i, double dz_i,
    double dx_j, double dy_j, double dz_j, double zeta, int n_k, int lx_i,
    int ly_i, int lz_i, int lx_j, int ly_j, int lz_j, int radial_count,
    int polar_count, int azimuthal_count)
{
    return Local_Primitive_Derivative_Quadrature_Order(
               alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta,
               n_k, lx_i, ly_i, lz_i, lx_j, ly_j, lz_j, radial_count,
               polar_count, azimuthal_count)
        .value;
}

QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_DERIVATIVE_RESULT
Local_Primitive_Derivative_Quadrature(
    double alpha_i, double alpha_j, double Ax, double Ay, double Az,
    double Bx, double By, double Bz, double Cx, double Cy, double Cz,
    double zeta, int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j,
    int lz_j)
{
    QC_ECP_PRIMITIVE_DERIVATIVE_RESULT result;
    if (n_k < 0) return result;
    const double dx_i = Cx - Ax, dy_i = Cy - Ay, dz_i = Cz - Az;
    const double dx_j = Cx - Bx, dy_j = Cy - By, dz_j = Cz - Bz;
    const int angular_degree =
        lx_i + ly_i + lz_i + lx_j + ly_j + lz_j + 1;
    // The exponential-CDF polar transform removes the exp(-k(1+u)) boundary
    // layer, so resolution is governed by the finite Cartesian polynomial
    // degree rather than by k.  Refuse an unsupported degree instead of
    // capping an estimated order and letting two under-resolved grids agree.
    int fine_polar = 28 + 2 * angular_degree;
    if (fine_polar & 1) ++fine_polar;
    if (fine_polar > QC_ECP_MAX_QUADRATURE_POLAR)
    {
        result.estimated_error = DBL_MAX;
        return result;
    }
    const int coarse_polar = fine_polar - 8;
    const int fine_azimuthal = 2 * fine_polar;
    const int coarse_azimuthal = 2 * coarse_polar;
    const QC_ECP_PRIMITIVE_DERIVATIVE_RESULT coarse =
        Local_Primitive_Derivative_Quadrature_Order(
            alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
            lx_i, ly_i, lz_i, lx_j, ly_j, lz_j,
            QC_ECP_LOCAL_QUADRATURE_RADIAL_COARSE, coarse_polar,
            coarse_azimuthal);
    const QC_ECP_PRIMITIVE_DERIVATIVE_RESULT fine =
        Local_Primitive_Derivative_Quadrature_Order(
            alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
            lx_i, ly_i, lz_i, lx_j, ly_j, lz_j,
            QC_ECP_MAX_QUADRATURE_RADIAL, fine_polar, fine_azimuthal);
    result = fine;
    // See the semi-local analogue above: this routine certifies the six
    // derivative observables, while Local_Primitive_Quadrature performs the
    // independent value-only certification used by the ECP matrix.
    result.converged = true;
    result.estimated_error = 0.0;
    for (int direction = 0; direction < 3; ++direction)
    {
        double difference =
            Abs(fine.derivative_i[direction] - coarse.derivative_i[direction]);
        result.derivative_i_error[direction] = difference;
        if (difference > result.estimated_error)
            result.estimated_error = difference;
        if (difference >
            2.0e-9 * (1.0 + Abs(fine.derivative_i[direction])))
            result.converged = false;
        difference =
            Abs(fine.derivative_j[direction] - coarse.derivative_j[direction]);
        result.derivative_j_error[direction] = difference;
        if (difference > result.estimated_error)
            result.estimated_error = difference;
        if (difference >
            2.0e-9 * (1.0 + Abs(fine.derivative_j[direction])))
            result.converged = false;
    }
    return result;
}

QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_RESULT Local_Primitive_Quadrature(
    double alpha_i, double alpha_j, double Ax, double Ay, double Az,
    double Bx, double By, double Bz, double Cx, double Cy, double Cz,
    double zeta, int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j,
    int lz_j)
{
    QC_ECP_PRIMITIVE_RESULT result;
    if (n_k < 0) return result;
    const double dx_i = Cx - Ax, dy_i = Cy - Ay, dz_i = Cz - Az;
    const double dx_j = Cx - Bx, dy_j = Cy - By, dz_j = Cz - Bz;
    const int angular_degree =
        lx_i + ly_i + lz_i + lx_j + ly_j + lz_j;
    int fine_polar = 28 + 2 * angular_degree;
    if (fine_polar & 1) ++fine_polar;
    if (fine_polar > QC_ECP_MAX_QUADRATURE_POLAR)
    {
        result.last_layer_bound = DBL_MAX;
        return result;
    }
    const int coarse_polar = fine_polar - 8;
    const int fine_azimuthal = 2 * fine_polar;
    const int coarse_azimuthal = 2 * coarse_polar;
    const double coarse = Local_Primitive_Quadrature_Order(
        alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
        lx_i, ly_i, lz_i, lx_j, ly_j, lz_j,
        QC_ECP_LOCAL_QUADRATURE_RADIAL_COARSE, coarse_polar,
        coarse_azimuthal);
    const double fine = Local_Primitive_Quadrature_Order(
        alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
        lx_i, ly_i, lz_i, lx_j, ly_j, lz_j,
        QC_ECP_MAX_QUADRATURE_RADIAL, fine_polar, fine_azimuthal);
    result.value = fine;
    result.last_layer_bound = Abs(fine - coarse);
    result.converged =
        result.last_layer_bound <= 2.0e-9 * (1.0 + Abs(fine));
    return result;
}

QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_RESULT
Semilocal_Primitive_Quadrature(
    double alpha_i, double alpha_j, double Ax, double Ay, double Az,
    double Bx, double By, double Bz, double Cx, double Cy, double Cz,
    double zeta, int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j,
    int lz_j, int projector_l)
{
    QC_ECP_PRIMITIVE_RESULT result;
    if (n_k < 0 || projector_l < 0 ||
        projector_l > QC_ECP_MAX_SEMILOCAL_L)
        return result;
    const double dx_i = Cx - Ax, dy_i = Cy - Ay, dz_i = Cz - Az;
    const double dx_j = Cx - Bx, dy_j = Cy - By, dz_j = Cz - Bz;
    const int angular_i = lx_i + ly_i + lz_i + projector_l;
    const int angular_j = lx_j + ly_j + lz_j + projector_l;
    const int angular_degree = angular_i > angular_j ? angular_i : angular_j;
    int fine_polar = 28 + 2 * angular_degree;
    if (fine_polar & 1) ++fine_polar;
    if (fine_polar > QC_ECP_MAX_QUADRATURE_POLAR)
    {
        result.last_layer_bound = DBL_MAX;
        return result;
    }
    const int coarse_polar = fine_polar - 8;
    const int fine_azimuthal = 2 * fine_polar;
    const int coarse_azimuthal = 2 * coarse_polar;
    const double coarse = Semilocal_Primitive_Quadrature_Order(
        alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
        lx_i, ly_i, lz_i, lx_j, ly_j, lz_j, projector_l,
        QC_ECP_SEMILOCAL_QUADRATURE_RADIAL_COARSE, coarse_polar,
        coarse_azimuthal);
    const double fine = Semilocal_Primitive_Quadrature_Order(
        alpha_i, alpha_j, dx_i, dy_i, dz_i, dx_j, dy_j, dz_j, zeta, n_k,
        lx_i, ly_i, lz_i, lx_j, ly_j, lz_j, projector_l,
        QC_ECP_SEMILOCAL_QUADRATURE_RADIAL, fine_polar, fine_azimuthal);
    result.value = fine;
    result.last_layer_bound = Abs(fine - coarse);
    result.converged =
        result.last_layer_bound <= 2.0e-9 * (1.0 + Abs(fine));
    return result;
}

QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_RESULT Semilocal_Primitive(
    double alpha_i, double alpha_j, double Ax, double Ay, double Az,
    double Bx, double By, double Bz, double Cx, double Cy, double Cz,
    double zeta, int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j,
    int lz_j, int projector_l,
    int max_order = QC_ECP_MAX_SERIES_ORDER)
{
    QC_ECP_PRIMITIVE_RESULT result;
    if (n_k < 0 || projector_l < 0 ||
        projector_l > QC_ECP_MAX_SEMILOCAL_L || max_order < projector_l ||
        max_order > QC_ECP_MAX_SERIES_ORDER)
        return result;
    if (((max_order - projector_l) & 1) != 0) --max_order;
    const double dx_i = Cx - Ax, dy_i = Cy - Ay, dz_i = Cz - Az;
    const double dx_j = Cx - Bx, dy_j = Cy - By, dz_j = Cz - Bz;
    const double distance_i2 = dx_i * dx_i + dy_i * dy_i + dz_i * dz_i;
    const double distance_j2 = dx_j * dx_j + dy_j * dy_j + dz_j * dz_j;
    const double scale = exp(-alpha_i * distance_i2 -
                             alpha_j * distance_j2);
    const double eta = alpha_i + alpha_j + zeta;
    double bra[QC_ECP_MAX_PROJECTOR_COMPONENTS]
               [QC_ECP_MAX_SERIES_LAYERS];
    double ket[QC_ECP_MAX_PROJECTOR_COMPONENTS]
              [QC_ECP_MAX_SERIES_LAYERS];
    Build_Projection_Series(projector_l, lx_i, ly_i, lz_i, alpha_i, dx_i,
                            dy_i, dz_i, max_order, bra);
    Build_Projection_Series(projector_l, lx_j, ly_j, lz_j, alpha_j, dx_j,
                            dy_j, dz_j, max_order, ket);
    const int component_count = 2 * projector_l + 1;
    const int layer_count = (max_order - projector_l) / 2 + 1;
    double layers[QC_ECP_MAX_SERIES_LAYERS] = {};
    double value = 0.0, correction = 0.0;
    for (int outer = 0; outer < layer_count; ++outer)
    {
        double layer_value = 0.0, layer_correction = 0.0;
        for (int left = 0; left <= outer; ++left)
            for (int right = 0; right <= outer; ++right)
            {
                if (left != outer && right != outer) continue;
                const int q = projector_l + 2 * left;
                const int s = projector_l + 2 * right;
                double angular = 0.0;
                for (int m = 0; m < component_count; ++m)
                    angular += bra[m][left] * ket[m][right];
                Kahan_Add(scale * angular *
                              Gaussian_Radial_Moment(n_k + q + s, eta),
                          layer_value, layer_correction);
            }
        layers[outer] = layer_value;
        Kahan_Add(layer_value, value, correction);
    }
    result.value = value;
    result.converged =
        Tail_Is_Converged(layers, layer_count, value,
                          result.last_layer_bound);
    return result;
}

QC_ECP_HOST_DEVICE inline QC_ECP_PRIMITIVE_RESULT Local_Primitive_Series(
    double alpha_i, double alpha_j, double Ax, double Ay, double Az,
    double Bx, double By, double Bz, double Cx, double Cy, double Cz,
    double zeta, int n_k, int lx_i, int ly_i, int lz_i, int lx_j, int ly_j,
    int lz_j, int max_order = QC_ECP_MAX_SERIES_ORDER)
{
    QC_ECP_PRIMITIVE_RESULT result;
    if (n_k < 0 || max_order < 0 || max_order > QC_ECP_MAX_SERIES_ORDER)
        return result;
    if (max_order & 1) --max_order;
    const double dx_i = Cx - Ax, dy_i = Cy - Ay, dz_i = Cz - Az;
    const double dx_j = Cx - Bx, dy_j = Cy - By, dz_j = Cz - Bz;
    double angular[QC_ECP_MAX_SERIES_LAYERS] = {};
    if (!Build_Local_Angular_Series(
            lx_i, ly_i, lz_i, alpha_i, dx_i, dy_i, dz_i, lx_j, ly_j, lz_j,
            alpha_j, dx_j, dy_j, dz_j, max_order, angular))
        return result;
    const double scale =
        exp(-alpha_i * (dx_i * dx_i + dy_i * dy_i + dz_i * dz_i) -
            alpha_j * (dx_j * dx_j + dy_j * dy_j + dz_j * dz_j));
    const double eta = alpha_i + alpha_j + zeta;
    const int layer_count = max_order / 2 + 1;
    double layers[QC_ECP_MAX_SERIES_LAYERS] = {};
    double value = 0.0, correction = 0.0;
    for (int layer = 0; layer < layer_count; ++layer)
    {
        layers[layer] =
            scale * angular[layer] *
            Gaussian_Radial_Moment(n_k + 2 * layer, eta);
        Kahan_Add(layers[layer], value, correction);
    }
    result.value = value;
    result.converged =
        Tail_Is_Converged(layers, layer_count, value,
                          result.last_layer_bound);
    return result;
}
}  // namespace qc_ecp_math

// ECP 参数集 (类比 QC_BASIS_SET)
using QC_ECP_MAP = std::map<std::string, QC_ECP_ATOM_DATA>;

struct QC_ECP_SET
{
    virtual ~QC_ECP_SET() = default;
    virtual void Initialize() = 0;
    const char* name = "";
    QC_ECP_MAP data;
};

// Helper to build ECP channels from inline data
inline QC_ECP_CHANNEL make_channel(int l,
                                   std::initializer_list<QC_ECP_TERM> terms)
{
    QC_ECP_CHANNEL ch;
    ch.l = l;
    ch.terms.assign(terms.begin(), terms.end());
    return ch;
}
