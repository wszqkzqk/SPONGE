#pragma once

#include <cstdint>
#include <cstring>
#include <limits>

#include "dft.hpp"

// The project is compiled with -ffast-math.  std::isfinite and comparisons
// against NaN are therefore not reliable guards: the optimizer is allowed to
// assume that non-finite values do not exist.  Inspect the IEEE-754 exponent
// bits instead.  This helper is used both on the host and in device kernels.
static inline __host__ __device__ std::uint64_t QC_XC_Double_Bits(double value)
{
    std::uint64_t bits = 0;
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    bits = static_cast<std::uint64_t>(__double_as_longlong(value));
#elif defined(__GNUC__) || defined(__clang__)
    static_assert(sizeof(bits) == sizeof(value) &&
                      std::numeric_limits<double>::is_iec559,
                  "SPONGE requires 64-bit IEEE-754 doubles");
    std::memcpy(&bits, &value, sizeof(bits));
    __asm__ __volatile__("" : "+r"(bits));
#else
    std::memcpy(&bits, &value, sizeof(bits));
#endif
    return bits;
}

static inline __host__ __device__ bool QC_XC_Double_Is_Finite(double value)
{
    const std::uint64_t bits = QC_XC_Double_Bits(value);
    return (bits & UINT64_C(0x7ff0000000000000)) !=
           UINT64_C(0x7ff0000000000000);
}

static inline __host__ __device__ double QC_XC_Double_From_Bits(
    std::uint64_t bits)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
    return __longlong_as_double(static_cast<long long>(bits));
#else
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
#if defined(__GNUC__) || defined(__clang__)
    // Keep -ffast-math from replacing a bit-constructed special value with a
    // value that satisfies its global finite-only assumption.
    __asm__ __volatile__("" : "+m"(value));
#endif
    return value;
#endif
}

static inline __host__ __device__ double QC_XC_Quiet_NaN()
{
    return QC_XC_Double_From_Bits(UINT64_C(0x7ff8000000000000));
}

static inline __host__ __device__ double QC_XC_Positive_Infinity()
{
    return QC_XC_Double_From_Bits(UINT64_C(0x7ff0000000000000));
}

static inline __host__ __device__ bool QC_XC_Method_Is_Supported(
    QC_METHOD method)
{
    return method == QC_METHOD::LDA || method == QC_METHOD::PBE ||
           method == QC_METHOD::BLYP || method == QC_METHOD::PBE0 ||
           method == QC_METHOD::B3LYP;
}

// XC is defined on non-negative densities and gradient invariants.  Vacuum is
// a single boundary point: a zero density cannot carry a non-zero gradient.
// No positive density, however small, is discarded.
static inline __host__ __device__ bool QC_XC_Inputs_Valid_RKS(double rho,
                                                              double sigma)
{
    if (!QC_XC_Double_Is_Finite(rho) ||
        !QC_XC_Double_Is_Finite(sigma))
        return false;
    if (rho < 0.0 || sigma < 0.0) return false;
    return rho != 0.0 || sigma == 0.0;
}

struct QC_XC_Normalized_IEEE_Double
{
    std::uint64_t significand;
    int exponent;
};

struct QC_XC_UInt106
{
    std::uint64_t high;
    std::uint64_t low;
    int exponent;
};

// Decompose a finite, positive IEEE-754 double as significand*2^exponent,
// with the significand normalized to exactly 53 bits.  Subnormals are handled
// by integer shifts, so no floating operation can flush them to zero.
static inline __host__ __device__ QC_XC_Normalized_IEEE_Double
QC_XC_Normalize_Positive_Double_Bits(std::uint64_t bits)
{
    const std::uint64_t fraction =
        bits & UINT64_C(0x000fffffffffffff);
    const unsigned exponent_field =
        static_cast<unsigned>((bits >> 52) & UINT64_C(0x7ff));
    QC_XC_Normalized_IEEE_Double result;
    if (exponent_field != 0u)
    {
        result.significand = fraction | UINT64_C(0x0010000000000000);
        result.exponent = static_cast<int>(exponent_field) - 1075;
        return result;
    }

    result.significand = fraction;
    result.exponent = -1074;
    while ((result.significand & UINT64_C(0x0010000000000000)) == 0u)
    {
        result.significand <<= 1;
        --result.exponent;
    }
    return result;
}

// Exact 53x53 -> 106-bit multiplication using only portable 64-bit integer
// operations.  The two 32-bit cross products are accumulated separately so
// neither limb can overflow before its carry is propagated.
static inline __host__ __device__ QC_XC_UInt106 QC_XC_Multiply_53x53(
    const QC_XC_Normalized_IEEE_Double& a,
    const QC_XC_Normalized_IEEE_Double& b)
{
    const std::uint64_t mask32 = UINT64_C(0xffffffff);
    const std::uint64_t a_low = a.significand & mask32;
    const std::uint64_t a_high = a.significand >> 32;
    const std::uint64_t b_low = b.significand & mask32;
    const std::uint64_t b_high = b.significand >> 32;

    QC_XC_UInt106 result;
    result.low = a_low * b_low;
    result.high = a_high * b_high;

    const std::uint64_t cross0 = a_low * b_high;
    const std::uint64_t add0 = cross0 << 32;
    const std::uint64_t old0 = result.low;
    result.low += add0;
    result.high += (cross0 >> 32) + (result.low < old0 ? 1u : 0u);

    const std::uint64_t cross1 = a_high * b_low;
    const std::uint64_t add1 = cross1 << 32;
    const std::uint64_t old1 = result.low;
    result.low += add1;
    result.high += (cross1 >> 32) + (result.low < old1 ? 1u : 0u);

    result.exponent = a.exponent + b.exponent;
    // A product of two normalized 53-bit integers has bit 104 or 105 set.
    // Put bit 105 at a fixed position so exponent comparison is sufficient.
    if ((result.high & (UINT64_C(1) << 41)) == 0u)
    {
        result.high = (result.high << 1) | (result.low >> 63);
        result.low <<= 1;
        --result.exponent;
    }
    return result;
}

static inline __host__ __device__ bool QC_XC_UInt106_Less_Or_Equal(
    const QC_XC_UInt106& a, const QC_XC_UInt106& b)
{
    if (a.exponent != b.exponent) return a.exponent < b.exponent;
    if (a.high != b.high) return a.high < b.high;
    return a.low <= b.low;
}

static inline __host__ __device__ bool QC_XC_Gradient_Gram_Is_PSD(
    double sigma_aa, double sigma_ab, double sigma_bb)
{
    const std::uint64_t aa_bits =
        QC_XC_Double_Bits(sigma_aa) & UINT64_C(0x7fffffffffffffff);
    const std::uint64_t ab_bits =
        QC_XC_Double_Bits(sigma_ab) & UINT64_C(0x7fffffffffffffff);
    const std::uint64_t bb_bits =
        QC_XC_Double_Bits(sigma_bb) & UINT64_C(0x7fffffffffffffff);
    if (ab_bits == 0u) return true;
    if (aa_bits == 0u || bb_bits == 0u) return false;

    const QC_XC_Normalized_IEEE_Double aa =
        QC_XC_Normalize_Positive_Double_Bits(aa_bits);
    const QC_XC_Normalized_IEEE_Double ab =
        QC_XC_Normalize_Positive_Double_Bits(ab_bits);
    const QC_XC_Normalized_IEEE_Double bb =
        QC_XC_Normalize_Positive_Double_Bits(bb_bits);
    const QC_XC_UInt106 ab_squared = QC_XC_Multiply_53x53(ab, ab);
    const QC_XC_UInt106 aa_times_bb = QC_XC_Multiply_53x53(aa, bb);
    return QC_XC_UInt106_Less_Or_Equal(ab_squared, aa_times_bb);
}

struct QC_XC_Interval
{
    double lower;
    double upper;
};

static inline __host__ __device__ bool QC_XC_Double_Is_Zero(double value)
{
    return (QC_XC_Double_Bits(value) & UINT64_C(0x7fffffffffffffff)) == 0u;
}

static inline __host__ __device__ double QC_XC_Next_Up(double value)
{
    std::uint64_t bits = QC_XC_Double_Bits(value);
    const std::uint64_t magnitude = bits & UINT64_C(0x7fffffffffffffff);
    if (magnitude > UINT64_C(0x7ff0000000000000)) return value;
    if (bits == UINT64_C(0x7ff0000000000000)) return value;
    if (bits == UINT64_C(0xfff0000000000000))
        return QC_XC_Double_From_Bits(UINT64_C(0xffefffffffffffff));
    if (magnitude == 0u) return QC_XC_Double_From_Bits(UINT64_C(1));
    bits += (bits >> 63) != 0u ? UINT64_C(-1) : UINT64_C(1);
    return QC_XC_Double_From_Bits(bits);
}

static inline __host__ __device__ double QC_XC_Next_Down(double value)
{
    std::uint64_t bits = QC_XC_Double_Bits(value);
    const std::uint64_t magnitude = bits & UINT64_C(0x7fffffffffffffff);
    if (magnitude > UINT64_C(0x7ff0000000000000)) return value;
    if (bits == UINT64_C(0xfff0000000000000)) return value;
    if (bits == UINT64_C(0x7ff0000000000000))
        return QC_XC_Double_From_Bits(UINT64_C(0x7fefffffffffffff));
    if (magnitude == 0u)
        return QC_XC_Double_From_Bits(UINT64_C(0x8000000000000001));
    bits += (bits >> 63) != 0u ? UINT64_C(1) : UINT64_C(-1);
    return QC_XC_Double_From_Bits(bits);
}

static inline __host__ __device__ QC_XC_Interval QC_XC_Product_Interval(
    double a, double b)
{
    if (QC_XC_Double_Is_Zero(a) || QC_XC_Double_Is_Zero(b))
        return {0.0, 0.0};
    volatile double rounded = fma(a, b, 0.0);
    if ((QC_XC_Double_Bits(rounded) & UINT64_C(0x7fffffffffffffff)) >
        UINT64_C(0x7ff0000000000000))
    {
        return {QC_XC_Double_From_Bits(UINT64_C(0xfff0000000000000)),
                QC_XC_Double_From_Bits(UINT64_C(0x7ff0000000000000))};
    }
    return {QC_XC_Next_Down(rounded), QC_XC_Next_Up(rounded)};
}

static inline __host__ __device__ bool QC_XC_Interval_Is_Exact_Zero(
    const QC_XC_Interval& interval)
{
    return QC_XC_Double_Is_Zero(interval.lower) &&
           QC_XC_Double_Is_Zero(interval.upper);
}

static inline __host__ __device__ QC_XC_Interval QC_XC_Add_Intervals(
    const QC_XC_Interval& a, const QC_XC_Interval& b)
{
    if (QC_XC_Interval_Is_Exact_Zero(a)) return b;
    if (QC_XC_Interval_Is_Exact_Zero(b)) return a;
    volatile double lower = fma(1.0, a.lower, b.lower);
    volatile double upper = fma(1.0, a.upper, b.upper);
    const std::uint64_t lower_magnitude =
        QC_XC_Double_Bits(lower) & UINT64_C(0x7fffffffffffffff);
    const std::uint64_t upper_magnitude =
        QC_XC_Double_Bits(upper) & UINT64_C(0x7fffffffffffffff);
    if (lower_magnitude > UINT64_C(0x7ff0000000000000) ||
        upper_magnitude > UINT64_C(0x7ff0000000000000))
    {
        return {QC_XC_Double_From_Bits(UINT64_C(0xfff0000000000000)),
                QC_XC_Double_From_Bits(UINT64_C(0x7ff0000000000000))};
    }
    return {QC_XC_Next_Down(lower), QC_XC_Next_Up(upper)};
}

static inline __host__ __device__ QC_XC_Interval QC_XC_Dot3_Interval(
    double ax, double ay, double az, double bx, double by, double bz)
{
    const QC_XC_Interval x = QC_XC_Product_Interval(ax, bx);
    const QC_XC_Interval y = QC_XC_Product_Interval(ay, by);
    const QC_XC_Interval z = QC_XC_Product_Interval(az, bz);
    return QC_XC_Add_Intervals(QC_XC_Add_Intervals(x, y), z);
}

// Construct representable spin-gradient invariants with a proof of PSD.
// Diagonal entries are outward upper bounds on the exact squared norms.  The
// cross entry is the endpoint of an exact-dot interval closest to zero (or
// zero if the interval straddles it).  Cauchy-Schwarz for the exact vectors
// then proves sigma_ab^2 <= sigma_aa*sigma_bb for the published doubles.
static inline __host__ __device__ bool QC_XC_Build_Gradient_Gram(
    double gax, double gay, double gaz, double gbx, double gby, double gbz,
    double& sigma_aa, double& sigma_ab, double& sigma_bb)
{
    if (!QC_XC_Double_Is_Finite(gax) || !QC_XC_Double_Is_Finite(gay) ||
        !QC_XC_Double_Is_Finite(gaz) || !QC_XC_Double_Is_Finite(gbx) ||
        !QC_XC_Double_Is_Finite(gby) || !QC_XC_Double_Is_Finite(gbz))
    {
        sigma_aa = sigma_ab = sigma_bb = QC_XC_Quiet_NaN();
        return false;
    }
    const QC_XC_Interval aa =
        QC_XC_Dot3_Interval(gax, gay, gaz, gax, gay, gaz);
    const QC_XC_Interval ab =
        QC_XC_Dot3_Interval(gax, gay, gaz, gbx, gby, gbz);
    const QC_XC_Interval bb =
        QC_XC_Dot3_Interval(gbx, gby, gbz, gbx, gby, gbz);
    sigma_aa = aa.upper;
    sigma_bb = bb.upper;
    sigma_ab = ab.lower > 0.0 ? ab.lower : (ab.upper < 0.0 ? ab.upper : 0.0);
    return QC_XC_Double_Is_Finite(sigma_aa) &&
           QC_XC_Double_Is_Finite(sigma_ab) &&
           QC_XC_Double_Is_Finite(sigma_bb) &&
           QC_XC_Gradient_Gram_Is_PSD(sigma_aa, sigma_ab, sigma_bb);
}

static inline __host__ __device__ bool QC_XC_Inputs_Valid_UKS(
    double rho_a, double rho_b, double sigma_aa, double sigma_ab,
    double sigma_bb)
{
    if (!QC_XC_Double_Is_Finite(rho_a) ||
        !QC_XC_Double_Is_Finite(rho_b) ||
        !QC_XC_Double_Is_Finite(sigma_aa) ||
        !QC_XC_Double_Is_Finite(sigma_ab) ||
        !QC_XC_Double_Is_Finite(sigma_bb))
        return false;
    if (rho_a < 0.0 || rho_b < 0.0 || sigma_aa < 0.0 || sigma_bb < 0.0)
        return false;

    const double rho = rho_a + rho_b;
    const double sigma = sigma_aa + 2.0 * sigma_ab + sigma_bb;
    if (!QC_XC_Double_Is_Finite(rho) ||
        !QC_XC_Double_Is_Finite(sigma) || sigma < 0.0)
        return false;

    if (rho == 0.0)
        return sigma_aa == 0.0 && sigma_ab == 0.0 && sigma_bb == 0.0;

    // The three invariants are a 2x2 Gram matrix.  Validate positive
    // semidefiniteness exactly in integer arithmetic: floating ratios can
    // round an out-of-domain matrix back onto the boundary.
    if (!QC_XC_Gradient_Gram_Is_PSD(sigma_aa, sigma_ab, sigma_bb)) return false;

    // A non-negative differentiable spin density has zero gradient wherever
    // that spin density itself is zero.  Enforcing this boundary identity
    // avoids asking for partial derivatives outside the admissible domain.
    if (rho_a == 0.0 && (sigma_aa != 0.0 || sigma_ab != 0.0)) return false;
    if (rho_b == 0.0 && (sigma_bb != 0.0 || sigma_ab != 0.0)) return false;
    return true;
}

// Output bits used by the host-side failure plumbing.  At a fully polarized
// PBE/PBE0 point with a non-zero total gradient, dE_c/d(rho_minority) has the
// genuine +infinity limit caused by d[(rho_minority/rho)^(2/3)]/d rho.
// It must not be confused with an accidental floating-point failure.
enum QC_XC_UKS_Output_Bit : unsigned
{
    QC_XC_UKS_VRHO_A_BIT = 1u << 0,
    QC_XC_UKS_VRHO_B_BIT = 1u << 1,
    QC_XC_UKS_VSIGMA_AA_BIT = 1u << 2,
    QC_XC_UKS_VSIGMA_AB_BIT = 1u << 3,
    QC_XC_UKS_VSIGMA_BB_BIT = 1u << 4
};

static inline __host__ __device__ unsigned
QC_XC_UKS_Expected_Infinite_Output_Mask(QC_METHOD method, double rho_a,
                                        double rho_b, double sigma_aa,
                                        double sigma_ab, double sigma_bb)
{
    if (method != QC_METHOD::PBE && method != QC_METHOD::PBE0) return 0u;
    const double sigma = sigma_aa + 2.0 * sigma_ab + sigma_bb;
    if (sigma == 0.0) return 0u;
    if (rho_a == 0.0 && rho_b > 0.0) return QC_XC_UKS_VRHO_A_BIT;
    if (rho_b == 0.0 && rho_a > 0.0) return QC_XC_UKS_VRHO_B_BIT;
    return 0u;
}

// Public energy-only interfaces.  Values and first derivatives deliberately
// share the same implementation so that the SCF energy and potential cannot
// drift apart algebraically.
static inline __host__ __device__ double QC_Local_Exc_Density(
    QC_METHOD method, double rho, double sigma)
{
    double energy, vrho, vsigma;
    QC_VXC_Analytical_RKS(method, rho, sigma, energy, vrho, vsigma);
    return energy;
}

static inline __host__ __device__ double QC_Local_Exc_Density_UKS(
    QC_METHOD method, double rho_a, double rho_b, double sigma_aa,
    double sigma_ab, double sigma_bb)
{
    double energy, vrho_a, vrho_b, vsigma_aa, vsigma_ab, vsigma_bb;
    QC_VXC_Analytical_UKS(method, rho_a, rho_b, sigma_aa, sigma_ab, sigma_bb,
                          energy, vrho_a, vrho_b, vsigma_aa, vsigma_ab,
                          vsigma_bb);
    return energy;
}

#include "xc_deriv.hpp"
