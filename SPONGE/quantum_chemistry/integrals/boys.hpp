#pragma once

// Boys function
//
//     F_m(t) = integral_0^1 x^(2m) exp(-t x^2) dx,  t >= 0.
//
// Upward recursion is inexpensive, but loses accuracy when its subtraction
// is ill-conditioned.  Rather than selecting it at a fixed value of t, track
// a conservative forward-error bound through every requested order.  If that
// bound exceeds the target, evaluate the highest order from the positive
// series
//
//   F_m(t) = exp(-t) sum_k (2t)^k /
//            ((2m+1)(2m+3)...(2m+2k+1))
//
// and recurse downward.  Once the series term ratios are below one they are
// monotonically decreasing, so term*r/(1-r) is an explicit upper bound on
// the uncomputed tail.  The downward recurrence is a positive sum and cannot
// amplify relative error from the highest-order seed.
static __device__ __forceinline__ void QC_Compute_Boys_Double(double* F,
                                                              float t,
                                                              int max_m)
{
    const double td = (double)t;

    // Invalid physical arguments are made visible to the caller as NaNs.
    if (!(td >= 0.0))
    {
        const double invalid = sqrt(td);
        for (int m = 0; m <= max_m; m++) F[m] = invalid;
        return;
    }
    if (td == 0.0)
    {
        for (int m = 0; m <= max_m; m++) F[m] = 1.0 / (2.0 * m + 1.0);
        return;
    }

    // t arrives as a float, so a value beyond the largest finite float can
    // only be +infinity.  The exact limiting value is zero for every order.
    if (td > 3.40282346638528859812e38)
    {
        for (int m = 0; m <= max_m; m++) F[m] = 0.0;
        return;
    }

    const double epsilon = 2.22044604925031308085e-16;
    const double target_relative_error = 64.0 * epsilon;
    const double series_relative_tolerance = 8.0 * epsilon;
    const double exp_t = exp(-td);
    const double sqrt_t = sqrt(td);
    const double f0 = 0.88622692545275801365 * erf(sqrt_t) / sqrt_t;

    // First try the exact upward recurrence while carrying a conservative
    // relative roundoff bound.  Four epsilons cover erf, multiplication and
    // division in the F0 expression.
    F[0] = f0;
    double previous = f0;
    double relative_error_bound = 4.0 * epsilon;
    bool upward_is_accurate = true;
    for (int m = 0; m < max_m; m++)
    {
        const double scaled_previous = (2.0 * m + 1.0) * previous;
        const double numerator = scaled_previous - exp_t;

        // With exp(-t) already underflowed, a zero previous value means this
        // and every higher order are below the double-precision range.
        if (numerator == 0.0 && exp_t == 0.0 && previous == 0.0)
        {
            for (int order = m + 1; order <= max_m; order++) F[order] = 0.0;
            return;
        }
        if (!(numerator > 0.0))
        {
            upward_is_accurate = false;
            break;
        }

        const double next = numerator / (2.0 * td);
        F[m + 1] = next;

        // Absolute errors in the product and exp(-t), plus one rounding for
        // the subtraction and one for the division, give this propagated
        // relative bound.  Stop immediately once the requested target cannot
        // be certified; the positive series below then supplies the seed.
        relative_error_bound =
            (fabs(scaled_previous) * relative_error_bound +
             2.0 * epsilon * (fabs(scaled_previous) + exp_t)) /
                fabs(numerator) +
            epsilon;
        if (relative_error_bound > target_relative_error)
        {
            upward_is_accurate = false;
            break;
        }
        previous = next;
    }
    if (upward_is_accurate) return;

    // Evaluate only F_max from the positive series.  Kahan summation keeps
    // accumulation error below the analytically bounded truncation error.
    const int order = max_m;
    double term = exp_t / (2.0 * order + 1.0);
    double sum = term;
    double compensation = 0.0;
    for (int k = 1;; k++)
    {
        const double ratio = (2.0 * td) / (2.0 * order + 2.0 * k + 1.0);
        term *= ratio;

        const double corrected_term = term - compensation;
        const double updated_sum = sum + corrected_term;
        compensation = (updated_sum - sum) - corrected_term;
        sum = updated_sum;

        const double next_ratio = (2.0 * td) / (2.0 * order + 2.0 * k + 3.0);
        if (next_ratio < 1.0)
        {
            const double tail_bound = term * next_ratio / (1.0 - next_ratio);
            if (tail_bound <= series_relative_tolerance * sum) break;
        }
    }

    F[order] = sum;
    for (int m = order - 1; m >= 0; m--)
        F[m] = (2.0 * td * F[m + 1] + exp_t) / (2.0 * m + 1.0);
}
