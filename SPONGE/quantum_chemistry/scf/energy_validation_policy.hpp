#pragma once

#include <cstdint>
#include <cstring>
#include <limits>

enum QC_SCF_Energy_Quantity
{
    QC_SCF_CURRENT_ENERGY = 0,
    QC_SCF_DELTA_ENERGY = 1,
};

struct QC_SCF_Energy_Validation_Failure
{
    QC_SCF_Energy_Quantity quantity;
    const char* quantity_name;
    int md_step;
    int iteration;
    double value;
};

static inline const char* QC_SCF_Energy_Quantity_Name(
    QC_SCF_Energy_Quantity quantity)
{
    switch (quantity)
    {
        case QC_SCF_CURRENT_ENERGY:
            return "current energy";
        case QC_SCF_DELTA_ENERGY:
            return "energy delta";
    }
    return "unknown energy quantity";
}

// Inspect the IEEE-754 exponent bits directly.  SCF translation units are
// built with finite-math optimizations, under which std::isfinite and ordinary
// comparisons are not a reliable validation boundary for injected or
// propagated NaN/Inf values.
static inline bool QC_SCF_Energy_Double_Is_Finite(double value)
{
    std::uint64_t bits = 0;
    static_assert(
        sizeof(bits) == sizeof(value) && std::numeric_limits<double>::is_iec559,
        "SPONGE requires 64-bit IEEE-754 doubles");
    std::memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return ((bits >> 52U) & UINT64_C(0x7ff)) != UINT64_C(0x7ff);
}

// Validate the host observation made immediately after SCF energy
// accumulation.  The delta has no meaning before a previous energy exists,
// so callers explicitly state whether it is part of this observation.
//
// The policy is controller-independent by design: production supplies a
// handler which raises its formatted fatal error, while tests inject exact
// IEEE-754 faults and record the same structured failure.  At most one failure
// is reported per observation; current energy takes precedence when both
// fields are invalid.
template <typename FailureHandler>
static inline bool QC_SCF_Require_Finite_Energy_Observation(
    int md_step, int iteration, bool has_previous_energy, double current_energy,
    double delta_energy, FailureHandler failure_handler)
{
    auto report = [&](QC_SCF_Energy_Quantity quantity, double value)
    {
        const QC_SCF_Energy_Validation_Failure failure = {
            quantity, QC_SCF_Energy_Quantity_Name(quantity), md_step, iteration,
            value,
        };
        failure_handler(failure);
    };

    if (!QC_SCF_Energy_Double_Is_Finite(current_energy))
    {
        report(QC_SCF_CURRENT_ENERGY, current_energy);
        return false;
    }
    if (has_previous_energy && !QC_SCF_Energy_Double_Is_Finite(delta_energy))
    {
        report(QC_SCF_DELTA_ENERGY, delta_energy);
        return false;
    }
    return true;
}
