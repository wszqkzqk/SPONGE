#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "quantum_chemistry/scf/energy_validation_policy.hpp"

namespace
{

double Double_From_Bits(std::uint64_t bit_pattern)
{
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bit_pattern));
#endif
    double value = 0.0;
    std::memcpy(&value, &bit_pattern, sizeof(value));
    return value;
}

std::uint64_t Double_To_Bits(double value)
{
    std::uint64_t bit_pattern = 0;
    std::memcpy(&bit_pattern, &value, sizeof(bit_pattern));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bit_pattern));
#endif
    return bit_pattern;
}

struct Failure_Recorder
{
    int calls = 0;
    QC_SCF_Energy_Validation_Failure last = {
        QC_SCF_CURRENT_ENERGY, nullptr, -1, -1, 0.0,
    };

    void operator()(const QC_SCF_Energy_Validation_Failure& failure)
    {
        ++calls;
        last = failure;
    }
};

bool Validate_Fault(std::uint64_t fault_bits, QC_SCF_Energy_Quantity quantity)
{
    const double finite_energy = Double_From_Bits(UINT64_C(0x3ff4000000000000));
    const double finite_delta = Double_From_Bits(UINT64_C(0xbeb0c6f7a0b5ed8d));
    const double fault = Double_From_Bits(fault_bits);
    Failure_Recorder recorder;
    const bool accepted = QC_SCF_Require_Finite_Energy_Observation(
        17, 23, true, quantity == QC_SCF_CURRENT_ENERGY ? fault : finite_energy,
        quantity == QC_SCF_DELTA_ENERGY ? fault : finite_delta,
        [&](const QC_SCF_Energy_Validation_Failure& failure)
        { recorder(failure); });

    return !accepted && recorder.calls == 1 &&
           recorder.last.quantity == quantity &&
           recorder.last.quantity_name != nullptr &&
           std::strcmp(recorder.last.quantity_name,
                       QC_SCF_Energy_Quantity_Name(quantity)) == 0 &&
           recorder.last.md_step == 17 && recorder.last.iteration == 23 &&
           Double_To_Bits(recorder.last.value) == fault_bits;
}

}  // namespace

int main()
{
    const double quiet_nan = Double_From_Bits(UINT64_C(0x7ff8000000000042));
    const double finite_max = Double_From_Bits(UINT64_C(0x7fefffffffffffff));
    const double finite_negative_max =
        Double_From_Bits(UINT64_C(0xffefffffffffffff));

    Failure_Recorder recorder;
    if (!QC_SCF_Require_Finite_Energy_Observation(
            0, 1, true, finite_max, finite_negative_max,
            [&](const QC_SCF_Energy_Validation_Failure& failure)
            { recorder(failure); }) ||
        recorder.calls != 0)
    {
        std::fprintf(stderr, "finite energy observation was rejected\n");
        return EXIT_FAILURE;
    }

    // The first SCF iteration has no previous energy.  Its delta storage may
    // therefore contain any sentinel without being consumed by the policy.
    recorder = Failure_Recorder{};
    if (!QC_SCF_Require_Finite_Energy_Observation(
            3, 1, false, 0.0, quiet_nan,
            [&](const QC_SCF_Energy_Validation_Failure& failure)
            { recorder(failure); }) ||
        recorder.calls != 0)
    {
        std::fprintf(stderr, "first-iteration delta sentinel was consumed\n");
        return EXIT_FAILURE;
    }

    const std::uint64_t faults[] = {
        UINT64_C(0x7ff8000000000042),  // quiet NaN with a non-zero payload
        UINT64_C(0x7ff0000000000000),  // positive infinity
        UINT64_C(0xfff0000000000000),  // negative infinity
    };
    for (std::uint64_t fault : faults)
    {
        if (!Validate_Fault(fault, QC_SCF_CURRENT_ENERGY) ||
            !Validate_Fault(fault, QC_SCF_DELTA_ENERGY))
        {
            std::fprintf(stderr,
                         "non-finite energy fault was not classified exactly "
                         "once\n");
            return EXIT_FAILURE;
        }
    }

    // If both fields are invalid, current energy is the primary failure and
    // the handler must still run exactly once.
    recorder = Failure_Recorder{};
    if (QC_SCF_Require_Finite_Energy_Observation(
            5, 9, true, quiet_nan,
            Double_From_Bits(UINT64_C(0x7ff0000000000000)),
            [&](const QC_SCF_Energy_Validation_Failure& failure)
            { recorder(failure); }) ||
        recorder.calls != 1 || recorder.last.quantity != QC_SCF_CURRENT_ENERGY)
    {
        std::fprintf(stderr,
                     "multiple energy faults were reported ambiguously\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
