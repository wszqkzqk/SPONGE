#pragma once

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace sponge_qc_input
{

inline bool Float_Is_Finite(float value) noexcept
{
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires a 32-bit float representation");
    std::memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return (bits & UINT32_C(0x7f800000)) != UINT32_C(0x7f800000);
}

inline int Parse_Exact_Int(const char* token)
{
    if (token == nullptr) throw std::invalid_argument("integer token is null");
    errno = 0;
    char* end = nullptr;
    const long long parsed = std::strtoll(token, &end, 10);
    if (end == token || end[0] != '\0')
        throw std::invalid_argument("value is not an exact integer token");
    if (errno == ERANGE || parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max())
        throw std::overflow_error("integer value is outside the int range");
    return static_cast<int>(parsed);
}

inline float Parse_Finite_Nonnegative_Float(const char* token)
{
    if (token == nullptr)
        throw std::invalid_argument("floating-point token is null");
    errno = 0;
    char* end = nullptr;
    const float parsed = std::strtof(token, &end);
    if (end == token || end[0] != '\0')
        throw std::invalid_argument(
            "value is not an exact floating-point token");
    if (errno == ERANGE)
        throw std::overflow_error(
            "floating-point value is outside the float range");
    if (!Float_Is_Finite(parsed) || parsed < 0.0f)
        throw std::domain_error(
            "floating-point value must be finite and nonnegative");
    return parsed;
}

}  // namespace sponge_qc_input
