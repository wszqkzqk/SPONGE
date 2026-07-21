#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "common.h"

namespace
{

struct Case
{
    const char* name;
    unsigned int bits;
    bool finite;
    bool normal;
    bool zero_or_normal;
};

float Float_From_Bits(unsigned int bits)
{
    float value = 0.0f;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

}  // namespace

int main()
{
    const Case cases[] = {
        {"positive zero", 0x00000000U, true, false, true},
        {"negative zero", 0x80000000U, true, false, true},
        {"positive subnormal", 0x00000001U, true, false, false},
        {"negative subnormal", 0x80000001U, true, false, false},
        {"minimum normal", 0x00800000U, true, true, true},
        {"maximum normal", 0x7f7fffffU, true, true, true},
        {"positive infinity", 0x7f800000U, false, false, false},
        {"negative infinity", 0xff800000U, false, false, false},
        {"quiet NaN", 0x7fc00000U, false, false, false},
    };
    for (const Case& test : cases)
    {
        const float value = Float_From_Bits(test.bits);
        const bool finite = Float_Memory_Is_Finite(&value);
        const bool normal = Float_Memory_Is_Normal(&value);
        const bool zero_or_normal = Float_Memory_Is_Zero_Or_Normal(&value);
        if (finite != test.finite || normal != test.normal ||
            zero_or_normal != test.zero_or_normal)
        {
            std::fprintf(
                stderr,
                "%s (0x%08x): finite/normal/zero-or-normal=%d/%d/%d, "
                "expected %d/%d/%d\n",
                test.name, test.bits, static_cast<int>(finite),
                static_cast<int>(normal), static_cast<int>(zero_or_normal),
                static_cast<int>(test.finite), static_cast<int>(test.normal),
                static_cast<int>(test.zero_or_normal));
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
