#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <type_traits>
#include <utility>

#include "Domain_decomposition/Domain_decomposition.h"

using DomainCornerStorage =
    decltype(std::declval<DOMAIN_INFORMATION&>().min_corner_set);

static_assert(!std::is_array_v<DomainCornerStorage>,
              "domain corners must not have a compile-time rank limit");
static_assert(
    std::is_same_v<decltype(std::declval<DOMAIN_INFORMATION&>().max_corner_set),
                   DomainCornerStorage>,
    "minimum and maximum domain corners must use the same storage model");

int CONTROLLER::MPI_rank = 0;
int CONTROLLER::PP_MPI_size = 1;

bool CONTROLLER::Command_Exist(const char* prefix, const char* key)
{
    return commands.count(std::string(prefix) + "_" + key) != 0;
}

const char* CONTROLLER::Command(const char* prefix, const char* key)
{
    return commands[std::string(prefix) + "_" + key].c_str();
}

bool Float_Memory_Is_Finite(const void* address)
{
    return std::isfinite(*static_cast<const float*>(address));
}

bool Float_Memory_Is_Normal(const void* address)
{
    return std::isnormal(*static_cast<const float*>(address));
}

bool Float_Memory_Is_Zero_Or_Normal(const void* address)
{
    const float value = *static_cast<const float*>(address);
    return value == 0.0f || std::isnormal(value);
}

bool Double_Memory_Is_Finite(const void* address)
{
    return std::isfinite(*static_cast<const double*>(address));
}

void deviceMemset(void* destination, int value, std::size_t size)
{
    std::memset(destination, value, size);
}

namespace
{

bool Check_Rank_Count(int rank_count)
{
    CONTROLLER controller{};
    CONTROLLER::PP_MPI_size = rank_count;
    controller.commands["DOM_DEC_split_nx"] = std::to_string(rank_count);
    controller.commands["DOM_DEC_split_ny"] = "1";
    controller.commands["DOM_DEC_split_nz"] = "1";

    MD_INFORMATION md_info{};
    constexpr float local_extent = 30.0f;
    const float box_x = local_extent * static_cast<float>(rank_count);
    md_info.sys.box_length = {box_x, local_extent, local_extent};
    md_info.sys.box_angle = {90.0f, 90.0f, 90.0f};
    md_info.pbc.cell = {box_x, 0.0f, local_extent, 0.0f, 0.0f, local_extent};
    md_info.nb.cutoff = 10.0f;
    md_info.nb.skin = 2.0f;

    DOMAIN_INFORMATION domain{};
    domain.Domain_Decomposition(&controller, &md_info);

    if (domain.min_corner_set.size() != static_cast<std::size_t>(rank_count) ||
        domain.max_corner_set.size() != static_cast<std::size_t>(rank_count))
    {
        std::fprintf(stderr, "wrong corner capacity for %d ranks\n",
                     rank_count);
        return false;
    }

    for (int rank = 0; rank < rank_count; ++rank)
    {
        const float expected_min = local_extent * static_cast<float>(rank);
        const float expected_max = expected_min + local_extent;
        const VECTOR& minimum = domain.min_corner_set[rank];
        const VECTOR& maximum = domain.max_corner_set[rank];
        if (minimum.x != expected_min || minimum.y != 0.0f ||
            minimum.z != 0.0f || maximum.x != expected_max ||
            maximum.y != local_extent || maximum.z != local_extent)
        {
            std::fprintf(stderr, "domain corner mismatch at rank %d\n", rank);
            return false;
        }
    }

    std::free(domain.d_sum_ene_local);
    std::free(domain.d_sum_ene_total);
    std::free(domain.d_ek_local);
    std::free(domain.d_ek_total);
    return true;
}

}  // namespace

int main()
{
    // 17 is the first rank count that overflowed the former fixed arrays.
    // A larger non-power-of-two count also exercises normal vector growth.
    if (!Check_Rank_Count(17) || !Check_Rank_Count(257)) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
