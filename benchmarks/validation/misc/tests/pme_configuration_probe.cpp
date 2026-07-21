#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

// Include the implementation in this isolated translation unit so the probe
// exercises the exact private validation and layout helpers used by
// Particle_Mesh::Initial, Domain_Decomposition, and Get_Atoms.  Dead-section
// elimination discards the unrelated force kernels at link time.
#include "PM_force/PM_force.cpp"

using PMToPPStorage = decltype(std::declval<Particle_Mesh&>().pm_pp_corres);
using PMCornerStorage = decltype(std::declval<Particle_Mesh&>().min_corner_set);
using PMCurrentPPStorage =
    decltype(std::declval<Particle_Mesh&>().pm_corres_pp_rank_set);

static_assert(!std::is_array_v<PMToPPStorage>,
              "PM/PP correspondence must not have a fixed rank capacity");
static_assert(!std::is_array_v<PMCornerStorage>,
              "PM corners must not have a fixed rank capacity");
static_assert(!std::is_array_v<PMCurrentPPStorage>,
              "current-PM PP ranks must not have a fixed rank capacity");
static_assert(PME_MPI_TAG_DOMAIN_MIN_CORNER >= 0 && PME_MPI_TAG_FORCES <= 32767,
              "PME protocol tags must fit MPI's guaranteed tag range");

int CONTROLLER::MPI_rank = 0;
int CONTROLLER::MPI_size = 1;
int CONTROLLER::PP_MPI_size = 1;
int CONTROLLER::PM_MPI_size = 1;
int CONTROLLER::CC_MPI_size = 0;

namespace
{

bool Fail(const char* message)
{
    std::fprintf(stderr, "%s\n", message);
    return false;
}

int Run_Device_Atom_ID_Validation(const std::vector<int>& local_to_global,
                                  std::vector<int>* candidate_inverse)
{
    candidate_inverse->assign(local_to_global.size(), -1);
    int error[3] = {PME_ATOM_IDS_VALID, -1, -1};
    Build_PME_Inverse_Atom_ID_Candidate(
        local_to_global.data(), candidate_inverse->data(), error,
        static_cast<int>(local_to_global.size()));
    Check_PME_Inverse_Atom_ID_Candidate(
        candidate_inverse->data(), error,
        static_cast<int>(local_to_global.size()));
    return error[0];
}

bool Check_Single_PM_Domain_Layout_Capacity(int pp_ranks)
{
    CONTROLLER controller{};
    CONTROLLER::MPI_rank = 0;
    CONTROLLER::MPI_size = pp_ranks + 1;
    CONTROLLER::PP_MPI_size = pp_ranks;
    CONTROLLER::PM_MPI_size = 1;
    CONTROLLER::CC_MPI_size = 0;

    Particle_Mesh pme{};
    pme.PM_MPI_size = 1;
    const float box_x = 10.0f * static_cast<float>(pp_ranks);
    pme.Domain_Decomposition(&controller, {box_x, 10.0f, 10.0f},
                             {pp_ranks, 1, 1});

    if (pme.pm_pp_corres.size() != 1 || pme.pm_pp_num.size() != 1 ||
        pme.min_corner_set.size() != 1 || pme.max_corner_set.size() != 1)
        return Fail("single-PM dynamic outer storage has the wrong size");
    if (pme.pm_pp_num[0] != pp_ranks ||
        pme.pm_pp_corres[0].size() != static_cast<std::size_t>(pp_ranks))
        return Fail("single PM did not retain every PP rank");
    for (int rank = 0; rank < pp_ranks; ++rank)
    {
        if (pme.pm_pp_corres[0][rank] != rank)
            return Fail("single-PM correspondence is incomplete or reordered");
    }
    if (pme.pm_dom_dec_split_num.int_x != 1 ||
        pme.pm_dom_dec_split_num.int_y != 1 ||
        pme.pm_dom_dec_split_num.int_z != 1)
        return Fail("single-PM split is not (1, 1, 1)");
    if (pme.min_corner_set[0].x != 0.0f || pme.max_corner_set[0].x != box_x ||
        pme.max_corner_set[0].y != 10.0f || pme.max_corner_set[0].z != 10.0f)
        return Fail("single-PM corner bounds are wrong");
    return true;
}

bool Check_Dormant_Multi_PM_Helper_Layout()
{
    // Multi-PM execution remains deliberately rejected by Initial.  This
    // helper-only case verifies that the dormant layout data structure itself
    // no longer truncates either 101 PM ranks or their 202 PP assignments.
    PME_Domain_Layout layout;
    std::string error;
    if (!Try_Build_PME_Domain_Layout(101, 202, {202, 1, 1},
                                     {2020.0f, 10.0f, 10.0f}, &layout, &error))
    {
        std::fprintf(stderr, "101-PM helper layout failed: %s\n",
                     error.c_str());
        return false;
    }
    if (layout.pm_to_pp.size() != 101 || layout.minimum_corners.size() != 101 ||
        layout.maximum_corners.size() != 101)
        return Fail("101-PM helper storage has the wrong size");

    std::vector<unsigned char> seen(202, 0);
    for (const std::vector<int>& ranks : layout.pm_to_pp)
    {
        if (ranks.size() != 2)
            return Fail("101-PM helper layout is not uniform");
        for (int rank : ranks)
        {
            if (rank < 0 || rank >= 202 || seen[rank])
                return Fail("101-PM helper layout has an invalid PP rank");
            seen[rank] = 1;
        }
    }
    for (unsigned char was_seen : seen)
        if (!was_seen) return Fail("101-PM helper layout lost a PP rank");
    return true;
}

bool Check_Strict_Configuration_Validation()
{
    std::string error;
    float parsed_float = 0.0f;
    int parsed_integer = 0;
    const char* invalid_positive_floats[] = {"",      "0",     "-1",   "nan",
                                             "NaN",   "inf",   "+inf", "1tail",
                                             "1e999", "1e-50", " 1",   "0x1p0"};
    for (const char* token : invalid_positive_floats)
    {
        if (Try_Parse_PME_Positive_Normal_Float(token, &parsed_float, &error))
        {
            std::fprintf(stderr, "invalid positive float was accepted: %s\n",
                         token);
            return false;
        }
    }
    if (!Try_Parse_PME_Positive_Normal_Float("1e-5", &parsed_float, &error) ||
        parsed_float != 1.0e-5f)
        return Fail("valid Direct_Tolerance token was rejected");
    if (!Try_Parse_PME_Positive_Normal_Float("1.25", &parsed_float, &error) ||
        parsed_float != 1.25f)
        return Fail("valid grid_spacing token was rejected");

    const char* invalid_integers[] = {"", "4x", "1.0", " 4",
                                      "999999999999999999999"};
    for (const char* token : invalid_integers)
    {
        if (Try_Parse_PME_Integer(token, &parsed_integer, &error))
        {
            std::fprintf(stderr, "invalid integer was accepted: %s\n", token);
            return false;
        }
    }

    float beta = 0.0f;
    if (!Try_Compute_PME_Beta(10.0, 1.0e-5, &beta, &error) ||
        !Float_Memory_Is_Normal(&beta) || !(beta > 0.0f))
        return Fail("valid PME beta could not be solved");
    const double invalid_tolerances[] = {
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
        0.1,
    };
    for (double tolerance : invalid_tolerances)
    {
        if (Try_Compute_PME_Beta(10.0, tolerance, &beta, &error))
            return Fail("invalid Direct_Tolerance produced a beta");
    }
    if (Try_Compute_PME_Beta(0.0, 1.0e-5, &beta, &error) ||
        Try_Compute_PME_Beta(-1.0, 1.0e-5, &beta, &error) ||
        Try_Compute_PME_Beta(std::numeric_limits<double>::infinity(), 1.0e-5,
                             &beta, &error))
        return Fail("invalid cutoff produced a beta");
    if (Try_Compute_PME_Beta(static_cast<double>(FLT_MAX),
                             0.5 / static_cast<double>(FLT_MAX), &beta, &error))
        return Fail("subnormal float beta was accepted under FTZ rules");

    int automatic_dimension = 0;
    if (!Try_Get_PME_Fft_Parameter(30.0, &automatic_dimension, &error) ||
        automatic_dimension != 32)
        return Fail("valid automatic FFT dimension was not reproduced");
    const double invalid_mesh_lengths[] = {
        0.0,
        -1.0,
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
        static_cast<double>(INT_MAX),
    };
    for (double mesh_length : invalid_mesh_lengths)
    {
        if (Try_Get_PME_Fft_Parameter(mesh_length, &automatic_dimension,
                                      &error))
            return Fail("invalid automatic FFT ratio was accepted");
    }

    PME_Grid_Shape shape;
    if (!Try_Build_PME_Grid_Shape(64, 64, 64, &shape, &error) ||
        shape.all != 64 * 64 * 64)
        return Fail("valid FFT grid shape was rejected");
    if (Try_Build_PME_Grid_Shape(0, 64, 64, &shape, &error) ||
        Try_Build_PME_Grid_Shape(-4, 64, 64, &shape, &error) ||
        Try_Build_PME_Grid_Shape(3, 64, 64, &shape, &error) ||
        Try_Build_PME_Grid_Shape(50000, 50000, 4, &shape, &error) ||
        Try_Build_PME_Grid_Shape(INT_MAX, 4, 4, &shape, &error))
        return Fail("invalid or overflowing FFT grid shape was accepted");

    if (!Try_Validate_PME_Process_Count(258, 0, 1, &error) ||
        !Try_Validate_PME_Process_Count(1, 0, 0, &error) ||
        Try_Validate_PME_Process_Count(0, 0, 0, &error) ||
        Try_Validate_PME_Process_Count(8, 0, -1, &error) ||
        Try_Validate_PME_Process_Count(8, 0, 8, &error) ||
        Try_Validate_PME_Process_Count(8, 8, 0, &error))
        return Fail("PM rank-count validation contract is wrong");
    if (!Try_Validate_PME_Implemented_Process_Roles(0, &error) ||
        Try_Validate_PME_Implemented_Process_Roles(1, &error))
        return Fail("unsupported CC-rank layouts are not rejected explicitly");

    PME_Domain_Layout invalid_layout;
    if (Try_Build_PME_Domain_Layout(2, 3, {3, 1, 1}, {30.0f, 10.0f, 10.0f},
                                    &invalid_layout, &error) ||
        Try_Build_PME_Domain_Layout(1, 3, {0, 3, 1}, {30.0f, 10.0f, 10.0f},
                                    &invalid_layout, &error) ||
        Try_Build_PME_Domain_Layout(1, 3, {2, 1, 1}, {30.0f, 10.0f, 10.0f},
                                    &invalid_layout, &error))
        return Fail("invalid PM/PP domain layout was accepted");
    return true;
}

bool Check_Atom_Transfer_Validation()
{
    std::string error;

    std::vector<int> prefixes = {91};
    if (!Try_Build_PME_Atom_Count_Prefix({2, 0, 3}, 5, &prefixes, &error) ||
        prefixes != std::vector<int>({0, 2, 2}))
        return Fail("a complete PP atom-count partition was rejected");

    const std::vector<std::vector<int>> invalid_count_partitions = {
        {2, 0, 2},  // short total
        {2, -1, 4},
        {6},
        {INT_MAX},  // MPI byte count does not fit for atom records
    };
    const int invalid_global_counts[] = {5, 5, 5, INT_MAX};
    for (std::size_t i = 0; i < invalid_count_partitions.size(); ++i)
    {
        std::vector<int> unchanged = {73, 74};
        if (Try_Build_PME_Atom_Count_Prefix(invalid_count_partitions[i],
                                            invalid_global_counts[i],
                                            &unchanged, &error))
            return Fail("an invalid PP atom-count partition was accepted");
        if (unchanged != std::vector<int>({73, 74}))
            return Fail("failed count validation committed a partial prefix");
    }

    std::vector<int> inverse;
    if (Run_Device_Atom_ID_Validation({2, 0, 3, 1}, &inverse) !=
            PME_ATOM_IDS_VALID ||
        inverse != std::vector<int>({1, 3, 0, 2}))
        return Fail("a complete atom-ID permutation was rejected");

    const std::vector<std::vector<int>> invalid_id_sets = {
        {0, 1, 1, 3},
        {0, 1, -1, 3},
        {0, 1, 4, 3},
    };
    const int expected_errors[] = {PME_ATOM_ID_DUPLICATE,
                                   PME_ATOM_ID_OUT_OF_RANGE,
                                   PME_ATOM_ID_OUT_OF_RANGE};
    for (std::size_t i = 0; i < invalid_id_sets.size(); ++i)
    {
        std::vector<int> candidate;
        std::vector<int> committed_inverse = {81, 82};
        const int validation =
            Run_Device_Atom_ID_Validation(invalid_id_sets[i], &candidate);
        if (validation != expected_errors[i])
            return Fail(
                "a duplicate or out-of-range atom-ID mapping was accepted");
        if (validation == PME_ATOM_IDS_VALID) committed_inverse = candidate;
        if (committed_inverse != std::vector<int>({81, 82}))
            return Fail("failed device ID validation committed its candidate");
    }

    const int incomplete_candidate[] = {0, -1, 2, 3};
    int missing_error[3] = {PME_ATOM_IDS_VALID, -1, -1};
    Check_PME_Inverse_Atom_ID_Candidate(incomplete_candidate, missing_error, 4);
    if (missing_error[0] != PME_ATOM_ID_MISSING || missing_error[1] != 1)
        return Fail("an incomplete inverse atom-ID mapping was not detected");

    const int protocol_tags[] = {
        PME_MPI_TAG_DOMAIN_MIN_CORNER,
        PME_MPI_TAG_DOMAIN_MAX_CORNER,
        PME_MPI_TAG_DOMAIN_PP_COUNT,
        PME_MPI_TAG_DOMAIN_PP_RANKS,
        PME_MPI_TAG_DOMAIN_SPLIT,
        PME_MPI_TAG_DOMAIN_PM_ASSIGNMENT,
        PME_MPI_TAG_ATOM_COUNT,
        PME_MPI_TAG_COORDINATES,
        PME_MPI_TAG_CHARGES,
        PME_MPI_TAG_ATOM_IDS,
        PME_MPI_TAG_FORCES,
    };
    for (std::size_t i = 0; i < std::size(protocol_tags); ++i)
    {
        if (protocol_tags[i] < 0 || protocol_tags[i] > 32767)
            return Fail("a PME protocol tag exceeds MPI's guaranteed range");
        for (std::size_t j = i + 1; j < std::size(protocol_tags); ++j)
            if (protocol_tags[i] == protocol_tags[j])
                return Fail("PME protocol tags are not unique");
    }
    return true;
}

}  // namespace

int main()
{
    // This isolated, non-MPI probe covers the production domain-layout helper,
    // not an end-to-end MPI exchange.  The two capacities correspond to
    // production-reachable totals of 102 and 258 ranks and both exceed the
    // former fixed PP storage limit of 100.
    if (!Check_Single_PM_Domain_Layout_Capacity(101) ||
        !Check_Single_PM_Domain_Layout_Capacity(257) ||
        !Check_Dormant_Multi_PM_Helper_Layout() ||
        !Check_Strict_Configuration_Validation() ||
        !Check_Atom_Transfer_Validation())
        return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
