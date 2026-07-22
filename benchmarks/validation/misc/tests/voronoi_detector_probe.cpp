#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

// Exercise the production parser and state machine directly.  The engine
// integration entry points are discarded by section-garbage collection.
#include "bias/voronoi_detector.cpp"

namespace
{

bool Fail(const std::string& message)
{
    std::fprintf(stderr, "%s\n", message.c_str());
    return false;
}

bool Check(bool condition, const std::string& message)
{
    return condition ? true : Fail(message);
}

bool Write_File(const std::string& path, const std::string& contents)
{
    std::ofstream output(path);
    output << contents;
    return output.good();
}

const char* const kBranchingManifest =
    "4\n"
    "M_0 0 0\n"
    "M_1 2 0\n"
    "M_2 0 2\n"
    "M_3 2 2\n"
    "3\n"
    "S_0_1 0 1\n"
    "S_0_2 0 2\n"
    "S_1_3 1 3\n";

bool Check_Source_Recrossing_State_Machine(const std::string& manifest)
{
    VORONOI_DETECTOR detector;
    std::string error;
    if (!detector.Load_Milestone_File(manifest, "S_0_1", 2, &error))
        return Fail("valid branching manifest was rejected: " + error);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (!detector.Observe_CV_Values({nan, nan}, -1, false, &error))
        return Fail("uncommitted initial observation was validated: " + error);
    if (!Check(!detector.initial_state_observed &&
                   detector.current_milestone == -1 &&
                   detector.source_recrossing_count == 0,
               "uncommitted initial observation changed detector state"))
        return false;

    if (!detector.Observe_CV_Values({0.9, 0.0}, 7, true, &error))
        return Fail("valid source-side initial state was rejected: " + error);
    if (!Check(detector.initial_state_observed &&
                   detector.current_milestone == 0 &&
                   detector.initial_step == 7 && !detector.terminal_hit,
               "initial source-side classification is wrong"))
        return false;

    // Exact contact with the source interface is not strict entry into its
    // other cell and must not depend on anchor ordering.
    if (!detector.Observe_CV_Values({1.0, 0.0}, 8, true, &error))
        return Fail("current-cell source tie was rejected: " + error);
    if (!Check(detector.current_milestone == 0 &&
                   detector.source_recrossing_count == 0,
               "source tie was misclassified as a crossing"))
        return false;

    if (!detector.Observe_CV_Values({0.9, 1.2}, 9, false, &error))
        return Fail("uncommitted destination trial was rejected: " + error);
    if (!Check(detector.current_milestone == 0 && !detector.terminal_hit &&
                   detector.initial_step == 7,
               "uncommitted destination trial changed detector state"))
        return false;

    if (!detector.Observe_CV_Values({1.1, 0.0}, 10, true, &error))
        return Fail("source recrossing 0 -> 1 was rejected: " + error);
    if (!Check(detector.current_milestone == 1 && !detector.terminal_hit &&
                   detector.source_recrossing_count == 1 &&
                   detector.initial_step == 7,
               "source recrossing 0 -> 1 became terminal or reset time"))
        return false;

    if (!detector.Observe_CV_Values({0.9, 0.0}, 12, true, &error))
        return Fail("source recrossing 1 -> 0 was rejected: " + error);
    if (!Check(detector.current_milestone == 0 && !detector.terminal_hit &&
                   detector.source_recrossing_count == 2 &&
                   detector.initial_step == 7,
               "source recrossing 1 -> 0 became terminal or reset time"))
        return false;

    if (!detector.Observe_CV_Values({0.9, 1.2}, 15, true, &error))
        return Fail("valid non-source first hit was rejected: " + error);
    if (!Check(detector.terminal_hit && detector.hit_step == 15 &&
                   detector.initial_step == 7 &&
                   detector.hit_from_milestone == 0 &&
                   detector.destination_milestone == 2 &&
                   detector.Hit_Interface().name == "S_0_2" &&
                   detector.Hit_Restart_Basename() == "voronoi_hit_S_0_2",
               "first non-source interface was not the terminal hit"))
        return false;

    const int stable_hit_step = detector.hit_step;
    if (!detector.Observe_CV_Values({nan, nan}, 16, true, &error))
        return Fail("terminal detector did not remain terminal: " + error);
    return Check(detector.hit_step == stable_hit_step &&
                     detector.destination_milestone == 2,
                 "observation after terminal hit changed the result");
}

bool Check_Either_Source_Side_Is_Valid(const std::string& manifest)
{
    VORONOI_DETECTOR detector;
    std::string error;
    if (!detector.Load_Milestone_File(manifest, "S_0_1", 2, &error) ||
        !detector.Observe_CV_Values({1.1, 0.0}, 0, true, &error))
        return Fail("second source side was rejected: " + error);
    if (!Check(detector.current_milestone == 1,
               "initial source side was not inferred from coordinates"))
        return false;
    if (!detector.Observe_CV_Values({1.1, 1.2}, 3, true, &error))
        return Fail("target from second source side was rejected: " + error);
    return Check(detector.terminal_hit &&
                     detector.Hit_Interface().name == "S_1_3" &&
                     detector.hit_from_milestone == 1 &&
                     detector.destination_milestone == 3,
                 "branching graph was treated as a hard-coded linear graph");
}

bool Check_Exact_Other_Interface_Contact(const std::string& directory,
                                         const std::string& manifest)
{
    VORONOI_DETECTOR detector;
    std::string error;
    if (!detector.Load_Milestone_File(manifest, "S_0_1", 2, &error) ||
        !detector.Observe_CV_Values({0.1, 0.0}, 0, true, &error) ||
        !detector.Observe_CV_Values({0.0, 1.0}, 1, true, &error))
        return Fail("exact non-source interface contact was rejected: " +
                    error);
    if (!Check(detector.terminal_hit && detector.hit_step == 1 &&
                   detector.hit_from_milestone == 0 &&
                   detector.destination_milestone == 2 &&
                   detector.Hit_Interface().name == "S_0_2",
               "exact non-source interface contact was not terminal"))
        return false;

    const std::string ambiguous_manifest = directory + "/ambiguous_contact.txt";
    if (!Write_File(ambiguous_manifest,
                    "4\n"
                    "M_0 0 0\n"
                    "M_1 -2 0\n"
                    "M_2 2 0\n"
                    "M_3 0 2\n"
                    "3\n"
                    "S_0_1 0 1\n"
                    "S_0_2 0 2\n"
                    "S_0_3 0 3\n"))
        return Fail("could not write ambiguous-contact manifest");

    VORONOI_DETECTOR ambiguous;
    if (!ambiguous.Load_Milestone_File(ambiguous_manifest, "S_0_1", 2,
                                       &error) ||
        !ambiguous.Observe_CV_Values({0.1, 0.0}, 0, true, &error))
        return Fail("could not initialize ambiguous-contact case: " + error);
    if (ambiguous.Observe_CV_Values({1.0, 1.0}, 1, true, &error))
        return Fail("simultaneous contact with two destinations was accepted");
    return Check(
        error.find("multiple non-source interfaces") != std::string::npos &&
            !ambiguous.terminal_hit && ambiguous.current_milestone == 0,
        "ambiguous interface contact did not fail closed");
}

bool Check_Initial_And_Runtime_Ambiguity(const std::string& directory)
{
    const std::string manifest = directory + "/ambiguous.txt";
    if (!Write_File(manifest,
                    "4\n"
                    "M_0 0 0\n"
                    "M_1 -2 0\n"
                    "M_2 2 1\n"
                    "M_3 2 -1\n"
                    "3\n"
                    "S_0_1 0 1\n"
                    "S_0_2 0 2\n"
                    "S_0_3 0 3\n"))
        return Fail("could not write ambiguity manifest");

    std::string error;
    VORONOI_DETECTOR initial_tie;
    if (!initial_tie.Load_Milestone_File(manifest, "S_0_1", 2, &error))
        return Fail("ambiguity manifest was rejected: " + error);
    if (initial_tie.Observe_CV_Values({-1.0, 0.0}, 0, true, &error))
        return Fail("initial source-interface tie was accepted");
    if (!Check(error.find("initial state is tied") != std::string::npos &&
                   !initial_tie.initial_state_observed,
               "initial tie did not fail without publishing state"))
        return false;

    VORONOI_DETECTOR outside_source;
    if (!outside_source.Load_Milestone_File(manifest, "S_0_1", 2, &error))
        return Fail("ambiguity manifest reload failed: " + error);
    if (outside_source.Observe_CV_Values({2.0, 1.0}, 0, true, &error))
        return Fail("initial state outside both source cells was accepted");
    if (!Check(error.find("not incident") != std::string::npos &&
                   !outside_source.initial_state_observed,
               "outside-source initial state published detector state"))
        return false;

    VORONOI_DETECTOR destination_tie;
    if (!destination_tie.Load_Milestone_File(manifest, "S_0_1", 2, &error) ||
        !destination_tie.Observe_CV_Values({0.1, 0.0}, 0, true, &error))
        return Fail("could not initialize destination-tie case: " + error);
    if (destination_tie.Observe_CV_Values({2.0, 0.0}, 1, true, &error))
        return Fail("ambiguous destination tie was accepted");
    return Check(error.find("ambiguous destination") != std::string::npos &&
                     destination_tie.current_milestone == 0 &&
                     !destination_tie.terminal_hit,
                 "ambiguous destination changed current or terminal state");
}

bool Check_Nonadjacent_Jump_Fails_Closed(const std::string& directory)
{
    const std::string manifest = directory + "/nonadjacent.txt";
    if (!Write_File(manifest,
                    "3\n"
                    "M_0 0\n"
                    "M_1 2\n"
                    "M_2 4\n"
                    "2\n"
                    "S_0_1 0 1\n"
                    "S_1_2 1 2\n"))
        return Fail("could not write nonadjacent manifest");
    VORONOI_DETECTOR detector;
    std::string error;
    if (!detector.Load_Milestone_File(manifest, "S_0_1", 1, &error) ||
        !detector.Observe_CV_Values({0.1}, 0, true, &error))
        return Fail("could not initialize nonadjacent case: " + error);
    if (detector.Observe_CV_Values({3.9}, 1, true, &error))
        return Fail("nonadjacent jump was accepted");
    return Check(error.find("nonadjacent Voronoi jump") != std::string::npos &&
                     detector.current_milestone == 0 && !detector.terminal_hit,
                 "nonadjacent jump did not fail closed");
}

bool Check_Strict_Parser_And_Finite_Runtime(const std::string& directory,
                                            const std::string& valid_manifest)
{
    const std::vector<std::pair<std::string, std::string>> invalid_cases = {
        {"truncated.txt", "2\nM_0 0 0\nM_1 2\n"},
        {"trailing.txt", std::string(kBranchingManifest) + "extra\n"},
        {"duplicate_edge.txt", "2\nM_0 0\nM_1 2\n2\nS_A 0 1\nS_B 1 0\n"},
        {"duplicate_anchor.txt", "2\nM_0 0\nM_1 0\n1\nS_0_1 0 1\n"},
        {"invalid_endpoint.txt", "2\nM_0 0\nM_1 2\n1\nS_0_1 0 2\n"},
        {"unsafe_name.txt", "2\nM_0 0\nM_1 2\n1\n../S 0 1\n"},
    };
    for (const auto& invalid : invalid_cases)
    {
        const std::string path = directory + "/" + invalid.first;
        if (!Write_File(path, invalid.second))
            return Fail("could not write " + invalid.first);
        VORONOI_DETECTOR detector;
        std::string error;
        if (detector.Load_Milestone_File(path, "S_0_1", 1, &error))
            return Fail("invalid manifest was accepted: " + invalid.first);
        if (!Check(detector.milestone_count == 0 &&
                       detector.source_interface == -1,
                   "failed parse published a partial graph: " + invalid.first))
            return false;
    }

    VORONOI_DETECTOR missing_source;
    std::string error;
    if (missing_source.Load_Milestone_File(valid_manifest, "S_missing", 2,
                                           &error))
        return Fail("missing explicit source_interface was accepted");
    if (!Check(error.find("not present") != std::string::npos,
               "missing source_interface produced the wrong error"))
        return false;

    VORONOI_DETECTOR detector;
    if (!detector.Load_Milestone_File(valid_manifest, "S_0_1", 2, &error))
        return Fail("valid manifest reload failed: " + error);
    const int retained_count = detector.milestone_count;
    const std::string bad_path = directory + "/trailing.txt";
    if (detector.Load_Milestone_File(bad_path, "S_0_1", 2, &error))
        return Fail("transactional reload accepted an invalid manifest");
    if (!Check(
            detector.milestone_count == retained_count &&
                detector.interfaces[detector.source_interface].name == "S_0_1",
            "failed reload destroyed the previously valid graph"))
        return false;

    const double infinity = std::numeric_limits<double>::infinity();
    if (detector.Observe_CV_Values({infinity, 0.0}, 0, true, &error))
        return Fail("non-finite runtime CV was accepted");
    return Check(error.find("non-finite") != std::string::npos &&
                     !detector.initial_state_observed,
                 "non-finite runtime CV published initial state");
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 2)
        return Fail("usage: voronoi_detector_probe DIRECTORY") ? 0 : 2;
    const std::string directory = argv[1];
    const std::string valid_manifest = directory + "/branching.txt";
    if (!Write_File(valid_manifest, kBranchingManifest)) return 2;

    if (!Check_Source_Recrossing_State_Machine(valid_manifest) ||
        !Check_Either_Source_Side_Is_Valid(valid_manifest) ||
        !Check_Exact_Other_Interface_Contact(directory, valid_manifest) ||
        !Check_Initial_And_Runtime_Ambiguity(directory) ||
        !Check_Nonadjacent_Jump_Fails_Closed(directory) ||
        !Check_Strict_Parser_And_Finite_Runtime(directory, valid_manifest))
        return 1;
    return 0;
}
