#include "voronoi_detector.h"

#include <cstring>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <utility>

namespace
{

bool Set_Voronoi_Error(std::string* error, const std::string& message)
{
    if (error != NULL) *error = message;
    return false;
}

bool Is_Safe_Interface_Name(const std::string& name)
{
    if (name.empty()) return false;
    for (char c : name)
    {
        if (!((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
              (c >= '0' && c <= '9') || c == '_' || c == '-'))
            return false;
    }
    return true;
}

std::pair<int, int> Canonical_Edge(int first, int second)
{
    if (first > second) std::swap(first, second);
    return {first, second};
}

#if defined(_MSC_VER)
__declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
bool Double_Bits_Are_Finite(const void* address)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(double), "unexpected double size");
    static_assert(std::numeric_limits<double>::is_iec559,
                  "IEEE-754 double precision is required");
    std::memcpy(&bits, address, sizeof(bits));
    return (bits & UINT64_C(0x7ff0000000000000)) !=
           UINT64_C(0x7ff0000000000000);
}

bool Is_Finite_Double(const double value)
{
    return Double_Bits_Are_Finite(&value);
}

}  // namespace

bool VORONOI_DETECTOR::Load_Milestone_File(
    const std::string& filename, const std::string& source_interface_name,
    int cv_numbers, std::string* error)
{
    if (error != NULL) error->clear();
    if (cv_numbers < 1)
        return Set_Voronoi_Error(error, "at least one CV is required");
    if (!Is_Safe_Interface_Name(source_interface_name))
        return Set_Voronoi_Error(
            error,
            "source_interface must be a nonempty filename-safe identifier");

    std::ifstream input(filename);
    if (!input.is_open())
        return Set_Voronoi_Error(
            error, "cannot open milestone file '" + filename + "'");

    int parsed_milestone_count = 0;
    if (!(input >> parsed_milestone_count))
        return Set_Voronoi_Error(error, "missing milestone count");
    if (parsed_milestone_count < 2)
        return Set_Voronoi_Error(error, "milestone count must be at least two");
    if (static_cast<std::size_t>(parsed_milestone_count) >
        std::numeric_limits<std::size_t>::max() /
            static_cast<std::size_t>(cv_numbers))
        return Set_Voronoi_Error(error, "milestone table size overflows");

    std::vector<std::string> parsed_names(parsed_milestone_count);
    std::vector<double> parsed_cvs(
        static_cast<std::size_t>(parsed_milestone_count) * cv_numbers);
    std::set<std::string> seen_milestone_names;
    for (int milestone = 0; milestone < parsed_milestone_count; ++milestone)
    {
        if (!(input >> parsed_names[milestone]))
        {
            return Set_Voronoi_Error(error,
                                     "truncated milestone record at index " +
                                         std::to_string(milestone));
        }
        if (!seen_milestone_names.insert(parsed_names[milestone]).second)
        {
            return Set_Voronoi_Error(error, "duplicate milestone name '" +
                                                parsed_names[milestone] + "'");
        }
        for (int cv = 0; cv < cv_numbers; ++cv)
        {
            double value = 0.0;
            if (!(input >> value))
            {
                return Set_Voronoi_Error(
                    error, "truncated or invalid CV value for milestone '" +
                               parsed_names[milestone] + "'");
            }
            if (!Is_Finite_Double(value))
            {
                return Set_Voronoi_Error(error,
                                         "non-finite CV value for milestone '" +
                                             parsed_names[milestone] + "'");
            }
            parsed_cvs[static_cast<std::size_t>(milestone) * cv_numbers + cv] =
                value;
        }
    }

    for (int first = 0; first < parsed_milestone_count; ++first)
    {
        for (int second = first + 1; second < parsed_milestone_count; ++second)
        {
            bool identical = true;
            for (int cv = 0; cv < cv_numbers; ++cv)
            {
                if (parsed_cvs[static_cast<std::size_t>(first) * cv_numbers +
                               cv] !=
                    parsed_cvs[static_cast<std::size_t>(second) * cv_numbers +
                               cv])
                {
                    identical = false;
                    break;
                }
            }
            if (identical)
            {
                return Set_Voronoi_Error(
                    error, "milestones '" + parsed_names[first] + "' and '" +
                               parsed_names[second] +
                               "' have identical anchor coordinates");
            }
        }
    }

    int interface_count = 0;
    if (!(input >> interface_count))
        return Set_Voronoi_Error(error, "missing interface count");
    if (interface_count < 1)
        return Set_Voronoi_Error(error, "interface count must be at least one");

    std::vector<VORONOI_INTERFACE_RECORD> parsed_interfaces;
    parsed_interfaces.reserve(static_cast<std::size_t>(interface_count));
    std::vector<std::vector<int>> parsed_incident(parsed_milestone_count);
    std::set<std::string> seen_interface_names;
    std::set<std::pair<int, int>> seen_edges;
    int parsed_source_interface = -1;
    for (int index = 0; index < interface_count; ++index)
    {
        VORONOI_INTERFACE_RECORD record;
        if (!(input >> record.name >> record.first >> record.second))
        {
            return Set_Voronoi_Error(
                error,
                "truncated interface record at index " + std::to_string(index));
        }
        if (!Is_Safe_Interface_Name(record.name))
        {
            return Set_Voronoi_Error(error,
                                     "interface name '" + record.name +
                                         "' is not a filename-safe identifier");
        }
        if (!seen_interface_names.insert(record.name).second)
        {
            return Set_Voronoi_Error(
                error, "duplicate interface name '" + record.name + "'");
        }
        if (record.first < 0 || record.first >= parsed_milestone_count ||
            record.second < 0 || record.second >= parsed_milestone_count ||
            record.first == record.second)
        {
            return Set_Voronoi_Error(
                error, "invalid endpoints for interface '" + record.name + "'");
        }
        const std::pair<int, int> edge =
            Canonical_Edge(record.first, record.second);
        if (!seen_edges.insert(edge).second)
        {
            return Set_Voronoi_Error(
                error, "duplicate undirected interface for milestones " +
                           std::to_string(edge.first) + " and " +
                           std::to_string(edge.second));
        }
        parsed_incident[record.first].push_back(index);
        parsed_incident[record.second].push_back(index);
        if (record.name == source_interface_name)
            parsed_source_interface = index;
        parsed_interfaces.push_back(std::move(record));
    }

    std::string trailing_token;
    if (input >> trailing_token)
    {
        return Set_Voronoi_Error(
            error, "unexpected trailing token '" + trailing_token + "'");
    }
    if (input.bad())
        return Set_Voronoi_Error(error,
                                 "I/O error while reading milestone file");
    if (parsed_source_interface < 0)
    {
        return Set_Voronoi_Error(error,
                                 "source_interface '" + source_interface_name +
                                     "' is not present in the milestone file");
    }

    CV_numbers = cv_numbers;
    milestone_count = parsed_milestone_count;
    milestone_cvs = std::move(parsed_cvs);
    interfaces = std::move(parsed_interfaces);
    incident_interfaces = std::move(parsed_incident);
    source_interface = parsed_source_interface;
    current_milestone = -1;
    initial_state_observed = false;
    terminal_hit = false;
    initial_step = -1;
    hit_step = -1;
    hit_from_milestone = -1;
    destination_milestone = -1;
    destination_interface = -1;
    source_recrossing_count = 0;
    return true;
}

bool VORONOI_DETECTOR::Observe_CV_Values(const std::vector<double>& values,
                                         int step, bool commit_sampling_state,
                                         std::string* error)
{
    if (error != NULL) error->clear();
    if (!commit_sampling_state || terminal_hit) return true;
    if (CV_numbers < 1 || milestone_count < 2 || source_interface < 0 ||
        source_interface >= static_cast<int>(interfaces.size()))
        return Set_Voronoi_Error(error, "detector geometry is not initialized");
    if (values.size() != static_cast<std::size_t>(CV_numbers))
        return Set_Voronoi_Error(error,
                                 "runtime CV count does not match the "
                                 "milestone manifest");
    if (step < 0)
        return Set_Voronoi_Error(error, "observation step must be nonnegative");
    for (int cv = 0; cv < CV_numbers; ++cv)
    {
        if (!Is_Finite_Double(values[cv]))
        {
            return Set_Voronoi_Error(error, "runtime CV " + std::to_string(cv) +
                                                " is non-finite at step " +
                                                std::to_string(step));
        }
    }

    std::vector<double> distances(milestone_count, 0.0);
    double best_distance = 0.0;
    int best_milestone = -1;
    int best_count = 0;
    for (int milestone = 0; milestone < milestone_count; ++milestone)
    {
        double distance = 0.0;
        for (int cv = 0; cv < CV_numbers; ++cv)
        {
            const double delta =
                values[cv] -
                milestone_cvs[static_cast<std::size_t>(milestone) * CV_numbers +
                              cv];
            distance += delta * delta;
        }
        if (!Is_Finite_Double(distance))
        {
            return Set_Voronoi_Error(
                error, "non-finite squared distance to milestone " +
                           std::to_string(milestone) + " at step " +
                           std::to_string(step));
        }
        distances[milestone] = distance;
        if (best_milestone < 0 || distance < best_distance)
        {
            best_distance = distance;
            best_milestone = milestone;
            best_count = 1;
        }
        else if (distance == best_distance)
        {
            ++best_count;
        }
    }

    if (!initial_state_observed)
    {
        if (best_count != 1)
        {
            return Set_Voronoi_Error(
                error, "initial state is tied between Voronoi cells");
        }
        const VORONOI_INTERFACE_RECORD& source = interfaces[source_interface];
        if (best_milestone != source.first && best_milestone != source.second)
        {
            return Set_Voronoi_Error(
                error, "initial state belongs to milestone " +
                           std::to_string(best_milestone) +
                           ", which is not incident on source_interface '" +
                           source.name + "'");
        }
        current_milestone = best_milestone;
        initial_state_observed = true;
        initial_step = step;
        return true;
    }

    const auto find_interface = [this](int first, int second)
    {
        for (int interface_index : incident_interfaces[first])
        {
            const VORONOI_INTERFACE_RECORD& record =
                interfaces[interface_index];
            if (record.first == second || record.second == second)
                return interface_index;
        }
        return -1;
    };

    // Contact with the source interface is a recrossing only after a strict
    // move to its other side.  Exact contact with one different interface is
    // already a first hit; contact with several different interfaces has no
    // unique destination and must fail closed.
    if (distances[current_milestone] == best_distance)
    {
        int tied_interface = -1;
        int tied_milestone = -1;
        for (int milestone = 0; milestone < milestone_count; ++milestone)
        {
            if (milestone == current_milestone ||
                distances[milestone] != best_distance)
                continue;
            const int interface_index =
                find_interface(current_milestone, milestone);
            if (interface_index < 0)
            {
                return Set_Voronoi_Error(
                    error, "nonadjacent Voronoi tie between milestones " +
                               std::to_string(current_milestone) + " and " +
                               std::to_string(milestone) + " at step " +
                               std::to_string(step));
            }
            if (interface_index == source_interface) continue;
            if (tied_interface >= 0)
            {
                return Set_Voronoi_Error(
                    error,
                    "ambiguous destination: multiple non-source "
                    "interfaces are tied at step " +
                        std::to_string(step));
            }
            tied_interface = interface_index;
            tied_milestone = milestone;
        }
        if (tied_interface < 0) return true;
        terminal_hit = true;
        hit_step = step;
        hit_from_milestone = current_milestone;
        destination_milestone = tied_milestone;
        destination_interface = tied_interface;
        return true;
    }
    if (best_count != 1)
    {
        return Set_Voronoi_Error(
            error,
            "ambiguous destination: multiple Voronoi cells are tied "
            "at step " +
                std::to_string(step));
    }

    const int crossed_interface =
        find_interface(current_milestone, best_milestone);
    if (crossed_interface < 0)
    {
        return Set_Voronoi_Error(
            error, "nonadjacent Voronoi jump from milestone " +
                       std::to_string(current_milestone) + " to " +
                       std::to_string(best_milestone) + " at step " +
                       std::to_string(step));
    }

    const int previous_milestone = current_milestone;
    if (crossed_interface == source_interface)
    {
        current_milestone = best_milestone;
        ++source_recrossing_count;
        return true;
    }

    terminal_hit = true;
    hit_step = step;
    hit_from_milestone = previous_milestone;
    destination_milestone = best_milestone;
    destination_interface = crossed_interface;
    return true;
}

void VORONOI_DETECTOR::Initial(CONTROLLER* controller,
                               COLLECTIVE_VARIABLE_CONTROLLER* cv_controller)
{
    strcpy(module_name, "voronoi_detector");
    if (!cv_controller->Command_Exist("voronoi_detector_CV")) return;

    controller->printf("START INITIALIZING VORONOI DETECTOR:\n");
    if (CONTROLLER::MPI_size != 1)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorNotImplemented, "VORONOI_DETECTOR::Initial",
            "Reason:\n\tvoronoi_detector requires exactly one MPI process; "
            "got %d\n",
            CONTROLLER::MPI_size);
    }

    cv_list = cv_controller->Ask_For_CV("voronoi_detector", -1);
    const std::vector<std::string> milestone_file =
        cv_controller->Ask_For_String_Parameter(
            "voronoi_detector", "milestone_file", 1, 1, true, "", 0);
    const std::vector<std::string> source_name =
        cv_controller->Ask_For_String_Parameter(
            "voronoi_detector", "source_interface", 1, 1, true, "", 0);
    std::string error;
    if (!Load_Milestone_File(milestone_file[0], source_name[0],
                             static_cast<int>(cv_list.size()), &error))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat, "VORONOI_DETECTOR::Initial",
            "Reason:\n\tinvalid Voronoi milestone configuration: %s\n",
            error.c_str());
    }

    controller->Step_Print_Initial(module_name, "%d");
    is_controller_printf_initialized = 1;
    is_initialized = 1;
    const VORONOI_INTERFACE_RECORD& source = interfaces[source_interface];
    controller->printf("    source interface %s (%d, %d)\n",
                       source.name.c_str(), source.first, source.second);
    controller->printf("END INITIALIZING VORONOI DETECTOR\n\n");
}

void VORONOI_DETECTOR::Observe(int atom_numbers, VECTOR* crd,
                               const LTMatrix3 cell, const LTMatrix3 rcell,
                               int step, bool commit_sampling_state,
                               CONTROLLER* controller)
{
    if (!is_initialized || !commit_sampling_state || terminal_hit) return;

    std::vector<double> values(CV_numbers, 0.0);
    for (int cv = 0; cv < CV_numbers; ++cv)
    {
        cv_list[cv]->Compute(atom_numbers, crd, cell, rcell, CV_NEED_CPU_VALUE,
                             step);
        values[cv] = static_cast<double>(cv_list[cv]->value);
    }

    std::string error;
    if (!Observe_CV_Values(values, step, true, &error))
    {
        std::ostringstream values_text;
        for (int cv = 0; cv < CV_numbers; ++cv)
        {
            if (cv != 0) values_text << ',';
            values_text << values[cv];
        }
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "VORONOI_DETECTOR::Observe",
            "Reason:\n\t%s; CV=[%s]\n", error.c_str(),
            values_text.str().c_str());
    }
}

void VORONOI_DETECTOR::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized || !is_controller_printf_initialized) return;
    controller->Step_Print(module_name, current_milestone);
}

std::string VORONOI_DETECTOR::Hit_Restart_Basename() const
{
    if (!terminal_hit || destination_interface < 0 ||
        destination_interface >= static_cast<int>(interfaces.size()))
        return std::string();
    return "voronoi_hit_" + interfaces[destination_interface].name;
}
