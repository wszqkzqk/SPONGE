#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "../collective_variable/collective_variable.h"

struct VORONOI_INTERFACE_RECORD
{
    std::string name;
    int first = -1;
    int second = -1;
};

// Detect first passage from a declared Voronoi interface to any other
// interface.  The two cells incident on source_interface form the source
// region: crossings of their shared edge are recrossings, not terminal hits.
struct VORONOI_DETECTOR
{
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260723;

    int CV_numbers = 0;
    CV_LIST cv_list;

    int milestone_count = 0;
    std::vector<double> milestone_cvs;
    std::vector<VORONOI_INTERFACE_RECORD> interfaces;
    std::vector<std::vector<int>> incident_interfaces;

    int source_interface = -1;
    int current_milestone = -1;
    bool initial_state_observed = false;

    bool terminal_hit = false;
    int initial_step = -1;
    int hit_step = -1;
    int hit_from_milestone = -1;
    int destination_milestone = -1;
    int destination_interface = -1;
    std::uint64_t source_recrossing_count = 0;

    // These two methods are separated from CONTROLLER I/O so the parser and
    // state machine can be tested directly.  On failure they do not publish a
    // partial graph or transition.
    bool Load_Milestone_File(const std::string& filename,
                             const std::string& source_interface_name,
                             int cv_numbers, std::string* error);
    bool Observe_CV_Values(const std::vector<double>& values, int step,
                           bool commit_sampling_state, std::string* error);

    void Initial(CONTROLLER* controller,
                 COLLECTIVE_VARIABLE_CONTROLLER* cv_controller);
    void Observe(int atom_numbers, VECTOR* crd, const LTMatrix3 cell,
                 const LTMatrix3 rcell, const LTMatrix3 reference_cell,
                 int step, bool commit_sampling_state, CONTROLLER* controller);
    void Step_Print(CONTROLLER* controller);

    bool Has_Terminal_Hit() const { return terminal_hit; }
    const VORONOI_INTERFACE_RECORD& Hit_Interface() const
    {
        return interfaces[destination_interface];
    }
    std::string Hit_Restart_Basename() const;
};
