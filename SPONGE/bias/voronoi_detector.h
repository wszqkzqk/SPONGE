#pragma once
#include <vector>

#include "../collective_variable/collective_variable.h"
#include "../MD_core/MD_core.h"

struct VORONOI_DETECTOR
{
    char module_name[CHAR_LENGTH_MAX];
    int is_initialized = 0;
    int is_controller_printf_initialized = 0;
    int last_modify_date = 20260629;

    int CV_numbers = 0;
    CV_LIST cv_list;

    int milestone_count = 0;
    std::vector<float> milestone_cvs;  // flattened [milestone_count][CV_numbers]

    std::vector<std::vector<int>> neighbors;  // adjacency list

    int start_milestone = -1;
    int current_milestone = -1;
    int crossing_detected = 0;
    int crossing_step = -1;
    int destination_milestone = -1;
    char crossing_interface[CHAR_LENGTH_MAX];
    char restart_name[CHAR_LENGTH_MAX];

    void Initial(CONTROLLER* controller,
                 COLLECTIVE_VARIABLE_CONTROLLER* cv_controller);
    void Detect(int atom_numbers, VECTOR* crd, LTMatrix3 cell,
                LTMatrix3 rcell, int step, MD_INFORMATION* md_info,
                CONTROLLER* controller);
    void Step_Print(CONTROLLER* controller);
};
