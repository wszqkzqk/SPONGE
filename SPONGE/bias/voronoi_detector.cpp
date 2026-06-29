#include "voronoi_detector.h"

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

void VORONOI_DETECTOR::Initial(CONTROLLER* controller,
                               COLLECTIVE_VARIABLE_CONTROLLER* cv_controller)
{
    strcpy(this->module_name, "voronoi_detector");
    if (!cv_controller->Command_Exist("voronoi_detector_CV"))
    {
        return;
    }
    controller->printf("START INITIALIZING VORONOI DETECTOR:\n");

    cv_list = cv_controller->Ask_For_CV("voronoi_detector", 0);
    CV_numbers = static_cast<int>(cv_list.size());
    if (CV_numbers < 1)
    {
        controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                       "VORONOI_DETECTOR::Initial",
                                       "Reason:\n\tat least one CV is required\n");
    }

    std::vector<std::string> milestone_file =
        cv_controller->Ask_For_String_Parameter(
            "voronoi_detector", "milestone_file", 1, 1, true);

    std::ifstream infile(milestone_file[0]);
    if (!infile.is_open())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "VORONOI_DETECTOR::Initial",
            (std::string("Reason:\n\tCannot open milestone file ") +
             milestone_file[0])
                .c_str());
    }

    infile >> milestone_count;
    if (milestone_count < 2)
    {
        controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                       "VORONOI_DETECTOR::Initial",
                                       "Reason:\n\tat least two milestones are "
                                       "required\n");
    }

    milestone_cvs.resize(static_cast<size_t>(milestone_count) * CV_numbers);
    neighbors.resize(milestone_count);

    std::string name;
    for (int i = 0; i < milestone_count; ++i)
    {
        infile >> name;
        for (int j = 0; j < CV_numbers; ++j)
        {
            infile >> milestone_cvs[static_cast<size_t>(i) * CV_numbers + j];
        }
    }

    int interface_count = 0;
    infile >> interface_count;
    for (int k = 0; k < interface_count; ++k)
    {
        std::string iface_name;
        int i, j;
        infile >> iface_name >> i >> j;
        if (i < 0 || i >= milestone_count || j < 0 || j >= milestone_count ||
            i == j)
        {
            controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                           "VORONOI_DETECTOR::Initial",
                                           "Reason:\n\tinvalid interface pair\n");
        }
        neighbors[i].push_back(j);
        neighbors[j].push_back(i);
    }
    infile.close();

    start_milestone = cv_controller->Ask_For_Int_Parameter(
        "voronoi_detector", "start_milestone", 1, 1, false, 0, 0)[0];
    if (start_milestone < 0 || start_milestone >= milestone_count)
    {
        controller->Throw_SPONGE_Error(spongeErrorValueErrorCommand,
                                       "VORONOI_DETECTOR::Initial",
                                       "Reason:\n\tstart_milestone out of range\n");
    }

    controller->Step_Print_Initial(module_name, "%d");
    is_controller_printf_initialized = 1;
    is_initialized = 1;
    controller->printf("END INITIALIZING VORONOI DETECTOR\n\n");
}

void VORONOI_DETECTOR::Detect(int atom_numbers, VECTOR* crd, LTMatrix3 cell,
                              LTMatrix3 rcell, int step,
                              MD_INFORMATION* md_info,
                              CONTROLLER* controller)
{
    if (!is_initialized || crossing_detected)
    {
        return;
    }

    std::vector<float> values(CV_numbers);
    for (int i = 0; i < CV_numbers; ++i)
    {
        cv_list[i]->Compute(atom_numbers, crd, cell, rcell,
                            CV_NEED_CPU_VALUE, step);
        values[i] = cv_list[i]->value;
    }

    int nearest = -1;
    float best_d2 = FLT_MAX;
    for (int m = 0; m < milestone_count; ++m)
    {
        float d2 = 0.0f;
        for (int j = 0; j < CV_numbers; ++j)
        {
            float dx = values[j] -
                       milestone_cvs[static_cast<size_t>(m) * CV_numbers + j];
            d2 += dx * dx;
        }
        if (d2 < best_d2)
        {
            best_d2 = d2;
            nearest = m;
        }
    }

    if (current_milestone < 0)
    {
        if (start_milestone >= 0 && nearest != start_milestone)
        {
            // Initial configuration is already in a different cell; treat it as
            // a crossing from the requested start milestone.
            controller->printf(
                "VORONOI_DETECTOR: initial coords are in milestone %d, "
                "crossing from requested start_milestone %d\n",
                nearest, start_milestone);
            current_milestone = start_milestone;
            // Fall through to the crossing logic below.
        }
        else
        {
            current_milestone = nearest;
            return;
        }
    }

    if (nearest != current_milestone)
    {
        bool adjacent = false;
        for (int nb : neighbors[current_milestone])
        {
            if (nb == nearest)
            {
                adjacent = true;
                break;
            }
        }
        if (!adjacent)
        {
            controller->printf(
                "VORONOI_DETECTOR WARNING: jumped from milestone %d to non-"
                "adjacent milestone %d at step %d; ignoring\n",
                current_milestone, nearest, step);
            current_milestone = nearest;
            return;
        }

        crossing_detected = 1;
        crossing_step = step;
        destination_milestone = nearest;

        int i_low = current_milestone < nearest ? current_milestone : nearest;
        int i_high = current_milestone < nearest ? nearest : current_milestone;
        snprintf(crossing_interface, CHAR_LENGTH_MAX, "S_%02d_%02d",
                 i_low + 1, i_high + 1);
        snprintf(restart_name, CHAR_LENGTH_MAX, "voronoi_hit_%s",
                 crossing_interface);

        controller->printf(
            "VORONOI_DETECTOR: crossing detected at step %d, "
            "interface %s (%d -> %d)\n",
            step, crossing_interface, current_milestone, nearest);

        md_info->sys.step_limit = step;
    }
}

void VORONOI_DETECTOR::Step_Print(CONTROLLER* controller)
{
    if (!is_initialized || !is_controller_printf_initialized)
    {
        return;
    }
    controller->Step_Print(module_name, current_milestone);
}
