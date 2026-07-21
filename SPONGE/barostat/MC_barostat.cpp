#include "MC_barostat.h"

void MC_BAROSTAT_INFORMATION::Volume_Change_Attempt(VECTOR boxlength, float dt)
{
    if (CONTROLLER::MPI_rank == 0)
    {
        if (!Float_Memory_Is_Finite(&dt) || !(dt > 0.0f))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MC_BAROSTAT_INFORMATION::Volume_Change_Attempt",
                "Reason:\n\tthe MC barostat timestep must be finite and "
                "positive\n");
        }
        const float box_components[3] = {boxlength.x, boxlength.y, boxlength.z};
        for (int axis = 0; axis < 3; ++axis)
        {
            if (!Float_Memory_Is_Finite(box_components + axis) ||
                !(box_components[axis] > 0.0f) ||
                !Float_Memory_Is_Finite(Delta_Box_Length_Max + axis) ||
                Delta_Box_Length_Max[axis] < 0.0f)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "MC_BAROSTAT_INFORMATION::Volume_Change_Attempt",
                    "Reason:\n\tinvalid MC box/proposal extent on axis %d\n",
                    axis);
            }
        }
        double nrand = ((double)2.0 * rand() / RAND_MAX - 1.0);

        Delta_Box_Length = {0.0f, 0.0f, 0.0f};
        switch (couple_dimension)
        {
            case NO:
                if (only_direction > 0)
                    xyz = only_direction - 1;
                else
                    xyz = rand() % 3;
                if (xyz == 0)
                {
                    Delta_Box_Length.x = nrand * Delta_Box_Length_Max[xyz];
                }
                else if (xyz == 1)
                {
                    Delta_Box_Length.y = nrand * Delta_Box_Length_Max[xyz];
                }
                else
                {
                    Delta_Box_Length.z = nrand * Delta_Box_Length_Max[xyz];
                }
                break;
            case XY:
                if (only_direction > 0)
                    xyz = only_direction - 1;
                else
                    xyz = rand() % 2;
                if (xyz == 0)
                {
                    Delta_Box_Length.z = nrand * Delta_Box_Length_Max[xyz];
                }
                else
                {
                    Delta_Box_Length.y = nrand * Delta_Box_Length_Max[xyz];
                    Delta_Box_Length.x = nrand * Delta_Box_Length_Max[xyz];
                }
                break;
            case XZ:
                if (only_direction > 0)
                    xyz = only_direction - 1;
                else
                    xyz = rand() % 2;
                if (xyz == 0)
                {
                    Delta_Box_Length.y = nrand * Delta_Box_Length_Max[xyz];
                }
                else
                {
                    Delta_Box_Length.z = nrand * Delta_Box_Length_Max[xyz];
                    Delta_Box_Length.x = nrand * Delta_Box_Length_Max[xyz];
                }
                break;
            case YZ:
                if (only_direction > 0)
                    xyz = only_direction - 1;
                else
                    xyz = rand() % 2;
                if (xyz == 0)
                {
                    Delta_Box_Length.x = nrand * Delta_Box_Length_Max[xyz];
                }
                else
                {
                    Delta_Box_Length.z = nrand * Delta_Box_Length_Max[xyz];
                    Delta_Box_Length.y = nrand * Delta_Box_Length_Max[xyz];
                }
                break;
            case XYZ:
                xyz = 0;
                Delta_Box_Length.x = nrand * Delta_Box_Length_Max[xyz];
                Delta_Box_Length.y = nrand * Delta_Box_Length_Max[xyz];
                Delta_Box_Length.z = nrand * Delta_Box_Length_Max[xyz];
                break;
        }

        New_Box_Length = boxlength + Delta_Box_Length;
        const float new_box_components[3] = {New_Box_Length.x, New_Box_Length.y,
                                             New_Box_Length.z};
        for (int axis = 0; axis < 3; ++axis)
        {
            if (!Float_Memory_Is_Finite(new_box_components + axis) ||
                !(new_box_components[axis] > 0.0f))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "MC_BAROSTAT_INFORMATION::Volume_Change_Attempt",
                    "Reason:\n\tthe proposed MC box length on axis %d is "
                    "non-finite or non-positive\n",
                    axis);
            }
        }
        DeltaS = 0.0f;
        switch (couple_dimension)
        {
            case NO:
                break;
            case XY:
                if (xyz == 1)
                {
                    DeltaS = New_Box_Length.x * New_Box_Length.y -
                             boxlength.x * boxlength.y;
                }
                break;
            case XZ:
                if (xyz == 1)
                {
                    DeltaS = New_Box_Length.x * New_Box_Length.z -
                             boxlength.x * boxlength.z;
                }
                break;
            case YZ:
                if (xyz == 1)
                {
                    DeltaS = New_Box_Length.z * New_Box_Length.y -
                             boxlength.z * boxlength.y;
                }
                break;
            case XYZ:
                break;
        }
        double V = boxlength.x * boxlength.y * boxlength.z;
        newV = New_Box_Length.x * New_Box_Length.y * New_Box_Length.z;
        DeltaV = newV - V;
        VDevided = newV / V;
        if (!isfinite(V) || !(V > 0.0) || !Float_Memory_Is_Finite(&newV) ||
            !(newV > 0.0f) || !isfinite(DeltaV) || !isfinite(VDevided) ||
            !(VDevided > 0.0))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MC_BAROSTAT_INFORMATION::Volume_Change_Attempt",
                "Reason:\n\tthe proposed MC volume or volume ratio is "
                "invalid\n");
        }
        VECTOR crd_scale_factor = New_Box_Length / boxlength;
        g = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        g.a11 = (crd_scale_factor.x - 1.0f) / dt;
        g.a22 = (crd_scale_factor.y - 1.0f) / dt;
        g.a33 = (crd_scale_factor.z - 1.0f) / dt;
    }
#ifdef USE_MPI
    MPI_Bcast(&g, 6, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
}

bool MC_BAROSTAT_INFORMATION::Will_Attempt(int steps) const
{
    return is_initialized && update_interval > 0 && steps >= 0 &&
           steps % update_interval == 0;
}

LTMatrix3 MC_BAROSTAT_INFORMATION::Get_Exact_Reverse_G(float dt) const
{
    LTMatrix3 reverse = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    const float forward[3] = {g.a11, g.a22, g.a33};
    float* reverse_data[3] = {&reverse.a11, &reverse.a22, &reverse.a33};
    for (int axis = 0; axis < 3; ++axis)
    {
        const double scale =
            1.0 + static_cast<double>(dt) * static_cast<double>(forward[axis]);
        const double value = -static_cast<double>(forward[axis]) / scale;
        const float stored = static_cast<float>(value);
        if (!isfinite(scale) || !(scale > 0.0) || !isfinite(value) ||
            !Float_Memory_Is_Finite(&stored))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MC_BAROSTAT_INFORMATION::Get_Exact_Reverse_G",
                "Reason:\n\tthe MC box scaling on axis %d has no finite "
                "positive inverse\n",
                axis);
        }
        *reverse_data[axis] = stored;
    }
    return reverse;
}

int MC_BAROSTAT_INFORMATION::Check_MC_Barostat_Accept()
{
    if (CONTROLLER::MPI_rank == 0)
    {
        total_count[xyz] += 1;
        const float tmp_rand = (float)rand() / RAND_MAX;
        accept = tmp_rand < accept_possibility;
        if (accept)
        {
            accep_count[xyz] += 1;
        }
    }
#ifdef USE_MPI
    MPI_Bcast(&accept, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    return accept;
}

void MC_BAROSTAT_INFORMATION::Initial(CONTROLLER* controller, int atom_numbers,
                                      float target_pressure, VECTOR boxlength,
                                      LTMatrix3 cell, const char* module_name)
{
    this->controller = controller;
    controller->printf("START INITIALIZING MC BAROSTAT:\n");
    if (module_name == NULL)
    {
        strcpy(this->module_name, "monte_carlo_barostat");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }
    if (!Float_Memory_Is_Finite(&target_pressure))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\tthe MC barostat target pressure must be finite\n");
    }
    const float box_components[3] = {boxlength.x, boxlength.y, boxlength.z};
    for (int axis = 0; axis < 3; ++axis)
    {
        if (!Float_Memory_Is_Finite(box_components + axis) ||
            !(box_components[axis] > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "MC_BAROSTAT_INFORMATION::Initial",
                "Reason:\n\tthe MC barostat box length on axis %d must be "
                "finite and positive\n",
                axis);
        }
    }
    if (cell.a21 != 0.0f || cell.a31 != 0.0f || cell.a32 != 0.0f)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\t MC barostat only supports orthogonal box now.\n");
    }
    controller->printf("    The target pressure is %.2f bar\n",
                       target_pressure * CONSTANT_PRES_CONVERTION);
    float mc_baro_initial_ratio = 0.001;
    if (controller[0].Command_Exist(this->module_name, "initial_ratio"))
    {
        controller->Check_Float(this->module_name, "initial_ratio",
                                "MC_BAROSTAT_INFORMATION::Initial");
        mc_baro_initial_ratio =
            atof(controller[0].Command(this->module_name, "initial_ratio"));
    }
    if (!Float_Memory_Is_Finite(&mc_baro_initial_ratio) ||
        mc_baro_initial_ratio < 0.0f || mc_baro_initial_ratio >= 1.0f)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\tmonte_carlo_barostat_initial_ratio must be finite "
            "and in [0, 1)\n");
    }
    Delta_Box_Length_Max[0] = mc_baro_initial_ratio * boxlength.x;
    Delta_Box_Length_Max[1] = mc_baro_initial_ratio * boxlength.y;
    Delta_Box_Length_Max[2] = mc_baro_initial_ratio * boxlength.z;
    controller->printf(
        "    The initial max box length to change is %f %f %f Angstrom for x y "
        "z\n",
        Delta_Box_Length_Max[0], Delta_Box_Length_Max[1],
        Delta_Box_Length_Max[2]);

    update_interval = 100;
    if (controller[0].Command_Exist(this->module_name, "update_interval"))
    {
        controller->Check_Int(this->module_name, "update_interval",
                              "MC_BAROSTAT_INFORMATION::Initial");
        update_interval =
            atoi(controller[0].Command(this->module_name, "update_interval"));
    }
    if (update_interval <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\tmonte_carlo_barostat_update_interval must be "
            "positive\n");
    }
    controller->printf("    The update_interval is %d\n", update_interval);

    check_interval = 10;
    if (controller[0].Command_Exist(this->module_name, "check_interval"))
    {
        controller->Check_Int(this->module_name, "check_interval",
                              "MC_BAROSTAT_INFORMATION::Initial");
        check_interval =
            atoi(controller[0].Command(this->module_name, "check_interval"));
    }
    if (check_interval <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\tmonte_carlo_barostat_check_interval must be "
            "positive\n");
    }
    controller->printf("    The check_interval is %d\n", check_interval);

    accept_rate_low = 30;
    if (controller[0].Command_Exist(this->module_name, "accept_rate_low"))
    {
        controller->Check_Float(this->module_name, "accept_rate_low",
                                "MC_BAROSTAT_INFORMATION::Initial");
        accept_rate_low =
            atof(controller[0].Command(this->module_name, "accept_rate_low"));
    }
    if (!Float_Memory_Is_Finite(&accept_rate_low) || accept_rate_low < 0.0f ||
        accept_rate_low > 100.0f)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\tmonte_carlo_barostat_accept_rate_low must be finite "
            "and in [0, 100]\n");
    }
    controller->printf("    The lowest accept rate is %.2f%%\n",
                       accept_rate_low);

    accept_rate_high = 40;
    if (controller[0].Command_Exist(this->module_name, "accept_rate_high"))
    {
        controller->Check_Float(this->module_name, "accept_rate_high",
                                "MC_BAROSTAT_INFORMATION::Initial");
        accept_rate_high =
            atof(controller[0].Command(this->module_name, "accept_rate_high"));
    }
    if (!Float_Memory_Is_Finite(&accept_rate_high) ||
        accept_rate_high < accept_rate_low || accept_rate_high > 100.0f)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\tmonte_carlo_barostat_accept_rate_high must be finite, "
            "at least accept_rate_low, and at most 100\n");
    }
    controller->printf("    The highest accept rate is %.2f%%\n",
                       accept_rate_high);

    if (!controller->Command_Exist(this->module_name, "couple_dimension") ||
        controller->Command_Choice(this->module_name, "couple_dimension",
                                   "XYZ"))
    {
        couple_dimension = XYZ;
    }
    else if (controller->Command_Choice(this->module_name, "couple_dimension",
                                        "NO"))
    {
        couple_dimension = NO;
    }
    else if (controller->Command_Choice(this->module_name, "couple_dimension",
                                        "XY"))
    {
        couple_dimension = XY;
    }
    else if (controller->Command_Choice(this->module_name, "couple_dimension",
                                        "XZ"))
    {
        couple_dimension = XZ;
    }
    else if (controller->Command_Choice(this->module_name, "couple_dimension",
                                        "YZ"))
    {
        couple_dimension = YZ;
    }
    else
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "MC_BAROSTAT_INFORMATION::Initial",
            "Reason:\n\tmonte_carlo_barostat_couple_dimension must be one of "
            "NO, XY, XZ, YZ, or XYZ\n");
    }
    if (!controller->Command_Exist(this->module_name, "couple_dimension"))
        controller->printf("    The couple dimension is %s (index %d)\n", "XYZ",
                           couple_dimension);
    else
        controller->printf(
            "    The couple dimension is %s (index %d)\n",
            controller->Command(this->module_name, "couple_dimension"),
            couple_dimension);
    if (controller->Command_Exist(this->module_name, "only_direction"))
    {
        bool direction_is_valid = false;
        if (couple_dimension == NO)
        {
            if (controller->Command_Choice(this->module_name, "only_direction",
                                           "x"))
            {
                only_direction = 1;
                direction_is_valid = true;
            }
            else if (controller->Command_Choice(this->module_name,
                                                "only_direction", "y"))
            {
                only_direction = 2;
                direction_is_valid = true;
            }
            else if (controller->Command_Choice(this->module_name,
                                                "only_direction", "z"))
            {
                only_direction = 3;
                direction_is_valid = true;
            }
        }
        else if (couple_dimension == XYZ)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "MC_BAROSTAT_INFORMATION::Initial",
                "Reason:\n\tonly_direction is not valid for isotropic pressure "
                "regulation\n");
        }
        else if (couple_dimension == XY)
        {
            if (controller->Command_Choice(this->module_name, "only_direction",
                                           "z"))
            {
                only_direction = 1;
                direction_is_valid = true;
            }
            else if (controller->Command_Choice(this->module_name,
                                                "only_direction", "xy"))
            {
                only_direction = 2;
                direction_is_valid = true;
            }
        }
        else if (couple_dimension == XZ)
        {
            if (controller->Command_Choice(this->module_name, "only_direction",
                                           "y"))
            {
                only_direction = 1;
                direction_is_valid = true;
            }
            else if (controller->Command_Choice(this->module_name,
                                                "only_direction", "xz"))
            {
                only_direction = 2;
                direction_is_valid = true;
            }
        }
        else if (couple_dimension == YZ)
        {
            if (controller->Command_Choice(this->module_name, "only_direction",
                                           "x"))
            {
                only_direction = 1;
                direction_is_valid = true;
            }
            else if (controller->Command_Choice(this->module_name,
                                                "only_direction", "yz"))
            {
                only_direction = 2;
                direction_is_valid = true;
            }
        }
        if (!direction_is_valid && couple_dimension != XYZ)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "MC_BAROSTAT_INFORMATION::Initial",
                "Reason:\n\tmonte_carlo_barostat_only_direction is invalid "
                "for the selected coupling mode\n");
        }
    }
    if (couple_dimension != NO && couple_dimension != XYZ)
    {
        surface_number = 0;
        if (controller->Command_Exist(this->module_name, "surface_number"))
        {
            controller->Check_Int(this->module_name, "surface_number",
                                  "MC_BAROSTAT_INFORMATION::Initial");
            surface_number =
                atoi(controller->Command(this->module_name, "surface_number"));
        }
        if (surface_number < 0)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "MC_BAROSTAT_INFORMATION::Initial",
                "Reason:\n\tmonte_carlo_barostat_surface_number must be "
                "non-negative\n");
        }
        surface_tension = 0.0f;
        if (controller->Command_Exist(this->module_name, "surface_tension"))
        {
            controller->Check_Float(this->module_name, "surface_tension",
                                    "MC_BAROSTAT_INFORMATION::Initial");
            surface_tension =
                atof(controller->Command(this->module_name, "surface_tension"));
        }
        if (!Float_Memory_Is_Finite(&surface_tension))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "MC_BAROSTAT_INFORMATION::Initial",
                "Reason:\n\tmonte_carlo_barostat_surface_tension must be "
                "finite\n");
        }
        surface_tension *= TENSION_UNIT_FACTOR;
        controller->printf("        The surface number is %d\n",
                           surface_number);
        controller->printf("        The surface tension is %f\n",
                           surface_tension);
    }
    Device_Malloc_Safely((void**)&crd_backup, sizeof(VECTOR) * atom_numbers);
    is_initialized = 1;
    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial("density", "%.4f");
        is_controller_printf_initialized = 1;
        controller[0].printf("    structure last modify date is %d\n",
                             last_modify_date);
    }
    controller[0].printf("END INITIALIZING MC BAROSTAT\n\n");
}

void MC_BAROSTAT_INFORMATION::Delta_Box_Length_Max_Update()
{
    if (CONTROLLER::MPI_rank != 0) return;
    if (total_count[xyz] % check_interval == 0)
    {
        accept_rate[xyz] = 100.0 * accep_count[xyz] / total_count[xyz];

        if (accept_rate[xyz] < accept_rate_low)
        {
            total_count[xyz] = 0;
            accep_count[xyz] = 0;
            Delta_Box_Length_Max[xyz] *= 0.9;
        }
        if (accept_rate[xyz] > accept_rate_high)
        {
            total_count[xyz] = 0;
            accep_count[xyz] = 0;
            Delta_Box_Length_Max[xyz] *= 1.1;
        }
    }
}

void MC_BAROSTAT_INFORMATION::Ask_For_Calculate_Potential(int steps,
                                                          int* need_potential)
{
    if (is_initialized && steps % update_interval == 0)
    {
        *need_potential = 1;
    }
}
