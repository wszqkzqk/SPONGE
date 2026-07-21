#pragma once

static __host__ __device__ __forceinline__ bool
Minimization_Float_Is_Zero_Or_Normal(float value)
{
#if defined(__CUDA_ARCH__) || \
    (defined(__HIP_DEVICE_COMPILE__) && __HIP_DEVICE_COMPILE__)
    const unsigned int bits = __float_as_uint(value);
#else
    unsigned int bits = 0;
    memcpy(&bits, &value, sizeof(bits));
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : "+r"(bits));
#endif
#endif
    const unsigned int magnitude = bits & 0x7fffffffU;
    const unsigned int exponent = magnitude & 0x7f800000U;
    return magnitude == 0U || (exponent != 0U && exponent != 0x7f800000U);
}

static __host__ __device__ __forceinline__ bool
Minimization_Vector_Is_Zero_Or_Normal(const VECTOR& value)
{
    return Minimization_Float_Is_Zero_Or_Normal(value.x) &&
           Minimization_Float_Is_Zero_Or_Normal(value.y) &&
           Minimization_Float_Is_Zero_Or_Normal(value.z);
}

static __host__ __device__ __forceinline__ VECTOR
Minimization_Limit_Move(const VECTOR& move, float max_move)
{
    if (max_move <= 0.0f) return move;
    const double move_x = static_cast<double>(move.x);
    const double move_y = static_cast<double>(move.y);
    const double move_z = static_cast<double>(move.z);
    const double norm =
        sqrt(move_x * move_x + move_y * move_y + move_z * move_z);
    if (!(norm > static_cast<double>(max_move))) return move;
    const double scale = static_cast<double>(max_move) / norm;
    return {static_cast<float>(move_x * scale),
            static_cast<float>(move_y * scale),
            static_cast<float>(move_z * scale)};
}

static __global__ void MD_Iteration_Gradient_Descent(
    const int atom_numbers, const int* atom_local, VECTOR* crd,
    const VECTOR* frc, const float* mass_inverse, const float dt, VECTOR* vel,
    const float momentum_keep, const float max_move, int* invalid_atom)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        const VECTOR old_crd = crd[i];
        const VECTOR old_vel = vel[i];
        const VECTOR force = frc[i];
        const float inverse_mass = mass_inverse[i];
        bool valid = Minimization_Vector_Is_Zero_Or_Normal(old_crd) &&
                     Minimization_Vector_Is_Zero_Or_Normal(old_vel) &&
                     Minimization_Vector_Is_Zero_Or_Normal(force) &&
                     Minimization_Float_Is_Zero_Or_Normal(inverse_mass) &&
                     inverse_mass >= 0.0f;

        const VECTOR optimizer_velocity = old_vel + dt * inverse_mass * force;
        const VECTOR raw_move = dt * optimizer_velocity;
        const VECTOR move = Minimization_Limit_Move(raw_move, max_move);
        const VECTOR new_crd = old_crd + move;
        const VECTOR retained_velocity = momentum_keep * optimizer_velocity;
        valid = valid &&
                Minimization_Vector_Is_Zero_Or_Normal(optimizer_velocity) &&
                Minimization_Vector_Is_Zero_Or_Normal(raw_move) &&
                Minimization_Vector_Is_Zero_Or_Normal(move) &&
                Minimization_Vector_Is_Zero_Or_Normal(new_crd) &&
                Minimization_Vector_Is_Zero_Or_Normal(retained_velocity);
        if (!valid)
        {
            atomicExch(invalid_atom, atom_local[i]);
        }
        else
        {
            crd[i] = new_crd;
            vel[i] = retained_velocity;
        }
    }
}

static __host__ __device__ __forceinline__ bool Get_Adam_Move_Component(
    float force, float mass_inverse, float beta1, float beta2, float epsilon,
    double first_bias, double second_bias_sqrt, double learning_rate,
    float* first_moment, float* root_second_moment, double* move)
{
    if (!Minimization_Float_Is_Zero_Or_Normal(force) ||
        !Minimization_Float_Is_Zero_Or_Normal(mass_inverse) ||
        mass_inverse < 0.0f ||
        !Minimization_Float_Is_Zero_Or_Normal(*first_moment) ||
        !Minimization_Float_Is_Zero_Or_Normal(*root_second_moment) ||
        *root_second_moment < 0.0f)
    {
        return false;
    }

    // Zero inverse mass marks a derived/non-independent site.  It must not be
    // moved directly, but positive inverse masses are deliberately not used
    // as a preconditioner: coordinate-space Adam is independent of the
    // arbitrary atom masses.
    if (mass_inverse == 0.0f)
    {
        *move = 0.0;
        return true;
    }

    // Store sqrt(v), rather than v itself.  The exact recurrence is evaluated
    // in double precision, so every accepted zero-or-normal float force has a
    // representable state instead of overflowing at force * force.
    const double force_double = static_cast<double>(force);
    const double moment_double =
        static_cast<double>(beta1) * static_cast<double>(*first_moment) +
        (1.0 - static_cast<double>(beta1)) * force_double;
    const double old_root = static_cast<double>(*root_second_moment);
    const double root_double =
        sqrt(static_cast<double>(beta2) * old_root * old_root +
             (1.0 - static_cast<double>(beta2)) * force_double * force_double);
    const float stored_moment = static_cast<float>(moment_double);
    const float stored_root = static_cast<float>(root_double);
    if (!Minimization_Float_Is_Zero_Or_Normal(stored_moment) ||
        !Minimization_Float_Is_Zero_Or_Normal(stored_root) ||
        stored_root < 0.0f)
    {
        return false;
    }

    *first_moment = stored_moment;
    *root_second_moment = stored_root;
    const double corrected_moment = moment_double / first_bias;
    const double corrected_root = root_double / second_bias_sqrt;
    *move = learning_rate * corrected_moment /
            (corrected_root + static_cast<double>(epsilon));
    return true;
}

static __global__ void Get_Adam_Move(int atom_numbers, const int* atom_local,
                                     const float* mass_inverse,
                                     const VECTOR* frc, VECTOR* first_moment,
                                     VECTOR* root_second_moment, float beta1,
                                     float beta2, float epsilon,
                                     float learning_rate, float max_move,
                                     double first_bias, double second_bias_sqrt,
                                     VECTOR* move, int* invalid_atom)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        VECTOR moment = first_moment[i];
        VECTOR root = root_second_moment[i];
        double move_x = 0.0, move_y = 0.0, move_z = 0.0;
        bool valid = first_bias > 0.0 && second_bias_sqrt > 0.0;
        valid = valid && Get_Adam_Move_Component(
                             frc[i].x, mass_inverse[i], beta1, beta2, epsilon,
                             first_bias, second_bias_sqrt, learning_rate,
                             &moment.x, &root.x, &move_x);
        valid = valid && Get_Adam_Move_Component(
                             frc[i].y, mass_inverse[i], beta1, beta2, epsilon,
                             first_bias, second_bias_sqrt, learning_rate,
                             &moment.y, &root.y, &move_y);
        valid = valid && Get_Adam_Move_Component(
                             frc[i].z, mass_inverse[i], beta1, beta2, epsilon,
                             first_bias, second_bias_sqrt, learning_rate,
                             &moment.z, &root.z, &move_z);

        const double move_norm =
            sqrt(move_x * move_x + move_y * move_y + move_z * move_z);
        if (valid && max_move > 0.0f && move_norm > max_move)
        {
            const double scale = static_cast<double>(max_move) / move_norm;
            move_x *= scale;
            move_y *= scale;
            move_z *= scale;
        }
        const VECTOR stored_move = {static_cast<float>(move_x),
                                    static_cast<float>(move_y),
                                    static_cast<float>(move_z)};
        valid = valid && Minimization_Vector_Is_Zero_Or_Normal(stored_move);
        if (!valid)
        {
            atomicExch(invalid_atom, atom_local[i]);
            move[i] = {0.0f, 0.0f, 0.0f};
        }
        else
        {
            first_moment[i] = moment;
            root_second_moment[i] = root;
            move[i] = stored_move;
        }
    }
}

static __global__ void MD_Iteration_Adam_Move(const int atom_numbers,
                                              const int* atom_local,
                                              VECTOR* crd, const VECTOR* move,
                                              int* invalid_atom)
{
    SIMPLE_DEVICE_FOR(i, atom_numbers)
    {
        const VECTOR old_crd = crd[i];
        const VECTOR atom_move = move[i];
        const VECTOR new_crd = old_crd + atom_move;
        if (!Minimization_Vector_Is_Zero_Or_Normal(old_crd) ||
            !Minimization_Vector_Is_Zero_Or_Normal(atom_move) ||
            !Minimization_Vector_Is_Zero_Or_Normal(new_crd))
        {
            atomicExch(invalid_atom, atom_local[i]);
        }
        else
        {
            crd[i] = new_crd;
        }
    }
}

void MD_INFORMATION::MINIMIZATION_iteration::Initial(CONTROLLER* controller,
                                                     MD_INFORMATION* md_info)
{
    this->md_info = md_info;
    this->controller = controller;
    if (md_info->mode == MINIMIZATION)
    {
        controller->printf("    Start initializing minimization:\n");
        max_move = 0.1f;
        if (controller[0].Command_Exist("minimization_max_move"))
        {
            controller->Check_Float(
                "minimization", "max_move",
                "MD_INFORMATION::MINIMIZATION_iteration::Initial");
            max_move = atof(controller[0].Command("minimization_max_move"));
        }
        controller->printf("        minimization max move is %f A\n", max_move);

        momentum_keep = 0;
        if (controller[0].Command_Exist("minimization_momentum_keep"))
        {
            controller->Check_Float(
                "minimization", "momentum_keep",
                "MD_INFORMATION::MINIMIZATION_iteration::Initial");
            momentum_keep =
                atof(controller[0].Command("minimization_momentum_keep"));
        }
        controller->printf("        minimization momentum keep is %f\n",
                           momentum_keep);

        dynamic_dt = 1;
        if (controller[0].Command_Exist("minimization_dynamic_dt"))
        {
            controller->Check_Int(
                "minimization", "dynamic_dt",
                "MD_INFORMATION::MINIMIZATION_iteration::Initial");
            dynamic_dt = atoi(controller[0].Command("minimization_dynamic_dt"));
        }
        controller->printf("        minimization dynamic dt is %d\n",
                           dynamic_dt);

        if ((dynamic_dt != 0 && dynamic_dt != 1) ||
            !Float_Memory_Is_Finite(&max_move) || max_move < 0.0f)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "MD_INFORMATION::MINIMIZATION_iteration::Initial",
                "Reason:\n\tminimization_dynamic_dt must be 0 or 1, and "
                "minimization_max_move must be finite and non-negative\n");
            return;
        }
        if (!Float_Memory_Is_Zero_Or_Normal(&max_move))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand,
                "MD_INFORMATION::MINIMIZATION_iteration::Initial",
                "Reason:\n\tminimization_max_move must be zero or a normal "
                "float so its meaning is identical on FTZ backends\n");
            return;
        }

        if (dynamic_dt)
        {
            md_info->dt = 3e-4f;
            // Adam owns dedicated first/second-moment state.  Physical
            // velocity and acceleration remain available to constraints,
            // kinetic observables, and trajectory/restart output.
            momentum_keep = 0;
            beta1 = 0.9;
            if (controller->Command_Exist("minimization", "beta1"))
            {
                controller->Check_Float(
                    "minimization", "beta1",
                    "MD_INFORMATION::MINIMIZATION_iteration::Initial");
                beta1 = atof(controller->Command("minimization", "beta1"));
            }
            controller->printf("        minimization beta1 is %f\n", beta1);

            beta2 = 0.9;
            if (controller->Command_Exist("minimization", "beta2"))
            {
                controller->Check_Float(
                    "minimization", "beta2",
                    "MD_INFORMATION::MINIMIZATION_iteration::Initial");
                beta2 = atof(controller->Command("minimization", "beta2"));
            }
            controller->printf("        minimization beta2 is %f\n", beta2);

            epsilon = 1e-4f;
            if (controller->Command_Exist("minimization", "epsilon"))
            {
                controller->Check_Float(
                    "minimization", "epsilon",
                    "MD_INFORMATION::MINIMIZATION_iteration::Initial");
                epsilon = atof(controller->Command("minimization", "epsilon"));
            }
            controller->printf("        minimization epsilon is %e\n", epsilon);

            learning_rate = 3e-4f;
            if (controller->Command_Exist("minimization", "learning_rate"))
            {
                controller->Check_Float(
                    "minimization", "learning_rate",
                    "MD_INFORMATION::MINIMIZATION_iteration::Initial");
                learning_rate =
                    atof(controller->Command("minimization", "learning_rate"));
            }
            controller->printf("        minimization learning rate is %e A\n",
                               learning_rate);

            if (!Float_Memory_Is_Finite(&beta1) || beta1 < 0.0f ||
                !(beta1 < 1.0f) || !Float_Memory_Is_Finite(&beta2) ||
                beta2 < 0.0f || !(beta2 < 1.0f) ||
                !Float_Memory_Is_Finite(&epsilon) || !(epsilon > 0.0f) ||
                !Float_Memory_Is_Finite(&learning_rate) ||
                !(learning_rate > 0.0f) ||
                !Float_Memory_Is_Zero_Or_Normal(&beta1) ||
                !Float_Memory_Is_Zero_Or_Normal(&beta2) ||
                !Float_Memory_Is_Normal(&epsilon) ||
                !Float_Memory_Is_Normal(&learning_rate))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "MD_INFORMATION::MINIMIZATION_iteration::Initial",
                    "Reason:\n\tAdam minimization requires 0 <= beta1,beta2 "
                    "< 1, epsilon > 0, and learning_rate > 0; beta values "
                    "must be zero "
                    "or normal floats and positive values must be normal "
                    "floats so FTZ backends use the same parameters\n");
                return;
            }
        }
        else
        {
            md_info->dt = 1e-8f;
            if (!Float_Memory_Is_Finite(&momentum_keep) ||
                momentum_keep < 0.0f || momentum_keep > 1.0f ||
                !Float_Memory_Is_Zero_Or_Normal(&momentum_keep))
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorValueErrorCommand,
                    "MD_INFORMATION::MINIMIZATION_iteration::Initial",
                    "Reason:\n\tminimization_momentum_keep must be finite "
                    "and within [0, 1], and any nonzero value must be a "
                    "normal float\n");
                return;
            }
        }
        Device_Malloc_Safely((void**)&d_invalid_atom, sizeof(int));
        controller->printf("    End initializing minimization\n\n");
    }
}

void MD_INFORMATION::MINIMIZATION_iteration::Gradient_Descent(
    int atom_numbers, const int* atom_local, VECTOR* crd, const VECTOR* frc,
    VECTOR* vel, const float* d_mass_inverse, const VECTOR* adam_move)
{
    if (atom_numbers == 0) return;
    if (atom_numbers < 0 || atom_local == NULL || crd == NULL || frc == NULL ||
        vel == NULL || d_mass_inverse == NULL || d_invalid_atom == NULL ||
        (dynamic_dt && adam_move == NULL))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "MD_INFORMATION::MINIMIZATION_iteration::Gradient_Descent",
            "Reason:\n\tminimization requires valid local-atom, coordinate, "
            "force, velocity, inverse-mass, numeric-error, and (for Adam) "
            "coordinate-move buffers\n");
        return;
    }

    int invalid_atom = -1;
    deviceMemcpy(d_invalid_atom, &invalid_atom, sizeof(int),
                 deviceMemcpyHostToDevice);
    if (dynamic_dt)
    {
        Launch_Device_Kernel(
            MD_Iteration_Adam_Move,
            (atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers, atom_local,
            crd, adam_move, d_invalid_atom);
    }
    else
    {
        Launch_Device_Kernel(
            MD_Iteration_Gradient_Descent,
            (atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers, atom_local,
            crd, frc, d_mass_inverse, md_info->dt, vel, momentum_keep, max_move,
            d_invalid_atom);
    }
    deviceMemcpy(&invalid_atom, d_invalid_atom, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (invalid_atom >= 0)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "MD_INFORMATION::MINIMIZATION_iteration::Gradient_Descent",
            "Reason:\n\tminimization coordinate update input, optimizer "
            "state, displacement, retained velocity, or resulting coordinate "
            "is not a finite zero or normal float for global atom %d at step "
            "%d\n",
            invalid_atom, md_info->sys.steps);
        return;
    }
}

void MD_INFORMATION::MINIMIZATION_iteration::Validate_Final_Time_Step()
{
    if (md_info->mode != MINIMIZATION) return;
    if (!Float_Memory_Is_Finite(&md_info->dt) || !(md_info->dt > 0.0f) ||
        !Float_Memory_Is_Normal(&md_info->dt))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand,
            "MD_INFORMATION::MINIMIZATION_iteration::"
            "Validate_Final_Time_Step",
            "Reason:\n\tminimization dt must be a finite positive normal "
            "float after unit conversion\n");
    }
}

void MD_INFORMATION::MINIMIZATION_iteration::Scale_Force_For_Dynamic_Dt(
    int atom_numbers, const int* atom_local, const float* d_mass_inverse,
    const VECTOR* frc, VECTOR* first_moment, VECTOR* root_second_moment,
    VECTOR* adam_move)
{
    if (md_info->mode == MINIMIZATION && dynamic_dt)
    {
        if (atom_numbers == 0) return;
        if (atom_local == NULL || d_mass_inverse == NULL || frc == NULL ||
            first_moment == NULL || root_second_moment == NULL ||
            adam_move == NULL || d_invalid_atom == NULL)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MD_INFORMATION::MINIMIZATION_iteration::Scale_Force_For_"
                "Dynamic_Dt",
                "Reason:\n\tdynamic minimization requires valid local atom, "
                "mass, force, optimizer-state, and coordinate-move buffers\n");
            return;
        }
        int invalid_atom = -1;
        const double bias_step = static_cast<double>(md_info->sys.steps) + 1.0;
        const double first_bias =
            1.0 - pow(static_cast<double>(beta1), bias_step);
        const double second_bias =
            1.0 - pow(static_cast<double>(beta2), bias_step);
        const double second_bias_sqrt = sqrt(second_bias);
        if (!(first_bias > 0.0) || !(second_bias > 0.0) ||
            !(second_bias_sqrt > 0.0))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MD_INFORMATION::MINIMIZATION_iteration::Scale_Force_For_"
                "Dynamic_Dt",
                "Reason:\n\tAdam bias correction is invalid at step %d\n",
                md_info->sys.steps);
            return;
        }
        deviceMemcpy(d_invalid_atom, &invalid_atom, sizeof(int),
                     deviceMemcpyHostToDevice);
        Launch_Device_Kernel(
            Get_Adam_Move,
            (atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, atom_numbers, atom_local,
            d_mass_inverse, frc, first_moment, root_second_moment, beta1, beta2,
            epsilon, learning_rate, max_move, first_bias, second_bias_sqrt,
            adam_move, d_invalid_atom);
        deviceMemcpy(&invalid_atom, d_invalid_atom, sizeof(int),
                     deviceMemcpyDeviceToHost);
        if (invalid_atom >= 0)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MD_INFORMATION::MINIMIZATION_iteration::Scale_Force_For_"
                "Dynamic_Dt",
                "Reason:\n\tAdam optimizer input, state, or coordinate move "
                "is not a finite zero or normal float for global atom %d at "
                "step %d\n",
                invalid_atom, md_info->sys.steps);
            return;
        }
    }
}

void MD_INFORMATION::MINIMIZATION_iteration::Clear()
{
    Free_Single_Device_Pointer((void**)&d_invalid_atom);
}
