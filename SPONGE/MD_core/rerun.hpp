#pragma once

#include <cerrno>

namespace
{

enum class RerunBinaryReadState
{
    complete,
    clean_eof,
    partial,
    io_error,
};

struct RerunBinaryReadResult
{
    RerunBinaryReadState state;
    std::size_t count;
};

static RerunBinaryReadResult Rerun_Read_Binary_Frame(FILE* input, void* output,
                                                     std::size_t element_size,
                                                     std::size_t count)
{
    if (input == NULL || output == NULL || element_size == 0 || count == 0)
    {
        return {RerunBinaryReadState::io_error, 0};
    }
    const std::size_t read_count = fread(output, element_size, count, input);
    if (read_count == count)
    {
        return {RerunBinaryReadState::complete, read_count};
    }
    if (ferror(input))
    {
        return {RerunBinaryReadState::io_error, read_count};
    }
    if (read_count == 0 && feof(input))
    {
        return {RerunBinaryReadState::clean_eof, 0};
    }
    return {RerunBinaryReadState::partial, read_count};
}

static const char* Rerun_Binary_State_Name(RerunBinaryReadState state)
{
    switch (state)
    {
        case RerunBinaryReadState::complete:
            return "complete";
        case RerunBinaryReadState::clean_eof:
            return "clean EOF";
        case RerunBinaryReadState::partial:
            return "partial";
        case RerunBinaryReadState::io_error:
            return "I/O error";
    }
    return "unknown";
}

static int Rerun_Parse_Nonnegative_Int(CONTROLLER* controller,
                                       const char* command,
                                       const char* error_by)
{
    const char* token = controller->Command(command);
    errno = 0;
    char* end = NULL;
    const long parsed = strtol(token, &end, 10);
    if (end == token || end == NULL || *end != '\0' || errno == ERANGE ||
        parsed < 0 || parsed > INT_MAX)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorValueErrorCommand, error_by,
            "Reason:\n\t%s must be a nonnegative int in [0, %d]; got "
            "\"%s\"\n",
            command, INT_MAX, token);
    }
    return static_cast<int>(parsed);
}

static int Rerun_Read_Box_Record(FILE* input, std::size_t* line_number,
                                 VECTOR* box_length, VECTOR* box_angle,
                                 CONTROLLER* controller)
{
    if (input == NULL || line_number == NULL || box_length == NULL ||
        box_angle == NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "MD_INFORMATION::RERUN_information::Read_Frame",
            "Reason:\n\tthe rerun box stream or destination is null\n");
    }

    while (true)
    {
        std::string line;
        bool read_character = false;
        int character = EOF;
        while ((character = getc(input)) != EOF)
        {
            read_character = true;
            if (character == '\n') break;
            try
            {
                line.push_back(static_cast<char>(character));
            }
            catch (const std::length_error&)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow,
                    "MD_INFORMATION::RERUN_information::Read_Frame",
                    "Reason:\n\trerun box line %zu exceeds the host string "
                    "container limit\n",
                    *line_number + 1);
            }
            catch (const std::bad_alloc&)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorMallocFailed,
                    "MD_INFORMATION::RERUN_information::Read_Frame",
                    "Reason:\n\tfailed to grow host storage while reading "
                    "rerun box line %zu\n",
                    *line_number + 1);
            }
        }
        if (character == EOF && ferror(input))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat,
                "MD_INFORMATION::RERUN_information::Read_Frame",
                "Reason:\n\tI/O error while reading rerun box line %zu\n",
                *line_number + 1);
        }
        if (!read_character && line.empty()) return 0;

        ++*line_number;
        if (line.find_first_not_of(" \t\r") == std::string::npos)
        {
            if (character == EOF) return 0;
            continue;
        }

        VECTOR candidate_length;
        VECTOR candidate_angle;
        std::istringstream parser(line);
        if (!(parser >> candidate_length.x >> candidate_length.y >>
              candidate_length.z >> candidate_angle.x >> candidate_angle.y >>
              candidate_angle.z))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat,
                "MD_INFORMATION::RERUN_information::Read_Frame",
                "Reason:\n\trerun box line %zu must contain exactly six "
                "finite float tokens\n",
                *line_number);
        }
        const float components[6] = {candidate_length.x, candidate_length.y,
                                     candidate_length.z, candidate_angle.x,
                                     candidate_angle.y,  candidate_angle.z};
        for (const float component : components)
        {
            if (!Float_Memory_Is_Finite(&component))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorBadFileFormat,
                    "MD_INFORMATION::RERUN_information::Read_Frame",
                    "Reason:\n\trerun box line %zu contains a non-finite "
                    "float\n",
                    *line_number);
            }
        }
        std::string trailing;
        if (parser >> trailing)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat,
                "MD_INFORMATION::RERUN_information::Read_Frame",
                "Reason:\n\trerun box line %zu has trailing token \"%s\"; "
                "each frame must contain exactly six floats\n",
                *line_number, trailing.c_str());
        }
        *box_length = candidate_length;
        *box_angle = candidate_angle;
        return 1;
    }
}

static bool Rerun_Box_Is_Exactly_Equal(VECTOR lhs_length, VECTOR lhs_angle,
                                       VECTOR rhs_length, VECTOR rhs_angle)
{
    return lhs_length.x == rhs_length.x && lhs_length.y == rhs_length.y &&
           lhs_length.z == rhs_length.z && lhs_angle.x == rhs_angle.x &&
           lhs_angle.y == rhs_angle.y && lhs_angle.z == rhs_angle.z;
}

static bool Rerun_Deformation_Is_Representable(const LTMatrix3& deformation)
{
    const float components[6] = {deformation.a11, deformation.a21,
                                 deformation.a22, deformation.a31,
                                 deformation.a32, deformation.a33};
    for (const float component : components)
    {
        if (!Float_Memory_Is_Zero_Or_Normal(&component)) return false;
    }
    return true;
}

}  // namespace

void MD_INFORMATION::RERUN_information::Initial(CONTROLLER* controller,
                                                MD_INFORMATION* md_info)
{
    this->controller = controller;
    this->md_info = md_info;
    if (md_info->mode != RERUN) return;

    controller->printf("    Start initializing rerun:\n");
    if (!controller->Command_Exist(TRAJ_COMMAND))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand,
            "MD_INFORMATION::RERUN_information::Initial",
            "Reason:\n\tno trajectory information found (command 'crd' "
            "required)");
    }
    if (!controller->Command_Exist(BOX_TRAJ_COMMAND))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand,
            "MD_INFORMATION::RERUN_information::Initial",
            "Reason:\n\tno box information found (command 'box' required)");
    }

    has_velocity_stream = controller->Command_Exist(VEL_TRAJ_COMMAND);
    if (controller->MPI_rank == 0)
    {
        const std::string trajectory_name =
            controller->Original_Command(TRAJ_COMMAND);
        Open_File_Safely(&traj_file, trajectory_name.c_str(), "rb");
        controller->printf("        Open rerun coordinate trajectory '%s'\n",
                           trajectory_name.c_str());
        controller->Set_File_Buffer(traj_file,
                                    sizeof(VECTOR) * md_info->atom_numbers);
        Malloc_Safely((void**)&coordinate_staging,
                      sizeof(VECTOR) * md_info->atom_numbers);

        const std::string box_name =
            controller->Original_Command(BOX_TRAJ_COMMAND);
        Open_File_Safely(&box_file, box_name.c_str(), "r");
        controller->printf("        Open rerun box trajectory '%s'\n",
                           box_name.c_str());
        controller->Set_File_Buffer(box_file, sizeof(char) * 50);

        if (has_velocity_stream)
        {
            const std::string velocity_name =
                controller->Original_Command(VEL_TRAJ_COMMAND);
            Open_File_Safely(&vel_file, velocity_name.c_str(), "rb");
            controller->printf("        Open rerun velocity trajectory '%s'\n",
                               velocity_name.c_str());
            controller->Set_File_Buffer(vel_file,
                                        sizeof(VECTOR) * md_info->atom_numbers);
            Malloc_Safely((void**)&velocity_staging,
                          sizeof(VECTOR) * md_info->atom_numbers);
        }
    }

    start_frame = 0;
    if (controller->Command_Exist("rerun_start"))
    {
        start_frame = Rerun_Parse_Nonnegative_Int(
            controller, "rerun_start",
            "MD_INFORMATION::RERUN_information::Initial");
    }
    strip_frame = 0;
    if (controller->Command_Exist("rerun_strip"))
    {
        strip_frame = Rerun_Parse_Nonnegative_Int(
            controller, "rerun_strip",
            "MD_INFORMATION::RERUN_information::Initial");
    }
    need_box_update = 0;
    if (controller->Command_Exist("rerun_need_box_update"))
    {
        need_box_update =
            controller->Get_Bool("rerun_need_box_update",
                                 "MD_INFORMATION::RERUN_information::Initial");
    }
    md_info->sys.step_limit = INT_MAX;
    controller->printf("    End initializing rerun\n\n");
}

int MD_INFORMATION::RERUN_information::Read_Frame(VECTOR* box_length,
                                                  VECTOR* box_angle)
{
    const RerunBinaryReadResult coordinate = Rerun_Read_Binary_Frame(
        traj_file, coordinate_staging, sizeof(VECTOR), md_info->atom_numbers);
    VECTOR candidate_length;
    VECTOR candidate_angle;
    const int box_complete =
        Rerun_Read_Box_Record(box_file, &box_line_number, &candidate_length,
                              &candidate_angle, controller);

    RerunBinaryReadResult velocity = {
        RerunBinaryReadState::complete,
        static_cast<std::size_t>(md_info->atom_numbers)};
    if (has_velocity_stream)
    {
        velocity = Rerun_Read_Binary_Frame(
            vel_file, velocity_staging, sizeof(VECTOR), md_info->atom_numbers);
    }

    if (coordinate.state == RerunBinaryReadState::io_error ||
        velocity.state == RerunBinaryReadState::io_error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat,
            "MD_INFORMATION::RERUN_information::Read_Frame",
            "Reason:\n\tI/O error while reading rerun frame %llu "
            "(coordinate=%s, velocity=%s)\n",
            static_cast<unsigned long long>(complete_frames_read),
            Rerun_Binary_State_Name(coordinate.state),
            has_velocity_stream ? Rerun_Binary_State_Name(velocity.state)
                                : "not configured");
    }
    if (coordinate.state == RerunBinaryReadState::partial ||
        velocity.state == RerunBinaryReadState::partial)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat,
            "MD_INFORMATION::RERUN_information::Read_Frame",
            "Reason:\n\ttruncated binary rerun frame %llu: coordinate "
            "%zu/%d atoms, velocity %zu/%d atoms\n",
            static_cast<unsigned long long>(complete_frames_read),
            coordinate.count, md_info->atom_numbers,
            has_velocity_stream ? velocity.count : 0,
            has_velocity_stream ? md_info->atom_numbers : 0);
    }

    const bool coordinate_eof =
        coordinate.state == RerunBinaryReadState::clean_eof;
    const bool velocity_eof = !has_velocity_stream ||
                              velocity.state == RerunBinaryReadState::clean_eof;
    const bool all_eof = coordinate_eof && box_complete == 0 && velocity_eof;
    if (all_eof) return 0;

    const bool all_complete =
        coordinate.state == RerunBinaryReadState::complete &&
        box_complete == 1 &&
        (!has_velocity_stream ||
         velocity.state == RerunBinaryReadState::complete);
    if (!all_complete)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat,
            "MD_INFORMATION::RERUN_information::Read_Frame",
            "Reason:\n\trerun streams are out of sync at frame %llu: "
            "coordinate=%s, box=%s, velocity=%s\n",
            static_cast<unsigned long long>(complete_frames_read),
            Rerun_Binary_State_Name(coordinate.state),
            box_complete ? "complete" : "clean EOF",
            has_velocity_stream ? Rerun_Binary_State_Name(velocity.state)
                                : "not configured");
    }
    if (complete_frames_read == std::numeric_limits<std::uint64_t>::max())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow,
            "MD_INFORMATION::RERUN_information::Read_Frame",
            "Reason:\n\trerun frame counter overflow\n");
    }
    ++complete_frames_read;
    *box_length = candidate_length;
    *box_angle = candidate_angle;
    return 1;
}

void MD_INFORMATION::RERUN_information::Close_Files()
{
    std::string close_errors;
    auto close_one = [&](FILE** input, const char* description)
    {
        if (*input == NULL) return;
        errno = 0;
        if (fclose(*input) != 0)
        {
            if (!close_errors.empty()) close_errors += "; ";
            close_errors += description;
            close_errors += ": ";
            close_errors += strerror(errno);
        }
        *input = NULL;
    };
    close_one(&traj_file, "coordinate trajectory");
    close_one(&box_file, "box trajectory");
    close_one(&vel_file, "velocity trajectory");
    if (!close_errors.empty() && controller != NULL)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorBadFileFormat,
            "MD_INFORMATION::RERUN_information::Close_Files",
            "Reason:\n\tI/O error while closing rerun input: %s\n",
            close_errors.c_str());
    }
}

void MD_INFORMATION::RERUN_information::Clear()
{
    Close_Files();
    free(coordinate_staging);
    free(velocity_staging);
    coordinate_staging = NULL;
    velocity_staging = NULL;
    controller = NULL;
    md_info = NULL;
    has_velocity_stream = false;
    has_frame = false;
    frame_box_length = VECTOR();
    frame_box_angle = VECTOR();
    box_line_number = 0;
    complete_frames_read = 0;
    g = LTMatrix3();
}

bool MD_INFORMATION::RERUN_information::Iteration(int strip)
{
    if (strip < 0) strip = strip_frame;
    VECTOR candidate_box_length;
    VECTOR candidate_box_angle;
    int frame_available = 0;
    int candidate_steps = md_info->sys.steps;

    if (CONTROLLER::MPI_rank == 0)
    {
        for (int i = 0; i < strip; ++i)
        {
            VECTOR skipped_length;
            VECTOR skipped_angle;
            if (!Read_Frame(&skipped_length, &skipped_angle))
            {
                if (!has_frame)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorBadFileFormat,
                        "MD_INFORMATION::RERUN_information::Iteration",
                        "Reason:\n\trerun input ended before the first "
                        "selected frame after skipping %d frame(s)\n",
                        strip);
                }
                md_info->sys.step_limit = md_info->sys.steps;
                Close_Files();
                break;
            }
        }
        if (traj_file != NULL)
        {
            frame_available =
                Read_Frame(&candidate_box_length, &candidate_box_angle);
            if (!frame_available)
            {
                if (!has_frame)
                {
                    controller->Throw_SPONGE_Error(
                        spongeErrorBadFileFormat,
                        "MD_INFORMATION::RERUN_information::Iteration",
                        "Reason:\n\trerun input contains no complete selected "
                        "frame\n");
                }
                md_info->sys.step_limit = md_info->sys.steps;
                Close_Files();
            }
            else
            {
                if (strip > INT_MAX - md_info->sys.steps)
                {
                    controller->Throw_SPONGE_Error(
                        spongeErrorOverflow,
                        "MD_INFORMATION::RERUN_information::Iteration",
                        "Reason:\n\trerun step counter overflows while "
                        "skipping frames\n");
                }
                candidate_steps = md_info->sys.steps + strip;
            }
        }
    }

#ifdef USE_MPI
    MPI_Bcast(&frame_available, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&md_info->sys.step_limit, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if (!frame_available)
    {
        g = LTMatrix3();
        return false;
    }

#ifdef USE_MPI
    MPI_Bcast(&candidate_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&candidate_box_length, sizeof(VECTOR), MPI_BYTE, 0,
              MPI_COMM_WORLD);
    MPI_Bcast(&candidate_box_angle, sizeof(VECTOR), MPI_BYTE, 0,
              MPI_COMM_WORLD);
#endif

    const bool box_changed =
        !has_frame ||
        !Rerun_Box_Is_Exactly_Equal(frame_box_length, frame_box_angle,
                                    candidate_box_length, candidate_box_angle);
    LTMatrix3 candidate_g;
    if (box_changed && md_info->pbc.is_initialized)
    {
        if (!Float_Memory_Is_Normal(&md_info->dt) || !(md_info->dt > 0.0f))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MD_INFORMATION::RERUN_information::Iteration",
                "Reason:\n\trerun box deformation requires a finite positive "
                "normal internal timestep; got %.9g\n",
                md_info->dt);
        }
        const LTMatrix3 new_cell =
            md_info->pbc.Get_Cell(candidate_box_length, candidate_box_angle);
        const LTMatrix3 old_cell =
            md_info->pbc.pbc ? md_info->pbc.cell : md_info->pbc.reference_cell;
        candidate_g = (1.0f / md_info->dt) *
                      (new_cell * inv(old_cell) - LTMatrix3(1, 0, 1, 0, 0, 1));
        if (!Rerun_Deformation_Is_Representable(candidate_g))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "MD_INFORMATION::RERUN_information::Iteration",
                "Reason:\n\tthe rerun box deformation matrix is not "
                "representable as finite zero-or-normal floats\n");
        }
    }

    // Publication begins only after every input stream, counter update and
    // semantic box/deformation check has succeeded.  Until this point the
    // selected frame exists solely in staging buffers and local candidates.
    if (CONTROLLER::MPI_rank == 0)
    {
        md_info->sys.steps = candidate_steps;
        memcpy(md_info->coordinate, coordinate_staging,
               sizeof(VECTOR) * md_info->atom_numbers);
        if (has_velocity_stream)
        {
            memcpy(md_info->velocity, velocity_staging,
                   sizeof(VECTOR) * md_info->atom_numbers);
        }
    }
#ifdef USE_MPI
    MPI_Bcast(&md_info->sys.steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(md_info->coordinate, sizeof(VECTOR) * md_info->atom_numbers,
              MPI_BYTE, 0, MPI_COMM_WORLD);
    if (has_velocity_stream)
    {
        MPI_Bcast(md_info->velocity, sizeof(VECTOR) * md_info->atom_numbers,
                  MPI_BYTE, 0, MPI_COMM_WORLD);
    }
#endif

    g = candidate_g;
    frame_box_length = candidate_box_length;
    frame_box_angle = candidate_box_angle;
    md_info->sys.box_length = candidate_box_length;
    md_info->sys.box_angle = candidate_box_angle;
    has_frame = true;

    deviceMemcpy(md_info->crd, md_info->coordinate,
                 sizeof(VECTOR) * md_info->atom_numbers,
                 deviceMemcpyHostToDevice);
    if (has_velocity_stream)
    {
        deviceMemcpy(md_info->vel, md_info->velocity,
                     sizeof(VECTOR) * md_info->atom_numbers,
                     deviceMemcpyHostToDevice);
    }
    return box_changed;
}
