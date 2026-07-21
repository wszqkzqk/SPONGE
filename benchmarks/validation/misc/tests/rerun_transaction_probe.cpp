#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "MD_core/MD_core.h"
#include "MD_core/pbc.hpp"
#include "common.h"
#include "control.h"
#define BOX_TRAJ_COMMAND "box"
#define TRAJ_COMMAND "crd"
#define VEL_TRAJ_COMMAND "vel"
#include "MD_core/rerun.hpp"
#undef BOX_TRAJ_COMMAND
#undef TRAJ_COMMAND
#undef VEL_TRAJ_COMMAND

int CONTROLLER::MPI_rank = 0;
int CONTROLLER::MPI_size = 1;
CONTROLLER controller{};

bool CONTROLLER::Command_Exist(const char*) { return false; }

bool CONTROLLER::Command_Exist(const char*, const char*) { return false; }

const char* CONTROLLER::Command(const char*) { return ""; }

const char* CONTROLLER::Command(const char*, const char*) { return ""; }

const char* CONTROLLER::Original_Command(const char*) { return ""; }

const char* CONTROLLER::Original_Command(const char*, const char*)
{
    return "";
}

void deviceMemcpy(void* destination, const void* source, std::size_t size,
                  deviceMemcpyKind)
{
    std::memcpy(destination, source, size);
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

static bool Same_Vector(const VECTOR& lhs, const VECTOR& rhs)
{
    return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
}

static bool Same_Matrix(const LTMatrix3& lhs, const LTMatrix3& rhs)
{
    return lhs.a11 == rhs.a11 && lhs.a21 == rhs.a21 && lhs.a22 == rhs.a22 &&
           lhs.a31 == rhs.a31 && lhs.a32 == rhs.a32 && lhs.a33 == rhs.a33;
}

struct Transaction_Check
{
    FILE* result;
    MD_INFORMATION* md;
    MD_INFORMATION::RERUN_information* rerun;
    VECTOR* host_coordinate;
    VECTOR* host_velocity;
    VECTOR* device_coordinate;
    VECTOR* device_velocity;
};

static void Record_Transaction_State(void* opaque) noexcept
{
    Transaction_Check* check = static_cast<Transaction_Check*>(opaque);
    const bool state_preserved =
        check->md->sys.steps == 7 &&
        Same_Vector(check->md->sys.box_length, VECTOR(1.0f, 1.0f, 1.0f)) &&
        Same_Vector(check->md->sys.box_angle, VECTOR(90.0f, 90.0f, 90.0f)) &&
        Same_Vector(*check->host_coordinate, VECTOR(1.0f, 2.0f, 3.0f)) &&
        Same_Vector(*check->host_velocity, VECTOR(4.0f, 5.0f, 6.0f)) &&
        Same_Vector(*check->device_coordinate, VECTOR(1.0f, 2.0f, 3.0f)) &&
        Same_Vector(*check->device_velocity, VECTOR(4.0f, 5.0f, 6.0f)) &&
        check->rerun->has_frame &&
        Same_Vector(check->rerun->frame_box_length, VECTOR(1.0f, 1.0f, 1.0f)) &&
        Same_Vector(check->rerun->frame_box_angle,
                    VECTOR(90.0f, 90.0f, 90.0f)) &&
        Same_Matrix(check->rerun->g, LTMatrix3(2.0f));
    if (check->result != nullptr)
    {
        std::fputc(state_preserved ? 'P' : 'M', check->result);
        std::fclose(check->result);
        check->result = nullptr;
    }
}

int main(int argc, char** argv)
{
    if (argc != 2) return EXIT_FAILURE;
    CONTROLLER controller{};
    MD_INFORMATION md{};
    md.atom_numbers = 1;
    md.dt = 1.0e-6f;
    md.sys.steps = 7;
    md.sys.box_length = VECTOR(1.0f, 1.0f, 1.0f);
    md.sys.box_angle = VECTOR(90.0f, 90.0f, 90.0f);

    VECTOR host_coordinate(1.0f, 2.0f, 3.0f);
    VECTOR host_velocity(4.0f, 5.0f, 6.0f);
    VECTOR device_coordinate = host_coordinate;
    VECTOR device_velocity = host_velocity;
    md.coordinate = &host_coordinate;
    md.velocity = &host_velocity;
    md.crd = &device_coordinate;
    md.vel = &device_velocity;

    md.pbc.controller = &controller;
    md.pbc.md_info = &md;
    md.pbc.is_initialized = true;
    md.pbc.pbc = true;
    const LTMatrix3 unit_cell(1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f);
    md.pbc.cell = unit_cell;
    md.pbc.reference_cell = unit_cell;
    md.pbc.rcell = unit_cell;

    auto& rerun = md.rerun;
    rerun.controller = &controller;
    rerun.md_info = &md;
    rerun.has_frame = true;
    rerun.frame_box_length = md.sys.box_length;
    rerun.frame_box_angle = md.sys.box_angle;
    rerun.g = LTMatrix3(2.0f);
    VECTOR coordinate_staging;
    rerun.coordinate_staging = &coordinate_staging;

    rerun.traj_file = std::tmpfile();
    rerun.box_file = std::tmpfile();
    if (rerun.traj_file == nullptr || rerun.box_file == nullptr)
    {
        return EXIT_FAILURE;
    }
    const VECTOR candidate_coordinate(9.0f, 8.0f, 7.0f);
    if (std::fwrite(&candidate_coordinate, sizeof(VECTOR), 1,
                    rerun.traj_file) != 1 ||
        std::fputs("1e37 1 1 90 90 90\n", rerun.box_file) < 0)
    {
        return EXIT_FAILURE;
    }
    std::rewind(rerun.traj_file);
    std::rewind(rerun.box_file);

    FILE* result = std::fopen(argv[1], "wb");
    if (result == nullptr) return EXIT_FAILURE;
    Transaction_Check check{result,          &md,
                            &rerun,          &host_coordinate,
                            &host_velocity,  &device_coordinate,
                            &device_velocity};
    if (!controller.Register_Fatal_Cleanup(Record_Transaction_State, &check))
    {
        std::fclose(result);
        return EXIT_FAILURE;
    }
    (void)rerun.Iteration(0);
    controller.Unregister_Fatal_Cleanup(Record_Transaction_State, &check);
    if (check.result != nullptr) std::fclose(check.result);
    return EXIT_FAILURE;
}
