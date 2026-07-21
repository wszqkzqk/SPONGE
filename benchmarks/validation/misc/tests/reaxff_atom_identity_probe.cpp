#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "manybody/reaxff/reaxff.h"

unsigned int CONTROLLER::device_max_thread = 8;
int CONTROLLER::MPI_rank = 0;
int CONTROLLER::PP_MPI_size = 1;

void deviceMemcpy(void* destination, const void* source, size_t size,
                  deviceMemcpyKind)
{
    std::memmove(destination, source, size);
}

namespace
{

bool Close(float lhs, float rhs) { return std::fabs(lhs - rhs) <= 1.0e-6f; }

bool Check_Local_View(const char* name, const int* global_values,
                      const int* local_values, const int* local_to_global,
                      int atom_numbers)
{
    for (int local_id = 0; local_id < atom_numbers; local_id++)
    {
        if (local_values[local_id] != global_values[local_to_global[local_id]])
        {
            std::fprintf(stderr, "%s gather failed at local atom %d\n", name,
                         local_id);
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc == 2 && std::strcmp(argv[1], "--multi-pp") == 0)
    {
        REAXFF reaxff;
        CONTROLLER controller;
        reaxff.is_initialized = 1;
        CONTROLLER::PP_MPI_size = 2;
        reaxff.Validate_Parallel_Layout(&controller);
        std::fprintf(stderr, "multi-PP validation unexpectedly returned\n");
        return EXIT_FAILURE;
    }

    constexpr int atom_numbers = 5;
    const int local_to_global[atom_numbers] = {3, 0, 4, 1, 2};
    int global_types[8][atom_numbers] = {};
    int local_types[8][atom_numbers] = {};
    for (int module = 0; module < 8; module++)
        for (int global_id = 0; global_id < atom_numbers; global_id++)
            global_types[module][global_id] = 100 * (module + 1) + global_id;
    int global_is_hydrogen[atom_numbers] = {0, 1, 0, 1, 1};
    int local_is_hydrogen[atom_numbers] = {};

    REAXFF reaxff;
    CONTROLLER controller;
    reaxff.is_initialized = 1;
    reaxff.atom_numbers = atom_numbers;
    reaxff.eeq.is_initialized = 1;
    reaxff.bond_order.is_initialized = 1;
    reaxff.bond.is_initialized = 1;
    reaxff.vdw.is_initialized = 1;
    reaxff.ovun.is_initialized = 1;
    reaxff.angle.is_initialized = 1;
    reaxff.torsion.is_initialized = 1;
    reaxff.hb.is_initialized = 1;
    reaxff.eeq.d_atom_type_global = global_types[0];
    reaxff.eeq.d_atom_type = local_types[0];
    reaxff.bond_order.d_atom_type_global = global_types[1];
    reaxff.bond_order.d_atom_type = local_types[1];
    reaxff.bond.d_atom_type_global = global_types[2];
    reaxff.bond.d_atom_type = local_types[2];
    reaxff.vdw.d_atom_type_global = global_types[3];
    reaxff.vdw.d_atom_type = local_types[3];
    reaxff.ovun.d_atom_type_global = global_types[4];
    reaxff.ovun.d_atom_type = local_types[4];
    reaxff.angle.d_atom_type_global = global_types[5];
    reaxff.angle.d_atom_type = local_types[5];
    reaxff.torsion.d_atom_type_global = global_types[6];
    reaxff.torsion.d_atom_type = local_types[6];
    reaxff.hb.d_atom_type_global = global_types[7];
    reaxff.hb.d_atom_type = local_types[7];
    reaxff.hb.d_is_hydrogen_global = global_is_hydrogen;
    reaxff.hb.d_is_hydrogen = local_is_hydrogen;

    reaxff.Get_Local(&controller, local_to_global, atom_numbers, 0);
    const char* module_names[8] = {"EEQ",  "bond order", "bond",    "vdW",
                                   "ovun", "angle",      "torsion", "HB"};
    for (int module = 0; module < 8; module++)
    {
        if (!Check_Local_View(module_names[module], global_types[module],
                              local_types[module], local_to_global,
                              atom_numbers))
            return EXIT_FAILURE;
    }
    if (!Check_Local_View("HB hydrogen flag", global_is_hydrogen,
                          local_is_hydrogen, local_to_global, atom_numbers))
        return EXIT_FAILURE;

    const float local_charge[atom_numbers] = {0.4f, 0.1f, 0.5f, 0.2f, 0.3f};
    float global_charge[atom_numbers] = {};
    ReaxFFAtomIdentity::Scatter_Float_By_Global_Id_Kernel(
        atom_numbers, local_to_global, local_charge, global_charge);
    for (int global_id = 0; global_id < atom_numbers; global_id++)
    {
        const float expected = 0.1f * static_cast<float>(global_id + 1);
        if (!Close(global_charge[global_id], expected))
        {
            std::fprintf(stderr, "charge scatter failed at global atom %d\n",
                         global_id);
            return EXIT_FAILURE;
        }
    }

    constexpr int history_count = 3;
    float global_history[history_count * atom_numbers] = {};
    for (int frame = 0; frame < history_count; frame++)
    {
        for (int global_id = 0; global_id < atom_numbers; global_id++)
        {
            global_history[frame * atom_numbers + global_id] =
                100.0f * static_cast<float>(frame) +
                static_cast<float>(global_id);
        }
    }
    const float coefficients[5] = {1.0f, -2.0f, 0.5f, 0.0f, 0.0f};
    float local_prediction[atom_numbers] = {};
    ReaxFFAtomIdentity::Gather_Float_History_By_Global_Id_Kernel(
        atom_numbers, local_to_global, local_prediction, global_history,
        atom_numbers, history_count, coefficients[0], coefficients[1],
        coefficients[2], coefficients[3], coefficients[4]);
    for (int local_id = 0; local_id < atom_numbers; local_id++)
    {
        const int global_id = local_to_global[local_id];
        const float expected =
            global_history[global_id] -
            2.0f * global_history[atom_numbers + global_id] +
            0.5f * global_history[2 * atom_numbers + global_id];
        if (!Close(local_prediction[local_id], expected))
        {
            std::fprintf(stderr,
                         "history gather failed at local/global atom %d/%d\n",
                         local_id, global_id);
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
