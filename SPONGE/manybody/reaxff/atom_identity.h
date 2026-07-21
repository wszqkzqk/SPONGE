#ifndef REAXFF_ATOM_IDENTITY_H
#define REAXFF_ATOM_IDENTITY_H

#include "../../control.h"

namespace ReaxFFAtomIdentity
{

static __global__ void Gather_Int_By_Global_Id_Kernel(
    int local_atom_numbers, const int* local_to_global,
    const int* global_values, int* local_values)
{
    SIMPLE_DEVICE_FOR(local_id, local_atom_numbers)
    {
        local_values[local_id] = global_values[local_to_global[local_id]];
    }
}

static __global__ void Scatter_Float_By_Global_Id_Kernel(
    int local_atom_numbers, const int* local_to_global,
    const float* local_values, float* global_values)
{
    SIMPLE_DEVICE_FOR(local_id, local_atom_numbers)
    {
        global_values[local_to_global[local_id]] = local_values[local_id];
    }
}

static __global__ void Gather_Float_History_By_Global_Id_Kernel(
    int local_atom_numbers, const int* local_to_global, float* local_values,
    const float* global_history, int global_stride, int history_count, float c0,
    float c1, float c2, float c3, float c4)
{
    SIMPLE_DEVICE_FOR(local_id, local_atom_numbers)
    {
        const int global_id = local_to_global[local_id];
        float value = c0 * global_history[global_id];
        if (history_count >= 2)
            value += c1 * global_history[global_stride + global_id];
        if (history_count >= 3)
            value += c2 * global_history[2 * global_stride + global_id];
        if (history_count >= 4)
            value += c3 * global_history[3 * global_stride + global_id];
        if (history_count >= 5)
            value += c4 * global_history[4 * global_stride + global_id];
        local_values[local_id] = value;
    }
}

inline void Gather_Int_By_Global_Id(int local_atom_numbers,
                                    const int* local_to_global,
                                    const int* global_values, int* local_values)
{
    if (local_atom_numbers <= 0) return;
    const dim3 block_size = {CONTROLLER::device_max_thread};
    const dim3 grid_size = {(local_atom_numbers + block_size.x - 1) /
                            block_size.x};
    Launch_Device_Kernel(Gather_Int_By_Global_Id_Kernel, grid_size, block_size,
                         0, NULL, local_atom_numbers, local_to_global,
                         global_values, local_values);
}

inline void Scatter_Float_By_Global_Id(int local_atom_numbers,
                                       const int* local_to_global,
                                       const float* local_values,
                                       float* global_values)
{
    if (local_atom_numbers <= 0) return;
    const dim3 block_size = {CONTROLLER::device_max_thread};
    const dim3 grid_size = {(local_atom_numbers + block_size.x - 1) /
                            block_size.x};
    Launch_Device_Kernel(Scatter_Float_By_Global_Id_Kernel, grid_size,
                         block_size, 0, NULL, local_atom_numbers,
                         local_to_global, local_values, global_values);
}

inline void Gather_Float_History_By_Global_Id(
    int local_atom_numbers, const int* local_to_global, float* local_values,
    const float* global_history, int global_stride, int history_count,
    const float* coefficients)
{
    if (local_atom_numbers <= 0) return;
    const dim3 block_size = {CONTROLLER::device_max_thread};
    const dim3 grid_size = {(local_atom_numbers + block_size.x - 1) /
                            block_size.x};
    Launch_Device_Kernel(
        Gather_Float_History_By_Global_Id_Kernel, grid_size, block_size, 0,
        NULL, local_atom_numbers, local_to_global, local_values, global_history,
        global_stride, history_count, coefficients[0], coefficients[1],
        coefficients[2], coefficients[3], coefficients[4]);
}

}  // namespace ReaxFFAtomIdentity

#endif
