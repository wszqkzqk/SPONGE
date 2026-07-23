#include "full_neighbor_list.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
void Reset_Active_Prefix(FULL_NEIGHBOR_LIST* list, int next_active)
{
    const int clear_count =
        std::max(list->active_owned_atom_numbers, next_active);
    if (list->h_nl != NULL && list->d_nl != NULL && clear_count > 0)
    {
        for (int i = 0; i < clear_count; ++i)
        {
            list->h_nl[i].atom_numbers = 0;
            list->h_nl[i].ghost_numbers = 0;
        }
        deviceMemcpy(list->d_nl, list->h_nl,
                     sizeof(ATOM_GROUP) * (size_t)clear_count,
                     deviceMemcpyHostToDevice);
    }
    list->active_owned_atom_numbers = next_active;
}

bool Fail_Build(FULL_NEIGHBOR_LIST* list, int error, int atom = -1,
                int value = -1, int required_neighbor_capacity = 0)
{
    Reset_Active_Prefix(list, 0);
    list->last_build_error = error;
    list->last_error_atom = atom;
    list->last_error_value = value;
    list->last_required_neighbor_capacity = required_neighbor_capacity;
    if (required_neighbor_capacity == 0 && list->d_overflow != NULL)
        deviceMemset(list->d_overflow, 0, sizeof(int));
    return false;
}

bool Prepare_Build(FULL_NEIGHBOR_LIST* list, const ATOM_GROUP* half_nl,
                   int owned_atom_numbers, int coordinate_numbers)
{
    if (!list->is_initialized)
        return Fail_Build(list, FULL_NEIGHBOR_LIST::BUILD_NOT_INITIALIZED);
    if (owned_atom_numbers < 0 || coordinate_numbers < owned_atom_numbers ||
        owned_atom_numbers > list->atom_capacity ||
        (owned_atom_numbers > 0 && half_nl == NULL))
        return Fail_Build(list, FULL_NEIGHBOR_LIST::BUILD_INVALID_ARGUMENT);

    list->last_build_error = FULL_NEIGHBOR_LIST::BUILD_OK;
    list->last_error_atom = -1;
    list->last_error_value = -1;
    list->last_required_neighbor_capacity = 0;
    deviceMemset(list->d_overflow, 0, sizeof(int));
    deviceMemset(list->d_build_error, 0, 3 * sizeof(int));
    Reset_Active_Prefix(list, owned_atom_numbers);
    return true;
}

static __device__ __forceinline__ int Claim_Build_Error(int* build_error,
                                                        int error)
{
#ifdef GPU_ARCH_NAME
    return atomicCAS(build_error, FULL_NEIGHBOR_LIST::BUILD_OK, error);
#else
    int observed = FULL_NEIGHBOR_LIST::BUILD_OK;
#pragma omp critical(sponge_full_neighbor_list_build_error)
    {
        observed = build_error[0];
        if (observed == FULL_NEIGHBOR_LIST::BUILD_OK) build_error[0] = error;
    }
    return observed;
#endif
}

static __device__ __forceinline__ void Record_Build_Error(int* build_error,
                                                          int error, int atom,
                                                          int value)
{
    if (Claim_Build_Error(build_error, error) == FULL_NEIGHBOR_LIST::BUILD_OK)
    {
        build_error[1] = atom;
        build_error[2] = value;
    }
}

static __device__ __forceinline__ void Record_Required_Capacity(
    int* required_capacity, int required)
{
    if (required <= 0) return;
#ifdef GPU_ARCH_NAME
    atomicMax(required_capacity, required);
#else
#pragma omp critical(sponge_full_neighbor_list_required_capacity)
    {
        if (required_capacity[0] < required) required_capacity[0] = required;
    }
#endif
}

static __device__ __forceinline__ int Reserve_One_Count(int* counter,
                                                        int* required_capacity)
{
#ifdef GPU_ARCH_NAME
    int observed = atomicAdd(counter, 0);
    while (true)
    {
        // INT_MAX: std::numeric_limits<int>::max() is a host-only constexpr
        // under nvcc and cannot be called from this device function.
        if (observed == INT_MAX)
        {
            Record_Required_Capacity(required_capacity, observed);
            return observed;
        }
        const int prior = atomicCAS(counter, observed, observed + 1);
        if (prior == observed) return observed;
        observed = prior;
    }
#else
    int slot = 0;
#pragma omp critical(sponge_full_neighbor_list_counter_reservation)
    {
        slot = counter[0];
        if (slot != std::numeric_limits<int>::max()) counter[0] = slot + 1;
    }
    if (slot == std::numeric_limits<int>::max())
        Record_Required_Capacity(required_capacity, slot);
    return slot;
#endif
}

static __device__ __forceinline__ void Append_Full_Neighbor(
    ATOM_GROUP* full_nl, int atom_i, int atom_j, bool neighbor_is_ghost,
    int max_neighbor_numbers, int* overflow_flag)
{
    const int slot =
        Reserve_One_Count(&full_nl[atom_i].atom_numbers, overflow_flag);
    if (neighbor_is_ghost)
        Reserve_One_Count(&full_nl[atom_i].ghost_numbers, overflow_flag);
    if (slot < max_neighbor_numbers)
        full_nl[atom_i].atom_serial[slot] = atom_j;
    else
        // INT_MAX: see Reserve_One_Count above.
        Record_Required_Capacity(overflow_flag,
                                 slot == INT_MAX ? slot : slot + 1);
}

static __device__ __forceinline__ bool Validate_Half_Group(
    const ATOM_GROUP& group, int atom_i, int max_neighbor_numbers,
    int* build_error)
{
    if (group.atom_numbers < 0 || group.atom_numbers > max_neighbor_numbers ||
        group.ghost_numbers < 0 || group.ghost_numbers > group.atom_numbers)
    {
        const int invalid_value =
            group.atom_numbers < 0 || group.atom_numbers > max_neighbor_numbers
                ? group.atom_numbers
                : group.ghost_numbers;
        Record_Build_Error(build_error,
                           FULL_NEIGHBOR_LIST::BUILD_INVALID_HALF_COUNT, atom_i,
                           invalid_value);
        return false;
    }
    if (group.atom_numbers > 0 && group.atom_serial == NULL)
    {
        Record_Build_Error(build_error,
                           FULL_NEIGHBOR_LIST::BUILD_INVALID_HALF_POINTER,
                           atom_i, group.atom_numbers);
        return false;
    }
    return true;
}

bool Matrix_Is_Finite(const LTMatrix3& matrix)
{
    return Float_Memory_Is_Finite(&matrix.a11) &&
           Float_Memory_Is_Finite(&matrix.a21) &&
           Float_Memory_Is_Finite(&matrix.a22) &&
           Float_Memory_Is_Finite(&matrix.a31) &&
           Float_Memory_Is_Finite(&matrix.a32) &&
           Float_Memory_Is_Finite(&matrix.a33);
}
}  // namespace

bool FULL_NEIGHBOR_LIST::Initial(int requested_atom_capacity,
                                 int requested_max_neighbor_numbers)
{
    if (is_initialized)
    {
        if (requested_atom_capacity == atom_capacity &&
            requested_max_neighbor_numbers == max_neighbor_numbers)
            return true;
        last_build_error = BUILD_INVALID_ARGUMENT;
        last_error_atom = -1;
        last_error_value = -1;
        return false;
    }
    if (requested_atom_capacity <= 0 || requested_max_neighbor_numbers <= 0)
    {
        last_build_error = BUILD_INVALID_ARGUMENT;
        last_error_atom = -1;
        last_error_value = -1;
        return false;
    }
    const size_t capacity = (size_t)requested_atom_capacity;
    const size_t neighbor_capacity = (size_t)requested_max_neighbor_numbers;
    if (capacity > std::numeric_limits<size_t>::max() / sizeof(ATOM_GROUP) ||
        neighbor_capacity > std::numeric_limits<size_t>::max() / capacity ||
        capacity * neighbor_capacity >
            std::numeric_limits<size_t>::max() / sizeof(int))
    {
        last_build_error = BUILD_INVALID_ARGUMENT;
        last_error_atom = -1;
        last_error_value = -1;
        return false;
    }

    atom_capacity = requested_atom_capacity;
    max_neighbor_numbers = requested_max_neighbor_numbers;
    const size_t neighbor_slots = capacity * neighbor_capacity;
    if (!Malloc_Safely((void**)&h_nl, sizeof(ATOM_GROUP) * capacity) ||
        !Device_Malloc_Safely((void**)&d_temp, sizeof(int) * neighbor_slots))
    {
        Clear();
        last_build_error = BUILD_ALLOCATION_FAILED;
        return false;
    }
    for (int i = 0; i < atom_capacity; ++i)
    {
        h_nl[i].atom_numbers = 0;
        h_nl[i].ghost_numbers = 0;
        h_nl[i].atom_serial = d_temp + (size_t)max_neighbor_numbers * i;
    }
    if (!Device_Malloc_And_Copy_Safely((void**)&d_nl, h_nl,
                                       sizeof(ATOM_GROUP) * capacity) ||
        !Device_Malloc_Safely((void**)&d_overflow, sizeof(int)) ||
        !Device_Malloc_Safely((void**)&d_build_error, 3 * sizeof(int)))
    {
        Clear();
        last_build_error = BUILD_ALLOCATION_FAILED;
        return false;
    }
    deviceMemset(d_overflow, 0, sizeof(int));
    deviceMemset(d_build_error, 0, 3 * sizeof(int));
    active_owned_atom_numbers = 0;
    last_build_error = BUILD_OK;
    last_required_neighbor_capacity = 0;
    is_initialized = true;
    return true;
}

static __global__ void Build_Full_Neighbor_List_Kernel(
    const ATOM_GROUP* half_nl, ATOM_GROUP* full_nl, int owned_atom_numbers,
    int coordinate_numbers, int max_neighbor_numbers, int* overflow_flag,
    int* build_error)
{
    SIMPLE_DEVICE_FOR(atom_i, owned_atom_numbers)
    {
        const ATOM_GROUP half_group = half_nl[atom_i];
        // Guard the loop instead of `continue`, which is illegal inside
        // CUDA's SIMPLE_DEVICE_FOR (an `if` guard, not a loop).
        const bool half_group_valid = Validate_Half_Group(
            half_group, atom_i, max_neighbor_numbers, build_error);
        for (int k = 0; half_group_valid && k < half_group.atom_numbers; ++k)
        {
            const int atom_j = half_group.atom_serial[k];
            if (atom_j < 0 || atom_j >= coordinate_numbers)
            {
                Record_Build_Error(
                    build_error,
                    FULL_NEIGHBOR_LIST::BUILD_INVALID_NEIGHBOR_INDEX, atom_i,
                    atom_j);
                continue;
            }
            Append_Full_Neighbor(full_nl, atom_i, atom_j,
                                 atom_j >= owned_atom_numbers,
                                 max_neighbor_numbers, overflow_flag);
            if (atom_j < owned_atom_numbers)
                Append_Full_Neighbor(full_nl, atom_j, atom_i, false,
                                     max_neighbor_numbers, overflow_flag);
        }
    }
}

bool FULL_NEIGHBOR_LIST::Build_From_Half(const ATOM_GROUP* half_nl,
                                         int owned_atom_numbers,
                                         int coordinate_numbers)
{
    if (!Prepare_Build(this, half_nl, owned_atom_numbers, coordinate_numbers))
        return false;
    if (owned_atom_numbers == 0) return true;
    Launch_Device_Kernel(
        Build_Full_Neighbor_List_Kernel,
        (owned_atom_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, half_nl, d_nl,
        owned_atom_numbers, coordinate_numbers, max_neighbor_numbers,
        d_overflow, d_build_error);
    int error[3] = {0, -1, -1};
    deviceMemcpy(error, d_build_error, sizeof(error), deviceMemcpyDeviceToHost);
    if (error[0] != BUILD_OK)
        return Fail_Build(this, error[0], error[1], error[2]);
    int required_neighbor_capacity = 0;
    deviceMemcpy(&required_neighbor_capacity, d_overflow, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (required_neighbor_capacity > 0)
        return Fail_Build(this, BUILD_CAPACITY_EXCEEDED, -1,
                          required_neighbor_capacity,
                          required_neighbor_capacity);
    return true;
}

static __global__ void Validate_Full_Neighbor_Coordinates(
    const VECTOR* crd, int coordinate_numbers, int* build_error)
{
    SIMPLE_DEVICE_FOR(atom_i, coordinate_numbers)
    {
        const VECTOR coordinate = crd[atom_i];
        if (!isfinite(coordinate.x) || !isfinite(coordinate.y) ||
            !isfinite(coordinate.z))
        {
            Record_Build_Error(build_error,
                               FULL_NEIGHBOR_LIST::BUILD_INVALID_GEOMETRY,
                               atom_i, atom_i);
        }
    }
}

static __global__ void Build_Full_Neighbor_List_With_Cutoff_Kernel(
    const ATOM_GROUP* half_nl, ATOM_GROUP* full_nl, int owned_atom_numbers,
    int coordinate_numbers, int max_neighbor_numbers, int* overflow_flag,
    int* build_error, const VECTOR* crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, float cutoff_squared)
{
    SIMPLE_DEVICE_FOR(atom_i, owned_atom_numbers)
    {
        const ATOM_GROUP half_group = half_nl[atom_i];
        // Guard the loop instead of `continue`, which is illegal inside
        // CUDA's SIMPLE_DEVICE_FOR (an `if` guard, not a loop).
        const bool half_group_valid = Validate_Half_Group(
            half_group, atom_i, max_neighbor_numbers, build_error);
        const VECTOR ri = crd[atom_i];
        for (int k = 0; half_group_valid && k < half_group.atom_numbers; ++k)
        {
            const int atom_j = half_group.atom_serial[k];
            if (atom_j < 0 || atom_j >= coordinate_numbers)
            {
                Record_Build_Error(
                    build_error,
                    FULL_NEIGHBOR_LIST::BUILD_INVALID_NEIGHBOR_INDEX, atom_i,
                    atom_j);
                continue;
            }
            const VECTOR dr =
                Get_Periodic_Displacement(ri, crd[atom_j], cell, rcell);
            const float distance_squared = dr * dr;
            if (!isfinite(distance_squared))
            {
                Record_Build_Error(build_error,
                                   FULL_NEIGHBOR_LIST::BUILD_INVALID_GEOMETRY,
                                   atom_i, atom_j);
                continue;
            }
            if (!(distance_squared <= cutoff_squared)) continue;
            Append_Full_Neighbor(full_nl, atom_i, atom_j,
                                 atom_j >= owned_atom_numbers,
                                 max_neighbor_numbers, overflow_flag);
            if (atom_j < owned_atom_numbers)
                Append_Full_Neighbor(full_nl, atom_j, atom_i, false,
                                     max_neighbor_numbers, overflow_flag);
        }
    }
}

bool FULL_NEIGHBOR_LIST::Build_From_Half_With_Cutoff(
    const ATOM_GROUP* half_nl, int owned_atom_numbers, int coordinate_numbers,
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell,
    float cutoff)
{
    if (!Prepare_Build(this, half_nl, owned_atom_numbers, coordinate_numbers))
        return false;
    const float cutoff_squared = cutoff * cutoff;
    if ((coordinate_numbers > 0 && crd == NULL) || !(cutoff > 0.0f) ||
        !Float_Memory_Is_Normal(&cutoff) ||
        !Float_Memory_Is_Normal(&cutoff_squared))
        return Fail_Build(this, BUILD_INVALID_ARGUMENT);
    if (!Matrix_Is_Finite(cell) || !Matrix_Is_Finite(rcell))
        return Fail_Build(this, BUILD_INVALID_GEOMETRY);
    if (coordinate_numbers > 0)
    {
        Launch_Device_Kernel(
            Validate_Full_Neighbor_Coordinates,
            (coordinate_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, crd, coordinate_numbers,
            d_build_error);
        int coordinate_error[3] = {0, -1, -1};
        deviceMemcpy(coordinate_error, d_build_error, sizeof(coordinate_error),
                     deviceMemcpyDeviceToHost);
        if (coordinate_error[0] != BUILD_OK)
            return Fail_Build(this, coordinate_error[0], coordinate_error[1],
                              coordinate_error[2]);
    }
    if (owned_atom_numbers == 0) return true;
    Launch_Device_Kernel(
        Build_Full_Neighbor_List_With_Cutoff_Kernel,
        (owned_atom_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, half_nl, d_nl,
        owned_atom_numbers, coordinate_numbers, max_neighbor_numbers,
        d_overflow, d_build_error, crd, cell, rcell, cutoff_squared);
    int error[3] = {0, -1, -1};
    deviceMemcpy(error, d_build_error, sizeof(error), deviceMemcpyDeviceToHost);
    if (error[0] != BUILD_OK)
        return Fail_Build(this, error[0], error[1], error[2]);
    int required_neighbor_capacity = 0;
    deviceMemcpy(&required_neighbor_capacity, d_overflow, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (required_neighbor_capacity > 0)
        return Fail_Build(this, BUILD_CAPACITY_EXCEEDED, -1,
                          required_neighbor_capacity,
                          required_neighbor_capacity);
    return true;
}

const char* FULL_NEIGHBOR_LIST::Last_Error_Message() const
{
    switch (last_build_error)
    {
        case BUILD_OK:
            return "no error";
        case BUILD_NOT_INITIALIZED:
            return "full neighbor list is not initialized";
        case BUILD_INVALID_ARGUMENT:
            return "invalid full-neighbor-list count, pointer, or cutoff";
        case BUILD_INVALID_HALF_COUNT:
            return "invalid half-neighbor count";
        case BUILD_INVALID_HALF_POINTER:
            return "null half-neighbor storage for a nonempty row";
        case BUILD_INVALID_NEIGHBOR_INDEX:
            return "half-neighbor index is outside local-plus-ghost storage";
        case BUILD_ALLOCATION_FAILED:
            return "failed to allocate full-neighbor-list storage";
        case BUILD_CAPACITY_EXCEEDED:
            return "full-neighbor-list row capacity was exceeded";
        case BUILD_INVALID_GEOMETRY:
            return "non-finite coordinate or periodic geometry";
        default:
            return "unknown full-neighbor-list build error";
    }
}

void FULL_NEIGHBOR_LIST::Invalidate_Active()
{
    Reset_Active_Prefix(this, 0);
    last_build_error = BUILD_OK;
    last_error_atom = -1;
    last_error_value = -1;
    last_required_neighbor_capacity = 0;
    if (d_overflow != NULL) deviceMemset(d_overflow, 0, sizeof(int));
    if (d_build_error != NULL) deviceMemset(d_build_error, 0, 3 * sizeof(int));
}

void FULL_NEIGHBOR_LIST::Clear()
{
    if (d_temp != NULL) Free_Single_Device_Pointer((void**)&d_temp);
    if (h_nl != NULL || d_nl != NULL)
        Free_Host_And_Device_Pointer((void**)&h_nl, (void**)&d_nl);
    if (d_overflow != NULL) Free_Single_Device_Pointer((void**)&d_overflow);
    if (d_build_error != NULL)
        Free_Single_Device_Pointer((void**)&d_build_error);
    is_initialized = false;
    atom_capacity = 0;
    active_owned_atom_numbers = 0;
    max_neighbor_numbers = 0;
    last_build_error = BUILD_OK;
    last_error_atom = -1;
    last_error_value = -1;
    last_required_neighbor_capacity = 0;
}
