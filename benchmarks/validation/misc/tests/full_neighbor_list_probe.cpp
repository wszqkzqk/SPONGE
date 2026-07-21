#define NO_GLOBAL_CONTROLLER

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

#include "neighbor_list/full_neighbor_list.h"

unsigned int CONTROLLER::device_max_thread = 64;

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

int atomicAdd(int* address, int value)
{
    int previous;
#pragma omp atomic capture
    {
        previous = *address;
        *address += value;
    }
    return previous;
}

int atomicExch(int* address, int value)
{
    int previous;
#pragma omp atomic capture
    {
        previous = *address;
        *address = value;
    }
    return previous;
}

void deviceMemcpy(void* destination, const void* source, std::size_t size,
                  deviceMemcpyKind)
{
    if (destination != source) std::memcpy(destination, source, size);
}

void deviceMemset(void* destination, int value, std::size_t size)
{
    int* values = static_cast<int*>(destination);
    for (std::size_t i = 0; i < size / sizeof(int); ++i) values[i] = value;
}

void deviceFree(void* pointer) { std::free(pointer); }

void Free_Single_Device_Pointer(void** pointer)
{
    std::free(pointer[0]);
    pointer[0] = NULL;
}

void Free_Host_And_Device_Pointer(void** host_pointer, void** device_pointer)
{
    if (host_pointer != NULL)
    {
        std::free(host_pointer[0]);
        host_pointer[0] = NULL;
    }
    device_pointer[0] = NULL;
}

#include "neighbor_list/full_neighbor_list.cpp"

namespace
{
struct HALF_LIST
{
    std::vector<ATOM_GROUP> rows;
    std::vector<std::vector<int>> storage;

    explicit HALF_LIST(int row_count) : rows(row_count), storage(row_count) {}

    void Set(int row, std::initializer_list<int> neighbors, int ghost_count)
    {
        storage[row].assign(neighbors);
        rows[row].atom_numbers = static_cast<int>(storage[row].size());
        rows[row].ghost_numbers = ghost_count;
        rows[row].atom_serial =
            storage[row].empty() ? NULL : storage[row].data();
    }
};

bool Expect_Row(const FULL_NEIGHBOR_LIST& list, int row,
                std::initializer_list<int> expected, int expected_ghosts)
{
    const ATOM_GROUP& group = list.h_nl[row];
    if (group.atom_numbers != static_cast<int>(expected.size()) ||
        group.ghost_numbers != expected_ghosts)
        return false;
    std::vector<int> actual(group.atom_serial,
                            group.atom_serial + group.atom_numbers);
    std::vector<int> wanted(expected);
    std::sort(actual.begin(), actual.end());
    std::sort(wanted.begin(), wanted.end());
    return actual == wanted;
}

bool Check_Domain_Decomposition_Semantics()
{
    FULL_NEIGHBOR_LIST list;
    if (!list.Initial(5, 6)) return false;

    HALF_LIST initial(3);
    initial.Set(0, {1, 3}, 1);
    initial.Set(1, {2, 4}, 1);
    initial.Set(2, {}, 0);
    if (!list.Build_From_Half(initial.rows.data(), 3, 5) ||
        list.atom_capacity != 5 || list.active_owned_atom_numbers != 3 ||
        !Expect_Row(list, 0, {1, 3}, 1) || !Expect_Row(list, 1, {0, 2, 4}, 1) ||
        !Expect_Row(list, 2, {1}, 0) || !Expect_Row(list, 3, {}, 0) ||
        !Expect_Row(list, 4, {}, 0))
        return false;

    HALF_LIST shrunk(2);
    shrunk.Set(0, {2}, 1);
    shrunk.Set(1, {}, 0);
    if (!list.Build_From_Half(shrunk.rows.data(), 2, 3) ||
        list.active_owned_atom_numbers != 2 || !Expect_Row(list, 0, {2}, 1) ||
        !Expect_Row(list, 1, {}, 0) || !Expect_Row(list, 2, {}, 0))
        return false;

    HALF_LIST grown(4);
    grown.Set(0, {}, 0);
    grown.Set(1, {}, 0);
    grown.Set(2, {3}, 0);
    grown.Set(3, {}, 0);
    if (!list.Build_From_Half(grown.rows.data(), 4, 4) ||
        list.active_owned_atom_numbers != 4 || !Expect_Row(list, 0, {}, 0) ||
        !Expect_Row(list, 1, {}, 0) || !Expect_Row(list, 2, {3}, 0) ||
        !Expect_Row(list, 3, {2}, 0))
        return false;

    if (!list.Build_From_Half(NULL, 0, 0) ||
        list.active_owned_atom_numbers != 0)
        return false;
    for (int row = 0; row < 4; ++row)
        if (!Expect_Row(list, row, {}, 0)) return false;

    if (list.Initial(6, 6) || !list.is_initialized || list.atom_capacity != 5)
        return false;
    list.Clear();
    list.Clear();
    return !list.is_initialized && list.h_nl == NULL && list.d_nl == NULL &&
           list.d_temp == NULL && list.d_overflow == NULL &&
           list.d_build_error == NULL && list.atom_capacity == 0 &&
           list.active_owned_atom_numbers == 0;
}

bool Check_Cutoff_Semantics()
{
    FULL_NEIGHBOR_LIST list;
    if (!list.Initial(5, 6)) return false;

    HALF_LIST half(3);
    half.Set(0, {1, 3}, 1);
    half.Set(1, {2, 4}, 1);
    half.Set(2, {4}, 1);
    const VECTOR coordinates[5] = {{0.0f, 0.0f, 0.0f},
                                   {1.0f, 0.0f, 0.0f},
                                   {4.0f, 0.0f, 0.0f},
                                   {2.0f, 0.0f, 0.0f},
                                   {8.0f, 0.0f, 0.0f}};
    const LTMatrix3 cell = {100.0f, 0.0f, 100.0f, 0.0f, 0.0f, 100.0f};
    const LTMatrix3 reciprocal = {0.01f, 0.0f, 0.01f, 0.0f, 0.0f, 0.01f};
    if (!list.Build_From_Half_With_Cutoff(half.rows.data(), 3, 5, coordinates,
                                          cell, reciprocal, 2.0f) ||
        !Expect_Row(list, 0, {1, 3}, 1) || !Expect_Row(list, 1, {0}, 0) ||
        !Expect_Row(list, 2, {}, 0))
        return false;

    VECTOR invalid_coordinates[5];
    std::copy(std::begin(coordinates), std::end(coordinates),
              invalid_coordinates);
    invalid_coordinates[4].x = std::numeric_limits<float>::quiet_NaN();
    if (list.Build_From_Half_With_Cutoff(half.rows.data(), 3, 5,
                                         invalid_coordinates, cell, reciprocal,
                                         2.0f) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_GEOMETRY ||
        list.last_error_atom != 4 || list.active_owned_atom_numbers != 0)
        return false;
    for (int row = 0; row < 3; ++row)
        if (!Expect_Row(list, row, {}, 0)) return false;

    LTMatrix3 invalid_cell = cell;
    invalid_cell.a21 = std::numeric_limits<float>::infinity();
    if (list.Build_From_Half_With_Cutoff(half.rows.data(), 3, 5, coordinates,
                                         invalid_cell, reciprocal, 2.0f) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_GEOMETRY ||
        list.active_owned_atom_numbers != 0)
        return false;

    const float invalid_cutoffs[] = {
        0.0f,
        -1.0f,
        std::numeric_limits<float>::infinity(),
        std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::denorm_min(),
        std::numeric_limits<float>::max(),
    };
    for (float cutoff : invalid_cutoffs)
    {
        if (list.Build_From_Half_With_Cutoff(half.rows.data(), 3, 5,
                                             coordinates, cell, reciprocal,
                                             cutoff) ||
            list.last_build_error !=
                FULL_NEIGHBOR_LIST::BUILD_INVALID_ARGUMENT ||
            list.active_owned_atom_numbers != 0)
            return false;
    }
    if (list.Build_From_Half_With_Cutoff(half.rows.data(), 3, 5, NULL, cell,
                                         reciprocal, 2.0f))
        return false;
    list.Clear();
    return true;
}

bool Check_Invalid_Half_Lists()
{
    FULL_NEIGHBOR_LIST list;
    if (!list.Initial(5, 3)) return false;

    if (list.Build_From_Half(NULL, 1, 1) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_ARGUMENT)
        return false;
    if (list.Build_From_Half(NULL, -1, 0) || list.Build_From_Half(NULL, 6, 6) ||
        list.Build_From_Half(NULL, 2, 1))
        return false;

    HALF_LIST half(2);
    half.Set(0, {}, 0);
    half.Set(1, {}, 0);
    half.Set(0, {1}, 0);
    if (!list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.active_owned_atom_numbers != 2)
        return false;
    half.rows[0].atom_numbers = -1;
    if (list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_HALF_COUNT ||
        list.last_error_atom != 0 || list.active_owned_atom_numbers != 0 ||
        !Expect_Row(list, 0, {}, 0) || !Expect_Row(list, 1, {}, 0))
        return false;

    half.rows[0].atom_numbers = 4;
    if (list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_HALF_COUNT)
        return false;

    half.rows[0].atom_numbers = 0;
    half.rows[0].ghost_numbers = -1;
    if (list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_HALF_COUNT)
        return false;
    half.rows[0].ghost_numbers = 1;
    if (list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_HALF_COUNT)
        return false;

    half.rows[0].atom_numbers = 1;
    half.rows[0].ghost_numbers = 0;
    half.rows[0].atom_serial = NULL;
    if (list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_INVALID_HALF_POINTER)
        return false;

    int invalid_neighbor = -1;
    half.rows[0].atom_serial = &invalid_neighbor;
    if (list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.last_build_error !=
            FULL_NEIGHBOR_LIST::BUILD_INVALID_NEIGHBOR_INDEX ||
        list.last_error_value != -1)
        return false;
    invalid_neighbor = 2;
    if (list.Build_From_Half(half.rows.data(), 2, 2) ||
        list.last_error_value != 2)
        return false;

    list.Clear();
    return true;
}

bool Check_Overflow_Is_Bounded()
{
    FULL_NEIGHBOR_LIST list;
    if (!list.Initial(4, 2)) return false;
    HALF_LIST half(4);
    half.Set(0, {1}, 0);
    half.Set(1, {}, 0);
    half.Set(2, {1}, 0);
    half.Set(3, {1}, 0);
    if (list.Build_From_Half(half.rows.data(), 4, 4) ||
        list.last_build_error != FULL_NEIGHBOR_LIST::BUILD_CAPACITY_EXCEEDED ||
        list.last_required_neighbor_capacity != 3 ||
        list.active_owned_atom_numbers != 0)
        return false;
    int overflow = 0;
    deviceMemcpy(&overflow, list.d_overflow, sizeof(overflow),
                 deviceMemcpyDeviceToHost);
    const bool rejected_safely =
        overflow == 3 && list.h_nl[1].atom_numbers == 0;
    list.Clear();
    if (!rejected_safely || !list.Initial(4, overflow) ||
        !list.Build_From_Half(half.rows.data(), 4, 4) ||
        !Expect_Row(list, 1, {0, 2, 3}, 0))
        return false;
    list.Clear();
    return true;
}
}  // namespace

int main()
{
    if (!Check_Domain_Decomposition_Semantics()) return 1;
    if (!Check_Cutoff_Semantics()) return 2;
    if (!Check_Invalid_Half_Lists()) return 3;
    if (!Check_Overflow_Is_Bounded()) return 4;
    std::puts("full-neighbor-list probe passed");
    return 0;
}
