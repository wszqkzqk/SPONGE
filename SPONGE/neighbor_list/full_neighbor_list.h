#ifndef FULL_NEIGHBOR_LIST_H
#define FULL_NEIGHBOR_LIST_H
#include "../common.h"
#include "../control.h"

struct FULL_NEIGHBOR_LIST
{
    enum BUILD_ERROR
    {
        BUILD_OK = 0,
        BUILD_NOT_INITIALIZED = 1,
        BUILD_INVALID_ARGUMENT = 2,
        BUILD_INVALID_HALF_COUNT = 3,
        BUILD_INVALID_HALF_POINTER = 4,
        BUILD_INVALID_NEIGHBOR_INDEX = 5,
        BUILD_ALLOCATION_FAILED = 6,
        BUILD_CAPACITY_EXCEEDED = 7,
        BUILD_INVALID_GEOMETRY = 8
    };

    bool is_initialized = false;
    int atom_capacity = 0;
    int active_owned_atom_numbers = 0;
    int max_neighbor_numbers = 0;
    int last_build_error = BUILD_OK;
    int last_error_atom = -1;
    int last_error_value = -1;
    int last_required_neighbor_capacity = 0;

    ATOM_GROUP* d_nl = NULL;
    ATOM_GROUP* h_nl = NULL;
    int* d_temp = NULL;
    int* d_overflow = NULL;
    int* d_build_error = NULL;

    bool Initial(int atom_capacity, int max_neighbor_numbers);

    bool Build_From_Half(const ATOM_GROUP* half_nl, int owned_atom_numbers,
                         int coordinate_numbers);

    bool Build_From_Half_With_Cutoff(const ATOM_GROUP* half_nl,
                                     int owned_atom_numbers,
                                     int coordinate_numbers, const VECTOR* crd,
                                     const LTMatrix3 cell,
                                     const LTMatrix3 rcell, float cutoff);

    const char* Last_Error_Message() const;

    // Make every previously built row unavailable without releasing storage.
    // This is used when an upstream half-list build fails before a new full
    // list can be constructed.
    void Invalidate_Active();

    void Clear();
};
#endif
