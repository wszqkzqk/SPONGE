#pragma once

#include "../common.h"

// Coordinates and the current charge of an exclusion partner that is owned
// by another PP domain and is not part of the spatial neighbor-list halo.
// PME exclusion corrections are long ranged, so these dependencies must be
// communicated independently of the cutoff-based ghost layout.
struct PME_EXCLUSION_DEPENDENCY_STATE
{
    VECTOR crd;
    float charge;
};
