#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#define QC_CARTESIAN_HOST_DEVICE __host__ __device__
#else
#define QC_CARTESIAN_HOST_DEVICE
#endif

namespace qc_cartesian
{
// Cartesian components are ordered by descending lx and then descending ly.
// Return false for an invalid shell/component instead of fabricating a tuple;
// host initialization validates all production inputs before kernel launch.
QC_CARTESIAN_HOST_DEVICE inline bool Component(int l, int index, int& lx,
                                               int& ly, int& lz)
{
    if (l < 0 || index < 0) return false;
    for (lx = l; lx >= 0; --lx)
    {
        const int remaining = l - lx;
        const int count = remaining + 1;
        if (index < count)
        {
            ly = remaining - index;
            lz = index;
            return true;
        }
        index -= count;
    }
    return false;
}
}  // namespace qc_cartesian

#undef QC_CARTESIAN_HOST_DEVICE
