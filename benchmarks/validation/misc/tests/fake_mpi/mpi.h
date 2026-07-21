#pragma once

#include <cstdio>
#include <cstdlib>

// Minimal compile-time MPI surface for the isolated domain-decomposition
// probe.  The probe exercises layout construction only; it never performs
// communication.
using MPI_Comm = int;
using MPI_Datatype = int;
using MPI_Op = int;
using MPI_Request = int;
struct MPI_Status
{
    int unused;
};

#define MPI_COMM_WORLD 0
#define MPI_BYTE 0
#define MPI_CHAR 0
#define MPI_DOUBLE 0
#define MPI_FLOAT 0
#define MPI_INT 0
#define MPI_SUM 0
#define MPI_SUCCESS 0
#define MPI_THREAD_SINGLE 0
#define MPI_THREAD_FUNNELED 1
#define MPI_THREAD_SERIALIZED 2
#define MPI_THREAD_MULTIPLE 3
#define MPI_IN_PLACE nullptr
#define MPI_STATUS_IGNORE nullptr

#ifndef FAKE_MPI_INIT_STATUS
#define FAKE_MPI_INIT_STATUS MPI_SUCCESS
#endif
#ifndef FAKE_MPI_PROVIDED
#define FAKE_MPI_PROVIDED MPI_THREAD_MULTIPLE
#endif
#ifndef FAKE_MPI_INITIALIZED
#define FAKE_MPI_INITIALIZED 1
#endif

inline int MPI_Abort(MPI_Comm, int)
{
    std::fputs("FAKE_MPI_ABORT_RETURNED\n", stderr);
    std::fflush(stderr);
    return MPI_SUCCESS;
}
inline int MPI_Allgather(...) { return MPI_SUCCESS; }
inline int MPI_Allgatherv(...) { return MPI_SUCCESS; }
inline int MPI_Allreduce(...) { return MPI_SUCCESS; }
inline int MPI_Alltoall(...) { return MPI_SUCCESS; }
inline int MPI_Alltoallv(...) { return MPI_SUCCESS; }
inline int MPI_Barrier(...) { return MPI_SUCCESS; }
inline int MPI_Bcast(...) { return MPI_SUCCESS; }
inline int MPI_Comm_rank(MPI_Comm, int* rank)
{
    *rank = 0;
    return MPI_SUCCESS;
}
inline int MPI_Comm_size(MPI_Comm, int* size)
{
    *size = 1;
    return MPI_SUCCESS;
}
inline int MPI_Comm_split(...) { return MPI_SUCCESS; }
inline int MPI_Finalize(...) { return MPI_SUCCESS; }
inline int MPI_Finalized(int* finalized)
{
    *finalized = 0;
    return MPI_SUCCESS;
}
inline int MPI_Gather(...) { return MPI_SUCCESS; }
inline int MPI_Gatherv(...) { return MPI_SUCCESS; }
inline int MPI_Init(...) { return MPI_SUCCESS; }
inline int MPI_Init_thread(int*, char***, int, int* provided)
{
    *provided = FAKE_MPI_PROVIDED;
    return FAKE_MPI_INIT_STATUS;
}
inline int MPI_Initialized(int* initialized)
{
    *initialized = FAKE_MPI_INITIALIZED;
    return MPI_SUCCESS;
}
inline int MPI_Irecv(...) { return MPI_SUCCESS; }
inline int MPI_Isend(...) { return MPI_SUCCESS; }
inline int MPI_Recv(...) { return MPI_SUCCESS; }
inline int MPI_Reduce(...) { return MPI_SUCCESS; }
inline int MPI_Send(...) { return MPI_SUCCESS; }
inline int MPI_Waitall(...) { return MPI_SUCCESS; }
