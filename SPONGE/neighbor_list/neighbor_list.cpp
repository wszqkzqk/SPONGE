#include "neighbor_list.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>

#define MAX_GRID_NEIGHBORS 192

namespace
{
bool Checked_Product(size_t lhs, size_t rhs, size_t* product)
{
    if (product == NULL ||
        (lhs != 0 && rhs > std::numeric_limits<size_t>::max() / lhs))
        return false;
    *product = lhs * rhs;
    return true;
}

bool Checked_Bytes(size_t count, size_t element_size, size_t* bytes)
{
    return Checked_Product(count, element_size, bytes);
}

bool Finite_Matrix(const LTMatrix3& matrix)
{
    return Float_Memory_Is_Finite(&matrix.a11) &&
           Float_Memory_Is_Finite(&matrix.a21) &&
           Float_Memory_Is_Finite(&matrix.a22) &&
           Float_Memory_Is_Finite(&matrix.a31) &&
           Float_Memory_Is_Finite(&matrix.a32) &&
           Float_Memory_Is_Finite(&matrix.a33);
}

bool Grid_Dimension(float box_extent, float grid_length, int* dimension)
{
    const double ratio =
        static_cast<double>(box_extent) / static_cast<double>(grid_length);
    const double floored = std::floor(ratio);
    if (!std::isfinite(floored) || floored < 1.0 ||
        floored > static_cast<double>(std::numeric_limits<int>::max()))
        return false;
    *dimension = static_cast<int>(floored);
    return true;
}

bool Grid_Shear_Shift(float shear, float grid_length, int64_t* shift)
{
    const double ratio =
        static_cast<double>(shear) / static_cast<double>(grid_length);
    if (!std::isfinite(ratio) ||
        ratio < static_cast<double>(std::numeric_limits<int>::min()) ||
        ratio > static_cast<double>(std::numeric_limits<int>::max()))
        return false;
    *shift = static_cast<int64_t>(ratio);
    return true;
}

bool Parse_Int_Exactly(const char* token, int* value)
{
    if (token == NULL || token[0] == '\0' || value == NULL) return false;
    errno = 0;
    char* end = NULL;
    const long long parsed = std::strtoll(token, &end, 10);
    if (errno == ERANGE || end == token || end[0] != '\0' ||
        parsed < std::numeric_limits<int>::min() ||
        parsed > std::numeric_limits<int>::max())
        return false;
    *value = static_cast<int>(parsed);
    return true;
}
}  // namespace

static __global__ void Find_Neighor_Grids_Device(
    int grid_numbers, int* neighbor_grid_numbers, int* neighbor_grids, int Nx,
    int Ny, int Nz, int64_t dxy, int64_t dxz, int64_t dyz, LTMatrix3 cell)
{
    SIMPLE_DEVICE_FOR(grid_i, grid_numbers)
    {
        const int64_t NyNz = static_cast<int64_t>(Ny) * Nz;
        int local = 0;
        int dx, dy, dz;
        int* neighbor_i = neighbor_grids + (size_t)MAX_GRID_NEIGHBORS * grid_i;
        const int64_t grid_i_x = grid_i / NyNz;
        const int64_t grid_i_y = (grid_i % NyNz) / Nz;
        const int64_t grid_i_z = grid_i % Nz;
        int cx, cy;
        for (int grid_j = 0; grid_j < grid_numbers; grid_j += 1)
        {
            const int64_t grid_j_x = grid_j / NyNz;
            const int64_t grid_j_y = (grid_j % NyNz) / Nz;
            const int64_t grid_j_z = grid_j % Nz;
            for (int k = 0; k < 27; k += 1)
            {
                dx = k / 9 - 1;
                dy = (k % 9) / 3 - 1;
                dz = k % 3 - 1;
                const int64_t temp_x = grid_i_x +
                                       static_cast<int64_t>(dx) * Nx +
                                       dy * dxy + dz * dxz - grid_j_x;
                const int64_t temp_y = grid_i_y +
                                       static_cast<int64_t>(dy) * Ny +
                                       dz * dyz - grid_j_y;
                const int64_t temp_z =
                    grid_i_z + static_cast<int64_t>(dz) * Nz - grid_j_z;
                cx = static_cast<double>(dy) * cell.a21 +
                                 static_cast<double>(dz) * cell.a31 ==
                             0.0
                         ? -2
                         : -3;
                cy = static_cast<double>(dz) * cell.a32 == 0.0 ? -2 : -3;
                if (temp_x <= 2 && temp_x >= cx && temp_y <= 2 &&
                    temp_y >= cy && temp_z <= 2 && temp_z >= -2)
                {
                    if (local >= MAX_GRID_NEIGHBORS)
                    {
                        local = -1;
                    }
                    else
                    {
                        neighbor_i[local++] = grid_j;
                    }
                    break;
                }
            }
            if (local < 0) break;
        }
        neighbor_grid_numbers[grid_i] = local;
    }
}

bool NEIGHBOR_LIST::GRIDS::Initial(CONTROLLER* controller,
                                   int max_atom_in_grid_numbers,
                                   int max_ghost_in_grid_numbers,
                                   LTMatrix3 cell, LTMatrix3 rcell,
                                   float grid_length)
{
    if (controller == NULL || max_atom_in_grid_numbers <= 0 ||
        max_ghost_in_grid_numbers <= 0 || !(grid_length > 0.0f) ||
        !Float_Memory_Is_Normal(&grid_length) || !Finite_Matrix(cell) ||
        !Finite_Matrix(rcell) || !(cell.a11 > 0.0f) || !(cell.a22 > 0.0f) ||
        !(cell.a33 > 0.0f) || !Grid_Dimension(cell.a11, grid_length, &Nx) ||
        !Grid_Dimension(cell.a22, grid_length, &Ny) ||
        !Grid_Dimension(cell.a33, grid_length, &Nz))
    {
        if (controller != NULL)
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "NEIGHBOR_LIST::GRIDS::Initial",
                "Reason:\n\tinvalid finite grid geometry or capacity\n");
        return false;
    }

    int64_t dxy = 0, dxz = 0, dyz = 0;
    size_t xy = 0, grid_count = 0;
    if (!Grid_Shear_Shift(cell.a21, grid_length, &dxy) ||
        !Grid_Shear_Shift(cell.a31, grid_length, &dxz) ||
        !Grid_Shear_Shift(cell.a32, grid_length, &dyz) ||
        !Checked_Product((size_t)Nx, (size_t)Ny, &xy) ||
        !Checked_Product(xy, (size_t)Nz, &grid_count) || grid_count == 0 ||
        grid_count > (size_t)std::numeric_limits<int>::max())
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "NEIGHBOR_LIST::GRIDS::Initial",
            "Reason:\n\tneighbor-grid dimensions exceed the supported "
            "signed-index representation\n");
        return false;
    }
    grid_numbers = static_cast<int>(grid_count);

    size_t neighbor_slots = 0, atom_slots = 0, ghost_slots = 0;
    size_t grid_int_bytes = 0, neighbor_bytes = 0, atom_bytes = 0;
    size_t ghost_bytes = 0, atom_crd_bytes = 0, ghost_crd_bytes = 0;
    if (!Checked_Product(grid_count, (size_t)MAX_GRID_NEIGHBORS,
                         &neighbor_slots) ||
        !Checked_Product(grid_count, (size_t)max_atom_in_grid_numbers,
                         &atom_slots) ||
        !Checked_Product(grid_count, (size_t)max_ghost_in_grid_numbers,
                         &ghost_slots) ||
        !Checked_Bytes(grid_count, sizeof(int), &grid_int_bytes) ||
        !Checked_Bytes(neighbor_slots, sizeof(int), &neighbor_bytes) ||
        !Checked_Bytes(atom_slots, sizeof(int), &atom_bytes) ||
        !Checked_Bytes(ghost_slots, sizeof(int), &ghost_bytes) ||
        !Checked_Bytes(atom_slots, sizeof(VECTOR), &atom_crd_bytes) ||
        !Checked_Bytes(ghost_slots, sizeof(VECTOR), &ghost_crd_bytes))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "NEIGHBOR_LIST::GRIDS::Initial",
            "Reason:\n\tneighbor-grid allocation size overflow\n");
        return false;
    }

    controller->printf("    initializing grids\n");
    controller->printf("        Nx: %d        Ny: %d        Nz: %d\n", Nx, Ny,
                       Nz);
    controller->printf("        Max number of atoms in one grid: %d\n",
                       max_atom_in_grid_numbers);
    controller->printf("        Max number of ghosts in one grid: %d\n",
                       max_ghost_in_grid_numbers);
    if (!Malloc_Safely((void**)&h_neighbor_grid_numbers, grid_int_bytes) ||
        !Malloc_Safely((void**)&h_neighbor_grids, neighbor_bytes) ||
        !Malloc_Safely((void**)&h_grid_atoms, atom_bytes) ||
        !Malloc_Safely((void**)&h_grid_atom_numbers, grid_int_bytes) ||
        !Malloc_Safely((void**)&h_grid_ghosts, ghost_bytes) ||
        !Malloc_Safely((void**)&h_grid_ghost_numbers, grid_int_bytes))
    {
        Clear();
        return false;
    }
    memset(h_grid_atom_numbers, 0, grid_int_bytes);
    memset(h_grid_ghost_numbers, 0, grid_int_bytes);
    if (!Device_Malloc_And_Copy_Safely((void**)&d_grid_ghosts, h_grid_ghosts,
                                       ghost_bytes) ||
        !Device_Malloc_And_Copy_Safely((void**)&d_grid_ghost_numbers,
                                       h_grid_ghost_numbers, grid_int_bytes) ||
        !Device_Malloc_Safely((void**)&d_grid_ghost_crd, ghost_crd_bytes) ||
        !Device_Malloc_And_Copy_Safely((void**)&d_neighbor_grid_numbers,
                                       h_neighbor_grid_numbers,
                                       grid_int_bytes) ||
        !Device_Malloc_And_Copy_Safely((void**)&d_neighbor_grids,
                                       h_neighbor_grids, neighbor_bytes) ||
        !Device_Malloc_And_Copy_Safely((void**)&d_grid_atoms, h_grid_atoms,
                                       atom_bytes) ||
        !Device_Malloc_Safely((void**)&d_grid_atom_crd, atom_crd_bytes) ||
        !Device_Malloc_And_Copy_Safely((void**)&d_grid_atom_numbers,
                                       h_grid_atom_numbers, grid_int_bytes))
    {
        Clear();
        return false;
    }

    Launch_Device_Kernel(Find_Neighor_Grids_Device,
                         (grid_numbers + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, grid_numbers,
                         d_neighbor_grid_numbers, d_neighbor_grids, Nx, Ny, Nz,
                         dxy, dxz, dyz, cell);
    deviceMemcpy(h_neighbor_grid_numbers, d_neighbor_grid_numbers,
                 grid_int_bytes, deviceMemcpyDeviceToHost);
    for (int grid = 0; grid < grid_numbers; ++grid)
    {
        if (h_neighbor_grid_numbers[grid] < 0 ||
            h_neighbor_grid_numbers[grid] > MAX_GRID_NEIGHBORS)
        {
            Clear();
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "NEIGHBOR_LIST::GRIDS::Initial",
                "Reason:\n\tthe fixed neighbor-grid topology bound was "
                "exceeded\n");
            return false;
        }
    }
    return true;
}

void NEIGHBOR_LIST::GRIDS::Clear()
{
    if (h_neighbor_grid_numbers != NULL || d_neighbor_grid_numbers != NULL)
        Free_Host_And_Device_Pointer((void**)&h_neighbor_grid_numbers,
                                     (void**)&d_neighbor_grid_numbers);
    if (h_neighbor_grids != NULL || d_neighbor_grids != NULL)
        Free_Host_And_Device_Pointer((void**)&h_neighbor_grids,
                                     (void**)&d_neighbor_grids);
    if (h_grid_atoms != NULL || d_grid_atoms != NULL)
        Free_Host_And_Device_Pointer((void**)&h_grid_atoms,
                                     (void**)&d_grid_atoms);
    if (h_grid_atom_numbers != NULL || d_grid_atom_numbers != NULL)
        Free_Host_And_Device_Pointer((void**)&h_grid_atom_numbers,
                                     (void**)&d_grid_atom_numbers);
    if (h_grid_ghosts != NULL || d_grid_ghosts != NULL)
        Free_Host_And_Device_Pointer((void**)&h_grid_ghosts,
                                     (void**)&d_grid_ghosts);
    if (h_grid_ghost_numbers != NULL || d_grid_ghost_numbers != NULL)
        Free_Host_And_Device_Pointer((void**)&h_grid_ghost_numbers,
                                     (void**)&d_grid_ghost_numbers);
    if (d_grid_ghost_crd != NULL)
        Free_Single_Device_Pointer((void**)&d_grid_ghost_crd);
    if (d_grid_atom_crd != NULL)
        Free_Single_Device_Pointer((void**)&d_grid_atom_crd);
    grid_numbers = 0;
    Nx = Ny = Nz = 0;
}

bool NEIGHBOR_LIST::UPDATOR::Initial(CONTROLLER* controller, int atom_numbers)
{
    if (controller == NULL || atom_numbers <= 0) return false;
    controller->printf("    initializing updator\n");
    time_recorder = controller->Get_Time_Recorder("neighbor searching");
    refresh_interval = 0;
    if (controller->Command_Exist("neighbor_list", "refresh_interval"))
    {
        controller->Check_Int("neighbor_list", "refresh_interval",
                              "NEIGHBOR_LIST::Initial");
        refresh_interval =
            atoi(controller->Command("neighbor_list", "refresh_interval"));
    }
    controller->printf(
        "        The interval to refresh the neighbor list: %d\n",
        refresh_interval);

    h_need_update = 1;
    if (refresh_interval <= 0)
    {
        skin_permit = 0.5f;
        if (controller->Command_Exist("neighbor_list", "skin_permit"))
        {
            controller->Check_Float("neighbor_list", "skin_permit",
                                    "NEIGHBOR_LIST::Initial");
            skin_permit =
                atof(controller->Command("neighbor_list", "skin_permit"));
        }
        controller->printf(
            "        The permit of the skin to refresh the neighbor list: %f\n",
            skin_permit);
        size_t coordinate_bytes = 0;
        if (!(skin_permit > 0.0f) || !Float_Memory_Is_Normal(&skin_permit) ||
            !Checked_Bytes((size_t)atom_numbers, sizeof(VECTOR),
                           &coordinate_bytes) ||
            !Device_Malloc_Safely((void**)&old_crd, coordinate_bytes))
        {
            Clear();
            return false;
        }
    }
    if (!Device_Malloc_And_Copy_Safely((void**)&d_need_update, &h_need_update,
                                       sizeof(int)))
    {
        Clear();
        return false;
    }
    return true;
}

static __global__ void Check_Refresh(int h_need_update, int atom_numbers,
                                     VECTOR* crd, VECTOR* crd_old,
                                     LTMatrix3 cell, LTMatrix3 rcell,
                                     int* d_need_refresh, float permit_square)
{
    SIMPLE_DEVICE_FOR(tid, atom_numbers)
    {
        VECTOR dr =
            Get_Periodic_Displacement(crd[tid], crd_old[tid], cell, rcell);
        if (dr * dr > permit_square)
        {
            d_need_refresh[0] = 1;
        }
    }
}

void NEIGHBOR_LIST::UPDATOR::Check(int atom_numbers, float skin, VECTOR* crd,
                                   LTMatrix3 cell, LTMatrix3 rcell)
{
    if (atom_numbers <= 0) return;
    Launch_Device_Kernel(Check_Refresh,
                         (atom_numbers + CONTROLLER::device_max_thread - 1) /
                             CONTROLLER::device_max_thread,
                         CONTROLLER::device_max_thread, 0, NULL, h_need_update,
                         atom_numbers, crd, old_crd, cell, rcell, d_need_update,
                         skin * skin * skin_permit * skin_permit);
}

static __device__ __forceinline__ void Record_Required_Capacity(
    int* required_capacity, int required)
{
    if (required <= 0) return;
#ifdef GPU_ARCH_NAME
    atomicMax(required_capacity, required);
#else
#pragma omp critical(sponge_neighbor_list_required_capacity)
    {
        if (required_capacity[0] < required) required_capacity[0] = required;
    }
#endif
}

static __device__ __forceinline__ int Required_After_Reservation(int slot,
                                                                 int count)
{
    const int64_t required = static_cast<int64_t>(slot) + count;
    // INT_MAX: std::numeric_limits<int>::max() is a host-only constexpr
    // under nvcc and cannot be called from this device function.
    return required >= INT_MAX ? INT_MAX : static_cast<int>(required);
}

static __device__ __forceinline__ unsigned int Neighbor_Float_Bits(float value)
{
#ifdef GPU_ARCH_NAME
    return __float_as_uint(value);
#else
    unsigned int bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires 32-bit IEEE-754 floats");
    memcpy(&bits, &value, sizeof(value));
#if defined(__GNUC__) || defined(__clang__)
    // Keep -ffast-math from deriving finite-only facts about the original
    // float and folding the integer exponent test below.
    __asm__ __volatile__("" : "+r"(bits));
#endif
    return bits;
#endif
}

static __device__ __forceinline__ UNSIGNED_INT_VECTOR
Neighbor_Vector_Bits(const VECTOR& value)
{
    return {Neighbor_Float_Bits(value.x), Neighbor_Float_Bits(value.y),
            Neighbor_Float_Bits(value.z)};
}

static __device__ __forceinline__ bool Neighbor_Vector_Bits_Are_Finite(
    const UNSIGNED_INT_VECTOR& bits)
{
    return (bits.uint_x & 0x7f800000U) != 0x7f800000U &&
           (bits.uint_y & 0x7f800000U) != 0x7f800000U &&
           (bits.uint_z & 0x7f800000U) != 0x7f800000U;
}

static __device__ __forceinline__ void Record_Neighbor_Update_Error(
    NEIGHBOR_LIST::UPDATE_ERROR_RECORD* update_error, int error, int local_atom,
    int global_atom, UNSIGNED_INT_VECTOR raw_bits,
    UNSIGNED_INT_VECTOR wrapped_bits, INT_VECTOR grid_index, int grid_id)
{
#ifdef GPU_ARCH_NAME
    if (atomicCAS(&update_error->header.code, NEIGHBOR_LIST::UPDATE_OK,
                  error) == NEIGHBOR_LIST::UPDATE_OK)
    {
        update_error->header.local_atom = local_atom;
        update_error->global_atom = global_atom;
        update_error->raw_bits = raw_bits;
        update_error->wrapped_bits = wrapped_bits;
        update_error->grid_index = grid_index;
        update_error->grid_id = grid_id;
    }
#else
#pragma omp critical(sponge_neighbor_list_update_error)
    {
        if (update_error->header.code == NEIGHBOR_LIST::UPDATE_OK)
        {
            update_error->header.code = error;
            update_error->header.local_atom = local_atom;
            update_error->global_atom = global_atom;
            update_error->raw_bits = raw_bits;
            update_error->wrapped_bits = wrapped_bits;
            update_error->grid_index = grid_index;
            update_error->grid_id = grid_id;
        }
    }
#endif
}

static __global__ void Clear_Bucket(const int* need, int grid_numbers,
                                    int* grid_atom_numbers,
                                    int* grid_ghost_numbers)
{
    if (need[0] == 0) return;
    SIMPLE_DEVICE_FOR(tid, grid_numbers)
    {
        grid_atom_numbers[tid] = 0;
        grid_ghost_numbers[tid] = 0;
    }
}

static __global__ void Put_Atom_In_Grids(
    const int* need, const int need_copy, const int* atom_local,
    const int atom_numbers, const int ghost_numbers, const int grid_numbers,
    const VECTOR* crd, VECTOR* old_crd, const LTMatrix3 cell,
    const LTMatrix3 rcell, const float grid_length, const int Nx, const int Ny,
    const int Nz, int* grid_atoms, int* grid_atom_numbers, VECTOR* grid_crd,
    ATOM_GROUP* nl, const int max_grid_atoms, int* neighbor_grid_overflow,
    int* grid_ghosts, int* grid_ghost_numbers, VECTOR* grid_ghost_crd,
    const int max_grid_ghosts, int* neighbor_grid_ghost_overflow,
    NEIGHBOR_LIST::UPDATE_ERROR_RECORD* update_error)
{
    if (need[0] == 0) return;
    SIMPLE_DEVICE_FOR(tid, atom_numbers + ghost_numbers)
    {
        const VECTOR raw_crd = crd[tid];
        VECTOR local_crd = raw_crd;
        const UNSIGNED_INT_VECTOR raw_bits = Neighbor_Vector_Bits(raw_crd);
        bool coordinate_is_finite = Neighbor_Vector_Bits_Are_Finite(raw_bits);
        const bool raw_coordinate_is_finite = coordinate_is_finite;
        UNSIGNED_INT_VECTOR wrapped_bits = {0U, 0U, 0U};
        if (raw_coordinate_is_finite && need_copy && tid < atom_numbers)
        {
            old_crd[tid] = raw_crd;
        }
        if (coordinate_is_finite)
        {
            float k3 = floorf(local_crd.z / cell.a33);
            local_crd.z -= k3 * cell.a33;
            local_crd.y -= k3 * cell.a32;
            float k2 = floorf(local_crd.y / cell.a22);
            local_crd.y -= k2 * cell.a22;
            local_crd.x -= k3 * cell.a31 + k2 * cell.a21;
            local_crd.x -= floorf(local_crd.x / cell.a11) * cell.a11;
            wrapped_bits = Neighbor_Vector_Bits(local_crd);
            coordinate_is_finite =
                Neighbor_Vector_Bits_Are_Finite(wrapped_bits);
        }
        if (!raw_coordinate_is_finite)
        {
            Record_Neighbor_Update_Error(
                update_error, NEIGHBOR_LIST::UPDATE_RAW_COORDINATE_NONFINITE,
                tid, atom_local == NULL ? -1 : atom_local[tid], raw_bits,
                wrapped_bits, {-1, -1, -1}, -1);
        }
        else if (!coordinate_is_finite)
        {
            Record_Neighbor_Update_Error(
                update_error,
                NEIGHBOR_LIST::UPDATE_WRAPPED_COORDINATE_NONFINITE, tid,
                atom_local == NULL ? -1 : atom_local[tid], raw_bits,
                wrapped_bits, {-1, -1, -1}, -1);
        }
        else
        {
            int nx = local_crd.x / grid_length;
            int ny = local_crd.y / grid_length;
            int nz = local_crd.z / grid_length;

            nx = nx < 0 ? 0 : (nx < Nx ? nx : Nx - 1);
            ny = ny < 0 ? 0 : (ny < Ny ? ny : Ny - 1);
            nz = nz < 0 ? 0 : (nz < Nz ? nz : Nz - 1);
            int grid_id = nx * Ny * Nz + ny * Nz + nz;
            if (grid_id < 0 || grid_id >= grid_numbers)
            {
                Record_Neighbor_Update_Error(
                    update_error, NEIGHBOR_LIST::UPDATE_GRID_INDEX_INVALID, tid,
                    atom_local == NULL ? -1 : atom_local[tid], raw_bits,
                    wrapped_bits, {nx, ny, nz}, grid_id);
            }
            else if (tid < atom_numbers)
            {
                int k1 = atomicAdd(grid_atom_numbers + grid_id, 1);
                if (k1 >= max_grid_atoms)
                {
                    Record_Required_Capacity(neighbor_grid_overflow,
                                             Required_After_Reservation(k1, 1));
                }
                else
                {
                    const size_t slot =
                        (size_t)max_grid_atoms * grid_id + (size_t)k1;
                    grid_atoms[slot] = tid;
                    grid_crd[slot] = raw_crd;
                    nl[tid].atom_numbers = 0;
                    nl[tid].ghost_numbers = 0;
                }
            }
            else
            {
                int k1 = atomicAdd(grid_ghost_numbers + grid_id, 1);
                if (k1 >= max_grid_ghosts)
                {
                    Record_Required_Capacity(neighbor_grid_ghost_overflow,
                                             Required_After_Reservation(k1, 1));
                }
                else
                {
                    const size_t slot =
                        (size_t)max_grid_ghosts * grid_id + (size_t)k1;
                    grid_ghosts[slot] = tid;
                    grid_ghost_crd[slot] = raw_crd;
                }
            }
        }
    }
}

#ifdef USE_GPU

static __global__ void Find_Neighbors_Gridly(
    int* atom_local, int atom_numbers, const int* need, int grid_numbers,
    int* grid_neighbor_numbers, int* grid_neighbors, VECTOR* grid_crd,
    LTMatrix3 cell, LTMatrix3 rcell, int max_atom_numbers_in_grid,
    ATOM_GROUP* nl, float cutoff_skin_square, int* grid_atom_numbers,
    int* grid_atoms, int max_neighbor_numbers, int* neighbor_list_overflow,
    VECTOR* grid_ghost_crd, int max_ghost_numbers_in_grid,
    int* grid_ghost_numbers, int* grid_ghosts,
    const int* neighbor_grid_overflow, const int* neighbor_grid_ghost_overflow,
    const NEIGHBOR_LIST::UPDATE_ERROR_RECORD* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error->header.code != NEIGHBOR_LIST::UPDATE_OK)
        return;
    extern __shared__ unsigned char shared_mem[];
    VECTOR* sh_crd = reinterpret_cast<VECTOR*>(shared_mem);
    int* sh_atoms = reinterpret_cast<int*>(sh_crd + max_atom_numbers_in_grid);

    const int lane = threadIdx.x & (warpSize - 1);
    int warps_per_block = blockDim.x / warpSize;
    if (warps_per_block == 0)
    {
        warps_per_block = 1;
    }
    const int warp_id = threadIdx.x / warpSize;
    const int lane_stride = blockDim.x < warpSize ? blockDim.x : warpSize;
    const int lane_index = blockDim.x < warpSize ? threadIdx.x : lane;

    device_mask_t warp_mask = FULL_MASK;
    if (blockDim.x < warpSize)
    {
        warp_mask = deviceLowerLaneMask(blockDim.x);
    }

    for (int grid_i = blockIdx.x; grid_i < grid_numbers; grid_i += gridDim.x)
    {
        int atom_numbers_in_grid_i = grid_atom_numbers[grid_i];
        if (atom_numbers_in_grid_i == 0)
        {
            __syncthreads();
            continue;
        }

        int* bucket_i = grid_atoms + (size_t)grid_i * max_atom_numbers_in_grid;
        VECTOR* grid_crd_i =
            grid_crd + (size_t)grid_i * max_atom_numbers_in_grid;

        for (int idx = threadIdx.x; idx < atom_numbers_in_grid_i;
             idx += blockDim.x)
        {
            sh_atoms[idx] = bucket_i[idx];
            sh_crd[idx] = grid_crd_i[idx];
        }
        __syncthreads();

        int neighbor_count = grid_neighbor_numbers[grid_i];
        if (neighbor_count == 0)
        {
            __syncthreads();
            continue;
        }
        for (int jj = warp_id; jj < neighbor_count; jj += warps_per_block)
        {
            int grid_j =
                grid_neighbors[(size_t)grid_i * MAX_GRID_NEIGHBORS + jj];
            int atom_numbers_in_grid_j = grid_atom_numbers[grid_j];
            if (atom_numbers_in_grid_j == 0)
            {
                continue;
            }

            int* bucket_j =
                grid_atoms + (size_t)grid_j * max_atom_numbers_in_grid;
            VECTOR* grid_crd_j =
                grid_crd + (size_t)grid_j * max_atom_numbers_in_grid;

            for (int j_base = 0; j_base < atom_numbers_in_grid_j;
                 j_base += lane_stride)
            {
                int j = j_base + lane_index;
                bool active = j < atom_numbers_in_grid_j;
                int atom_j = 0;
                int global_j = 0;
                VECTOR crd_j = {0, 0, 0};
                if (active)
                {
                    atom_j = bucket_j[j];
                    global_j = atom_local[atom_j];
                    crd_j = grid_crd_j[j];
                }

                for (int i = 0; i < atom_numbers_in_grid_i; ++i)
                {
                    int atom_i = sh_atoms[i];
                    int global_i = atom_local[atom_i];
                    bool is_neighbor = false;
                    if (active && global_j > global_i)
                    {
                        VECTOR dr = Get_Periodic_Displacement(sh_crd[i], crd_j,
                                                              cell, rcell);
                        float dr2 = dr * dr;
                        if (dr2 < cutoff_skin_square)
                        {
                            is_neighbor = true;
                        }
                    }

                    LaneMask mask = LaneGroup::And(
                        LaneGroup::Ballot(is_neighbor), LaneMask(warp_mask));
                    if (LaneGroup::Any(mask))
                    {
                        int count = LaneGroup::Count(mask);
                        int base_slot = 0;
                        int leader_lane = LaneGroup::First_Lane(mask);
                        if (lane == leader_lane)
                        {
                            base_slot =
                                atomicAdd(&nl[atom_i].atom_numbers, count);
                            if (static_cast<int64_t>(base_slot) + count >
                                max_neighbor_numbers)
                            {
                                Record_Required_Capacity(
                                    neighbor_list_overflow,
                                    Required_After_Reservation(base_slot,
                                                               count));
                            }
                        }
                        base_slot = deviceShfl(warp_mask, base_slot,
                                               leader_lane, lane_stride);

                        if (is_neighbor)
                        {
                            int rank = LaneGroup::Count(LaneGroup::And(
                                mask, LaneGroup::Lower_Lane_Mask()));
                            if (base_slot < max_neighbor_numbers &&
                                rank < max_neighbor_numbers - base_slot)
                            {
                                nl[atom_i].atom_serial[base_slot + rank] =
                                    atom_j;
                            }
                        }
                    }
                }
            }
        }

        __syncthreads();

        for (int jj = warp_id; jj < neighbor_count; jj += warps_per_block)
        {
            int grid_j =
                grid_neighbors[(size_t)grid_i * MAX_GRID_NEIGHBORS + jj];
            int ghost_numbers_in_grid_j = grid_ghost_numbers[grid_j];
            if (ghost_numbers_in_grid_j == 0)
            {
                continue;
            }

            int* bucket_j =
                grid_ghosts + (size_t)grid_j * max_ghost_numbers_in_grid;
            VECTOR* grid_ghost_crd_j =
                grid_ghost_crd + (size_t)grid_j * max_ghost_numbers_in_grid;

            for (int j_base = 0; j_base < ghost_numbers_in_grid_j;
                 j_base += lane_stride)
            {
                int j = j_base + lane_index;
                bool active = j < ghost_numbers_in_grid_j;
                int atom_j = 0;
                VECTOR crd_j = {0, 0, 0};
                if (active)
                {
                    atom_j = bucket_j[j];
                    crd_j = grid_ghost_crd_j[j];
                }

                for (int i = 0; i < atom_numbers_in_grid_i; ++i)
                {
                    int atom_i = sh_atoms[i];
                    bool is_neighbor = false;
                    if (active)
                    {
                        VECTOR dr = Get_Periodic_Displacement(sh_crd[i], crd_j,
                                                              cell, rcell);
                        float dr2 = dr * dr;
                        if (dr2 < cutoff_skin_square)
                        {
                            is_neighbor = true;
                        }
                    }

                    LaneMask mask = LaneGroup::And(
                        LaneGroup::Ballot(is_neighbor), LaneMask(warp_mask));
                    if (LaneGroup::Any(mask))
                    {
                        int count = LaneGroup::Count(mask);
                        int base_slot = 0;
                        int leader_lane = LaneGroup::First_Lane(mask);
                        if (lane == leader_lane)
                        {
                            base_slot =
                                atomicAdd(&nl[atom_i].atom_numbers, count);
                            atomicAdd(&nl[atom_i].ghost_numbers, count);
                            if (static_cast<int64_t>(base_slot) + count >
                                max_neighbor_numbers)
                            {
                                Record_Required_Capacity(
                                    neighbor_list_overflow,
                                    Required_After_Reservation(base_slot,
                                                               count));
                            }
                        }
                        base_slot = deviceShfl(warp_mask, base_slot,
                                               leader_lane, lane_stride);

                        if (is_neighbor)
                        {
                            int rank = LaneGroup::Count(LaneGroup::And(
                                mask, LaneGroup::Lower_Lane_Mask()));
                            if (base_slot < max_neighbor_numbers &&
                                rank < max_neighbor_numbers - base_slot)
                            {
                                nl[atom_i].atom_serial[base_slot + rank] =
                                    atom_j;
                            }
                        }
                    }
                }
            }
        }

        __syncthreads();
    }
}

#else
static __global__ void Find_Neighbors_Gridly(
    int* atom_local, int atom_numbers, const int* need, int grid_numbers,
    int* grid_neighbor_numbers, int* grid_neighbors, VECTOR* grid_crd,
    LTMatrix3 cell, LTMatrix3 rcell, int max_atom_numbers_in_grid,
    ATOM_GROUP* nl, float cutoff_skin_square, int* grid_atom_numbers,
    int* grid_atoms, int max_neighbor_numbers, int* neighbor_list_overflow,
    VECTOR* grid_ghost_crd, int max_ghost_numbers_in_grid,
    int* grid_ghost_numbers, int* grid_ghosts,
    const int* neighbor_grid_overflow, const int* neighbor_grid_ghost_overflow,
    const NEIGHBOR_LIST::UPDATE_ERROR_RECORD* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error->header.code != NEIGHBOR_LIST::UPDATE_OK)
        return;
#pragma omp parallel for schedule(dynamic)
    for (int grid_i = 0; grid_i < grid_numbers; grid_i++)

    {
        int* bucket_i = grid_atoms + (size_t)grid_i * max_atom_numbers_in_grid;
        int atom_numbers_in_grid_i = grid_atom_numbers[grid_i];
        VECTOR* grid_crd_i =
            grid_crd + (size_t)grid_i * max_atom_numbers_in_grid;
        for (int jj = 0; jj < grid_neighbor_numbers[grid_i]; ++jj)
        {
            int grid_j =
                grid_neighbors[jj + (size_t)MAX_GRID_NEIGHBORS * grid_i];
            int* bucket_j =
                grid_atoms + (size_t)grid_j * max_atom_numbers_in_grid;
            int atom_numbers_in_grid_j = grid_atom_numbers[grid_j];

            VECTOR* grid_crd_j =
                grid_crd + (size_t)grid_j * max_atom_numbers_in_grid;
            for (int i = 0; i < atom_numbers_in_grid_i; i++)
            {
                int atom_i = bucket_i[i];
                int global_i = atom_local[atom_i];
                VECTOR crd_i = grid_crd_i[i];

                int* nl_atom_numbers_ptr = &nl[atom_i].atom_numbers;
                int* nl_atom_serial_ptr = nl[atom_i].atom_serial;
                for (int j = 0; j < atom_numbers_in_grid_j; ++j)
                {
                    int atom_j = bucket_j[j];
                    if (atom_local[atom_j] <= global_i) continue;
                    VECTOR crd_j = grid_crd_j[j];
                    VECTOR dr =
                        Get_Periodic_Displacement(crd_i, crd_j, cell, rcell);
                    float dr2 = dr * dr;
                    if (dr2 < cutoff_skin_square)
                    {
                        const int slot = atomicAdd(nl_atom_numbers_ptr, 1);
                        if (slot < max_neighbor_numbers)
                        {
                            nl_atom_serial_ptr[slot] = atom_j;
                        }
                        else
                        {
                            Record_Required_Capacity(
                                neighbor_list_overflow,
                                Required_After_Reservation(slot, 1));
                        }
                    }
                }
            }
        }

        for (int jj = 0; jj < grid_neighbor_numbers[grid_i]; ++jj)
        {
            int grid_j =
                grid_neighbors[jj + (size_t)MAX_GRID_NEIGHBORS * grid_i];
            int ghost_numbers_in_grid_j = grid_ghost_numbers[grid_j];
            if (ghost_numbers_in_grid_j == 0) continue;
            int* bucket_j =
                grid_ghosts + (size_t)grid_j * max_ghost_numbers_in_grid;
            VECTOR* grid_ghost_crd_j =
                grid_ghost_crd + (size_t)grid_j * max_ghost_numbers_in_grid;
            for (int i = 0; i < atom_numbers_in_grid_i; i++)
            {
                int atom_i = bucket_i[i];
                VECTOR crd_i = grid_crd_i[i];

                int* nl_atom_numbers_ptr = &nl[atom_i].atom_numbers;
                int* nl_ghost_numbers_ptr = &nl[atom_i].ghost_numbers;
                int* nl_atom_serial_ptr = nl[atom_i].atom_serial;
                for (int j = 0; j < ghost_numbers_in_grid_j; ++j)
                {
                    int atom_j = bucket_j[j];
                    VECTOR crd_j = grid_ghost_crd_j[j];
                    VECTOR dr =
                        Get_Periodic_Displacement(crd_i, crd_j, cell, rcell);
                    float dr2 = dr * dr;
                    if (dr2 < cutoff_skin_square)
                    {
                        const int slot = atomicAdd(nl_atom_numbers_ptr, 1);
                        atomicAdd(nl_ghost_numbers_ptr, 1);
                        if (slot < max_neighbor_numbers)
                        {
                            nl_atom_serial_ptr[slot] = atom_j;
                        }
                        else
                        {
                            Record_Required_Capacity(
                                neighbor_list_overflow,
                                Required_After_Reservation(slot, 1));
                        }
                    }
                }
            }
        }
    }
}
#endif

static __global__ void Delete_Excluded_Atoms_Serial_In_Neighbor_List(
    const int* need, const int local_atom_numbers, int* atom_local,
    ATOM_GROUP* nl, const int* excluded_list_start, const int* excluded_list,
    const int* excluded_atom_numbers, const int* neighbor_grid_overflow,
    const int* neighbor_grid_ghost_overflow, const int* neighbor_list_overflow,
    const NEIGHBOR_LIST::UPDATE_ERROR_RECORD* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        neighbor_list_overflow[0] != 0 ||
        update_error->header.code != NEIGHBOR_LIST::UPDATE_OK)
        return;
    SIMPLE_DEVICE_FOR(atom_i_local, local_atom_numbers)
    {
        int atom_i = atom_local[atom_i_local];
        int excluded_number = excluded_atom_numbers[atom_i];
        int list_start_i =
            excluded_number > 0 ? excluded_list_start[atom_i] : 0;
        int list_end_i = list_start_i + excluded_number;
        int atom_min_i = excluded_number > 0 ? excluded_list[list_start_i] : 0;
        int atom_max_i =
            excluded_number > 0 ? excluded_list[list_end_i - 1] : -1;

        ATOM_GROUP nl_i = nl[atom_i_local];
        int atomnumbers_in_nl_lin = nl_i.atom_numbers;
        int atom_j_local, atom_j;
        for (int i = 0; i < atomnumbers_in_nl_lin; ++i)
        {
            atom_j_local = nl_i.atom_serial[i];
            atom_j = atom_local[atom_j_local];
            bool is_excluded = false;

            if (excluded_number > 0 && atom_j >= atom_min_i &&
                atom_j <= atom_max_i)
            {
                for (int j = list_start_i; j < list_end_i; ++j)
                {
                    if (atom_j == excluded_list[j])
                    {
                        is_excluded = true;
                        break;
                    }
                }
            }

            if (!is_excluded && atom_j_local >= local_atom_numbers)
            {
                int excluded_number_j = excluded_atom_numbers[atom_j];
                if (excluded_number_j > 0)
                {
                    int list_start_j = excluded_list_start[atom_j];
                    int list_end_j = list_start_j + excluded_number_j;
                    int atom_min_j = excluded_list[list_start_j];
                    int atom_max_j = excluded_list[list_end_j - 1];
                    if (atom_i >= atom_min_j && atom_i <= atom_max_j)
                    {
                        for (int j = list_start_j; j < list_end_j; ++j)
                        {
                            if (atom_i == excluded_list[j])
                            {
                                is_excluded = true;
                                break;
                            }
                        }
                    }
                }
            }

            if (is_excluded)
            {
                if (atom_j_local >= local_atom_numbers)
                {
                    nl_i.ghost_numbers--;
                }
                atomnumbers_in_nl_lin = atomnumbers_in_nl_lin - 1;
                nl_i.atom_serial[i] = nl_i.atom_serial[atomnumbers_in_nl_lin];
                i--;
            }
        }
        nl[atom_i_local].atom_numbers = atomnumbers_in_nl_lin;
        nl[atom_i_local].ghost_numbers = nl_i.ghost_numbers;
    }
}

void NEIGHBOR_LIST::UPDATOR::Update(
    int* atom_local, int local_atom_numbers, int ghost_numbers, int need_copy,
    VECTOR* crd, LTMatrix3 cell, LTMatrix3 rcell, NEIGHBOR_LIST::GRIDS* grids,
    int max_atom_in_grid_numbers, int max_ghost_in_grid_numbers,
    int max_neighbor_numbers, float grid_length, int* d_neighbor_grid_overflow,
    int* d_neighbor_grid_ghost_overflow, int* d_neighbor_list_overflow,
    NEIGHBOR_LIST::UPDATE_ERROR_RECORD* d_update_error, ATOM_GROUP* d_nl,
    int* excluded_list_start, int* excluded_list, int* excluded_numbers)
{
    int total_atom_numbers = local_atom_numbers + ghost_numbers;
    if (total_atom_numbers <= 0) return;
    Launch_Device_Kernel(
        Clear_Bucket,
        (grids->grid_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, d_need_update,
        grids->grid_numbers, grids->d_grid_atom_numbers,
        grids->d_grid_ghost_numbers);

    Launch_Device_Kernel(
        Put_Atom_In_Grids,
        (total_atom_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, d_need_update, need_copy,
        atom_local, local_atom_numbers, ghost_numbers, grids->grid_numbers, crd,
        old_crd, cell, rcell, grid_length, grids->Nx, grids->Ny, grids->Nz,
        grids->d_grid_atoms, grids->d_grid_atom_numbers, grids->d_grid_atom_crd,
        d_nl, max_atom_in_grid_numbers, d_neighbor_grid_overflow,
        grids->d_grid_ghosts, grids->d_grid_ghost_numbers,
        grids->d_grid_ghost_crd, max_ghost_in_grid_numbers,
        d_neighbor_grid_ghost_overflow, d_update_error);
    Launch_Device_Kernel(
        Find_Neighbors_Gridly, grids->grid_numbers,
        CONTROLLER::device_max_thread,
        (size_t)(max_atom_in_grid_numbers * (sizeof(VECTOR) + sizeof(int))),
        NULL, atom_local, local_atom_numbers, d_need_update,
        grids->grid_numbers, grids->d_neighbor_grid_numbers,
        grids->d_neighbor_grids, grids->d_grid_atom_crd, cell, rcell,
        max_atom_in_grid_numbers, d_nl, grid_length * grid_length * 4.0f,
        grids->d_grid_atom_numbers, grids->d_grid_atoms, max_neighbor_numbers,
        d_neighbor_list_overflow, grids->d_grid_ghost_crd,
        max_ghost_in_grid_numbers, grids->d_grid_ghost_numbers,
        grids->d_grid_ghosts, d_neighbor_grid_overflow,
        d_neighbor_grid_ghost_overflow, d_update_error);

    if (local_atom_numbers > 0)
    {
        Launch_Device_Kernel(
            Delete_Excluded_Atoms_Serial_In_Neighbor_List,
            (local_atom_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, d_need_update,
            local_atom_numbers, atom_local, d_nl, excluded_list_start,
            excluded_list, excluded_numbers, d_neighbor_grid_overflow,
            d_neighbor_grid_ghost_overflow, d_neighbor_list_overflow,
            d_update_error);
    }
}

void NEIGHBOR_LIST::UPDATOR::Clear()
{
    if (d_need_update != NULL)
        Free_Host_And_Device_Pointer(NULL, (void**)&d_need_update);
    if (old_crd != NULL) Free_Single_Device_Pointer((void**)&old_crd);
    time_recorder = NULL;
    refresh_interval = 0;
    h_need_update = 0;
}

void NEIGHBOR_LIST::Initial(CONTROLLER* controller, int atom_numbers,
                            float cutoff, float skin, LTMatrix3 cell,
                            LTMatrix3 rcell)
{
    if (controller == NULL) return;
    if (is_initialized || h_nl != NULL || d_nl != NULL || d_temp != NULL)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "NEIGHBOR_LIST::Initial",
            "Reason:\n\tneighbor list was initialized without first clearing "
            "its previous storage\n");
        return;
    }
    const double grid_length_double =
        0.5 * (static_cast<double>(cutoff) + static_cast<double>(skin));
    const float grid_length = static_cast<float>(grid_length_double);
    if (atom_numbers <= 0 || !(cutoff > 0.0f) ||
        !Float_Memory_Is_Normal(&cutoff) || !(skin >= 0.0f) ||
        !Float_Memory_Is_Zero_Or_Normal(&skin) ||
        !std::isfinite(grid_length_double) || !(grid_length > 0.0f) ||
        !Float_Memory_Is_Normal(&grid_length))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "NEIGHBOR_LIST::Initial",
            "Reason:\n\tatom count, cutoff, and skin must define a finite "
            "positive neighbor grid\n");
        return;
    }

    this->atom_numbers = atom_numbers;
    this->cutoff = cutoff;
    this->skin = skin;
    controller->printf("START INITIALIZING NEIGHBOR LIST:\n");

    throw_error_when_overflow = 0;
    if (controller->Command_Exist("neighbor_list", "throw_error_when_overflow"))
        throw_error_when_overflow =
            controller->Get_Bool("neighbor_list", "throw_error_when_overflow",
                                 "NEIGHBOR_LIST::Initial");

    max_neighbor_numbers = 1200;
    if (controller->Command_Exist("neighbor_list", "max_neighbor_numbers"))
    {
        if (!Parse_Int_Exactly(
                controller->Command("neighbor_list", "max_neighbor_numbers"),
                &max_neighbor_numbers))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "NEIGHBOR_LIST::Initial",
                "Reason:\n\tneighbor_list_max_neighbor_numbers is not a "
                "representable integer\n");
            return;
        }
    }
    if (max_neighbor_numbers <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "NEIGHBOR_LIST::Initial",
            "Reason:\n\tneighbor_list_max_neighbor_numbers must be "
            "positive\n");
    }
    controller->printf("    Max number of neighbors for one atom: %d\n",
                       max_neighbor_numbers);

    max_atom_in_grid_numbers = 150;
    if (controller->Command_Exist("neighbor_list", "max_atom_in_grid_numbers"))
    {
        if (!Parse_Int_Exactly(controller->Command("neighbor_list",
                                                   "max_atom_in_grid_numbers"),
                               &max_atom_in_grid_numbers))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "NEIGHBOR_LIST::Initial",
                "Reason:\n\tneighbor_list_max_atom_in_grid_numbers is not a "
                "representable integer\n");
            return;
        }
    }
    max_ghost_in_grid_numbers = 150;
    if (controller->Command_Exist("neighbor_list", "max_ghost_in_grid_numbers"))
    {
        if (!Parse_Int_Exactly(controller->Command("neighbor_list",
                                                   "max_ghost_in_grid_numbers"),
                               &max_ghost_in_grid_numbers))
        {
            controller->Throw_SPONGE_Error(
                spongeErrorValueErrorCommand, "NEIGHBOR_LIST::Initial",
                "Reason:\n\tneighbor_list_max_ghost_in_grid_numbers is not a "
                "representable integer\n");
            return;
        }
    }
    if (max_atom_in_grid_numbers <= 0 || max_ghost_in_grid_numbers <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorValueErrorCommand, "NEIGHBOR_LIST::Initial",
            "Reason:\n\tneighbor-list atom and ghost grid capacities must "
            "be positive\n");
    }

    size_t group_bytes = 0, neighbor_slots = 0, neighbor_bytes = 0;
    if (!Checked_Bytes((size_t)atom_numbers, sizeof(ATOM_GROUP),
                       &group_bytes) ||
        !Checked_Product((size_t)atom_numbers, (size_t)max_neighbor_numbers,
                         &neighbor_slots) ||
        !Checked_Bytes(neighbor_slots, sizeof(int), &neighbor_bytes))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorOverflow, "NEIGHBOR_LIST::Initial",
            "Reason:\n\tneighbor-list allocation size overflow\n");
        return;
    }

    h_neighbor_grid_overflow = 0;
    h_neighbor_list_overflow = 0;
    h_neighbor_grid_ghost_overflow = 0;
    if (!Device_Malloc_And_Copy_Safely((void**)&d_neighbor_grid_overflow,
                                       &h_neighbor_grid_overflow,
                                       sizeof(int)) ||
        !Device_Malloc_And_Copy_Safely((void**)&d_neighbor_list_overflow,
                                       &h_neighbor_list_overflow,
                                       sizeof(int)) ||
        !Device_Malloc_And_Copy_Safely((void**)&d_neighbor_grid_ghost_overflow,
                                       &h_neighbor_grid_ghost_overflow,
                                       sizeof(int)) ||
        !Device_Malloc_Safely((void**)&d_update_error,
                              sizeof(UPDATE_ERROR_RECORD)) ||
        !Malloc_Safely((void**)&h_nl, group_bytes) ||
        !Device_Malloc_Safely((void**)&d_temp, neighbor_bytes))
    {
        Clear();
        controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                       "NEIGHBOR_LIST::Initial");
        return;
    }
    deviceMemset(d_update_error, 0, sizeof(UPDATE_ERROR_RECORD));
    for (int i = 0; i < atom_numbers; ++i)
    {
        h_nl[i].atom_numbers = 0;
        h_nl[i].ghost_numbers = 0;
        h_nl[i].atom_serial = d_temp + (size_t)max_neighbor_numbers * (size_t)i;
    }
    if (!Device_Malloc_And_Copy_Safely((void**)&d_nl, h_nl, group_bytes) ||
        !grids.Initial(controller, max_atom_in_grid_numbers,
                       max_ghost_in_grid_numbers, cell, rcell, grid_length) ||
        !updator.Initial(controller, atom_numbers))
    {
        Clear();
        controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                       "NEIGHBOR_LIST::Initial");
        return;
    }

    if (grids.Nx <= 2 || grids.Ny <= 2 || grids.Nz <= 2)
    {
        controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                       "NEIGHBOR_LIST::Initial",
                                       "the box is too small.");
    }
    if (this->cutoff_full > 0.0f)
    {
        controller->printf("    cutoff_full (from module): %f\n",
                           this->cutoff_full);
    }
    controller->printf("    is_needed_half: %s\n",
                       is_needed_half ? "true" : "false");
    controller->printf("    is_needed_full: %s\n",
                       is_needed_full ? "true" : "false");

    if (this->is_needed_full)
    {
        controller->printf("    Initializing full neighbor list...\n");
        if (!full_neighbor_list.Initial(atom_numbers, max_neighbor_numbers))
        {
            controller->Throw_Formatted_SPONGE_Error(
                full_neighbor_list.last_build_error ==
                        FULL_NEIGHBOR_LIST::BUILD_ALLOCATION_FAILED
                    ? spongeErrorMallocFailed
                    : spongeErrorSimulationBreakDown,
                "NEIGHBOR_LIST::Initial",
                "Reason:\n\tfailed to initialize full neighbor list: %s "
                "(capacity=%d, max_neighbors=%d)\n",
                full_neighbor_list.Last_Error_Message(), atom_numbers,
                max_neighbor_numbers);
        }
    }

    active_local_atom_numbers = 0;
    is_initialized = true;

    controller->printf("END INITIALIZING NEIGHBOR LIST\n\n");
}

static __global__ void Clear_Neighbor_List_Rows(ATOM_GROUP* nl, int row_numbers)
{
    SIMPLE_DEVICE_FOR(row, row_numbers)
    {
        nl[row].atom_numbers = 0;
        nl[row].ghost_numbers = 0;
    }
}

static void Invalidate_Half_Neighbor_List(NEIGHBOR_LIST* list,
                                          int requested_active_rows)
{
    int clear_rows = list->active_local_atom_numbers;
    if (requested_active_rows >= 0 &&
        requested_active_rows <= list->atom_numbers)
        clear_rows = std::max(clear_rows, requested_active_rows);
    if (list->d_nl != NULL && clear_rows > 0)
    {
        Launch_Device_Kernel(Clear_Neighbor_List_Rows,
                             (clear_rows + CONTROLLER::device_max_thread - 1) /
                                 CONTROLLER::device_max_thread,
                             CONTROLLER::device_max_thread, 0, NULL, list->d_nl,
                             clear_rows);
    }
    list->active_local_atom_numbers = 0;
}

bool NEIGHBOR_LIST::Update(int* atom_local, int local_atom_numbers,
                           int ghost_numbers, VECTOR* crd, LTMatrix3 cell,
                           LTMatrix3 rcell, int step, int update,
                           int* excluded_list_start, int* excluded_list,
                           int* excluded_numbers)
{
    const bool invalid_geometry = !Finite_Matrix(cell) ||
                                  !Finite_Matrix(rcell) || !(cell.a11 > 0.0f) ||
                                  !(cell.a22 > 0.0f) || !(cell.a33 > 0.0f);
    if (!is_initialized || local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > atom_numbers ||
        local_atom_numbers > std::numeric_limits<int>::max() - ghost_numbers ||
        (local_atom_numbers > 0 && atom_local == NULL) ||
        (local_atom_numbers + ghost_numbers > 0 && crd == NULL) ||
        (is_needed_full && !full_neighbor_list.is_initialized) ||
        invalid_geometry ||
        (update != NEIGHBOR_LIST_UPDATE_PARAMETER::FORCED_UPDATE &&
         update != NEIGHBOR_LIST_UPDATE_PARAMETER::CONDITIONAL_UPDATE))
    {
        last_update_error =
            invalid_geometry ? UPDATE_INVALID_PERIODIC_GEOMETRY : UPDATE_OK;
        last_error_atom = -1;
        if (is_initialized)
        {
            Invalidate_Half_Neighbor_List(this, local_atom_numbers);
            full_neighbor_list.Invalidate_Active();
        }
        return false;
    }
    updator.time_recorder->Start();
    // Overflow flags describe this build only.  Leaving them sticky makes a
    // later capacity check unable to distinguish a successful retry from the
    // truncated build that caused it.
    deviceMemset(d_neighbor_grid_overflow, 0, sizeof(int));
    deviceMemset(d_neighbor_grid_ghost_overflow, 0, sizeof(int));
    deviceMemset(d_neighbor_list_overflow, 0, sizeof(int));
    h_neighbor_grid_overflow = 0;
    h_neighbor_grid_ghost_overflow = 0;
    h_neighbor_list_overflow = 0;
    last_update_error = UPDATE_OK;
    last_error_atom = -1;
    deviceMemset(d_update_error, 0, sizeof(UPDATE_ERROR_RECORD));
    if (full_neighbor_list.is_initialized)
    {
        deviceMemset(full_neighbor_list.d_overflow, 0, sizeof(int));
    }
    if (update == NEIGHBOR_LIST_UPDATE_PARAMETER::FORCED_UPDATE)
    {
        deviceMemset(updator.d_need_update, -1, sizeof(int));
    }
    else if (updator.refresh_interval <= 0)
    {
        deviceMemset(updator.d_need_update, 0, sizeof(int));
        updator.Check(local_atom_numbers, skin, crd, cell, rcell);
    }
    else if (Next_Step_Is_Interval_Boundary(step,
                                            updator.refresh_interval))
    {
        deviceMemset(updator.d_need_update, -1, sizeof(int));
    }
    else
    {
        deviceMemset(updator.d_need_update, 0, sizeof(int));
    }
    if (this->is_needed_half && local_atom_numbers + ghost_numbers > 0)
    {
        updator.Update(atom_local, local_atom_numbers, ghost_numbers,
                       updator.refresh_interval <= 0, crd, cell, rcell, &grids,
                       max_atom_in_grid_numbers, max_ghost_in_grid_numbers,
                       max_neighbor_numbers, 0.5f * (cutoff + skin),
                       d_neighbor_grid_overflow, d_neighbor_grid_ghost_overflow,
                       d_neighbor_list_overflow, d_update_error, this->d_nl,
                       excluded_list_start, excluded_list, excluded_numbers);
    }

    deviceMemcpy(&h_neighbor_grid_overflow, d_neighbor_grid_overflow,
                 sizeof(int), deviceMemcpyDeviceToHost);
    deviceMemcpy(&h_neighbor_grid_ghost_overflow,
                 d_neighbor_grid_ghost_overflow, sizeof(int),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&h_neighbor_list_overflow, d_neighbor_list_overflow,
                 sizeof(int), deviceMemcpyDeviceToHost);
    UPDATE_ERROR_HEADER update_error = {UPDATE_OK, -1};
    deviceMemcpy(&update_error, d_update_error, sizeof(update_error),
                 deviceMemcpyDeviceToHost);
    last_update_error = update_error.code;
    last_error_atom =
        last_update_error == UPDATE_OK ? -1 : update_error.local_atom;
    if (last_update_error != UPDATE_OK || h_neighbor_grid_overflow != 0 ||
        h_neighbor_grid_ghost_overflow != 0 || h_neighbor_list_overflow != 0)
    {
        Invalidate_Half_Neighbor_List(this, local_atom_numbers);
        full_neighbor_list.Invalidate_Active();
        updator.time_recorder->Stop();
        return false;
    }

    if (this->is_needed_full && full_neighbor_list.is_initialized)
    {
        const int coordinate_numbers = local_atom_numbers + ghost_numbers;
        if (this->cutoff_full > 0.0f)
        {
            if (!full_neighbor_list.Build_From_Half_With_Cutoff(
                    this->d_nl, local_atom_numbers, coordinate_numbers, crd,
                    cell, rcell, this->cutoff_full + skin))
            {
                Invalidate_Half_Neighbor_List(this, local_atom_numbers);
                updator.time_recorder->Stop();
                return false;
            }
        }
        else
        {
            if (!full_neighbor_list.Build_From_Half(
                    this->d_nl, local_atom_numbers, coordinate_numbers))
            {
                Invalidate_Half_Neighbor_List(this, local_atom_numbers);
                updator.time_recorder->Stop();
                return false;
            }
        }
    }
    active_local_atom_numbers = is_needed_half ? local_atom_numbers : 0;
    updator.time_recorder->Stop();
    return true;
}

bool NEIGHBOR_LIST::Update_With_Overflow_Recovery(
    CONTROLLER* controller, int* atom_local, int local_atom_numbers,
    int ghost_numbers, VECTOR* crd, LTMatrix3 cell, LTMatrix3 rcell, int step,
    int update, int* excluded_list_start, int* excluded_list,
    int* excluded_numbers)
{
    if (controller == NULL) return false;
    if (!is_initialized)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
            "Reason:\n\tneighbor list was updated before initialization\n");
        return false;
    }
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > atom_numbers ||
        local_atom_numbers > std::numeric_limits<int>::max() - ghost_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
            "Reason:\n\tinvalid local/ghost atom counts %d/%d for full-list "
            "owned capacity %d\n",
            local_atom_numbers, ghost_numbers, atom_numbers);
        return false;
    }
    const int coordinate_numbers = local_atom_numbers + ghost_numbers;
    if ((local_atom_numbers > 0 && atom_local == NULL) ||
        (coordinate_numbers > 0 && crd == NULL))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
            "Reason:\n\tnull local-index or coordinate storage for a "
            "nonempty neighbor-list state\n");
        return false;
    }
    if (is_needed_full && !full_neighbor_list.is_initialized)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
            "Reason:\n\tfull neighbor list is required but not initialized\n");
        return false;
    }

    bool rebound = false;
    bool first_attempt = true;
    while (true)
    {
        if (Update(atom_local, local_atom_numbers, ghost_numbers, crd, cell,
                   rcell, step, first_attempt ? update : FORCED_UPDATE,
                   excluded_list_start, excluded_list, excluded_numbers))
            return rebound;
        first_attempt = false;

        const bool full_capacity_failure =
            is_needed_full && full_neighbor_list.last_build_error ==
                                  FULL_NEIGHBOR_LIST::BUILD_CAPACITY_EXCEEDED;
        const int required_full_neighbors =
            full_capacity_failure
                ? full_neighbor_list.last_required_neighbor_capacity
                : 0;
        const bool grid_atom_capacity_failure =
            h_neighbor_grid_overflow > max_atom_in_grid_numbers;
        const bool grid_ghost_capacity_failure =
            h_neighbor_grid_ghost_overflow > max_ghost_in_grid_numbers;
        const bool half_capacity_failure =
            h_neighbor_list_overflow > max_neighbor_numbers;
        const bool capacity_failure =
            grid_atom_capacity_failure || grid_ghost_capacity_failure ||
            half_capacity_failure || full_capacity_failure;

        if (!capacity_failure)
        {
            if (last_update_error != UPDATE_OK)
            {
                if (last_update_error == UPDATE_INVALID_PERIODIC_GEOMETRY)
                {
                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorSimulationBreakDown,
                        "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
                        "Reason:\n\tneighbor-list build rejected invalid "
                        "periodic geometry\n\t(step=%d, MPI/PP rank=%d/%d, "
                        "owned/ghost=%d/%d, grid_shape=(%d,%d,%d), "
                        "grid_count=%d)\n",
                        step, CONTROLLER::MPI_rank, CONTROLLER::PP_MPI_rank,
                        local_atom_numbers, ghost_numbers, grids.Nx, grids.Ny,
                        grids.Nz, grids.grid_numbers);
                }
                else
                {
                    UPDATE_ERROR_RECORD error_record = {};
                    deviceMemcpy(&error_record, d_update_error,
                                 sizeof(error_record),
                                 deviceMemcpyDeviceToHost);
                    const char* ownership =
                        error_record.header.local_atom >= 0 &&
                                error_record.header.local_atom <
                                    local_atom_numbers
                            ? "owned"
                            : (error_record.header.local_atom >=
                                           local_atom_numbers &&
                                       error_record.header.local_atom <
                                           coordinate_numbers
                                   ? "ghost"
                                   : "invalid");
                    const char* failure_stage = "unknown update error";
                    if (last_update_error == UPDATE_RAW_COORDINATE_NONFINITE)
                        failure_stage = "raw coordinate non-finite";
                    else if (last_update_error ==
                             UPDATE_WRAPPED_COORDINATE_NONFINITE)
                        failure_stage = "wrapped coordinate non-finite";
                    else if (last_update_error == UPDATE_GRID_INDEX_INVALID)
                        failure_stage = "grid index invalid";
                    const char* wrapped_available =
                        last_update_error == UPDATE_RAW_COORDINATE_NONFINITE
                            ? "no"
                            : "yes";
                    const char* grid_available =
                        last_update_error == UPDATE_GRID_INDEX_INVALID ? "yes"
                                                                       : "no";

                    controller->Throw_Formatted_SPONGE_Error(
                        spongeErrorSimulationBreakDown,
                        "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
                        "Reason:\n\tneighbor-list build rejected: %s\n"
                        "\t(error=%d, step=%d, MPI/PP rank=%d/%d, "
                        "local/global atom (0-based)=%d/%d, ownership=%s, "
                        "owned/ghost=%d/%d, raw_bits=(0x%08x,0x%08x,"
                        "0x%08x), wrapped_available=%s, "
                        "wrapped_bits=(0x%08x,0x%08x,0x%08x), "
                        "grid_available=%s, grid_index=(%d,%d,%d), "
                        "grid_id=%d, grid_shape=(%d,%d,%d), "
                        "grid_count=%d)\n",
                        failure_stage, last_update_error, step,
                        CONTROLLER::MPI_rank, CONTROLLER::PP_MPI_rank,
                        error_record.header.local_atom,
                        error_record.global_atom, ownership, local_atom_numbers,
                        ghost_numbers, error_record.raw_bits.uint_x,
                        error_record.raw_bits.uint_y,
                        error_record.raw_bits.uint_z, wrapped_available,
                        error_record.wrapped_bits.uint_x,
                        error_record.wrapped_bits.uint_y,
                        error_record.wrapped_bits.uint_z, grid_available,
                        error_record.grid_index.int_x,
                        error_record.grid_index.int_y,
                        error_record.grid_index.int_z, error_record.grid_id,
                        grids.Nx, grids.Ny, grids.Nz, grids.grid_numbers);
                }
            }
            else if (is_needed_full)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
                    "Reason:\n\tfull neighbor-list build failed: %s "
                    "(atom=%d, value=%d, owned=%d, coordinates=%d, "
                    "grid-required=%d/%d, half-required=%d)\n",
                    full_neighbor_list.Last_Error_Message(),
                    full_neighbor_list.last_error_atom,
                    full_neighbor_list.last_error_value, local_atom_numbers,
                    coordinate_numbers, h_neighbor_grid_overflow,
                    h_neighbor_grid_ghost_overflow, h_neighbor_list_overflow);
            }
            else
            {
                controller->Throw_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
                    "Reason:\n\tneighbor-list update rejected its current "
                    "state\n");
            }
            return rebound;
        }
        if (throw_error_when_overflow)
        {
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow,
                "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
                "Reason:\n\tneighbor-list storage overflowed while building "
                "the current force-evaluation state\n");
            return rebound;
        }

        const auto required_capacity =
            [&](int current, int required, const char* name)
        {
            if (current <= 0 || required <= current)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorOverflow,
                    "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
                    "Reason:\n\t%s build reported required capacity %d from "
                    "current capacity %d; the rejected build cannot be "
                    "recovered safely\n",
                    name, required, current);
            }
            return required;
        };
        const int grown_grid_atoms =
            grid_atom_capacity_failure
                ? required_capacity(max_atom_in_grid_numbers,
                                    h_neighbor_grid_overflow, "grid-atom")
                : max_atom_in_grid_numbers;
        const int grown_grid_ghosts =
            grid_ghost_capacity_failure
                ? required_capacity(max_ghost_in_grid_numbers,
                                    h_neighbor_grid_ghost_overflow,
                                    "grid-ghost")
                : max_ghost_in_grid_numbers;
        const int required_neighbors =
            std::max(h_neighbor_list_overflow, required_full_neighbors);
        const int grown_neighbors =
            (half_capacity_failure || full_capacity_failure)
                ? required_capacity(max_neighbor_numbers, required_neighbors,
                                    "per-atom neighbor")
                : max_neighbor_numbers;

        controller->commands["neighbor_list_max_atom_in_grid_numbers"] =
            std::to_string(grown_grid_atoms);
        controller->commands["neighbor_list_max_ghost_in_grid_numbers"] =
            std::to_string(grown_grid_ghosts);
        controller->commands["neighbor_list_max_neighbor_numbers"] =
            std::to_string(grown_neighbors);
        controller->printf(
            "Neighbor-list capacity was insufficient for the current state; "
            "rebuilding exactly with grid atoms=%d, grid ghosts=%d, "
            "neighbors=%d (grid_atoms=%d, grid_ghosts=%d, half=%d, "
            "full=%d).\n",
            grown_grid_atoms, grown_grid_ghosts, grown_neighbors,
            h_neighbor_grid_overflow, h_neighbor_grid_ghost_overflow,
            h_neighbor_list_overflow, required_full_neighbors);

        const int configured_atom_capacity = atom_numbers;
        const float configured_cutoff = cutoff;
        const float configured_skin = skin;
        Clear();
        Initial(controller, configured_atom_capacity, configured_cutoff,
                configured_skin, cell, rcell);
        rebound = true;
    }
}

void NEIGHBOR_LIST::Clear()
{
    is_initialized = false;
    active_local_atom_numbers = 0;
    if (d_temp != NULL) Free_Single_Device_Pointer((void**)&d_temp);
    if (h_nl != NULL || d_nl != NULL)
        Free_Host_And_Device_Pointer((void**)&h_nl, (void**)&d_nl);
    if (d_neighbor_grid_overflow != NULL)
        Free_Host_And_Device_Pointer(NULL, (void**)&d_neighbor_grid_overflow);
    if (d_neighbor_list_overflow != NULL)
        Free_Host_And_Device_Pointer(NULL, (void**)&d_neighbor_list_overflow);
    if (d_neighbor_grid_ghost_overflow != NULL)
        Free_Host_And_Device_Pointer(NULL,
                                     (void**)&d_neighbor_grid_ghost_overflow);
    if (d_update_error != NULL)
        Free_Single_Device_Pointer((void**)&d_update_error);
    full_neighbor_list.Clear();
    grids.Clear();
    updator.Clear();
    atom_numbers = 0;
    ghost_numbers = 0;
    cutoff = 0.0f;
    skin = 0.0f;
    max_neighbor_numbers = 0;
    max_atom_in_grid_numbers = 0;
    max_ghost_in_grid_numbers = 0;
    h_neighbor_grid_overflow = 0;
    h_neighbor_grid_ghost_overflow = 0;
    h_neighbor_list_overflow = 0;
    last_update_error = UPDATE_OK;
    last_error_atom = -1;
    throw_error_when_overflow = false;
}
