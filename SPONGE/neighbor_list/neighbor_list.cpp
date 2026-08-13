#include "neighbor_list.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>

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

// 按 grid_j 升序逐个检查 27 种周期镜像位移，命中即记录，输出为升序邻居表。
static __host__ __device__ __forceinline__ int Find_Neighbor_Grids_Scan(
    int grid_i, int Nx, int Ny, int Nz, int64_t dxy, int64_t dxz, int64_t dyz,
    LTMatrix3 cell, int* neighbor_i)
{
    const int64_t NyNz = static_cast<int64_t>(Ny) * Nz;
    int local = 0;
    int dx, dy, dz;
    const int64_t grid_i_x = grid_i / NyNz;
    const int64_t grid_i_y = (grid_i % NyNz) / Nz;
    const int64_t grid_i_z = grid_i % Nz;
    int cx, cy;
    int grid_j = 0;
    for (int64_t grid_j_x = 0; grid_j_x < Nx; grid_j_x += 1)
    {
        const int64_t base_x = grid_i_x - grid_j_x;
        for (int64_t grid_j_y = 0; grid_j_y < Ny; grid_j_y += 1)
        {
            const int64_t base_y = grid_i_y - grid_j_y;
            for (int64_t grid_j_z = 0; grid_j_z < Nz;
                 grid_j_z += 1, grid_j += 1)
            {
                const int64_t base_z = grid_i_z - grid_j_z;
                for (int k = 0; k < 27; k += 1)
                {
                    dx = k / 9 - 1;
                    dy = (k % 9) / 3 - 1;
                    dz = k % 3 - 1;
                    const int64_t temp_x = base_x +
                                           static_cast<int64_t>(dx) * Nx +
                                           dy * dxy + dz * dxz;
                    const int64_t temp_y = base_y +
                                           static_cast<int64_t>(dy) * Ny +
                                           dz * dyz;
                    const int64_t temp_z = base_z + static_cast<int64_t>(dz) * Nz;
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
            if (local < 0) break;
        }
        if (local < 0) break;
    }
    return local;
}

#ifdef GPU_ARCH_NAME
// 枚举路径单线程候选上限；超出即回退到全扫描，语义不变
#define GRID_ENUM_CAPACITY 2048
#endif

static __global__ void Find_Neighor_Grids_Device(
    int grid_numbers, int* neighbor_grid_numbers, int* neighbor_grids, int Nx,
    int Ny, int Nz, int64_t dxy, int64_t dxz, int64_t dyz, LTMatrix3 cell)
{
    SIMPLE_DEVICE_FOR(grid_i, grid_numbers)
    {
        int* neighbor_i = neighbor_grids + (size_t)MAX_GRID_NEIGHBORS * grid_i;
        int local = 0;
#ifdef GPU_ARCH_NAME
        // 全扫描是 O(grid_numbers^2)。命中条件对每种周期镜像 (dx,dy,dz)
        // 都是关于 grid_j 三个分量的区间约束，因此可以按镜像直接枚举候选
        // 小盒，再排序去重得到与全扫描逐位一致的升序邻居表。
        const int64_t grid_i_x = grid_i / (static_cast<int64_t>(Ny) * Nz);
        const int64_t grid_i_y =
            (grid_i % (static_cast<int64_t>(Ny) * Nz)) / Nz;
        const int64_t grid_i_z = grid_i % Nz;
        int candidates[GRID_ENUM_CAPACITY];
        int candidate_numbers = 0;
        bool enumeration_ok = true;
        for (int k = 0; k < 27 && enumeration_ok; k += 1)
        {
            const int dx = k / 9 - 1;
            const int dy = (k % 9) / 3 - 1;
            const int dz = k % 3 - 1;
            const int cx = static_cast<double>(dy) * cell.a21 +
                                       static_cast<double>(dz) * cell.a31 ==
                                   0.0
                               ? -2
                               : -3;
            const int cy =
                static_cast<double>(dz) * cell.a32 == 0.0 ? -2 : -3;
            const int64_t shift_x =
                static_cast<int64_t>(dx) * Nx + dy * dxy + dz * dxz;
            const int64_t shift_y = static_cast<int64_t>(dy) * Ny + dz * dyz;
            const int64_t shift_z = static_cast<int64_t>(dz) * Nz;
            const int64_t lo_x = grid_i_x + shift_x - 2 > 0
                                     ? grid_i_x + shift_x - 2
                                     : 0;
            const int64_t hi_x = grid_i_x + shift_x - cx < Nx - 1
                                     ? grid_i_x + shift_x - cx
                                     : Nx - 1;
            const int64_t lo_y = grid_i_y + shift_y - 2 > 0
                                     ? grid_i_y + shift_y - 2
                                     : 0;
            const int64_t hi_y = grid_i_y + shift_y - cy < Ny - 1
                                     ? grid_i_y + shift_y - cy
                                     : Ny - 1;
            const int64_t lo_z =
                grid_i_z + shift_z - 2 > 0 ? grid_i_z + shift_z - 2 : 0;
            const int64_t hi_z =
                grid_i_z + shift_z + 2 < Nz - 1 ? grid_i_z + shift_z + 2
                                                : Nz - 1;
            if (lo_x > hi_x || lo_y > hi_y || lo_z > hi_z) continue;
            const int64_t box =
                (hi_x - lo_x + 1) * (hi_y - lo_y + 1) * (hi_z - lo_z + 1);
            if (candidate_numbers + box > GRID_ENUM_CAPACITY)
            {
                enumeration_ok = false;
                break;
            }
            for (int64_t jx = lo_x; jx <= hi_x; jx += 1)
            {
                for (int64_t jy = lo_y; jy <= hi_y; jy += 1)
                {
                    for (int64_t jz = lo_z; jz <= hi_z; jz += 1)
                    {
                        candidates[candidate_numbers] =
                            static_cast<int>((jx * Ny + jy) * Nz + jz);
                        candidate_numbers += 1;
                    }
                }
            }
        }
        if (enumeration_ok)
        {
            for (int gap = candidate_numbers / 2; gap > 0; gap /= 2)
            {
                for (int i = gap; i < candidate_numbers; i += 1)
                {
                    const int value = candidates[i];
                    int j = i;
                    for (; j >= gap && candidates[j - gap] > value; j -= gap)
                    {
                        candidates[j] = candidates[j - gap];
                    }
                    candidates[j] = value;
                }
            }
            for (int i = 0; i < candidate_numbers; i += 1)
            {
                if (i > 0 && candidates[i] == candidates[i - 1]) continue;
                if (local >= MAX_GRID_NEIGHBORS)
                {
                    local = -1;
                    break;
                }
                neighbor_i[local++] = candidates[i];
            }
        }
        else
        {
            local = Find_Neighbor_Grids_Scan(grid_i, Nx, Ny, Nz, dxy, dxz, dyz,
                                             cell, neighbor_i);
        }
#else
        local = Find_Neighbor_Grids_Scan(grid_i, Nx, Ny, Nz, dxy, dxz, dyz,
                                         cell, neighbor_i);
#endif
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
    size_t prefix_slots = 0, prefix_bytes = 0;
    size_t grid_int_bytes = 0, neighbor_bytes = 0, atom_bytes = 0;
    size_t ghost_bytes = 0, atom_crd_bytes = 0, ghost_crd_bytes = 0;
    if (!Checked_Product(grid_count, (size_t)MAX_GRID_NEIGHBORS,
                         &neighbor_slots) ||
        !Checked_Product(grid_count, (size_t)(MAX_GRID_NEIGHBORS + 1),
                         &prefix_slots) ||
        !Checked_Bytes(prefix_slots, sizeof(int), &prefix_bytes) ||
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
                                       h_grid_atom_numbers, grid_int_bytes) ||
        !Device_Malloc_Safely((void**)&d_grid_neighbor_prefix, prefix_bytes))
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
    if (d_grid_neighbor_prefix != NULL)
        Free_Single_Device_Pointer((void**)&d_grid_neighbor_prefix);
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

static __device__ __forceinline__ void Record_Neighbor_Update_Error(
    int* update_error, int error, int atom)
{
#ifdef GPU_ARCH_NAME
    if (atomicCAS(update_error, NEIGHBOR_LIST::UPDATE_OK, error) ==
        NEIGHBOR_LIST::UPDATE_OK)
        update_error[1] = atom;
#else
#pragma omp critical(sponge_neighbor_list_update_error)
    {
        if (update_error[0] == NEIGHBOR_LIST::UPDATE_OK)
        {
            update_error[0] = error;
            update_error[1] = atom;
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
    int* update_error)
{
    if (need[0] == 0) return;
    SIMPLE_DEVICE_FOR(tid, atom_numbers + ghost_numbers)
    {
        VECTOR local_crd = crd[tid];
        bool coordinate_is_finite = isfinite(local_crd.x) &&
                                    isfinite(local_crd.y) &&
                                    isfinite(local_crd.z);
        if (coordinate_is_finite && need_copy && tid < atom_numbers)
        {
            old_crd[tid] = local_crd;
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
            coordinate_is_finite = isfinite(local_crd.x) &&
                                   isfinite(local_crd.y) &&
                                   isfinite(local_crd.z);
        }
        if (!coordinate_is_finite)
        {
            Record_Neighbor_Update_Error(
                update_error, NEIGHBOR_LIST::UPDATE_INVALID_GEOMETRY, tid);
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
                    update_error, NEIGHBOR_LIST::UPDATE_INVALID_GEOMETRY, tid);
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
                    grid_crd[slot] = crd[tid];
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
                    grid_ghost_crd[slot] = crd[tid];
                }
            }
        }
    }
}

#ifdef USE_GPU

// 为每个 grid_i 计算其邻居格原子数的独占前缀和（行长 MAX_GRID_NEIGHBORS+1）。
// 邻居格拓扑在盒子不变时是静态的，但每格原子数每次重建都会变，
// 因此前缀和在每次重建时由本 kernel 重算。早退条件与
// Find_Neighbors_Gridly 一致：不重建或构建已失败时不做任何事。
static __global__ void Build_Grid_Neighbor_Prefix(
    const int* need, int grid_numbers, const int* grid_neighbor_numbers,
    const int* grid_neighbors, const int* grid_atom_numbers,
    int* grid_neighbor_prefix, const int* neighbor_grid_overflow,
    const int* neighbor_grid_ghost_overflow, const int* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error[0] != NEIGHBOR_LIST::UPDATE_OK)
        return;
    SIMPLE_DEVICE_FOR(grid_i, grid_numbers)
    {
        const int neighbor_count = grid_neighbor_numbers[grid_i];
        const int* neighbor_row =
            grid_neighbors + (size_t)grid_i * MAX_GRID_NEIGHBORS;
        int* prefix_row =
            grid_neighbor_prefix + (size_t)grid_i * (MAX_GRID_NEIGHBORS + 1);
        int sum = 0;
        prefix_row[0] = 0;
        for (int k = 0; k < neighbor_count; ++k)
        {
            sum += grid_atom_numbers[neighbor_row[k]];
            prefix_row[k + 1] = sum;
        }
    }
}

static __global__ void Find_Neighbors_Gridly(
    int* atom_local, int atom_numbers, const int* need, int grid_numbers,
    int* grid_neighbor_numbers, int* grid_neighbors, VECTOR* grid_crd,
    LTMatrix3 cell, LTMatrix3 rcell, int max_atom_numbers_in_grid,
    ATOM_GROUP* nl, float cutoff_skin_square, int* grid_atom_numbers,
    int* grid_atoms, int max_neighbor_numbers, int* neighbor_list_overflow,
    const int* grid_neighbor_prefix, VECTOR* grid_ghost_crd,
    int max_ghost_numbers_in_grid, int* grid_ghost_numbers, int* grid_ghosts,
    const int* neighbor_grid_overflow, const int* neighbor_grid_ghost_overflow,
    const int* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error[0] != NEIGHBOR_LIST::UPDATE_OK)
        return;
    extern __shared__ unsigned char shared_mem[];
    VECTOR* sh_crd = reinterpret_cast<VECTOR*>(shared_mem);
    int* sh_atoms = reinterpret_cast<int*>(sh_crd + max_atom_numbers_in_grid);
    int* sh_globals = sh_atoms + max_atom_numbers_in_grid;
    int* sh_prefix = sh_globals + max_atom_numbers_in_grid;

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

        // 预载 grid_i 的原子 id、全局 id（旧实现内层循环对每个
        // (j-chunk, i) 重复从全局内存 gather atom_local，这里顺便预载
        // 到 shared）与坐标
        for (int idx = threadIdx.x; idx < atom_numbers_in_grid_i;
             idx += blockDim.x)
        {
            const int atom_i = bucket_i[idx];
            sh_atoms[idx] = atom_i;
            sh_globals[idx] = atom_local[atom_i];
            sh_crd[idx] = grid_crd_i[idx];
        }

        const int neighbor_count = grid_neighbor_numbers[grid_i];
        const int* prefix_row =
            grid_neighbor_prefix + (size_t)grid_i * (MAX_GRID_NEIGHBORS + 1);
        for (int idx = threadIdx.x; idx <= neighbor_count; idx += blockDim.x)
        {
            sh_prefix[idx] = prefix_row[idx];
        }
        __syncthreads();

        if (neighbor_count == 0)
        {
            __syncthreads();
            continue;
        }

        // 把该 grid_i 全部邻居格的原子看成展平的一维候选空间（总长
        // sh_prefix[neighbor_count]）。warp 以 lane_stride 个连续展平 j
        // 为一组认领，lane 用 shared 前缀和二分定位 (grid_j, 局部 j)，
        // lane 利用率不再受单格平均原子数（~12 < 32）的限制。
        const int total_j = sh_prefix[neighbor_count];
        const int* neighbor_row =
            grid_neighbors + (size_t)grid_i * MAX_GRID_NEIGHBORS;
        for (int j_base = warp_id * lane_stride; j_base < total_j;
             j_base += warps_per_block * lane_stride)
        {
            const int flat_j = j_base + lane_index;
            const bool active = flat_j < total_j;
            int atom_j = 0;
            int global_j = 0;
            VECTOR crd_j = {0, 0, 0};
            if (active)
            {
                // 二分：sh_prefix[1..neighbor_count] 中首个 > flat_j 的
                // 位置减一，即最大的满足 sh_prefix[k] <= flat_j 的 k
                // （空邻居格的前缀和与前一格相同，天然不会被命中）
                int lo = 1;
                int hi = neighbor_count;
                while (lo < hi)
                {
                    const int mid = (lo + hi) >> 1;
                    if (sh_prefix[mid] <= flat_j)
                        lo = mid + 1;
                    else
                        hi = mid;
                }
                const int k = lo - 1;
                const int grid_j = neighbor_row[k];
                const int local_j = flat_j - sh_prefix[k];
                const size_t slot_j =
                    (size_t)grid_j * max_atom_numbers_in_grid + local_j;
                atom_j = grid_atoms[slot_j];
                global_j = atom_local[atom_j];
                crd_j = grid_crd[slot_j];
            }

            for (int i = 0; i < atom_numbers_in_grid_i; ++i)
            {
                int atom_i = sh_atoms[i];
                int global_i = sh_globals[i];
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

// ---- LJ cluster-pair tile 表构建（S2，设计文档 §1.3；暂无消费者）----

// 单 block 扫描：为每格计算 local/ghost cluster 起始编号（独占前缀）与总数。
// ghost 段统一编号接在 local 段之后。早退条件与 Find_Neighbors_Gridly 一致。
static __global__ void Build_LJ_Cluster_Prefix(
    const int* need, int grid_numbers, const int* grid_atom_numbers,
    const int* grid_ghost_numbers, int* grid_cluster_base,
    int* grid_ghost_cluster_base, int* cluster_numbers,
    const int* neighbor_grid_overflow, const int* neighbor_grid_ghost_overflow,
    const int* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error[0] != NEIGHBOR_LIST::UPDATE_OK)
        return;
    extern __shared__ int sh_scan[];
    __shared__ int sh_carry;
    if (threadIdx.x == 0) sh_carry = 0;
    __syncthreads();
    for (int segment = 0; segment < 2; ++segment)
    {
        const int* counts = segment == 0 ? grid_atom_numbers : grid_ghost_numbers;
        int* bases = segment == 0 ? grid_cluster_base : grid_ghost_cluster_base;
        int segment_total = 0;
        for (int base = 0; base < grid_numbers; base += blockDim.x)
        {
            const int g = base + threadIdx.x;
            const int v = g < grid_numbers
                              ? (counts[g] + LJ_TILE_CLUSTER_SIZE - 1) /
                                    LJ_TILE_CLUSTER_SIZE
                              : 0;
            sh_scan[threadIdx.x] = v;
            __syncthreads();
            for (int offset = 1; offset < blockDim.x; offset <<= 1)
            {
                const int t =
                    threadIdx.x >= offset ? sh_scan[threadIdx.x - offset] : 0;
                __syncthreads();
                sh_scan[threadIdx.x] += t;
                __syncthreads();
            }
            if (g < grid_numbers)
                bases[g] = sh_carry + sh_scan[threadIdx.x] - v;
            __syncthreads();
            if (threadIdx.x == 0) sh_carry += sh_scan[blockDim.x - 1];
            __syncthreads();
        }
        segment_total = sh_carry;
        if (threadIdx.x == 0)
        {
            cluster_numbers[segment] =
                segment == 0 ? segment_total : segment_total - cluster_numbers[0];
        }
        __syncthreads();
    }
}

// 把格桶顺序的原子按 8 个一组填进 cluster 表（local 段 + ghost 段），
// 尾部 padding 写 -1，并写 ghost 标志。
static __global__ void Fill_LJ_Cluster_Atoms(
    const int* need, int grid_numbers, const int* grid_atom_numbers,
    const int* grid_atoms, const int* grid_cluster_base,
    int max_atom_in_grid, const int* grid_ghost_numbers, const int* grid_ghosts,
    const int* grid_ghost_cluster_base, int max_ghost_in_grid,
    int* cluster_atoms, int* cluster_flags, const int* neighbor_grid_overflow,
    const int* neighbor_grid_ghost_overflow, const int* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error[0] != NEIGHBOR_LIST::UPDATE_OK)
        return;
    SIMPLE_DEVICE_FOR(grid, grid_numbers)
    {
        const int count = grid_atom_numbers[grid];
        const int base = grid_cluster_base[grid];
        const int nc = (count + LJ_TILE_CLUSTER_SIZE - 1) / LJ_TILE_CLUSTER_SIZE;
        const int* bucket = grid_atoms + (size_t)grid * max_atom_in_grid;
        for (int k = 0; k < nc * LJ_TILE_CLUSTER_SIZE; ++k)
        {
            cluster_atoms[base * LJ_TILE_CLUSTER_SIZE + k] =
                k < count ? bucket[k] : -1;
        }
        for (int c = 0; c < nc; ++c) cluster_flags[base + c] = 0;

        const int ghost_count = grid_ghost_numbers[grid];
        const int ghost_base = grid_ghost_cluster_base[grid];
        const int gnc =
            (ghost_count + LJ_TILE_CLUSTER_SIZE - 1) / LJ_TILE_CLUSTER_SIZE;
        const int* ghost_bucket = grid_ghosts + (size_t)grid * max_ghost_in_grid;
        for (int k = 0; k < gnc * LJ_TILE_CLUSTER_SIZE; ++k)
        {
            cluster_atoms[ghost_base * LJ_TILE_CLUSTER_SIZE + k] =
                k < ghost_count ? ghost_bucket[k] : -1;
        }
        for (int c = 0; c < gnc; ++c) cluster_flags[ghost_base + c] = 1;
    }
}

// 排除表查询：调用方预载 (list_start, n, min, max)，先 min/max 剪枝再
// 二分。语义与 Delete_Excluded_Atoms_Serial_In_Neighbor_List 的单向
// 检查逐语句对应（排除表每原子段升序）。
static __device__ __forceinline__ bool LJ_Tile_Excluded_Range_Contains(
    LJ_EXCL_RANGE range, int global_other, const int* excluded_list)
{
    if (range.n <= 0 || global_other < range.min_atom ||
        global_other > range.max_atom)
        return false;
    int lo = range.list_start;
    int hi = range.list_start + range.n;
    while (lo < hi)
    {
        const int mid = (lo + hi) >> 1;
        if (excluded_list[mid] < global_other)
            lo = mid + 1;
        else
            hi = mid;
    }
    return lo < range.list_start + range.n &&
           excluded_list[lo] == global_other;
}

// 从排除表数组预载一个 global 原子的排除范围（无排除时返回 n=0）
static __device__ __forceinline__ LJ_EXCL_RANGE LJ_Tile_Load_Excluded_Range(
    int global_atom, const int* excluded_list_start, const int* excluded_list,
    const int* excluded_numbers)
{
    LJ_EXCL_RANGE range = {0, 0, 0, -1};
    if (excluded_numbers == NULL || excluded_list == NULL) return range;
    const int n = excluded_numbers[global_atom];
    if (n > 0)
    {
        const int start = excluded_list_start[global_atom];
        range.list_start = start;
        range.n = n;
        range.min_atom = excluded_list[start];
        range.max_atom = excluded_list[start + n - 1];
    }
    return range;
}

// 一次性预计算全部 global 原子的排除范围缓存（排除表运行期静态）
static __global__ void Build_LJ_Excluded_Range(
    int atom_numbers, const int* excluded_list_start,
    const int* excluded_list, const int* excluded_numbers,
    LJ_EXCL_RANGE* excl_range)
{
    SIMPLE_DEVICE_FOR(g, atom_numbers)
    {
        excl_range[g] = LJ_Tile_Load_Excluded_Range(g, excluded_list_start,
                                                    excluded_list,
                                                    excluded_numbers);
    }
}

// 对格对 (grid_i, grid_j) 展开 8x8 tile：
// - local×local 段只保留 grid_j >= grid_i（格对去重）；同格内 cluster 对
//   只保留 ci <= cj；自对 tile 靠 global_j > global_i 出下半三角。
//   注意 ci != cj 时无序对 {x,y} 在整个 tile 表中只出现一次（格对去重
//   保证），因此不再需要 pair 级方向判据——与半表"每对只出现一次"等价。
// - local×ghost 段不判方向、不做格对去重（ghost 只在它所在的格桶里出现
//   一次，天然无重复），与 Find_Neighbors_Gridly ghost 段一致。
// - 排除剔除：local×local 查 global 较小者（半表行主）的排除表；
//   local×ghost 正向查 i 的排除表、反向查 j 的排除表（镜像逻辑）。
//   排除检查按需做——只对距离存活的 pair 触发 gather（min/max 剪枝 +
//   罕发二分），不占用预载与 j-cluster 载入的关键路径。
// 任务 = (邻居格, j-cluster)，按 shared 前缀和展平后由 warp 认领，下一
// 任务的 j 数据软件预取（与当前任务评估重叠）；warp 内 lane 分摊 64 个
// pair（lane 的 j_local = lane&7 恒定，j 端数据 lanes 0-7 载入后 shuffle
// 广播），两轮 32-lane ballot 聚成 64 位掩码，非零才认领槽位写出
//（块级 shared 缓冲聚合，每 grid_i 一次全局 atomicAdd）；溢出走
// Record_Required_Capacity。
// 块级 tile 发射缓冲容量（1024×16B=16KB shared；典型 grid_i 产出
// ~300-500 tile，溢出时回退到逐 tile 直接全局认领，正确性不受影响）
#define LJ_TILE_BLOCK_BUFFER 1024

template <int ABLATE>
static __global__ void Build_LJ_Tile_List(
    const int* need, int grid_numbers, const int* grid_neighbor_numbers,
    const int* grid_neighbors, const VECTOR* grid_crd, LTMatrix3 cell,
    LTMatrix3 rcell, int max_atom_in_grid, const int* grid_atom_numbers,
    const int* grid_atoms, const int* grid_cluster_base,
    const VECTOR* grid_ghost_crd, int max_ghost_in_grid,
    const int* grid_ghost_numbers, const int* grid_ghosts,
    const int* grid_ghost_cluster_base, const int* atom_local,
    float cutoff_skin_square, const int* excluded_list,
    const LJ_EXCL_RANGE* excl_range, LJ_TILE* tiles, int tile_capacity,
    int* tile_count, int* tile_overflow,
    const int* neighbor_grid_overflow, const int* neighbor_grid_ghost_overflow,
    const int* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error[0] != NEIGHBOR_LIST::UPDATE_OK)
        return;
    extern __shared__ unsigned char shared_mem[];
    float4* sh_crd = reinterpret_cast<float4*>(shared_mem);
    // i 端排除范围 (list_start, n, min, max)，格级预载避免逐对依赖 gather
    LJ_EXCL_RANGE* sh_excl =
        reinterpret_cast<LJ_EXCL_RANGE*>(sh_crd + max_atom_in_grid);
    int* sh_globals = reinterpret_cast<int*>(sh_excl + max_atom_in_grid);
    // 邻居格任务数（j-cluster 对数）扫描缓冲，inclusive 前缀和
    int* sh_scan = reinterpret_cast<int*>(sh_globals + max_atom_in_grid);
    // 块级 tile 发射缓冲：块内先攒进 shared（shared 原子加），每个 grid_i
    // 结束只做一次全局 atomicAdd 认领整段槽位再协作写出——逐 tile 全局
    // 原子加在同一计数器上串行，实测是构建 kernel 的主要开销
    LJ_TILE* sh_tiles =
        reinterpret_cast<LJ_TILE*>(sh_scan + blockDim.x);
    __shared__ int sh_tile_fill;
    __shared__ int sh_tile_base;

    const int lane = threadIdx.x & (warpSize - 1);
    int warps_per_block = blockDim.x / warpSize;
    if (warps_per_block == 0)
    {
        warps_per_block = 1;
    }
    const int warp_id = threadIdx.x / warpSize;

    device_mask_t warp_mask = FULL_MASK;
    if (blockDim.x < warpSize)
    {
        warp_mask = deviceLowerLaneMask(blockDim.x);
    }

    for (int grid_i = blockIdx.x; grid_i < grid_numbers; grid_i += gridDim.x)
    {
        const int count_i = grid_atom_numbers[grid_i];
        if (count_i == 0)
        {
            __syncthreads();
            continue;
        }
        if (threadIdx.x == 0) sh_tile_fill = 0;
        // 预载 i 端 global id、坐标（float4 对齐，单发 16B shared 读）与
        // 排除范围
        for (int idx = threadIdx.x; idx < count_i; idx += blockDim.x)
        {
            const int atom_i =
                grid_atoms[(size_t)grid_i * max_atom_in_grid + idx];
            const int global_i = atom_local[atom_i];
            sh_globals[idx] = global_i;
            const VECTOR c = grid_crd[(size_t)grid_i * max_atom_in_grid + idx];
            sh_crd[idx] = make_float4(c.x, c.y, c.z, 0.0f);
            sh_excl[idx] =
                ABLATE & 1 ? LJ_EXCL_RANGE{0, 0, 0, -1} : excl_range[global_i];
        }

        const int neighbor_count = grid_neighbor_numbers[grid_i];
        const int* neighbor_row =
            grid_neighbors + (size_t)grid_i * MAX_GRID_NEIGHBORS;
        const int nci =
            (count_i + LJ_TILE_CLUSTER_SIZE - 1) / LJ_TILE_CLUSTER_SIZE;
        const int base_ci = grid_cluster_base[grid_i];

        // local / ghost 两段合并成一次扫描：任务 = (邻居格, j-cluster 对)，
        // 每个任务覆盖相邻两个 j-cluster，i 端 shared 读与任务管理开销
        // 减半，且每个 warp 同时评估 4 条独立距离链（2 cluster × 2 半区）
        // 提高 ILP。local 段只收 grid_j >= grid_i（格对去重）；ghost 段
        // 全收（ghost 只在它所在的格桶里出现一次，天然无重复）
        {
            int v = 0;
            if (threadIdx.x < (unsigned int)neighbor_count)
            {
                const int grid_j = neighbor_row[threadIdx.x];
                if (grid_j >= grid_i)
                    v += ((grid_atom_numbers[grid_j] +
                           LJ_TILE_CLUSTER_SIZE - 1) /
                              LJ_TILE_CLUSTER_SIZE +
                          1) >>
                         1;
                v += ((grid_ghost_numbers[grid_j] + LJ_TILE_CLUSTER_SIZE - 1) /
                          LJ_TILE_CLUSTER_SIZE +
                      1) >>
                     1;
            }
            sh_scan[threadIdx.x] = v;
        }
        __syncthreads();
        for (int offset = 1; offset < blockDim.x; offset <<= 1)
        {
            const int t = threadIdx.x >= (unsigned int)offset
                              ? sh_scan[threadIdx.x - offset]
                              : 0;
            __syncthreads();
            sh_scan[threadIdx.x] += t;
            __syncthreads();
        }
        const int total_tasks =
            neighbor_count > 0 ? sh_scan[neighbor_count - 1] : 0;

            // 任务循环（带下一任务 j-cluster 数据的软件预取：进入迭代时
            // 当前任务的 j 数据已在寄存器，评估当前任务前先发射下一任务的
            // 载入，隐藏依赖 gather 的延迟）
            auto load_task = [&](int t, int& out_grid_j, int& out_segment,
                                 int& out_cj, int& out_base_cj, int& out_rem_a,
                                 int& out_rem_b, int& out_j0_global,
                                 float& out_j0x, float& out_j0y,
                                 float& out_j0z, LJ_EXCL_RANGE& out_j0_range,
                                 int& out_j1_global, float& out_j1x,
                                 float& out_j1y, float& out_j1z,
                                 LJ_EXCL_RANGE& out_j1_range)
            {
                // 二分：首个 sh_scan[jj] > t 的邻居格
                int blo = 0;
                int bhi = neighbor_count - 1;
                while (blo < bhi)
                {
                    const int mid = (blo + bhi) >> 1;
                    if (sh_scan[mid] <= t)
                        blo = mid + 1;
                    else
                        bhi = mid;
                }
                out_grid_j = neighbor_row[blo];
                const int t_in = t - (blo == 0 ? 0 : sh_scan[blo - 1]);
                const int ncj_local =
                    (grid_atom_numbers[out_grid_j] + LJ_TILE_CLUSTER_SIZE - 1) /
                    LJ_TILE_CLUSTER_SIZE;
                const int local_tasks =
                    out_grid_j >= grid_i ? (ncj_local + 1) >> 1 : 0;
                const int* seg_atoms;
                const VECTOR* seg_crd;
                int seg_max, seg_count;
                if (t_in < local_tasks)
                {
                    out_segment = 0;
                    out_cj = t_in << 1;
                    seg_atoms = grid_atoms;
                    seg_crd = grid_crd;
                    seg_max = max_atom_in_grid;
                    seg_count = grid_atom_numbers[out_grid_j];
                    out_base_cj = grid_cluster_base[out_grid_j];
                }
                else
                {
                    out_segment = 1;
                    out_cj = (t_in - local_tasks) << 1;
                    seg_atoms = grid_ghosts;
                    seg_crd = grid_ghost_crd;
                    seg_max = max_ghost_in_grid;
                    seg_count = grid_ghost_numbers[out_grid_j];
                    out_base_cj = grid_ghost_cluster_base[out_grid_j];
                }
                out_rem_a = seg_count - (out_cj << 3);
                out_rem_b = out_rem_a - LJ_TILE_CLUSTER_SIZE;
                out_j0_global = -1;
                out_j1_global = -1;
                out_j0x = 0.0f;
                out_j0y = 0.0f;
                out_j0z = 0.0f;
                out_j1x = 0.0f;
                out_j1y = 0.0f;
                out_j1z = 0.0f;
                out_j0_range = LJ_EXCL_RANGE{0, 0, 0, -1};
                out_j1_range = LJ_EXCL_RANGE{0, 0, 0, -1};
                if (!(ABLATE & 16) && lane < LJ_TILE_CLUSTER_SIZE)
                {
                    const size_t slot = (size_t)out_grid_j * seg_max +
                                        (size_t)(out_cj << 3) + lane;
                    if (lane < out_rem_a)
                    {
                        out_j0_global = atom_local[seg_atoms[slot]];
                        const VECTOR c = seg_crd[slot];
                        out_j0x = c.x;
                        out_j0y = c.y;
                        out_j0z = c.z;
                        if (!(ABLATE & 1))
                            out_j0_range = excl_range[out_j0_global];
                    }
                    if (lane < out_rem_b)
                    {
                        out_j1_global =
                            atom_local[seg_atoms[slot + LJ_TILE_CLUSTER_SIZE]];
                        const VECTOR c = seg_crd[slot + LJ_TILE_CLUSTER_SIZE];
                        out_j1x = c.x;
                        out_j1y = c.y;
                        out_j1z = c.z;
                        if (!(ABLATE & 1))
                            out_j1_range = excl_range[out_j1_global];
                    }
                }
            };
            for (int task = warp_id; !(ABLATE & 4) && task < total_tasks;
                 task += warps_per_block)
            {
                int cur_grid_j = -1, cur_segment = 0, cur_cj = 0,
                    cur_base_cj = 0, cur_rem_a = 0, cur_rem_b = 0;
                int cur_j0_global = -1, cur_j1_global = -1;
                float cur_j0x = 0.0f, cur_j0y = 0.0f, cur_j0z = 0.0f;
                float cur_j1x = 0.0f, cur_j1y = 0.0f, cur_j1z = 0.0f;
                LJ_EXCL_RANGE cur_j0_range = {0, 0, 0, -1};
                LJ_EXCL_RANGE cur_j1_range = {0, 0, 0, -1};
                load_task(task, cur_grid_j, cur_segment, cur_cj, cur_base_cj,
                          cur_rem_a, cur_rem_b, cur_j0_global, cur_j0x,
                          cur_j0y, cur_j0z, cur_j0_range, cur_j1_global,
                          cur_j1x, cur_j1y, cur_j1z, cur_j1_range);

                // lane 的 j_local = lane&7 在两个 j-cluster 间复用
                const int jl = lane & (LJ_TILE_CLUSTER_SIZE - 1);
                const int my_j0_global =
                    deviceShfl(FULL_MASK, cur_j0_global, jl, 32);
                const float my_j0x = deviceShfl(FULL_MASK, cur_j0x, jl, 32);
                const float my_j0y = deviceShfl(FULL_MASK, cur_j0y, jl, 32);
                const float my_j0z = deviceShfl(FULL_MASK, cur_j0z, jl, 32);
                const int my_j1_global =
                    deviceShfl(FULL_MASK, cur_j1_global, jl, 32);
                const float my_j1x = deviceShfl(FULL_MASK, cur_j1x, jl, 32);
                const float my_j1y = deviceShfl(FULL_MASK, cur_j1y, jl, 32);
                const float my_j1z = deviceShfl(FULL_MASK, cur_j1z, jl, 32);
                const LJ_EXCL_RANGE my_j0_range = {
                    deviceShfl(FULL_MASK, cur_j0_range.list_start, jl, 32),
                    deviceShfl(FULL_MASK, cur_j0_range.n, jl, 32),
                    deviceShfl(FULL_MASK, cur_j0_range.min_atom, jl, 32),
                    deviceShfl(FULL_MASK, cur_j0_range.max_atom, jl, 32)};
                const LJ_EXCL_RANGE my_j1_range = {
                    deviceShfl(FULL_MASK, cur_j1_range.list_start, jl, 32),
                    deviceShfl(FULL_MASK, cur_j1_range.n, jl, 32),
                    deviceShfl(FULL_MASK, cur_j1_range.min_atom, jl, 32),
                    deviceShfl(FULL_MASK, cur_j1_range.max_atom, jl, 32)};
                const int base_cj = cur_base_cj;
                const bool jv_a = jl < cur_rem_a;
                const bool jv_b = jl < cur_rem_b;

                if (ABLATE & 8)
                {
                    // 消融：只跑任务管理 + j 载入；假使用防死代码消除
                    if (lane == 0 && my_j0x == 1.0e30f &&
                        my_j0_global == 0x7fffffff)
                        tiles[0].mask = 0;
                }
                else
                {
                    // local 段同格去重：cluster A 要求 ci <= cj、cluster B
                    //（cj+1）要求 ci <= cj+1；ci_end 取两者上界，A 越界在
                    // 评估时清掉；自对 tile 靠 global 判据出下半三角。
                    // ghost 段不判方向、不判大小
                    const bool same_grid =
                        cur_segment == 0 && cur_grid_j == grid_i;
                    const int ci_end =
                        same_grid ? min(cur_cj + 1, nci - 1) : nci - 1;
                    for (int ci = 0; ci <= ci_end; ++ci)
                    {
                        const int rem_i = count_i - (ci << 3);
                        const bool self_a = same_grid && ci == cur_cj;
                        const bool self_b = same_grid && ci == cur_cj + 1;
                        unsigned long long mask_a = 0, mask_b = 0;
                        for (int half = 0; half < 2; ++half)
                        {
                            // bit = i_local*8 + j_local；i_local = lane>>3
                            //（+4 第二轮），j_local = lane&7
                            const int il = (lane >> 3) + (half << 2);
                            const int ia = (ci << 3) + il;
                            bool pass_a = false, pass_b = false;
                            if (il < rem_i)
                            {
                                const float4 c4 = sh_crd[ia];
                                const VECTOR crd_i = {c4.x, c4.y, c4.z};
                                if (ABLATE & 2)
                                {
                                    pass_a = jv_a;
                                    pass_b = jv_b;
                                }
                                else
                                {
                                    if (jv_a)
                                    {
                                        const VECTOR crd_j = {my_j0x, my_j0y,
                                                              my_j0z};
                                        VECTOR dr =
                                            Get_Periodic_Displacement(
                                                crd_i, crd_j, cell, rcell);
                                        pass_a =
                                            dr * dr < cutoff_skin_square;
                                    }
                                    if (jv_b)
                                    {
                                        const VECTOR crd_j = {my_j1x, my_j1y,
                                                              my_j1z};
                                        VECTOR dr =
                                            Get_Periodic_Displacement(
                                                crd_i, crd_j, cell, rcell);
                                        pass_b =
                                            dr * dr < cutoff_skin_square;
                                    }
                                }
                                if (same_grid)
                                {
                                    if (ci > cur_cj)
                                        pass_a = false;
                                    else if (self_a && pass_a)
                                        pass_a =
                                            my_j0_global > sh_globals[ia];
                                    if (self_b && pass_b)
                                        pass_b =
                                            my_j1_global > sh_globals[ia];
                                }
                                // 排除检查：i 端范围在格级预载（shared），
                                // j 端范围随任务预取；逐对只做 min/max
                                // 剪枝 + 罕发二分，无依赖 gather
                                if (!(ABLATE & 1) && (pass_a || pass_b))
                                {
                                    const int global_i = sh_globals[ia];
                                    if (cur_segment == 0)
                                    {
                                        // local×local：查 global 较小者
                                        //（半表行主）的排除表
                                        if (pass_a)
                                        {
                                            pass_a =
                                                global_i < my_j0_global
                                                    ? !LJ_Tile_Excluded_Range_Contains(
                                                          sh_excl[ia],
                                                          my_j0_global,
                                                          excluded_list)
                                                    : !LJ_Tile_Excluded_Range_Contains(
                                                          my_j0_range,
                                                          global_i,
                                                          excluded_list);
                                        }
                                        if (pass_b)
                                        {
                                            pass_b =
                                                global_i < my_j1_global
                                                    ? !LJ_Tile_Excluded_Range_Contains(
                                                          sh_excl[ia],
                                                          my_j1_global,
                                                          excluded_list)
                                                    : !LJ_Tile_Excluded_Range_Contains(
                                                          my_j1_range,
                                                          global_i,
                                                          excluded_list);
                                        }
                                    }
                                    else
                                    {
                                        // local×ghost：正向查 i、反向查 j
                                        if (pass_a)
                                            pass_a =
                                                !LJ_Tile_Excluded_Range_Contains(
                                                    sh_excl[ia], my_j0_global,
                                                    excluded_list) &&
                                                !LJ_Tile_Excluded_Range_Contains(
                                                    my_j0_range, global_i,
                                                    excluded_list);
                                        if (pass_b)
                                            pass_b =
                                                !LJ_Tile_Excluded_Range_Contains(
                                                    sh_excl[ia], my_j1_global,
                                                    excluded_list) &&
                                                !LJ_Tile_Excluded_Range_Contains(
                                                    my_j1_range, global_i,
                                                    excluded_list);
                                    }
                                }
                            }
                            const unsigned int ballot_a =
                                LaneGroup::And(LaneGroup::Ballot(pass_a),
                                               LaneMask(warp_mask))
                                    .bits;
                            const unsigned int ballot_b =
                                LaneGroup::And(LaneGroup::Ballot(pass_b),
                                               LaneMask(warp_mask))
                                    .bits;
                            if (half == 0)
                            {
                                mask_a = ballot_a;
                                mask_b = ballot_b;
                            }
                            else
                            {
                                mask_a |= (unsigned long long)ballot_a << 32;
                                mask_b |= (unsigned long long)ballot_b << 32;
                            }
                        }
                        if (mask_a != 0 && lane == 0)
                        {
                            const int pos = atomicAdd(&sh_tile_fill, 1);
                            if (pos < LJ_TILE_BLOCK_BUFFER)
                            {
                                sh_tiles[pos].cluster_i = base_ci + ci;
                                sh_tiles[pos].cluster_j = base_cj + cur_cj;
                                sh_tiles[pos].mask = mask_a;
                            }
                            else
                            {
                                // 缓冲溢出（罕见）：回退逐 tile 直接全局认领
                                const int slot = atomicAdd(tile_count, 1);
                                if (slot < tile_capacity)
                                {
                                    tiles[slot].cluster_i = base_ci + ci;
                                    tiles[slot].cluster_j = base_cj + cur_cj;
                                    tiles[slot].mask = mask_a;
                                }
                                else
                                {
                                    Record_Required_Capacity(
                                        tile_overflow,
                                        Required_After_Reservation(slot, 1));
                                }
                            }
                        }
                        if (mask_b != 0 && lane == 0)
                        {
                            const int pos = atomicAdd(&sh_tile_fill, 1);
                            if (pos < LJ_TILE_BLOCK_BUFFER)
                            {
                                sh_tiles[pos].cluster_i = base_ci + ci;
                                sh_tiles[pos].cluster_j = base_cj + cur_cj + 1;
                                sh_tiles[pos].mask = mask_b;
                            }
                            else
                            {
                                const int slot = atomicAdd(tile_count, 1);
                                if (slot < tile_capacity)
                                {
                                    tiles[slot].cluster_i = base_ci + ci;
                                    tiles[slot].cluster_j =
                                        base_cj + cur_cj + 1;
                                    tiles[slot].mask = mask_b;
                                }
                                else
                                {
                                    Record_Required_Capacity(
                                        tile_overflow,
                                        Required_After_Reservation(slot, 1));
                                }
                            }
                        }
                    }
                }  // !(ABLATE & 8)
            }
            __syncthreads();  // 重写 sh_scan / 读 sh_tile_fill 前等所有 warp

        // 冲刷块级缓冲：一次全局 atomicAdd 认领整段槽位，协作写出
        const int buffered =
            sh_tile_fill < LJ_TILE_BLOCK_BUFFER ? sh_tile_fill
                                                : LJ_TILE_BLOCK_BUFFER;
        if (buffered > 0)
        {
            if (threadIdx.x == 0)
                sh_tile_base = atomicAdd(tile_count, buffered);
            __syncthreads();
            const int base = sh_tile_base;
            for (int k = threadIdx.x; k < buffered; k += blockDim.x)
            {
                if (base + k < tile_capacity) tiles[base + k] = sh_tiles[k];
            }
            if (base + buffered > tile_capacity && threadIdx.x == 0)
            {
                Record_Required_Capacity(
                    tile_overflow,
                    Required_After_Reservation(base, buffered));
            }
        }

        __syncthreads();
    }
}

// ---- S3：tile 按 cluster_i 计数排序（分组序下标表）----
// tile kernel 按行（同一 i-cluster 的连续 tile 段）消费时，f_i/能量/维里可
// 在寄存器累计、每段只归约写出一次，i 侧全局原子加降一个量级。行内顺序
// 任意（散射原子加序），消费端不依赖行界，靠 tile.cluster_i 变化检测分段。
static __global__ void Zero_LJ_Tile_Rows(const int* need, int* row_cursor,
                                         const int cluster_capacity,
                                         const int* update_error)
{
    if (need[0] == 0 || update_error[0] != NEIGHBOR_LIST::UPDATE_OK) return;
    SIMPLE_DEVICE_FOR(i, cluster_capacity)
    {
        row_cursor[i] = 0;
    }
}

static __global__ void Count_LJ_Tile_Rows(const int* need, const LJ_TILE* tiles,
                                          const int* tile_count,
                                          int* row_cursor,
                                          const int* update_error)
{
    if (need[0] == 0 || update_error[0] != NEIGHBOR_LIST::UPDATE_OK) return;
    const int n = tile_count[0];
    for (int t = blockIdx.x * blockDim.x + threadIdx.x; t < n;
         t += gridDim.x * blockDim.x)
    {
        atomicAdd(row_cursor + tiles[t].cluster_i, 1);
    }
}
// 单块独占前缀扫描：把 row_cursor 从行计数改写为各行散射游标（行首槽位）。
// local cluster 数从 d_lj_cluster_numbers[0] 读取（host 侧重建末才回读）
static __global__ void Scan_LJ_Tile_Rows(const int* need,
                                         const int* cluster_numbers,
                                         int* row_cursor,
                                         const int* update_error)
{
    if (need[0] == 0 || update_error[0] != NEIGHBOR_LIST::UPDATE_OK) return;
    const int n = cluster_numbers[0];
    __shared__ int sh_scan[1024];
    __shared__ int sh_carry;
    if (threadIdx.x == 0) sh_carry = 0;
    __syncthreads();
    for (int base = 0; base < n; base += blockDim.x)
    {
        const int idx = base + threadIdx.x;
        const int v = idx < n ? row_cursor[idx] : 0;
        sh_scan[threadIdx.x] = v;
        __syncthreads();
        for (int offset = 1; offset < blockDim.x; offset <<= 1)
        {
            const int t =
                threadIdx.x >= offset ? sh_scan[threadIdx.x - offset] : 0;
            __syncthreads();
            sh_scan[threadIdx.x] += t;
            __syncthreads();
        }
        if (idx < n)
        {
            row_cursor[idx] = sh_carry + sh_scan[threadIdx.x] - v;
        }
        __syncthreads();
        if (threadIdx.x == blockDim.x - 1)
        {
            sh_carry += sh_scan[threadIdx.x];
        }
        __syncthreads();
    }
}

static __global__ void Scatter_LJ_Tile_Rows(const int* need,
                                            const LJ_TILE* tiles,
                                            const int* tile_count,
                                            int* row_cursor, int* tile_sorted,
                                            const int* update_error)
{
    if (need[0] == 0 || update_error[0] != NEIGHBOR_LIST::UPDATE_OK) return;
    const int n = tile_count[0];
    for (int t = blockIdx.x * blockDim.x + threadIdx.x; t < n;
         t += gridDim.x * blockDim.x)
    {
        tile_sorted[atomicAdd(row_cursor + tiles[t].cluster_i, 1)] = t;
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
    const int* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        update_error[0] != NEIGHBOR_LIST::UPDATE_OK)
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
    const int* update_error)
{
    if (need[0] == 0 || neighbor_grid_overflow[0] != 0 ||
        neighbor_grid_ghost_overflow[0] != 0 ||
        neighbor_list_overflow[0] != 0 ||
        update_error[0] != NEIGHBOR_LIST::UPDATE_OK)
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

#ifdef USE_GPU
// S2 调试开关（诊断路径，默认关闭，对数值零影响）：
// SPONGE_LJ_TILE_DEBUG=1 打印 tile 构建各阶段耗时与规模；
// SPONGE_LJ_TILE_VALIDATE=1 在每次成功重建后把 tile 表展开的 (global_i,
// global_j) 无序对集合与半表 d_nl 展开的对集合做逐对集合相等比较。
static bool LJ_Tile_Debug_Enabled()
{
    static int enabled = -1;
    if (enabled < 0)
    {
        const char* v = std::getenv("SPONGE_LJ_TILE_DEBUG");
        enabled = v != NULL && v[0] != '\0' && v[0] != '0' ? 1 : 0;
    }
    return enabled != 0;
}

static bool LJ_Tile_Validate_Enabled()
{
    static int enabled = -1;
    if (enabled < 0)
    {
        const char* v = std::getenv("SPONGE_LJ_TILE_VALIDATE");
        enabled = v != NULL && v[0] != '\0' && v[0] != '0' ? 1 : 0;
    }
    return enabled != 0;
}

// SPONGE_LJ_TILE_ABLATE：消融测时（仅诊断，结果不正确）：
// bit0 跳过排除检查，bit1 跳过距离测试，bit2 跳过整个任务循环，
// bit3 跳过 tile 评估（只跑任务管理+j 载入），bit4 跳过 j 数据载入
static int LJ_Tile_Ablate_Mode()
{
    static int mode = -1;
    if (mode < 0)
    {
        const char* v = std::getenv("SPONGE_LJ_TILE_ABLATE");
        mode = v != NULL ? atoi(v) : 0;
    }
    return mode;
}

static __global__ void Count_Half_List_Pair_Total(const ATOM_GROUP* nl, int n,
                                                  unsigned long long* total)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        atomicAdd(total, (unsigned long long)nl[i].atom_numbers);
    }
}

static __global__ void Expand_Half_List_Pair_Keys(
    const int* atom_local, const ATOM_GROUP* nl, int n,
    unsigned long long* keys, unsigned long long* counter)
{
    SIMPLE_DEVICE_FOR(i, n)
    {
        const unsigned int global_i = (unsigned int)atom_local[i];
        const int count = nl[i].atom_numbers;
        const int* serial = nl[i].atom_serial;
        for (int k = 0; k < count; ++k)
        {
            const unsigned int global_j = (unsigned int)atom_local[serial[k]];
            const unsigned int lo = global_i < global_j ? global_i : global_j;
            const unsigned int hi = global_i < global_j ? global_j : global_i;
            keys[atomicAdd(counter, 1ULL)] =
                ((unsigned long long)lo << 32) | hi;
        }
    }
}

static __global__ void Count_LJ_Tile_Pair_Total(const LJ_TILE* tiles, int n,
                                                unsigned long long* total)
{
    SIMPLE_DEVICE_FOR(t, n)
    {
        atomicAdd(total, (unsigned long long)__popcll(tiles[t].mask));
    }
}

static __global__ void Expand_LJ_Tile_Pair_Keys(
    const LJ_TILE* tiles, int n, const int* cluster_atoms,
    const int* atom_local, unsigned long long* keys,
    unsigned long long* counter, int* bad_padding)
{
    SIMPLE_DEVICE_FOR(t, n)
    {
        const LJ_TILE tile = tiles[t];
        int atoms_i[LJ_TILE_CLUSTER_SIZE], atoms_j[LJ_TILE_CLUSTER_SIZE];
        for (int k = 0; k < LJ_TILE_CLUSTER_SIZE; ++k)
        {
            atoms_i[k] = cluster_atoms[tile.cluster_i * LJ_TILE_CLUSTER_SIZE + k];
            atoms_j[k] = cluster_atoms[tile.cluster_j * LJ_TILE_CLUSTER_SIZE + k];
        }
        unsigned long long m = tile.mask;
        while (m != 0)
        {
            const int b = __ffsll(m) - 1;
            m &= m - 1;
            const int atom_i = atoms_i[b >> 3];
            const int atom_j = atoms_j[b & 7];
            if (atom_i < 0 || atom_j < 0)
            {
                // 掩码引用了 padding 位：构建 bug，计入错误
                atomicAdd(bad_padding, 1);
                continue;
            }
            const unsigned int global_i = (unsigned int)atom_local[atom_i];
            const unsigned int global_j = (unsigned int)atom_local[atom_j];
            const unsigned int lo = global_i < global_j ? global_i : global_j;
            const unsigned int hi = global_i < global_j ? global_j : global_i;
            keys[atomicAdd(counter, 1ULL)] =
                ((unsigned long long)lo << 32) | hi;
        }
    }
}

// tile 展开对集 vs 半表展开对集的逐对集合相等校验（诊断，仅
// SPONGE_LJ_TILE_VALIDATE=1 时触发；主机端排序后比对）。
static void Validate_LJ_Tile_Pair_Set(NEIGHBOR_LIST* list, int* atom_local,
                                      int local_atom_numbers)
{
    unsigned long long* d_total = NULL;
    unsigned long long* d_counter = NULL;
    int* d_bad_padding = NULL;
    unsigned long long* d_keys_half = NULL;
    unsigned long long* d_keys_tile = NULL;
    unsigned long long h_half_total = 0, h_tile_total = 0;
    int h_bad_padding = 0;
    const int threads = CONTROLLER::device_max_thread;
    if (!Device_Malloc_Safely((void**)&d_total, sizeof(unsigned long long)) ||
        !Device_Malloc_Safely((void**)&d_counter, sizeof(unsigned long long)) ||
        !Device_Malloc_Safely((void**)&d_bad_padding, sizeof(int)))
    {
        fprintf(stderr, "[lj-tile-validate] device alloc failed, skipped\n");
        goto cleanup;
    }

    // 第一段：半表展开
    deviceMemset(d_total, 0, sizeof(unsigned long long));
    Launch_Device_Kernel(Count_Half_List_Pair_Total,
                         (local_atom_numbers + threads - 1) / threads, threads,
                         0, NULL, list->d_nl, local_atom_numbers, d_total);
    deviceMemcpy(&h_half_total, d_total, sizeof(unsigned long long),
                 deviceMemcpyDeviceToHost);
    if (h_half_total > 0 &&
        !Device_Malloc_Safely((void**)&d_keys_half,
                              h_half_total * sizeof(unsigned long long)))
    {
        fprintf(stderr, "[lj-tile-validate] half key alloc failed\n");
        goto cleanup;
    }
    deviceMemset(d_counter, 0, sizeof(unsigned long long));
    Launch_Device_Kernel(Expand_Half_List_Pair_Keys,
                         (local_atom_numbers + threads - 1) / threads, threads,
                         0, NULL, atom_local, list->d_nl, local_atom_numbers,
                         d_keys_half, d_counter);

    // 第二段：tile 表展开
    deviceMemset(d_total, 0, sizeof(unsigned long long));
    Launch_Device_Kernel(Count_LJ_Tile_Pair_Total,
                         (list->h_lj_tile_numbers + threads - 1) / threads,
                         threads, 0, NULL, list->d_lj_tiles,
                         list->h_lj_tile_numbers, d_total);
    deviceMemcpy(&h_tile_total, d_total, sizeof(unsigned long long),
                 deviceMemcpyDeviceToHost);
    if (h_tile_total > 0 &&
        !Device_Malloc_Safely((void**)&d_keys_tile,
                              h_tile_total * sizeof(unsigned long long)))
    {
        fprintf(stderr, "[lj-tile-validate] tile key alloc failed\n");
        goto cleanup;
    }
    deviceMemset(d_counter, 0, sizeof(unsigned long long));
    deviceMemset(d_bad_padding, 0, sizeof(int));
    Launch_Device_Kernel(Expand_LJ_Tile_Pair_Keys,
                         (list->h_lj_tile_numbers + threads - 1) / threads,
                         threads, 0, NULL, list->d_lj_tiles,
                         list->h_lj_tile_numbers, list->d_lj_cluster_atoms,
                         atom_local, d_keys_tile, d_counter, d_bad_padding);
    deviceMemcpy(&h_bad_padding, d_bad_padding, sizeof(int),
                 deviceMemcpyDeviceToHost);

    {
        std::vector<unsigned long long> h_half(h_half_total),
            h_tile(h_tile_total);
        if (h_half_total > 0)
            deviceMemcpy(h_half.data(), d_keys_half,
                         h_half_total * sizeof(unsigned long long),
                         deviceMemcpyDeviceToHost);
        if (h_tile_total > 0)
            deviceMemcpy(h_tile.data(), d_keys_tile,
                         h_tile_total * sizeof(unsigned long long),
                         deviceMemcpyDeviceToHost);
        std::sort(h_half.begin(), h_half.end());
        std::sort(h_tile.begin(), h_tile.end());

        size_t only_half = 0, only_tile = 0, common = 0;
        size_t a = 0, b = 0;
        while (a < h_half.size() || b < h_tile.size())
        {
            if (a < h_half.size() &&
                (b == h_tile.size() || h_half[a] < h_tile[b]))
            {
                if (only_half < 5)
                    fprintf(stderr,
                            "[lj-tile-validate] only in half-list: "
                            "(%u, %u)\n",
                            (unsigned int)(h_half[a] >> 32),
                            (unsigned int)(h_half[a] & 0xffffffffu));
                ++only_half;
                ++a;
            }
            else if (b < h_tile.size() &&
                     (a == h_half.size() || h_tile[b] < h_half[a]))
            {
                if (only_tile < 5)
                    fprintf(stderr,
                            "[lj-tile-validate] only in tile-list: "
                            "(%u, %u)\n",
                            (unsigned int)(h_tile[b] >> 32),
                            (unsigned int)(h_tile[b] & 0xffffffffu));
                ++only_tile;
                ++b;
            }
            else
            {
                ++common;
                ++a;
                ++b;
            }
        }
        fprintf(stderr,
                "[lj-tile-validate] tiles=%d, half-list pairs=%llu, "
                "tile pairs=%llu, common=%zu, only_half=%zu, only_tile=%zu, "
                "bad_padding=%d\n",
                list->h_lj_tile_numbers, h_half_total, h_tile_total, common,
                only_half, only_tile, h_bad_padding);
        if (only_half == 0 && only_tile == 0 && h_bad_padding == 0)
            fprintf(stderr, "[lj-tile-validate] PAIR SETS EQUAL\n");
        else
            fprintf(stderr, "[lj-tile-validate] PAIR SETS DIFFER\n");
    }

cleanup:
    if (d_total != NULL) Free_Single_Device_Pointer((void**)&d_total);
    if (d_counter != NULL) Free_Single_Device_Pointer((void**)&d_counter);
    if (d_bad_padding != NULL)
        Free_Single_Device_Pointer((void**)&d_bad_padding);
    if (d_keys_half != NULL) Free_Single_Device_Pointer((void**)&d_keys_half);
    if (d_keys_tile != NULL) Free_Single_Device_Pointer((void**)&d_keys_tile);
}
#endif  // USE_GPU

void NEIGHBOR_LIST::UPDATOR::Update(
    int* atom_local, int local_atom_numbers, int ghost_numbers, int need_copy,
    VECTOR* crd, LTMatrix3 cell, LTMatrix3 rcell, NEIGHBOR_LIST::GRIDS* grids,
    int max_atom_in_grid_numbers, int max_ghost_in_grid_numbers,
    int max_neighbor_numbers, float grid_length, int* d_neighbor_grid_overflow,
    int* d_neighbor_grid_ghost_overflow, int* d_neighbor_list_overflow,
    int* d_update_error, ATOM_GROUP* d_nl, int* excluded_list_start,
    int* excluded_list, int* excluded_numbers, NEIGHBOR_LIST* owner)
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
#ifdef USE_GPU
    Launch_Device_Kernel(
        Build_Grid_Neighbor_Prefix,
        (grids->grid_numbers + CONTROLLER::device_max_thread - 1) /
            CONTROLLER::device_max_thread,
        CONTROLLER::device_max_thread, 0, NULL, d_need_update,
        grids->grid_numbers, grids->d_neighbor_grid_numbers,
        grids->d_neighbor_grids, grids->d_grid_atom_numbers,
        grids->d_grid_neighbor_prefix, d_neighbor_grid_overflow,
        d_neighbor_grid_ghost_overflow, d_update_error);
    Launch_Device_Kernel(
        Find_Neighbors_Gridly, grids->grid_numbers, 256,
        (size_t)(max_atom_in_grid_numbers *
                     (sizeof(VECTOR) + 2 * sizeof(int)) +
                 (MAX_GRID_NEIGHBORS + 1) * sizeof(int)),
        NULL, atom_local, local_atom_numbers, d_need_update,
        grids->grid_numbers, grids->d_neighbor_grid_numbers,
        grids->d_neighbor_grids, grids->d_grid_atom_crd, cell, rcell,
        max_atom_in_grid_numbers, d_nl, grid_length * grid_length * 4.0f,
        grids->d_grid_atom_numbers, grids->d_grid_atoms, max_neighbor_numbers,
        d_neighbor_list_overflow, grids->d_grid_neighbor_prefix,
        grids->d_grid_ghost_crd, max_ghost_in_grid_numbers,
        grids->d_grid_ghost_numbers, grids->d_grid_ghosts,
        d_neighbor_grid_overflow, d_neighbor_grid_ghost_overflow,
        d_update_error);

    // S2：cluster-pair tile 表构建（暂无消费者，与半表并存）。
    // owner 为空（不应发生）或 tile 表未分配（CPU 路径）时跳过。
    if (owner != NULL && owner->d_lj_tiles != NULL)
    {
        // 排除范围缓存：运行期静态，首次重建时预计算（排除表为 NULL 时
        // 得到全零范围，kernel 无需再判空）
        if (!owner->lj_excl_range_built)
        {
            Launch_Device_Kernel(
                Build_LJ_Excluded_Range,
                (owner->atom_numbers + CONTROLLER::device_max_thread - 1) /
                    CONTROLLER::device_max_thread,
                CONTROLLER::device_max_thread, 0, NULL, owner->atom_numbers,
                excluded_list_start, excluded_list, excluded_numbers,
                owner->d_lj_excl_range);
            owner->lj_excl_range_built = 1;
        }
        const bool tile_debug = LJ_Tile_Debug_Enabled();
        deviceEvent_t tile_events[4] = {nullptr, nullptr, nullptr, nullptr};
        if (tile_debug)
        {
            for (int e = 0; e < 4; ++e) deviceEventCreate(&tile_events[e]);
            deviceEventRecord(tile_events[0], 0);
        }
        Launch_Device_Kernel(
            Build_LJ_Cluster_Prefix, 1, CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread * sizeof(int), NULL, d_need_update,
            grids->grid_numbers, grids->d_grid_atom_numbers,
            grids->d_grid_ghost_numbers, owner->d_grid_cluster_base,
            owner->d_grid_ghost_cluster_base, owner->d_lj_cluster_numbers,
            d_neighbor_grid_overflow, d_neighbor_grid_ghost_overflow,
            d_update_error);
        if (tile_debug) deviceEventRecord(tile_events[1], 0);
        Launch_Device_Kernel(
            Fill_LJ_Cluster_Atoms,
            (grids->grid_numbers + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, d_need_update,
            grids->grid_numbers, grids->d_grid_atom_numbers,
            grids->d_grid_atoms, owner->d_grid_cluster_base,
            max_atom_in_grid_numbers, grids->d_grid_ghost_numbers,
            grids->d_grid_ghosts, owner->d_grid_ghost_cluster_base,
            max_ghost_in_grid_numbers, owner->d_lj_cluster_atoms,
            owner->d_lj_cluster_flags, d_neighbor_grid_overflow,
            d_neighbor_grid_ghost_overflow, d_update_error);
        if (tile_debug) deviceEventRecord(tile_events[2], 0);
// 消融测时模板分发（仅诊断；生产路径 ABLATE=0，分支编译期消除）
#define SPONGE_LAUNCH_LJ_TILE_LIST(ABLATE_MODE)                             \
    Launch_Device_Kernel(                                                   \
        Build_LJ_Tile_List<ABLATE_MODE>, grids->grid_numbers, 256,          \
        (size_t)(max_atom_in_grid_numbers *                                 \
                     (sizeof(float4) + sizeof(LJ_EXCL_RANGE) + sizeof(int)) +        \
                 256 * sizeof(int) + LJ_TILE_BLOCK_BUFFER * sizeof(LJ_TILE)), \
        NULL, d_need_update, grids->grid_numbers,                           \
        grids->d_neighbor_grid_numbers, grids->d_neighbor_grids,            \
        grids->d_grid_atom_crd, cell, rcell, max_atom_in_grid_numbers,      \
        grids->d_grid_atom_numbers, grids->d_grid_atoms,                    \
        owner->d_grid_cluster_base, grids->d_grid_ghost_crd,                \
        max_ghost_in_grid_numbers, grids->d_grid_ghost_numbers,             \
        grids->d_grid_ghosts, owner->d_grid_ghost_cluster_base, atom_local, \
        grid_length * grid_length * 4.0f, excluded_list,                  \
        owner->d_lj_excl_range, owner->d_lj_tiles,                        \
        owner->lj_tile_capacity, owner->d_lj_tile_count,                    \
        owner->d_neighbor_tile_overflow, d_neighbor_grid_overflow,          \
        d_neighbor_grid_ghost_overflow, d_update_error)
        switch (LJ_Tile_Ablate_Mode())
        {
            case 1: SPONGE_LAUNCH_LJ_TILE_LIST(1); break;
            case 2: SPONGE_LAUNCH_LJ_TILE_LIST(2); break;
            case 3: SPONGE_LAUNCH_LJ_TILE_LIST(3); break;
            case 4: SPONGE_LAUNCH_LJ_TILE_LIST(4); break;
            case 8: SPONGE_LAUNCH_LJ_TILE_LIST(8); break;
            case 24: SPONGE_LAUNCH_LJ_TILE_LIST(24); break;
            default: SPONGE_LAUNCH_LJ_TILE_LIST(0); break;
        }
#undef SPONGE_LAUNCH_LJ_TILE_LIST
        if (tile_debug)
        {
            deviceEventRecord(tile_events[3], 0);
            deviceEventSynchronize(tile_events[3]);
            float ms[3] = {0, 0, 0};
            for (int e = 0; e < 3; ++e)
                deviceEventElapsedTime(&ms[e], tile_events[e],
                                       tile_events[e + 1]);
            int h_clusters[2] = {0, 0};
            deviceMemcpy(h_clusters, owner->d_lj_cluster_numbers,
                         sizeof(h_clusters), deviceMemcpyDeviceToHost);
            fprintf(stderr,
                    "[lj-tile] rebuild: cluster_prefix=%.3f ms, "
                    "cluster_fill=%.3f ms, tile_build=%.3f ms, "
                    "clusters=%d+%d\n",
                    ms[0], ms[1], ms[2], h_clusters[0], h_clusters[1]);
            for (int e = 0; e < 4; ++e) deviceEventDestroy(tile_events[e]);
        }
        // S3：tile 按 cluster_i 分组排序（tile kernel 按行消费，i 侧全局
        // 原子加降一个量级）；四个内核自带 need/错误短路，非重建步零开销
        Launch_Device_Kernel(
            Zero_LJ_Tile_Rows,
            (owner->lj_cluster_capacity + CONTROLLER::device_max_thread - 1) /
                CONTROLLER::device_max_thread,
            CONTROLLER::device_max_thread, 0, NULL, d_need_update,
            owner->d_lj_tile_row_cursor, owner->lj_cluster_capacity,
            d_update_error);
        Launch_Device_Kernel(Count_LJ_Tile_Rows, 2048, 256, 0, NULL,
                             d_need_update, owner->d_lj_tiles,
                             owner->d_lj_tile_count, owner->d_lj_tile_row_cursor,
                             d_update_error);
        Launch_Device_Kernel(Scan_LJ_Tile_Rows, 1, 1024, 0, NULL, d_need_update,
                             owner->d_lj_cluster_numbers,
                             owner->d_lj_tile_row_cursor, d_update_error);
        Launch_Device_Kernel(Scatter_LJ_Tile_Rows, 2048, 256, 0, NULL,
                             d_need_update, owner->d_lj_tiles,
                             owner->d_lj_tile_count, owner->d_lj_tile_row_cursor,
                             owner->d_lj_tile_sorted, d_update_error);
    }
#else
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
#endif

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
    if (!Device_Malloc_Safely((void**)&d_update_status,
                              UPDATE_STATUS_INTS * sizeof(int)) ||
        !Malloc_Safely((void**)&h_nl, group_bytes) ||
        !Device_Malloc_Safely((void**)&d_temp, neighbor_bytes))
    {
        Clear();
        controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                       "NEIGHBOR_LIST::Initial");
        return;
    }
    // 状态块内各字段的别名
    d_neighbor_grid_overflow = d_update_status + UPDATE_STATUS_GRID_OVERFLOW;
    d_neighbor_grid_ghost_overflow =
        d_update_status + UPDATE_STATUS_GRID_GHOST_OVERFLOW;
    d_neighbor_list_overflow = d_update_status + UPDATE_STATUS_LIST_OVERFLOW;
    d_update_error = d_update_status + UPDATE_STATUS_ERROR;
    deviceMemset(d_update_status, 0, UPDATE_STATUS_INTS * sizeof(int));
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

#ifdef USE_GPU
    // S2：LJ cluster-pair tile 表存储（暂无消费者）。cluster 槽数上界由
    // 格数 × 每格容量决定，构建时不会溢出；tile 表容量按原子数预估，
    // 不够用 Record_Required_Capacity + 溢出恢复扩容重建兜底。
    {
        const int clusters_per_grid =
            (max_atom_in_grid_numbers + LJ_TILE_CLUSTER_SIZE - 1) /
                LJ_TILE_CLUSTER_SIZE +
            (max_ghost_in_grid_numbers + LJ_TILE_CLUSTER_SIZE - 1) /
                LJ_TILE_CLUSTER_SIZE;
        size_t cluster_slots = 0, cluster_atom_slots = 0, tile_slots = 0;
        size_t base_bytes = 0, cluster_atom_bytes = 0, flag_bytes = 0;
        size_t tile_bytes = 0;
        // 初始容量按实测 ~17 tile/原子（skin=4.0, 55.6 万原子 ~940 万
        // tile）估到 20 tile/原子，溢出恢复兜底
        int64_t default_tiles = static_cast<int64_t>(atom_numbers) * 20;
        if (default_tiles > std::numeric_limits<int>::max())
            default_tiles = std::numeric_limits<int>::max();
        lj_tile_capacity = lj_tile_capacity_hint > default_tiles
                               ? lj_tile_capacity_hint
                               : static_cast<int>(default_tiles);
        if (!Checked_Product((size_t)grids.grid_numbers,
                             (size_t)clusters_per_grid, &cluster_slots) ||
            cluster_slots > (size_t)std::numeric_limits<int>::max() ||
            !Checked_Product(cluster_slots, (size_t)LJ_TILE_CLUSTER_SIZE,
                             &cluster_atom_slots) ||
            !Checked_Bytes((size_t)grids.grid_numbers, sizeof(int),
                           &base_bytes) ||
            !Checked_Bytes(cluster_atom_slots, sizeof(int),
                           &cluster_atom_bytes) ||
            !Checked_Bytes(cluster_slots, sizeof(int), &flag_bytes) ||
            !Checked_Bytes((size_t)lj_tile_capacity, sizeof(LJ_TILE),
                           &tile_bytes))
        {
            Clear();
            controller->Throw_SPONGE_Error(
                spongeErrorOverflow, "NEIGHBOR_LIST::Initial",
                "Reason:\n\tLJ tile-list allocation size overflow\n");
            return;
        }
        lj_cluster_capacity = static_cast<int>(cluster_slots);
        if (!Device_Malloc_Safely((void**)&d_grid_cluster_base, base_bytes) ||
            !Device_Malloc_Safely((void**)&d_grid_ghost_cluster_base,
                                  base_bytes) ||
            !Device_Malloc_Safely((void**)&d_lj_cluster_numbers,
                                  2 * sizeof(int)) ||
            !Device_Malloc_Safely((void**)&d_lj_cluster_atoms,
                                  cluster_atom_bytes) ||
            !Device_Malloc_Safely((void**)&d_lj_cluster_flags, flag_bytes) ||
            !Device_Malloc_Safely((void**)&d_lj_excl_range,
                                  (size_t)atom_numbers * sizeof(LJ_EXCL_RANGE)) ||
            !Device_Malloc_Safely((void**)&d_lj_tiles, tile_bytes) ||
            !Device_Malloc_Safely((void**)&d_lj_tile_row_cursor,
                                  flag_bytes) ||
            !Device_Malloc_Safely((void**)&d_lj_tile_sorted,
                                  (size_t)lj_tile_capacity * sizeof(int)))
        {
            Clear();
            controller->Throw_SPONGE_Error(spongeErrorMallocFailed,
                                           "NEIGHBOR_LIST::Initial");
            return;
        }
        d_neighbor_tile_overflow = d_update_status + UPDATE_STATUS_TILE_OVERFLOW;
        d_lj_tile_count = d_update_status + UPDATE_STATUS_TILE_COUNT;
        controller->printf(
            "    LJ tile list: cluster capacity %d, tile capacity %d "
            "(%zu MB)\n",
            lj_cluster_capacity, lj_tile_capacity, tile_bytes >> 20);
    }
#endif

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
            invalid_geometry ? UPDATE_INVALID_GEOMETRY : UPDATE_OK;
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
    h_neighbor_grid_overflow = 0;
    h_neighbor_grid_ghost_overflow = 0;
    h_neighbor_list_overflow = 0;
    h_neighbor_tile_overflow = 0;
    last_update_error = UPDATE_OK;
    last_error_atom = -1;
    if (full_neighbor_list.is_initialized)
    {
        deviceMemset(full_neighbor_list.d_overflow, 0, sizeof(int));
    }

    // 先在主机侧判定本步是否需要重建：强制/定间隔由参数与步数直接可知；
    // 动态模式跑 Check kernel 并取回标志（非重建步唯一的同步点）。
    // 只有重建步才清零状态块、启动构建 kernel 并取回构建状态。
    bool need_rebuild = false;
    if (update == NEIGHBOR_LIST_UPDATE_PARAMETER::FORCED_UPDATE)
    {
        need_rebuild = true;
    }
    else if (updator.refresh_interval <= 0)
    {
        deviceMemset(updator.d_need_update, 0, sizeof(int));
        updator.Check(local_atom_numbers, skin, crd, cell, rcell);
        int h_need_update = 0;
        deviceMemcpy(&h_need_update, updator.d_need_update, sizeof(int),
                     deviceMemcpyDeviceToHost);
        need_rebuild = (h_need_update != 0);
    }
    else
    {
        need_rebuild = Next_Step_Is_Interval_Boundary(
            step, updator.refresh_interval);
    }

    if (need_rebuild && this->is_needed_half &&
        local_atom_numbers + ghost_numbers > 0)
    {
        deviceMemset(d_update_status, 0, UPDATE_STATUS_INTS * sizeof(int));
        if (update == NEIGHBOR_LIST_UPDATE_PARAMETER::FORCED_UPDATE ||
            updator.refresh_interval > 0)
        {
            deviceMemset(updator.d_need_update, -1, sizeof(int));
        }
        // 动态模式下 Check kernel 已把 d_need_update 置 1
        updator.Update(atom_local, local_atom_numbers, ghost_numbers,
                       updator.refresh_interval <= 0, crd, cell, rcell, &grids,
                       max_atom_in_grid_numbers, max_ghost_in_grid_numbers,
                       max_neighbor_numbers, 0.5f * (cutoff + skin),
                       d_neighbor_grid_overflow, d_neighbor_grid_ghost_overflow,
                       d_neighbor_list_overflow, d_update_error, this->d_nl,
                       excluded_list_start, excluded_list, excluded_numbers,
                       this);
        // 阻塞式取回全部构建状态
        int status[UPDATE_STATUS_INTS] = {0};
        deviceMemcpy(status, d_update_status, sizeof(status),
                     deviceMemcpyDeviceToHost);
        h_neighbor_grid_overflow = status[UPDATE_STATUS_GRID_OVERFLOW];
        h_neighbor_grid_ghost_overflow =
            status[UPDATE_STATUS_GRID_GHOST_OVERFLOW];
        h_neighbor_list_overflow = status[UPDATE_STATUS_LIST_OVERFLOW];
        h_neighbor_tile_overflow = status[UPDATE_STATUS_TILE_OVERFLOW];
        h_lj_tile_numbers = status[UPDATE_STATUS_TILE_COUNT];
#ifdef USE_GPU
        if (d_lj_cluster_numbers != NULL)
        {
            int h_clusters[2] = {0, 0};
            deviceMemcpy(h_clusters, d_lj_cluster_numbers, sizeof(h_clusters),
                         deviceMemcpyDeviceToHost);
            h_lj_cluster_atom_slots =
                (h_clusters[0] + h_clusters[1]) * LJ_TILE_CLUSTER_SIZE;
        }
#endif
        last_update_error = status[UPDATE_STATUS_ERROR];
        last_error_atom = status[UPDATE_STATUS_ERROR + 1];
        if (last_update_error != UPDATE_OK || h_neighbor_grid_overflow != 0 ||
            h_neighbor_grid_ghost_overflow != 0 ||
            h_neighbor_list_overflow != 0 || h_neighbor_tile_overflow != 0)
        {
            Invalidate_Half_Neighbor_List(this, local_atom_numbers);
            full_neighbor_list.Invalidate_Active();
            updator.time_recorder->Stop();
            return false;
        }
#ifdef USE_GPU
        // S2 诊断路径（默认关闭）：tile 规模打印与对集校验
        if (d_lj_tiles != NULL && LJ_Tile_Debug_Enabled())
        {
            fprintf(stderr, "[lj-tile] rebuild done: tiles=%d (capacity %d)\n",
                    h_lj_tile_numbers, lj_tile_capacity);
        }
        if (d_lj_tiles != NULL && LJ_Tile_Validate_Enabled())
        {
            Validate_LJ_Tile_Pair_Set(this, atom_local, local_atom_numbers);
        }
#endif
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
        const bool tile_capacity_failure =
            d_lj_tiles != NULL && h_neighbor_tile_overflow > lj_tile_capacity;
        const bool capacity_failure =
            grid_atom_capacity_failure || grid_ghost_capacity_failure ||
            half_capacity_failure || full_capacity_failure ||
            tile_capacity_failure;

        if (!capacity_failure)
        {
            if (last_update_error != UPDATE_OK)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "NEIGHBOR_LIST::Update_With_Overflow_Recovery",
                    "Reason:\n\tneighbor-list build rejected non-finite or "
                    "invalid periodic geometry (atom=%d)\n",
                    last_error_atom);
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
        // tile 表容量不走 mdin 命令，直接抬高跨 Clear/Initial 保留的提示
        //（加 25% 余量避免紧贴最小需求反复触发重建）
        if (tile_capacity_failure)
        {
            lj_tile_capacity_hint =
                h_neighbor_tile_overflow + h_neighbor_tile_overflow / 4;
        }
        controller->printf(
            "Neighbor-list capacity was insufficient for the current state; "
            "rebuilding exactly with grid atoms=%d, grid ghosts=%d, "
            "neighbors=%d (grid_atoms=%d, grid_ghosts=%d, half=%d, "
            "full=%d, tiles=%d).\n",
            grown_grid_atoms, grown_grid_ghosts, grown_neighbors,
            h_neighbor_grid_overflow, h_neighbor_grid_ghost_overflow,
            h_neighbor_list_overflow, required_full_neighbors,
            h_neighbor_tile_overflow);

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
    // 状态块单独释放；四个字段指针只是块内别名，置空即可
    if (d_update_status != NULL)
        Free_Single_Device_Pointer((void**)&d_update_status);
    d_neighbor_grid_overflow = NULL;
    d_neighbor_list_overflow = NULL;
    d_neighbor_grid_ghost_overflow = NULL;
    d_update_error = NULL;
    // S2 tile 表存储释放；lj_tile_capacity_hint 刻意保留，供溢出恢复后的
    // Initial 按提示扩容
    if (d_grid_cluster_base != NULL)
        Free_Single_Device_Pointer((void**)&d_grid_cluster_base);
    if (d_grid_ghost_cluster_base != NULL)
        Free_Single_Device_Pointer((void**)&d_grid_ghost_cluster_base);
    if (d_lj_cluster_numbers != NULL)
        Free_Single_Device_Pointer((void**)&d_lj_cluster_numbers);
    if (d_lj_cluster_atoms != NULL)
        Free_Single_Device_Pointer((void**)&d_lj_cluster_atoms);
    if (d_lj_cluster_flags != NULL)
        Free_Single_Device_Pointer((void**)&d_lj_cluster_flags);
    if (d_lj_tiles != NULL) Free_Single_Device_Pointer((void**)&d_lj_tiles);
    if (d_lj_tile_row_cursor != NULL)
        Free_Single_Device_Pointer((void**)&d_lj_tile_row_cursor);
    if (d_lj_tile_sorted != NULL)
        Free_Single_Device_Pointer((void**)&d_lj_tile_sorted);
    if (d_lj_excl_range != NULL)
        Free_Single_Device_Pointer((void**)&d_lj_excl_range);
    lj_excl_range_built = 0;
    lj_cluster_capacity = 0;
    lj_tile_capacity = 0;
    h_lj_tile_numbers = 0;
    h_lj_cluster_atom_slots = 0;
    h_neighbor_tile_overflow = 0;
    d_neighbor_tile_overflow = NULL;
    d_lj_tile_count = NULL;
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
