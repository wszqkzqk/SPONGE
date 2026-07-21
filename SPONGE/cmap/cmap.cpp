#include "cmap.h"

#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "../xponge/load/native/cmap.hpp"
#include "../xponge/xponge.h"

// clang-format off
// 由于求导带来的系数矩阵的逆矩阵A_inv
static const float A_inv[16][16] =
{ { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, },
{ -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, },
{ 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, },
{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, },
{ 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, },
{ 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, },
{ -3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{ 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, },
{ 9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1, },
{ -6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1, },
{ 2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
{ 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, },
{ -6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1, },
{ 4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1 } };
// clang-format on

static bool CMap_Float_Is_Finite(float value)
{
    return Float_Memory_Is_Finite(&value);
}

static bool CMap_Float_Is_Zero_Or_Normal(float value)
{
    return Float_Memory_Is_Zero_Or_Normal(&value);
}

static std::string Validate_CMap(const Xponge::CMap& cmap,
                                 std::size_t atom_numbers)
{
    const std::size_t cmap_numbers = cmap.atom_a.size();
    if (cmap_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return "CMAP interaction count exceeds the kernel index range";
    }
    if (atom_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()))
    {
        return "system atom count exceeds the CMAP kernel index range";
    }
    if (cmap.atom_b.size() != cmap_numbers ||
        cmap.atom_c.size() != cmap_numbers ||
        cmap.atom_d.size() != cmap_numbers ||
        cmap.atom_e.size() != cmap_numbers ||
        cmap.cmap_type.size() != cmap_numbers)
    {
        return "CMAP interaction arrays have inconsistent lengths";
    }
    if (cmap.unique_type_numbers < 0)
    {
        return "CMAP unique type count is negative";
    }
    const std::size_t type_numbers =
        static_cast<std::size_t>(cmap.unique_type_numbers);
    if (cmap.resolution.size() != type_numbers ||
        cmap.type_offset.size() != type_numbers)
    {
        return "CMAP type arrays do not match the unique type count";
    }

    std::size_t gridpoint_numbers = 0;
    for (std::size_t type = 0; type < type_numbers; type++)
    {
        const int resolution = cmap.resolution[type];
        if (resolution <= 0)
        {
            return "CMAP type " + std::to_string(type) +
                   " has a non-positive resolution";
        }
        const std::size_t n = static_cast<std::size_t>(resolution);
        if (n > std::numeric_limits<std::size_t>::max() / n)
        {
            return "CMAP type " + std::to_string(type) +
                   " resolution overflows the grid size";
        }
        const std::size_t type_gridpoint_numbers = n * n;
        if (gridpoint_numbers >
            std::numeric_limits<std::size_t>::max() - type_gridpoint_numbers)
        {
            return "CMAP grid size overflows";
        }
        if (gridpoint_numbers >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) / 16)
        {
            return "CMAP interpolation coefficient offset overflows int";
        }
        const int expected_offset = static_cast<int>(16 * gridpoint_numbers);
        if (cmap.type_offset[type] != expected_offset)
        {
            return "CMAP type " + std::to_string(type) +
                   " has an inconsistent interpolation offset";
        }
        gridpoint_numbers += type_gridpoint_numbers;
    }
    if (gridpoint_numbers >
        static_cast<std::size_t>(std::numeric_limits<int>::max()) / 16)
    {
        return "CMAP interpolation table exceeds the kernel index range";
    }
    if (cmap.unique_gridpoint_numbers != static_cast<int>(gridpoint_numbers))
    {
        return "CMAP grid point count is inconsistent with its resolutions";
    }
    if (cmap.grid_value.size() != gridpoint_numbers)
    {
        return "CMAP grid value count is inconsistent with its resolutions";
    }
    for (std::size_t point = 0; point < cmap.grid_value.size(); point++)
    {
        if (!CMap_Float_Is_Finite(cmap.grid_value[point]) ||
            !CMap_Float_Is_Zero_Or_Normal(cmap.grid_value[point]))
        {
            return "CMAP grid point " + std::to_string(point) +
                   " is not a finite zero or normal float";
        }
    }

    for (std::size_t term = 0; term < cmap_numbers; term++)
    {
        const int type = cmap.cmap_type[term];
        if (type < 0 || static_cast<std::size_t>(type) >= type_numbers)
        {
            return "CMAP interaction " + std::to_string(term) +
                   " has an out-of-range type";
        }
        const int atoms[5] = {cmap.atom_a[term], cmap.atom_b[term],
                              cmap.atom_c[term], cmap.atom_d[term],
                              cmap.atom_e[term]};
        for (int atom : atoms)
        {
            if (atom < 0 || static_cast<std::size_t>(atom) >= atom_numbers)
            {
                return "CMAP interaction " + std::to_string(term) +
                       " has atom index " + std::to_string(atom) +
                       " outside the system";
            }
        }
        for (int torsion = 0; torsion < 2; torsion++)
        {
            for (int first = torsion; first < torsion + 4; first++)
            {
                for (int second = first + 1; second < torsion + 4; second++)
                {
                    if (atoms[first] == atoms[second])
                    {
                        return "CMAP interaction " + std::to_string(term) +
                               " repeats atom " + std::to_string(atoms[first]) +
                               " within its " +
                               (torsion == 0 ? "first" : "second") + " torsion";
                    }
                }
            }
        }
    }
    return "";
}

// Match the CMAP preprocessing used by GROMACS.  It first tiles the periodic
// grid to twice its size, then uses an ordinary natural cubic spline on the
// tiled rows and columns.  Taking derivatives only from the central copy
// makes the artificial zero-curvature boundaries remote from the map seam.
static void CMap_Natural_Spline_Second_Derivatives(
    const std::vector<double>& values, std::vector<double>* second_derivatives)
{
    const std::size_t n = values.size();
    second_derivatives->assign(n, 0.0);
    if (n <= 2)
    {
        return;
    }

    std::vector<double> work(n, 0.0);
    for (std::size_t i = 1; i + 1 < n; i++)
    {
        const double pivot = 0.5 * (*second_derivatives)[i - 1] + 2.0;
        (*second_derivatives)[i] = -0.5 / pivot;
        const double second_difference =
            values[i + 1] - 2.0 * values[i] + values[i - 1];
        work[i] = (3.0 * second_difference - 0.5 * work[i - 1]) / pivot;
    }
    for (std::size_t i = n - 1; i-- > 0;)
    {
        (*second_derivatives)[i] =
            (*second_derivatives)[i] * (*second_derivatives)[i + 1] + work[i];
    }
}

static void CMap_Natural_Spline_Interpolate(
    const std::vector<double>& values,
    const std::vector<double>& second_derivatives, double coordinate,
    double* value, double* derivative)
{
    const std::size_t left = static_cast<std::size_t>(coordinate);
    const double lower_weight = static_cast<double>(left + 1) - coordinate;
    const double upper_weight = coordinate - static_cast<double>(left);
    *value = lower_weight * values[left] + upper_weight * values[left + 1] +
             ((lower_weight * lower_weight * lower_weight - lower_weight) *
                  second_derivatives[left] +
              (upper_weight * upper_weight * upper_weight - upper_weight) *
                  second_derivatives[left + 1]) /
                 6.0;
    *derivative = values[left + 1] - values[left] -
                  (3.0 * lower_weight * lower_weight - 1.0) / 6.0 *
                      second_derivatives[left] +
                  (3.0 * upper_weight * upper_weight - 1.0) / 6.0 *
                      second_derivatives[left + 1];
}

void CMAP::Initial(CONTROLLER* controller, const char* module_name)
{
    this->controller = controller;
    if (module_name == NULL)
    {
        strcpy(this->module_name, "cmap");
    }
    else
    {
        strcpy(this->module_name, module_name);
    }

    const auto& cmap = Xponge::system.classical_force_field.cmap;
    Xponge::CMap local_cmap;
    const Xponge::CMap* cmap_to_use = NULL;
    const char* init_source = NULL;
    if (module_name == NULL)
    {
        cmap_to_use = &cmap;
        init_source = "Xponge::system";
    }
    else if (controller->Command_Exist(this->module_name, "in_file"))
    {
        Xponge::Native_Load_CMap(&local_cmap, controller, this->module_name);
        cmap_to_use = &local_cmap;
    }

    tot_cmap_num = 0;
    uniq_cmap_num = 0;
    uniq_gridpoint_num = 0;
    if (cmap_to_use != NULL)
    {
        const std::string validation_error =
            Validate_CMap(*cmap_to_use, Xponge::system.atoms.mass.size());
        if (!validation_error.empty())
        {
            const std::string reason =
                "Reason:\n\tinvalid CMAP data: " + validation_error + "\n";
            controller->Throw_SPONGE_Error(spongeErrorBadFileFormat,
                                           "CMAP::Initial", reason.c_str());
        }
        tot_cmap_num = static_cast<int>(cmap_to_use->atom_a.size());
        uniq_cmap_num = cmap_to_use->unique_type_numbers;
        uniq_gridpoint_num = cmap_to_use->unique_gridpoint_numbers;
    }
    if (tot_cmap_num > 0)
    {
        if (module_name == NULL)
        {
            controller->printf("START INITIALIZING CMAP (%s):\n", init_source);
        }
        else
        {
            controller->printf("START INITIALIZING CMAP (%s_in_file):\n",
                               this->module_name);
        }
        controller->printf(
            "    total CMAP number is %d\n    unique CMAP number is %d\n",
            tot_cmap_num, uniq_cmap_num);
        this->Memory_Allocate();
        Malloc_Safely((void**)&grid_value, sizeof(float) * uniq_gridpoint_num);
        Malloc_Safely((void**)&h_inter_coeff,
                      sizeof(float) * 16 * uniq_gridpoint_num);
        for (int i = 0; i < uniq_cmap_num; i++)
        {
            h_cmap_resolution[i] = cmap_to_use->resolution[i];
            type_offset[i] = cmap_to_use->type_offset[i];
        }
        memcpy(grid_value, cmap_to_use->grid_value.data(),
               sizeof(float) * uniq_gridpoint_num);
        for (int i = 0; i < tot_cmap_num; i++)
        {
            h_atom_a[i] = cmap_to_use->atom_a[i];
            h_atom_b[i] = cmap_to_use->atom_b[i];
            h_atom_c[i] = cmap_to_use->atom_c[i];
            h_atom_d[i] = cmap_to_use->atom_d[i];
            h_atom_e[i] = cmap_to_use->atom_e[i];
            h_cmap_type[i] = cmap_to_use->cmap_type[i];
        }
        is_initialized = 1;
    }

    if (is_initialized && !is_controller_printf_initialized)
    {
        controller->Step_Print_Initial(this->module_name, "%.6f");
        is_controller_printf_initialized = 1;
        controller->printf("    structure last modify date is %d\n",
                           last_modify_date);
    }
    if (is_initialized)
    {
        // 完成插值系数计算，完成初始化
        this->Interpolation(controller);
        Parameter_Host_to_Device();
        controller->printf("END INITIALIZING CMAP\n\n");
    }
    else
    {
        controller->printf("CMAP IS NOT INITIALIZED\n\n");
    }
}

void CMAP::Parameter_Host_to_Device()
{
    Device_Malloc_And_Copy_Safely((void**)&d_atom_a, h_atom_a,
                                  sizeof(int) * tot_cmap_num);
    Device_Malloc_And_Copy_Safely((void**)&d_atom_b, h_atom_b,
                                  sizeof(int) * tot_cmap_num);
    Device_Malloc_And_Copy_Safely((void**)&d_atom_c, h_atom_c,
                                  sizeof(int) * tot_cmap_num);
    Device_Malloc_And_Copy_Safely((void**)&d_atom_d, h_atom_d,
                                  sizeof(int) * tot_cmap_num);
    Device_Malloc_And_Copy_Safely((void**)&d_atom_e, h_atom_e,
                                  sizeof(int) * tot_cmap_num);
    Device_Malloc_And_Copy_Safely((void**)&d_cmap_type, h_cmap_type,
                                  sizeof(int) * tot_cmap_num);
    Device_Malloc_And_Copy_Safely((void**)&d_inter_coeff, h_inter_coeff,
                                  sizeof(float) * 16 * uniq_gridpoint_num);
    Device_Malloc_And_Copy_Safely((void**)&d_cmap_resolution, h_cmap_resolution,
                                  sizeof(int) * uniq_cmap_num);
    Device_Malloc_And_Copy_Safely((void**)&d_type_offset, type_offset,
                                  sizeof(int) * uniq_cmap_num);

    Device_Malloc_Safely((void**)&d_atom_a_local, sizeof(int) * tot_cmap_num);
    Device_Malloc_Safely((void**)&d_atom_b_local, sizeof(int) * tot_cmap_num);
    Device_Malloc_Safely((void**)&d_atom_c_local, sizeof(int) * tot_cmap_num);
    Device_Malloc_Safely((void**)&d_atom_d_local, sizeof(int) * tot_cmap_num);
    Device_Malloc_Safely((void**)&d_atom_e_local, sizeof(int) * tot_cmap_num);
    Device_Malloc_Safely((void**)&d_cmap_type_local,
                         sizeof(int) * tot_cmap_num);
    Device_Malloc_Safely((void**)&d_cmap_global_index_local,
                         sizeof(int) * tot_cmap_num);
    Device_Malloc_Safely((void**)&d_num_cmap_local, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_term, sizeof(int));
    Device_Malloc_Safely((void**)&d_invalid_local_atom, sizeof(int));
#ifndef GPU_ARCH_NAME
    Device_Malloc_Safely((void**)&d_invalid_geometry_term, sizeof(int));
#endif
    deviceMemset(d_num_cmap_local, 0, sizeof(int));
}

void CMAP::Memory_Allocate()
{
    Malloc_Safely((void**)&h_cmap_resolution, sizeof(int) * uniq_cmap_num);
    Malloc_Safely((void**)&h_cmap_type, sizeof(int) * tot_cmap_num);
    Malloc_Safely((void**)&h_atom_a, sizeof(int) * tot_cmap_num);
    Malloc_Safely((void**)&h_atom_b, sizeof(int) * tot_cmap_num);
    Malloc_Safely((void**)&h_atom_c, sizeof(int) * tot_cmap_num);
    Malloc_Safely((void**)&h_atom_d, sizeof(int) * tot_cmap_num);
    Malloc_Safely((void**)&h_atom_e, sizeof(int) * tot_cmap_num);
    Malloc_Safely((void**)&type_offset, sizeof(int) * uniq_cmap_num);
    Device_Malloc_Safely((void**)&d_sigma_of_cmap_ene, sizeof(float));
    Device_Malloc_Safely((void**)&d_cmap_ene, sizeof(float) * tot_cmap_num);
}

void CMAP::Interpolation(CONTROLLER* controller)
{
    controller->printf("    Start Interpolating the CMAP Grid Value\n");
    std::size_t grid_offset = 0;
    for (int type = 0; type < uniq_cmap_num; type++)
    {
        const std::size_t resolution =
            static_cast<std::size_t>(h_cmap_resolution[type]);
        const std::size_t grid_size = resolution * resolution;
        std::vector<double> values(grid_size, 0.0);
        for (std::size_t i = 0; i < grid_size; i++)
        {
            values[i] = grid_value[grid_offset + i];
        }

        const std::size_t extended_resolution = 2 * resolution;
        const std::size_t extended_grid_size =
            extended_resolution * extended_resolution;
        std::vector<double> extended_grid(extended_grid_size, 0.0);
        const std::size_t extension_shift = resolution - resolution / 2;
        for (std::size_t phi = 0; phi < extended_resolution; phi++)
        {
            const std::size_t source_phi = (phi + extension_shift) % resolution;
            for (std::size_t psi = 0; psi < extended_resolution; psi++)
            {
                const std::size_t source_psi =
                    (psi + extension_shift) % resolution;
                extended_grid[phi * extended_resolution + psi] =
                    values[source_phi * resolution + source_psi];
            }
        }

        // First interpolate every extended phi row at the psi coordinates of
        // the original map.  For odd resolutions these coordinates fall
        // halfway between tiled samples; retaining that behavior is required
        // to reproduce GROMACS' setup_cmap rather than an idealized periodic
        // spline.
        std::vector<double> values_at_psi(extended_resolution * resolution,
                                          0.0);
        std::vector<double> dpsi_at_psi(extended_resolution * resolution, 0.0);
        std::vector<double> line(extended_resolution, 0.0);
        std::vector<double> second_derivatives;
        for (std::size_t phi = 0; phi < extended_resolution; phi++)
        {
            for (std::size_t psi = 0; psi < extended_resolution; psi++)
            {
                line[psi] = extended_grid[phi * extended_resolution + psi];
            }
            CMap_Natural_Spline_Second_Derivatives(line, &second_derivatives);
            for (std::size_t psi = 0; psi < resolution; psi++)
            {
                const double coordinate = static_cast<double>(psi) +
                                          0.5 * static_cast<double>(resolution);
                CMap_Natural_Spline_Interpolate(
                    line, second_derivatives, coordinate,
                    &values_at_psi[phi * resolution + psi],
                    &dpsi_at_psi[phi * resolution + psi]);
            }
        }

        std::vector<double> dphi(grid_size, 0.0);
        std::vector<double> dpsi(grid_size, 0.0);
        std::vector<double> dphi_dpsi(grid_size, 0.0);
        for (std::size_t psi = 0; psi < resolution; psi++)
        {
            for (std::size_t phi = 0; phi < extended_resolution; phi++)
            {
                line[phi] = values_at_psi[phi * resolution + psi];
            }
            CMap_Natural_Spline_Second_Derivatives(line, &second_derivatives);
            for (std::size_t phi = 0; phi < resolution; phi++)
            {
                const double coordinate = static_cast<double>(phi) +
                                          0.5 * static_cast<double>(resolution);
                double interpolated_value = 0.0;
                CMap_Natural_Spline_Interpolate(line, second_derivatives,
                                                coordinate, &interpolated_value,
                                                &dphi[phi * resolution + psi]);
            }

            for (std::size_t phi = 0; phi < extended_resolution; phi++)
            {
                line[phi] = dpsi_at_psi[phi * resolution + psi];
            }
            CMap_Natural_Spline_Second_Derivatives(line, &second_derivatives);
            for (std::size_t phi = 0; phi < resolution; phi++)
            {
                const double coordinate = static_cast<double>(phi) +
                                          0.5 * static_cast<double>(resolution);
                CMap_Natural_Spline_Interpolate(
                    line, second_derivatives, coordinate,
                    &dpsi[phi * resolution + psi],
                    &dphi_dpsi[phi * resolution + psi]);
            }
        }

        for (std::size_t phi = 0; phi < resolution; phi++)
        {
            const std::size_t next_phi = (phi + 1) % resolution;
            for (std::size_t psi = 0; psi < resolution; psi++)
            {
                const std::size_t next_psi = (psi + 1) % resolution;
                const std::size_t p00 = phi * resolution + psi;
                const std::size_t p10 = next_phi * resolution + psi;
                const std::size_t p01 = phi * resolution + next_psi;
                const std::size_t p11 = next_phi * resolution + next_psi;
                const double parameters[16] = {
                    values[p00],    values[p10],    values[p01],
                    values[p11],    dphi[p00],      dphi[p10],
                    dphi[p01],      dphi[p11],      dpsi[p00],
                    dpsi[p10],      dpsi[p01],      dpsi[p11],
                    dphi_dpsi[p00], dphi_dpsi[p10], dphi_dpsi[p01],
                    dphi_dpsi[p11],
                };
                const std::size_t coefficient_base =
                    static_cast<std::size_t>(type_offset[type]) + 16 * p00;
                for (int coefficient = 0; coefficient < 16; coefficient++)
                {
                    double value = 0.0;
                    for (int parameter = 0; parameter < 16; parameter++)
                    {
                        value += A_inv[coefficient][parameter] *
                                 parameters[parameter];
                    }
                    const float stored_value = static_cast<float>(value);
                    const double float_max =
                        static_cast<double>(std::numeric_limits<float>::max());
                    if (!Double_Memory_Is_Finite(&value) || value > float_max ||
                        value < -float_max ||
                        !CMap_Float_Is_Finite(stored_value) ||
                        (value != 0.0 && stored_value == 0.0f))
                    {
                        controller->Throw_Formatted_SPONGE_Error(
                            spongeErrorBadFileFormat, "CMAP::Interpolation",
                            "Reason:\n\tCMAP type %d produces an "
                            "interpolation coefficient outside the finite "
                            "float range\n",
                            type);
                    }
                    if (!CMap_Float_Is_Zero_Or_Normal(stored_value))
                    {
                        controller->Throw_Formatted_SPONGE_Error(
                            spongeErrorBadFileFormat, "CMAP::Interpolation",
                            "Reason:\n\tCMAP type %d produces a subnormal "
                            "interpolation coefficient; derived CMAP values "
                            "must be finite zero or normal floats for "
                            "consistent FTZ behavior\n",
                            type);
                    }
                    h_inter_coeff[coefficient_base + coefficient] =
                        stored_value;
                }
            }
        }
        grid_offset += grid_size;
    }
    controller->printf("    End Interpolating CMAP Grid Value\n");
}

static __device__ __forceinline__ void CMap_Fail_Invalid_Geometry(
    int global_term, int* invalid_geometry_term,
    bool accumulator_failure = false)
{
#ifdef GPU_ARCH_NAME
    if (accumulator_failure)
    {
        printf(
            "Fatal SPONGE CMAP error: global term %d would produce a "
            "non-finite force/energy/virial accumulator.\n",
            global_term);
    }
    else
    {
        printf(
            "Fatal SPONGE CMAP error: global term %d has undefined or "
            "non-finite torsion geometry/energy/force/virial.\n",
            global_term);
    }
#if defined(USE_CUDA)
    asm volatile("trap;");
#elif defined(USE_HIP)
    __builtin_trap();
#endif
#else
    // -1 means no failure.  Encode accumulator failures as -term-2 so the
    // host can report the actual cause without another device-side buffer.
    atomicExch(invalid_geometry_term,
               accumulator_failure ? -global_term - 2 : global_term);
#endif
}

static __device__ __forceinline__ bool CMap_Runtime_Float_Is_Finite(float value)
{
#ifdef GPU_ARCH_NAME
    // Inspect the representation directly so CUDA/HIP fast-math cannot fold
    // validation of an arithmetic result away under a finite-value assumption.
    return (__float_as_uint(value) & 0x7f800000U) != 0x7f800000U;
#elif defined(__GNUC__) || defined(__clang__)
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value),
                  "SPONGE requires 32-bit IEEE-754 floats");
    std::memcpy(&bits, &value, sizeof(value));
    // Keep the hot-path test as an integer mask while preventing fast-math
    // from assuming away validation of intermediate results.
    __asm__ __volatile__("" : "+r"(bits));
    return (bits & 0x7f800000U) != 0x7f800000U;
#else
    return Float_Memory_Is_Finite(&value);
#endif
}

static __device__ __forceinline__ bool CMap_Finite_Atomic_Add(
    VECTOR* accumulator, const VECTOR& value)
{
    return Finite_Atomic_Add(&accumulator->x, value.x) &&
           Finite_Atomic_Add(&accumulator->y, value.y) &&
           Finite_Atomic_Add(&accumulator->z, value.z);
}

static __device__ __forceinline__ bool CMap_Finite_Atomic_Add(
    LTMatrix3* accumulator, const LTMatrix3& value)
{
    return Finite_Atomic_Add(&accumulator->a11, value.a11) &&
           Finite_Atomic_Add(&accumulator->a21, value.a21) &&
           Finite_Atomic_Add(&accumulator->a22, value.a22) &&
           Finite_Atomic_Add(&accumulator->a31, value.a31) &&
           Finite_Atomic_Add(&accumulator->a32, value.a32) &&
           Finite_Atomic_Add(&accumulator->a33, value.a33);
}

static __global__
    __launch_bounds__(1024) void CMAP_Force_With_Atom_Energy_And_Virial_Device(
        const int cmap_numbers, const int local_atom_numbers, const VECTOR* crd,
        const LTMatrix3 cell, const LTMatrix3 rcell, const int* atom_a,
        const int* atom_b, const int* atom_c, const int* atom_d,
        const int* atom_e, const int* cmap_type, const int* cmap_global_index,
        const int* resolution, const int* type_offset,
        const float* all_inter_coeff, VECTOR* frc, int need_potential,
        float* ene, float* cmap_ene, int need_pressure, LTMatrix3* virial,
        int* invalid_geometry_term)
{
#ifdef USE_GPU
    int cmap_i = blockDim.x * blockIdx.x + threadIdx.x;
    if (cmap_i < cmap_numbers)
#else
#pragma omp parallel for
    for (int cmap_i = 0; cmap_i < cmap_numbers; cmap_i++)
#endif
    {
        int atom_i = atom_a[cmap_i];
        int atom_j = atom_b[cmap_i];
        int atom_k = atom_c[cmap_i];
        int atom_l = atom_d[cmap_i];
        int atom_m = atom_e[cmap_i];
        VECTOR ri = crd[atom_i];
        VECTOR rj = crd[atom_j];
        VECTOR rk = crd[atom_k];
        VECTOR rl = crd[atom_l];
        VECTOR rm = crd[atom_m];
        // 计算phi
        VECTOR drij = Get_Periodic_Displacement(ri, rj, cell, rcell);
        VECTOR drkj = Get_Periodic_Displacement(rk, rj, cell, rcell);
        VECTOR drkl = Get_Periodic_Displacement(rk, rl, cell, rcell);
        // Reuse each minimum-image bond in reverse so both torsions share one
        // continuous five-atom image, including exact half-box ties.
        VECTOR drjk = -drkj;
        VECTOR drlk = -drkl;
        VECTOR drlm = Get_Periodic_Displacement(rl, rm, cell, rcell);

        // 法向量夹角。二面角在中心键或任一平面法向量为零时
        // 没有定义；在求倒数之前检测，避免让 NaN/Inf 进入力累加。
        VECTOR r1_phi = drij ^ drkj;
        VECTOR r2_phi = drkl ^ drkj;
        VECTOR r1_psi = drjk ^ drlk;
        VECTOR r2_psi = drlm ^ drlk;

        const float r1_phi_squared = r1_phi * r1_phi;
        const float r2_phi_squared = r2_phi * r2_phi;
        const float phi_bond_squared = drkj * drkj;
        const float r1_psi_squared = r1_psi * r1_psi;
        const float r2_psi_squared = r2_psi * r2_psi;
        const float psi_bond_squared = drlk * drlk;
        bool geometry_is_valid =
            r1_phi_squared > 0.0f && r2_phi_squared > 0.0f &&
            phi_bond_squared > 0.0f && r1_psi_squared > 0.0f &&
            r2_psi_squared > 0.0f && psi_bond_squared > 0.0f &&
            CMap_Runtime_Float_Is_Finite(r1_phi_squared) &&
            CMap_Runtime_Float_Is_Finite(r2_phi_squared) &&
            CMap_Runtime_Float_Is_Finite(phi_bond_squared) &&
            CMap_Runtime_Float_Is_Finite(r1_psi_squared) &&
            CMap_Runtime_Float_Is_Finite(r2_psi_squared) &&
            CMap_Runtime_Float_Is_Finite(psi_bond_squared);

        float r1_1_phi = 0.0f;
        float r2_1_phi = 0.0f;
        float rkj_1 = 0.0f;
        float r1_1_psi = 0.0f;
        float r2_1_psi = 0.0f;
        float rlk_1 = 0.0f;
        float inverse_r1_phi_squared = 0.0f;
        float inverse_r2_phi_squared = 0.0f;
        float inverse_r1_psi_squared = 0.0f;
        float inverse_r2_psi_squared = 0.0f;
        float phi_bond_length = 0.0f;
        float psi_bond_length = 0.0f;
        float cos_phi = 0.0f;
        float sin_phi = 0.0f;
        float cos_psi = 0.0f;
        float sin_psi = 0.0f;
        if (geometry_is_valid)
        {
            r1_1_phi = rnorm3df(r1_phi.x, r1_phi.y, r1_phi.z);
            r2_1_phi = rnorm3df(r2_phi.x, r2_phi.y, r2_phi.z);
            rkj_1 = rnorm3df(drkj.x, drkj.y, drkj.z);
            r1_1_psi = rnorm3df(r1_psi.x, r1_psi.y, r1_psi.z);
            r2_1_psi = rnorm3df(r2_psi.x, r2_psi.y, r2_psi.z);
            rlk_1 = rnorm3df(drlk.x, drlk.y, drlk.z);
            inverse_r1_phi_squared = 1.0f / r1_phi_squared;
            inverse_r2_phi_squared = 1.0f / r2_phi_squared;
            inverse_r1_psi_squared = 1.0f / r1_psi_squared;
            inverse_r2_psi_squared = 1.0f / r2_psi_squared;
            phi_bond_length = sqrtf(phi_bond_squared);
            psi_bond_length = sqrtf(psi_bond_squared);
            geometry_is_valid =
                CMap_Runtime_Float_Is_Finite(r1_1_phi) &&
                CMap_Runtime_Float_Is_Finite(r2_1_phi) &&
                CMap_Runtime_Float_Is_Finite(rkj_1) &&
                CMap_Runtime_Float_Is_Finite(r1_1_psi) &&
                CMap_Runtime_Float_Is_Finite(r2_1_psi) &&
                CMap_Runtime_Float_Is_Finite(rlk_1) &&
                CMap_Runtime_Float_Is_Finite(inverse_r1_phi_squared) &&
                CMap_Runtime_Float_Is_Finite(inverse_r2_phi_squared) &&
                CMap_Runtime_Float_Is_Finite(inverse_r1_psi_squared) &&
                CMap_Runtime_Float_Is_Finite(inverse_r2_psi_squared) &&
                CMap_Runtime_Float_Is_Finite(phi_bond_length) &&
                CMap_Runtime_Float_Is_Finite(psi_bond_length);
        }
        if (geometry_is_valid)
        {
            const float r1_1_r2_1_phi = r1_1_phi * r2_1_phi;
            const float r1_1_r2_1_psi = r1_1_psi * r2_1_psi;
            cos_phi = -(r1_phi * r2_phi) * r1_1_r2_1_phi;
            sin_phi = ((r2_phi ^ r1_phi) * drkj) * r1_1_r2_1_phi * rkj_1;
            cos_psi = -(r1_psi * r2_psi) * r1_1_r2_1_psi;
            sin_psi = ((r2_psi ^ r1_psi) * drlk) * r1_1_r2_1_psi * rlk_1;
            geometry_is_valid = CMap_Runtime_Float_Is_Finite(cos_phi) &&
                                CMap_Runtime_Float_Is_Finite(sin_phi) &&
                                CMap_Runtime_Float_Is_Finite(cos_psi) &&
                                CMap_Runtime_Float_Is_Finite(sin_psi) &&
                                (cos_phi != 0.0f || sin_phi != 0.0f) &&
                                (cos_psi != 0.0f || sin_psi != 0.0f);
        }
        if (!geometry_is_valid)
        {
            CMap_Fail_Invalid_Geometry(cmap_global_index[cmap_i],
                                       invalid_geometry_term);
            if (need_potential) cmap_ene[cmap_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }
        const float phi = atan2f(sin_phi, cos_phi);
        const float psi = atan2f(sin_psi, cos_psi);

        // Map [-pi, pi] to the periodic grid whose first knot is -pi.
        int temp_reso = resolution[cmap_type[cmap_i]];
        const float normalized_phi =
            (phi + CONSTANT_Pi) / (2.0f * CONSTANT_Pi / temp_reso);
        const float normalized_psi =
            (psi + CONSTANT_Pi) / (2.0f * CONSTANT_Pi / temp_reso);
        const float floor_phi = floorf(normalized_phi);
        const float floor_psi = floorf(normalized_psi);
        const float parm_phi = normalized_phi - floor_phi;
        const float parm_psi = normalized_psi - floor_psi;
        int locate_phi = static_cast<int>(floor_phi) % temp_reso;
        int locate_psi = static_cast<int>(floor_psi) % temp_reso;
        if (locate_phi < 0) locate_phi += temp_reso;
        if (locate_psi < 0) locate_psi += temp_reso;

        // 定义幂次
        float parm_phi_2 = parm_phi * parm_phi;
        float parm_phi_3 = parm_phi_2 * parm_phi;
        float parm_psi_2 = parm_psi * parm_psi;
        float parm_psi_3 = parm_psi_2 * parm_psi;

        // 用于定位的中间变量
        const float* inter_coeff =
            all_inter_coeff + type_offset[cmap_type[cmap_i]];
        int locate = 16 * (locate_phi * temp_reso + locate_psi);

        // 计算能量对有符号归一化二面角的偏微分
        float dE_dphi =
            (inter_coeff[locate + 4] + parm_psi * inter_coeff[locate + 5] +
             parm_psi_2 * inter_coeff[locate + 6] +
             parm_psi_3 * inter_coeff[locate + 7]) +
            2 * parm_phi *
                (inter_coeff[locate + 8] + parm_psi * inter_coeff[locate + 9] +
                 parm_psi_2 * inter_coeff[locate + 10] +
                 parm_psi_3 * inter_coeff[locate + 11]) +
            3 * parm_phi_2 *
                (inter_coeff[locate + 12] +
                 parm_psi * inter_coeff[locate + 13] +
                 parm_psi_2 * inter_coeff[locate + 14] +
                 parm_psi_3 * inter_coeff[locate + 15]);

        float dE_dpsi =
            inter_coeff[locate + 1] + 2 * parm_psi * inter_coeff[locate + 2] +
            3 * parm_psi_2 * inter_coeff[locate + 3] +
            parm_phi * (inter_coeff[locate + 5] +
                        2 * parm_psi * inter_coeff[locate + 6] +
                        3 * parm_psi_2 * inter_coeff[locate + 7]) +
            parm_phi_2 * (inter_coeff[locate + 9] +
                          2 * parm_psi * inter_coeff[locate + 10] +
                          3 * parm_psi_2 * inter_coeff[locate + 11]) +
            parm_phi_3 * (inter_coeff[locate + 13] +
                          2 * parm_psi * inter_coeff[locate + 14] +
                          3 * parm_psi_2 * inter_coeff[locate + 15]);

        // 将有符号归一化二面角映射回弧度制二面角
        dE_dphi = dE_dphi / (2.0 * CONSTANT_Pi / temp_reso);
        dE_dpsi = dE_dpsi / (2.0 * CONSTANT_Pi / temp_reso);

        float Energy = 0.0f;
        if (need_potential)
        {
            Energy = inter_coeff[locate] + parm_psi * inter_coeff[locate + 1] +
                     parm_psi_2 * inter_coeff[locate + 2] +
                     parm_psi_3 * inter_coeff[locate + 3] +
                     parm_phi * (inter_coeff[locate + 4] +
                                 parm_psi * inter_coeff[locate + 5] +
                                 parm_psi_2 * inter_coeff[locate + 6] +
                                 parm_psi_3 * inter_coeff[locate + 7]) +
                     parm_phi_2 * (inter_coeff[locate + 8] +
                                   parm_psi * inter_coeff[locate + 9] +
                                   parm_psi_2 * inter_coeff[locate + 10] +
                                   parm_psi_3 * inter_coeff[locate + 11]) +
                     parm_phi_3 * (inter_coeff[locate + 12] +
                                   parm_psi * inter_coeff[locate + 13] +
                                   parm_psi_2 * inter_coeff[locate + 14] +
                                   parm_psi_3 * inter_coeff[locate + 15]);
        }
        if (!CMap_Runtime_Float_Is_Finite(dE_dphi) ||
            !CMap_Runtime_Float_Is_Finite(dE_dpsi) ||
            (need_potential && !CMap_Runtime_Float_Is_Finite(Energy)))
        {
            CMap_Fail_Invalid_Geometry(cmap_global_index[cmap_i],
                                       invalid_geometry_term);
            if (need_potential) cmap_ene[cmap_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        // phi角部分
        VECTOR temp_phi_A = drij ^ drjk;
        VECTOR temp_phi_B = drlk ^ drjk;

        VECTOR dphi_dri =
            -phi_bond_length * inverse_r1_phi_squared * temp_phi_A;
        VECTOR dphi_drj =
            +phi_bond_length * inverse_r1_phi_squared * temp_phi_A +
            (drij * drjk) * inverse_r1_phi_squared * rkj_1 * temp_phi_A -
            (drlk * drjk) * inverse_r2_phi_squared * rkj_1 * temp_phi_B;
        VECTOR dphi_drk =
            -phi_bond_length * inverse_r2_phi_squared * temp_phi_B -
            (drij * drjk) * inverse_r1_phi_squared * rkj_1 * temp_phi_A +
            (drlk * drjk) * inverse_r2_phi_squared * rkj_1 * temp_phi_B;
        VECTOR dphi_drl =
            +phi_bond_length * inverse_r2_phi_squared * temp_phi_B;
        VECTOR dphi_drm = {0, 0, 0};

        // psi角部分
        VECTOR drml = -drlm;

        VECTOR temp_psi_A = drjk ^ drkl;
        VECTOR temp_psi_B = drml ^ drkl;

        VECTOR dpsi_dri = {0, 0, 0};
        VECTOR dpsi_drj =
            -psi_bond_length * inverse_r1_psi_squared * temp_psi_A;
        VECTOR dpsi_drk =
            psi_bond_length * inverse_r1_psi_squared * temp_psi_A +
            (drjk * drkl) * inverse_r1_psi_squared * rlk_1 * temp_psi_A -
            (drml * drkl) * inverse_r2_psi_squared * rlk_1 * temp_psi_B;
        VECTOR dpsi_drl =
            -psi_bond_length * inverse_r2_psi_squared * temp_psi_B -
            (drjk * drkl) * inverse_r1_psi_squared * rlk_1 * temp_psi_A +
            (drml * drkl) * inverse_r2_psi_squared * rlk_1 * temp_psi_B;
        VECTOR dpsi_drm = psi_bond_length * inverse_r2_psi_squared * temp_psi_B;

        // 计算力
        VECTOR fi = -(dE_dphi * dphi_dri + dE_dpsi * dpsi_dri);
        VECTOR fj = -(dE_dphi * dphi_drj + dE_dpsi * dpsi_drj);
        VECTOR fk = -(dE_dphi * dphi_drk + dE_dpsi * dpsi_drk);
        VECTOR fl = -(dE_dphi * dphi_drl + dE_dpsi * dpsi_drl);
        VECTOR fm = -(dE_dphi * dphi_drm + dE_dpsi * dpsi_drm);

        if (!CMap_Runtime_Float_Is_Finite(fi.x) ||
            !CMap_Runtime_Float_Is_Finite(fi.y) ||
            !CMap_Runtime_Float_Is_Finite(fi.z) ||
            !CMap_Runtime_Float_Is_Finite(fj.x) ||
            !CMap_Runtime_Float_Is_Finite(fj.y) ||
            !CMap_Runtime_Float_Is_Finite(fj.z) ||
            !CMap_Runtime_Float_Is_Finite(fk.x) ||
            !CMap_Runtime_Float_Is_Finite(fk.y) ||
            !CMap_Runtime_Float_Is_Finite(fk.z) ||
            !CMap_Runtime_Float_Is_Finite(fl.x) ||
            !CMap_Runtime_Float_Is_Finite(fl.y) ||
            !CMap_Runtime_Float_Is_Finite(fl.z) ||
            !CMap_Runtime_Float_Is_Finite(fm.x) ||
            !CMap_Runtime_Float_Is_Finite(fm.y) ||
            !CMap_Runtime_Float_Is_Finite(fm.z))
        {
            CMap_Fail_Invalid_Geometry(cmap_global_index[cmap_i],
                                       invalid_geometry_term);
            if (need_potential) cmap_ene[cmap_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }

        LTMatrix3 term_virial = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        if (need_pressure && atom_i < local_atom_numbers)
        {
            const VECTOR ri_unwrapped = drij;
            const VECTOR rj_unwrapped = {0.0f, 0.0f, 0.0f};
            const VECTOR rk_unwrapped = drkj;
            const VECTOR rl_unwrapped = drkj - drkl;
            const VECTOR rm_unwrapped = rl_unwrapped + drml;
            term_virial = Get_Virial_From_Force_Dis(fi, ri_unwrapped) +
                          Get_Virial_From_Force_Dis(fj, rj_unwrapped) +
                          Get_Virial_From_Force_Dis(fk, rk_unwrapped) +
                          Get_Virial_From_Force_Dis(fl, rl_unwrapped) +
                          Get_Virial_From_Force_Dis(fm, rm_unwrapped);
            if (!CMap_Runtime_Float_Is_Finite(term_virial.a11) ||
                !CMap_Runtime_Float_Is_Finite(term_virial.a21) ||
                !CMap_Runtime_Float_Is_Finite(term_virial.a22) ||
                !CMap_Runtime_Float_Is_Finite(term_virial.a31) ||
                !CMap_Runtime_Float_Is_Finite(term_virial.a32) ||
                !CMap_Runtime_Float_Is_Finite(term_virial.a33))
            {
                CMap_Fail_Invalid_Geometry(cmap_global_index[cmap_i],
                                           invalid_geometry_term);
                if (need_potential) cmap_ene[cmap_i] = 0.0f;
#ifdef USE_GPU
                return;
#else
                continue;
#endif
            }
        }

        bool accumulation_is_valid = true;
        if (atom_i < local_atom_numbers)
        {
            accumulation_is_valid = CMap_Finite_Atomic_Add(frc + atom_i, fi);
        }
        if (accumulation_is_valid && atom_j < local_atom_numbers)
        {
            accumulation_is_valid = CMap_Finite_Atomic_Add(frc + atom_j, fj);
        }
        if (accumulation_is_valid && atom_k < local_atom_numbers)
        {
            accumulation_is_valid = CMap_Finite_Atomic_Add(frc + atom_k, fk);
        }
        if (accumulation_is_valid && atom_l < local_atom_numbers)
        {
            accumulation_is_valid = CMap_Finite_Atomic_Add(frc + atom_l, fl);
        }
        if (accumulation_is_valid && atom_m < local_atom_numbers)
        {
            accumulation_is_valid = CMap_Finite_Atomic_Add(frc + atom_m, fm);
        }
        if (accumulation_is_valid && need_potential &&
            atom_i < local_atom_numbers)
        {
            accumulation_is_valid = Finite_Atomic_Add(ene + atom_i, Energy);
        }
        if (accumulation_is_valid && need_pressure &&
            atom_i < local_atom_numbers)
        {
            accumulation_is_valid =
                CMap_Finite_Atomic_Add(virial + atom_i, term_virial);
        }
        if (!accumulation_is_valid)
        {
            CMap_Fail_Invalid_Geometry(cmap_global_index[cmap_i],
                                       invalid_geometry_term, true);
            if (need_potential) cmap_ene[cmap_i] = 0.0f;
#ifdef USE_GPU
            return;
#else
            continue;
#endif
        }
        if (need_potential)
        {
            cmap_ene[cmap_i] = atom_i < local_atom_numbers ? Energy : 0.0f;
        }
    }
}

void CMAP::CMAP_Force_With_Atom_Energy_And_Virial(
    const VECTOR* crd, const LTMatrix3 cell, const LTMatrix3 rcell, VECTOR* frc,
    int need_potential, float* atom_energy, int need_pressure,
    LTMatrix3* atom_virial)
{
    if (is_initialized && num_cmap_local > 0)
    {
#ifndef GPU_ARCH_NAME
        deviceMemset(d_invalid_geometry_term, -1, sizeof(int));
#endif
        Launch_Device_Kernel(
            CMAP_Force_With_Atom_Energy_And_Virial_Device,
            (num_cmap_local - 1) / CONTROLLER::device_max_thread + 1,
            CONTROLLER::device_max_thread, 0, NULL, this->num_cmap_local,
            this->local_atom_numbers, crd, cell, rcell, this->d_atom_a_local,
            this->d_atom_b_local, this->d_atom_c_local, this->d_atom_d_local,
            this->d_atom_e_local, this->d_cmap_type_local,
            this->d_cmap_global_index_local, this->d_cmap_resolution,
            this->d_type_offset, this->d_inter_coeff, frc, need_potential,
            atom_energy, d_cmap_ene, need_pressure, atom_virial,
            d_invalid_geometry_term);
#ifndef GPU_ARCH_NAME
        int invalid_code = -1;
        deviceMemcpy(&invalid_code, d_invalid_geometry_term, sizeof(int),
                     deviceMemcpyDeviceToHost);
        if (invalid_code != -1)
        {
            const bool accumulator_failure = invalid_code < -1;
            const int invalid_term =
                accumulator_failure ? -(invalid_code + 2) : invalid_code;
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "CMAP::CMAP_Force_With_Atom_Energy_And_Virial",
                "Reason:\n\t%s CMAP term %d (global atoms %d %d %d %d %d) "
                "%s\n",
                module_name, invalid_term, h_atom_a[invalid_term],
                h_atom_b[invalid_term], h_atom_c[invalid_term],
                h_atom_d[invalid_term], h_atom_e[invalid_term],
                accumulator_failure
                    ? "would produce a non-finite force/energy/virial "
                      "accumulator"
                    : "has undefined or non-finite torsion "
                      "geometry/energy/force/virial");
            return;
        }
#endif
    }
}

void CMAP::Step_Print(CONTROLLER* controller, bool print_sum)
{
    if (is_initialized && CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        Sum_Of_List(d_cmap_ene, d_sigma_of_cmap_ene, num_cmap_local);
        deviceMemcpy(&h_sigma_of_cmap_ene, d_sigma_of_cmap_ene, sizeof(float),
                     deviceMemcpyDeviceToHost);
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &h_sigma_of_cmap_ene, 1, MPI_FLOAT, MPI_SUM,
                      CONTROLLER::pp_comm);
#endif
        if (!CMap_Float_Is_Finite(h_sigma_of_cmap_ene))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown, "CMAP::Step_Print",
                "Reason:\n\t%s CMAP total energy is non-finite after "
                "local/MPI reduction\n",
                module_name);
            return;
        }
        controller->Step_Print(this->module_name, h_sigma_of_cmap_ene,
                               print_sum);
    }
}

static __global__ void get_local_device(
    int tot_cmap_num, int local_coordinate_numbers, const int* d_atom_a,
    const int* d_atom_b, const int* d_atom_c, const int* d_atom_d,
    const int* d_atom_e, const int* d_cmap_type, const char* atom_local_label,
    const int* atom_local_id, int* d_atom_a_local, int* d_atom_b_local,
    int* d_atom_c_local, int* d_atom_d_local, int* d_atom_e_local,
    int* d_cmap_type_local, int* d_cmap_global_index_local,
    int* d_num_cmap_local, int* d_invalid_local_term, int* d_invalid_local_atom)
{
#ifdef USE_GPU
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx != 0) return;
#endif
    d_num_cmap_local[0] = 0;
    d_invalid_local_term[0] = -1;
    d_invalid_local_atom[0] = -1;
    for (int i = 0; i < tot_cmap_num; i++)
    {
        const int global_atoms[5] = {d_atom_a[i], d_atom_b[i], d_atom_c[i],
                                     d_atom_d[i], d_atom_e[i]};
        if (atom_local_label[global_atoms[0]] == 1 ||
            atom_local_label[global_atoms[1]] == 1 ||
            atom_local_label[global_atoms[2]] == 1 ||
            atom_local_label[global_atoms[3]] == 1 ||
            atom_local_label[global_atoms[4]] == 1)
        {
            int local_atoms[5];
            for (int atom = 0; atom < 5; atom++)
            {
                local_atoms[atom] = atom_local_id[global_atoms[atom]];
                if (local_atoms[atom] < 0 ||
                    local_atoms[atom] >= local_coordinate_numbers)
                {
                    d_invalid_local_term[0] = i;
                    d_invalid_local_atom[0] = global_atoms[atom];
                    return;
                }
            }
            const int local_term = d_num_cmap_local[0];
            d_atom_a_local[local_term] = local_atoms[0];
            d_atom_b_local[local_term] = local_atoms[1];
            d_atom_c_local[local_term] = local_atoms[2];
            d_atom_d_local[local_term] = local_atoms[3];
            d_atom_e_local[local_term] = local_atoms[4];
            d_cmap_type_local[local_term] = d_cmap_type[i];
            d_cmap_global_index_local[local_term] = i;
            d_num_cmap_local[0]++;
        }
    }
}

void CMAP::Get_Local(int* atom_local, int local_atom_numbers, int ghost_numbers,
                     char* atom_local_label, int* atom_local_id)
{
    if (!is_initialized) return;
    if (local_atom_numbers < 0 || ghost_numbers < 0 ||
        local_atom_numbers > std::numeric_limits<int>::max() - ghost_numbers)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "CMAP::Get_Local",
            "Reason:\n\t%s received invalid local/ghost atom counts %d/%d\n",
            module_name, local_atom_numbers, ghost_numbers);
        return;
    }
    const int local_coordinate_numbers = local_atom_numbers + ghost_numbers;
    num_cmap_local = 0;
    this->local_atom_numbers = local_atom_numbers;
    Launch_Device_Kernel(
        get_local_device, 1, 1, 0, NULL, this->tot_cmap_num,
        local_coordinate_numbers, this->d_atom_a, this->d_atom_b,
        this->d_atom_c, this->d_atom_d, this->d_atom_e, this->d_cmap_type,
        atom_local_label, atom_local_id, this->d_atom_a_local,
        this->d_atom_b_local, this->d_atom_c_local, this->d_atom_d_local,
        this->d_atom_e_local, this->d_cmap_type_local,
        this->d_cmap_global_index_local, this->d_num_cmap_local,
        this->d_invalid_local_term, this->d_invalid_local_atom);
    int invalid_term = -1;
    int invalid_atom = -1;
    deviceMemcpy(&invalid_term, d_invalid_local_term, sizeof(int),
                 deviceMemcpyDeviceToHost);
    deviceMemcpy(&invalid_atom, d_invalid_local_atom, sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (invalid_term >= 0)
    {
        int invalid_local_id = -1;
        deviceMemcpy(&invalid_local_id, atom_local_id + invalid_atom,
                     sizeof(int), deviceMemcpyDeviceToHost);
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown, "CMAP::Get_Local",
            "Reason:\n\t%s CMAP term %d (global atoms %d %d %d %d %d) "
            "maps global atom %d to local index %d outside the valid "
            "owned/ghost range [0, %d) on this domain\n",
            module_name, invalid_term, h_atom_a[invalid_term],
            h_atom_b[invalid_term], h_atom_c[invalid_term],
            h_atom_d[invalid_term], h_atom_e[invalid_term], invalid_atom,
            invalid_local_id, local_coordinate_numbers);
        return;
    }
    deviceMemcpy(&num_cmap_local, d_num_cmap_local, sizeof(int),
                 deviceMemcpyDeviceToHost);
}
