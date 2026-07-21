#pragma once

#include "becke.hpp"
#include "dft.hpp"
#include "lebedev.hpp"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

struct QC_GRID_POINT
{
    float x;
    float y;
    float z;
    float w;
};

static inline bool QC_Grid_Double_Is_Finite(double value)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x7ff0000000000000ULL) != 0x7ff0000000000000ULL;
}

static inline bool QC_Grid_Float_Is_Finite(float value)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return (bits & 0x7f800000U) != 0x7f800000U;
}

static inline double QC_Get_Covalent_Radius_Bohr(int Z)
{
    return sponge_qc_elements::Cordero_Covalent_Radius_Bohr(Z);
}

static void QC_Gauss_Legendre_01(int n, std::vector<double>& nodes,
                                 std::vector<double>& weights)
{
    nodes.assign(n, 0.0f);
    weights.assign(n, 0.0f);
    const double eps = 1e-14;
    const int m = (n + 1) / 2;
    for (int i = 0; i < m; i++)
    {
        const double i1 = (double)i + 1.0;
        double z = cos(CONSTANT_Pi * (i1 - 0.25) / ((double)n + 0.5));
        double z1 = 0.0;
        double p1 = 0.0, p2 = 0.0, pp = 0.0;
        do
        {
            p1 = 1.0;
            p2 = 0.0;
            for (int j = 1; j <= n; j++)
            {
                double p3 = p2;
                p2 = p1;
                p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / (double)j;
            }
            pp = (double)n * (z * p1 - p2) / (z * z - 1.0);
            z1 = z;
            z = z1 - p1 / pp;
        } while (fabs(z - z1) > eps);
        const double x_low = -z;
        const double x_high = z;
        const double w = 2.0 / ((1.0 - z * z) * pp * pp);
        // Map [-1,1] -> [0,1]
        nodes[i] = 0.5 * (x_low + 1.0);
        nodes[n - 1 - i] = 0.5 * (x_high + 1.0);
        weights[i] = 0.5 * w;
        weights[n - 1 - i] = 0.5 * w;
    }
}

static void QC_Build_Fibonacci_Angular_Grid(int n_ang, std::vector<float>& dirs,
                                            std::vector<float>& w_ang)
{
    if (sponge_qc_lebedev::Load_Lebedev_Angular_Grid(n_ang, dirs, w_ang))
        return;
    dirs.assign(static_cast<std::size_t>(n_ang) * 3, 0.0f);
    w_ang.assign(static_cast<std::size_t>(n_ang), 0.0f);
    const float golden = (sqrtf(5.0f) - 1.0f) * 0.5f;
    const float w = 4.0f * CONSTANT_Pi / (float)n_ang;
    for (int i = 0; i < n_ang; i++)
    {
        const float z = 1.0f - 2.0f * ((float)i + 0.5f) / (float)n_ang;
        const float r = sqrtf(fmaxf(0.0f, 1.0f - z * z));
        const float phi = 2.0f * CONSTANT_Pi *
                          ((float)i * golden - floorf((float)i * golden));
        dirs[(int)i * 3 + 0] = r * cosf(phi);
        dirs[(int)i * 3 + 1] = r * sinf(phi);
        dirs[(int)i * 3 + 2] = z;
        w_ang[i] = w;
    }
}

static inline void QC_Get_Atom_Coord_From_Env(const QC_MOLECULE& mol,
                                              const float* env, int atom_i,
                                              float& x, float& y, float& z)
{
    const int ptr = mol.h_atm[(int)atom_i * 6 + 1];
    x = env[ptr + 0];
    y = env[ptr + 1];
    z = env[ptr + 2];
}

static double QC_Atom_Partition_Weight(
    const QC_MOLECULE& mol, const float* env,
    const std::vector<double>& covalent_radii, int atom_i, double x, double y,
    double z)
{
    const int natm = mol.natm;
    std::vector<double> log_p(static_cast<std::size_t>(natm), 0.0);
    // Do not encode an exactly-zero Becke product as -infinity.  The project
    // is compiled with finite-math optimizations, under which non-finite
    // sentinels and std::isfinite are not reliable control flow.
    std::vector<unsigned char> owns_point(static_cast<std::size_t>(natm), 1U);
    for (int a = 0; a < natm; a++)
    {
        float ax = 0.0f, ay = 0.0f, az = 0.0f;
        QC_Get_Atom_Coord_From_Env(mol, env, a, ax, ay, az);
        const double dax = x - static_cast<double>(ax);
        const double day = y - static_cast<double>(ay);
        const double daz = z - static_cast<double>(az);
        const double ra = sqrt(dax * dax + day * day + daz * daz);
        for (int b = 0; b < natm; b++)
        {
            if (a == b) continue;
            float bx = 0.0f, by = 0.0f, bz = 0.0f;
            QC_Get_Atom_Coord_From_Env(mol, env, b, bx, by, bz);
            const double dbx = x - static_cast<double>(bx);
            const double dby = y - static_cast<double>(by);
            const double dbz = z - static_cast<double>(bz);
            const double rb = sqrt(dbx * dbx + dby * dby + dbz * dbz);
            const double abx = static_cast<double>(ax) - bx;
            const double aby = static_cast<double>(ay) - by;
            const double abz = static_cast<double>(az) - bz;
            const double rab = sqrt(abx * abx + aby * aby + abz * abz);
            if (rab == 0.0)
                throw std::domain_error(
                    "DFT atom partition contains exactly coincident nuclei");
            const double mu = (ra - rb) / rab;
            const double adjusted_mu = QC_Becke_Size_Adjusted_Mu(
                mu, covalent_radii[static_cast<std::size_t>(a)],
                covalent_radii[static_cast<std::size_t>(b)]);
            const double shape = QC_Becke_Shape(adjusted_mu);
            if (shape == 0.0)
            {
                owns_point[static_cast<std::size_t>(a)] = 0U;
                break;
            }
            if (!(shape > 0.0) || !QC_Grid_Double_Is_Finite(shape))
                throw std::domain_error(
                    "DFT atom partition produced an invalid Becke factor");
            log_p[a] += log(shape);
            if (!QC_Grid_Double_Is_Finite(log_p[a]))
                throw std::domain_error(
                    "DFT atom partition produced a non-finite logarithmic "
                    "product");
        }
    }
    bool has_owner = false;
    double max_log = -std::numeric_limits<double>::max();
    for (int a = 0; a < natm; ++a)
    {
        if (owns_point[static_cast<std::size_t>(a)] == 0U) continue;
        has_owner = true;
        max_log = std::max(max_log, log_p[static_cast<std::size_t>(a)]);
    }
    if (!has_owner || !QC_Grid_Double_Is_Finite(max_log))
        throw std::domain_error(
            "DFT atom partition has no finite owning weight");
    double denominator = 0.0;
    for (int a = 0; a < natm; ++a)
        if (owns_point[static_cast<std::size_t>(a)] != 0U)
            denominator += exp(log_p[static_cast<std::size_t>(a)] - max_log);
    if (!(denominator > 0.0) || !QC_Grid_Double_Is_Finite(denominator))
        throw std::domain_error(
            "DFT atom partition normalization is not finite and positive");
    return owns_point[static_cast<std::size_t>(atom_i)] != 0U
               ? exp(log_p[atom_i] - max_log) / denominator
               : 0.0;
}

static std::vector<QC_GRID_POINT> QC_Build_Molecular_Grid(
    const QC_MOLECULE& mol, const float* env, const int* atomic_numbers,
    int atomic_number_stride, int n_radial, int n_angular)
{
    if (env == NULL || atomic_numbers == NULL || mol.natm <= 0) return {};

    std::vector<double> r_nodes, r_weights;
    QC_Gauss_Legendre_01(n_radial, r_nodes, r_weights);
    std::vector<float> dirs, w_ang;
    QC_Build_Fibonacci_Angular_Grid(n_angular, dirs, w_ang);
    std::vector<double> covalent_radii(static_cast<std::size_t>(mol.natm));
    for (int ia = 0; ia < mol.natm; ++ia)
        covalent_radii[static_cast<std::size_t>(ia)] =
            QC_Get_Covalent_Radius_Bohr(
                atomic_numbers[ia * atomic_number_stride]);

    std::vector<QC_GRID_POINT> grid;
    grid.reserve(static_cast<std::size_t>(mol.natm) *
                 static_cast<std::size_t>(n_radial) *
                 static_cast<std::size_t>(n_angular));
    for (int ia = 0; ia < mol.natm; ia++)
    {
        float cx = 0.0f, cy = 0.0f, cz = 0.0f;
        QC_Get_Atom_Coord_From_Env(mol, env, ia, cx, cy, cz);
        const double scale = covalent_radii[static_cast<std::size_t>(ia)];
        for (int ir = 0; ir < n_radial; ir++)
        {
            const double x = r_nodes[ir];
            const double one_m_x = 1.0 - x;
            if (!(one_m_x > 0.0))
                throw std::domain_error(
                    "DFT radial quadrature produced a node outside [0, 1)");
            const double r = scale * x / one_m_x;
            const double drdx = scale / (one_m_x * one_m_x);
            const double wr = r_weights[ir] * r * r * drdx;
            for (int iang = 0; iang < n_angular; iang++)
            {
                const double gx = cx + r * dirs[(int)iang * 3 + 0];
                const double gy = cy + r * dirs[(int)iang * 3 + 1];
                const double gz = cz + r * dirs[(int)iang * 3 + 2];
                const float gx_f = static_cast<float>(gx);
                const float gy_f = static_cast<float>(gy);
                const float gz_f = static_cast<float>(gz);
                const double w_part =
                    QC_Atom_Partition_Weight(
                        mol, env, covalent_radii, ia, (double)gx_f,
                        (double)gy_f, (double)gz_f);
                const double w = wr * w_ang[iang] * w_part;
                const float w_f = static_cast<float>(w);
                if (!QC_Grid_Float_Is_Finite(gx_f) ||
                    !QC_Grid_Float_Is_Finite(gy_f) ||
                    !QC_Grid_Float_Is_Finite(gz_f) ||
                    !QC_Grid_Float_Is_Finite(w_f) ||
                    !QC_Grid_Double_Is_Finite(w) || w < 0.0 ||
                    (w != 0.0 && w_f == 0.0f))
                {
                    throw std::overflow_error(
                        "DFT grid point or weight is not representable as "
                        "a finite float");
                }
                grid.push_back({gx_f, gy_f, gz_f, w_f});
            }
        }
    }
    return grid;
}

void QUANTUM_CHEMISTRY::Update_DFT_Grid()
{
    std::vector<float> env_host;
    std::vector<QC_GRID_POINT> grid;
    std::vector<float> staged_coords;
    std::vector<float> staged_weights;
    std::vector<double> staged_covalent_radii;
    try
    {
        env_host.resize(mol.h_env.size());
        if (mol.d_env != NULL)
            deviceMemcpy(env_host.data(), mol.d_env,
                         sizeof(float) * env_host.size(),
                         deviceMemcpyDeviceToHost);
        else
            std::copy(mol.h_env.begin(), mol.h_env.end(), env_host.begin());
        if (mol.h_atomic_numbers.size() !=
            static_cast<std::size_t>(mol.natm))
            throw std::logic_error(
                "DFT grid atomic identities do not match the molecule");
        grid = QC_Build_Molecular_Grid(
            mol, env_host.data(), mol.h_atomic_numbers.data(), 1,
            dft.dft_radial_points, dft.dft_angular_points);
        staged_covalent_radii.resize(static_cast<std::size_t>(mol.natm));
        for (int ia = 0; ia < mol.natm; ++ia)
            staged_covalent_radii[static_cast<std::size_t>(ia)] =
                QC_Get_Covalent_Radius_Bohr(mol.Atomic_Number(ia));
        if (grid.size() > static_cast<std::size_t>(dft.max_grid_capacity) ||
            grid.size() >
                static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            throw std::length_error(
                "generated grid exceeds its validated capacity");
        }
        staged_coords.resize(grid.size() * 3);
        staged_weights.resize(grid.size());
        for (std::size_t i = 0; i < grid.size(); i++)
        {
            staged_coords[i * 3 + 0] = grid[i].x;
            staged_coords[i * 3 + 1] = grid[i].y;
            staged_coords[i * 3 + 2] = grid[i].z;
            staged_weights[i] = grid[i].w;
        }
    }
    catch (const std::length_error& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "QUANTUM_CHEMISTRY::Update_DFT_Grid",
            "Reason:\n    DFT grid working storage is not representable "
            "(natm=%d, radial=%d, angular=%d): %s\n",
            mol.natm, dft.dft_radial_points, dft.dft_angular_points,
            error.what());
    }
    catch (const std::bad_alloc& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorMallocFailed, "QUANTUM_CHEMISTRY::Update_DFT_Grid",
            "Reason:\n    failed to allocate DFT grid working storage "
            "(natm=%d, radial=%d, angular=%d): %s\n",
            mol.natm, dft.dft_radial_points, dft.dft_angular_points,
            error.what());
    }
    catch (const std::exception& error)
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorSimulationBreakDown,
            "QUANTUM_CHEMISTRY::Update_DFT_Grid",
            "Reason:\n    failed to construct a valid DFT grid "
            "(natm=%d, radial=%d, angular=%d): %s\n",
            mol.natm, dft.dft_radial_points, dft.dft_angular_points,
            error.what());
    }

    deviceMemcpy(dft.d_covalent_radii, staged_covalent_radii.data(),
                 sizeof(double) * staged_covalent_radii.size(),
                 deviceMemcpyHostToDevice);
#ifdef CPU_ARCH_NAME
    // CPU "device" pointers alias their host storage.  Publish the staged
    // vectors first, then refresh both aliases; otherwise the old vectors are
    // destroyed at function exit and leave the DFT hot path with dangling
    // pointers.
    dft.h_grid_coords.swap(staged_coords);
    dft.h_grid_weights.swap(staged_weights);
    dft.d_grid_coords = dft.h_grid_coords.data();
    dft.d_grid_weights = dft.h_grid_weights.data();
#else
    deviceMemcpy(dft.d_grid_coords, staged_coords.data(),
                 sizeof(float) * staged_coords.size(),
                 deviceMemcpyHostToDevice);
    deviceMemcpy(dft.d_grid_weights, staged_weights.data(),
                 sizeof(float) * staged_weights.size(),
                 deviceMemcpyHostToDevice);
    dft.h_grid_coords.swap(staged_coords);
    dft.h_grid_weights.swap(staged_weights);
#endif
    dft.max_grid_size = static_cast<int>(grid.size());
}
