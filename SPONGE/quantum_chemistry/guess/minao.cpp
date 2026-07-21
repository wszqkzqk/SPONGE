#include "minao.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "../quantum_chemistry.h"
#include "minao_occupancy.hpp"

// 设备核函数：将 AO 占据数写入密度矩阵对角线
static __global__ void QC_Set_Diagonal_Kernel(const int count, const int nao,
                                              const int* ao_indices,
                                              const float* ao_occs, float* P)
{
    SIMPLE_DEVICE_FOR(i, count)
    {
        P[ao_indices[i] * nao + ao_indices[i]] = ao_occs[i];
    }
}

void QC_Build_Minao_Guess(const QC_MOLECULE& mol,
                          const QC_SCF_Runtime_State& runtime, float* d_P,
                          float* d_P_beta)
{
    const int nao = mol.nao;
    const bool unrestricted = runtime.unrestricted;
    const bool is_spherical = mol.is_spherical;
    if (mol.natm <= 0 || mol.nbas < 0 || nao <= 0 || mol.nao2 <= 0 ||
        d_P == nullptr || (unrestricted && d_P_beta == nullptr))
        throw std::domain_error("invalid MINAO molecule or density storage");
    if (nao > std::numeric_limits<int>::max() / nao || mol.nao2 != nao * nao)
        throw std::logic_error("MINAO AO square dimension is inconsistent");
    const std::size_t atom_count = static_cast<std::size_t>(mol.natm);
    const std::size_t shell_count = static_cast<std::size_t>(mol.nbas);
    if (shell_count > std::numeric_limits<std::size_t>::max() / 8U ||
        mol.h_bas.size() < 8U * shell_count ||
        mol.h_l_list.size() < shell_count ||
        mol.h_ao_offsets.size() < shell_count ||
        mol.h_ao_offsets_sph.size() < shell_count ||
        mol.h_atomic_numbers.size() != atom_count ||
        mol.h_Z.size() != atom_count || mol.h_ecp_n_core.size() != atom_count)
        throw std::logic_error("MINAO molecule metadata is incomplete");

    auto shell_ao_dimension = [&](int angular_momentum)
    {
        if (angular_momentum < 0)
            throw std::domain_error(
                "MINAO shell angular momentum cannot be negative");
        const long long dimension =
            is_spherical
                ? 2LL * angular_momentum + 1LL
                : (static_cast<long long>(angular_momentum) + 1LL) *
                      (static_cast<long long>(angular_momentum) + 2LL) / 2LL;
        if (dimension <= 0 || dimension > std::numeric_limits<int>::max())
            throw std::overflow_error("MINAO shell AO dimension overflows int");
        return static_cast<int>(dimension);
    };

    // 按原子分组壳层
    std::vector<std::vector<int>> atom_shells(mol.natm);
    for (int ish = 0; ish < mol.nbas; ish++)
    {
        const int iatm = mol.h_bas[(std::size_t)ish * 8U];
        if (iatm < 0 || iatm >= mol.natm)
            throw std::logic_error("MINAO shell atom index is out of range");
        atom_shells[iatm].push_back(ish);
    }

    // Build the neutral-atom/ECP-core-removed preference first. Molecular
    // charge and multiplicity are global constraints and must not be imposed
    // independently atom by atom.
    std::vector<double> preferred_total(static_cast<std::size_t>(nao), 0.0);
    std::vector<unsigned char> ao_seen(static_cast<std::size_t>(nao), 0U);

    for (int iatm = 0; iatm < mol.natm; iatm++)
    {
        const int atomic_number = mol.Atomic_Number(iatm);
        const std::array<int, 4> electrons_by_l =
            sponge_qc_minao::Explicit_Electrons_By_Angular_Momentum(
                atomic_number, mol.h_ecp_n_core[iatm]);
        int explicit_electron_count = 0;
        for (int count : electrons_by_l) explicit_electron_count += count;
        if (explicit_electron_count != mol.Effective_Nuclear_Charge(iatm))
            throw std::logic_error(
                "MINAO explicit-electron count differs from the effective "
                "nuclear charge");

        // 按角动量统计该原子在当前基组中的 AO 总数
        const int l_max = 3;
        int nao_per_l[4] = {0, 0, 0, 0};
        for (int ish : atom_shells[iatm])
        {
            const int l = mol.h_l_list[ish];
            const int dimension = shell_ao_dimension(l);
            if (l > l_max) continue;
            if (nao_per_l[l] > std::numeric_limits<int>::max() - dimension)
                throw std::overflow_error(
                    "MINAO per-angular-momentum AO count overflows int");
            nao_per_l[l] += dimension;
        }
        for (int l = 0; l <= l_max; ++l)
        {
            if (electrons_by_l[l] != 0 && nao_per_l[l] == 0)
                throw std::domain_error(
                    "MINAO cannot represent " +
                    std::to_string(electrons_by_l[l]) +
                    " explicit electron(s) with angular momentum l=" +
                    std::to_string(l) + " for atomic number " +
                    std::to_string(atomic_number));
            if (static_cast<long long>(electrons_by_l[l]) > 2LL * nao_per_l[l])
                throw std::domain_error(
                    "MINAO cannot place " + std::to_string(electrons_by_l[l]) +
                    " explicit electron(s) into the available l=" +
                    std::to_string(l) + " AO capacity for atomic number " +
                    std::to_string(atomic_number));
        }

        // Spread each neutral explicit-atom angular population over the AOs
        // of that l. Higher-l polarization AOs carry zero preference but are
        // retained as legitimate capacity for an anionic molecular target.
        for (int ish : atom_shells[iatm])
        {
            const int l = mol.h_l_list[ish];
            const int nao_shell = shell_ao_dimension(l);
            const int ao_start = is_spherical ? mol.h_ao_offsets_sph[ish]
                                              : mol.h_ao_offsets[ish];
            if (ao_start < 0 || ao_start > nao || nao_shell > nao - ao_start)
                throw std::logic_error("MINAO shell AO range is out of bounds");
            const double total_occ =
                l <= l_max ? static_cast<double>(electrons_by_l[l]) /
                                 static_cast<double>(nao_per_l[l])
                           : 0.0;

            for (int k = 0; k < nao_shell; k++)
            {
                const int ao = ao_start + k;
                if (ao < 0 || ao >= nao || ao_seen[(std::size_t)ao] != 0U)
                    throw std::logic_error(
                        "MINAO AO layout is overlapping or out of range");
                ao_seen[(std::size_t)ao] = 1U;
                preferred_total[(std::size_t)ao] = total_occ;
            }
        }
    }

    for (unsigned char seen : ao_seen)
        if (seen == 0U)
            throw std::logic_error("MINAO AO layout does not cover every AO");

    const long long runtime_electron_target =
        unrestricted ? static_cast<long long>(runtime.n_alpha) + runtime.n_beta
                     : 2LL * runtime.n_alpha;
    if (runtime_electron_target != mol.nelectron)
        throw std::logic_error(
            "MINAO runtime occupations do not match the molecular electron "
            "count");

    const sponge_qc_minao::Global_AO_Occupancies occupations =
        sponge_qc_minao::Allocate_Global_AO_Occupancies(
            preferred_total, unrestricted, runtime.n_alpha, runtime.n_beta,
            runtime.occ_factor);
    if (occupations.alpha.size() != static_cast<std::size_t>(nao) ||
        (unrestricted &&
         occupations.beta.size() != static_cast<std::size_t>(nao)))
        throw std::logic_error("MINAO global occupation size mismatch");

    std::vector<int> h_ao_indices(static_cast<std::size_t>(nao));
    std::vector<float> h_occ_alpha(static_cast<std::size_t>(nao), 0.0f);
    std::vector<float> h_occ_beta(
        unrestricted ? static_cast<std::size_t>(nao) : 0U, 0.0f);
    double staged_alpha_sum = 0.0;
    double staged_beta_sum = 0.0;
    for (int ao = 0; ao < nao; ++ao)
    {
        h_ao_indices[(std::size_t)ao] = ao;
        h_occ_alpha[(std::size_t)ao] =
            static_cast<float>(occupations.alpha[(std::size_t)ao]);
        staged_alpha_sum += h_occ_alpha[(std::size_t)ao];
        if (unrestricted)
        {
            h_occ_beta[(std::size_t)ao] =
                static_cast<float>(occupations.beta[(std::size_t)ao]);
            staged_beta_sum += h_occ_beta[(std::size_t)ao];
        }
    }
    const double target_alpha =
        unrestricted ? runtime.n_alpha : runtime.occ_factor * runtime.n_alpha;
    const double target_beta = unrestricted ? runtime.n_beta : 0.0;
    const double float_rounding_bound =
        8.0 * std::numeric_limits<float>::epsilon() *
        std::max({1.0, static_cast<double>(nao), target_alpha, target_beta});
    if (std::fabs(staged_alpha_sum - target_alpha) > float_rounding_bound ||
        std::fabs(staged_beta_sum - target_beta) > float_rounding_bound)
        throw std::logic_error(
            "MINAO float occupations do not preserve molecular electrons");

    // 用设备核函数写入密度矩阵对角线
    const int threads = 256;
    deviceMemset(d_P, 0, sizeof(float) * (std::size_t)mol.nao2);
    if (unrestricted && d_P_beta != nullptr)
        deviceMemset(d_P_beta, 0, sizeof(float) * (std::size_t)mol.nao2);

    if (!h_ao_indices.empty())
    {
        const int n = (int)h_ao_indices.size();
        int* d_idx = nullptr;
        float* d_occ = nullptr;
        Device_Malloc_Safely((void**)&d_idx, sizeof(int) * n);
        Device_Malloc_Safely((void**)&d_occ, sizeof(float) * n);
        deviceMemcpy(d_idx, h_ao_indices.data(), sizeof(int) * n,
                     deviceMemcpyHostToDevice);
        deviceMemcpy(d_occ, h_occ_alpha.data(), sizeof(float) * n,
                     deviceMemcpyHostToDevice);
        Launch_Device_Kernel(QC_Set_Diagonal_Kernel,
                             (n + threads - 1) / threads, threads, 0, 0, n, nao,
                             d_idx, d_occ, d_P);
        deviceFree(d_idx);
        deviceFree(d_occ);
    }

    if (unrestricted && d_P_beta && !h_ao_indices.empty())
    {
        const int n = (int)h_ao_indices.size();
        int* d_idx = nullptr;
        float* d_occ = nullptr;
        Device_Malloc_Safely((void**)&d_idx, sizeof(int) * n);
        Device_Malloc_Safely((void**)&d_occ, sizeof(float) * n);
        deviceMemcpy(d_idx, h_ao_indices.data(), sizeof(int) * n,
                     deviceMemcpyHostToDevice);
        deviceMemcpy(d_occ, h_occ_beta.data(), sizeof(float) * n,
                     deviceMemcpyHostToDevice);
        Launch_Device_Kernel(QC_Set_Diagonal_Kernel,
                             (n + threads - 1) / threads, threads, 0, 0, n, nao,
                             d_idx, d_occ, d_P_beta);
        deviceFree(d_idx);
        deviceFree(d_occ);
    }
}
