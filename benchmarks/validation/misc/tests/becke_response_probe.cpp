#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>

#include "quantum_chemistry/dft/becke.hpp"

struct Probe_Vector
{
    double x;
    double y;
    double z;
};

static Probe_Vector Subtract(const Probe_Vector& a, const Probe_Vector& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

static double Norm(const Probe_Vector& value)
{
    return std::sqrt(value.x * value.x + value.y * value.y +
                     value.z * value.z);
}

static std::uint64_t Double_Bits(double value)
{
    std::uint64_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static double Double_From_Bits(std::uint64_t bits)
{
    double value = 0.0;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
}

template <std::size_t N>
static std::array<double, N> Partition_Weights(
    const std::array<Probe_Vector, N>& atoms,
    const std::array<double, N>& radii, const Probe_Vector& point)
{
    std::array<double, N> products{};
    double denominator = 0.0;
    for (std::size_t a = 0; a < N; ++a)
    {
        products[a] = 1.0;
        const Probe_Vector point_a = Subtract(point, atoms[a]);
        const double ra = Norm(point_a);
        for (std::size_t b = 0; b < N; ++b)
        {
            if (a == b) continue;
            const Probe_Vector point_b = Subtract(point, atoms[b]);
            const Probe_Vector atom_ab = Subtract(atoms[a], atoms[b]);
            const double rb = Norm(point_b);
            const double rab = Norm(atom_ab);
            const double mu = (ra - rb) / rab;
            products[a] *= QC_Becke_Shape(
                QC_Becke_Size_Adjusted_Mu(mu, radii[a], radii[b]));
        }
        denominator += products[a];
    }
    for (double& product : products) product /= denominator;
    return products;
}

template <std::size_t N>
static std::array<Probe_Vector, N> Analytic_Weight_Response(
    const std::array<Probe_Vector, N>& atoms,
    const std::array<double, N>& radii, const Probe_Vector& point,
    std::size_t grid_owner, std::size_t target)
{
    const std::array<double, N> weights =
        Partition_Weights(atoms, radii, point);
    std::array<Probe_Vector, N> response{};
    for (std::size_t a = 0; a < N; ++a)
    {
        const Probe_Vector point_a = Subtract(point, atoms[a]);
        const double ra = Norm(point_a);
        for (std::size_t b = 0; b < N; ++b)
        {
            if (a == b) continue;
            const Probe_Vector point_b = Subtract(point, atoms[b]);
            const Probe_Vector atom_ab = Subtract(atoms[a], atoms[b]);
            const double rb = Norm(point_b);
            const double rab = Norm(atom_ab);
            const double inverse_rab = 1.0 / rab;
            const double mu = (ra - rb) * inverse_rab;
            const double adjustment =
                QC_Becke_Size_Adjustment_Coefficient(radii[a], radii[b]);
            const double adjusted_mu =
                mu + adjustment * (1.0 - mu * mu);
            double shape = 0.0, shape_derivative = 0.0;
            QC_Becke_Shape_And_Derivative(adjusted_mu, shape,
                                          shape_derivative);
            const double coefficient =
                (a == target ? 1.0 : 0.0) - weights[a];
            const double scale = coefficient * (shape_derivative / shape) *
                                 (1.0 - 2.0 * adjustment * mu);

            const Probe_Vector ua = {point_a.x / ra, point_a.y / ra,
                                     point_a.z / ra};
            const Probe_Vector ub = {point_b.x / rb, point_b.y / rb,
                                     point_b.z / rb};
            const Probe_Vector h = {atom_ab.x * inverse_rab,
                                    atom_ab.y * inverse_rab,
                                    atom_ab.z * inverse_rab};
            double pair_a[3], pair_b[3], pair_grid[3];
            QC_Becke_Pair_Response(scale, mu, inverse_rab, ua.x, ua.y, ua.z,
                                   ub.x, ub.y, ub.z, h.x, h.y, h.z, pair_a,
                                   pair_b, pair_grid);
            double* response_a = &response[a].x;
            double* response_b = &response[b].x;
            double* response_owner = &response[grid_owner].x;
            for (int axis = 0; axis < 3; ++axis)
            {
                if (grid_owner == a)
                {
                    response_a[axis] -= pair_b[axis];
                    response_b[axis] += pair_b[axis];
                }
                else if (grid_owner == b)
                {
                    response_a[axis] += pair_a[axis];
                    response_b[axis] -= pair_a[axis];
                }
                else
                {
                    response_a[axis] += pair_a[axis];
                    response_b[axis] += pair_b[axis];
                    response_owner[axis] += pair_grid[axis];
                }
            }
        }
    }
    for (Probe_Vector& atom_response : response)
    {
        atom_response.x *= weights[target];
        atom_response.y *= weights[target];
        atom_response.z *= weights[target];
    }
    return response;
}

int main()
{
    const std::uint64_t mu_bits[] = {
        UINT64_C(0xbff0000000000000), UINT64_C(0xbfefffffffffffff),
        UINT64_C(0xbfe0000000000000), UINT64_C(0xbfb999999999999a),
        UINT64_C(0x0000000000000000), UINT64_C(0x3fb999999999999a),
        UINT64_C(0x3fe0000000000000), UINT64_C(0x3fefffffffffffff),
        UINT64_C(0x3ff0000000000000)};
    for (std::uint64_t bits : mu_bits)
    {
        const double mu = Double_From_Bits(bits);
        double shape = 0.0, derivative = 0.0;
        QC_Becke_Shape_And_Derivative(mu, shape, derivative);
        std::printf("shape %016llx %.17e %.17e %.17e\n",
                    static_cast<unsigned long long>(bits), QC_Becke_Shape(mu),
                    shape, derivative);
    }

    const double maximum = std::numeric_limits<double>::max();
    const double minimum = std::numeric_limits<double>::min();
    const std::array<std::array<double, 2>, 7> radius_pairs = {{
        {{1.0, 1.0}},
        {{2.0, 1.0}},
        {{1.0, 2.0}},
        {{2.4, 1.0}},
        {{1.0, 2.4}},
        {{maximum, minimum}},
        {{minimum, maximum}},
    }};
    for (const auto& pair : radius_pairs)
        std::printf("radius %.17e %.17e %.17e\n", pair[0], pair[1],
                    QC_Becke_Size_Adjustment_Coefficient(pair[0], pair[1]));

    const std::array<Probe_Vector, 3> atoms = {{
        {-0.7, 0.2, 0.1},
        {0.4, -0.3, 0.5},
        {0.9, 0.8, -0.6},
    }};
    const std::array<double, 3> radii = {{0.8, 1.1, 0.65}};
    constexpr std::size_t owner = 1;
    const Probe_Vector local_point = {0.31, -0.22, 0.17};
    const Probe_Vector point = {atoms[owner].x + local_point.x,
                                atoms[owner].y + local_point.y,
                                atoms[owner].z + local_point.z};
    constexpr double step = 2.0e-6;

    for (std::size_t target : {std::size_t(0), owner})
    {
        const auto analytic =
            Analytic_Weight_Response(atoms, radii, point, owner, target);
        double translation[3] = {0.0, 0.0, 0.0};
        for (std::size_t nucleus = 0; nucleus < atoms.size(); ++nucleus)
        {
            const double* derivative = &analytic[nucleus].x;
            for (int axis = 0; axis < 3; ++axis)
            {
                auto plus_atoms = atoms;
                auto minus_atoms = atoms;
                Probe_Vector plus_point = point;
                Probe_Vector minus_point = point;
                (&plus_atoms[nucleus].x)[axis] += step;
                (&minus_atoms[nucleus].x)[axis] -= step;
                if (nucleus == owner)
                {
                    (&plus_point.x)[axis] += step;
                    (&minus_point.x)[axis] -= step;
                }
                const double plus =
                    Partition_Weights(plus_atoms, radii, plus_point)[target];
                const double minus =
                    Partition_Weights(minus_atoms, radii, minus_point)[target];
                std::printf(
                    "weight-derivative %zu %zu %d %.17e %.17e %.17e %.17e\n",
                    target, nucleus, axis, derivative[axis], plus, minus,
                    step);
                translation[axis] += derivative[axis];
            }
        }
        std::printf("translation %zu %.17e %.17e %.17e\n", target,
                    translation[0], translation[1], translation[2]);
    }

    const Probe_Vector point_a = Subtract(point, atoms[0]);
    const Probe_Vector point_b = Subtract(point, atoms[2]);
    const Probe_Vector atom_ab = Subtract(atoms[0], atoms[2]);
    const double ra = Norm(point_a), rb = Norm(point_b), rab = Norm(atom_ab);
    const double mu = (ra - rb) / rab;
    double pair_a[3], pair_b[3], pair_grid[3];
    QC_Becke_Pair_Response(
        1.234567890123, mu, 1.0 / rab, point_a.x / ra, point_a.y / ra,
        point_a.z / ra, point_b.x / rb, point_b.y / rb, point_b.z / rb,
        atom_ab.x / rab, atom_ab.y / rab, atom_ab.z / rab, pair_a, pair_b,
        pair_grid);
    for (int axis = 0; axis < 3; ++axis)
        std::printf("pair %d %.17e %.17e %.17e\n", axis, pair_a[axis],
                    pair_b[axis], pair_grid[axis]);
    return 0;
}
