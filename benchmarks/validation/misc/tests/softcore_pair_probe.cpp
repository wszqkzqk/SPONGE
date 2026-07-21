#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

#include "Lennard_Jones_force/LJ_soft_core.h"

namespace
{

struct Scenario
{
    const char* name;
    float AA;
    float AB;
    float BA;
    float BB;
    float q1_a;
    float q1_b;
    float q2_a;
    float q2_b;
};

VECTOR_LJ_SOFT_TYPE Make_Atom(float q_a, float q_b, float lambda)
{
    VECTOR_LJ_SOFT_TYPE atom{};
    atom.charge_A = q_a;
    atom.charge_B = q_b;
    atom.charge_BA = q_b - q_a;
    atom.charge = PairwiseInteraction::Interpolate_Charge(q_a, q_b, lambda);
    return atom;
}

LJ_SOFT_CORE_PAIR_RESULT Evaluate(const Scenario& scenario, float lambda,
                                  float distance, bool force, bool energy,
                                  bool derivative)
{
    const VECTOR_LJ_SOFT_TYPE r1 =
        Make_Atom(scenario.q1_a, scenario.q1_b, lambda);
    const VECTOR_LJ_SOFT_TYPE r2 =
        Make_Atom(scenario.q2_a, scenario.q2_b, lambda);
    return Evaluate_LJ_Soft_Core_Pair(
        r1, r2, distance, scenario.AA, scenario.AB, scenario.BA, scenario.BB,
        0.27f, lambda, 0.5f, 1.0f, 729.0f, 0.0f, force, energy, true,
        derivative);
}

double Energy(const Scenario& scenario, float lambda, float distance = 2.0f)
{
    const auto result =
        Evaluate(scenario, lambda, distance, false, true, false);
    return static_cast<double>(result.lj_energy) + result.coulomb_energy;
}

bool Close(double actual, double expected, double relative, double absolute)
{
    return std::fabs(actual - expected) <=
           absolute + relative * std::fmax(std::fabs(actual),
                                           std::fabs(expected));
}

bool Check_Derivative(const Scenario& scenario, float lambda)
{
    const float h = 2.0e-4f;
    const double numeric =
        (Energy(scenario, lambda + h) - Energy(scenario, lambda - h)) /
        (2.0 * h);
    const auto result =
        Evaluate(scenario, lambda, 2.0f, false, false, true);
    const double analytic = static_cast<double>(result.du_dlambda_lj) +
                            result.du_dlambda_coulomb;
    if (!std::isfinite(analytic) || !std::isfinite(numeric) ||
        !Close(analytic, numeric, 8.0e-3, 3.0e-4))
    {
        std::fprintf(stderr,
                     "%s: analytic dH/dlambda %.9g != finite difference %.9g\n",
                     scenario.name, analytic, numeric);
        return false;
    }
    return true;
}

bool Check_Force(const Scenario& scenario, float lambda)
{
    const float distance = 2.0f;
    const float h = 2.0e-4f;
    const double numeric =
        (Energy(scenario, lambda, distance + h) -
         Energy(scenario, lambda, distance - h)) /
        (2.0 * h);
    const auto result =
        Evaluate(scenario, lambda, distance, true, false, false);
    // Production multiplies the returned radial coefficient by the pair
    // displacement.  For a displacement along +x its magnitude is F_scalar*r
    // and must equal dU/dr.
    const double analytic = static_cast<double>(result.force) * distance;
    if (!std::isfinite(analytic) || !std::isfinite(numeric) ||
        !Close(analytic, numeric, 8.0e-3, 3.0e-4))
    {
        std::fprintf(stderr,
                     "%s: analytic radial force %.9g != finite difference "
                     "%.9g\n",
                     scenario.name, analytic, numeric);
        return false;
    }
    return true;
}

double Hard_Endpoint_Energy(const Scenario& scenario, bool endpoint_b)
{
    const float lj_a = endpoint_b ? scenario.BA : scenario.AA;
    const float lj_b = endpoint_b ? scenario.BB : scenario.AB;
    const float q1 = endpoint_b ? scenario.q1_b : scenario.q1_a;
    const float q2 = endpoint_b ? scenario.q2_b : scenario.q2_a;
    const VECTOR_LJ_SOFT_TYPE r1 = Make_Atom(q1, q1, 0.0f);
    const VECTOR_LJ_SOFT_TYPE r2 = Make_Atom(q2, q2, 0.0f);
    double energy = 0.0;
    if (PairwiseInteraction::Lennard_Jones_Is_Active(lj_a, lj_b))
    {
        energy += Get_LJ_Energy(r1, r2, 2.0f, lj_a, lj_b);
    }
    if (PairwiseInteraction::Coulomb_Is_Active(q1, q2))
    {
        energy += Get_Direct_Coulomb_Energy(r1, r2, 2.0f, 0.27f);
    }
    return energy;
}

bool Check_Lambda_Endpoint(const Scenario& scenario, bool endpoint_b)
{
    const float lambda = endpoint_b ? 1.0f : 0.0f;
    const auto result =
        Evaluate(scenario, lambda, 2.0f, false, true, true);
    const double energy =
        static_cast<double>(result.lj_energy) + result.coulomb_energy;
    const double expected_energy =
        Hard_Endpoint_Energy(scenario, endpoint_b);
    if (!std::isfinite(energy) ||
        !Close(energy, expected_energy, 2.0e-6, 2.0e-6))
    {
        std::fprintf(stderr,
                     "%s lambda=%g: endpoint energy %.9g != hard endpoint "
                     "%.9g\n",
                     scenario.name, lambda, energy, expected_energy);
        return false;
    }

    const float h = 5.0e-4f;
    const double numeric = endpoint_b
                               ? (3.0 * Energy(scenario, 1.0f) -
                                  4.0 * Energy(scenario, 1.0f - h) +
                                  Energy(scenario, 1.0f - 2.0f * h)) /
                                     (2.0 * h)
                               : (-3.0 * Energy(scenario, 0.0f) +
                                  4.0 * Energy(scenario, h) -
                                  Energy(scenario, 2.0f * h)) /
                                     (2.0 * h);
    const double analytic = static_cast<double>(result.du_dlambda_lj) +
                            result.du_dlambda_coulomb;
    if (!std::isfinite(analytic) || !std::isfinite(numeric) ||
        !Close(analytic, numeric, 2.0e-2, 1.0e-3))
    {
        std::fprintf(stderr,
                     "%s lambda=%g: endpoint dH/dlambda %.9g != one-sided "
                     "finite difference %.9g\n",
                     scenario.name, lambda, analytic, numeric);
        return false;
    }
    return true;
}

bool Check_Endpoint_Contract()
{
    const float lambda = 0.25f;
    const float charge_a = 1.0f;
    const float charge_b = -0.5f;
    const float current = PairwiseInteraction::Interpolate_Charge(
        charge_a, charge_b, lambda);
    float derivative = charge_b - charge_a;
    auto validation = PairwiseInteraction::Validate_Charge_Endpoints(
        current, charge_a, charge_b, &derivative, lambda);
    if (validation.error !=
            PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONE ||
        validation.current != current || validation.derivative != derivative)
    {
        std::fprintf(stderr, "valid explicit charge endpoints were rejected\n");
        return false;
    }

    validation = PairwiseInteraction::Validate_Charge_Endpoints(
        current + 0.125f, charge_a, charge_b, &derivative, lambda);
    if (validation.error !=
        PairwiseInteraction::CHARGE_ENDPOINT_ERROR_CURRENT_MISMATCH)
    {
        std::fprintf(stderr, "mismatched current charge was accepted\n");
        return false;
    }

    float wrong_derivative = derivative + 0.125f;
    validation = PairwiseInteraction::Validate_Charge_Endpoints(
        current, charge_a, charge_b, &wrong_derivative, lambda);
    if (validation.error !=
        PairwiseInteraction::CHARGE_ENDPOINT_ERROR_DERIVATIVE_MISMATCH)
    {
        std::fprintf(stderr, "mismatched charge derivative was accepted\n");
        return false;
    }

    const float nan = std::numeric_limits<float>::quiet_NaN();
    validation = PairwiseInteraction::Validate_Charge_Endpoints(
        current, nan, charge_b, NULL, lambda);
    if (validation.error !=
        PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONFINITE)
    {
        std::fprintf(stderr,
                     "non-finite charge endpoint survived -ffast-math check\n");
        return false;
    }

    validation = PairwiseInteraction::Validate_Charge_Endpoints(
        current, std::numeric_limits<float>::max(),
        -std::numeric_limits<float>::max(), NULL, lambda);
    if (validation.error !=
        PairwiseInteraction::CHARGE_ENDPOINT_ERROR_NONFINITE)
    {
        std::fprintf(stderr, "overflowing endpoint derivative was accepted\n");
        return false;
    }
    return true;
}

}  // namespace

int main()
{
    if (Get_LJ_Type(65534, 65534) != 2147450879 ||
        Get_LJ_Type(65535, 0) != -1 || Get_LJ_Type(0, 65535) != -1)
    {
        std::fprintf(stderr,
                     "LJ triangular indexing overflow boundary is wrong\n");
        return EXIT_FAILURE;
    }
    if (!Check_Endpoint_Contract())
    {
        return EXIT_FAILURE;
    }
    const auto false_active_false =
        PairwiseInteraction::Classify_Coulomb_Endpoints(1.0f, 0.0f, 0.0f,
                                                        -1.0f);
    if (false_active_false.state_a || false_active_false.state_b ||
        !false_active_false.ever_active || !false_active_false.Changes())
    {
        std::fprintf(stderr,
                     "false-active-false endpoint classification is wrong\n");
        return EXIT_FAILURE;
    }
    if (!PairwiseInteraction::Coulomb_Is_Active(1.0e-30f, 1.0e-30f))
    {
        std::fprintf(stderr, "nonzero underflowing charges were classified off\n");
        return EXIT_FAILURE;
    }

    const Scenario scenarios[] = {
        {"LJ-only", 12000.0f, 600.0f, 0.0f, 0.0f, 0.0f, 0.0f,
         0.0f, 0.0f},
        {"charge-only", 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,
         -1.0f, 0.0f},
        {"LJ-and-charge", 12000.0f, 600.0f, 0.0f, 0.0f, 1.0f,
         0.0f, -1.0f, 0.0f},
        {"false-active-false", 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,
         0.0f, 0.0f, -1.0f},
    };
    for (const Scenario& scenario : scenarios)
    {
        if (!Check_Derivative(scenario, 0.37f) ||
            !Check_Force(scenario, 0.37f) ||
            !Check_Lambda_Endpoint(scenario, false) ||
            !Check_Lambda_Endpoint(scenario, true))
        {
            return EXIT_FAILURE;
        }
    }

    // p=1 has a well-defined, nonzero derivative of lambda^p at lambda=0;
    // this explicitly exercises the endpoint branch without an FLT_MIN offset.
    const float endpoint_term = Get_Soft_Core_dU_dlambda(
        2.0f, 729.0f, 3.0f, 0.5f, 1.0f, 0.0f);
    if (!std::isfinite(endpoint_term) || endpoint_term == 0.0f)
    {
        std::fprintf(stderr, "p=1 lambda=0 derivative boundary is wrong\n");
        return EXIT_FAILURE;
    }
    const float smooth_endpoint_term = Get_Soft_Core_dU_dlambda(
        2.0f, 729.0f, 3.0f, 0.5f, 2.0f, 0.0f);
    if (smooth_endpoint_term != 0.0f)
    {
        std::fprintf(stderr, "p>1 lambda=0 derivative must be exact zero\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
