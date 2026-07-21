#pragma once

#include <cerrno>

#include "../../../Lennard_Jones_force/pair_activity.h"
#include "md_core_parse.hpp"

namespace Xponge
{

static float Native_LJ_Soft_Core_Read_Lambda(CONTROLLER* controller)
{
    const char* error_by = "Xponge::Native_Load_LJ_Soft_Core";
    if (!controller->Command_Exist("lambda_lj"))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand, error_by,
            "Reason:\n\t'lambda_lj' is required when charge endpoint files "
            "are provided\n");
    }
    const std::string token = controller->Command("lambda_lj");
    if (!Native_Core_Is_Strict_Decimal(token))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand, error_by,
            "Reason:\n\t'lambda_lj' must be a strict finite decimal in "
            "[0,1], received '%s'\n",
            token.c_str());
    }
    char* end = NULL;
    errno = 0;
    const float lambda = strtof(token.c_str(), &end);
    if (end != token.c_str() + token.size() || errno == ERANGE ||
        !Float_Memory_Is_Finite(&lambda) || lambda < 0.0f || lambda > 1.0f ||
        !Float_Memory_Is_Zero_Or_Normal(&lambda))
    {
        controller->Throw_Formatted_SPONGE_Error(
            spongeErrorConflictingCommand, error_by,
            "Reason:\n\t'lambda_lj' must be representable as a finite zero "
            "or normal float in [0,1], received '%s'\n",
            token.c_str());
    }
    return lambda;
}

static std::vector<float> Native_Load_LJ_Soft_Core_Charge_Endpoint(
    CONTROLLER* controller, const char* command, const char* endpoint_name,
    const int expected_atom_numbers)
{
    Native_Core_Parser parser(controller->Original_Command(command), command,
                              "Xponge::Native_Load_LJ_Soft_Core", controller);
    const int atom_numbers =
        parser.Validate_Atom_Count(parser.Read_Int("atom count"), "atom count");
    if (atom_numbers != expected_atom_numbers)
    {
        parser.Fail(spongeErrorConflictingCommand,
                    std::string(command) + " atom count " +
                        std::to_string(atom_numbers) +
                        " differs from LJ_soft_core_in_file atom count " +
                        std::to_string(expected_atom_numbers));
    }
    std::vector<float> endpoint(static_cast<std::size_t>(atom_numbers));
    for (int atom = 0; atom < atom_numbers; atom++)
    {
        const std::string field = Native_Core_Entry_Field(endpoint_name, atom);
        endpoint[atom] = parser.Checked_Float(parser.Read_Double(field), field);
    }
    parser.Ensure_End();
    parser.Close();
    return endpoint;
}

static void Native_Load_LJ_Soft_Core_Charge_Endpoints(System* system,
                                                      CONTROLLER* controller)
{
    const bool has_endpoint_a = controller->Command_Exist("charge_A_in_file");
    const bool has_endpoint_b = controller->Command_Exist("charge_B_in_file");
    if (has_endpoint_a != has_endpoint_b)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMissingCommand, "Xponge::Native_Load_LJ_Soft_Core",
            "Reason:\n\t'charge_A_in_file' and 'charge_B_in_file' must be "
            "provided together\n");
        return;
    }
    if (!has_endpoint_a)
    {
        return;
    }

    LJSoftCore& soft_core = system->classical_force_field.lj_soft_core;
    if (soft_core.atom_numbers <= 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "Xponge::Native_Load_LJ_Soft_Core",
            "Reason:\n\tcharge endpoint files require "
            "LJ_soft_core_in_file\n");
        return;
    }
    std::vector<float> endpoint_a = Native_Load_LJ_Soft_Core_Charge_Endpoint(
        controller, "charge_A_in_file", "endpoint A charge",
        soft_core.atom_numbers);
    std::vector<float> endpoint_b = Native_Load_LJ_Soft_Core_Charge_Endpoint(
        controller, "charge_B_in_file", "endpoint B charge",
        soft_core.atom_numbers);
    const float lambda = Native_LJ_Soft_Core_Read_Lambda(controller);

    if (system->atoms.charge.size() !=
        static_cast<std::size_t>(soft_core.atom_numbers))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "Xponge::Native_Load_LJ_Soft_Core",
            "Reason:\n\tthe current charge array and soft-core charge "
            "endpoints have different atom counts\n");
        return;
    }
    const bool current_was_supplied =
        controller->Command_Exist("charge_in_file");
    std::vector<float> current_charge(
        static_cast<std::size_t>(soft_core.atom_numbers));
    for (int atom = 0; atom < soft_core.atom_numbers; atom++)
    {
        const float charge_a = endpoint_a[atom];
        const float charge_b = endpoint_b[atom];
        const float derivative = charge_b - charge_a;
        const float current =
            PairwiseInteraction::Interpolate_Charge(charge_a, charge_b, lambda);
        if (!Float_Memory_Is_Finite(&derivative) ||
            !Float_Memory_Is_Zero_Or_Normal(&derivative) ||
            !Float_Memory_Is_Finite(&current) ||
            !Float_Memory_Is_Zero_Or_Normal(&current))
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorBadFileFormat, "Xponge::Native_Load_LJ_Soft_Core",
                "Reason:\n\tcharge endpoints for atom %d produce a "
                "non-finite or subnormal current/difference at lambda=%g "
                "(qA=%g, qB=%g)\n",
                atom, lambda, charge_a, charge_b);
            return;
        }
        if (current_was_supplied && system->atoms.charge[atom] != current)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorConflictingCommand,
                "Xponge::Native_Load_LJ_Soft_Core",
                "Reason:\n\tcharge_in_file entry %d (%g) does not equal "
                "fma(lambda_lj, qB-qA, qA) (%g) for qA=%g, qB=%g, "
                "lambda=%g\n",
                atom, system->atoms.charge[atom], current, charge_a, charge_b,
                lambda);
            return;
        }
        current_charge[atom] = current;
    }
    soft_core.charge_A = std::move(endpoint_a);
    soft_core.charge_B = std::move(endpoint_b);
    system->atoms.charge = std::move(current_charge);
}

static void Native_Load_LJ_Soft_Core(LJSoftCore* lj_soft,
                                     CONTROLLER* controller,
                                     const char* module_name = "LJ_soft_core",
                                     int atom_numbers_hint = 0)
{
    if (!controller->Command_Exist(module_name, "in_file"))
    {
        return;
    }

    const char* input_path =
        controller->Original_Command(module_name, "in_file");
    Native_Core_Parser parser(input_path, "LJ_soft_core_in_file",
                              "Xponge::Native_Load_LJ_Soft_Core", controller);
    LJSoftCore parsed;
    parsed.atom_numbers =
        parser.Validate_Atom_Count(parser.Read_Int("atom count"), "atom count");
    parsed.atom_type_numbers_A = parser.Read_Int("endpoint A atom type count");
    parsed.atom_type_numbers_B = parser.Read_Int("endpoint B atom type count");
    const int pair_type_numbers_A = parser.Validate_Triangular_Type_Count(
        parsed.atom_type_numbers_A, "endpoint A atom type count");
    const int pair_type_numbers_B = parser.Validate_Triangular_Type_Count(
        parsed.atom_type_numbers_B, "endpoint B atom type count");
    if (parsed.atom_numbers > 0 &&
        (parsed.atom_type_numbers_A == 0 || parsed.atom_type_numbers_B == 0))
    {
        parser.Fail(spongeErrorBadFileFormat,
                    "LJ_soft_core_in_file has atoms but an endpoint has no "
                    "LJ atom types");
    }
    if (atom_numbers_hint > 0 && parsed.atom_numbers != atom_numbers_hint)
    {
        parser.Fail(spongeErrorConflictingCommand,
                    "LJ_soft_core_in_file atom count " +
                        std::to_string(parsed.atom_numbers) +
                        " differs from the previously loaded atom count " +
                        std::to_string(atom_numbers_hint));
    }
    for (int i = 0; i < pair_type_numbers_A; i++)
    {
        const std::string field = Native_Core_Entry_Field("endpoint A LJ A", i);
        parsed.LJ_AA.push_back(parser.Checked_Float(
            parser.Read_Double(field) * 12.0, field + " scaled by 12"));
    }
    for (int i = 0; i < pair_type_numbers_A; i++)
    {
        const std::string field = Native_Core_Entry_Field("endpoint A LJ B", i);
        parsed.LJ_AB.push_back(parser.Checked_Float(
            parser.Read_Double(field) * 6.0, field + " scaled by 6"));
    }
    for (int i = 0; i < pair_type_numbers_B; i++)
    {
        const std::string field = Native_Core_Entry_Field("endpoint B LJ A", i);
        parsed.LJ_BA.push_back(parser.Checked_Float(
            parser.Read_Double(field) * 12.0, field + " scaled by 12"));
    }
    for (int i = 0; i < pair_type_numbers_B; i++)
    {
        const std::string field = Native_Core_Entry_Field("endpoint B LJ B", i);
        parsed.LJ_BB.push_back(parser.Checked_Float(
            parser.Read_Double(field) * 6.0, field + " scaled by 6"));
    }
    for (int i = 0; i < parsed.atom_numbers; i++)
    {
        const std::string field_A =
            Native_Core_Entry_Field("endpoint A atom LJ type", i);
        const std::string field_B =
            Native_Core_Entry_Field("endpoint B atom LJ type", i);
        const int atom_type_A = parser.Read_Int(field_A);
        const int atom_type_B = parser.Read_Int(field_B);
        if (atom_type_A < 0 || atom_type_A >= parsed.atom_type_numbers_A)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "LJ_soft_core_in_file " + field_A + " " +
                            std::to_string(atom_type_A) + " is outside [0, " +
                            std::to_string(parsed.atom_type_numbers_A) + ")");
        }
        if (atom_type_B < 0 || atom_type_B >= parsed.atom_type_numbers_B)
        {
            parser.Fail(spongeErrorBadFileFormat,
                        "LJ_soft_core_in_file " + field_B + " " +
                            std::to_string(atom_type_B) + " is outside [0, " +
                            std::to_string(parsed.atom_type_numbers_B) + ")");
        }
        parsed.atom_LJ_type_A.push_back(atom_type_A);
        parsed.atom_LJ_type_B.push_back(atom_type_B);
    }
    parser.Ensure_End();
    parser.Close();

    if (controller->Command_Exist("subsys_division_in_file"))
    {
        const char* division_path =
            controller->Original_Command("subsys_division_in_file");
        Native_Core_Parser division_parser(
            division_path, "subsys_division_in_file",
            "Xponge::Native_Load_LJ_Soft_Core", controller);
        const int atom_numbers = division_parser.Validate_Atom_Count(
            division_parser.Read_Int("atom count"), "atom count");
        if (atom_numbers != parsed.atom_numbers)
        {
            division_parser.Fail(
                spongeErrorConflictingCommand,
                "subsys_division_in_file atom count " +
                    std::to_string(atom_numbers) +
                    " differs from LJ_soft_core_in_file atom count " +
                    std::to_string(parsed.atom_numbers));
        }
        for (int i = 0; i < atom_numbers; i++)
        {
            parsed.subsystem_division.push_back(division_parser.Read_Int(
                Native_Core_Entry_Field("subsystem mask", i)));
        }
        division_parser.Ensure_End();
        division_parser.Close();
    }
    *lj_soft = std::move(parsed);
}

static void Native_Load_LJ_Soft_Core(System* system, CONTROLLER* controller)
{
    Native_Load_LJ_Soft_Core(&system->classical_force_field.lj_soft_core,
                             controller, "LJ_soft_core",
                             Load_Get_Atom_Numbers(system));
    Native_Load_LJ_Soft_Core_Charge_Endpoints(system, controller);
}

}  // namespace Xponge
