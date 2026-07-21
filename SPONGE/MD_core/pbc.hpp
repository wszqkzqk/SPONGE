#pragma once

#include <cstdint>

namespace
{

struct PBC_Build_Context
{
    CONTROLLER* controller;
    int error_number;
    const char* error_by;
    const char* operation;
};

struct PBC_Raw_Cell
{
    double a11;
    double a21;
    double a22;
    double a31;
    double a32;
    double a33;
};

struct PBC_Checked_Cell
{
    LTMatrix3 cell;
    LTMatrix3 rcell;
};

static std::uint32_t PBC_Float_Bits(float value)
{
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value));
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static std::uint64_t PBC_Double_Bits(double value)
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value));
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static bool PBC_Double_Is_Zero(double value)
{
    return (PBC_Double_Bits(value) & UINT64_C(0x7fffffffffffffff)) == 0;
}

static bool PBC_Double_Is_Strictly_Positive(double value)
{
    const std::uint64_t bits = PBC_Double_Bits(value);
    return (bits & UINT64_C(0x8000000000000000)) == 0 &&
           (bits & UINT64_C(0x7fffffffffffffff)) != 0;
}

static bool PBC_Fail(const PBC_Build_Context& context,
                     const std::string& reason)
{
    if (context.controller != NULL)
    {
        const std::string message = "Reason:\n\t" + reason +
                                    "\n\tCell operation: " + context.operation +
                                    "\n";
        context.controller->Throw_SPONGE_Error(
            context.error_number, context.error_by, message.c_str());
    }
    return false;
}

static bool PBC_Read_Positive_Box_Length(float value, const char* axis,
                                         const PBC_Build_Context& context,
                                         double* result)
{
    const std::string field = std::string("box length ") + axis;
    if (!Float_Memory_Is_Finite(&value))
    {
        return PBC_Fail(context, field + " is non-finite");
    }
    const std::uint32_t bits = PBC_Float_Bits(value);
    const std::uint32_t magnitude = bits & UINT32_C(0x7fffffff);
    if ((bits & UINT32_C(0x80000000)) != 0 || magnitude == 0)
    {
        return PBC_Fail(context,
                        field + " must be finite and strictly positive");
    }
    if ((bits & UINT32_C(0x7f800000)) == 0)
    {
        return PBC_Fail(context,
                        field +
                            " is subnormal and would be flushed to zero by an "
                            "accelerator before reciprocal-cell construction");
    }
    *result = static_cast<double>(value);
    return true;
}

static bool PBC_Read_Box_Angle(float value, const char* name,
                               const PBC_Build_Context& context,
                               double* degrees)
{
    const std::string field = std::string("box angle ") + name;
    if (!Float_Memory_Is_Finite(&value))
    {
        return PBC_Fail(context, field + " is non-finite");
    }
    const std::uint32_t bits = PBC_Float_Bits(value);
    const std::uint32_t magnitude = bits & UINT32_C(0x7fffffff);
    if ((bits & UINT32_C(0x80000000)) != 0 || magnitude == 0 ||
        (bits & UINT32_C(0x7f800000)) == 0 ||
        !(static_cast<double>(value) < 180.0))
    {
        return PBC_Fail(context,
                        field +
                            " must be finite and strictly between 0 and 180 "
                            "degrees");
    }
    *degrees = static_cast<double>(value);
    return true;
}

static bool PBC_Checked_Float(double value, const std::string& field,
                              bool allow_zero, const PBC_Build_Context& context,
                              float* result)
{
    if (!Double_Memory_Is_Finite(&value))
    {
        return PBC_Fail(context, field + " is non-finite");
    }
    if (PBC_Double_Is_Zero(value))
    {
        if (!allow_zero)
        {
            return PBC_Fail(context, field + " must be nonzero");
        }
        *result = 0.0f;
        return true;
    }
    const double float_max =
        static_cast<double>(std::numeric_limits<float>::max());
    if (value > float_max || value < -float_max)
    {
        return PBC_Fail(context, field + " is outside the finite float range");
    }
    const float stored = static_cast<float>(value);
    if (!Float_Memory_Is_Finite(&stored))
    {
        return PBC_Fail(context, field + " is outside the finite float range");
    }
    const std::uint32_t bits = PBC_Float_Bits(stored);
    if ((bits & UINT32_C(0x7f800000)) == 0)
    {
        return PBC_Fail(
            context,
            field +
                " cannot be represented as a nonzero normal float and would "
                "be flushed to zero by an accelerator");
    }
    *result = stored;
    return true;
}

static bool PBC_Checked_Positive_Float(double value, const std::string& field,
                                       const PBC_Build_Context& context,
                                       float* result)
{
    if (!Double_Memory_Is_Finite(&value) ||
        !PBC_Double_Is_Strictly_Positive(value))
    {
        return PBC_Fail(context,
                        field + " must be finite and strictly positive");
    }
    if (!PBC_Checked_Float(value, field, false, context, result)) return false;
    const std::uint32_t bits = PBC_Float_Bits(*result);
    if ((bits & UINT32_C(0x80000000)) != 0)
    {
        return PBC_Fail(context,
                        field + " must be finite and strictly positive");
    }
    return true;
}

static bool PBC_Raw_Cell_From_Box(VECTOR box_length, VECTOR box_angle,
                                  const PBC_Build_Context& context,
                                  PBC_Raw_Cell* raw_cell)
{
    const char* axes[3] = {"x", "y", "z"};
    const char* angle_names[3] = {"alpha", "beta", "gamma"};
    const float length_values[3] = {box_length.x, box_length.y, box_length.z};
    const float angle_values[3] = {box_angle.x, box_angle.y, box_angle.z};
    double lengths[3] = {};
    double degrees[3] = {};
    for (int axis = 0; axis < 3; axis++)
    {
        if (!PBC_Read_Positive_Box_Length(length_values[axis], axes[axis],
                                          context, &lengths[axis]))
        {
            return false;
        }
        if (!PBC_Read_Box_Angle(angle_values[axis], angle_names[axis], context,
                                &degrees[axis]))
        {
            return false;
        }
    }

    double radians[3] = {};
    double cosines[3] = {};
    for (int angle = 0; angle < 3; angle++)
    {
        radians[angle] = CONSTANT_DEG_TO_RAD_DOUBLE * degrees[angle];
        cosines[angle] =
            degrees[angle] == 90.0 ? 0.0 : std::cos(radians[angle]);
        if (!Double_Memory_Is_Finite(&radians[angle]) ||
            !Double_Memory_Is_Finite(&cosines[angle]))
        {
            return PBC_Fail(context,
                            std::string("box angle ") + angle_names[angle] +
                                " cannot be converted to finite trigonometric "
                                "cell parameters");
        }
    }

    const double sin_gamma = degrees[2] == 90.0 ? 1.0 : std::sin(radians[2]);
    float checked_sin_gamma = 0.0f;
    if (!PBC_Checked_Positive_Float(sin_gamma, "sin(gamma)", context,
                                    &checked_sin_gamma))
    {
        return false;
    }

    const double cos_alpha = cosines[0];
    const double cos_beta = cosines[1];
    const double cos_gamma = cosines[2];
    const double gram_determinant =
        1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta -
        cos_gamma * cos_gamma + 2.0 * cos_alpha * cos_beta * cos_gamma;
    if (!Double_Memory_Is_Finite(&gram_determinant) ||
        !PBC_Double_Is_Strictly_Positive(gram_determinant))
    {
        return PBC_Fail(
            context,
            "box angles have a non-finite or non-positive triclinic Gram "
            "determinant (1-ca^2-cb^2-cg^2+2*ca*cb*cg)");
    }
    const double gram_sqrt = std::sqrt(gram_determinant);
    if (!Double_Memory_Is_Finite(&gram_sqrt) ||
        !PBC_Double_Is_Strictly_Positive(gram_sqrt))
    {
        return PBC_Fail(
            context,
            "the triclinic Gram determinant has no finite positive square "
            "root");
    }

    raw_cell->a11 = lengths[0];
    raw_cell->a21 = lengths[1] * cos_gamma;
    raw_cell->a22 = lengths[1] * sin_gamma;
    raw_cell->a31 = lengths[2] * cos_beta;
    raw_cell->a32 = lengths[2] / sin_gamma * (cos_alpha - cos_beta * cos_gamma);
    raw_cell->a33 = lengths[2] / sin_gamma * gram_sqrt;
    return true;
}

static bool PBC_Check_Volume(const LTMatrix3& cell,
                             const PBC_Build_Context& context)
{
    const double diagonals[3] = {static_cast<double>(cell.a11),
                                 static_cast<double>(cell.a22),
                                 static_cast<double>(cell.a33)};
    const char* pair_names[3] = {"a11*a22", "a11*a33", "a22*a33"};
    const int first[3] = {0, 0, 1};
    const int second[3] = {1, 2, 2};
    const int remaining[3] = {2, 1, 0};
    for (int association = 0; association < 3; association++)
    {
        float pair = 0.0f;
        if (!PBC_Checked_Positive_Float(
                diagonals[first[association]] * diagonals[second[association]],
                std::string("cell volume intermediate ") +
                    pair_names[association],
                context, &pair))
        {
            return false;
        }
        float volume = 0.0f;
        if (!PBC_Checked_Positive_Float(
                static_cast<double>(pair) * diagonals[remaining[association]],
                "cell volume", context, &volume))
        {
            return false;
        }
        float volume_inverse = 0.0f;
        if (!PBC_Checked_Positive_Float(1.0 / static_cast<double>(volume),
                                        "reciprocal cell volume", context,
                                        &volume_inverse))
        {
            return false;
        }
    }
    return true;
}

static bool PBC_Finalize_Cell(const PBC_Raw_Cell& raw_cell,
                              const PBC_Build_Context& context,
                              PBC_Checked_Cell* result)
{
    const double raw_components[6] = {raw_cell.a11, raw_cell.a21, raw_cell.a22,
                                      raw_cell.a31, raw_cell.a32, raw_cell.a33};
    const char* component_names[6] = {"cell a11", "cell a21", "cell a22",
                                      "cell a31", "cell a32", "cell a33"};
    const bool is_diagonal[6] = {true, false, true, false, false, true};
    float components[6] = {};
    for (int component = 0; component < 6; component++)
    {
        if (is_diagonal[component])
        {
            if (!PBC_Checked_Positive_Float(raw_components[component],
                                            component_names[component], context,
                                            &components[component]))
            {
                return false;
            }
        }
        else if (!PBC_Checked_Float(raw_components[component],
                                    component_names[component], true, context,
                                    &components[component]))
        {
            return false;
        }
    }
    const LTMatrix3 cell(components[0], components[1], components[2],
                         components[3], components[4], components[5]);

    const double a11 = static_cast<double>(cell.a11);
    const double a21 = static_cast<double>(cell.a21);
    const double a22 = static_cast<double>(cell.a22);
    const double a31 = static_cast<double>(cell.a31);
    const double a32 = static_cast<double>(cell.a32);
    const double a33 = static_cast<double>(cell.a33);
    const double expected_inverse[6] = {
        1.0 / a11,          -a21 / (a11 * a22),
        1.0 / a22,          (a21 * a32 - a22 * a31) / (a11 * a22 * a33),
        -a32 / (a22 * a33), 1.0 / a33};
    const char* inverse_names[6] = {
        "reciprocal cell a11", "reciprocal cell a21", "reciprocal cell a22",
        "reciprocal cell a31", "reciprocal cell a32", "reciprocal cell a33"};
    for (int component = 0; component < 6; component++)
    {
        float checked = 0.0f;
        if (is_diagonal[component])
        {
            if (!PBC_Checked_Positive_Float(expected_inverse[component],
                                            inverse_names[component], context,
                                            &checked))
            {
                return false;
            }
        }
        else if (!PBC_Checked_Float(expected_inverse[component],
                                    inverse_names[component], true, context,
                                    &checked))
        {
            return false;
        }
    }

    const LTMatrix3 reciprocal = inv(cell);
    const float reciprocal_components[6] = {reciprocal.a11, reciprocal.a21,
                                            reciprocal.a22, reciprocal.a31,
                                            reciprocal.a32, reciprocal.a33};
    for (int component = 0; component < 6; component++)
    {
        if (!Float_Memory_Is_Finite(&reciprocal_components[component]))
        {
            return PBC_Fail(context,
                            std::string(inverse_names[component]) +
                                " is non-finite after shared inversion");
        }
        const std::uint32_t bits =
            PBC_Float_Bits(reciprocal_components[component]);
        const bool expected_nonzero =
            !PBC_Double_Is_Zero(expected_inverse[component]);
        if (expected_nonzero && (bits & UINT32_C(0x7f800000)) == 0)
        {
            return PBC_Fail(context,
                            std::string(inverse_names[component]) +
                                " became zero or subnormal during shared "
                                "single-precision inversion");
        }
        if (is_diagonal[component] && ((bits & UINT32_C(0x80000000)) != 0 ||
                                       (bits & UINT32_C(0x7fffffff)) == 0))
        {
            return PBC_Fail(context,
                            std::string(inverse_names[component]) +
                                " must remain strictly positive and nonzero");
        }
    }
    if (!PBC_Check_Volume(cell, context)) return false;

    result->cell = cell;
    result->rcell = reciprocal;
    return true;
}

static bool PBC_Checked_Build_From_Box(VECTOR box_length, VECTOR box_angle,
                                       const PBC_Build_Context& context,
                                       PBC_Checked_Cell* result)
{
    PBC_Raw_Cell raw_cell = {};
    return PBC_Raw_Cell_From_Box(box_length, box_angle, context, &raw_cell) &&
           PBC_Finalize_Cell(raw_cell, context, result);
}

static bool PBC_Checked_Box_From_Cell(const LTMatrix3& cell,
                                      const PBC_Build_Context& context,
                                      VECTOR* box_length, VECTOR* box_angle)
{
    const double a11 = static_cast<double>(cell.a11);
    const double a21 = static_cast<double>(cell.a21);
    const double a22 = static_cast<double>(cell.a22);
    const double a31 = static_cast<double>(cell.a31);
    const double a32 = static_cast<double>(cell.a32);
    const double a33 = static_cast<double>(cell.a33);
    const double length_a = a11;
    const double length_b = std::hypot(a21, a22);
    const double length_c = std::hypot(std::hypot(a31, a32), a33);

    const double gamma = std::atan2(a22, a21);
    const double beta = std::atan2(std::hypot(a32, a33), a31);
    const double alpha_cross_x = a22 * a33;
    const double alpha_cross_y = -a21 * a33;
    const double alpha_cross_z = a21 * a32 - a22 * a31;
    const double alpha_cross =
        std::hypot(std::hypot(alpha_cross_x, alpha_cross_y), alpha_cross_z);
    const double alpha_dot = a21 * a31 + a22 * a32;
    const double alpha = std::atan2(alpha_cross, alpha_dot);

    float lengths[3] = {};
    float angles[3] = {};
    const double raw_lengths[3] = {length_a, length_b, length_c};
    const double radians[3] = {alpha, beta, gamma};
    const char* axes[3] = {"x", "y", "z"};
    const char* angle_names[3] = {"alpha", "beta", "gamma"};
    for (int axis = 0; axis < 3; axis++)
    {
        if (!PBC_Checked_Positive_Float(raw_lengths[axis],
                                        std::string("box length ") + axes[axis],
                                        context, &lengths[axis]))
        {
            return false;
        }
        const double degrees = radians[axis] / CONSTANT_DEG_TO_RAD_DOUBLE;
        if (!PBC_Checked_Positive_Float(
                degrees, std::string("box angle ") + angle_names[axis], context,
                &angles[axis]) ||
            !(static_cast<double>(angles[axis]) < 180.0))
        {
            return PBC_Fail(
                context, std::string("box angle ") + angle_names[axis] +
                             " derived from the updated cell is not strictly "
                             "between 0 and 180 degrees");
        }
    }

    const VECTOR candidate_length = {lengths[0], lengths[1], lengths[2]};
    const VECTOR candidate_angle = {angles[0], angles[1], angles[2]};
    PBC_Checked_Cell reconstructed;
    if (!PBC_Checked_Build_From_Box(candidate_length, candidate_angle, context,
                                    &reconstructed))
    {
        return false;
    }
    *box_length = candidate_length;
    *box_angle = candidate_angle;
    return true;
}

}  // namespace

void MD_INFORMATION::periodic_box_condition_information::Initial(
    CONTROLLER* controller, MD_INFORMATION* md_info)
{
    this->controller = controller;
    this->md_info = md_info;
    this->pbc = true;
    if (controller->Command_Exist("pbc"))
    {
        this->pbc = controller->Get_Bool(
            "pbc",
            "MD_INFORMATION::periodic_box_condition_information::Initial");
    }
    this->PBC_Check();
    this->No_PBC_Check(controller);
    this->cell0 = cell;
    this->is_initialized = true;
}

void MD_INFORMATION::periodic_box_condition_information::No_PBC_Check(
    CONTROLLER* controller)
{
    if (this->pbc) return;

    if (controller->MPI_size > 1)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "MD_INFORMATION::periodic_box_condition_information::No_PBC_Check",
            "NOPBC can not be used in Multi-Process mode");
    }

    if (md_info->mode == md_info->NPT)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "MD_INFORMATION::periodic_box_condition_information::No_PBC_Check",
            "NPT mode can not be used for NOPBC");
    }
    if (!(controller->Command_Exist("SITS", "atom_numbers") &&
          (strcmp(controller->Command("SITS", "atom_numbers"), "ITS") == 0 ||
           strcmp(controller->Command("SITS", "atom_numbers"), "ALL") == 0)) &&
        controller->Command_Exist("SITS", "mode"))
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand,
            "MD_INFORMATION::periodic_box_condition_information::No_PBC_Check",
            "SITS can not be used for NOPBC now");
    }

    // All force/CV kernels that take the cell and reciprocal cell compute a
    // displacement as dr - image * cell.  A zero interaction-cell matrix
    // therefore expresses NOPBC geometry exactly: no image can alter dr.
    // reference_cell and rcell retain the validated input cell for volume,
    // stress normalization and explicitly scaled-position CVs.  This removes
    // the historical dependence on arbitrary 100/900-Angstrom thresholds
    // without pretending that the interaction geometry is periodic.
    cell = LTMatrix3();
}

void MD_INFORMATION::periodic_box_condition_information::PBC_Check()
{
    const PBC_Build_Context context = {
        controller, spongeErrorBadFileFormat,
        "MD_INFORMATION::periodic_box_condition_information::Initial",
        "initial periodic-cell construction"};
    PBC_Checked_Cell candidate;
    if (!PBC_Checked_Build_From_Box(md_info->sys.box_length,
                                    md_info->sys.box_angle, context,
                                    &candidate))
    {
        return;
    }
    reference_cell = candidate.cell;
    cell = candidate.cell;
    rcell = candidate.rcell;
}

void MD_INFORMATION::periodic_box_condition_information::Update_Box(LTMatrix3 g)
{
    const double dt = static_cast<double>(md_info->dt);
    const LTMatrix3 source_cell = reference_cell;
    const PBC_Raw_Cell raw_candidate = {
        static_cast<double>(source_cell.a11) +
            dt * static_cast<double>(source_cell.a11) *
                static_cast<double>(g.a11),
        static_cast<double>(source_cell.a21) +
            dt * (static_cast<double>(source_cell.a21) *
                      static_cast<double>(g.a11) +
                  static_cast<double>(source_cell.a22) *
                      static_cast<double>(g.a21)),
        static_cast<double>(source_cell.a22) +
            dt * static_cast<double>(source_cell.a22) *
                static_cast<double>(g.a22),
        static_cast<double>(source_cell.a31) +
            dt * (static_cast<double>(source_cell.a31) *
                      static_cast<double>(g.a11) +
                  static_cast<double>(source_cell.a32) *
                      static_cast<double>(g.a21) +
                  static_cast<double>(source_cell.a33) *
                      static_cast<double>(g.a31)),
        static_cast<double>(source_cell.a32) +
            dt * (static_cast<double>(source_cell.a32) *
                      static_cast<double>(g.a22) +
                  static_cast<double>(source_cell.a33) *
                      static_cast<double>(g.a32)),
        static_cast<double>(source_cell.a33) +
            dt * static_cast<double>(source_cell.a33) *
                static_cast<double>(g.a33)};
    const PBC_Build_Context context = {
        controller, spongeErrorSimulationBreakDown,
        "MD_INFORMATION::periodic_box_condition_information::Update_Box",
        "runtime periodic-cell update"};
    PBC_Checked_Cell candidate;
    if (!PBC_Finalize_Cell(raw_candidate, context, &candidate)) return;

    VECTOR candidate_length;
    VECTOR candidate_angle;
    if (!PBC_Checked_Box_From_Cell(candidate.cell, context, &candidate_length,
                                   &candidate_angle))
    {
        return;
    }

    reference_cell = candidate.cell;
    cell = pbc ? candidate.cell : LTMatrix3();
    rcell = candidate.rcell;
    md_info->sys.box_length = candidate_length;
    md_info->sys.box_angle = candidate_angle;
}

void MD_INFORMATION::periodic_box_condition_information::Update_Box_From_Input(
    VECTOR box_length, VECTOR box_angle)
{
    const PBC_Build_Context context = {
        controller, spongeErrorSimulationBreakDown,
        "MD_INFORMATION::periodic_box_condition_information::"
        "Update_Box_From_Input",
        "runtime trajectory-cell publication"};
    PBC_Checked_Cell candidate;
    if (!PBC_Checked_Build_From_Box(box_length, box_angle, context, &candidate))
    {
        return;
    }

    // The trajectory record is the authoritative public float state.  Do not
    // replace it with lengths/angles reconstructed from the float-valued cell:
    // a valid triclinic round trip can differ by one ulp and would make an
    // unchanged following frame look like another deformation.
    reference_cell = candidate.cell;
    cell = pbc ? candidate.cell : LTMatrix3();
    rcell = candidate.rcell;
    md_info->sys.box_length = box_length;
    md_info->sys.box_angle = box_angle;
}

bool MD_INFORMATION::periodic_box_condition_information::Check_Change_Large()
{
    bool result = false;
    float grid_length = 0.5f * (md_info->nb.cutoff + md_info->nb.skin);
    float* cell = (float*)&this->cell;
    float* cell0 = (float*)&this->cell0;
    int i1, i0;
    float f1, f0;
    for (int i = 0; i < 6; i += 1)
    {
        i1 = cell[i] / grid_length;
        i0 = cell0[i] / grid_length;
        f1 = cell[i];
        f0 = cell0[i];
        if (fabsf(f1 - f0) > 0.5f * md_info->nb.skin && i1 != i0)
        {
            result = true;
        }
    }
    if (result)
    {
        this->cell0 = this->cell;
    }
    return result;
}

LTMatrix3 MD_INFORMATION::periodic_box_condition_information::Get_Cell(
    VECTOR box_length, VECTOR box_angle)
{
    const PBC_Build_Context context = {
        controller, spongeErrorSimulationBreakDown,
        "MD_INFORMATION::periodic_box_condition_information::Get_Cell",
        "runtime periodic-cell rebuild"};
    PBC_Checked_Cell candidate;
    if (!PBC_Checked_Build_From_Box(box_length, box_angle, context, &candidate))
    {
        return LTMatrix3();
    }
    return candidate.cell;
}
