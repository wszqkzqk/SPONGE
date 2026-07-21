#pragma once

#include <type_traits>

#include "../xponge.h"

namespace Xponge
{

// Input discovery is not observationally const: CONTROLLER::Command_Exist and
// the Command/Original_Command accessors update command_check, and the default
// input-prefix fallback may also materialize entries in the command maps.
// Stage those maps together with System so an exception from a loader cannot
// leave command bookkeeping (or an inferred command such as AMBER's pbc=false)
// committed while the System itself is rolled back.
class Load_Controller_Command_Transaction
{
   public:
    explicit Load_Controller_Command_Transaction(CONTROLLER* controller)
        : controller_(controller)
    {
        if (controller_ == nullptr) return;

        // Complete every potentially throwing copy before changing the live
        // controller.  The following swaps use equal default allocators and
        // are noexcept.
        StringMap staged_original_commands = controller_->original_commands;
        StringMap staged_commands = controller_->commands;
        CheckMap staged_command_check = controller_->command_check;
        CheckMap staged_choice_check = controller_->choice_check;

        saved_original_commands_.swap(controller_->original_commands);
        saved_commands_.swap(controller_->commands);
        saved_command_check_.swap(controller_->command_check);
        saved_choice_check_.swap(controller_->choice_check);

        controller_->original_commands.swap(staged_original_commands);
        controller_->commands.swap(staged_commands);
        controller_->command_check.swap(staged_command_check);
        controller_->choice_check.swap(staged_choice_check);
        active_ = true;
    }

    Load_Controller_Command_Transaction(
        const Load_Controller_Command_Transaction&) = delete;
    Load_Controller_Command_Transaction& operator=(
        const Load_Controller_Command_Transaction&) = delete;

    ~Load_Controller_Command_Transaction()
    {
        if (active_) Rollback();
    }

    void Commit() noexcept { active_ = false; }

   private:
    void Rollback() noexcept
    {
        controller_->original_commands.swap(saved_original_commands_);
        controller_->commands.swap(saved_commands_);
        controller_->command_check.swap(saved_command_check_);
        controller_->choice_check.swap(saved_choice_check_);
        active_ = false;
    }

    CONTROLLER* controller_ = nullptr;
    bool active_ = false;
    StringMap saved_original_commands_;
    StringMap saved_commands_;
    CheckMap saved_command_check_;
    CheckMap saved_choice_check_;
};

static bool Load_Float_Is_Finite(float value)
{
    return Float_Memory_Is_Finite(&value);
}

static int Load_Get_Atom_Numbers(const System* system)
{
    std::size_t atom_numbers = 0;
    bool has_atom_numbers = false;
    auto merge_count = [&](std::size_t count)
    {
        if (count >
            static_cast<std::size_t>(std::numeric_limits<int>::max()))
        {
            return false;
        }
        if (has_atom_numbers && atom_numbers != count) return false;
        atom_numbers = count;
        has_atom_numbers = true;
        return true;
    };

    if (!system->atoms.mass.empty() &&
        !merge_count(system->atoms.mass.size()))
    {
        return -1;
    }
    if (!system->atoms.charge.empty() &&
        !merge_count(system->atoms.charge.size()))
    {
        return -1;
    }
    if (!system->atoms.coordinate.empty())
    {
        if (system->atoms.coordinate.size() % 3 != 0 ||
            !merge_count(system->atoms.coordinate.size() / 3))
        {
            return -1;
        }
    }
    if (!system->atoms.velocity.empty())
    {
        if (system->atoms.velocity.size() % 3 != 0 ||
            !merge_count(system->atoms.velocity.size() / 3))
        {
            return -1;
        }
    }
    return has_atom_numbers ? static_cast<int>(atom_numbers) : 0;
}

static int Load_Ensure_Atom_Numbers(System* system, int atom_numbers,
                                    CONTROLLER* controller,
                                    const char* error_by)
{
    int current_atom_numbers = Load_Get_Atom_Numbers(system);
    if (current_atom_numbers < 0)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe retained atom arrays have inconsistent, "
            "misaligned, or unsupported sizes\n");
    }
    if (current_atom_numbers > 0 && current_atom_numbers != atom_numbers)
    {
        controller->Throw_SPONGE_Error(spongeErrorConflictingCommand, error_by,
                                       "Reason:\n\t'atom_numbers' is different "
                                       "in different input files\n");
    }
    return atom_numbers;
}

static void Load_Reset_Classical_Force_Field(ClassicalForceField* ff)
{
    ff->bonds = Bonds{};
    ff->constraints = DistanceConstraints{};
    ff->angles = Angles{};
    ff->dihedrals = Torsions{};
    ff->impropers = Torsions{};
    ff->nb14 = NB14{};
    ff->lj = LennardJones{};
    ff->cmap = CMap{};
    ff->urey_bradley = UreyBradley{};
    ff->lj_soft_core = LJSoftCore{};
}

enum class Load_System_Seed
{
    kEmpty,
    kCurrent,
};

template <typename Loader>
static void Load_System_Transaction(System* system, CONTROLLER* controller,
                                    const char* error_by, Load_System_Seed seed,
                                    Loader&& loader)
{
    static_assert(std::is_nothrow_move_assignable<System>::value,
                  "transactional System publication must not throw");
    try
    {
        Load_Controller_Command_Transaction command_transaction(controller);
        System staged;
        // Complete replacements must not copy (or allocate in proportion to)
        // an old System that they will discard.  AMBER alone uses kCurrent so
        // an omitted parm7/rst7 half is retained intentionally.
        if (seed == Load_System_Seed::kCurrent) staged = *system;
        loader(&staged);
        *system = std::move(staged);
        command_transaction.Commit();
    }
    catch (const std::length_error&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorBadFileFormat, error_by,
            "Reason:\n\tthe transactional input state exceeds a host "
            "container limit\n");
        return;
    }
    catch (const std::bad_alloc&)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorMallocFailed, error_by,
            "Reason:\n\tcould not allocate the transactional input state\n");
        return;
    }
}

template <typename Loader>
static void Load_System_Transaction(System* system, CONTROLLER* controller,
                                    const char* error_by, Loader&& loader)
{
    Load_System_Transaction(system, controller, error_by,
                            Load_System_Seed::kCurrent,
                            std::forward<Loader>(loader));
}

}  // namespace Xponge
