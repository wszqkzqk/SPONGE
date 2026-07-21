#include "xponge.h"

#include "load/amber.hpp"
#include "load/gromacs.hpp"
#include "load/native.hpp"

void Xponge::System::Load_Inputs(CONTROLLER* controller)
{
    // Source selection itself calls Command_Exist, which updates command
    // bookkeeping and may materialize default-prefix input commands.  Keep
    // those changes provisional until the selected source transaction has
    // published a complete System.
    Load_Controller_Command_Transaction command_transaction(controller);
    const bool use_gromacs = controller->Command_Exist("gromacs_top") ||
                             controller->Command_Exist("gromacs_gro");
    const bool use_amber = controller->Command_Exist("amber_parm7") ||
                           controller->Command_Exist("amber_rst7");
    const bool has_charge_endpoint =
        controller->Command_Exist("charge_A_in_file") ||
        controller->Command_Exist("charge_B_in_file");
    if (use_gromacs && use_amber)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "Xponge::System::Load_Inputs",
            "Reason:\n\tGROMACS and AMBER input sources cannot be selected "
            "together; provide exactly one source family\n");
        return;
    }
    if ((use_gromacs || use_amber) && has_charge_endpoint)
    {
        controller->Throw_SPONGE_Error(
            spongeErrorConflictingCommand, "Xponge::System::Load_Inputs",
            "Reason:\n\t'charge_A_in_file' and 'charge_B_in_file' belong to "
            "the native LJ_soft_core input contract and cannot be combined "
            "with GROMACS or AMBER input selection; provide the endpoints "
            "together with a native LJ_soft_core_in_file\n");
        return;
    }
    if (use_gromacs)
    {
        Load_Gromacs_Inputs(this, controller);
    }
    else if (use_amber)
    {
        Load_Amber_Inputs(this, controller);
    }
    else
    {
        Load_Native_Inputs(this, controller);
    }
    command_transaction.Commit();
}
