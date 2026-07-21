#ifndef MAIN_RUN_H
#define MAIN_RUN_H

#include "Domain_decomposition/Domain_decomposition.h"
#include "Lennard_Jones_force/LJ_soft_core.h"
#include "Lennard_Jones_force/Lennard_Jones_force.h"
#include "Lennard_Jones_force/solvent_LJ.h"
#include "MD_core/MD_core.h"
#include "NO_PBC/Coulomb_Force_No_PBC.h"
#include "NO_PBC/Lennard_Jones_force_No_PBC.h"
#include "NO_PBC/generalized_Born.h"
#include "PM_force/PM_force.h"
#include "SITS/SITS.h"
#include "angle/Urey_Bradley_force.h"
#include "angle/angle.h"
#include "barostat/MC_barostat.h"
#include "barostat/pressure_based_barostat.h"
#include "bias/restrain_cv.h"
#include "bias/sinkmeta.h"
#include "bias/steer.h"
#include "bond/bond.h"
#include "cmap/cmap.h"
#include "collective_variable/collective_variable.h"
#include "common.h"
#include "constrain/constrain.h"
#include "constrain/settle.h"
#include "constrain/shake.h"
#include "control.h"
#include "custom_force/listed_forces.h"
#include "custom_force/pairwise_force.h"
#include "dihedral/dihedral.h"
#include "dihedral/improper_dihedral.h"
#include "manybody/eam.h"
#include "manybody/edip.h"
#include "manybody/reaxff/reaxff.h"
#include "manybody/sw.h"
#include "manybody/tersoff.h"
#include "nb14/nb14.h"
#include "neighbor_list/full_neighbor_list.h"
#include "neighbor_list/neighbor_list.h"
#include "plugin/plugin.h"
#include "quantum_chemistry/quantum_chemistry.h"
#include "restrain/restrain.h"
#include "thermostat/Andersen_thermostat.h"
#include "thermostat/Berendsen_thermostat.h"
#include "thermostat/Bussi_thermostat.h"
#include "thermostat/Middle_Langevin_MD.h"
#include "thermostat/Nose_Hoover_Chain.h"
#include "virtual_atoms/virtual_atoms.h"
#include "wall/hard_wall.h"
#include "wall/soft_wall.h"
#include "xponge/xponge.h"

struct FORCE_EVALUATION_CONTEXT
{
    bool commit_sampling_state;
    bool exact_state;

    FORCE_EVALUATION_CONTEXT(bool commit_sampling_state = true,
                             bool exact_state = false)
        : commit_sampling_state(commit_sampling_state), exact_state(exact_state)
    {
    }
};

void Main_Initial(int argc, char* argv[]);
void Main_Process_Management();
void Main_Calculate_Force(
    const FORCE_EVALUATION_CONTEXT& evaluation = FORCE_EVALUATION_CONTEXT());
void Main_Iteration();
void Main_Print();
void Main_Clear();
void Main_Sync_Dynamic_Targets_To_Controllers();
void Main_Refresh_Local_State(bool rebuild_dd);

void Main_MC_Barostat();
float Main_Box_Change(LTMatrix3 g, int scale_box, int scale_crd, int scale_vel);
float Main_Box_Change_Transactional(
    LTMatrix3 g, int scale_box, int scale_crd, int scale_vel,
    bool commit_box_state, const VECTOR* authoritative_box_length = NULL,
    const VECTOR* authoritative_box_angle = NULL);
float Main_Rerun_Box_Change(LTMatrix3 g, VECTOR box_length, VECTOR box_angle);
void Main_Box_Change_Largely();

#endif
