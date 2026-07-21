#ifndef REAXFF_H
#define REAXFF_H

#include "../../Domain_decomposition/Domain_decomposition.h"
#include "../../control.h"
#include "../../neighbor_list/neighbor_list.h"
#include "bond.h"
#include "bond_order.h"
#include "eeq.h"
#include "hydrogen_bond.h"
#include "over_under.h"
#include "torsion.h"
#include "valence_angle.h"
#include "vdw.h"

struct REAXFF
{
    int is_initialized = 0;
    int atom_numbers = 0;
    CONTROLLER* controller = NULL;

    REAXFF_EEQ eeq;
    REAXFF_BOND_ORDER bond_order;
    REAXFF_BOND bond;
    REAXFF_VDW vdw;
    REAXFF_OVER_UNDER ovun;
    REAXFF_VALENCE_ANGLE angle;
    REAXFF_TORSION torsion;
    REAXFF_HYDROGEN_BOND hb;

    // A ReaxFF evaluation is transactional.  Every module accumulates into
    // these private buffers; DD force/energy/virial and local/global charge are
    // published only after the complete evaluation has passed validation.
    VECTOR* d_staged_frc = NULL;
    float* d_staged_energy = NULL;
    LTMatrix3* d_staged_virial = NULL;
    float* d_staged_charge = NULL;
    int staging_capacity = 0;
    int* d_seen_global = NULL;
    int* d_evaluation_error = NULL;

    void Initial(CONTROLLER* controller, int atom_numbers, float cutoff,
                 float* cutoff_full, bool* need_full_nl_flag);
    void Calculate_Force(DOMAIN_INFORMATION* dd, MD_INFORMATION* md_info,
                         NEIGHBOR_LIST* neighbor_list,
                         bool commit_sampling_state = true);
    void Validate_Parallel_Layout(CONTROLLER* controller) const;
    void Get_Local(CONTROLLER* controller, const int* atom_local,
                   int local_atom_numbers, int ghost_numbers);
    void Step_Print(CONTROLLER* controller, const float* d_charge);

   private:
    void Wire_Shared_State();
    void Ensure_Staging_Capacity(int coordinate_count);
    void Validate_Module_Consistency();
    void Preflight_Evaluation(DOMAIN_INFORMATION* dd, MD_INFORMATION* md_info,
                              NEIGHBOR_LIST* neighbor_list);
    void Validate_Staged_Result(DOMAIN_INFORMATION* dd,
                                MD_INFORMATION* md_info,
                                bool commit_sampling_state);
    void Commit_Staged_Result(DOMAIN_INFORMATION* dd, MD_INFORMATION* md_info,
                              bool commit_sampling_state);
    void Check_Evaluation_Error(const char* error_by);
};

#include "reaxff_local.hpp"

#endif
