#pragma once

#include "./common.hpp"
#include "./native/angle.hpp"
#include "./native/bond.hpp"
#include "./native/cmap.hpp"
#include "./native/dihedral.hpp"
#include "./native/gb.hpp"
#include "./native/improper_dihedral.hpp"
#include "./native/lj.hpp"
#include "./native/lj_soft.hpp"
#include "./native/md_core.hpp"
#include "./native/nb14.hpp"
#include "./native/urey_bradley.hpp"
#include "./native/virtual_atoms.hpp"

namespace Xponge
{

static void Native_Load_Classical_Force_Field(System* system,
                                              CONTROLLER* controller)
{
    Native_Load_Bonds(system, controller);
    Native_Load_Angles(system, controller);
    Native_Load_Dihedrals(system, controller);
    Native_Load_Impropers(system, controller);
    Native_Load_LJ(system, controller);
    Native_Load_NB14(system, controller);
    Native_Load_Urey_Bradley(system, controller);
    Native_Load_CMap(system, controller);
    Native_Load_LJ_Soft_Core(system, controller);
}

void Load_Native_Inputs(System* system, CONTROLLER* controller)
{
    Load_System_Transaction(
        system, controller, "Xponge::Load_Native_Inputs",
        Load_System_Seed::kEmpty,
        [&](System* staged)
        {
            // Native input is a complete source selection: coordinate input
            // is mandatory and omitted mass/charge/residue/exclusion fields
            // receive native defaults rather than inheriting prior values.
            // Build from an empty System so a previous source cannot impose a
            // stale atom count or leak fields that this load did not publish.
            staged->source = InputSource::kNative;
            Native_Load_Mass(staged, controller);
            Native_Load_Charge(staged, controller);
            Native_Load_Coordinate_And_Velocity(staged, controller);
            Native_Load_Residues(staged, controller);
            Native_Load_Exclusions(staged, controller);
            Native_Load_Classical_Force_Field(staged, controller);
            Native_Load_Generalized_Born(staged, controller);
            Native_Load_Virtual_Atoms(staged, controller);
        });
}

}  // namespace Xponge
