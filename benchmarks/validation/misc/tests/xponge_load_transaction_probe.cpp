#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include "xponge/load/common.hpp"

int CONTROLLER::MPI_rank = 0;

namespace
{

struct Abort_Transaction
{
};

bool abort_during_source_selection = false;

bool Has_Original_State(const Xponge::System& system)
{
    return system.source == Xponge::InputSource::kAmber &&
           system.start_time == 7.5 && system.atoms.mass.size() == 1 &&
           system.atoms.mass[0] == 12.0f &&
           system.atoms.charge.size() == 1 &&
           system.atoms.charge[0] == -0.25f &&
           system.box.box_length.size() == 3 &&
           system.box.box_length[0] == 10.0f &&
           system.box.box_length[1] == 11.0f &&
           system.box.box_length[2] == 12.0f;
}

bool Has_Same_Virtual_Atoms(const Xponge::VirtualAtoms& lhs,
                            const Xponge::VirtualAtoms& rhs)
{
    if (lhs.records.size() != rhs.records.size()) return false;
    for (std::size_t i = 0; i < lhs.records.size(); i++)
    {
        if (lhs.records[i].type != rhs.records[i].type ||
            lhs.records[i].virtual_atom != rhs.records[i].virtual_atom ||
            lhs.records[i].from != rhs.records[i].from ||
            lhs.records[i].parameter != rhs.records[i].parameter)
        {
            return false;
        }
    }
    return true;
}

bool Has_Same_Amber_Topology(const Xponge::System& lhs,
                             const Xponge::System& rhs)
{
    const Xponge::ClassicalForceField& left = lhs.classical_force_field;
    const Xponge::ClassicalForceField& right = rhs.classical_force_field;
    return lhs.atoms.mass == rhs.atoms.mass &&
           lhs.atoms.charge == rhs.atoms.charge &&
           lhs.residues.atom_numbers == rhs.residues.atom_numbers &&
           lhs.exclusions.excluded_atoms == rhs.exclusions.excluded_atoms &&
           lhs.generalized_born.radius == rhs.generalized_born.radius &&
           lhs.generalized_born.scale_factor ==
               rhs.generalized_born.scale_factor &&
           Has_Same_Virtual_Atoms(lhs.virtual_atoms, rhs.virtual_atoms) &&
           left.bonds.atom_a == right.bonds.atom_a &&
           left.bonds.atom_b == right.bonds.atom_b &&
           left.bonds.k == right.bonds.k && left.bonds.r0 == right.bonds.r0 &&
           left.constraints.atom_a == right.constraints.atom_a &&
           left.constraints.atom_b == right.constraints.atom_b &&
           left.constraints.r0 == right.constraints.r0 &&
           left.angles.atom_a == right.angles.atom_a &&
           left.angles.atom_b == right.angles.atom_b &&
           left.angles.atom_c == right.angles.atom_c &&
           left.angles.k == right.angles.k &&
           left.angles.theta0 == right.angles.theta0 &&
           left.dihedrals.atom_a == right.dihedrals.atom_a &&
           left.dihedrals.atom_b == right.dihedrals.atom_b &&
           left.dihedrals.atom_c == right.dihedrals.atom_c &&
           left.dihedrals.atom_d == right.dihedrals.atom_d &&
           left.dihedrals.pk == right.dihedrals.pk &&
           left.dihedrals.pn == right.dihedrals.pn &&
           left.dihedrals.ipn == right.dihedrals.ipn &&
           left.dihedrals.gamc == right.dihedrals.gamc &&
           left.dihedrals.gams == right.dihedrals.gams &&
           left.impropers.atom_a == right.impropers.atom_a &&
           left.impropers.atom_b == right.impropers.atom_b &&
           left.impropers.atom_c == right.impropers.atom_c &&
           left.impropers.atom_d == right.impropers.atom_d &&
           left.impropers.pk == right.impropers.pk &&
           left.impropers.pn == right.impropers.pn &&
           left.impropers.ipn == right.impropers.ipn &&
           left.impropers.gamc == right.impropers.gamc &&
           left.impropers.gams == right.impropers.gams &&
           left.nb14.atom_a == right.nb14.atom_a &&
           left.nb14.atom_b == right.nb14.atom_b &&
           left.nb14.A == right.nb14.A && left.nb14.B == right.nb14.B &&
           left.nb14.cf_scale_factor == right.nb14.cf_scale_factor &&
           left.lj.atom_type == right.lj.atom_type &&
           left.lj.pair_A == right.lj.pair_A &&
           left.lj.pair_B == right.lj.pair_B &&
           left.lj.atom_type_numbers == right.lj.atom_type_numbers &&
           left.cmap.atom_a == right.cmap.atom_a &&
           left.cmap.atom_b == right.cmap.atom_b &&
           left.cmap.atom_c == right.cmap.atom_c &&
           left.cmap.atom_d == right.cmap.atom_d &&
           left.cmap.atom_e == right.cmap.atom_e &&
           left.cmap.cmap_type == right.cmap.cmap_type &&
           left.cmap.resolution == right.cmap.resolution &&
           left.cmap.grid_value == right.cmap.grid_value &&
           left.cmap.interpolation_coeff == right.cmap.interpolation_coeff &&
           left.cmap.type_offset == right.cmap.type_offset &&
           left.cmap.unique_type_numbers == right.cmap.unique_type_numbers &&
           left.cmap.unique_gridpoint_numbers ==
               right.cmap.unique_gridpoint_numbers &&
           left.urey_bradley.atom_a == right.urey_bradley.atom_a &&
           left.urey_bradley.atom_b == right.urey_bradley.atom_b &&
           left.urey_bradley.atom_c == right.urey_bradley.atom_c &&
           left.urey_bradley.angle_k == right.urey_bradley.angle_k &&
           left.urey_bradley.angle_theta0 ==
               right.urey_bradley.angle_theta0 &&
           left.urey_bradley.bond_k == right.urey_bradley.bond_k &&
           left.urey_bradley.bond_r0 == right.urey_bradley.bond_r0 &&
           left.lj_soft_core.atom_numbers == right.lj_soft_core.atom_numbers &&
           left.lj_soft_core.LJ_AA == right.lj_soft_core.LJ_AA &&
           left.lj_soft_core.charge_A == right.lj_soft_core.charge_A;
}

}  // namespace

bool CONTROLLER::Command_Exist(const char* key)
{
    command_check[key] = 0;
    if (abort_during_source_selection &&
        std::string(key) == "amber_parm7")
    {
        original_commands["selection_temporary"] = "full value";
        commands["selection_temporary"] = "value";
        command_check["selection_temporary"] = 0;
        choice_check["selection_temporary"] = 3;
        throw Abort_Transaction{};
    }
    return commands.count(key) != 0;
}

bool CONTROLLER::Command_Exist(const char* prefix, const char* key)
{
    const std::string full_key = std::string(prefix) + "_" + key;
    return Command_Exist(full_key.c_str());
}

const char* CONTROLLER::Command(const char* key)
{
    command_check[key] = 0;
    return commands[key].c_str();
}

const char* CONTROLLER::Command(const char* prefix, const char* key)
{
    const std::string full_key = std::string(prefix) + "_" + key;
    return Command(full_key.c_str());
}

const char* CONTROLLER::Original_Command(const char* key)
{
    command_check[key] = 0;
    return original_commands[key].c_str();
}

const char* CONTROLLER::Original_Command(const char* prefix, const char* key)
{
    const std::string full_key = std::string(prefix) + "_" + key;
    return Original_Command(full_key.c_str());
}

void CONTROLLER::Set_Command(const char* flag, const char* value, int check,
                             const char* prefix, bool preserve_full_value)
{
    std::string full_key;
    if (prefix != nullptr && prefix[0] != '\0' &&
        std::string(prefix) != "main")
    {
        full_key = std::string(prefix) + "_";
    }
    full_key += flag;
    original_commands[full_key] = value;
    if (preserve_full_value)
    {
        commands[full_key] = value;
    }
    else
    {
        std::istringstream input(value);
        input >> commands[full_key];
    }
    command_check[full_key] = check;
}

int main(int argc, char* argv[])
{
    if (argc != 5) return EXIT_FAILURE;

    Xponge::System inconsistent_atom_arrays;
    inconsistent_atom_arrays.atoms.mass = {1.0f, 2.0f};
    inconsistent_atom_arrays.atoms.charge = {0.0f};
    if (Xponge::Load_Get_Atom_Numbers(&inconsistent_atom_arrays) != -1)
    {
        return EXIT_FAILURE;
    }

    Xponge::System system;
    system.source = Xponge::InputSource::kAmber;
    system.start_time = 7.5;
    system.atoms.mass = {12.0f};
    system.atoms.charge = {-0.25f};
    system.box.box_length = {10.0f, 11.0f, 12.0f};

    CONTROLLER controller{};
    controller.original_commands = {{"retained", "original value"}};
    controller.commands = {{"retained", "original"}};
    controller.command_check = {{"retained", 1}};
    controller.choice_check = {{"retained", 2}};
    const StringMap original_commands = controller.original_commands;
    const StringMap commands = controller.commands;
    const CheckMap command_check = controller.command_check;
    const CheckMap choice_check = controller.choice_check;

    try
    {
        Xponge::Load_System_Transaction(
            &system, &controller, "transaction probe",
            [&](Xponge::System* staged)
            {
                staged->source = Xponge::InputSource::kGromacs;
                staged->start_time = 99.0;
                staged->atoms.mass = {1.0f, 2.0f};
                staged->box.box_length = {1.0f, 1.0f, 1.0f};
                controller.original_commands["temporary"] = "full value";
                controller.commands["temporary"] = "value";
                controller.command_check["temporary"] = 0;
                controller.choice_check["temporary"] = 3;
                throw Abort_Transaction{};
            });
        return EXIT_FAILURE;
    }
    catch (const Abort_Transaction&)
    {
    }
    if (!Has_Original_State(system)) return EXIT_FAILURE;
    if (controller.original_commands != original_commands ||
        controller.commands != commands ||
        controller.command_check != command_check ||
        controller.choice_check != choice_check)
    {
        return EXIT_FAILURE;
    }

    Xponge::Load_System_Transaction(
        &system, &controller, "transaction probe",
        [&](Xponge::System* staged)
        {
            staged->source = Xponge::InputSource::kNative;
            staged->start_time = 3.0;
            staged->atoms.mass = {4.0f, 5.0f};
            controller.original_commands["committed"] = "full value";
            controller.commands["committed"] = "value";
            controller.command_check["committed"] = 0;
            controller.choice_check["committed"] = 1;
        });
    if (system.source != Xponge::InputSource::kNative ||
        system.start_time != 3.0 || system.atoms.mass.size() != 2 ||
        system.atoms.mass[0] != 4.0f || system.atoms.mass[1] != 5.0f)
    {
        return EXIT_FAILURE;
    }
    if (controller.original_commands.at("committed") != "full value" ||
        controller.commands.at("committed") != "value" ||
        controller.command_check.at("committed") != 0 ||
        controller.choice_check.at("committed") != 1)
    {
        return EXIT_FAILURE;
    }

    Xponge::System empty_seed_system = system;
    empty_seed_system.generalized_born.radius = {9.0f};
    empty_seed_system.virtual_atoms.records.push_back(
        Xponge::VirtualAtomRecord{});
    Xponge::Load_System_Transaction(
        &empty_seed_system, &controller, "empty-seed transaction probe",
        Xponge::Load_System_Seed::kEmpty,
        [](Xponge::System* staged)
        {
            if (staged->source != Xponge::InputSource::kUnknown ||
                !staged->atoms.mass.empty() ||
                !staged->generalized_born.radius.empty() ||
                !staged->virtual_atoms.records.empty())
            {
                throw Abort_Transaction{};
            }
            staged->source = Xponge::InputSource::kGromacs;
            staged->atoms.mass = {6.0f};
        });
    if (empty_seed_system.source != Xponge::InputSource::kGromacs ||
        empty_seed_system.atoms.mass != std::vector<float>({6.0f}) ||
        !empty_seed_system.generalized_born.radius.empty() ||
        !empty_seed_system.virtual_atoms.records.empty())
    {
        return EXIT_FAILURE;
    }

    // Exercise the real System::Load_Inputs entry boundary.  Source probing
    // mutates command bookkeeping before a source loader is selected; an
    // exception there must restore all maps even though no inner source
    // transaction has started yet.
    Xponge::System selection_system = system;
    CONTROLLER selection_controller{};
    selection_controller.original_commands = {{"retained", "full value"}};
    selection_controller.commands = {{"retained", "value"}};
    selection_controller.command_check = {{"retained", 1}};
    selection_controller.choice_check = {{"retained", 2}};
    const StringMap selection_original_commands =
        selection_controller.original_commands;
    const StringMap selection_commands = selection_controller.commands;
    const CheckMap selection_command_check = selection_controller.command_check;
    const CheckMap selection_choice_check = selection_controller.choice_check;
    abort_during_source_selection = true;
    try
    {
        selection_system.Load_Inputs(&selection_controller);
        return EXIT_FAILURE;
    }
    catch (const Abort_Transaction&)
    {
    }
    abort_during_source_selection = false;
    if (selection_system.source != system.source ||
        selection_system.start_time != system.start_time ||
        selection_system.atoms.mass != system.atoms.mass ||
        selection_controller.original_commands !=
            selection_original_commands ||
        selection_controller.commands != selection_commands ||
        selection_controller.command_check != selection_command_check ||
        selection_controller.choice_check != selection_choice_check)
    {
        return EXIT_FAILURE;
    }

    // Two successful Native loads with different atom counts must replace the
    // complete source state.  Optional data planted after the first load must
    // not survive merely because the second input omits its files.
    const fs::path input_directory = argv[1];
    CONTROLLER native_controller{};
    auto configure_native = [&](const char* suffix)
    {
        native_controller.original_commands.clear();
        native_controller.commands.clear();
        native_controller.command_check.clear();
        native_controller.choice_check.clear();
        for (const char* field : {"mass", "charge", "coordinate"})
        {
            const std::string key = std::string(field) + "_in_file";
            const std::string path =
                (input_directory /
                 (std::string(field) + "_" + suffix + ".txt"))
                    .string();
            native_controller.original_commands[key] = path;
            native_controller.commands[key] = path;
        }
    };

    Xponge::System native_system;
    configure_native("first");
    native_system.Load_Inputs(&native_controller);
    native_system.generalized_born.radius = {1.0f, 1.0f};
    native_system.generalized_born.scale_factor = {0.5f, 0.5f};
    Xponge::VirtualAtomRecord virtual_record;
    virtual_record.type = 0;
    virtual_record.virtual_atom = 1;
    virtual_record.from = {0};
    virtual_record.parameter = {1.0f};
    native_system.virtual_atoms.records.push_back(virtual_record);
    native_system.classical_force_field.bonds.atom_a = {0};
    native_system.classical_force_field.bonds.atom_b = {1};
    native_system.classical_force_field.nb14.atom_a = {0};
    native_system.classical_force_field.nb14.atom_b = {1};
    native_system.classical_force_field.nb14.A = {12.0f};
    native_system.classical_force_field.nb14.B = {6.0f};
    native_system.classical_force_field.nb14.cf_scale_factor = {0.5f};
    native_system.classical_force_field.lj_soft_core.atom_numbers = 2;
    native_system.classical_force_field.lj_soft_core.charge_A = {1.0f, 1.0f};

    configure_native("second");
    native_system.Load_Inputs(&native_controller);
    const Xponge::ClassicalForceField& ff =
        native_system.classical_force_field;
    if (native_system.source != Xponge::InputSource::kNative ||
        native_system.start_time != 8.5 ||
        native_system.atoms.mass !=
            std::vector<float>({2.0f, 3.0f, 4.0f}) ||
        native_system.atoms.charge !=
            std::vector<float>({-1.0f, 0.0f, 1.0f}) ||
        native_system.atoms.coordinate.size() != 9 ||
        native_system.atoms.velocity != std::vector<float>(9, 0.0f) ||
        native_system.box.box_length !=
            std::vector<float>({21.0f, 22.0f, 23.0f}) ||
        native_system.box.box_angle !=
            std::vector<float>({80.0f, 90.0f, 100.0f}) ||
        native_system.residues.atom_numbers != std::vector<int>({1, 1, 1}) ||
        native_system.exclusions.excluded_atoms.size() != 3 ||
        !native_system.generalized_born.radius.empty() ||
        !native_system.generalized_born.scale_factor.empty() ||
        !native_system.virtual_atoms.records.empty() || !ff.bonds.atom_a.empty() ||
        !ff.nb14.atom_a.empty() || ff.lj_soft_core.atom_numbers != 0 ||
        !ff.lj_soft_core.charge_A.empty())
    {
        return EXIT_FAILURE;
    }
    for (const auto& excluded : native_system.exclusions.excluded_atoms)
    {
        if (!excluded.empty()) return EXIT_FAILURE;
    }

    // GROMACS top+gro is always a complete source replacement.  A second,
    // smaller system must clear every topology and coordinate field from the
    // first materialization rather than merely resizing selected arrays.
    CONTROLLER gromacs_controller{};
    auto configure_gromacs = [&](const char* suffix)
    {
        gromacs_controller.original_commands.clear();
        gromacs_controller.commands.clear();
        gromacs_controller.command_check.clear();
        gromacs_controller.choice_check.clear();
        const std::string top_path =
            (input_directory /
             (std::string("gromacs_") + suffix + ".top"))
                .string();
        const std::string gro_path =
            (input_directory /
             (std::string("gromacs_") + suffix + ".gro"))
                .string();
        gromacs_controller.original_commands["gromacs_top"] = top_path;
        gromacs_controller.commands["gromacs_top"] = top_path;
        gromacs_controller.original_commands["gromacs_gro"] = gro_path;
        gromacs_controller.commands["gromacs_gro"] = gro_path;
    };

    Xponge::System gromacs_system;
    configure_gromacs("first");
    gromacs_system.Load_Inputs(&gromacs_controller);
    gromacs_system.generalized_born.radius = {1.0f, 1.0f};
    gromacs_system.virtual_atoms.records.push_back(virtual_record);
    gromacs_system.classical_force_field.lj_soft_core.atom_numbers = 2;
    gromacs_system.classical_force_field.lj_soft_core.charge_A = {1.0f, 1.0f};
    configure_gromacs("second");
    gromacs_system.Load_Inputs(&gromacs_controller);
    const Xponge::ClassicalForceField& gromacs_ff =
        gromacs_system.classical_force_field;
    if (gromacs_system.source != Xponge::InputSource::kGromacs ||
        gromacs_system.start_time != 4.25 ||
        gromacs_system.atoms.mass != std::vector<float>({20.0f}) ||
        gromacs_system.atoms.charge !=
            std::vector<float>({0.5f * CONSTANT_SPONGE_CHARGE_SCALE}) ||
        gromacs_system.atoms.coordinate !=
            std::vector<float>({2.0f, 3.0f, 4.0f}) ||
        gromacs_system.atoms.velocity != std::vector<float>(3, 0.0f) ||
        gromacs_system.box.box_length !=
            std::vector<float>({60.0f, 70.0f, 80.0f}) ||
        gromacs_system.box.box_angle !=
            std::vector<float>({90.0f, 90.0f, 90.0f}) ||
        gromacs_system.residues.atom_numbers != std::vector<int>({1}) ||
        gromacs_system.exclusions.excluded_atoms.size() != 1 ||
        !gromacs_system.exclusions.excluded_atoms[0].empty() ||
        !gromacs_system.generalized_born.radius.empty() ||
        !gromacs_system.virtual_atoms.records.empty() ||
        !gromacs_ff.bonds.atom_a.empty() || !gromacs_ff.nb14.atom_a.empty() ||
        gromacs_ff.lj.atom_type != std::vector<int>({0}) ||
        gromacs_ff.lj_soft_core.atom_numbers != 0 ||
        !gromacs_ff.lj_soft_core.charge_A.empty())
    {
        return EXIT_FAILURE;
    }

    // AMBER deliberately supports either half of its source pair.  A restart
    // replacement retains the complete parm7-owned topology, while a parm7
    // replacement retains the complete rst7-owned coordinate state.
    CONTROLLER amber_controller{};
    auto configure_amber = [&](const char* parm7, const char* rst7)
    {
        amber_controller.original_commands.clear();
        amber_controller.commands.clear();
        amber_controller.command_check.clear();
        amber_controller.choice_check.clear();
        if (parm7 != nullptr)
        {
            amber_controller.original_commands["amber_parm7"] = parm7;
            amber_controller.commands["amber_parm7"] = parm7;
        }
        if (rst7 != nullptr)
        {
            amber_controller.original_commands["amber_rst7"] = rst7;
            amber_controller.commands["amber_rst7"] = rst7;
        }
    };

    Xponge::System amber_system;
    configure_amber(argv[2], argv[3]);
    amber_system.Load_Inputs(&amber_controller);
    const Xponge::System amber_topology_snapshot = amber_system;
    const std::vector<float> original_amber_coordinate =
        amber_system.atoms.coordinate;

    configure_amber(nullptr, argv[4]);
    amber_system.Load_Inputs(&amber_controller);
    if (!Has_Same_Amber_Topology(amber_system, amber_topology_snapshot) ||
        amber_system.atoms.coordinate == original_amber_coordinate)
    {
        return EXIT_FAILURE;
    }

    const std::size_t amber_component_count =
        amber_system.atoms.mass.size() * 3;
    amber_system.start_time = 123.5;
    amber_system.atoms.coordinate.assign(amber_component_count, 7.0f);
    amber_system.atoms.velocity.assign(amber_component_count, -2.0f);
    amber_system.box.box_length = {31.0f, 32.0f, 33.0f};
    amber_system.box.box_angle = {70.0f, 80.0f, 90.0f};
    const double retained_start_time = amber_system.start_time;
    const std::vector<float> retained_coordinate =
        amber_system.atoms.coordinate;
    const std::vector<float> retained_velocity = amber_system.atoms.velocity;
    const Xponge::Box retained_box = amber_system.box;
    configure_amber(argv[2], nullptr);
    amber_system.Load_Inputs(&amber_controller);
    if (amber_system.source != Xponge::InputSource::kAmber ||
        amber_system.start_time != retained_start_time ||
        amber_system.atoms.coordinate != retained_coordinate ||
        amber_system.atoms.velocity != retained_velocity ||
        amber_system.box.box_length != retained_box.box_length ||
        amber_system.box.box_angle != retained_box.box_angle)
    {
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
