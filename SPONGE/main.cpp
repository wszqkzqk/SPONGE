#include "main.h"

#define SUBPACKAGE_HINT \
    "SPONGE, for general-purpose molecular dynamics simulations"
#define THERMOSTAT_IS(name)                              \
    (md_info.mode >= md_info.NVT &&                      \
     (controller.Command_Choice("thermostat", (name)) || \
      controller.Command_Choice("thermostat_mode", (name))))
#define BAROSTAT_IS(name)                              \
    (md_info.mode == md_info.NPT &&                    \
     (controller.Command_Choice("barostat", (name)) || \
      controller.Command_Choice("barostat_mode", (name))))

CONTROLLER controller;
Xponge::System Xponge::system;
MD_INFORMATION md_info;
INITIAL_VELOCITY_INFORMATION initial_velocity;
DOMAIN_INFORMATION dd;
MIDDLE_Langevin_INFORMATION middle_langevin;
ANDERSEN_THERMOSTAT_INFORMATION ad_thermo;
BERENDSEN_THERMOSTAT_INFORMATION bd_thermo;
BUSSI_THERMOSTAT_INFORMATION bussi_thermo;
NOSE_HOOVER_CHAIN_INFORMATION nhc;
PRESSURE_BASED_BAROSTAT_INFORMATION press_baro;
MC_BAROSTAT_INFORMATION mc_baro;
NEIGHBOR_LIST neighbor_list;
LENNARD_JONES_INFORMATION lj;
LJ_SOFT_CORE lj_soft;
SOLVENT_LENNARD_JONES solvent_lj;
Particle_Mesh pm;
ANGLE angle;
UREY_BRADLEY urey_bradley;
BOND bond;
CMAP cmap;
DIHEDRAL dihedral;
IMPROPER_DIHEDRAL improper;
NON_BOND_14 nb14;
RESTRAIN_INFORMATION restrain;
CONSTRAIN constrain;
SETTLE settle;
SHAKE shake;
VIRTUAL_INFORMATION vatom;
COLLECTIVE_VARIABLE_CONTROLLER cv_controller;
STEER_CV steer_cv;
RESTRAIN_CV restrain_cv;
META meta;
VORONOI_DETECTOR voronoi_detector;
LISTED_FORCES listed_forces;
PAIRWISE_FORCE pairwise_force;
HARD_WALL hard_wall;
SOFT_WALLS soft_walls;
LENNARD_JONES_NO_PBC_INFORMATION LJ_NOPBC;
COULOMB_FORCE_NO_PBC_INFORMATION CF_NOPBC;
GENERALIZED_BORN_INFORMATION gb;
SITS_INFORMATION sits;
DIHEDRAL sits_dihedral;
NON_BOND_14 sits_nb14;
CMAP sits_cmap;
STILLINGER_WEBER_INFORMATION sw;
EDIP_INFORMATION edip;
EAM_INFORMATION eam;
TERSOFF_INFORMATION tersoff;
REAXFF reaxff;
QUANTUM_CHEMISTRY qc;
SPONGE_PLUGIN plugin;

deviceStream_t main_stream;

static std::uint64_t coordinate_generation_counter = 0;
static std::uint64_t current_coordinate_generation = 0;

static void Main_Advance_Coordinate_Generation(const char* reason)
{
    if (coordinate_generation_counter ==
        std::numeric_limits<std::uint64_t>::max())
    {
        controller.Throw_Formatted_SPONGE_Error(
            spongeErrorOverflow, "Main_Advance_Coordinate_Generation",
            "Reason:\n    coordinate-state generation overflow while %s\n",
            reason);
        return;
    }
    current_coordinate_generation = ++coordinate_generation_counter;
}

static bool Main_Cell_Is_Exactly_Equal(const LTMatrix3& lhs,
                                       const LTMatrix3& rhs)
{
    return lhs.a11 == rhs.a11 && lhs.a21 == rhs.a21 && lhs.a22 == rhs.a22 &&
           lhs.a31 == rhs.a31 && lhs.a32 == rhs.a32 && lhs.a33 == rhs.a33;
}

static bool Main_Update_Neighbor_List(int update)
{
    // NOPBC nonbonded forces use their dedicated all-pairs kernels.  No
    // neighbor list is constructed for that execution mode, so it must not be
    // sent through the PBC neighbor-list lifecycle.  Keep the neighbor-list
    // API strict: an actual call before initialization remains a hard error.
    if (!md_info.pbc.pbc) return false;

    return neighbor_list.Update_With_Overflow_Recovery(
        &controller, dd.atom_local, dd.atom_numbers, dd.ghost_numbers, dd.crd,
        md_info.pbc.cell, md_info.pbc.rcell, md_info.sys.steps, update,
        md_info.nb.d_excluded_list_start, md_info.nb.d_excluded_list,
        md_info.nb.d_excluded_numbers);
}

static void Main_Export_Voronoi_Hit()
{
    const std::string basename = voronoi_detector.Hit_Restart_Basename();
    if (basename.empty())
    {
        controller.Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "Main_Export_Voronoi_Hit",
            "Reason:\n\tthe terminal Voronoi hit has no destination "
            "artifact name\n");
    }

    // A hit is observed on the committed x_n state, before Main_Iteration.
    // Gather that exact local state and force a host refresh; ordinary output
    // may already have synchronized another state under the same step label.
    md_info.Crd_Vel_dd_to_Device(dd.crd, dd.vel, dd.atom_local_label,
                                 dd.atom_local_id, main_stream);
    deviceStreamSynchronize(main_stream);
    md_info.Crd_Vel_Device_To_Host(true);

    const int completed_steps =
        voronoi_detector.hit_step - voronoi_detector.initial_step;
    const VORONOI_INTERFACE_RECORD& source =
        voronoi_detector.interfaces[voronoi_detector.source_interface];
    const VORONOI_INTERFACE_RECORD& destination =
        voronoi_detector.Hit_Interface();
    controller.printf(
        "VORONOI_HIT source=%s from=%d destination=%s to=%d "
        "completed_steps=%d hit_time_ps=%.10g source_recrossings=%llu "
        "artifact=%s\n",
        source.name.c_str(), voronoi_detector.hit_from_milestone,
        destination.name.c_str(), voronoi_detector.destination_milestone,
        completed_steps,
        static_cast<double>(completed_steps) * md_info.sys.dt_in_ps,
        static_cast<unsigned long long>(
            voronoi_detector.source_recrossing_count),
        basename.c_str());
    md_info.output.Export_Restart_File(basename.c_str(), false);
}

int main(int argc, char* argv[])
{
    Main_Initial(argc, argv);
    for (md_info.sys.steps = 0; md_info.sys.steps <= md_info.sys.step_limit;)
    {
        Main_Sync_Dynamic_Targets_To_Controllers();
        const bool mc_attempt = mc_baro.Will_Attempt(md_info.sys.steps);
        if (mc_attempt)
        {
            Main_MC_Barostat();
        }
        // MC old/trial evaluations never commit adaptive sampling history.
        // The accepted state is the sole committed sample for this physical
        // step and is evaluated exactly after any box transaction.
        Main_Calculate_Force(FORCE_EVALUATION_CONTEXT(true, mc_attempt));
        if (voronoi_detector.Has_Terminal_Hit())
        {
            Main_Export_Voronoi_Hit();
            break;
        }
        Main_Iteration();
        Main_Print();
        // Keep int-valued public/plugin step counters for compatibility, but
        // do not overflow the loop increment at the supported upper bound.
        if (md_info.sys.steps == std::numeric_limits<int>::max()) break;
        md_info.sys.steps++;
    }
    Main_Clear();
    return 0;
}

void Main_Initial(int argc, char* argv[])
{
    controller.Initial(argc, argv, SUBPACKAGE_HINT);
    Xponge::system.Load_Inputs(&controller);
    cv_controller.Initial(&controller,
                          &md_info.no_direct_interaction_virtual_atom_numbers);
    md_info.Initial(&controller);
    // A schedule point at step zero defines the initial thermodynamic state,
    // including values consumed while temperature-dependent modules load
    // persistent history.
    md_info.sys.Update_Targets_By_Schedule(&controller, 0);
    controller.Step_Print_Initial("potential", "%.2f");
    controller.Step_Print_Initial("eff_pot", "%.7e");
    qc.Initial(&controller, md_info.atom_numbers, md_info.crd,
               md_info.sys.box_length, NULL, current_coordinate_generation);
    cv_controller.atom_numbers = md_info.atom_numbers;
    plugin.Initial(&md_info, &controller, &cv_controller, &neighbor_list);

    if (md_info.mode >= md_info.NVT &&
        (!controller.Command_Exist("thermostat") &&
         !controller.Command_Exist("thermostat_mode")))
    {
        controller.Throw_SPONGE_Error(
            spongeErrorMissingCommand, "Main_Initial",
            "Reason:\n\tthermostat is required for NVT or NPT simulations\n");
    }
    if (THERMOSTAT_IS("middle_langevin") || THERMOSTAT_IS("langevin"))
    {
        middle_langevin.Initial(&controller, md_info.atom_numbers,
                                md_info.sys.target_temperature, md_info.h_mass);
    }
    else if (THERMOSTAT_IS("andersen"))
    {
        ad_thermo.Initial(&controller, md_info.sys.target_temperature,
                          md_info.atom_numbers, md_info.sys.dt_in_ps,
                          md_info.h_mass);
    }
    else if (THERMOSTAT_IS("bussi_thermostat"))
    {
        bussi_thermo.Initial(&controller, md_info.sys.target_temperature);
    }
    else if (THERMOSTAT_IS("berendsen_thermostat"))
    {
        bd_thermo.Initial(&controller, md_info.sys.target_temperature);
    }
    else if (THERMOSTAT_IS("nose_hoover_chain"))
    {
        nhc.Initial(&controller, md_info.atom_numbers,
                    md_info.sys.target_temperature, md_info.h_mass);
    }

    if (md_info.mode == md_info.NPT && !controller.Command_Exist("barostat") &&
        !controller.Command_Exist("barostat_mode"))
    {
        controller.Throw_SPONGE_Error(
            spongeErrorMissingCommand, "Main_Initial",
            "Reason:\n\tbarostat is required for NPT simulations\n");
    }
    if (BAROSTAT_IS("andersen_barostat") || BAROSTAT_IS("bussi_barostat") ||
        BAROSTAT_IS("berendsen_barostat"))
    {
        press_baro.Initial(&controller, md_info.sys.target_pressure,
                           md_info.pbc.cell, &Main_Box_Change);
    }
    if (BAROSTAT_IS("monte_carlo_barostat"))
    {
        mc_baro.Initial(&controller, md_info.atom_numbers,
                        md_info.sys.target_pressure, md_info.sys.box_length,
                        md_info.pbc.cell);
    }

    if (md_info.pbc.pbc)
    {
        lj.Initial(&controller, md_info.nb.cutoff);
        lj_soft.Initial(&controller, md_info.nb.cutoff);
        pm.Initial(&controller, md_info.atom_numbers, md_info.pbc.cell,
                   md_info.pbc.rcell, md_info.sys.box_length, md_info.nb.cutoff,
                   md_info.no_direct_interaction_virtual_atom_numbers);
        pairwise_force.Initial(&controller);
        nb14.Initial(&controller, lj.h_LJ_A, lj.h_LJ_B, lj.h_atom_LJ_type,
                     md_info.atom_numbers, lj.atom_type_numbers);

        sits.Initial(&controller, md_info.atom_numbers);
        if (sits.is_initialized && sits.selectively_applied)
        {
            sits_dihedral.Initial(&controller, "sits_dihedral");
            sits_nb14.Initial(&controller, lj.h_LJ_A, lj.h_LJ_B,
                              lj.h_atom_LJ_type, md_info.atom_numbers,
                              lj.atom_type_numbers, "sits_nb14");
            sits_cmap.Initial(&controller, "sits_cmap");
        }
        sits.Check_Solvent(&controller, md_info.atom_numbers,
                           solvent_lj.solvent_numbers);
    }
    else
    {
        LJ_NOPBC.Initial(&controller);
        CF_NOPBC.Initial(&controller, md_info.atom_numbers);
        if (controller.Command_Exist("gb", "in_file"))
        {
            gb.Initial(&controller, md_info.atom_numbers, md_info.nb.cutoff);
        }
        nb14.Initial(&controller, LJ_NOPBC.h_LJ_A, LJ_NOPBC.h_LJ_B,
                     LJ_NOPBC.h_atom_LJ_type, md_info.atom_numbers,
                     LJ_NOPBC.atom_type_numbers);
        sits.Initial(&controller, md_info.atom_numbers);
    }

    bond.Initial(&controller, &md_info.sys.connectivity,
                 &md_info.sys.connected_distance);
    angle.Initial(&controller);
    urey_bradley.Initial(&controller);
    cmap.Initial(&controller);
    dihedral.Initial(&controller);
    improper.Initial(&controller);
    listed_forces.Initial(&controller, &md_info.sys.connectivity,
                          &md_info.sys.connected_distance);

    sw.Initial(&controller, "SW", &neighbor_list.is_needed_full);
    edip.Initial(&controller, "EDIP", &neighbor_list.is_needed_full);
    eam.Initial(&controller, md_info.atom_numbers, "EAM",
                &neighbor_list.is_needed_full);
    tersoff.Initial(&controller, md_info.atom_numbers, "TERSOFF",
                    &neighbor_list.is_needed_full);
    reaxff.Initial(&controller, md_info.atom_numbers, md_info.nb.cutoff,
                   &neighbor_list.cutoff_full, &neighbor_list.is_needed_full);

    restrain.Initial(&controller, md_info.atom_numbers, md_info.crd);
    hard_wall.Initial(&controller, md_info.sys.target_temperature,
                      md_info.sys.target_pressure, md_info.mode == md_info.NPT);
    soft_walls.Initial(&controller, md_info.atom_numbers);

    if (controller.Command_Exist("constrain_mode"))
    {
        constrain.Initial_List(&controller, md_info.sys.connected_distance,
                               md_info.h_mass);
        constrain.Initial_Constrain(&controller, md_info.atom_numbers,
                                    md_info.dt, md_info.sys.box_length,
                                    md_info.h_mass, &md_info.sys.freedom);
        settle.Initial(&controller, &constrain, md_info.h_mass);
        if (controller.Command_Choice("constrain_mode", "SHAKE"))
        {
            shake.Initial_SHAKE(&controller, &constrain);
        }
        if (md_info.mode == md_info.MINIMIZATION)
        {
            constrain.v_factor = 0.0f;
        }
        if (middle_langevin.is_initialized)
        {
            constrain.v_factor = middle_langevin.exp_gamma;
            constrain.x_factor = 0.5 * (1. + middle_langevin.exp_gamma);
        }
    }
    vatom.Initial(&controller, &cv_controller, md_info.atom_numbers,
                  md_info.no_direct_interaction_virtual_atom_numbers,
                  cv_controller.cv_vatom_name, md_info.h_mass,
                  &md_info.sys.freedom, &md_info.sys.connectivity);
    vatom.Coordinate_Refresh(md_info.crd, md_info.pbc.cell, md_info.pbc.rcell);
    initial_velocity.Initial(&controller, &md_info);

    if (md_info.pbc.pbc)
    {
        neighbor_list.Initial(&controller, md_info.atom_numbers,
                              md_info.nb.cutoff, md_info.nb.skin,
                              md_info.pbc.cell, md_info.pbc.rcell);
    }
    steer_cv.Initial(&controller, &cv_controller);
    restrain_cv.Initial(&controller, &cv_controller);
    meta.Initial(&controller, &cv_controller, NULL,
                 md_info.sys.target_temperature);
    voronoi_detector.Initial(&controller, &cv_controller);

    cv_controller.Print_Initial();
    plugin.After_Initial();
    cv_controller.Input_Check();

    md_info.ug.Initial_Edge(md_info.atom_numbers);
    constrain.update_ug_connectivity(&md_info.ug.connectivity);
    settle.update_ug_connectivity(&md_info.ug.connectivity);
    vatom.update_ug_connectivity(&md_info.ug.connectivity);
    md_info.ug.Read_Update_Group(md_info.atom_numbers);
    md_info.mol.Initial(&controller);
    if (md_info.pbc.pbc)
    {
        solvent_lj.Initial(&controller, &lj, &lj_soft, &md_info,
                           md_info.mode >= md_info.NVT);
    }
    Main_Process_Management();
    reaxff.Validate_Parallel_Layout(&controller);

    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        Main_Refresh_Local_State(true);
        initial_velocity.Finalize(&controller, &md_info, &dd, &settle, &shake);
        plugin.Set_Domain_Information(&dd);
    }

    pm.Get_Atoms(&controller, md_info.crd, md_info.d_charge, dd.atom_numbers,
                 dd.crd, dd.d_charge, dd.atom_local, true, true, true, true);

    controller.Print_First_Line_To_Mdout();
}

void Main_Calculate_Force(const FORCE_EVALUATION_CONTEXT& evaluation)
{
    bool use_reaxff_eeq = reaxff.eeq.is_initialized;
    const int cv_atom_numbers =
        md_info.atom_numbers +
        md_info.no_direct_interaction_virtual_atom_numbers;
    cv_controller.Invalidate_Evaluation_Caches();
    md_info.MD_Reset_Atom_Energy_And_Virial_And_Force();
    if (qc.is_initialized && CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        // DD storage is local-order and can change after migration/remeshing.
        // QC atom IDs are global, so every PP rank solves from the same
        // explicitly gathered global-order coordinates.
        md_info.Crd_Vel_dd_to_Device(dd.crd, dd.vel, dd.atom_local_label,
                                     dd.atom_local_id, main_stream);
        qc.Solve_SCF(md_info.crd, md_info.sys.box_length, true,
                     md_info.sys.steps, evaluation.commit_sampling_state,
                     current_coordinate_generation);
    }
    if (md_info.mode == md_info.MINIMIZATION && md_info.min.dynamic_dt)
    {
        md_info.need_potential = 1;
    }
    mc_baro.Ask_For_Calculate_Potential(md_info.sys.steps,
                                        &md_info.need_potential);
    press_baro.Ask_For_Calculate_Pressure(md_info.sys.steps,
                                          &md_info.need_pressure);
    if (press_baro.is_initialized && md_info.output.Check_Mdout_Step())
    {
        md_info.need_pressure = 1;
    }
    if (bd_thermo.is_initialized || bussi_thermo.is_initialized ||
        nhc.is_initialized)
    {
        md_info.need_kinetic = 1;
    }
    sits.Reset_Force_Energy(&md_info.need_potential);

    controller.Get_Time_Recorder("Calculate_Force")->Start();
    pm.Get_Atoms(&controller, md_info.crd, md_info.d_charge, dd.atom_numbers,
                 dd.crd, dd.d_charge, dd.atom_local, false, false, true, false);
    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        dd.Reset_Force_and_Virial(&md_info);
        // QC 梯度必须在 dd.Reset_Force_and_Virial 之后调用
        if (qc.is_initialized && qc.need_gradient)
            qc.Compute_Gradient(dd.frc, md_info.crd, dd.atom_local_id,
                                dd.atom_numbers, md_info.sys.box_length,
                                md_info.need_pressure, dd.d_virial);
        if (qc.is_initialized && md_info.need_potential)
            qc.Accumulate_Energy(dd.d_energy, dd.atom_local_id,
                                 dd.atom_numbers);
        dd.Update_Ghost(&controller);
        const bool neighbor_storage_rebound = Main_Update_Neighbor_List(
            evaluation.exact_state ? neighbor_list.FORCED_UPDATE
                                   : neighbor_list.CONDITIONAL_UPDATE);
        if (neighbor_storage_rebound)
        {
            // Plugins may expose zero-copy views of neighbor storage.  Any
            // recovery reallocation therefore requires an explicit rebind
            // before their force callback is allowed to run.
            plugin.Set_Domain_Information(&dd);
        }

        reaxff.Calculate_Force(&dd, &md_info, &neighbor_list,
                               evaluation.commit_sampling_state);

        LJ_NOPBC.LJ_Force_With_Atom_Energy(
            dd.atom_numbers, dd.crd, dd.frc, md_info.need_potential,
            dd.d_energy, dd.d_excluded_list_start, dd.d_excluded_list,
            dd.d_excluded_numbers);
        CF_NOPBC.Coulomb_Force_With_Atom_Energy(
            dd.atom_numbers, dd.crd, dd.d_charge, dd.frc,
            md_info.need_potential, dd.d_energy, dd.d_excluded_list_start,
            dd.d_excluded_list, dd.d_excluded_numbers);
        gb.GB_Force_With_Atom_Energy(dd.atom_numbers, dd.crd, dd.d_charge,
                                     dd.frc, dd.d_energy);

        if (!use_reaxff_eeq)
        {
            pm.MPI_PME_Excluded_Force_With_Atom_Energy(
                dd.atom_numbers, dd.crd, md_info.pbc.cell, md_info.pbc.rcell,
                dd.d_charge, dd.d_excluded_list_start, dd.d_excluded_list,
                dd.d_excluded_numbers, dd.d_exclusion_dependency_state, dd.frc,
                md_info.need_potential, dd.d_energy, md_info.need_pressure,
                dd.d_virial);
        }

        // S3：tile kernel 覆盖全部 local i（含溶剂原子）的半表对，启用时
        // 必须跳过下方的溶剂 dispatch，否则水-水对双计（tile 路径等价于强制
        // solvent_numbers=0 由 tile kernel 统一覆盖）。SITS 选择性施加时主
        // kernel 不参与非键计算，tile 路径必须停用，溶剂 kernel 保持原分工。
        LJ_TILE_SET lj_tile_set;
        lj_tile_set.tiles = neighbor_list.d_lj_tiles;
        lj_tile_set.cluster_atoms = neighbor_list.d_lj_cluster_atoms;
        lj_tile_set.cluster_flags = neighbor_list.d_lj_cluster_flags;
        lj_tile_set.tile_sorted = neighbor_list.d_lj_tile_sorted;
        lj_tile_set.tile_numbers = neighbor_list.h_lj_tile_numbers;
        lj_tile_set.cluster_atom_slots = neighbor_list.h_lj_cluster_atom_slots;
        const bool lj_tile_active =
            lj.use_tile && !(sits.is_initialized && sits.selectively_applied) &&
            lj_tile_set.tiles != NULL && lj_tile_set.cluster_atoms != NULL;
        if (sits.is_initialized && sits.selectively_applied)
        {
            sits_dihedral.Dihedral_Force_With_Atom_Energy_And_Virial(
                dd.crd, md_info.pbc.cell, md_info.pbc.rcell,
                sits.pw_select.select_force[0], md_info.need_potential,
                sits.pw_select.select_atom_energy[0], md_info.need_pressure,
                sits.pw_select.select_atom_virial_tensor[0]);
            sits_nb14.Non_Bond_14_LJ_CF_Force_With_Atom_Energy_And_Virial(
                dd.crd, dd.d_charge, md_info.pbc.cell, md_info.pbc.rcell,
                sits.pw_select.select_force[0], md_info.need_potential,
                sits.pw_select.select_atom_energy[0], md_info.need_pressure,
                sits.pw_select.select_atom_virial_tensor[0]);
            sits_cmap.CMAP_Force_With_Atom_Energy_And_Virial(
                dd.crd, md_info.pbc.cell, md_info.pbc.rcell,
                sits.pw_select.select_force[0], md_info.need_potential,
                sits.pw_select.select_atom_energy[0], md_info.need_pressure,
                sits.pw_select.select_atom_virial_tensor[0]);
            sits.SITS_LJ_Direct_CF_Force_With_Atom_Energy_And_Virial(
                md_info.atom_numbers, dd.atom_numbers,
                solvent_lj.local_solvent_numbers, dd.ghost_numbers, dd.crd,
                dd.d_charge, &lj, dd.frc, md_info.pbc.cell, md_info.pbc.rcell,
                neighbor_list.d_nl, md_info.nb.cutoff, pm.beta,
                md_info.need_potential, dd.d_energy, md_info.need_pressure,
                dd.d_virial, pm.d_direct_atom_energy);
            sits.SITS_LJ_Soft_Core_Direct_CF_Force_With_Atom_Energy_And_Virial(
                md_info.atom_numbers, dd.atom_numbers,
                solvent_lj.local_solvent_numbers, dd.ghost_numbers, dd.crd,
                dd.d_charge, &lj_soft, dd.frc, md_info.pbc.cell,
                md_info.pbc.rcell, neighbor_list.d_nl, md_info.nb.cutoff,
                pm.beta, md_info.need_potential, dd.d_energy,
                md_info.need_pressure, dd.d_virial, pm.d_direct_atom_energy);
        }
        else
        {
            lj.LJ_PME_Direct_Force_With_Atom_Energy_And_Virial(
                md_info.atom_numbers, dd.atom_numbers,
                solvent_lj.local_solvent_numbers, dd.ghost_numbers, dd.crd,
                dd.d_charge, dd.frc, md_info.pbc.cell, md_info.pbc.rcell,
                neighbor_list.d_nl, pm.beta, md_info.need_potential,
                dd.d_energy, md_info.need_pressure, dd.d_virial,
                pm.d_direct_atom_energy,
                lj_tile_active ? &lj_tile_set : NULL);

            lj_soft.LJ_Soft_Core_PME_Direct_Force_With_Atom_Energy_And_Virial(
                md_info.atom_numbers, dd.atom_numbers,
                solvent_lj.local_solvent_numbers, dd.ghost_numbers, dd.crd,
                dd.d_charge, dd.frc, md_info.pbc.cell, md_info.pbc.rcell,
                neighbor_list.d_nl, pm.beta, md_info.need_potential,
                dd.d_energy, md_info.need_pressure, dd.d_virial,
                pm.d_direct_atom_energy);
        }
        if (!lj_tile_active)
        {
            solvent_lj.LJ_PME_Direct_Force_With_Atom_Energy_And_Virial(
                dd.atom_numbers, dd.res_numbers, dd.res_start, dd.crd,
                dd.d_charge, dd.frc, md_info.pbc.cell, md_info.pbc.rcell,
                neighbor_list.d_nl, pm.beta, md_info.need_potential,
                dd.d_energy, md_info.need_pressure, dd.d_virial,
                pm.d_direct_atom_energy);
        }

        // 单进程路径：LJ 直接空间 kernel 占满整机，实测与之重叠会因
        // SM/缓存争用把它拖慢得更多；把 PME 倒易链投在这里（LJ 之后、
        // 成键 kernel 之前），让它与下方的小 grid 成键 kernel 重叠。
        // force_backup 的累加与能量尾部仍由后面的 Join 在原调用点按原
        // 顺序完成
        if (CONTROLLER::MPI_size == 1 && !use_reaxff_eeq)
        {
            pm.PME_Reciprocal_Force_Async_Start(
                dd.crd, md_info.pbc.cell, md_info.pbc.rcell, dd.d_charge,
                md_info.need_pressure, md_info.need_potential, dd.d_virial,
                md_info.sys.steps, evaluation.exact_state);
        }

        lj.Long_Range_Correction(
            md_info.need_pressure, dd.d_virial, md_info.need_potential,
            dd.d_energy,
            md_info.pbc.cell.a11 * md_info.pbc.cell.a22 * md_info.pbc.cell.a33);

        lj_soft.Long_Range_Correction(
            md_info.need_pressure, dd.d_virial, md_info.need_potential,
            dd.d_energy,
            md_info.pbc.cell.a11 * md_info.pbc.cell.a22 * md_info.pbc.cell.a33);
        sw.SW_Force_With_Atom_Energy_And_Virial_Full_NL(
            dd.atom_numbers, dd.crd, dd.frc, md_info.pbc.cell,
            md_info.pbc.rcell, neighbor_list.full_neighbor_list.d_nl,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        edip.EDIP_Force_With_Atom_Energy_And_Virial_Full_NL(
            dd.atom_numbers, dd.crd, dd.frc, md_info.pbc.cell,
            md_info.pbc.rcell, neighbor_list.full_neighbor_list.d_nl,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        eam.EAM_Force_With_Atom_Energy_And_Virial(
            dd.atom_numbers, dd.crd, dd.frc, md_info.pbc.cell,
            md_info.pbc.rcell, neighbor_list.full_neighbor_list.d_nl,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        tersoff.TERSOFF_Force_With_Atom_Energy_And_Virial(
            dd.atom_numbers, dd.crd, dd.frc, md_info.pbc.cell,
            md_info.pbc.rcell, neighbor_list.full_neighbor_list.d_nl,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        listed_forces.Compute_Force(dd.atom_numbers, dd.crd, md_info.pbc.cell,
                                    md_info.pbc.rcell, dd.frc,
                                    md_info.need_potential, dd.d_energy,
                                    md_info.need_pressure, dd.d_virial);
        pairwise_force.Compute_Force(
            neighbor_list.d_nl, dd.crd, md_info.pbc.cell, md_info.pbc.rcell,
            md_info.nb.cutoff, pm.beta, dd.d_charge, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial, pm.d_direct_atom_energy);
        angle.Angle_Force_With_Atom_Energy_And_Virial(
            dd.crd, md_info.pbc.cell, md_info.pbc.rcell, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        urey_bradley.Urey_Bradley_Force_With_Atom_Energy_And_Virial(
            dd.crd, md_info.pbc.cell, md_info.pbc.rcell, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        bond.Bond_Force_With_Atom_Energy_And_Virial(
            dd.crd, md_info.pbc.cell, md_info.pbc.rcell, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        cmap.CMAP_Force_With_Atom_Energy_And_Virial(
            dd.crd, md_info.pbc.cell, md_info.pbc.rcell, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        dihedral.Dihedral_Force_With_Atom_Energy_And_Virial(
            dd.crd, md_info.pbc.cell, md_info.pbc.rcell, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        improper.Dihedral_Force_With_Atom_Energy_And_Virial(
            dd.crd, md_info.pbc.cell, md_info.pbc.rcell, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        nb14.Non_Bond_14_LJ_CF_Force_With_Atom_Energy_And_Virial(
            dd.crd, dd.d_charge, md_info.pbc.cell, md_info.pbc.rcell, dd.frc,
            md_info.need_potential, dd.d_energy, md_info.need_pressure,
            dd.d_virial);
        soft_walls.Compute_Force(dd.atom_numbers, dd.crd, dd.frc,
                                 md_info.need_potential, dd.d_energy);
        plugin.Calculate_Force(
            evaluation.commit_sampling_state, evaluation.exact_state,
            md_info.need_potential != 0, md_info.need_pressure != 0);

        restrain.Restraint(dd.crd, md_info.pbc.cell, md_info.pbc.rcell,
                           md_info.need_potential, dd.d_energy,
                           md_info.need_pressure, dd.d_virial, dd.frc, &md_info,
                           &dd);

        // With one process the PP owner also evaluates CV work locally,
        // whether PME is enabled (PM_MPI_size=1) or the NOPBC execution mode
        // has no PM owner (PM_MPI_size=0).  PM ownership only selects a
        // remote CV evaluator when the process topology is actually split.
        if (CONTROLLER::MPI_size == 1)
        {
            vatom.Coordinate_Refresh_CV(dd.crd, md_info.pbc.cell,
                                        md_info.pbc.rcell);
            if (!use_reaxff_eeq)
            {
                pm.PME_Reciprocal_Force_Async_Join(dd.d_charge, dd.frc,
                                                   md_info.need_potential,
                                                   dd.d_energy);
            }

            if (md_info.output.Check_Mdout_Step())
            {
                cv_controller.Compute_CV_For_Print(
                    cv_atom_numbers, dd.crd, md_info.pbc.cell,
                    md_info.pbc.rcell, md_info.pbc.reference_cell,
                    md_info.sys.steps);
            }

            steer_cv.Steer(cv_atom_numbers, dd.crd, md_info.pbc.cell,
                           md_info.pbc.rcell, md_info.pbc.reference_cell,
                           md_info.sys.steps, dd.d_energy, dd.d_virial, dd.frc,
                           md_info.need_potential, md_info.need_pressure);
            restrain_cv.Restraint(cv_atom_numbers, dd.crd, md_info.pbc.cell,
                                  md_info.pbc.rcell, md_info.pbc.reference_cell,
                                  md_info.sys.steps, dd.d_energy, dd.d_virial,
                                  dd.frc, md_info.need_potential,
                                  md_info.need_pressure);
            meta.Do_Metadynamics(cv_atom_numbers, dd.crd, md_info.pbc.cell,
                                 md_info.pbc.rcell, md_info.pbc.reference_cell,
                                 md_info.sys.steps, md_info.need_potential,
                                 md_info.need_pressure, dd.frc, dd.d_energy,
                                 dd.d_virial, md_info.sys.target_temperature,
                                 evaluation.commit_sampling_state);
            voronoi_detector.Observe(
                cv_atom_numbers, dd.crd, md_info.pbc.cell, md_info.pbc.rcell,
                md_info.pbc.reference_cell, md_info.sys.steps,
                evaluation.commit_sampling_state, &controller);
            vatom.Force_Redistribute_CV(dd.crd, md_info.pbc.cell,
                                        md_info.pbc.rcell, dd.frc);
        }
        else
        {
            if (!use_reaxff_eeq)
            {
                pm.Send_Recv_Force(&controller, md_info.frc, dd.frc,
                                   dd.atom_numbers);
            }
        }
    }
    else
    {
        if (!use_reaxff_eeq)
        {
            pm.reset_global_force(
                md_info.no_direct_interaction_virtual_atom_numbers);
            vatom.Coordinate_Refresh_CV(pm.g_crd, md_info.pbc.cell,
                                        md_info.pbc.rcell);
            pm.PME_Reciprocal_Force_With_Energy_And_Virial(
                md_info.crd, md_info.pbc.cell, md_info.pbc.rcell,
                md_info.d_charge, md_info.frc, md_info.need_pressure,
                md_info.need_potential, md_info.d_atom_virial_tensor,
                md_info.d_atom_energy, md_info.sys.steps,
                evaluation.exact_state);
            if (md_info.output.Check_Mdout_Step())
            {
                cv_controller.Compute_CV_For_Print(
                    cv_atom_numbers, pm.g_crd, md_info.pbc.cell,
                    md_info.pbc.rcell, md_info.pbc.reference_cell,
                    md_info.sys.steps);
            }
            steer_cv.Steer(cv_atom_numbers, pm.g_crd, md_info.pbc.cell,
                           md_info.pbc.rcell, md_info.pbc.reference_cell,
                           md_info.sys.steps, md_info.d_atom_energy,
                           md_info.d_atom_virial_tensor, pm.g_frc,
                           md_info.need_potential, md_info.need_pressure);
            restrain_cv.Restraint(
                cv_atom_numbers, pm.g_crd, md_info.pbc.cell, md_info.pbc.rcell,
                md_info.pbc.reference_cell, md_info.sys.steps,
                md_info.d_atom_energy, md_info.d_atom_virial_tensor, pm.g_frc,
                md_info.need_potential, md_info.need_pressure);
            meta.Do_Metadynamics(
                cv_atom_numbers, pm.g_crd, md_info.pbc.cell, md_info.pbc.rcell,
                md_info.pbc.reference_cell, md_info.sys.steps,
                md_info.need_potential, md_info.need_pressure, pm.g_frc,
                md_info.d_atom_energy, md_info.d_atom_virial_tensor,
                md_info.sys.target_temperature,
                evaluation.commit_sampling_state);

            vatom.Force_Redistribute_CV(pm.g_crd, md_info.pbc.cell,
                                        md_info.pbc.rcell, pm.g_frc);
            pm.add_force_g_to_l(md_info.frc);
            pm.Send_Recv_Force(&controller, md_info.frc, dd.frc,
                               dd.atom_numbers);
        }
    }
    // SITS in ALL/ITS mode enhances the complete Hamiltonian, so it must see
    // the potential assembled for this force evaluation, including energy
    // produced on dedicated PM ranks.  Reducing in Main_Iteration is too late:
    // sys.d_potential was cleared at the beginning of this routine and the
    // delayed reduction would also overwrite the bias added by SITS.
    dd.Get_Potential(&controller, &md_info);
    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        sits.Update_And_Enhance(
            md_info.sys.steps, md_info.sys.d_potential, md_info.need_pressure,
            dd.d_virial, dd.frc,
            1.0f / (CONSTANT_kB * md_info.sys.target_temperature),
            evaluation.commit_sampling_state);
        vatom.Force_Redistribute(dd.crd, md_info.pbc.cell, md_info.pbc.rcell,
                                 dd.frc);
    }
    if (sits.is_initialized && md_info.need_potential)
    {
        float effective_potential = 0.0f;
        if (CONTROLLER::MPI_rank == 0)
        {
            deviceMemcpy(&effective_potential, md_info.sys.d_potential,
                         sizeof(float), deviceMemcpyDeviceToHost);
        }
#ifdef USE_MPI
        MPI_Bcast(&effective_potential, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
        md_info.sys.h_potential = effective_potential;
        deviceMemcpy(md_info.sys.d_potential, &effective_potential,
                     sizeof(float), deviceMemcpyHostToDevice);
    }
    md_info.min.Scale_Force_For_Dynamic_Dt(
        dd.atom_numbers, dd.atom_local, dd.d_mass_inverse, dd.frc,
        dd.min_first_moment, dd.min_root_second_moment, dd.min_move);
    controller.Get_Time_Recorder("Calculate_Force")->Stop();
}

void Main_Refresh_Local_State(bool rebuild_dd)
{
    if (rebuild_dd)
    {
        dd.Send_Recv_Dom_Dec(&controller);
        dd.Find_Neighbor_Domain(&controller, &md_info);
        dd.Get_Atoms(&controller, &md_info);
    }
    dd.Get_Ghost(&controller, &md_info);
    dd.Get_Excluded(&controller, &md_info);

    if (Main_Update_Neighbor_List(neighbor_list.FORCED_UPDATE))
    {
        plugin.Set_Domain_Information(&dd);
    }

    middle_langevin.Get_Local(dd.atom_local, dd.atom_numbers);
    ad_thermo.Get_Local(dd.atom_local, dd.atom_numbers);
    nhc.Get_Local(dd.atom_local, dd.atom_numbers);

    lj.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                 dd.d_charge);
    lj_soft.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers);
    solvent_lj.Get_Local(dd.res_numbers, dd.res_len, dd.atom_numbers,
                         dd.d_mass);
    listed_forces.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                            dd.atom_local_label, dd.atom_local_id);
    pairwise_force.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                             dd.atom_local_label, dd.atom_local_id);

    angle.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                    dd.atom_local_label, dd.atom_local_id);
    urey_bradley.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                           dd.atom_local_label, dd.atom_local_id);
    bond.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                   dd.atom_local_label, dd.atom_local_id);
    cmap.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                   dd.atom_local_label, dd.atom_local_id);
    dihedral.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                       dd.atom_local_label, dd.atom_local_id);
    improper.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                       dd.atom_local_label, dd.atom_local_id);
    nb14.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                   dd.atom_local_label, dd.atom_local_id);
    restrain.Get_Local(dd.atom_local, dd.atom_numbers, dd.atom_local_label,
                       dd.atom_local_id);
    constrain.Get_Local(dd.atom_local_id, dd.atom_local_label, dd.atom_numbers);
    settle.Get_Local(dd.atom_local_id, dd.atom_local_label, dd.atom_numbers);
    vatom.Get_Local(dd.atom_local_id, dd.atom_local_label, dd.atom_numbers);
    sits.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers);
    if (sits.is_initialized && sits.selectively_applied)
    {
        sits_dihedral.Get_Local(dd.atom_local, dd.atom_numbers,
                                dd.ghost_numbers, dd.atom_local_label,
                                dd.atom_local_id);
        sits_nb14.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                            dd.atom_local_label, dd.atom_local_id);
        sits_cmap.Get_Local(dd.atom_local, dd.atom_numbers, dd.ghost_numbers,
                            dd.atom_local_label, dd.atom_local_id);
    }
    reaxff.Get_Local(&controller, dd.atom_local, dd.atom_numbers,
                     dd.ghost_numbers);
}

void Main_Iteration()
{
    controller.Get_Time_Recorder("Iteration")->Start();
    if (md_info.need_pressure || md_info.need_kinetic)
    {
        dd.Get_Ek_and_Temperature(&controller, &md_info);
    }
    if (md_info.mode != md_info.RERUN)
    {
        if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
        {
            settle.Remember_Last_Coordinates(dd.crd, md_info.pbc.cell,
                                             md_info.pbc.rcell);
            shake.Remember_Last_Coordinates(dd.crd, md_info.pbc.cell,
                                            md_info.pbc.rcell);

            if (md_info.mode == md_info.NVE)
            {
                md_info.nve.Leap_Frog(dd.atom_numbers, dd.vel, dd.crd, dd.frc,
                                      dd.d_mass_inverse, md_info.dt);
            }
            else if (md_info.mode == md_info.MINIMIZATION)
            {
                md_info.min.Gradient_Descent(dd.atom_numbers, dd.atom_local,
                                             dd.crd, dd.frc, dd.vel,
                                             dd.d_mass_inverse, dd.min_move);
                constrain.v_factor = fmaxf(FLT_MIN, md_info.min.momentum_keep);
            }
            else if (middle_langevin.is_initialized)
            {
                middle_langevin.MD_Iteration_Leap_Frog(dd.frc, dd.vel, dd.acc,
                                                       dd.crd);
                constrain.v_factor = middle_langevin.exp_gamma;
                constrain.x_factor = 0.5f * middle_langevin.exp_gamma + 0.5f;
            }
            else if (bd_thermo.is_initialized)
            {
                bd_thermo.Record_Temperature(dd.temperature,
                                             md_info.sys.freedom);
                md_info.nve.Leap_Frog(dd.atom_numbers, dd.vel, dd.crd, dd.frc,
                                      dd.d_mass_inverse, md_info.dt);
                bd_thermo.Scale_Velocity(dd.atom_numbers, dd.vel);
            }
            else if (bussi_thermo.is_initialized)
            {
                bussi_thermo.Record_Temperature(dd.temperature,
                                                md_info.sys.freedom);
                md_info.nve.Leap_Frog(dd.atom_numbers, dd.vel, dd.crd, dd.frc,
                                      dd.d_mass_inverse, md_info.dt);
                bussi_thermo.Scale_Velocity(dd.atom_numbers, dd.vel);
            }
            else if (ad_thermo.is_initialized)
            {
                if ((md_info.sys.steps - 1) % ad_thermo.update_interval == 0)
                {
                    ad_thermo.MD_Iteration_Leap_Frog(dd.vel, dd.crd, dd.frc,
                                                     dd.acc, md_info.dt);
                    settle.Project_Velocity_To_Constraint_Manifold(
                        dd.vel, dd.crd, dd.d_mass_inverse, md_info.pbc.cell,
                        md_info.pbc.rcell);
                    shake.Project_Velocity_To_Constraint_Manifold(
                        dd.vel, dd.crd, dd.d_mass_inverse, md_info.pbc.cell,
                        md_info.pbc.rcell, dd.atom_numbers);
                    constrain.v_factor = FLT_MIN;
                    constrain.x_factor = 0.5;
                }
                else
                {
                    md_info.nve.Leap_Frog(dd.atom_numbers, dd.vel, dd.crd,
                                          dd.frc, dd.d_mass_inverse,
                                          md_info.dt);
                    constrain.v_factor = 1.0;
                    constrain.x_factor = 1.0;
                }
            }
            else if (nhc.is_initialized)
            {
                nhc.MD_Iteration_Leap_Frog(dd.vel, dd.crd, dd.frc, dd.acc,
                                           md_info.dt, dd.h_ek_total,
                                           md_info.sys.freedom);
            }

            settle.Do_SETTLE(&controller, dd.atom_local, dd.d_mass, dd.crd,
                             md_info.pbc.cell, md_info.pbc.rcell, dd.vel,
                             md_info.need_pressure, md_info.sys.d_stress);
            shake.Constrain(dd.atom_numbers, dd.crd, dd.vel, dd.d_mass_inverse,
                            dd.d_mass, md_info.pbc.cell, md_info.pbc.rcell,
                            md_info.need_pressure, md_info.sys.d_stress);
            hard_wall.Reflect(dd.atom_numbers, dd.crd, dd.vel);
        }
        if (md_info.need_pressure && !mc_baro.is_initialized)
        {
            md_info.Get_pressure(&controller, dd.atom_numbers, dd.vel,
                                 dd.d_mass, dd.d_virial, main_stream);
            md_info.sys.Get_Density();
            press_baro.Regulate_Pressure(
                md_info.sys.steps, md_info.sys.h_stress, md_info.pbc.cell,
                md_info.dt, md_info.sys.target_pressure,
                md_info.sys.target_temperature);
        }
    }
    else
    {
        const bool rerun_box_changed = md_info.rerun.Iteration();
        if (md_info.rerun.need_box_update && rerun_box_changed)
        {
            Main_Rerun_Box_Change(md_info.rerun.g,
                                  md_info.rerun.frame_box_length,
                                  md_info.rerun.frame_box_angle);
        }
        md_info.Crd_Vel_Device_to_dd(dd.crd, dd.vel, dd.atom_local_label,
                                     dd.atom_local_id, main_stream);
    }

    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        vatom.Coordinate_Refresh(dd.crd, md_info.pbc.cell, md_info.pbc.rcell);
        if (Next_Step_Is_Interval_Boundary(md_info.sys.steps,
                                           dd.update_interval) ||
            md_info.mode == md_info.RERUN)
        {
            if (CONTROLLER::PP_MPI_size != 1)
            {
                controller.Get_Time_Recorder("Communication")->Start();
                dd.Exchange_Particles(&controller, &md_info);
                controller.Get_Time_Recorder("Communication")->Stop();
                Main_Refresh_Local_State(false);
            }
            else
            {
                if (Main_Update_Neighbor_List(neighbor_list.FORCED_UPDATE))
                {
                    plugin.Set_Domain_Information(&dd);
                }
            }
        }
    }
    if (Next_Step_Is_Interval_Boundary(md_info.sys.steps, dd.update_interval) ||
        md_info.mode == md_info.RERUN)
    {
        controller.Get_Time_Recorder("Communication")->Start();
        pm.Get_Atoms(&controller, md_info.crd, md_info.d_charge,
                     dd.atom_numbers, dd.crd, dd.d_charge, dd.atom_local, true,
                     true, true, true);
        controller.Get_Time_Recorder("Communication")->Stop();
    }
    const bool coordinate_state_may_have_changed =
        md_info.mode == md_info.RERUN || md_info.mode == md_info.MINIMIZATION ||
        md_info.dt != 0.0f || constrain.is_initialized ||
        settle.is_initialized || shake.is_initialized ||
        hard_wall.is_initialized || vatom.is_initialized ||
        press_baro.is_initialized;
    if (coordinate_state_may_have_changed)
        Main_Advance_Coordinate_Generation("advancing an MD coordinate state");
    controller.Get_Time_Recorder("Iteration")->Stop();
}

void Main_Print()
{
    if (md_info.output.Check_Mdout_Step())
    {
        md_info.Step_Print(&controller);
        if (!md_info.pbc.pbc)
        {
            CF_NOPBC.Step_Print(&controller);
            LJ_NOPBC.Step_Print(&controller);
            gb.Step_Print(&controller);
        }
        else
        {
            lj.Step_Print(&controller);
            lj_soft.Step_Print(&controller);
            pm.Step_Print(&controller);
            sits.Step_Print(&controller, 1.0f / md_info.sys.target_temperature /
                                             CONSTANT_kB);
        }
        sits_dihedral.Step_Print(&controller, false);
        sits_nb14.Step_Print(&controller, false);
        sits_cmap.Step_Print(&controller, false);

        sw.Step_Print(&controller);
        eam.Step_Print(&controller);
        tersoff.Step_Print(&controller);
        reaxff.Step_Print(&controller, md_info.d_charge);
        pairwise_force.Step_Print(&controller);
        angle.Step_Print(&controller);
        urey_bradley.Step_Print(&controller);
        bond.Step_Print(&controller);
        cmap.Step_Print(&controller);
        listed_forces.Step_Print(&controller);
        dihedral.Step_Print(&controller);
        improper.Step_Print(&controller);
        nb14.Step_Print(&controller);

        controller.Step_Print("potential", dd.h_sum_ene_total);

        restrain.Step_Print(&controller);
        if (qc.is_initialized)
        {
            qc.Step_Print(&controller);
        }
        cv_controller.Step_Print();
        plugin.Mdout_Print();
        steer_cv.Step_Print(&controller);
        restrain_cv.Step_Print(&controller);
        meta.Step_Print(&controller);
        voronoi_detector.Step_Print(&controller);
        soft_walls.Step_Print(&controller);
        controller.Print_To_Screen_And_Mdout();
    }

    if (md_info.output.Check_Trajectory_Step())
    {
        md_info.Crd_Vel_dd_to_Device(dd.crd, dd.vel, dd.atom_local_label,
                                     dd.atom_local_id, main_stream);
        if (md_info.pbc.pbc)
        {
            md_info.mol.Molecule_Crd_Map();
            md_info.Crd_Vel_Device_to_dd(dd.crd, dd.vel, dd.atom_local_label,
                                         dd.atom_local_id, main_stream);
        }
        md_info.output.Append_Crd_Traj_File();
        md_info.output.Append_Vel_Traj_File();
        md_info.output.Append_Box_Traj_File();
        meta.Write_Potential();
        nhc.Save_Trajectory_File();
    }

    if (md_info.output.is_frc_traj && md_info.output.Check_Force_Step())
    {
        md_info.Frc_dd_to_Host(dd.frc, dd.atom_local_label, dd.atom_local_id,
                               main_stream);
        md_info.output.Append_Frc_Traj_File();
    }

    if (md_info.output.Check_Restart_Step())
    {
        md_info.output.Export_Restart_File();
        nhc.Save_Restart_File();
    }
}

void Main_Clear()
{
    md_info.rerun.Clear();
    md_info.min.Clear();
    lj_soft.Clear();
    LJ_NOPBC.Clear();
    CF_NOPBC.Clear();
    gb.Clear();
    controller.Final_Time_Summary(
        md_info.sys.steps, md_info.sys.speed_time_factor,
        md_info.sys.speed_unit_name.c_str(), md_info.mode);

    controller.Clear();
}

float Main_Box_Change(LTMatrix3 g, int scale_box, int scale_crd, int scale_vel)
{
    return Main_Box_Change_Transactional(g, scale_box, scale_crd, scale_vel,
                                         true);
}

float Main_Box_Change_Transactional(LTMatrix3 g, int scale_box, int scale_crd,
                                    int scale_vel, bool commit_box_state,
                                    const VECTOR* authoritative_box_length,
                                    const VECTOR* authoritative_box_angle)
{
    if ((authoritative_box_length == NULL) != (authoritative_box_angle == NULL))
    {
        controller.Throw_SPONGE_Error(
            spongeErrorSimulationBreakDown, "Main_Box_Change_Transactional",
            "Reason:\n\tauthoritative box length and angle must be supplied "
            "together\n");
    }
    if (scale_box)
    {
        if (authoritative_box_length != NULL)
        {
            md_info.pbc.Update_Box_From_Input(*authoritative_box_length,
                                              *authoritative_box_angle);
        }
        else
        {
            md_info.pbc.Update_Box(g);
        }
    }
    // 放缩坐标与速度
    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        md_info.Scale_Positions_And_Velocities(
            g, scale_crd, scale_vel, dd.crd,
            dd.vel);  // rescale dd进程原子坐标与速度
        if (commit_box_state)
        {
            restrain.Update_Refcoord_Scaling(
                &md_info, g, md_info.dt, dd.atom_local, dd.atom_numbers,
                dd.atom_local_label, dd.atom_local_id);
        }
    }

    // 大幅度放缩盒子时，重新初始化相关模块
    if (scale_box && commit_box_state && md_info.pbc.Check_Change_Large())
    {
        Main_Box_Change_Largely();
    }
    else  // 更新域分解盒子
    {
        if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
        {
            dd.Update_Box(g, md_info.dt);
        }
        if (CONTROLLER::PM_MPI_rank < CONTROLLER::PM_MPI_size &&
            CONTROLLER::PM_MPI_rank != -1)
        {
            pm.Update_Box(md_info.pbc.cell, md_info.pbc.rcell, g, md_info.dt);
        }
    }
    return md_info.sys.Get_Volume();
}

float Main_Rerun_Box_Change(LTMatrix3 g, VECTOR box_length, VECTOR box_angle)
{
    return Main_Box_Change_Transactional(g, 1, 0, 0, true, &box_length,
                                         &box_angle);
}

void Main_Box_Change_Largely()
{
    controller.printf(
        "Some modules are based on the meshing methods, and it is more "
        "precise "
        "to re-initialize these modules now for a large box change.\n");

    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        md_info.Crd_Vel_dd_to_Device(dd.crd, dd.vel, dd.atom_local_label,
                                     dd.atom_local_id, main_stream);
    }
    neighbor_list.Clear();
    neighbor_list.Initial(&controller, md_info.atom_numbers, md_info.nb.cutoff,
                          md_info.nb.skin, md_info.pbc.cell, md_info.pbc.rcell);
    pm.Clear();
    pm.Initial(&controller, md_info.atom_numbers, md_info.pbc.cell,
               md_info.pbc.rcell, md_info.sys.box_length, md_info.nb.cutoff,
               md_info.no_direct_interaction_virtual_atom_numbers);
    dd.Free_Buffer();
    dd.Domain_Decomposition(&controller, &md_info);
    pm.Domain_Decomposition(&controller, md_info.sys.box_length,
                            dd.dom_dec_split_num);
    pm.Send_Recv_Dom_Dec(&controller);
    pm.Find_Neighbor_Domain(&controller);
    if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
    {
        Main_Refresh_Local_State(true);
        plugin.Set_Domain_Information(&dd);
    }
    pm.Get_Atoms(&controller, md_info.crd, md_info.d_charge, dd.atom_numbers,
                 dd.crd, dd.d_charge, dd.atom_local, true, true, true, true);
    MPI_Barrier(MPI_COMM_WORLD);
    controller.printf(
        "------------------------------------------------------------------"
        "----"
        "--------------------------------------\n");
}

void Main_Process_Management()
{
    CONTROLLER::PM_MPI_size = pm.PM_MPI_size;
    CONTROLLER::PP_MPI_size =
        (CONTROLLER::MPI_size - CONTROLLER::PM_MPI_size -
         CONTROLLER::CC_MPI_size) <= 0
            ? 1
            : (CONTROLLER::MPI_size - CONTROLLER::PM_MPI_size -
               CONTROLLER::CC_MPI_size);

    if (CONTROLLER::MPI_size == 1)
    {
        CONTROLLER::pp_comm = MPI_COMM_WORLD;
        CONTROLLER::pm_comm = MPI_COMM_WORLD;
        CONTROLLER::PP_MPI_rank = 0;
        dd.pp_rank = 0;
        if (CONTROLLER::PM_MPI_size != 0)
        {
            CONTROLLER::PM_MPI_rank = 0;
            pm.pm_rank = 0;
        }
        else
        {
            CONTROLLER::PM_MPI_rank = -1;
            pm.pm_rank = -1;
        }
    }
    else if (CONTROLLER::PM_MPI_size == 0)
    {
        CONTROLLER::pp_comm = MPI_COMM_WORLD;
        CONTROLLER::PP_MPI_rank = CONTROLLER::MPI_rank;
        dd.pp_rank = CONTROLLER::PP_MPI_rank;
        pm.pm_rank = -1;
#ifdef USE_XCCL
        xcclUniqueId pp_id;
        if (CONTROLLER::PP_MPI_rank == 0)
        {
            xcclGetUniqueId(&pp_id);
        }
        MPI_Bcast(&pp_id, sizeof(pp_id), MPI_BYTE, 0, CONTROLLER::pp_comm);
        xcclCommInitRank(&CONTROLLER::d_pp_comm, CONTROLLER::PP_MPI_size, pp_id,
                         CONTROLLER::PP_MPI_rank);
#else
        CONTROLLER::d_pp_comm = CONTROLLER::pp_comm;
#endif
    }
    else
    {
        if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
        {
            MPI_Comm_split(MPI_COMM_WORLD, 0, CONTROLLER::MPI_rank,
                           &CONTROLLER::pp_comm);
            MPI_Comm_rank(CONTROLLER::pp_comm, &dd.pp_rank);
            CONTROLLER::PP_MPI_rank = dd.pp_rank;
#ifdef USE_XCCL
            xcclUniqueId pp_id;
            if (CONTROLLER::PP_MPI_rank == 0)
            {
                xcclGetUniqueId(&pp_id);
            }
            MPI_Bcast(&pp_id, sizeof(pp_id), MPI_BYTE, 0, CONTROLLER::pp_comm);
            xcclCommInitRank(&CONTROLLER::d_pp_comm, CONTROLLER::PP_MPI_size,
                             pp_id, CONTROLLER::PP_MPI_rank);
#else
            CONTROLLER::d_pp_comm = CONTROLLER::pp_comm;
#endif
        }
        else
        {
            CONTROLLER::PP_MPI_rank =
                CONTROLLER::PP_MPI_size;  // PP_MPI_rank 设置>=
                                          // PP_MPI_size，表示非PP进程
            MPI_Comm_split(MPI_COMM_WORLD, 1, CONTROLLER::MPI_rank,
                           &CONTROLLER::pm_comm);
            MPI_Comm_rank(CONTROLLER::pm_comm, &pm.pm_rank);
            CONTROLLER::PM_MPI_rank = pm.pm_rank;
#ifdef USE_XCCL
            xcclUniqueId pm_id;
            if (CONTROLLER::PM_MPI_rank == 0)
            {
                xcclGetUniqueId(&pm_id);
            }
            MPI_Bcast(&pm_id, sizeof(pm_id), MPI_BYTE, 0, CONTROLLER::pm_comm);
            xcclCommInitRank(&CONTROLLER::d_pm_comm, CONTROLLER::PM_MPI_size,
                             pm_id, CONTROLLER::PM_MPI_rank);
#else
            CONTROLLER::d_pm_comm = CONTROLLER::pm_comm;
#endif
        }
    }

    controller.printf(
        "MPI process total: MPI_size=%d, PP_MPI_size=%d, PM_MPI_size=%d\n",
        CONTROLLER::MPI_size, CONTROLLER::PP_MPI_size, CONTROLLER::PM_MPI_size);
    controller.MPI_printf(
        "MPI process partition: MPI_rank=%d, PP_MPI_rank=%d, "
        "PM_MPI_rank=%d\n",
        CONTROLLER::MPI_rank, CONTROLLER::PP_MPI_rank, CONTROLLER::PM_MPI_rank);

    if (CONTROLLER::PP_MPI_size > 1)
    {
        md_info.nb.Excluded_List_Reform(&controller, md_info.atom_numbers);
    }
    pm.exclude_factor = CONTROLLER::PP_MPI_size == 1 ? 1.0f : 0.5f;

    deviceStreamCreate(&main_stream);
    dd.Create_Stream();
    pm.Create_Stream();

    dd.Domain_Decomposition(&controller, &md_info);
    pm.Domain_Decomposition(&controller, md_info.sys.box_length,
                            dd.dom_dec_split_num);
    pm.Send_Recv_Dom_Dec(&controller);
    pm.Find_Neighbor_Domain(&controller);
}

void Main_MC_Barostat()
{
    if (mc_baro.Will_Attempt(md_info.sys.steps))
    {
        const FORCE_EVALUATION_CONTEXT transactional_evaluation(false, true);

        // Establish the old Hamiltonian from the current coordinates and the
        // already-committed adaptive history.  MC runs before this step's sole
        // committed force evaluation, so old and trial necessarily see the
        // same SITS statistics and metadynamics hills.
        Main_Calculate_Force(transactional_evaluation);

        SITS_STATE_SNAPSHOT sits_state_old;
        mc_baro.energy_old = md_info.sys.h_potential;
        const LTMatrix3 accepted_cell = md_info.pbc.cell;
        const LTMatrix3 accepted_rcell = md_info.pbc.rcell;
        const LTMatrix3 accepted_reference_cell = md_info.pbc.reference_cell;
        const LTMatrix3 accepted_cell0 = md_info.pbc.cell0;
        const VECTOR accepted_box_length = md_info.sys.box_length;
        const VECTOR accepted_box_angle = md_info.sys.box_angle;
        const std::uint64_t accepted_coordinate_generation =
            current_coordinate_generation;
        if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
        {
            sits.Save_State(&sits_state_old, md_info.sys.d_potential);
            restrain.Save_Refcoord_Transaction_State();
            // The backup is in global atom order.  A local-layout backup is
            // invalid if an accepted large box change later remeshes DD.
            md_info.Crd_Vel_dd_to_Device(dd.crd, dd.vel, dd.atom_local_label,
                                         dd.atom_local_id, main_stream);
            deviceMemcpy(mc_baro.crd_backup, md_info.crd,
                         sizeof(VECTOR) * md_info.atom_numbers,
                         deviceMemcpyDeviceToDevice);
        }
        mc_baro.Volume_Change_Attempt(md_info.sys.box_length, md_info.dt);
        // Suppress the ordinary one-way box side effects here.  The trial's
        // reference state and full remesh are applied explicitly below under
        // transaction control.
        Main_Box_Change_Transactional(mc_baro.g, 1, 0, 0, false);
        if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
        {
            dd.Res_Crd_Map(mc_baro.g, md_info.dt);
            // Reference scaling is part of the trial Hamiltonian.  Its global
            // state was snapshotted above, so rejection restores it exactly.
            restrain.Update_Refcoord_Scaling(
                &md_info, mc_baro.g, md_info.dt, dd.atom_local, dd.atom_numbers,
                dd.atom_local_label, dd.atom_local_id);
        }

        // A trial must use the same DD halo membership and PME discretization
        // that it would have after acceptance.  Rebuild transactionally now;
        // the global-order coordinate backup makes rejection safe even if the
        // local layout changes completely.
        Main_Box_Change_Largely();
        md_info.pbc.cell0 = md_info.pbc.cell;

        if (!Main_Cell_Is_Exactly_Equal(md_info.pbc.cell, accepted_cell))
            Main_Advance_Coordinate_Generation(
                "constructing an MC barostat trial state");

        Main_Calculate_Force(transactional_evaluation);
        mc_baro.energy_new = md_info.sys.h_potential;
        if (CONTROLLER::MPI_rank == 0)
        {
            mc_baro.extra_term =
                static_cast<double>(md_info.sys.target_pressure) *
                    mc_baro.DeltaV -
                static_cast<double>(md_info.ug.ug_numbers) * CONSTANT_kB *
                    md_info.sys.target_temperature * log(mc_baro.VDevided);
            if (mc_baro.couple_dimension != mc_baro.NO &&
                mc_baro.couple_dimension != mc_baro.XYZ)
            {
                mc_baro.extra_term -=
                    static_cast<double>(mc_baro.surface_number) *
                    mc_baro.surface_tension * mc_baro.DeltaS;
            }
            const double acceptance_delta =
                static_cast<double>(mc_baro.energy_new) -
                static_cast<double>(mc_baro.energy_old) + mc_baro.extra_term;
            const double thermal_energy = static_cast<double>(CONSTANT_kB) *
                                          md_info.sys.target_temperature;
            if (!isfinite(mc_baro.energy_old) ||
                !isfinite(mc_baro.energy_new) ||
                !isfinite(mc_baro.extra_term) || !isfinite(acceptance_delta) ||
                !isfinite(thermal_energy) || !(thermal_energy > 0.0))
            {
                controller.Throw_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "Main_MC_Barostat",
                    "Reason:\n\tthe MC acceptance Hamiltonian is non-finite "
                    "or has a non-positive thermal energy\n");
            }
            mc_baro.accept_possibility =
                acceptance_delta <= 0.0
                    ? 1.0
                    : exp(-acceptance_delta / thermal_energy);
            if (!isfinite(mc_baro.accept_possibility) ||
                mc_baro.accept_possibility < 0.0 ||
                mc_baro.accept_possibility > 1.0)
            {
                controller.Throw_SPONGE_Error(
                    spongeErrorSimulationBreakDown, "Main_MC_Barostat",
                    "Reason:\n\tthe MC acceptance probability is invalid\n");
            }
        }

        const bool mc_accepted = mc_baro.Check_MC_Barostat_Accept() != 0;
        if (!mc_accepted)  // 如果不接受
        {
            const LTMatrix3 reverse_g = mc_baro.Get_Exact_Reverse_G(md_info.dt);
            if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
            {
                restrain.Restore_Refcoord_Transaction_State(
                    dd.atom_local, dd.atom_numbers, dd.atom_local_label,
                    dd.atom_local_id);
                // Restore accepted coordinates into the current (possibly
                // remeshed) trial layout before rebuilding the old box.
                deviceMemcpy(md_info.crd, mc_baro.crd_backup,
                             sizeof(VECTOR) * md_info.atom_numbers,
                             deviceMemcpyDeviceToDevice);
                md_info.Crd_Vel_Device_to_dd(dd.crd, dd.vel,
                                             dd.atom_local_label,
                                             dd.atom_local_id, main_stream);
            }
            Main_Box_Change_Transactional(reverse_g, 1, 0, 0, false);

            // Remove round-trip floating-point drift from the authoritative
            // box, then rebuild every local/mesh representation from it.
            md_info.pbc.cell = accepted_cell;
            md_info.pbc.rcell = accepted_rcell;
            md_info.pbc.reference_cell = accepted_reference_cell;
            md_info.pbc.cell0 = accepted_cell0;
            md_info.sys.box_length = accepted_box_length;
            md_info.sys.box_angle = accepted_box_angle;
            current_coordinate_generation = accepted_coordinate_generation;
            Main_Box_Change_Largely();
            if (CONTROLLER::MPI_rank < CONTROLLER::PP_MPI_size)
            {
                sits.Restore_State(sits_state_old, md_info.sys.d_potential);
            }
        }
        mc_baro.Delta_Box_Length_Max_Update();
    }
}

void Main_Sync_Dynamic_Targets_To_Controllers()
{
    md_info.sys.Update_Targets_By_Schedule(&controller, md_info.sys.steps);
    const float target_temperature = md_info.sys.target_temperature;
    bd_thermo.Set_Target_Temperature(target_temperature);
    bussi_thermo.Set_Target_Temperature(target_temperature);
    ad_thermo.Set_Target_Temperature(target_temperature);
    middle_langevin.Set_Target_Temperature(target_temperature);
    nhc.Set_Target_Temperature(target_temperature);
}
