import os
from pathlib import Path

import numpy as np
import pytest

from benchmarks.utils import Outputer, Runner
from benchmarks.validation.utils import parse_mdout_column

MDIN = """md_name = "validation tip3p neighbor overflow exact recovery"
mode = "nvt"
step_limit = 50
dt = 0.001
cutoff = 8.0
skin = 1.0
default_in_file_prefix = "tip3p"
constrain_mode = "SETTLE"
print_zeroth_frame = 1
write_mdout_interval = 1
write_information_interval = 10
write_restart_file_interval = 0
thermostat = "middle_langevin"
thermostat_tau = 0.01
thermostat_seed = 2026
target_temperature = 300.0
neighbor_list.max_atom_in_grid_numbers = 2
neighbor_list.max_ghost_in_grid_numbers = 2
neighbor_list.max_neighbor_numbers = 32
"""


def test_neighbor_overflow_exact_recovery_keeps_md_stable(
    statics_path, outputs_path, mpi_np
):
    # Regression coverage for NEIGHBOR_LIST::Update_With_Overflow_Recovery:
    # tiny grid/neighbor capacities force exact-size rebuilds during the run,
    # which is the path that exercised the odd-capacity GPU shared-memory
    # alignment bug (fixed in 17c7448).
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p",
        mpi_np=mpi_np,
        run_name="tip3p_overflow_exact_recovery",
    )
    Path(case_dir, "mdin.spg.toml").write_text(MDIN)

    output = Runner.run_sponge(
        case_dir,
        timeout=600,
        mpi_np=mpi_np,
        sponge_cmd=os.environ.get("SPONGE_BIN"),
    )

    recovery_count = output.count("rebuilding exactly")
    assert recovery_count > 0, (
        "neighbor-list exact overflow recovery was not triggered\n" + output
    )

    potentials = np.asarray(
        parse_mdout_column(case_dir / "mdout.txt", "potential"),
        dtype=float,
    )
    assert potentials.size > 1
    Outputer.print_table(
        ["Metric", "Value"],
        [
            ["Case", "tip3p_overflow_exact_recovery"],
            ["Recoveries", recovery_count],
            ["Frames", potentials.size],
            ["PotentialMean", f"{np.mean(potentials):.4f}"],
            ["PotentialMin", f"{np.min(potentials):.4f}"],
            ["PotentialMax", f"{np.max(potentials):.4f}"],
            ["Status", "PASS"],
        ],
        title="Misc Validation: Neighbor Overflow Exact Recovery",
    )
    assert np.all(np.isfinite(potentials)), potentials
    median = float(np.median(potentials))
    assert median < 0.0
    spread = float(np.max(np.abs(potentials - median))) / abs(median)
    assert spread < 0.05, f"potential energy exploded: {potentials}"
