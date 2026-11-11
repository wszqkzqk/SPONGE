import os

import pytest

from benchmarks.utils import Outputer, Runner
from benchmarks.validation.cv.tests.utils import (
    DISTANCE_ATOMS,
    load_box_lengths,
    load_coordinates,
    write_validation_mdin,
)
from benchmarks.validation.utils import parse_mdout_column


def test_distance_path_uses_continuous_coordinates(
    statics_path, outputs_path, mpi_np
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="alanine_dipeptide_phi_psi",
        mpi_np=mpi_np,
        run_name="distance_path",
    )
    coordinate_path = case_dir / "sys_flexible_coordinate.txt"
    coordinates = load_coordinates(coordinate_path)
    box = load_box_lengths(coordinate_path)
    coordinates[DISTANCE_ATOMS[0]] = (0.1 * box[0], 7.0, 7.0)
    coordinates[DISTANCE_ATOMS[1]] = (0.9 * box[0], 9.0, 9.0)
    coordinates[5] = (0.5 * box[0], 20.0, 15.0)
    tokens = coordinate_path.read_text(encoding="utf-8").split()
    box_tokens = tokens[1 + len(coordinates) * 3 :]
    coordinate_path.write_text(
        f"{len(coordinates)}\n"
        + "\n".join(
            " ".join(f"{value:.9f}" for value in point) for point in coordinates
        )
        + "\n"
        + " ".join(box_tokens)
        + "\n",
        encoding="utf-8",
    )
    (case_dir / "cv.txt").write_text(
        "print\n{\n    CV = distance distance_path\n}\n"
        "distance\n{\n    CV_type = distance\n    atom = 4 6\n}\n"
        "distance_path\n{\n    CV_type = distance_path\n    path_atom = 4 5 6\n}\n",
        encoding="utf-8",
    )
    write_validation_mdin(case_dir, step_limit=1)
    Runner.run_sponge(
        case_dir,
        timeout=1200,
        mpi_np=mpi_np,
        sponge_cmd=os.environ.get("SPONGE_BIN", "SPONGE"),
    )

    unwrapped = ((0.8 * box[0]) ** 2 + 2.0**2 + 2.0**2) ** 0.5
    minimum_image = ((0.2 * box[0]) ** 2 + 2.0**2 + 2.0**2) ** 0.5
    actual_unwrapped = float(
        parse_mdout_column(case_dir / "mdout.txt", "distance_path")[-1]
    )
    actual_minimum_image = float(
        parse_mdout_column(case_dir / "mdout.txt", "distance")[-1]
    )
    assert actual_unwrapped == pytest.approx(unwrapped, abs=5e-4)
    assert actual_minimum_image == pytest.approx(minimum_image, abs=5e-4)
