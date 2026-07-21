import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

PROBE_SOURCE = Path(__file__).with_name("plugin_api_negotiation_probe.cpp")


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = (
        shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the plugin ABI probe")
    return [compiler]


def _sponge_binary():
    configured = os.environ.get("SPONGE_BIN", "SPONGE")
    path = shutil.which(configured)
    if path is None and Path(configured).is_file():
        path = str(Path(configured).resolve())
    if path is None:
        pytest.skip("SPONGE executable is required for the plugin ABI probe")
    return path


def _compile_probe(tmp_path, version):
    library = tmp_path / f"plugin_api_v{version}.so"
    result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-shared",
            "-fPIC",
            f"-DPROBE_API_VERSION={version}",
            str(PROBE_SOURCE),
            "-o",
            str(library),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, (
        f"failed to compile API v{version} probe\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    return library


def _write_case(case_dir, plugins):
    (case_dir / "mass.txt").write_text("2\n12\n12\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("2\n0\n0\n", encoding="utf-8")
    (case_dir / "coordinate.txt").write_text(
        "2\n5 10 10\n12 10 10\n20 20 20\n90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "velocity.txt").write_text(
        "2\n0 0 0\n0 0 0\n", encoding="utf-8"
    )
    (case_dir / "lj.txt").write_text(
        "2 1\n0\n0\n0\n0\n", encoding="utf-8"
    )
    settings = {
        "md_name": "stable plugin API negotiation probe",
        "mode": "nvt",
        "step_limit": 0,
        "dt": 0.0,
        "cutoff": 4.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "velocity_in_file": "velocity.txt",
        "LJ_in_file": "lj.txt",
        "thermostat": "berendsen_thermostat",
        "thermostat_tau": 1.0,
        "target_temperature": 300.0,
        "plugin": " ".join(path.as_posix() for path in plugins),
        "mdout": "mdout.txt",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )


def _read_fields(path):
    return dict(field.split("=", 1) for field in path.read_text().split())


@pytest.mark.skipif(sys.platform == "win32", reason="probe builds POSIX .so files")
def test_stable_plugin_api_negotiates_append_only_versions(tmp_path):
    # Exercise both the complete plugin-list value and each individual path
    # beyond the historical CHAR_LENGTH_MAX stack buffer.  Shortening these
    # paths would hide command-parser and loader truncation/overwrite bugs.
    plugin_dir = tmp_path.joinpath(
        *(
            f"plugin_path_segment_{index}{'_hash#' if index == 2 else '_'}"
            + "x" * 110
            for index in range(4)
        )
    )
    plugin_dir.mkdir(parents=True)
    plugins = [
        _compile_probe(plugin_dir, version) for version in (2, 3, 4, 5)
    ]
    plugins.append(_compile_probe(plugin_dir, 0))
    assert all(len(path.as_posix()) > 512 for path in plugins)
    assert len(" ".join(path.as_posix() for path in plugins)) > 5 * 512
    _write_case(tmp_path, plugins)

    result = subprocess.run(
        [_sponge_binary(), "-mdin", "mdin.spg.toml"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output

    for version in (2, 3, 4, 5):
        fields = _read_fields(tmp_path / f"stable_api_v{version}.log")
        assert fields["initial"] == str(version)
        assert fields["after"] == str(version)
        assert fields["callbacks"] == "1"
        assert fields["copy"] == "1"
        if version == 5:
            if fields["type"] == "1":
                assert fields["device"] == "0"
            else:
                assert int(fields["device"]) >= 0

    legacy = _read_fields(tmp_path / "legacy_api.log")
    assert legacy["initialized"] == "1"
    assert int(legacy["argument"]) not in {2, 3, 4, 5}
