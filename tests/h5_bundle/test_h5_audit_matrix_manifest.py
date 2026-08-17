#!/usr/bin/env python3
import argparse
import re
from pathlib import Path


def fail(message):
    raise AssertionError(message)


def normalize(text):
    return re.sub(r"\s+", " ", text)


def require_tokens(text, label, tokens):
    missing = [token for token in tokens if token not in text]
    if missing:
        fail(f"{label} missing required audit tokens: {missing}")


def parse_plan_runtime_exceptions(plan_text):
    exceptions = []
    in_table = False
    for line in plan_text.splitlines():
        stripped = line.strip()
        if stripped == "Current runtime native-parity exceptions:":
            in_table = True
            continue
        if not in_table:
            continue
        if not stripped:
            if exceptions:
                break
            continue
        if not stripped.startswith("|"):
            continue
        if stripped.startswith("|---") or stripped.startswith("| Exception "):
            continue
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if cells:
            exception = cells[0].strip("`")
            if exception:
                exceptions.append(exception)
    if not exceptions:
        fail("failed to parse runtime native-parity exceptions from plan")
    return exceptions


def parse_plan_pure_bundled_core_columns(plan_text):
    match = re.search(
        r"Pure bundled rerun native core mdout columns:\s+```text\n"
        r"(?P<body>.*?)```",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse pure bundled rerun native core mdout columns")
    columns = [
        line.strip()
        for line in match.group("body").splitlines()
        if line.strip()
    ]
    if not columns:
        fail("pure bundled rerun native core mdout column list is empty")
    return columns


def parse_plan_phases(plan_text):
    phases = []
    for match in re.finditer(
        r"^## Phase\s+(\d+):\s+(.+)$", plan_text, re.MULTILINE
    ):
        phases.append(f"Phase {match.group(1)}")
    if not phases:
        fail("failed to parse phase headings from plan")
    return phases


def parse_implementation_file_map(plan_text):
    match = re.search(
        r"## Implementation File Map\s+"
        r".*?"
        r"\| Phase \| Primary files \| Verification boundary \|\n"
        r"\|---\|---\|---\|\n"
        r"(?P<body>.*?)(?:\n\n|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Implementation File Map from plan")

    rows = []
    for raw_line in match.group("body").splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue
        if not stripped.startswith("|"):
            break
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if len(cells) != 3:
            fail(f"Implementation File Map row must have 3 cells: {stripped}")
        rows.append(tuple(cells))
    if not rows:
        fail("Implementation File Map contains no rows")
    return rows


def parse_input_fixture_generation(plan_text):
    match = re.search(
        r"## Input Fixture Generation\s+"
        r"(?P<body>.*?)(?:\n## Fixture Roots|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Input Fixture Generation from plan")
    body = match.group("body")
    required_tokens = [
        "legacy_input/",
        "bundled_input/",
        "bundled_input_with_legacy_sidecar/",
        "manifest.json",
        "core_structural",
        "full_contract_rerun",
        "/parameters/sponge/files/legacy_sidecars",
        "legacy_sidecars/<key>/<basename>",
        "byte-equivalent",
        "byte-match",
        "relocation-safe",
        "Typed native datasets",
        "qc_type.txt",
        "/qc/type",
    ]
    require_tokens(body, "Input Fixture Generation plan", required_tokens)
    return body


def parse_fixture_roots(plan_text):
    match = re.search(
        r"## Fixture Roots\s+"
        r"(?P<body>.*?)(?:\n## Phase 0: Fixture Helper|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Fixture Roots from plan")
    body = match.group("body")
    required_tokens = [
        "tests/h5_bundle/fixtures/input_matrix",
        "Primary fixture groups",
        "core_structural",
        "Minimal normal-mode case",
        "legacy_input",
        "bundled_input",
        "bundled_input_with_legacy_sidecar",
        "first contract tests",
        "first smoke matrix",
        "full_contract_rerun",
        "Broad rerun-mode case",
        "XPONGE's full converter contract",
        "rerun smoke",
        "later contract coverage",
    ]
    require_tokens(body, "Fixture Roots plan", required_tokens)
    return body


def parse_phase0_fixture_helper(plan_text):
    match = re.search(
        r"## Phase 0: Fixture Helper\s+"
        r"(?P<body>.*?)(?:\n## Phase 1: Input Contract Tests|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase 0 Fixture Helper from plan")
    body = match.group("body")
    required_tokens = [
        "tests/h5_bundle/h5_input_matrix_fixture.hpp",
        "Resolve the repository root",
        "tests/h5_bundle/fixtures/input_matrix",
        "Copy a selected fixture case",
        "per-test temporary working directory",
        "Return paths to mdin files",
        "bundle directories",
        "expected output",
        "Normalize path handling",
        "CTest",
        "direct binaries",
        "local scripts",
        "should not run SPONGE",
        "stable test plumbing",
        "copy `core_structural` and `full_contract_rerun`",
        "isolated temporary workspaces",
        "legacy mdin",
        "bundled mdin",
        "topology/protocol/restart H5",
        "rerun trajectory H5",
        "Missing fixture files fail with explicit messages",
    ]
    require_tokens(body, "Phase 0 Fixture Helper plan", required_tokens)
    return body


def parse_phase1_input_contract_tests(plan_text):
    match = re.search(
        r"## Phase 1: Input Contract Tests\s+"
        r"(?P<body>.*?)(?:\n## Phase 2: Minimal Normal-Mode Smoke Matrix|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase 1 Input Contract Tests from plan")
    body = match.group("body")
    required_tokens = [
        "tests/h5_bundle/test_h5_input_matrix_contract.cpp",
        "inspect resolved plans",
        "validation behavior",
        "sidecar injection",
        "without launching a full SPONGE run",
        "Legacy input",
        "No H5 input binding exists",
        "legacy_input_allowed == true",
        "Pure bundled input",
        "input_h5_topology_path",
        "input_h5_protocol_path",
        "input_h5_restart_path",
        "/parameters/sponge/files/legacy_sidecars",
        "No legacy sidecar key is injected",
        "Bundled input with legacy sidecars",
        "/parameters/sponge/files/legacy_sidecars/key",
        "/parameters/sponge/files/legacy_sidecars/path",
        "Sidecar paths resolve relative to the H5 container directory",
        "Allowed topology/protocol/restart sidecar keys",
        "Legacy override/conflict behavior",
        "absent from mdin succeeds",
        "same path is idempotent",
        "different path fails with a conflict error",
        "mdin.override_conflict.spg.toml",
        "H5 and legacy restart input mutual exclusion",
        "coordinate_in_file",
        "velocity_in_file",
        "rst7",
        "first diagnostic layer",
        "cheap to debug",
        "internal state and validation results",
        "not enough to claim legacy alignment",
        "H5 input binding discovery",
        "Legacy fallback and mutual exclusion",
        "Sidecar key allowlists",
        "Sidecar path resolution relative to the containing H5 file",
        "same-key different-path is rejected before simulation",
        "test_h5_input_matrix_contract",
        "test_h5_legacy_sidecar",
        "test_h5_input_fixture_equivalence",
        "test_h5_io_contract_manifest",
        "Pure bundled fixtures are proven not to carry legacy sidecar tables",
        "Bundled-with-sidecar fixtures are proven to carry sidecar tables",
        "materialized files",
    ]
    require_tokens(body, "Phase 1 Input Contract Tests plan", required_tokens)
    return body


def parse_phase2_minimal_normal_smoke_matrix(plan_text):
    match = re.search(
        r"## Phase 2: Minimal Normal-Mode Smoke Matrix\s+"
        r"(?P<body>.*?)(?:\n## Phase 3: Legacy Sidecar Smoke Tests|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail(
            "failed to parse Phase 2 Minimal Normal-Mode Smoke Matrix from plan"
        )
    body = match.group("body")
    required_tokens = [
        "tests/h5_bundle/test_h5_input_output_smoke_matrix.cpp",
        "Use `core_structural`",
        "legacy in",
        "legacy out",
        "bundled out",
        "bundled in",
        "bundled with sidecar in",
        "Do not include rerun or VDS",
        "SPONGE exits successfully",
        "Expected output files exist",
        "`mdout` or observable core fields align",
        "legacy baseline",
        "tolerance",
        "explicit `mdout` and `mdinfo` comparison",
        "does not keep default legacy trajectory/restart sidecar outputs",
        "Bundled restart output matches expected",
        "position",
        "velocity",
        "box",
        "step/time",
        "Bundled trajectory output matches expected",
        "frame count",
        "position/box",
        "actual SPONGE initialization",
        "output paths",
        "Run this phase only after Phase 1 passes",
        "small and deterministic",
        "one or two MD steps",
        "Treat process success without output comparison as insufficient",
        "Every input family has at least one successful SPONGE run",
        "`legacy`, `bundled`, and `bundled with sidecar`",
        "Both output families have at least one successful SPONGE run",
        "`legacy` and `bundled`",
        "Normal-mode bundled output is readable by the H5 readers used by tests",
    ]
    require_tokens(
        body, "Phase 2 Minimal Normal-Mode Smoke Matrix plan", required_tokens
    )
    return body


def parse_phase3_legacy_sidecar_smoke_tests(plan_text):
    match = re.search(
        r"## Phase 3: Legacy Sidecar Smoke Tests\s+"
        r"(?P<body>.*?)(?:\n## Phase 4: Rerun Smoke|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase 3 Legacy Sidecar Smoke Tests from plan")
    body = match.group("body")
    required_tokens = [
        "focused on sidecar and override semantics",
        "core_structural/bundled_input_with_legacy_sidecar",
        "Bundled mdin without explicit legacy sidecar keys",
        "Sidecars are injected from H5",
        "SPONGE run succeeds",
        "explicit same-key same-path legacy value",
        "idempotent and succeeds",
        "explicit same-key different-path legacy value",
        "Run fails before simulation",
        "conflicting key",
        "indicates a conflict",
        "Pure bundled input",
        "Deleting or hiding legacy sidecar files",
        "must not affect the run",
        "accidentally depend on",
        "compatibility sidecars",
        "without explicit legacy topology/protocol sidecar",
        "same-key same-path explicit override remains idempotent",
        "same-key different-path explicit override fails before simulation",
        "reports the conflicting key",
        "`legacy_sidecars` directory is absent",
    ]
    require_tokens(
        body, "Phase 3 Legacy Sidecar Smoke Tests plan", required_tokens
    )
    return body


def parse_phase4_rerun_smoke(plan_text):
    match = re.search(
        r"## Phase 4: Rerun Smoke\s+"
        r"(?P<body>.*?)(?:\n## Phase 5: Parameterized Matrix|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase 4 Rerun Smoke from plan")
    body = match.group("body")
    required_tokens = [
        "full_contract_rerun",
        "legacy rerun in",
        "pure bundled rerun in",
        "bundled with sidecar rerun in",
        "legacy out",
        "bundled out, VDS off",
        "bundled out, VDS on",
        "second-frame selection",
        "rerun_start",
        "rerun_strip",
        "rerun_frame_limit",
        "Trajectory input step/time and frame count",
        "Position and box values align",
        "Rerun output observables align",
        "Bundled trajectory output is readable",
        "expected frame metadata",
        "VDS output creates both the wrapper and shard files",
        "completed-frame sequence",
        "Pure bundled rerun native core mdout columns",
        "temperature",
        "urey_bradley",
        "dihedral",
        "Current runtime native-parity exceptions",
        "input_h5_restart_load",
        "Force H5 rerun restart loading to `structural`",
        "test_h5_restart_load_runtime_closure.cpp",
        "test_h5_vds_terminal_resume_smoke.cpp",
        "Legacy rerun input establishes the baseline",
        "Pure bundled rerun input can read H5 trajectory frames",
        "Bundled-with-sidecar rerun input can exercise the broader full-contract",
        "`input_h5_restart_load` dynamic, protocol, and full policies have targeted",
        "initialized metadynamics restart",
        "pure-bundled native custom pairwise/listed-force payloads",
        "Rerun frame selection is tested separately",
        "Bundled rerun output writes trajectory and observable H5 artifacts",
        "VDS-on rerun output writes both the wrapper and shard directory",
        "VDS terminal tail-shard failure is repaired to the complete prefix",
        "complete-prefix resume policy is a no-op",
    ]
    require_tokens(body, "Phase 4 Rerun Smoke plan", required_tokens)
    return body


def parse_phase5_parameterized_matrix(plan_text):
    match = re.search(
        r"## Phase 5: Parameterized Matrix\s+"
        r"(?P<body>.*?)(?:\n## Phase 6: Contract Coverage|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase 5 Parameterized Matrix from plan")
    body = match.group("body")
    body_flat = normalize(body)
    required_tokens = [
        "After Phases 1-4 are stable",
        "mode:",
        "normal",
        "rerun",
        "input:",
        "legacy",
        "bundled",
        "bundled_with_sidecar",
        "output:",
        "vds:",
        "Use an explicit table of legal cases",
        "instead of blindly generating a Cartesian product",
        "`vds` only applies to bundled trajectory output",
        "Legacy output uses `vds = N/A`",
        "Normal-mode cases that do not write trajectory output should not test VDS",
        "Rerun plus bundled trajectory output is the main VDS on/off target",
        "PR CI should run only the smoke subset",
        "Nightly or explicit CI should run the full matrix",
        "Minimum legal case set",
        "complete support matrix for behavior",
        "not necessarily the default runtime matrix",
        "test_h5_io_matrix_spec.cpp",
        "test_h5_matrix_plan_manifest.py",
        "static matrix spec",
        "runtime smoke case names",
        "represented either by runtime execution or by a static legality check",
        "Avoid adding illegal Cartesian rows",
        "`vds=on` is meaningful only for bundled trajectory output",
    ]
    require_tokens(
        body_flat, "Phase 5 Parameterized Matrix plan", required_tokens
    )
    return body


def parse_phase5_legal_cases(plan_text):
    body = parse_phase5_parameterized_matrix(plan_text)
    rows = []
    in_table = False
    for raw_line in body.splitlines():
        stripped = raw_line.strip()
        if stripped == "Minimum legal case set:":
            in_table = True
            continue
        if not in_table:
            continue
        if stripped.startswith("This is the complete support matrix"):
            break
        if not stripped.startswith("|"):
            continue
        if stripped.startswith("|---") or stripped.startswith("| Case "):
            continue
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if len(cells) != 5:
            fail(f"Phase 5 legal matrix row must have 5 cells: {stripped}")
        rows.append(tuple(cell.strip("`") for cell in cells))
    expected = [
        ("normal_legacy_in_legacy_out", "normal", "legacy", "legacy", "N/A"),
        ("normal_legacy_in_bundled_out", "normal", "legacy", "bundled", "off"),
        ("normal_bundled_in_legacy_out", "normal", "bundled", "legacy", "N/A"),
        (
            "normal_sidecar_in_legacy_out",
            "normal",
            "bundled with sidecar",
            "legacy",
            "N/A",
        ),
        (
            "normal_bundled_in_bundled_out",
            "normal",
            "bundled",
            "bundled",
            "off",
        ),
        (
            "normal_sidecar_in_bundled_out",
            "normal",
            "bundled with sidecar",
            "bundled",
            "off",
        ),
        ("rerun_legacy_in_legacy_out", "rerun", "legacy", "legacy", "N/A"),
        (
            "rerun_legacy_in_bundled_out_vds_off",
            "rerun",
            "legacy",
            "bundled",
            "off",
        ),
        (
            "rerun_legacy_in_bundled_out_vds_on",
            "rerun",
            "legacy",
            "bundled",
            "on",
        ),
        ("rerun_bundled_in_legacy_out", "rerun", "bundled", "legacy", "N/A"),
        (
            "rerun_sidecar_in_legacy_out",
            "rerun",
            "bundled with sidecar",
            "legacy",
            "N/A",
        ),
        (
            "rerun_bundled_in_bundled_out_vds_off",
            "rerun",
            "bundled",
            "bundled",
            "off",
        ),
        (
            "rerun_bundled_in_bundled_out_vds_on",
            "rerun",
            "bundled",
            "bundled",
            "on",
        ),
        (
            "rerun_sidecar_in_bundled_out_vds_off",
            "rerun",
            "bundled with sidecar",
            "bundled",
            "off",
        ),
        (
            "rerun_sidecar_in_bundled_out_vds_on",
            "rerun",
            "bundled with sidecar",
            "bundled",
            "on",
        ),
    ]
    if rows != expected:
        fail(
            f"Phase 5 legal matrix changed without guard update:\nactual={rows}"
        )
    return rows


def parse_audit_phase_coverage_rows(audit_text):
    match = re.search(
        r"## Phase Coverage\s+"
        r"\| Plan item \| Evidence \| Default status \|\n"
        r"\|---\|---\|---\|\n"
        r"(?P<body>.*?)(?:\n## |\Z)",
        audit_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase Coverage table from audit matrix")

    rows = []
    for line in match.group("body").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not stripped.startswith("|"):
            break
        cells = [cell.strip() for cell in stripped.strip("|").split("|")]
        if len(cells) != 3:
            fail(f"Phase Coverage row must have 3 cells: {stripped}")
        rows.append(tuple(cells))
    if not rows:
        fail("Phase Coverage table contains no rows")
    return rows


def parse_acceptance_criteria(plan_text):
    match = re.search(
        r"## Acceptance Criteria\s+"
        r"The new I/O path is considered aligned with legacy only when:\s+"
        r"(?P<body>.*?)(?:\n## |\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Acceptance Criteria from plan")

    criteria = []
    current = None
    for raw_line in match.group("body").splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue
        if stripped.startswith("- "):
            if current:
                criteria.append(current)
            current = stripped[2:].strip()
        elif current and raw_line[:1].isspace():
            current += " " + stripped
        elif current:
            criteria.append(current)
            current = None
    if current:
        criteria.append(current)
    if not criteria:
        fail("Acceptance Criteria block contains no bullets")
    return criteria


def parse_current_runtime_boundary(plan_text):
    match = re.search(
        r"## Current Runtime Boundary\s+"
        r"(?P<body>.*?)(?:\n## |\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Current Runtime Boundary from plan")
    return match.group("body")


def parse_phase6_required_buckets(plan_text):
    match = re.search(
        r"## Phase 6: Contract Coverage\s+"
        r".*?Required coverage buckets:\s+"
        r"(?P<body>.*?)(?:\n\n|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase 6 Required coverage buckets from plan")

    buckets = []
    for raw_line in match.group("body").splitlines():
        stripped = raw_line.strip()
        if stripped.startswith("- "):
            buckets.append(stripped[2:].rstrip("."))
    if not buckets:
        fail("Phase 6 Required coverage buckets block contains no bullets")
    return buckets


def parse_phase7_native_parity_expansion(plan_text):
    match = re.search(
        r"## Phase 7: Full-Contract Native Parity Expansion\s+"
        r"(?P<body>.*?)(?:\n## CI Labels|\Z)",
        plan_text,
        re.DOTALL,
    )
    if not match:
        fail("failed to parse Phase 7 Full-Contract Native Parity Expansion")
    body = match.group("body")
    body_flat = normalize(body)
    required_tokens = [
        "closes the gap between",
        "pure bundled input runs",
        "fully equivalent to the legacy path",
        "sidecar-backed full-contract fixture",
        "complete legacy feature set through existing loaders",
        "replace one sidecar class at a time with native H5 readers",
        "extend pure-bundled parity assertions",
        "Topology core terms",
        "mass",
        "charge",
        "LJ",
        "bonds",
        "angles",
        "dihedrals",
        "NB14",
        "exclusions",
        "QC type mapping",
        "qc_type.txt",
        "/qc/type",
        "Custom pairwise/listed force payloads",
        "Enhanced sampling and protocol sidecars",
        "CV",
        "restraints",
        "SITS",
        "metadynamics",
        "Many-body and specialized force fields",
        "EAM",
        "SW",
        "EDIP",
        "TERSOFF",
        "ReaxFF",
        "Dynamic/protocol restart state",
        "Nose-Hoover chain",
        "hills",
        "restart-embedded sidecar text",
        "test_h5_restart_load_runtime_closure.cpp",
        "supported dynamic NHC",
        "protocol SITS sidecar",
        "full NHC+SITS",
        "Add or extend a reader-level unit test",
        "Add fixture semantic equivalence checks against the legacy source",
        "Add a pure-bundled smoke assertion for the observable terms affected",
        "Remove the sidecar dependency from the pure-bundled fixture",
        "Keep the bundled-with-sidecar case as compatibility coverage",
        "Pure bundled rerun compares all stable full-contract mdout/observable columns",
        "Remaining sidecar-backed payloads",
        "listed explicitly with the reason",
        "plan, manifest tests, static matrix, and runtime smoke case names agree",
    ]
    require_tokens(
        body_flat, "Phase 7 Native Parity Expansion plan", required_tokens
    )
    return body


def parse_coverage_manifest_buckets(coverage_text):
    buckets = set(
        re.findall(
            r'\{\s*"([^"]+)"\s*,\s*"[^"]+"\s*,\s*"[^"]+"\s*\}', coverage_text
        )
    )
    if not buckets:
        fail("failed to parse ManifestRequirement buckets from coverage test")
    return buckets


def parse_coverage_phase6_bucket_evidence(coverage_text):
    match = re.search(
        r"std::vector<Phase6BucketEvidence>\s+"
        r"Phase6_Coverage_Bucket_Evidence\(\)\s*\{"
        r".*?return\s*\{(?P<body>.*?)\};\s*\}",
        coverage_text,
        re.DOTALL,
    )
    if not match:
        fail(
            "failed to parse Phase6_Coverage_Bucket_Evidence from coverage test"
        )
    entries = re.findall(
        r'\{\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\}', match.group("body")
    )
    if not entries:
        fail("Phase6_Coverage_Bucket_Evidence contains no buckets")
    return entries


def require_plan_phases_are_audited(plan_text, audit_text):
    rows = parse_audit_phase_coverage_rows(audit_text)
    missing = []
    weak_rows = []
    for phase in parse_plan_phases(plan_text):
        phase_rows = [row for row in rows if phase in row[0]]
        if not phase_rows:
            missing.append(phase)
            continue
        for plan_item, evidence, default_status in phase_rows:
            if not evidence or not default_status:
                weak_rows.append(plan_item)
            if "test_" not in evidence and "tests/" not in evidence:
                weak_rows.append(plan_item)
            if not (
                default_status.startswith("Runs under")
                or "runtime smoke" in default_status
                or "Runtime smoke" in default_status
                or "Static" in default_status
                or "Contract runs" in default_status
                or "Registered by default" in default_status
            ):
                weak_rows.append(plan_item)
    if missing:
        fail(f"audit matrix does not mention plan phases: {missing}")
    if weak_rows:
        fail(
            f"audit matrix has phase rows without strong evidence/status: {weak_rows}"
        )


def require_implementation_file_map_is_bound(
    plan_text, audit_text, evidence_text
):
    expected = [
        (
            "0 fixture plumbing",
            [
                "h5_input_matrix_fixture.hpp",
                "test_h5_io_matrix_spec.cpp",
                "test_h5_io_contract_manifest.py",
            ],
            ["isolated workspaces", "relocatable"],
        ),
        (
            "1 input contracts",
            [
                "test_h5_input_matrix_contract.cpp",
                "test_h5_legacy_sidecar.cpp",
                "test_h5_input_validation.cpp",
            ],
            ["Legacy fallback", "sidecar injection", "mutual exclusion"],
        ),
        (
            "2 normal runtime smoke",
            [
                "test_h5_input_output_smoke_matrix.cpp",
                "test_h5_matrix_plan_manifest.py",
            ],
            [
                "normal input family",
                "default legacy trajectory/restart leakage",
            ],
        ),
        (
            "3 sidecar runtime smoke",
            [
                "test_h5_input_output_smoke_matrix.cpp",
                "test_h5_io_matrix_spec.cpp",
            ],
            ["without explicit legacy keys", "different-path override"],
        ),
        (
            "4 rerun and VDS smoke",
            [
                "test_h5_input_output_smoke_matrix.cpp",
                "test_vds_trajectory_writer_with_mock_backend.cpp",
                "test_highfive_backend_io.cpp",
            ],
            ["Rerun frame selection", "wrapper-relative source paths"],
        ),
        (
            "5 legal matrix guard",
            [
                "test_h5_io_matrix_spec.cpp",
                "test_h5_matrix_plan_manifest.py",
                "test_h5_ci_plan_manifest.py",
            ],
            ["15 legal case names", "runtime gate"],
        ),
        (
            "6 full-contract coverage",
            [
                "test_h5_io_contract_coverage.cpp",
                "test_h5_io_contract_manifest.py",
                "test_h5_input_fixture_equivalence.py",
            ],
            ["Manifest buckets", "semantic equivalence"],
        ),
        (
            "7 native parity expansion",
            [
                "Reader-specific tests",
                "test_h5_input_fixture_equivalence.py",
                "runtime smoke comparisons",
            ],
            ["native H5 reader coverage", "pure-bundled runtime parity"],
        ),
    ]
    rows = parse_implementation_file_map(plan_text)
    if len(rows) != len(expected):
        fail(
            "Implementation File Map row count changed without guard update: "
            f"actual={len(rows)} expected={len(expected)}"
        )

    for row, (phase, files, boundary_tokens) in zip(rows, expected):
        actual_phase, actual_files, actual_boundary = row
        if actual_phase != phase:
            fail(
                "Implementation File Map phase order changed without guard "
                f"update: actual={actual_phase!r} expected={phase!r}"
            )
        require_tokens(
            actual_files, f"Implementation File Map files: {phase}", files
        )
        require_tokens(
            actual_boundary,
            f"Implementation File Map boundary: {phase}",
            boundary_tokens,
        )

    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)
    require_tokens(
        audit_flat,
        "Implementation File Map audit coverage",
        [
            "Implementation File Map",
            "phase-to-file ownership",
            "test_h5_audit_matrix_manifest.py",
        ],
    )
    require_tokens(
        evidence_flat,
        "Implementation File Map source coverage",
        [
            "parse_implementation_file_map",
            "require_implementation_file_map_is_bound",
            "Implementation File Map row count changed without guard update",
        ],
    )


def require_input_fixture_generation_is_bound(
    plan_text, audit_text, evidence_text
):
    parse_input_fixture_generation(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Input Fixture Generation audit coverage",
        [
            "Input fixture generation emits stable",
            "legacy_input",
            "bundled_input",
            "bundled_input_with_legacy_sidecar",
            "manifest schema/version",
            "sidecar_key",
            "sidecar_path",
            "legacy_sidecars/<key>/<basename>",
            "byte-match",
            "relocation-safe",
        ],
    )
    require_tokens(
        evidence_flat,
        "Input Fixture Generation source coverage",
        [
            "require_fixture_group_structure",
            '"core_structural"',
            '"full_contract_rerun"',
            "require_no_legacy_input_refs",
            "require_sidecar_mdin_refs_are_materialized",
            "require_manifest_lists_h5_sidecar_tables",
            "require_manifest_lists_materialized_sidecar_files",
            "require_materialized_sidecars_match_legacy",
            "require_manifest_path_relocates",
            "require_required_inputs_have_semantic_equivalence_evidence",
            "compare_qc_type_to_legacy",
        ],
    )


def require_fixture_roots_are_bound(plan_text, audit_text, evidence_text):
    parse_fixture_roots(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Fixture Roots audit coverage",
        [
            "Phase 0 fixture helper",
            "tests/h5_bundle/h5_input_matrix_fixture.hpp",
            "Test_Fixture_Helper_Copies_Case_To_Temp",
            "Test_Fixture_Helper_Describes_Case_Paths",
            "Test_Fixture_Helper_Missing_Path_Has_Explicit_Message",
        ],
    )
    require_tokens(
        evidence_flat,
        "Fixture Roots source coverage",
        [
            "SPONGE_H5_INPUT_MATRIX_FIXTURE_ROOT",
            "Fixture_Root()",
            "Core_Structural_Path()",
            "Full_Contract_Rerun_Path()",
            '"core_structural"',
            '"full_contract_rerun"',
            "Copy_Case_To_Temp(",
            "Describe_Case(",
            "Require_Path_Exists(",
            "missing H5 input matrix fixture path",
            "Test_Fixture_Helper_Copies_Case_To_Temp",
            "Test_Fixture_Helper_Describes_Case_Paths",
            "Test_Fixture_Helper_Missing_Path_Has_Explicit_Message",
        ],
    )


def require_phase0_fixture_helper_is_bound(
    plan_text, audit_text, evidence_text
):
    parse_phase0_fixture_helper(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Phase 0 Fixture Helper audit coverage",
        [
            "Phase 0 fixture helper",
            "resolves fixture roots",
            "copies fixture cases into temp workspaces",
            "describes mdin, bundle input, and expected normal/rerun output paths",
            "reports missing fixture paths with explicit path-bearing errors",
        ],
    )
    require_tokens(
        evidence_flat,
        "Phase 0 Fixture Helper source coverage",
        [
            "Fixture_Root()",
            "Group_Path(",
            "Core_Structural_Path()",
            "Full_Contract_Rerun_Path()",
            "Normal_Output_Paths(",
            "Rerun_Output_Paths(",
            "Describe_Case(",
            "Copy_Case_To_Temp(",
            "Require_Path_Exists(",
            "std::filesystem::copy_options::recursive",
            "normal_output.h5_restart",
            "normal_output.h5_trajectory",
            "normal_output.h5_observable",
            "rerun_output.h5_restart.empty()",
            "rerun_output.h5_trajectory",
            "rerun_output.h5_observable",
            "Test_Fixture_Helper_Copies_Case_To_Temp",
            "Test_Fixture_Helper_Describes_Case_Paths",
            "Test_Fixture_Helper_Missing_Path_Has_Explicit_Message",
            "missing H5 input matrix fixture path",
        ],
    )


def require_phase1_input_contract_tests_are_bound(
    plan_text, audit_text, evidence_text
):
    parse_phase1_input_contract_tests(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Phase 1 Input Contract Tests audit coverage",
        [
            "Phase 1 legacy input",
            "Phase 1 pure bundled input",
            "Phase 1 any H5 input requires",
            "Phase 1 bundled input with legacy sidecars",
            "Phase 1 legacy override behavior",
            "Phase 1 H5 restart input rejects legacy restart inputs",
            "Phase 1 H5 input validation",
            "Test_Legacy_Input_Allows_Legacy_Path",
            "Test_Pure_Bundled_Input_Has_No_Legacy_Sidecars",
            "Test_Bundled_Input_With_Sidecar_Injects_Allowed_Keys",
            "Test_Legacy_Sidecar_Override_Conflict",
            "Test_H5_Restart_Rejects_Legacy_Restart_Inputs",
            "Test_Validates_Controller_Input_Bindings",
        ],
    )
    require_tokens(
        evidence_flat,
        "Phase 1 Input Contract Tests source coverage",
        [
            "Test_Legacy_Input_Allows_Legacy_Path",
            "Test_Empty_H5_Input_Paths_Do_Not_Disable_Legacy_Fallback",
            "Test_Pure_Bundled_Input_Has_No_Legacy_Sidecars",
            "Test_Bundled_Input_With_Sidecar_Injects_Allowed_Keys",
            "Test_Legacy_Sidecar_Override_Conflict",
            "Test_H5_Restart_Rejects_Legacy_Restart_Inputs",
            "Test_H5_Input_Requires_Topology_And_Protocol_Bindings",
            "Test_Normal_H5_Input_Requires_Restart_Binding",
            "Test_H5_Restart_Load_Policy_Validation",
            "Test_Validates_Controller_Input_Bindings",
            "Test_Controller_Input_Bindings_Propagate_Resolver_Errors",
            "Test_Accepts_Relative_Existing_Path_For_Same_Sidecar",
            "Test_Reads_And_Resolves_Relative_Sidecar_Paths",
            "Test_Rejects_Conflicting_Sidecar_Command",
            "Resolve_Input_Plan(&controller)",
            "plan.legacy_input_allowed",
            "plan.any_h5_input_enabled",
            "Has_H5_Input_Binding",
            "legacy_sidecars",
            "H5_Topology_Sidecar_Command_Keys",
            "H5_Protocol_Sidecar_Command_Keys",
            "coordinate_in_file",
            "velocity_in_file",
            "rst7",
            "conflicts with existing command",
            "missing required H5 input binding",
            "Validate_Input_Bindings(&controller)",
            "require_no_legacy_input_refs",
            "require_manifest_lists_h5_sidecar_tables",
            "require_materialized_sidecars_match_legacy",
        ],
    )


def require_phase2_minimal_normal_smoke_matrix_is_bound(
    plan_text, audit_text, evidence_text
):
    parse_phase2_minimal_normal_smoke_matrix(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Phase 2 Minimal Normal-Mode Smoke Matrix audit coverage",
        [
            "Phase 2 minimal normal-mode smoke",
            "legacy, bundled, and bundled-with-sidecar inputs",
            "legacy and bundled outputs",
            "Run_Normal_Mode_Matrix",
            "Require_Normal_Prepared_Mdin",
            "Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs",
            "Require_H5_Restart_Matches_Core_State",
            "Require_H5_Trajectory_Has_Frames",
            "Require_H5_Trajectory_First_Frame_Matches_Core_State",
            "Require_H5_Observable_Stream_Matches_Mdout",
        ],
    )
    require_tokens(
        evidence_flat,
        "Phase 2 Minimal Normal-Mode Smoke Matrix source coverage",
        [
            "Normal_Smoke_Cases(",
            '"normal_legacy_in_bundled_out"',
            '"normal_bundled_in_legacy_out"',
            '"normal_sidecar_in_legacy_out"',
            '"normal_bundled_in_bundled_out"',
            '"normal_sidecar_in_bundled_out"',
            '"normal_legacy_in_legacy_out"',
            "Run_Normal_Mode_Matrix(",
            "Prepare_Normal_Case(",
            "Require_Normal_Prepared_Mdin(",
            "Run_SPONGE(",
            "Require_Core_Mdout_Equivalent(baseline.mdout, test_case.mdout)",
            "Require_Text_Equivalent(baseline.mdout, test_case.mdout)",
            "if (spec.bundled_output)",
            "Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs(test_case)",
            "Require_H5_Restart_Matches_Core_State(test_case.h5_restart)",
            "Require_H5_Trajectory_Has_Frames(test_case.h5_trajectory, {1}",
            "Require_H5_Trajectory_First_Frame_Matches_Core_State(",
            "Require_H5_Observable_Stream_Matches_Mdout(",
            "Require_Normal_Legacy_Restart_Output(",
            "Test_Resolve_Legacy_Output_Plan_Matrix",
            "Test_Explicit_Legacy_Sidecar_Collection",
            '"mdcrd.dat"',
            '"mdbox.txt"',
            '"restart"',
        ],
    )


def require_phase3_legacy_sidecar_smoke_tests_are_bound(
    plan_text, audit_text, evidence_text
):
    parse_phase3_legacy_sidecar_smoke_tests(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Phase 3 Legacy Sidecar Smoke Tests audit coverage",
        [
            "Phase 3 sidecar smoke",
            "injected sidecars",
            "same-path override",
            "different-path conflict",
            "pure bundled input without sidecar files",
            "Run_Legacy_Sidecar_Smoke_Cases",
            "Sidecar_Smoke_Cases",
            "Require_Materialized_Sidecars_Are_Exactly_H5_Referenced",
            "Require_Normal_Legacy_Restart_Output",
        ],
    )
    require_tokens(
        evidence_flat,
        "Phase 3 Legacy Sidecar Smoke Tests source coverage",
        [
            "Sidecar_Smoke_Cases(",
            "Run_Legacy_Sidecar_Smoke_Cases(",
            "SidecarSmokeKind::injected_without_explicit_keys",
            "SidecarSmokeKind::same_key_same_path",
            "SidecarSmokeKind::same_key_different_path",
            "SidecarSmokeKind::pure_bundled_without_sidecar_files",
            '"sidecar_injected_without_explicit_legacy_keys"',
            '"sidecar_same_key_same_path"',
            '"sidecar_same_key_different_path"',
            '"pure_bundled_without_sidecar_files"',
            '"mdin.override_same_path.spg.toml"',
            '"mdin.override_conflict.spg.toml"',
            "Run_SPONGE(sponge_executable, test_case)",
            "Run_SPONGE_Expect_Failure(",
            "Require_Text_Equivalent(baseline_mdout, test_case.mdout)",
            "Require_Core_Mdout_Equivalent(baseline_mdout, test_case.mdout)",
            "Require_Normal_Legacy_Restart_Output(test_case)",
            "Remove_Legacy_Sidecar_Directories(test_case.root)",
            "Require_No_Legacy_Sidecar_Directories(prepared.root)",
            "Require_Materialized_Sidecars_Are_Exactly_H5_Referenced(",
            'mass_in_file = \\"override_mass.txt\\"',
            "legacy_sidecars/qc_type_in_file/qc_type.txt",
            "Native H5 topology data and legacy text topology input cannot both own atom masses",
        ],
    )


def require_phase4_rerun_smoke_is_bound(plan_text, audit_text, evidence_text):
    parse_phase4_rerun_smoke(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Phase 4 Rerun Smoke audit coverage",
        [
            "Phase 4 rerun smoke covers bundled and bundled-with-sidecar rerun input to legacy output",
            "Phase 4 rerun smoke covers legacy, bundled, and bundled-with-sidecar rerun input to bundled output",
            "Phase 4 full-contract rerun H5 input plan",
            "Phase 4 full-contract bundled mdin preserves rerun mode",
            "Phase 4 rerun frame selection obeys",
            "Phase 4 rerun trajectory input step/time",
            "Run_Rerun_Mode_Matrix",
            "Rerun_Bundled_Output_Cases",
            "Run_Rerun_Frame_Selection_Smoke_Cases",
            "Require_VDS_Shards_Are_Complete",
            "Scrub_Runtime_Unstable_Rerun_Features",
        ],
    )
    require_tokens(
        evidence_flat,
        "Phase 4 Rerun Smoke source coverage",
        [
            "Run_Rerun_Mode_Matrix(",
            '"rerun_legacy_in_legacy_out"',
            "Rerun_Legacy_Output_Cases(",
            '"rerun_bundled_in_legacy_out"',
            '"rerun_sidecar_in_legacy_out"',
            "Rerun_Bundled_Output_Cases(",
            '"rerun_legacy_in_bundled_out_vds_off"',
            '"rerun_legacy_in_bundled_out_vds_on"',
            '"rerun_bundled_in_bundled_out_vds_off"',
            '"rerun_bundled_in_bundled_out_vds_on"',
            '"rerun_sidecar_in_bundled_out_vds_off"',
            '"rerun_sidecar_in_bundled_out_vds_on"',
            "Run_Rerun_Frame_Selection_Smoke_Cases(",
            "Rerun_Selection_Cases(",
            '"rerun_sidecar_second_frame_only_legacy_out"',
            "const RerunSelection second_frame_only = {1, 0, 1}",
            "Require_Rerun_Selection_Mdout_Equivalent(",
            "Require_Pure_Bundled_Rerun_Mdout_Core_Equivalent(",
            "Require_Rerun_Mdout_Equivalent(",
            "Require_H5_Trajectory_Has_Frames(",
            "Require_H5_Trajectory_Frame_Matches_Rerun_Runtime_State(",
            "Require_H5_Observable_Stream_Matches_Mdout(",
            "Require_H5_Observable_Stream_Has_Frames(",
            "Require_VDS_Shards_Are_Complete(",
            "VDS shard position",
            "VDS shard box",
            "Scrub_Runtime_Unstable_Rerun_Features(",
            "Require_Runtime_Unstable_Rerun_Features_Scrubbed(",
            'input_h5_restart_load = \\"structural\\"',
            "Test_Rerun_H5_Trajectory_Input_Does_Not_Require_Restart_Binding",
            "Test_H5_Trajectory_Rejects_Legacy_Rerun_Inputs",
        ],
    )


def require_phase5_parameterized_matrix_is_bound(
    plan_text, audit_text, evidence_text
):
    parse_phase5_legal_cases(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Phase 5 Parameterized Matrix audit coverage",
        [
            "Phase 5 legal matrix is explicit instead of blind Cartesian product",
            "full plan-table case-name set and mode/input/output/VDS combination counts",
            "Test_Legal_Matrix_Cases_Are_Explicit_And_Fixture_Backed",
            "test_h5_matrix_plan_manifest.py",
            "runtime matrix dimensions",
            "Phase 5 fixture-backed bundled rerun mdin aligns",
            "Phase 5 VDS on/off is covered for rerun bundled trajectory output",
            "Rerun_Bundled_Output_Cases",
            "Require_VDS_Shards_Are_Complete",
            "test_h5_ci_plan_manifest.py",
            "runtime smoke gate",
        ],
    )
    require_tokens(
        evidence_flat,
        "Phase 5 Parameterized Matrix source coverage",
        [
            "parse_plan_cases",
            "Minimum legal case set",
            "static_matrix_case_matches",
            "parse_static_matrix_cases",
            "parse_static_matrix_evidence",
            "require_static_matrix_evidence_is_explicit",
            "current legal matrix cases must be backed by runtime smoke",
            "parse_runtime_expected_names",
            "parse_runtime_matrix_cases",
            "runtime matrix should expose 15 legal cases",
            "Normal_Smoke_Cases",
            "Rerun_Legacy_Output_Cases",
            "Rerun_Bundled_Output_Cases",
            "require_runtime_smoke_execution_links",
            "require_runtime_smoke_output_assertions",
            "Legal_Matrix_Cases(",
            "Test_Legal_Matrix_Cases_Are_Explicit_And_Fixture_Backed",
            "REQUIRE_EQ(cases.size(), static_cast<std::size_t>(15))",
            "REQUIRE_EQ(normal_cases, 6)",
            "REQUIRE_EQ(rerun_cases, 9)",
            "REQUIRE_EQ(vds_applicable_cases, 6)",
            "REQUIRE_EQ(vds_on_cases, 3)",
            "REQUIRE_EQ(vds_off_cases, 3)",
            "REQUIRE_EQ(runtime_smoke_evidence_cases, 15)",
            "expected_combination_counts",
            "Require_Input_Fixture_Exists",
            '"normal_legacy_in_legacy_out"',
            '"normal_sidecar_in_bundled_out"',
            '"rerun_sidecar_in_bundled_out_vds_on"',
            "parse_plan_labels",
            "require_plan_labels_registered",
            "require_runtime_smoke_gate",
            "require_smoke_preparation_before_gate",
            "SPONGE_H5_ENABLE_RUNTIME_SMOKE=1",
            "SKIP_RETURN_CODE 77",
        ],
    )


def require_phase6_buckets_are_bound(plan_text, audit_text, coverage_text):
    expected = [
        "Topology typed datasets",
        "Restart structural state",
        "Restart dynamic state",
        "Protocol sidecars",
        "SITS state and sidecars",
        "Metadynamics state and sidecars",
        "Custom pairwise/listed force payloads",
        "QC/ReaxFF sidecars",
        "Rerun trajectory input",
        "Legacy sidecar key/path tables",
        "Bundled output trajectory/restart/observable paths",
    ]
    actual = parse_phase6_required_buckets(plan_text)
    if actual != expected:
        fail(
            f"Phase 6 coverage buckets changed without test update:\nactual={actual}"
        )

    coverage_phase6_evidence = parse_coverage_phase6_bucket_evidence(
        coverage_text
    )
    coverage_phase6_buckets = [bucket for bucket, _ in coverage_phase6_evidence]
    if coverage_phase6_buckets != actual:
        fail(
            "Phase 6 coverage bucket evidence does not match plan:\n"
            f"actual={coverage_phase6_buckets}\nexpected={actual}"
        )

    coverage_buckets = parse_coverage_manifest_buckets(coverage_text)
    coverage_flat = normalize(coverage_text)
    audit_flat = normalize(audit_text)
    for bucket, evidence_token in coverage_phase6_evidence:
        if evidence_token not in coverage_flat:
            fail(
                "Phase 6 coverage bucket evidence token is not backed by "
                f"coverage source: bucket={bucket!r} evidence={evidence_token!r}"
            )
    require_tokens(
        coverage_text,
        "Phase 6 coverage bucket explicitness test",
        [
            "Test_Full_Contract_Rerun_Manifest_Buckets_Are_Explicit",
            "Test_Phase6_Plan_Buckets_Are_Represented",
            "Phase6_Coverage_Bucket_Evidence",
            "expected_buckets",
            "actual_buckets == expected_buckets",
        ],
    )
    evidence = {
        "Topology typed datasets": ["topology typed datasets"],
        "Restart structural state": ["restart structural state"],
        "Restart dynamic state": ["restart dynamic state"],
        "Protocol sidecars": ["protocol sidecars"],
        "SITS state and sidecars": ["SITS state and sidecars"],
        "Metadynamics state and sidecars": ["metadynamics state and sidecars"],
        "Custom pairwise/listed force payloads": ["custom force payloads"],
        "QC/ReaxFF sidecars": ["QC/ReaxFF sidecars"],
        "Rerun trajectory input": ["rerun trajectory input"],
        "Legacy sidecar key/path tables": [
            "Test_Full_Contract_Rerun_H5_Files_Cover_Sidecar_Tables",
            "Require_Legacy_Sidecar_Table",
        ],
        "Bundled output trajectory/restart/observable paths": [
            "bundled output paths",
            "output.h5.output_h5_trajectory_path",
            "output.h5.output_h5_restart_path",
            "output.h5.output_h5_observable_path",
        ],
    }
    for bucket, tokens in evidence.items():
        coverage_ok = False
        for token in tokens:
            if token in coverage_buckets or token in coverage_flat:
                coverage_ok = True
                break
        if not coverage_ok:
            fail(f"Phase 6 bucket lacks coverage-test evidence: {bucket}")
        if bucket not in audit_flat and not any(
            token in audit_flat for token in tokens
        ):
            fail(f"Phase 6 bucket lacks audit evidence: {bucket}")


def require_phase7_native_parity_expansion_is_bound(
    plan_text, audit_text, evidence_text
):
    parse_phase7_native_parity_expansion(plan_text)
    audit_flat = normalize(audit_text)
    evidence_flat = normalize(evidence_text)

    require_tokens(
        audit_flat,
        "Phase 7 Native Parity Expansion audit coverage",
        [
            "Phase 7 full-contract native parity expansion",
            "Current runtime native-parity exceptions",
            "pure bundled rerun native core mdout columns",
            "REAXFF/EDIP pure bundled closure",
            "`test_h5_matrix_plan_manifest.py` requires the exception table",
            "Scrub_Runtime_Unstable_Rerun_Features",
            "runtime preparation builds `rerun_runtime_scrub_prepare_check`",
            "Require_Runtime_Unstable_Rerun_Features_Scrubbed",
            "test_h5_audit_matrix_manifest.py",
        ],
    )
    require_tokens(
        evidence_flat,
        "Phase 7 Native Parity Expansion source coverage",
        [
            "parse_plan_runtime_exceptions",
            "parse_plan_pure_bundled_core_columns",
            "parse_smoke_runtime_exceptions",
            "parse_smoke_pure_bundled_core_columns",
            "require_pure_bundled_core_columns_match_plan",
            "require_runtime_exception_scrub_behavior",
            "Scrub_Runtime_Unstable_Rerun_Features(",
            "Require_Runtime_Unstable_Rerun_Features_Scrubbed(",
            "rerun_runtime_scrub_prepare_check",
            "Require_Pure_Bundled_Rerun_Mdout_Core_Equivalent(",
            "Require_Mdout_Columns_Equivalent(",
            '"temperature"',
            '"urey_bradley"',
            '"dihedral"',
            "compare_qc_type_to_legacy",
            "compare_custom_pair_to_h5",
            "compare_custom_bond_to_h5",
            "compare_manybody_pair_triple_to_h5",
            "compare_full_contract_protocol_typed",
            "compare_full_contract_restart_dynamic_typed",
            "Test_Reads_Native_Mass_And_Charge",
            "TopologyNativeH5Reader",
            "Test_Restart_Reader_Round_Trips_Structural_State",
            "Test_Restart_Reader_Round_Trips_Dynamic_And_Protocol_State",
            "RestartH5Reader",
            "Test_Trajectory_Reader_Reads_Metadata_And_Frame",
            "Test_Trajectory_Reader_Reads_Vds_Wrapper",
            "TrajectoryH5Reader",
        ],
    )


def require_drift_checks_are_audited(audit_text):
    rows = parse_audit_phase_coverage_rows(audit_text)
    coverage_text = normalize("\n".join(" | ".join(row) for row in rows))
    require_tokens(
        coverage_text,
        "Phase Coverage drift-check rows",
        [
            "test_h5_test_targets_manifest.py",
            "test_h5_matrix_plan_manifest.py",
            "test_h5_ci_plan_manifest.py",
            "test_h5_audit_matrix_manifest.py",
        ],
    )


def require_acceptance_criteria_are_audited(
    plan_text, audit_text, evidence_text
):
    expected = [
        "Internal plan and validation tests pass.",
        "Sidecar injection and conflict handling are explicitly covered.",
        "At least one normal-mode end-to-end smoke passes for each input family: legacy, bundled, and bundled with sidecar.",
        "At least one output smoke compares legacy output against bundled output.",
        "Rerun input is covered by an end-to-end smoke.",
        "Bundled trajectory output is covered with VDS off and on.",
        "The broad full-contract fixture has manifest-level coverage checks.",
        "The converter-generated sidecar files byte-match their legacy sources, and sidecar paths are relocation-safe.",
        "Legacy override behavior is covered for absent key, same-key same-path, and same-key different-path conflict.",
    ]
    actual = parse_acceptance_criteria(plan_text)
    if actual != expected:
        fail(
            f"Acceptance Criteria changed without audit update:\nactual={actual}"
        )

    audit_flat = normalize(audit_text)
    evidence = {
        expected[0]: [
            "test_h5_output_plan",
            "test_h5_input_matrix_contract",
            "test_h5_input_validation",
        ],
        expected[1]: [
            "Test_Bundled_Input_With_Sidecar_Injects_Allowed_Keys",
            "Test_Legacy_Sidecar_Override_Conflict",
        ],
        expected[2]: ["Run_Normal_Mode_Matrix"],
        expected[3]: [
            "Require_Core_Mdout_Equivalent",
            "Require_H5_Observable_Stream_Matches_Mdout",
        ],
        expected[4]: ["Run_Rerun_Mode_Matrix"],
        expected[5]: ["VDS on/off", "Rerun_Bundled_Output_Cases"],
        expected[6]: [
            "test_h5_io_contract_manifest.py",
            "test_h5_io_contract_coverage.cpp",
        ],
        expected[7]: ["byte-match", "relocation"],
        expected[8]: [
            "absent key",
            "same-path",
            "different-path conflict",
        ],
    }
    for criterion, tokens in evidence.items():
        require_tokens(
            audit_flat, f"Acceptance Criteria evidence: {criterion}", tokens
        )

    source_flat = normalize(evidence_text)
    source_evidence = {
        expected[0]: [
            "Test_Defaults_And_Legacy_Gating",
            "Test_Legacy_Input_Allows_Legacy_Path",
            "Test_Validates_Controller_Input_Bindings",
        ],
        expected[1]: [
            "Test_Bundled_Input_With_Sidecar_Injects_Allowed_Keys",
            "Test_Legacy_Sidecar_Override_Conflict",
            "Run_Legacy_Sidecar_Smoke_Cases",
        ],
        expected[2]: [
            "Run_Normal_Mode_Matrix",
            "normal_legacy_in_legacy_out",
            "normal_bundled_in_legacy_out",
            "normal_sidecar_in_legacy_out",
        ],
        expected[3]: [
            "Require_Core_Mdout_Equivalent",
            "Require_Text_Equivalent",
            "Require_H5_Observable_Stream_Matches_Mdout",
        ],
        expected[4]: [
            "Run_Rerun_Mode_Matrix",
            "rerun_legacy_in_legacy_out",
            "Rerun_Legacy_Output_Cases",
        ],
        expected[5]: [
            "Rerun_Bundled_Output_Cases",
            "Require_VDS_Shards_Are_Complete",
            "Require_H5_Trajectory_Has_Frames",
            '"VDS shard position"',
            '"VDS shard box"',
        ],
        expected[6]: [
            "Test_Phase6_Plan_Buckets_Are_Represented",
            "Require_Manifest_Entry_Fields",
            "Test_Full_Contract_Rerun_Legacy_Output_Sidecar_Plan_Is_Preserved",
            "require_required_inputs_have_semantic_equivalence_evidence",
            "REQUIRED_ENTRIES",
        ],
        expected[7]: [
            "Require_Materialized_Sidecars_Are_Exactly_H5_Referenced",
            "require_materialized_sidecars_match_legacy",
            "require_manifest_lists_materialized_sidecar_files",
            "require_manifest_path_relocates",
        ],
        expected[8]: [
            "injected_without_explicit_keys",
            "same_key_same_path",
            "same_key_different_path",
            "Run_SPONGE_Expect_Failure",
        ],
    }
    for criterion, tokens in source_evidence.items():
        require_tokens(
            source_flat,
            f"Acceptance Criteria source evidence: {criterion}",
            tokens,
        )


def require_current_runtime_boundary_is_audited(plan_text, audit_text):
    boundary = normalize(parse_current_runtime_boundary(plan_text))
    boundary_tokens = [
        "default suite can prove parser, contract, fixture, reader, writer, static matrix, sidecar provenance, and manifest behavior",
        "final legacy-alignment claim still requires the gated runtime smoke",
        "SPONGE_H5_ENABLE_RUNTIME_SMOKE=1",
        "test_h5_input_output_smoke_matrix",
        "tests/h5_bundle/run_h5_bundle_tests.sh test-smoke-runtime",
    ]
    require_tokens(boundary, "Current Runtime Boundary plan", boundary_tokens)

    require_tokens(
        normalize(audit_text),
        "Current Runtime Boundary audit",
        [
            "Runtime Gate",
            "parser, contract, fixture, reader, writer, static matrix, sidecar provenance, and manifest behavior",
            "End-to-end legacy alignment remains proven only when the runtime smoke gate succeeds",
            "SPONGE_H5_ENABLE_RUNTIME_SMOKE=1",
            "skip code 77",
            "tests/h5_bundle/run_h5_bundle_tests.sh test-smoke-runtime",
        ],
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--coverage-source", required=True, type=Path)
    parser.add_argument(
        "--evidence-source", action="append", default=[], type=Path
    )
    args = parser.parse_args()

    plan_text = args.plan.read_text()
    audit_text = args.audit.read_text()
    coverage_text = args.coverage_source.read_text()
    evidence_text = "\n".join(
        [coverage_text, Path(__file__).read_text()]
        + [source.read_text() for source in args.evidence_source]
    )
    audit_flat = normalize(audit_text)

    require_plan_phases_are_audited(plan_text, audit_text)
    require_implementation_file_map_is_bound(
        plan_text, audit_text, evidence_text
    )
    require_input_fixture_generation_is_bound(
        plan_text, audit_text, evidence_text
    )
    require_fixture_roots_are_bound(plan_text, audit_text, evidence_text)
    require_phase0_fixture_helper_is_bound(plan_text, audit_text, evidence_text)
    require_phase1_input_contract_tests_are_bound(
        plan_text, audit_text, evidence_text
    )
    require_phase2_minimal_normal_smoke_matrix_is_bound(
        plan_text, audit_text, evidence_text
    )
    require_phase3_legacy_sidecar_smoke_tests_are_bound(
        plan_text, audit_text, evidence_text
    )
    require_phase4_rerun_smoke_is_bound(plan_text, audit_text, evidence_text)
    require_phase5_parameterized_matrix_is_bound(
        plan_text, audit_text, evidence_text
    )
    require_phase6_buckets_are_bound(plan_text, audit_text, coverage_text)
    require_phase7_native_parity_expansion_is_bound(
        plan_text, audit_text, evidence_text
    )
    require_drift_checks_are_audited(audit_text)
    require_acceptance_criteria_are_audited(
        plan_text, audit_text, evidence_text
    )
    require_current_runtime_boundary_is_audited(plan_text, audit_text)

    require_tokens(
        audit_flat,
        "input family coverage",
        ["legacy_input", "bundled_input", "bundled_input_with_legacy_sidecar"],
    )
    require_tokens(
        audit_flat,
        "input resolver contract coverage",
        [
            "Test_Legacy_Input_Allows_Legacy_Path",
            "Test_Empty_H5_Input_Paths_Do_Not_Disable_Legacy_Fallback",
            "Test_H5_Restart_Load_Policy_Validation",
            "Test_H5_Input_Path_Suffix_Flags_Are_Non_Fatal",
            "Test_H5_Input_Requires_Topology_And_Protocol_Bindings",
            "Test_Normal_H5_Input_Requires_Restart_Binding",
            "Test_H5_Trajectory_Rejects_Legacy_Rerun_Inputs",
        ],
    )
    require_tokens(
        evidence_text,
        "input resolver contract source coverage",
        [
            "Test_Legacy_Input_Allows_Legacy_Path",
            "Test_Empty_H5_Input_Paths_Do_Not_Disable_Legacy_Fallback",
            "Test_H5_Input_Path_Suffix_Flags_Are_Non_Fatal",
            "Test_Normal_H5_Input_Requires_Restart_Binding",
            "Test_H5_Input_Requires_Topology_And_Protocol_Bindings",
            "Test_H5_Restart_Load_Policy_Validation",
            "Test_H5_Trajectory_Rejects_Legacy_Rerun_Inputs",
            "Test_H5_Trajectory_Is_Rerun_Only",
            "Test_Rerun_H5_Trajectory_Input_Does_Not_Require_Restart_Binding",
            "Test_Rerun_H5_Trajectory_Default_Stream_And_Mode_Case",
            "Resolve_Input_Plan(&controller)",
            "plan.legacy_input_allowed",
            "plan.any_h5_input_enabled",
            "Has_H5_Input_Binding",
            "missing required H5 input binding",
            "input_h5_topology_path",
            "input_h5_protocol_path",
            "input_h5_restart_path",
            "input_h5_restart_load",
            "RestartLoadPolicy::structural",
            "RestartLoadPolicy::dynamic",
            "RestartLoadPolicy::protocol",
            "RestartLoadPolicy::full",
            "input_h5_trajectory_path is currently only valid",
            "crd/box/vel",
            "coordinate/velocity restart inputs",
        ],
    )
    require_tokens(
        audit_flat,
        "sidecar and override coverage",
        [
            "legacy sidecar",
            "H5 sidecar-table keys",
            "sidecar_key",
            "sidecar_path",
            "Require_Materialized_Sidecars_Are_Exactly_H5_Referenced",
            "union of H5 sidecar-table references",
            "same-path",
            "different-path conflict",
            "Test_Legacy_Sidecar_Override_Conflict",
        ],
    )
    require_tokens(
        audit_flat,
        "v2 artifact identity and restart integrity coverage",
        [
            "v2 artifact identity and restart integrity",
            "sponge.input.v2",
            "/identity/uuid",
            "/run/topology_hash",
            "/run/atom_order_hash",
            "/run/producer_protocol_hash",
            "/run/state_hash",
            "cross-language integrity path",
        ],
    )
    require_tokens(
        evidence_text,
        "v2 artifact identity and restart integrity source coverage",
        [
            "Test_Rejects_Invalid_Artifact_Schema_And_Identity",
            "Test_Rejects_Restart_Lineage_And_State_Hash_Mismatch",
            "Refresh_Restart_State_Hash",
            "Compute_State_Hash(&state_hash)",
            "restart state_hash does not match",
        ],
    )
    require_tokens(
        evidence_text,
        "runtime sidecar exact-set source coverage",
        [
            "Require_Materialized_Sidecars_Are_Exactly_H5_Referenced(",
            "H5_Referenced_Sidecar_File_Set",
            "Materialized_Sidecar_File_Set",
            "Relative_Materialized_Sidecar_Path",
            "Materialized_Sidecar_File_Set(bundle_root) ==",
            "H5_Referenced_Sidecar_File_Set(bundle_root)",
        ],
    )
    require_tokens(
        audit_flat,
        "runtime VDS shard payload coverage",
        [
            "shard step/time arrays",
            "shard position/box payload slices",
            "VDS shard position",
            "VDS shard box",
        ],
    )
    require_tokens(
        evidence_text,
        "runtime VDS shard payload source coverage",
        [
            "Require_VDS_Shards_Are_Complete(",
            'Read_Float_Vector(file, "/particles/all/position/value")',
            'Read_Float_Vector(shard_file, "/particles/all/position/value")',
            'Read_Float_Vector(shard_file, "/particles/all/box/edges/value")',
            '"VDS shard position"',
            '"VDS shard box"',
        ],
    )
    require_tokens(
        audit_flat,
        "VDS wrapper-relative source path coverage",
        [
            "wrapper-relative VDS source paths",
            "particle, observable, and module virtual datasets",
        ],
    )
    require_tokens(
        evidence_text,
        "VDS wrapper-relative source path source coverage",
        [
            "Test_Vds_Source_Path_Relativization",
            "Test_Vds_Source_Path_Without_Wrapper_Parent",
            "Vds_Source_Path(",
            'plan.trajectory.path = "/tmp/sponge_h5_vds_case/wrappers/prod.spg.h5md"',
            "plan.trajectory.derived_shard_root =",
            '"../shards/prod.spg.shards/segment_000000.spg.h5md"',
            'plan.trajectory.path = "prod.spg.h5md"',
            '"prod.spg.shards/segment_000000.spg.h5md"',
            "sources[0].file_path",
        ],
    )
    require_tokens(
        audit_flat,
        "runtime bundled-output default legacy file coverage",
        [
            "Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs",
            "mdcrd.dat",
            "mdbox.txt",
            "restart",
            "default legacy",
        ],
    )
    require_tokens(
        audit_flat,
        "legacy output plan direct coverage",
        [
            "Test_Resolve_Legacy_Output_Plan_Matrix",
            "Test_Explicit_Legacy_Sidecar_Collection",
            "mdout",
            "mdinfo",
            "crd",
            "box",
            "vel",
            "frc",
            "rst",
            "qc_scf_output",
            "explicit sidecar collection order",
        ],
    )
    require_tokens(
        evidence_text,
        "legacy output plan direct source coverage",
        [
            "Test_Resolve_Legacy_Output_Plan_Matrix",
            "Test_Explicit_Legacy_Sidecar_Collection",
            "Resolve_Legacy_Output_Plan(",
            "Collect_Explicit_Legacy_Sidecars(",
            '"mdout"',
            '"mdinfo"',
            '"crd"',
            '"box"',
            '"vel"',
            '"frc"',
            '"rst"',
            '"qc_scf_output"',
            "legacy.default_enabled",
            "legacy.Enabled(legacy_keys[i])",
            "legacy.Explicitly_Requested(legacy_keys[i])",
            "keys[0]",
            "keys[1]",
            "keys[2]",
        ],
    )
    require_tokens(
        evidence_text,
        "runtime bundled-output default legacy file source coverage",
        [
            "Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs(",
            '"mdcrd.dat"',
            '"mdbox.txt"',
            '"restart"',
            "std::filesystem::exists(test_case.root / filename)",
            "Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs(test_case)",
            "Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs(sidecar_bundled)",
        ],
    )
    require_tokens(
        audit_flat,
        "runtime normal legacy-output restart coverage",
        [
            "Require_Normal_Legacy_Restart_Output",
            "normal legacy-output runtime rows",
            "Successful sidecar legacy-output smoke branches",
        ],
    )
    require_tokens(
        evidence_text,
        "runtime normal legacy-output restart source coverage",
        [
            "Require_Normal_Legacy_Restart_Output(",
            '"restart_coordinate.txt"',
            '"restart_velocity.txt"',
            "Require_Normal_Legacy_Restart_Output(baseline)",
            "Require_Normal_Legacy_Restart_Output(test_case)",
        ],
    )
    require_tokens(
        audit_flat,
        "phase 6 native semantic coverage",
        [
            "Require_Manifest_Entry_Fields",
            "Test_Full_Contract_Rerun_Legacy_Output_Sidecar_Plan_Is_Preserved",
            "exact topology/protocol/restart key sets",
            "legacy_sidecars/<key>/...",
            "topology typed datasets",
            "restart structural/dynamic/protocol payloads",
            "protocol typed datasets",
            "custom pairwise/listed force payloads",
            "many-body and specialized force fields",
            "ReaxFF",
        ],
    )
    require_tokens(
        audit_flat,
        "phase 6 legacy output sidecar coverage",
        [
            "legacy output sidecar preservation",
            "Test_Full_Contract_Rerun_Legacy_Output_Sidecar_Plan_Is_Preserved",
        ],
    )
    require_tokens(
        evidence_text,
        "phase 6 manifest field source coverage",
        [
            "Require_Manifest_Entry_Fields",
            '"component"',
            '"direction"',
            '"payload_kind"',
            '"override_policy"',
            '"bundle_file"',
            '"bundle_path"',
            '"source_key"',
            '"source_path"',
            '"sidecar_key"',
            '"sidecar_path"',
        ],
    )
    require_tokens(
        evidence_text,
        "phase 6 exact sidecar table source coverage",
        [
            "Require_Legacy_Sidecar_Table",
            "actual_keys == expected_keys",
            "sidecar_root",
            "relative.begin()->string() == sidecar.key",
            "mass_in_file",
            "REAXFF_in_file",
            "cv_in_file",
            "restrain_weight_in_file",
            "nose_hoover_chain_restart_input",
            "hills_in_file",
        ],
    )
    require_tokens(
        evidence_text,
        "phase 6 legacy output sidecar source coverage",
        [
            "Full_Contract_Legacy_Output_Sidecar_Requirements",
            "Test_Full_Contract_Rerun_Legacy_Output_Sidecar_Plan_Is_Preserved",
            "legacy_output_sidecar_preserved",
            "output.legacy_sidecar.mdout",
            "output.legacy_sidecar.crd",
            "output.legacy_sidecar.box",
            "output.legacy_sidecar.vel",
            "active_bundled_legacy_outputs",
            'std::set<std::string>{"mdout"}',
        ],
    )
    require_tokens(
        audit_flat,
        "typed restraint coverage",
        [
            "protocol.restrain_atom_id",
            "protocol.restrain_weight",
            "restart.restrain_coordinate",
        ],
    )
    require_tokens(
        audit_flat,
        "fixture semantic equivalence coverage",
        [
            "test_h5_input_fixture_equivalence.py",
            "mass, charge",
            "SITS",
            "metadynamics",
            "EAM, SW, EDIP, TERSOFF, and ReaxFF",
        ],
    )
    require_tokens(
        audit_flat,
        "runtime exception coverage",
        [
            "Current runtime native-parity exceptions",
            "Scrub_Runtime_Unstable_Rerun_Features",
            "rerun_runtime_scrub_prepare_check",
            "Require_Runtime_Unstable_Rerun_Features_Scrubbed",
            "test_h5_reaxff_edip_runtime_parity.cpp",
            "test_h5_restart_load_runtime_closure.cpp",
            "test_h5_vds_terminal_resume_smoke.cpp",
            "topology_custom_force_h5_materializer.hpp",
            "manybody_legacy_in_legacy_out_no_manybody_scrub",
            "manybody_sidecar_in_legacy_out_no_manybody_scrub",
            "manybody_sidecar_in_bundled_out_vds_*_no_manybody_scrub",
            "manybody_pure_bundled_in_legacy_out_no_sidecar",
            "manybody_pure_bundled_in_bundled_out_vds_*_no_sidecar",
            "protocol SITS sidecar",
            "SITS `nk`/`log_norm`/`log_nk` device-state restoration",
            "NHC dynamic state",
            "initialized metadynamics restart loading",
            "pure-bundled native custom-force direct initialization",
        ],
    )
    require_tokens(
        evidence_text,
        "runtime scrub preparation source coverage",
        [
            "rerun_runtime_scrub_prepare_check",
            "Require_H5_Legacy_Sidecar_Keys_Present(",
            "Scrub_Runtime_Unstable_Rerun_Features(",
            "Require_Runtime_Unstable_Rerun_Features_Scrubbed(",
            'input_h5_restart_load = \\"structural\\"',
            "Manybody_Columns()",
            "Require_Manybody_Not_Scrubbed(",
            "restart_load_protocol_sidecar_meta_initialized",
            "START INITIALIZING 1D-META",
            "restart_load_protocol_native_sits_full_state",
            "Install_Full_Sits_Restart_State(",
            "native_sits_log_norm",
            "native_sits_log_nk",
            "Read Nk from native H5 restart",
            ".sponge_h5_restart_protocol",
            "restart_load_protocol_pure_bundled_custom_force_native",
            "_assert_pure_bundle_has_no_native_bridge_artifacts(",
            ".sponge_h5_native_custom_force",
            ".sponge_h5_native_qc",
            ".sponge_h5_native_protocol",
            ".sponge_h5_native_manybody",
            "QC =",
            "START INITIALIZING STILLINGER WEBER FORCE FROM NATIVE H5",
            "START INITIALIZING TERSOFF FORCE FROM NATIVE H5",
            "START INITIALIZING EAM FORCE FROM NATIVE H5",
            "START INITIALIZING PAIRWISE FORCE FROM NATIVE H5",
            "START INITIALIZING LISTED FORCES FROM NATIVE H5",
            "Read_Pairwise(",
            "Read_Listed(",
            "Read_ReaxFF(",
            "Initial_ReaxFF_From_Native(",
            "START INITIALIZING REAXFF_EEQ (native H5)",
            "reaxff_type.txt",
        ],
    )
    require_tokens(
        audit_flat,
        "plan runtime exception entries",
        parse_plan_runtime_exceptions(plan_text),
    )
    require_tokens(
        audit_flat,
        "pure bundled rerun native core mdout column boundary",
        parse_plan_pure_bundled_core_columns(plan_text),
    )


if __name__ == "__main__":
    main()
