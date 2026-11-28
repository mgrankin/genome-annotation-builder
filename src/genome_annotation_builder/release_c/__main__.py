"""Unified dispatcher for the generic Release C command modules."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Callable, Sequence

from . import (
    aldy_adapter,
    binding_builder,
    bundle_extension,
    caller_lifecycle,
    caller_profile,
    caller_receipt,
    cnv_concordance,
    cnvpytor_adapter,
    cnvpytor_call_capture,
    cnvpytor_combined_adapter,
    cnvpytor_combined_hdf5_probe,
    cnvpytor_hdf5_probe,
    cyrius_adapter,
    delly_adapter,
    expansionhunter_adapter,
    extension_stage,
    hla_concordance,
    hla_la_adapter,
    input_gate,
    manifest,
    manta_adapter,
    matrix_runner,
    mt_adapter,
    mt_concordance,
    mt_workflow,
    oci_runner,
    pgx_concordance,
    profile_adapter,
    profile_catalog,
    project_release,
    readback_adapter,
    readback_capture_v2,
    readback_concordance,
    readback_targets,
    repeat_concordance,
    resource_manifest,
    resource_tree,
    spechla_adapter,
    subject_stage,
    sv_concordance,
)

Handler = Callable[[], int]

COMMANDS: dict[str, tuple[str, Handler]] = {
    "manifest": ("build or validate the private input manifest", manifest.main),
    "input-gate": ("rehash and validate one caller input", input_gate.main),
    "binding-builder": (
        "compile private mount bindings from a public profile",
        binding_builder.main,
    ),
    "caller-lifecycle": (
        "run or resume one complete caller lifecycle",
        caller_lifecycle.main,
    ),
    "caller-profile": ("prepare an immutable caller execution plan", caller_profile.main),
    "profile-catalog": ("validate and list shipped caller profiles", profile_catalog.main),
    "oci-runner": ("execute a prepared OCI caller plan", oci_runner.main),
    "caller-receipt": ("bind execution and raw-artifact provenance", caller_receipt.main),
    "profile-adapter": (
        "materialize receipt-bound typed caller evidence",
        profile_adapter.main,
    ),
    "resource-tree": ("freeze or validate a resource tree", resource_tree.main),
    "resource-manifest": (
        "bind validated resources to the exact profile catalog",
        resource_manifest.main,
    ),
    "manta-adapter": ("materialize Manta structural evidence", manta_adapter.main),
    "delly-adapter": ("materialize DELLY structural evidence", delly_adapter.main),
    "cnvpytor-adapter": ("materialize CNVpytor read-depth evidence", cnvpytor_adapter.main),
    "cnvpytor-combined-adapter": (
        "materialize CNVpytor combined RD+BAF evidence",
        cnvpytor_combined_adapter.main,
    ),
    "cnvpytor-call-capture": (
        "capture and audit CNVpytor combined calls",
        cnvpytor_call_capture.main,
    ),
    "cnvpytor-hdf5-probe": ("audit CNVpytor RD HDF5 rows", cnvpytor_hdf5_probe.main),
    "cnvpytor-combined-hdf5-probe": (
        "audit CNVpytor combined HDF5 rows",
        cnvpytor_combined_hdf5_probe.main,
    ),
    "expansionhunter-adapter": (
        "materialize ExpansionHunter repeat evidence",
        expansionhunter_adapter.main,
    ),
    "hla-la-adapter": ("materialize HLA*LA evidence", hla_la_adapter.main),
    "spechla-adapter": ("materialize SpecHLA evidence", spechla_adapter.main),
    "cyrius-adapter": ("materialize licensed Cyrius output", cyrius_adapter.main),
    "aldy-adapter": ("materialize licensed Aldy output", aldy_adapter.main),
    "mitochondrial-adapter": ("materialize mitochondrial evidence", mt_adapter.main),
    "mitochondrial-workflow": (
        "run pinned normal/shifted mitochondrial calling",
        mt_workflow.main,
    ),
    "readback-adapter": ("materialize bounded read-back evidence", readback_adapter.main),
    "readback-capture": (
        "capture bounded read evidence without identifiers",
        readback_capture_v2.main,
    ),
    "readback-targets": (
        "select and validate bounded read-back targets",
        readback_targets.main,
    ),
    "readback-concordance": (
        "annotate technical-replicate read-back concordance",
        readback_concordance.main,
    ),
    "sv-concordance": ("annotate structural-event concordance", sv_concordance.main),
    "cnv-concordance": ("annotate copy-number concordance", cnv_concordance.main),
    "repeat-concordance": ("annotate repeat-call concordance", repeat_concordance.main),
    "hla-concordance": ("annotate HLA concordance", hla_concordance.main),
    "pgx-concordance": ("annotate complex-PGx concordance", pgx_concordance.main),
    "mt-concordance": ("annotate mitochondrial concordance", mt_concordance.main),
    "extension-stage": ("stage or validate portable extension tables", extension_stage.main),
    "bundle-extension": (
        "atomically enrich or validate a .genome/1.1 bundle",
        bundle_extension.main,
    ),
    "matrix-runner": (
        "run or validate the sequential caller matrix",
        matrix_runner.main,
    ),
    "subject-stage": (
        "assemble one subject from a completed caller matrix",
        subject_stage.main,
    ),
    "project-release": (
        "build or validate deterministic multi-bundle reports",
        project_release.main,
    ),
}


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genome-annotation-builder release-c",
        description="Run checksum-bound, read-backed Release C operations.",
    )
    parser.add_argument("command", nargs="?", choices=sorted(COMMANDS))
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    lines = [f"  {name:<34} {description}" for name, (description, _) in COMMANDS.items()]
    parser.epilog = "commands:\n" + "\n".join(lines)
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    parser = _parser()
    if not raw or raw[0] in {"-h", "--help"}:
        parser.print_help()
        return 0
    command = raw[0]
    if command not in COMMANDS:
        parser.error(f"unknown Release C command: {command}")
    handler = COMMANDS[command][1]
    original = sys.argv
    try:
        sys.argv = [f"{parser.prog} {command}", *raw[1:]]
        return handler()
    finally:
        sys.argv = original


if __name__ == "__main__":
    raise SystemExit(main())
