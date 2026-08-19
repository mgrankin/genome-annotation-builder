"""Command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .acmg_resources import sync_acmg_sf, verify_acmg_sf
from .clingen_resources import sync_clingen_resources, verify_clingen_resources
from .clinvar_resources import sync_clinvar, verify_clinvar
from .config import load_config
from .constants import PIPELINE_VERSION
from .errors import GenomeBuilderError
from .gwas_resources import sync_gwas_catalog, verify_gwas_catalog
from .pgs import verify_pgs_resources
from .pgx_resources import sync_pgx_resources, verify_pgx_resources
from .pharmcat_runner import sync_pharmcat_runtime, verify_pharmcat_runtime
from .predictor_doctor import run_predictor_doctor, verify_predictor_doctor
from .predictor_resources import sync_predictor_resources, verify_predictor_resources
from .preflight import run_preflight
from .resources import sync_pinned_resources, verify_pinned_resources
from .vep_resources import sync_vep_cache, verify_vep_cache


def _json_print(value: Any) -> None:
    print(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="genome-annotation-builder",
        description="Build and validate evidence-preserving .genome/1.1 bundles.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {PIPELINE_VERSION}")
    subcommands = parser.add_subparsers(dest="command", required=True)

    resources = subcommands.add_parser("resources", help="manage checksum-pinned resources")
    resource_commands = resources.add_subparsers(dest="resource_command", required=True)
    for name in (
        "sync",
        "verify",
        "vep-sync",
        "vep-verify",
        "clinvar-sync",
        "clinvar-verify",
        "clingen-sync",
        "clingen-verify",
        "acmg-sf-sync",
        "acmg-sf-verify",
        "pgs-verify",
        "pgx-verify",
        "pgx-sync",
        "pharmcat-sync",
        "pharmcat-verify",
        "gwas-sync",
        "gwas-verify",
        "predictors-sync",
        "predictors-verify",
        "predictors-doctor",
        "predictors-doctor-verify",
    ):
        command = resource_commands.add_parser(name)
        command.add_argument("--config", required=True, type=Path)

    preflight = subcommands.add_parser("preflight", help="validate inputs and environment")
    preflight.add_argument("--config", required=True, type=Path)

    # Registered early so the public CLI is stable; implementation is added in
    # the same release before a bundle is advertised as buildable.
    build = subcommands.add_parser("build", help="build and atomically validate a bundle")
    build.add_argument("--config", required=True, type=Path)

    validate = subcommands.add_parser("validate", help="validate an existing bundle")
    validate.add_argument("bundle", type=Path)
    validate.add_argument("--reference", type=Path)

    query = subcommands.add_parser("query", help="run a bounded read-only bundle query")
    query.add_argument("bundle", type=Path)
    query.add_argument(
        "operation",
        choices=[
            "overview",
            "search",
            "variant",
            "gene",
            "clinical",
            "pgx",
            "prs",
            "gwas",
            "callability",
            "evidence",
        ],
    )
    query.add_argument("term", nargs="?")
    query.add_argument("--limit", type=int, default=50)

    release_c = subcommands.add_parser(
        "release-c",
        help="run checksum-bound read-backed Release C operations",
    )
    release_c.add_argument("release_c_arguments", nargs=argparse.REMAINDER)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        if arguments.command == "resources":
            config = load_config(arguments.config)
            handlers = {
                "sync": sync_pinned_resources,
                "verify": verify_pinned_resources,
                "vep-sync": sync_vep_cache,
                "vep-verify": verify_vep_cache,
                "clinvar-sync": sync_clinvar,
                "clinvar-verify": verify_clinvar,
                "clingen-sync": sync_clingen_resources,
                "clingen-verify": verify_clingen_resources,
                "acmg-sf-sync": sync_acmg_sf,
                "acmg-sf-verify": verify_acmg_sf,
                "pgs-verify": verify_pgs_resources,
                "pgx-verify": verify_pgx_resources,
                "pgx-sync": sync_pgx_resources,
                "pharmcat-sync": sync_pharmcat_runtime,
                "pharmcat-verify": verify_pharmcat_runtime,
                "gwas-sync": sync_gwas_catalog,
                "gwas-verify": verify_gwas_catalog,
                "predictors-sync": sync_predictor_resources,
                "predictors-verify": verify_predictor_resources,
                "predictors-doctor": run_predictor_doctor,
                "predictors-doctor-verify": verify_predictor_doctor,
            }
            result = handlers[arguments.resource_command](config)
            _json_print(result)
            return 0
        if arguments.command == "preflight":
            _json_print(run_preflight(load_config(arguments.config)))
            return 0
        if arguments.command == "build":
            from .pipeline import build_bundle

            _json_print(build_bundle(load_config(arguments.config)))
            return 0
        if arguments.command == "validate":
            from .validator import validate_bundle

            _json_print(validate_bundle(arguments.bundle, reference_fasta=arguments.reference))
            return 0
        if arguments.command == "query":
            from .query import query_bundle

            _json_print(
                query_bundle(
                    arguments.bundle,
                    arguments.operation,
                    term=arguments.term,
                    limit=arguments.limit,
                )
            )
            return 0
        if arguments.command == "release-c":
            from .release_c.__main__ import main as release_c_main

            return release_c_main(arguments.release_c_arguments)
    except (GenomeBuilderError, FileNotFoundError, FileExistsError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
