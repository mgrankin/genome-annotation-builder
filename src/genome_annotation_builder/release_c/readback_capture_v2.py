"""Capture v2 adds the canonical zero-target selection to the pinned v1 pileup.

The OCI profile mounts both files in /tools; the installed CLI imports them as a package.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pysam

if __package__:
    from . import readback_capture as v1
else:
    import readback_capture as v1

EMPTY_TARGET_MANIFEST = ("\t".join(v1.TARGET_COLUMNS) + "\n").encode()
EMPTY_TARGET_SHA256 = hashlib.sha256(EMPTY_TARGET_MANIFEST).hexdigest()


def build_capture(
    alignment_path: str | Path,
    reference_fasta: str | Path,
    target_manifest: str | Path,
    *,
    policy: v1.CapturePolicy = v1.DEFAULT_CAPTURE_POLICY,
) -> dict[str, Any]:
    try:
        return v1.build_capture(alignment_path, reference_fasta, target_manifest, policy=policy)
    except v1.CaptureError as exc:
        if str(exc) != "target_manifest_empty":
            raise

    # v1 validated the policy and TSV schema. Only the selector's exact empty
    # representation can produce a successful zero-record capture.
    if v1._sha256_file(Path(target_manifest)) != EMPTY_TARGET_SHA256:
        raise v1.CaptureError("empty_target_identity_invalid")
    if not Path(str(reference_fasta) + ".fai").is_file():
        raise v1.CaptureError("reference_index_unusable")
    with (
        pysam.FastaFile(str(reference_fasta)) as reference,
        pysam.AlignmentFile(
            str(alignment_path), reference_filename=str(reference_fasta)
        ) as alignment,
    ):
        if not reference.references or not alignment.has_index():
            raise v1.CaptureError("alignment_or_reference_index_unusable")
        bindings = [str(r["SM"]) for r in alignment.header.to_dict().get("RG", []) if r.get("SM")]
    payload = {
        "schema_version": v1.SCHEMA_VERSION,
        "status": "PASS",
        "assembly": "GRCh38",
        "target_manifest_sha256": EMPTY_TARGET_SHA256,
        "policy": asdict(policy),
        "records": [],
        "record_count": 0,
        "sample_identifier_persisted": False,
        "read_names_persisted": False,
        "read_sequences_persisted": False,
        "private_paths_persisted": False,
    }
    if any(binding.encode() in v1._canonical_bytes(payload) for binding in bindings):
        raise v1.CaptureError("sample_identifier_persisted")
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    for name in ("alignment", "reference-fasta", "target-manifest", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    policy_defaults = asdict(v1.DEFAULT_CAPTURE_POLICY)
    for name, default in policy_defaults.items():
        flag = "--" + name.replace("_", "-")
        if isinstance(default, bool):
            parser.add_argument(flag, action="store_true", default=default)
        else:
            parser.add_argument(flag, type=int, default=default)
    args = parser.parse_args()
    try:
        policy = v1.CapturePolicy(**{name: getattr(args, name) for name in policy_defaults})
        result = v1.write_capture(
            args.output,
            build_capture(
                args.alignment, args.reference_fasta, args.target_manifest, policy=policy
            ),
        )
        print(json.dumps(result, sort_keys=True))
    except Exception:
        print(
            "read-back capture failed; inspect the authorized private process log", file=sys.stderr
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
