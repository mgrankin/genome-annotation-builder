from __future__ import annotations

from pathlib import Path

import pytest

from genome_annotation_builder.release_c.caller_lifecycle import run_caller_lifecycle
from genome_annotation_builder.release_c.profile_adapter import (
    ProfileAdapterError,
    materialize_profile_adapter,
    validate_profile_adapter_output,
)

from .test_caller_lifecycle import (
    PRIVATE_DATASET,
    PRIVATE_LAB,
    PRIVATE_SUBJECT,
    install_manta_runtime,
    private_manifest,
)
from .test_caller_profile import make_fixture


def test_receipt_bound_profile_adapter_materializes_replays_and_detects_tamper(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    fixture = make_fixture(tmp_path / "inputs")
    manifest = private_manifest(tmp_path, fixture)
    lifecycle = tmp_path / "lifecycle"
    lifecycle.mkdir()
    (lifecycle / "input-gate.json").write_bytes(fixture["gate"].read_bytes())
    commands: list[list[str]] = []
    install_manta_runtime(monkeypatch, fixture, commands)
    monkeypatch.setattr(
        "genome_annotation_builder.release_c.profile_adapter.resolve_shipped_profile",
        lambda _profile_id: fixture["profile"],
    )
    run_caller_lifecycle(
        manifest,
        PRIVATE_DATASET,
        fixture["bindings"],
        fixture["resources"],
        lifecycle,
        profile_path=fixture["profile"],
        minimum_free_bytes=1,
    )

    output = tmp_path / "adapter-output"
    first = materialize_profile_adapter(
        lifecycle,
        fixture["bindings"],
        "primary",
        output,
        profile_id="manta-germline-1.6.0",
    )
    assert first["status"] == "PASS"
    assert first["resume"] == "created"
    assert first["typed_files"] == 1
    assert first["typed_rows"] == 3
    portable = (output / "adapter-summary.json").read_text(encoding="utf-8")
    for private in (
        PRIVATE_DATASET,
        PRIVATE_SUBJECT,
        PRIVATE_LAB,
        str(manifest),
        str(fixture["alignment"]),
        str(fixture["workspace"]),
    ):
        assert private not in portable

    replay = materialize_profile_adapter(
        lifecycle,
        fixture["bindings"],
        "primary",
        output,
        profile_id="manta-germline-1.6.0",
    )
    assert replay["resume"] == "reused"
    assert replay["adapter_identity_sha256"] == first["adapter_identity_sha256"]

    typed = output / "manta.parquet"
    typed.write_bytes(typed.read_bytes() + b"changed")
    with pytest.raises(ProfileAdapterError, match="adapter_output_file_identity_changed"):
        validate_profile_adapter_output(
            output,
            lifecycle,
            fixture["bindings"],
            "primary",
            profile_id="manta-germline-1.6.0",
        )


def test_only_version_two_readback_profile_accepts_canonical_empty_capture(tmp_path: Path) -> None:
    from genome_annotation_builder.release_c import profile_adapter, readback_capture_v2
    from genome_annotation_builder.release_c.caller_profile import load_caller_profile
    from genome_annotation_builder.release_c.portable_model import PortableEvidenceError
    from genome_annotation_builder.release_c.profile_catalog import (
        CURRENT_READBACK,
        LEGACY_READBACK,
        resolve_shipped_profile,
    )
    from genome_annotation_builder.release_c.readback_capture import write_capture
    from genome_annotation_builder.release_c.readback_targets import _render_targets

    from .test_readback_adapter import build_inputs

    alignment, reference, targets, provenance = build_inputs(tmp_path)
    targets.write_bytes(_render_targets([]))
    capture = tmp_path / "empty-capture.json"
    write_capture(capture, readback_capture_v2.build_capture(alignment, reference, targets))
    for profile_id in (LEGACY_READBACK, CURRENT_READBACK):
        context = {
            "profile": load_caller_profile(resolve_shipped_profile(profile_id)),
            "components": {"raw_capture_json": capture},
            "boundary": {
                field: getattr(provenance, field)
                for field in (
                    "reference_sha256",
                    "caller_run_identity_sha256",
                    "resource_manifest_sha256",
                    "source_alignment_sha256",
                )
            },
            "artifact": {"raw_artifact_id": provenance.raw_artifact_id},
        }
        args = (profile_id, tmp_path, context, {}, provenance.technical_source_role, {}, ())
        if profile_id == LEGACY_READBACK:
            with pytest.raises(PortableEvidenceError):
                profile_adapter._materialize(*args)
        else:
            [(domain, path, result)] = profile_adapter._materialize(*args)
            assert domain == "read_backed_evidence" and result["status"] == "PASS"
            assert result["parsing"]["typed_rows"] == 0 and result["output"] is None
            assert not path.exists()
