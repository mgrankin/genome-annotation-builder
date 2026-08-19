from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from datetime import date

import pytest

from genome_annotation_builder.cpic_freeze import (
    GITHUB_REF_URL,
    GITHUB_RELEASE_URL,
    OPENAPI_URL,
    PAIR_FIELDS,
    PAIR_VIEW_URL,
    JsonResponse,
    freeze_cpic_pairs,
)
from genome_annotation_builder.errors import ResourceError
from genome_annotation_builder.util import canonical_json_bytes, sha256_bytes


def _pair(pair_id: int, drug: str = "example-drug") -> dict[str, object]:
    values: dict[str, object] = {
        "pairid": pair_id,
        "drugid": f"drug-{pair_id}",
        "drugname": drug,
        "genesymbol": "GENE1",
        "guidelinename": "Example guideline",
        "guidelineurl": "https://example.test/guideline",
        "cpiclevel": "A",
        "clinpgxlevel": "1A",
        "pgxtesting": "Actionable PGx",
        "pmids": ["12345"],
        "usedforrecommendation": "Yes",
        "provisional": False,
    }
    assert tuple(values) == PAIR_FIELDS
    return values


def _openapi() -> dict[str, object]:
    return {
        "swagger": "2.0",
        "info": {"title": "PostgREST API", "version": "12.0.2"},
        "basePath": "/v1/",
        "definitions": {"pair_view": {"type": "object"}},
        "paths": {"/pair_view": {"get": {}}},
    }


def _fetcher(*, mutate_second_pass: bool = False, content_range: str = "0-1/2"):
    pair_calls = 0

    def fetch(url: str, headers: Mapping[str, str]) -> JsonResponse:
        nonlocal pair_calls
        if url == GITHUB_RELEASE_URL.format(version="1.2.3"):
            return JsonResponse(
                value={
                    "id": 42,
                    "tag_name": "v1.2.3",
                    "published_at": "2026-08-01T00:00:00Z",
                    "html_url": "https://github.com/cpicpgx/cpic-data/releases/tag/v1.2.3",
                    "draft": False,
                    "prerelease": False,
                },
                headers={},
                status=200,
            )
        if url == GITHUB_REF_URL.format(version="1.2.3"):
            return JsonResponse(
                value={"object": {"type": "commit", "sha": "a" * 40}},
                headers={},
                status=200,
            )
        if url == OPENAPI_URL:
            return JsonResponse(value=_openapi(), headers={}, status=200)
        if url.startswith(PAIR_VIEW_URL + "?"):
            assert headers["Prefer"] == "count=exact"
            assert headers["Range"] == "0-999"
            pair_calls += 1
            rows = [_pair(2), _pair(1)]
            if mutate_second_pass and pair_calls == 2:
                rows[0] = _pair(2, "changed-drug")
            return JsonResponse(
                value=rows,
                headers={"content-range": content_range},
                status=200,
            )
        raise AssertionError(f"unexpected URL: {url}")

    return fetch


def test_cpic_freeze_binds_release_api_and_stable_rowset() -> None:
    document = freeze_cpic_pairs("1.2.3", date(2026, 8, 15), fetcher=_fetcher())
    assert document["schema_version"] == 2
    assert document["source_release"]["commit"] == "a" * 40
    assert document["api_identity"]["postgrest_version"] == "12.0.2"
    assert document["stable_retrieval_passes"] == 2
    assert [row["pairid"] for row in document["rows"]] == [1, 2]
    assert document["rowset_canonical_sha256"] == sha256_bytes(
        canonical_json_bytes(document["rows"])
    )


def test_cpic_freeze_rejects_live_rowset_change_between_passes() -> None:
    with pytest.raises(ResourceError, match="changed between"):
        freeze_cpic_pairs(
            "1.2.3",
            date(2026, 8, 15),
            fetcher=_fetcher(mutate_second_pass=True),
        )


def test_cpic_freeze_rejects_false_content_range() -> None:
    with pytest.raises(ResourceError, match="Content-Range"):
        freeze_cpic_pairs(
            "1.2.3",
            date(2026, 8, 15),
            fetcher=_fetcher(content_range="0-0/2"),
        )


def test_cpic_freeze_rejects_moved_release_tag() -> None:
    base = _fetcher()

    def moved(url: str, headers: Mapping[str, str]) -> JsonResponse:
        response = base(url, headers)
        if url == GITHUB_REF_URL.format(version="1.2.3"):
            value = deepcopy(response.value)
            value["object"]["sha"] = "not-a-commit"
            return JsonResponse(value=value, headers=response.headers, status=response.status)
        return response

    with pytest.raises(ResourceError, match="commit identity"):
        freeze_cpic_pairs("1.2.3", date(2026, 8, 15), fetcher=moved)
