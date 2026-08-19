"""Deterministically freeze the small CPIC gene-drug pair view.

The live CPIC endpoint does not expose a database-release identifier in each
response.  A snapshot therefore records, without conflating them:

* the operator-selected official CPIC data release and its resolved Git tag;
* the OpenAPI schema actually served during retrieval;
* an exact, stable, twice-retrieved row set from ``pair_view``.

The second retrieval and before/after OpenAPI checks catch a changing live API
during a freeze.  The resulting file remains a checksum-pinned input; none of
these metadata fields are treated as a substitute for its SHA-256.
"""

from __future__ import annotations

import json
import re
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from datetime import date
from typing import Any

from .errors import ResourceError
from .util import canonical_json_bytes, sha256_bytes

PAIR_VIEW_URL = "https://api.cpicpgx.org/v1/pair_view"
OPENAPI_URL = "https://api.cpicpgx.org/v1/"
GITHUB_RELEASE_URL = "https://api.github.com/repos/cpicpgx/cpic-data/releases/tags/v{version}"
GITHUB_REF_URL = "https://api.github.com/repos/cpicpgx/cpic-data/git/ref/tags/v{version}"

PAIR_FIELDS = (
    "pairid",
    "drugid",
    "drugname",
    "genesymbol",
    "guidelinename",
    "guidelineurl",
    "cpiclevel",
    "clinpgxlevel",
    "pgxtesting",
    "pmids",
    "usedforrecommendation",
    "provisional",
)
_VERSION_RE = re.compile(r"[0-9]+\.[0-9]+\.[0-9]+")
_SHA1_RE = re.compile(r"[0-9a-f]{40}")
_CONTENT_RANGE_RE = re.compile(r"([0-9]+)-([0-9]+)/([0-9]+)")


@dataclass(frozen=True)
class JsonResponse:
    """Decoded response plus only the headers needed for validation."""

    value: Any
    headers: Mapping[str, str]
    status: int


JsonFetcher = Callable[[str, Mapping[str, str]], JsonResponse]


def fetch_json(url: str, headers: Mapping[str, str]) -> JsonResponse:
    """Retrieve JSON with bounded time and stable user-facing failures."""

    request_headers = {
        "User-Agent": "genome-annotation-builder-cpic-freeze/0.2",
        **headers,
    }
    request = urllib.request.Request(url, headers=request_headers)
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            status = int(response.status)
            payload = response.read()
            response_headers = {key.lower(): value for key, value in response.headers.items()}
    except (urllib.error.URLError, TimeoutError) as exc:
        raise ResourceError(f"failed to retrieve CPIC snapshot resource {url}: {exc}") from exc
    if status not in {200, 206}:
        raise ResourceError(f"unexpected HTTP status {status} for {url}")
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ResourceError(f"CPIC snapshot resource is not valid JSON: {url}") from exc
    return JsonResponse(value=value, headers=response_headers, status=status)


def _canonical_digest(value: Any) -> str:
    return sha256_bytes(canonical_json_bytes(value))


def _release_identity(source_version: str, fetcher: JsonFetcher) -> dict[str, Any]:
    release_url = GITHUB_RELEASE_URL.format(version=source_version)
    ref_url = GITHUB_REF_URL.format(version=source_version)
    release = fetcher(release_url, {"Accept": "application/vnd.github+json"}).value
    ref = fetcher(ref_url, {"Accept": "application/vnd.github+json"}).value
    expected_tag = f"v{source_version}"
    if not isinstance(release, dict) or release.get("tag_name") != expected_tag:
        raise ResourceError(f"official CPIC release metadata does not resolve {expected_tag}")
    if release.get("draft") is True or release.get("prerelease") is True:
        raise ResourceError(f"CPIC release {expected_tag} is not a final published release")
    published_at = release.get("published_at")
    html_url = release.get("html_url")
    release_id = release.get("id")
    if not isinstance(published_at, str) or not published_at:
        raise ResourceError(f"CPIC release {expected_tag} has no publication timestamp")
    if not isinstance(html_url, str) or not html_url.startswith(
        "https://github.com/cpicpgx/cpic-data/releases/"
    ):
        raise ResourceError(f"CPIC release {expected_tag} has an unexpected canonical URL")
    if not isinstance(release_id, int):
        raise ResourceError(f"CPIC release {expected_tag} has no numeric release ID")
    if not isinstance(ref, dict):
        raise ResourceError(f"official CPIC Git ref does not resolve {expected_tag}")
    ref_object = ref.get("object")
    if not isinstance(ref_object, dict) or ref_object.get("type") != "commit":
        raise ResourceError(f"CPIC tag {expected_tag} does not resolve directly to a commit")
    commit = ref_object.get("sha")
    if not isinstance(commit, str) or _SHA1_RE.fullmatch(commit) is None:
        raise ResourceError(f"CPIC tag {expected_tag} has an invalid commit identity")
    return {
        "tag": expected_tag,
        "commit": commit,
        "github_release_id": release_id,
        "published_at": published_at,
        "release_url": html_url,
    }


def _openapi_identity(fetcher: JsonFetcher) -> dict[str, Any]:
    response = fetcher(OPENAPI_URL, {"Accept": "application/openapi+json"})
    document = response.value
    if not isinstance(document, dict) or document.get("swagger") != "2.0":
        raise ResourceError("CPIC API did not return the expected OpenAPI 2.0 document")
    if document.get("basePath") != "/v1/":
        raise ResourceError("CPIC OpenAPI basePath is not /v1/")
    info = document.get("info")
    if not isinstance(info, dict) or not isinstance(info.get("version"), str):
        raise ResourceError("CPIC OpenAPI document lacks its serving software version")
    definitions = document.get("definitions")
    paths = document.get("paths")
    if not isinstance(definitions, dict) or "pair_view" not in definitions:
        raise ResourceError("CPIC OpenAPI document has no pair_view definition")
    if not isinstance(paths, dict) or "/pair_view" not in paths:
        raise ResourceError("CPIC OpenAPI document has no pair_view route")
    return {
        "openapi_canonical_sha256": _canonical_digest(document),
        "postgrest_version": info["version"],
        "swagger_version": document["swagger"],
        "base_path": document["basePath"],
    }


def _parse_content_range(value: str | None) -> tuple[int, int, int]:
    if value is None:
        raise ResourceError("CPIC pair_view response lacks exact Content-Range")
    match = _CONTENT_RANGE_RE.fullmatch(value.strip())
    if match is None:
        raise ResourceError(f"invalid CPIC pair_view Content-Range: {value!r}")
    return tuple(int(group) for group in match.groups())  # type: ignore[return-value]


def _retrieve_rows(fetcher: JsonFetcher, *, page_size: int = 1000) -> list[dict[str, Any]]:
    query = urllib.parse.urlencode(
        {"select": ",".join(PAIR_FIELDS), "order": "pairid.asc"}, safe=","
    )
    url = f"{PAIR_VIEW_URL}?{query}"
    rows: list[dict[str, Any]] = []
    total: int | None = None
    while total is None or len(rows) < total:
        start = len(rows)
        stop = start + page_size - 1
        response = fetcher(
            url,
            {
                "Accept": "application/json",
                "Range": f"{start}-{stop}",
                "Range-Unit": "items",
                "Prefer": "count=exact",
            },
        )
        page = response.value
        if not isinstance(page, list) or not all(isinstance(row, dict) for row in page):
            raise ResourceError("CPIC pair_view response is not an array of objects")
        if not page:
            raise ResourceError("CPIC pair_view ended before its declared exact row count")
        range_start, range_end, range_total = _parse_content_range(
            response.headers.get("content-range")
        )
        if range_start != start or range_end != start + len(page) - 1:
            raise ResourceError("CPIC pair_view Content-Range does not match response rows")
        if total is None:
            total = range_total
        elif total != range_total:
            raise ResourceError("CPIC pair_view exact row count changed during retrieval")
        if len(rows) + len(page) > total:
            raise ResourceError("CPIC pair_view returned more rows than its exact count")
        rows.extend(page)
    if not rows:
        raise ResourceError("CPIC pair_view returned no rows")
    return rows


def _validate_and_sort_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    expected_fields = set(PAIR_FIELDS)
    for index, row in enumerate(rows, start=1):
        if set(row) != expected_fields:
            raise ResourceError(f"CPIC pair row {index} has an unexpected field boundary")
        try:
            int(row["pairid"])
        except (TypeError, ValueError) as exc:
            raise ResourceError(f"CPIC pair row {index} has an invalid pairid") from exc
    ordered = sorted(rows, key=lambda row: int(row["pairid"]))
    ids = [int(row["pairid"]) for row in ordered]
    if len(ids) != len(set(ids)):
        raise ResourceError("CPIC pair_view contains duplicate pairid values")
    return ordered


def freeze_cpic_pairs(
    source_version: str,
    snapshot_date: date,
    *,
    fetcher: JsonFetcher = fetch_json,
) -> dict[str, Any]:
    """Return a deterministic, stability-checked CPIC pair snapshot."""

    if _VERSION_RE.fullmatch(source_version) is None:
        raise ResourceError("CPIC source version must be an exact X.Y.Z data release")
    release = _release_identity(source_version, fetcher)
    openapi_before = _openapi_identity(fetcher)
    first_rows = _validate_and_sort_rows(_retrieve_rows(fetcher))
    second_rows = _validate_and_sort_rows(_retrieve_rows(fetcher))
    first_digest = _canonical_digest(first_rows)
    second_digest = _canonical_digest(second_rows)
    if first_digest != second_digest:
        raise ResourceError("CPIC pair_view changed between the two snapshot passes")
    openapi_after = _openapi_identity(fetcher)
    if openapi_before != openapi_after:
        raise ResourceError("CPIC OpenAPI identity changed during snapshot retrieval")
    return {
        "schema_version": 2,
        "source": "CPIC PostgREST pair_view",
        "source_url": PAIR_VIEW_URL,
        "source_version": source_version,
        "snapshot_date": snapshot_date.isoformat(),
        "fields": list(PAIR_FIELDS),
        "row_count": len(first_rows),
        "rowset_canonical_sha256": first_digest,
        "stable_retrieval_passes": 2,
        "source_release": release,
        "api_identity": openapi_before,
        "rows": first_rows,
    }
