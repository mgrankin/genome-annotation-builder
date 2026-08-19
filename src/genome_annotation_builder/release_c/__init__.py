"""Read-backed evidence adapters and immutable Release C staging."""

from .portable_model import DOMAINS, PortableEvidenceError, schema_for

RELEASE_C_SCHEMA_VERSION = 1

__all__ = [
    "DOMAINS",
    "RELEASE_C_SCHEMA_VERSION",
    "PortableEvidenceError",
    "schema_for",
]
