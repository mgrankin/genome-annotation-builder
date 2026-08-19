"""Pipeline-specific exceptions with stable exit semantics."""


class GenomeBuilderError(RuntimeError):
    """Base class for user-facing build failures."""


class ConfigurationError(GenomeBuilderError):
    """Configuration is missing, inconsistent, or unsafe."""


class ResourceError(GenomeBuilderError):
    """A pinned resource is unavailable or has the wrong digest."""


class InputValidationError(GenomeBuilderError):
    """An immutable input failed structural or reference validation."""


class ReconciliationError(GenomeBuilderError):
    """Accepted source evidence and generated output are not one-to-one."""


class BundleValidationError(GenomeBuilderError):
    """A staged or existing bundle is not conformant."""
