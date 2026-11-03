"""
Validation toolkit for the H5MD ``biophysics`` metadata schema.

Validates against specification 1.0.0-alpha (files stamped ``[1, 0]``); other
versions are reported as ``MOD-VERSION-UNSUPPORTED``. A new specification
release that changes elements needs a matching release of this package.
"""

from .findings import Finding, Report, Severity
from .spec import MODULE_VERSION, SPEC_VERSION
from .validate import Validator, validate

__all__ = [
    "MODULE_VERSION",
    "SPEC_VERSION",
    "Finding",
    "Report",
    "Severity",
    "Validator",
    "validate",
]

#: Version of this package, not of the specification it implements.
__version__ = "1.0.0"
