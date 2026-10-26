"""
Finding types emitted by the validator.

One :class:`Finding` per problem, carrying a check code and a severity; a
:class:`Report` collects the findings of one file and decides its verdict --
conformant means no ERROR.

Used by:
    Validation 3: internal consistency
    Validation 4: negative tests
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from enum import Enum


class Severity(str, Enum):
    """How serious a finding is; only ERROR decides the verdict.

    The severity is chosen by the check that raises the finding, at its
    ``rep.add`` call in :mod:`h5md_validator.validate`.
    """

    #: The file violates a MUST of the specification.  A conformant file has none.
    ERROR = "ERROR"
    #: The file violates a SHOULD / RECOMMENDED of the specification, e.g. an
    #: unknown unit symbol or an algorithm name outside the predefined set.
    WARNING = "WARNING"
    #: Informational, no conformance consequence, e.g. an interaction type
    #: outside Table SI that carries a reference.
    INFO = "INFO"


@dataclass
class Finding:
    """One problem found in one file, as written to ``validation.json``."""

    #: Stable check identifier, e.g. ``TOPO-INVARIANT-PARTICLES``; the
    #: families are listed in docs/check-codes.md.
    code: str
    severity: Severity
    #: HDF5 path the finding is anchored to ("" for the file itself).
    path: str
    #: Human-readable explanation; not stable, do not parse it.
    message: str
    #: Free-form supporting values (expected/actual, indices, ...).
    detail: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """JSON-ready form, with the severity as its string value."""
        d = asdict(self)
        d["severity"] = self.severity.value
        return d

    def __str__(self) -> str:
        return f"[{self.severity.value}] {self.code} at {self.path}: {self.message}"


@dataclass
class Report:
    """Result of validating one file: every finding, in the order raised.

    Checks only record findings; the verdict is derived afterwards by
    :attr:`conformant`.
    """

    file: str
    findings: list = field(default_factory=list)

    def add(self, code, severity, path, message, **detail):
        """Record one finding; keyword arguments become its ``detail``."""
        self.findings.append(Finding(code, severity, path, message, detail))

    def by_severity(self, severity: Severity) -> list:
        """The findings raised at one severity."""
        return [f for f in self.findings if f.severity is severity]

    @property
    def errors(self) -> list:
        return self.by_severity(Severity.ERROR)

    @property
    def warnings(self) -> list:
        return self.by_severity(Severity.WARNING)

    @property
    def conformant(self) -> bool:
        """True when no ERROR-level finding was raised."""
        return not self.errors

    def to_dict(self) -> dict:
        """The ``conformance`` section of a ``validation.json`` entry."""
        return {
            "file": self.file,
            "conformant": self.conformant,
            "counts": {s.value: len(self.by_severity(s)) for s in Severity},
            "findings": [f.to_dict() for f in self.findings],
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def summary(self) -> str:
        """One line: PASS/FAIL, the file, and the counts per severity."""
        c = {s.value: len(self.by_severity(s)) for s in Severity}
        verdict = "PASS" if self.conformant else "FAIL"
        return f"{verdict}  {self.file}  errors={c['ERROR']} warnings={c['WARNING']}"
