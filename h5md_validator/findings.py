"""Finding types emitted by the validator."""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict, field
from enum import Enum


class Severity(str, Enum):
    #: The file violates a MUST of the specification.  A conformant file has none.
    ERROR = "ERROR"
    #: A known specification-vs-implementation mismatch, tracked in
    #: :mod:`h5md_validator.deviations`.  Reported separately so that it does
    #: not drown out genuine errors while the mismatch is unresolved.
    DEVIATION = "DEVIATION"
    #: The file violates a SHOULD / RECOMMENDED of the specification.
    WARNING = "WARNING"
    #: Informational, no conformance consequence.
    INFO = "INFO"


@dataclass
class Finding:
    #: Stable check identifier, e.g. ``TOPO-INVARIANT-PARTICLES``.
    code: str
    severity: Severity
    #: HDF5 path the finding is anchored to.
    path: str
    message: str
    #: Free-form supporting values (expected/actual, indices, ...).
    detail: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["severity"] = self.severity.value
        return d

    def __str__(self) -> str:
        return f"[{self.severity.value}] {self.code} at {self.path}: {self.message}"


@dataclass
class Report:
    """Result of validating one file."""

    file: str
    findings: list = field(default_factory=list)

    def add(self, code, severity, path, message, **detail):
        self.findings.append(Finding(code, severity, path, message, detail))

    def by_severity(self, severity: Severity) -> list:
        return [f for f in self.findings if f.severity is severity]

    @property
    def errors(self) -> list:
        return self.by_severity(Severity.ERROR)

    @property
    def deviations(self) -> list:
        return self.by_severity(Severity.DEVIATION)

    @property
    def warnings(self) -> list:
        return self.by_severity(Severity.WARNING)

    @property
    def conformant(self) -> bool:
        """True when no ERROR-level finding was raised.

        Deviations are deliberately excluded: they are mismatches between the
        specification text and the reference implementation that are pending an
        editorial decision, not defects in the file under test.
        """
        return not self.errors

    def to_dict(self) -> dict:
        return {
            "file": self.file,
            "conformant": self.conformant,
            "counts": {
                s.value: len(self.by_severity(s)) for s in Severity
            },
            "findings": [f.to_dict() for f in self.findings],
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def summary(self) -> str:
        c = {s.value: len(self.by_severity(s)) for s in Severity}
        verdict = "PASS" if self.conformant else "FAIL"
        return (
            f"{verdict}  {self.file}  "
            f"errors={c['ERROR']} deviations={c['DEVIATION']} "
            f"warnings={c['WARNING']}"
        )
