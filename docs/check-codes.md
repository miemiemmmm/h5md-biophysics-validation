# Check codes

Findings carry a stable code so they can be cited and tracked.

| Prefix | Family |
|:--|:--|
| `MOD` | module declaration and version |
| `STRUCT` | presence, type and shape of specified elements |
| `ENUM` | values of enumerated elements |
| `UNIT`, `UNITS` | unit attributes and their notation |
| `COND` | conditional subgroup instantiation |
| `TIME` | step/time/value alignment |
| `TOPO` | topology invariants |
| `REF` | referential integrity of indices |
| `SEL`, `CONN` | selection and connectivity bounds |
| `FILE` | the file cannot be read at all |
| `DEV-*` | registered specification deviation |

## Findings

Each finding is one record in the `conformance.findings` list of
`validation.json`:

| Field | Meaning |
|:--|:--|
| `code` | the stable check code, from the families above |
| `severity` | one of the four below |
| `path` | the HDF5 path it is anchored to (empty for the file itself) |
| `message` | a human-readable explanation; it may change, so do not parse it |
| `detail` | supporting values, such as expected and actual |

`conformance.counts` gives the number per severity, and `conformance.conformant`
the verdict.

## Severities

| Severity | Meaning | Effect on the verdict |
|:--|:--|:--|
| `ERROR` | violates a MUST of the specification | the file is not conformant |
| `DEVIATION` | a registered mismatch between the specification and the writer, pending an editorial decision (`h5md_validator/deviations.py`) | none |
| `WARNING` | violates a SHOULD, e.g. an unknown unit symbol or an algorithm outside the predefined set | none |
| `INFO` | informational, e.g. an interaction type outside Table SI that carries a reference | none |

A file is conformant when it has no `ERROR` finding. Each check sets its own
severity when it raises the finding.
