# Findings

The validator (Validation 3) reports each problem it finds in a file as a
finding; Validation 4 checks that each corruption raises its expected code.
Findings are written to the `findings` list of each file's entry in
`conformance.json`:

| Field | Meaning |
|:--|:--|
| `code`     | the stable check code, from the list below |
| `severity` | one of the three levels below |
| `path`     | the HDF5 path it is anchored to (empty for the file itself) |
| `message`  | a human-readable explanation; it may change, so do not parse it |
| `detail`   | supporting values, such as expected and actual |

Each entry's `counts` gives the number per severity, and `conformant` the
verdict.

## Check codes

Findings carry a stable code so they can be cited and tracked.

| Prefix | Family |
|:--|:--|
| `MOD`    | module declaration and version |
| `STRUCT` | presence, type and shape of specified elements |
| `ENUM`   | values of enumerated elements |
| `UNIT`, `UNITS` | unit attributes and their notation |
| `COND`   | conditional subgroup instantiation |
| `TIME`   | step/time/value alignment |
| `TOPO`   | topology invariants |
| `REF`    | referential integrity of indices and selection names |
| `SEL`, `CONN` | selection and connectivity bounds, coupling-group coverage |
| `EXT`    | extensions of the specification where permitted |
| `FILE`   | the file cannot be read at all |

## Severities

| Severity | Meaning | Effect on the verdict |
|:--|:--|:--|
| `ERROR`     | violates a MUST of the specification | the file is not conformant |
| `WARNING`   | violates a SHOULD, e.g. an unknown unit symbol or an algorithm outside the predefined set | none |
| `INFO`      | informational, e.g. an interaction type outside Table SI that carries a reference | none |

A file is conformant when it has no `ERROR` finding. Each check sets its own
severity when it raises the finding.
