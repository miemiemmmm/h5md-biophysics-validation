# Findings

The validator (Validation 3) reports each problem it finds in a file as a
finding; Validation 4 checks that each corruption raises its expected code.
Findings are written to the `conformance.findings` list of `validation.json`:

| Field | Meaning |
|:--|:--|
| `code`     | the stable check code, from the list below |
| `severity` | one of the four levels below |
| `path`     | the HDF5 path it is anchored to (empty for the file itself) |
| `message`  | a human-readable explanation; it may change, so do not parse it |
| `detail`   | supporting values, such as expected and actual |

`conformance.counts` gives the number per severity, and `conformance.conformant`
the verdict.

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
| `REF`    | referential integrity of indices |
| `SEL`, `CONN` | selection and connectivity bounds |
| `EXT`    | extensions of the specification where permitted |
| `FILE`   | the file cannot be read at all |
| `DEV-*`  | registered deviations from the specification |

## Severities

| Severity | Meaning | Effect on the verdict |
|:--|:--|:--|
| `ERROR`     | violates a MUST of the specification | the file is not conformant |
| `DEVIATION` | a registered mismatch between the specification and the writer, pending an editorial decision (`h5md_validator/deviations.py`) | none |
| `WARNING`   | violates a SHOULD, e.g. an unknown unit symbol or an algorithm outside the predefined set | none |
| `INFO`      | informational, e.g. an interaction type outside Table SI that carries a reference | none |

A file is conformant when it has no `ERROR` finding. Each check sets its own
severity when it raises the finding.

## Deviations

A deviation is a registered mismatch between the specification and the
reference implementation, awaiting an editorial decision. Each is kept in
`h5md_validator/deviations.py` with what the implementation does, why it
matters and a suggested resolution. Four are raised as `DEVIATION` findings
when a file shows them; the other six are gaps in the specification text that
no single file can show, and are recorded for the decision only. All are open.

| Code | Mismatch | Raised per file |
|:--|:--|:--|
| `DEV-FF-BONDED-LAYER` | `force_field` has an undocumented `bonded`/`nonbonded` layer | yes |
| `DEV-OBS-GROUP-NESTING` | observables may nest per-selection-group subgroups | yes |
| `DEV-CONSTRAINTS-MODE-UNWRITTEN` | `constraints/mode` is specified but never written | yes |
| `DEV-CONSTRAINTS-ALGORITHM-UNCONDITIONAL` | `constraints/algorithm` is written even when nothing is constrained | yes |
| `DEV-REFCOORD-NO-VOCABULARY` | `reference_coordinate_scaling` has no defined value set | no |
| `DEV-FF-PROVENANCE-UNWRITTEN` | no file records which force field it used | no |
| `DEV-FF-NO-EXTENSION-CONVENTION` | S8 gives no convention for interaction types outside the table | no |
| `DEV-FF-NO-FUNCTIONAL-FORM` | interaction types do not state their functional form | no |
| `DEV-SPEC-GAP-ARITY` | the particle arity of interactions is not specified | no |
| `DEV-NO-ATOM-NAME-CONVENTION` | atom and component names have no canonical convention | no |
