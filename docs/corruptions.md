# Corruptions used by the negative tests

Validation 4 applies 65 corruptions, each injecting exactly one defect into a
copy of a valid reference file. The validator must raise the expected finding
code (see [check codes](check-codes.md)) for the corrupted file to count as
detected. A corruption that cannot apply to a given file, for example a
selection defect on a file without selections, is skipped rather than counted.

The list mirrors `CORRUPTIONS` in `h5md_validator/corrupt.py` of the validator in [gromacs/h5md-biophysics](https://gitlab.com/gromacs/h5md-biophysics/-/tree/main/validator), in the same order.

## Module (3)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-MOD-VERSION-UNSUPPORTED` | declare module version 9.9 | `MOD-VERSION-UNSUPPORTED` |
| `C-MOD-VERSION-MISSING` | remove the module version attribute | `MOD-VERSION-MISSING` |
| `C-MOD-MISSING` | remove the biophysics module declaration | `MOD-MISSING` |

## Topology invariant (8)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-TOPO-TOTAL-PARTICLES` | increment total_particle_count by one | `TOPO-INVARIANT-PARTICLES` |
| `C-TOPO-TOTAL-RESIDUES` | increment total_residue_count by one | `TOPO-INVARIANT-RESIDUES` |
| `C-TOPO-BLOCK-COUNT` | increment the molecule count of the first block by one | `TOPO-INVARIANT-PARTICLES` |
| `C-TOPO-PARTICLE-DATA-SHORT` | drop the last row of particle_data of the largest molecule type | `TOPO-PARTICLE-DATA-LENGTH` |
| `C-TOPO-TYPE-ID-SHORT` | drop the last type_id entry of the largest molecule type | `TOPO-TYPE-ID-LENGTH` |
| `C-TOPO-COULOMB14-COUNT` | drop the last coulomb_14 entry of the largest molecule type | `TOPO-COULOMB14-PAIRING` |
| `C-TOPO-PARTICLE-COUNT-ATTR` | increment particle_count of a molecule type by one | `TOPO-PARTICLE-DATA-LENGTH` |
| `C-TOPO-SEQUENCE-SHORT` | drop the last row of sequence | `TOPO-SEQUENCE-LENGTH` |

## Referential integrity (15)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-TOPO-BLOCK-DANGLING` | point a molecule block at a molecule type that does not exist | `REF-BLOCK-DANGLING` |
| `C-TOPO-SEQUENCE-INDEX` | set a sequence_index beyond residue_count | `REF-SEQUENCE-INDEX` |
| `C-TOPO-NAME-INDEX` | set a particle name index beyond particle_name_table | `REF-PARTICLE-NAME-INDEX` |
| `C-TOPO-COMPONENT-NAME-INDEX` | set a component name index beyond component_name_table | `REF-COMPONENT-NAME-INDEX` |
| `C-FF-PARAM-ID-DANGLING` | reference a force-field parameter index that does not exist | `REF-PARAMETER-ID` |
| `C-FF-PARAM-ID-NEGATIVE` | reference a negative force-field parameter index | `REF-PARAMETER-ID` |
| `C-FF-PARAM-TABLE-MISSING` | remove a parameters dataset that an interaction list references | `REF-PARAMETER-TABLE-MISSING` |
| `C-FF-PARAM-TABLE-SHORT` | drop the parameter records from the highest referenced index onward | `REF-PARAMETER-ID` |
| `C-TOPO-TYPE-ID-INDEX` | set a type_id to the number of nonbonded types, one past the matrix | `REF-TYPE-ID` |
| `C-INT-PARTICLE-ID` | set an interaction particle_id beyond particle_count | `REF-PARTICLE-ID` |
| `C-EXCL-ORDER` | store an exclusion pair higher index first | `REF-EXCL-ORDER` |
| `C-EXCL-INDEX` | set an exclusion index beyond particle_count | `REF-EXCL-INDEX` |
| `C-CONN-INDEX` | set a bond particle index beyond the system size | `CONN-PARTICLE-ID` |
| `C-SEL-INDEX` | set a selection index beyond the system size | `REF-SEL-INDEX` |
| `C-SEL-RANGE` | set a selection range end beyond the system size | `REF-SEL-RANGE` |

## Required element (8)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-TOPO-BLOCKS-MISSING` | remove the molecule_blocks dataset | `STRUCT-MISSING-DATASET` |
| `C-TOPO-PARTICLE-DATA-MISSING` | remove particle_data from a molecule type | `STRUCT-MISSING-DATASET` |
| `C-TOPO-NAME-TABLE-MISSING` | remove particle_name_table | `STRUCT-MISSING-DATASET` |
| `C-PARAM-REQUIRED-MISSING` | remove run_control/simulation_type | `STRUCT-MISSING-DATASET` |
| `C-TOPO-SEQUENCE-MISSING` | remove sequence from the largest molecule type | `STRUCT-MISSING-DATASET` |
| `C-INT-PARTICLE-IDS-MISSING` | remove particle_ids from a bonded interaction list | `STRUCT-MISSING-DATASET` |
| `C-TOPO-TYPE-ID-MISSING` | remove the nonbonded type_id of the largest molecule type | `STRUCT-MISSING-DATASET` |
| `C-TOPO-TOTAL-PARTICLES-MISSING` | remove the total_particle_count attribute of topology | `STRUCT-MISSING-ATTR` |

## Physical validity (2)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-TOPO-NEGATIVE-MASS` | give a particle a negative mass | `TOPO-NEGATIVE-MASS` |
| `C-CONN-SELF` | bond a particle to itself | `CONN-SELF-BOND` |

## Structure (14)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-TOPO-FIELD-RENAMED` | rename the charge field of particle_data | `STRUCT-COMPOUND-FIELDS` |
| `C-FF-FIELD-RENAMED` | rename a force-field compound parameter field | `STRUCT-COMPOUND-FIELDS` |
| `C-INT-ARITY` | drop a column from particle_ids | `STRUCT-INTERACTION-ARITY` |
| `C-INT-LENGTH` | make parameter_id and particle_ids differ in length | `STRUCT-INTERACTION-LENGTH` |
| `C-EXCL-DUPLICATE` | list an exclusion pair twice | `REF-EXCL-DUPLICATE` |
| `C-CONN-SHAPE` | store bonds as Integer[N][3] | `STRUCT-SHAPE` |
| `C-TIME-STEP-DTYPE` | store the step dataset of a time-dependent observable as Float | `STRUCT-DTYPE` |
| `C-OBS-SHAPE` | store the pressure tensor as 3x2 | `STRUCT-SHAPE` |
| `C-OBS-UNKNOWN` | add an observable that the specification does not define | `STRUCT-UNDOCUMENTED` |
| `C-OBS-SUBGROUP-NAME` | rename a selection subgroup of observables to a predefined observable | `STRUCT-OBS-SUBGROUP-NAME` |
| `C-PARAM-PER-GROUP` | give reference_temperature one more entry than there are groups | `STRUCT-PER-GROUP-LENGTH` |
| `C-PARAM-PERIODICITY-SHAPE` | store periodicity as a 2-vector | `STRUCT-SHAPE` |
| `C-PARAM-DTYPE` | store the timestep as an Integer | `STRUCT-DTYPE` |
| `C-SEL-EMPTY` | add a selection group that defines nothing | `SEL-EMPTY` |

## Units (6)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-FF-UNITS-MISSING` | remove the units attribute of a compound parameters dataset | `UNITS-MISSING` |
| `C-FF-UNITS-LENGTH` | shorten the units attribute to one entry | `UNITS-LENGTH` |
| `C-FF-UNIT-MALFORMED` | write a unit as 'nm^3' instead of 'nm+3' | `UNIT-SYNTAX` |
| `C-TIME-UNIT-MISSING` | remove the unit attribute of a time dataset | `UNIT-MISSING` |
| `C-OBS-UNIT-MISSING` | remove the unit attribute of an observable value | `UNIT-MISSING` |
| `C-PARAM-UNIT-MISSING` | remove the unit attribute of the timestep | `UNIT-MISSING` |

## Time alignment (4)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-TIME-LENGTH` | truncate the step dataset of a time-dependent observable | `TIME-LENGTH-MISMATCH` |
| `C-TIME-MONOTONIC` | swap two time values so the series is not monotonic | `TIME-NOT-MONOTONIC` |
| `C-TIME-STEP-MISSING` | remove the step dataset of a time-dependent observable | `TIME-STEP-MISSING` |
| `C-TIME-VALUE-MISSING` | remove the value dataset of a time-dependent observable | `TIME-VALUE-MISSING` |

## Conditional structure (3)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-PARAM-CONDITIONAL` | name a thermostat whose parameter subgroup is absent | `COND-SUBGROUP-MISSING` |
| `C-PARAM-COND-ORPHAN` | add a subgroup for a simulation type that is not selected | `COND-SUBGROUP-ORPHAN` |
| `C-PARAM-SD-THERMOSTAT` | add temperature_coupling to a stochastic dynamics run | `COND-SD-TEMPERATURE-COUPLING` |

## Cross-group consistency (2)

| Corruption | Defect injected | Expected finding |
|:--|:--|:--|
| `C-TRAJ-MISMATCH` | drop a particle from the trajectory but not from the topology | `TOPO-TRAJECTORY-MISMATCH` |
| `C-SEL-REFERENCE-MISSING` | name a parameter group that is not defined in selections | `REF-SELECTION-MISSING` |
