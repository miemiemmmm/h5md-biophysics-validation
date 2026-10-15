# Timing data

## Benchmark files

Ten repeats of the GROMACS benchmark
`H5mdTopologyBenchmarkTest.BenchmarkDecomposeWritingProcess`
(`src/gromacs/fileio/h5md/tests/h5md_benchmark.cpp`), one pair of files per
repeat:

- `decompose_runN.json`: Google Test JSON output. Each of the 58 test cases is
  one system (`value_param` names its `.tpr`) and records the elapsed time of
  each phase of the metadata write, in µs, as `metadata_elapsed_us`,
  `particle_block_elapsed_us`, `topology_elapsed_us`,
  `bond_connectivity_elapsed_us`, `force_field_elapsed_us`,
  `md_parameters_elapsed_us` and `selection_groups_elapsed_us`, plus `natoms`.
- `decompose_runN.log`: the console output of the same run.

The input is the `nvt-vrescale` run input of each of the 58 systems, as written
by `scripts/build_reference_files.py` to `_data/tprs/`; the benchmark reads them
from a directory fixed in its source (`kTprBaseDir`).

## Generating the benchmark data

The benchmark is not part of GROMACS; `h5md_benchmark.patch` adds it. Apply it
to the GROMACS source shipped with the paper, set `kTprBaseDir` in
`h5md_benchmark.cpp` to your `_data/tprs` directory, and build the `h5md-test`
target. Besides adding the test, the patch comments out one `private:` in
`h5md.h` so that the benchmark can time the writer's internal phases.

```bash
git apply /path/to/validation/data/h5md_benchmark.patch   # in the GROMACS source tree
```

The runs used a Release build on one AMD Ryzen 9 7950X workstation:

```bash
h5md-test --gtest_filter='BenchmarkTprs/H5mdTopologyBenchmarkTest.BenchmarkDecomposeWritingProcess/*' \
          --gtest_output=json:decompose_run1.json > decompose_run1.log
```

`scripts/validation5_cost_figure.py` reads every `decompose_run*.json`, takes the
median of each phase over the repeats, and plots each phase's share of the write
time as panel C.

The patch is not maintained. It is kept only to show how these data were
generated.
