
# 06_cerebellar_eQTL_coloc

## Purpose

Cell-type-resolved cerebellar cis-eQTL extraction for carbonic anhydrase genes and prior-sensitivity colocalization of CA3 with ET association statistics.

## Contents

- `scripts/`: archived analysis or preparation scripts.
- `derived_results/`: nonrestricted outputs used to verify manuscript results and figures.
- `input_schema.md`: required fields and access boundary for source inputs.
- `software_versions.txt`: captured or documented software information.
- `random_seed.txt`: stochasticity and seed status.
- `session_info/`: captured runtime/session records; when no historical session record exists, this directory contains an explicit note.

## Use

Review `input_schema.md` before re-execution. Archived scripts preserve the analysis implementation and may require local path configuration for authorized source files. Do not substitute a different ET outcome dataset or redistribute access-controlled inputs under this release.

## Restricted-data boundary

The variant-level ET/eQTL colocalization input and variant posterior table are not redistributed because they contain or encode access-controlled ET summary statistics; aggregate colocalization results and input counts are included.
