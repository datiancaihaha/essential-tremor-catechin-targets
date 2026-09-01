
# 01_microbial_MR

## Purpose

Genus-level instrument selection, harmonization, two-sample Mendelian randomization, sensitivity analyses, and nonduplicate taxonomic summaries.

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

Variant-level ET effect estimates are access-controlled and are not included. The archived instrument table retains exposure-side fields and harmonization status only.
