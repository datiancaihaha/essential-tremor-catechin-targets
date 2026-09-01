# Reproducibility archive

## Manuscript

Catechin-Related Carbonic Anhydrase Genes in Essential Tremor: Microbial Genetic and Cerebellar Evidence

## Scope

This release organizes the analysis code, input schemas, captured software environments, and nonrestricted derived results used for the reported microbial genetic, chemical-target, essential-tremor genetic, and cerebellar analyses. It is divided into ten analysis modules that match the manuscript workflow.

The archive does not redistribute access-controlled deCODE essential-tremor summary statistics, controlled source datasets, or third-party files whose terms prohibit redistribution. Each module documents the required input fields and the source from which authorized users can obtain those inputs.

## Modules

1. `01_microbial_MR`: genus-level instruments, harmonization, Mendelian-randomization estimators, and diagnostics.
2. `02_species_analysis`: independent shotgun-metagenomic species-resolution analyses.
3. `03_metabolite_curation`: taxon-to-metabolite evidence and exact chemical identities.
4. `04_target_prediction`: structure-resolved target prediction and experimental-database audit.
5. `05_MAGMA`: target-gene association and competitive gene-set analysis.
6. `06_cerebellar_eQTL_coloc`: cerebellar cell-type cis-eQTL extraction and CA3 colocalization.
7. `07_Purkinje_transcriptomics`: laser-captured Purkinje-cell preprocessing and case-control analysis.
8. `08_single_nucleus_pseudobulk`: donor-level single-nucleus pseudobulk and threshold sensitivity.
9. `09_cell_composition`: propeller, matched-donor, scCODA, and CAMERA analyses.
10. `10_figures`: figure-generation code, documentation, and final PDF figures.

## Reproduction boundary

The included derived results permit verification of the reported tables and figures without redistributing restricted source statistics. Re-execution from raw data requires the authorized source files described in each `input_schema.md` and local path configuration in the archived scripts.

## Release

Version 1.0.0 is the reproducibility archive accompanying the manuscript. The tagged source is available at https://github.com/datiancaihaha/essential-tremor-catechin-targets/releases/tag/v1.0.0.

The permanent Zenodo archive is available at https://doi.org/10.5281/zenodo.22236977.
