# Required data inputs

The scripts expect the following files in this directory. They are not redistributed because they include third-party source data or results derived from access-controlled summary statistics.

## Mendelian randomization

- `mendelian_randomization_nonduplicate_instrument_sets.csv`
- `mendelian_randomization_nominal_associations.csv`
- `mendelian_randomization_focal_estimators.csv`
- `mendelian_randomization_instruments.csv`
- `mendelian_randomization_leave_one_out.csv`
- `mendelian_randomization_all_estimators.csv`
- `mendelian_randomization_sensitivity_statistics.csv`

## Metabolite identity and target prediction

- `metabolite_identity.json`
- `swisstargetprediction_targets.csv`
- `similarity_ensemble_approach_targets.csv`
- `predicted_target_union.csv`
- `predicted_target_annotations.csv`
- `predicted_targets_with_human_evidence.csv`

## Gene association and cell-type evidence

- `magma_predicted_targets.csv`
- `magma_competitive_gene_set_test.csv`
- `carbonic_anhydrase_locus.csv`
- `carbonic_anhydrase_cell_type_cis_eqtl.csv`
- `ca3_colocalization_prior_grid.csv`
- `cerebellar_cell_type_expression.csv`
- `gene_identifier_mapping.csv`

## Purkinje-cell and cerebellar transcriptomic evidence

- `gse197345_raw_count_matrix.csv.gz`
- `gse197345_sample_metadata.csv`
- `gse197345_rlog_pca_coordinates.csv`
- `gse197345_rlog_pca_variance.csv`
- `gse197345_selected_gene_rlog_zscores.csv`
- `gse197345_heatmap_sample_order.csv`
- `purkinje_published_differential_expression.csv`
- `purkinje_differential_expression_comparison.csv`
- `predicted_target_expression_results.csv`
- `mroast_gene_set_results.csv`
- `fry_gene_set_results.csv`
- `prediction_probability_weighted_roast.csv`
- `scp3177_umap_coordinates.csv.gz`

The raw GSE197345 count files can be obtained from the Gene Expression Omnibus. The study-level source and access route for each remaining input are described in the manuscript Data Availability Statement.
