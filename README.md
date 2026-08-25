# Essential tremor catechin-related target prioritization

This repository accompanies the manuscript "Gut Microbial Genetics Prioritize Catechin-Related Cerebellar Targets in Essential Tremor."

The repository contains the code used to prepare the regularized-logarithm principal-component and heatmap inputs for the Purkinje-cell RNA-sequencing analysis and to assemble the five main figures and three supplementary figures. The accompanying figure documentation records the source paper, source figure, borrowed presentation structure, and study-specific data used for every panel.

## Contents

- `code/prepare_rlog_pca.R`: regularized-logarithm transformation, principal-component analysis, and gene-wise z-score preparation for GSE197345.
- `code/generate_figures.py`: generation of publication figures and supplementary tables from the analysis result tables.
- `code/figure_style.py`: shared typography, color, and panel-label settings.
- `figure_documentation/`: panel-level reference mapping, legends, specifications, and color definitions.
- `assets/`: PubChem structure depictions used in the chemical-identity panels.
- `data/README.md`: required inputs and access information.

## Data availability

Microbial genome-wide association summary data were obtained from the MiBioGen consortium. Essential-tremor summary statistics were obtained from deCODE genetics under the provider's access terms. Cerebellar transcriptomic data are available from the Gene Expression Omnibus under accession numbers GSE134878 and GSE197345. Single-cell cerebellar data are available from the Broad Institute Single Cell Portal under study SCP3177.

Third-party source data, access-controlled summary statistics, and derived variant-level tables are not redistributed in this repository. The required filenames and their roles are listed in `data/README.md`.

## Software

Python dependencies are listed in `code/requirements.txt`. The transcriptomic preparation script additionally requires R, DESeq2, and matrixStats.

From the repository root, run:

```text
Rscript code/prepare_rlog_pca.R
python code/generate_figures.py
```

The figure script writes PDF and TIFF files to `output/main_figures/` and `output/supplementary_figures/`.

## Citation

Please cite the accompanying manuscript. Citation metadata are provided in `CITATION.cff`.

## Contact

Feichi Chen  
Department of Geriatrics, the First Affiliated Hospital of Wenzhou Medical University  
chenfeichi@wmu.edu.cn
