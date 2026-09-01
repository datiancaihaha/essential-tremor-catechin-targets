import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
import torch
import gc
import glob
import sys
import os

gc.collect()
torch.cuda.empty_cache()

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas {pd.__version__}")

cell_type_id = sys.argv[1]

#LOAD GENOTYPES
#plink_prefix_path = "/lustre07/scratch/cec/TOPMed_imputation/chr*"
plink_prefix_path = "/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/genotyping/plink_files_HRC/chr*"
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

#LOAD PHENOTYPES
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(f"/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/pseudo_per_celltype/{cell_type_id}.bed")

#LOAD COVARIATES
covariates_df = pd.read_csv(f"/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/covariates/{cell_type_id}_cov.txt", sep='\t',index_col = 0)


#EMPIRICAL cis-eQTLs (gene-wide)
cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, maf_threshold = 0.05, nperm = 10000)
tensorqtl.calculate_qvalues(cis_df, fdr=0.05)#, qvalue_lambda=0.85)

cis_df = cis_df.sort_values('qval')

cis_df.to_csv(f"/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/results/{cell_type_id}_empirical_cis.csv")


#INDEPENDENT eQTL
indep_df = cis.map_independent(genotype_df, variant_df, cis_df,
                               phenotype_df, phenotype_pos_df, covariates_df)

cis_df.to_csv(f"/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/results/{cell_type_id}_ind_cis.csv")

#NOMINAL cis-eQTLs (all variants)
prefix = cell_type_id
cis_df_nom = cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df,
               prefix, covariates_df, maf_threshold = 0.05, output_dir='/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/results')


#merge parquet files together
files = glob.glob(f'/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/results/{cell_type_id}.cis_qtl_pairs.chr*.parquet')
data = [pd.read_parquet(f,engine='fastparquet') for f in files]
merged_data = pd.concat(data,ignore_index=True)
merged_data.to_csv(f"/lustre06/project/6001220/cec/ET_cereb_splitseq/CEREB_eQTL/results/{cell_type_id}_merged_nominal.txt")






