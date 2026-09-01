#Load
library(rtracklayer)
library(edgeR)
library(dplyr)
library(tibble)
library(foreach)
library(parallel)
library(doParallel)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)


registerDoParallel(cores=8)


###MAKE PHENOTYPE FILE
#load pseudobulk expression file
pb <- readRDS("~/projects/def-grouleau/cec/ET_cereb_splitseq/RDS_files/pseudo_mega/pseudo_bulk_lit_ALLcell_types.RData")

gtf <- rtracklayer::import('Homo_sapiens.GRCh38.105.gtf') %>%
  as.data.frame() %>%
  dplyr::filter(type=='gene') %>%
  mutate(TSS_start=ifelse(strand=='+',start,end),
         TSS_end=ifelse(strand=='+',start+1,end+1)) %>%
  mutate(gene=paste0(gene_name,'_',gene_id)) %>%
  mutate(gene_name=ifelse(gene_name %in% NA, gene_id, gene_name)) %>%
  dplyr::select(seqnames,TSS_start,TSS_end, gene_name) %>%
  filter(seqnames%in%c(1:22))
gtf$seqnames <- droplevels(gtf$seqnames)

#LOAD GENOTYPING INFO
#geno.samples <- read_tsv("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/genotyping/HRC/rename_chr1_plink_hg38.fam", col_names = F)
geno.samples <- read_tsv("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/genotyping/plink_files_HRC/chr1.fam", col_names = F)
geno.samples <- as.data.frame(geno.samples)
geno.samples_split <- as.data.frame(str_split_fixed(geno.samples$X2, "#", 3))
geno.samples_split$geno.id <- geno.samples$X2
colnames(geno.samples_split) <- c('V1', 'rou.id', 'V3', 'geno.id')
geno.samples_split <- select(geno.samples_split, 'rou.id', 'geno.id')


foreach (i=names(pb))  %do%{
       #extract counts for each cell-type
	counts <- pb[[paste0(i)]] %>% as.matrix() %>% t(.) %>% as.data.frame()
	
	#Rename with genotyping IDs and drop samples that didn't pass genoQC
	counts <- counts %>% rownames_to_column(., 'rou.id') %>% inner_join(geno.samples_split, by = 'rou.id')
	rownames(counts) <- counts$geno.id
	counts <- counts %>% dplyr::select(., -c(rou.id, geno.id)) %>% t(.)

	#Normalize with TMM from edgeR and the inverse normal transform
	#Filter for low count genes
	#counts.filtered <- counts[rowSums(counts > 1) >= 103*0.2,]
	#dge <- edgeR::DGEList(counts = counts.filtered)
	#counts.filtered.cpm <- counts.filtered[rowMeans(cpm(dge)) > 2,]
	#dge <- edgeR::DGEList(counts = counts.filtered.cpm)
	#dge <- calcNormFactors(dge, method = "TMM")
	#CPM <- cpm(dge, log = F)
	#keep <- rownames(CPM[rowMeans(CPM) > 6, ])
	#CPM <- CPM[rownames(CPM) %in% keep,]

	#Bryois filtering
	dge <- edgeR::DGEList(counts = counts)
	dge <- calcNormFactors(dge, method = "TMM")
	CPM <- cpm(dge, log = F)
	keep <- rownames(CPM[rowMeans(CPM) > 6, ])
	CPM <- CPM[rownames(CPM) %in% keep,]


	#Perform TMM on filtered counts
	#dge <- edgeR::DGEList(counts = counts_filtered)
	#dge <- calcNormFactors(dge, method = "TMM")
	#CPM <- cpm(dge, log = F)
	#keep <- rownames(CPM[rowMeans(CPM) > 6, ])
	#CPM <- CPM[rownames(CPM) %in% keep,]
	
	##Inverse normal transform
	n <- ncol(CPM)
	zvalues <- qnorm(ppoints(n))
	z <- CPM
	for (x in 1:nrow(z)) z[x,] <- zvalues[order(order(z[x,]))]
	bed <- z %>% as.data.frame(.) %>% 
        rownames_to_column('gene_name')


	#Quantile norm
	#CPM <- t(apply(CPM, 1, rank, ties.method = "average"))
	#bed <- qnorm(CPM / (ncol(CPM) + 1)) %>%
	#as.data.frame(.) %>% 
        #rownames_to_column('gene_name')

	#Merge with gtf and remove samples that didn't pass genotyping QC
	bed <- bed %>% left_join(., gtf, by = 'gene_name') %>% 
	select(seqnames,TSS_start,TSS_end,gene_name, everything()) %>% 
	arrange(seqnames,TSS_start) %>%
	drop_na(seqnames) %>%
	distinct(gene_name, .keep_all = T)
	
	bed$seqnames <- paste0("chr", bed$seqnames)
	colnames(bed)[c(1,2,3,4)] <- c('#chr','start','end','phenotype_id')
	

	write_tsv(bed, paste0("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/pseudo_per_celltype/", i, ".bed"))
}


###MAKE COVARIATES FILE
geno.pca <- read_delim("~/projects/def-grouleau/cec/ET_cereb_splitseq/data_for_CEC/CEC_FINAL_QC_PCs.eigenvec", col_names = F)

#keep only 3 PCs
geno.pca$X2 <- paste0(geno.pca$X1, "_", geno.pca$X2)
geno.pca <- as.data.frame(geno.pca[,2:5])

colnames(geno.pca) <- c("geno.id", 'gPC1', 'gPC2', 'gPC3')

#keep case control status
metadata <- read.table("~/projects/def-grouleau/cec/ET_cereb_splitseq/RDS_files/cereb.metadata.2023_08_17.txt", header = T)
condition <- select(metadata, ROU.id, condition, sex, age, seq_batch)

cov_common <- left_join(geno.samples_split, condition, by = c("rou.id" = "ROU.id"), multiple = 'first')
cov_common <- left_join(cov_common, geno.pca, by = 'geno.id')


#Wrrite covariate files
foreach (i=names(pb))  %do%{
	#read bed file for each cell_type
	bed <- as.data.frame(read_tsv(paste0("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/pseudo_per_celltype/", i ,".bed")))	
	bed <- bed[-c(1,2,3,4)] %>% t(.) 
	pca <- bed %>% prcomp(.,scale.=T)
	
	

	pca <- as.data.frame(pca$x[,1:30])
	pca <- rownames_to_column(pca, 'geno.id')
	

	cov_ct <- cov_common %>% right_join(., pca, by = 'geno.id')
	cov_ct <- select(cov_ct, -c(rou.id))
	cov_ct <- cov_ct[match(pca$geno.id, cov_ct$geno.id),]
	cov_ct$condition <- ifelse(cov_ct$condition=="case", 1, 2)
	cov_ct$sex <- ifelse(cov_ct$sex=="F", 1, 2)
	cov_ct$seq_batch <- ifelse(cov_ct$seq_batch=="WT", 1, 2)
	cov_ct$age <- scale(cov_ct$age, center = T) %>% as.numeric()
	

	write_tsv(cov_ct, paste0("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/covariates/", i, "_cov.txt"))
}





	



