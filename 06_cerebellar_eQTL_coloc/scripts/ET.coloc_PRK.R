library(coloc)
library(dplyr)
library(data.table)
library(sqldf)
library(arrow)
library(parallel)
library(doParallel)

setwd("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/coloc/")

snp_pos <- fread("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/snp_pos.txt.gz", header = T)
snp_pos$chr_post_alleles <- paste0(gsub('chr', '', snp_pos$SNP_id_hg38), ":", snp_pos$other_allele, ":", snp_pos$effect_allele)

head(snp_pos)

loci <- fread("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/PRK_loci_LDlinkR.r2.0.1.EUR.csv", header = T)

sumstats <- fread("~/projects/def-grouleau/cec/ET_cereb_splitseq/raw_gwas_files/GCST009325.tsv")

prepare_eqtl <- function(celltype,chrom_locus,sumstats_locus, vcf){
  eqtl <- read_parquet(paste0("~/scratch/CEREB_eQTL/results/", celltype,".cis_qtl_pairs.", chrom_locus, ".parquet"))
  eqtl.rs <- eqtl %>% left_join(vcf, c('variant_id' = 'ID')) %>%
  left_join(snp_pos, by = c("variant_id" = "chr_post_alleles"), keep = T) %>%  # Join on chr:number format
    mutate(variant_id = ifelse(is.na(SNP), variant_id, SNP)) %>%  # Replace if match found
    left_join(snp_pos, c('variant_id' = 'SNP'), keep = T) %>%
    select(phenotype_id, variant_id, SNP_id_hg19 = SNP_id_hg19.y, start_distance, af, ma_samples, ma_count, pval_nominal, slope, slope_se, eqtl_A2 = REF, eqtl_A1 = ALT)%>%
    filter(SNP_id_hg19%in%sumstats_locus$SNP_id_hg19) %>% 
    add_count(phenotype_id) %>% 
    filter(n>10) #Only keep genes with at least 10 SNPs
	return(eqtl.rs)
}



run_coloc <- function(eqtl_sumstats, celltype, sumstats_locus){
  
	
	if(nrow(eqtl_sumstats)==0){
    return (NULL)
  }


	out <- lapply(unique(eqtl_sumstats$phenotype_id),function(x){
    	message(x)
    		eqtl_sumstats_gene <- filter(eqtl_sumstats, phenotype_id==x)
    		sumstats_locus_gene <- sumstats_locus %>% inner_join(.,eqtl_sumstats_gene,by='SNP_id_hg19')

		sumstats_locus_gene <- sumstats_locus_gene %>%
      			mutate(eqtl_direction=case_when(
      			(effect_allele==eqtl_A1 & other_allele==eqtl_A2)  ~ sign(beta*slope),
      			(effect_allele==eqtl_A2 & other_allele==eqtl_A1)  ~ -sign(beta*slope),
    				TRUE ~ 0))
	


		if (nrow(sumstats_locus_gene)>0){

     	coloc_res_pval <- coloc.abf(
       	dataset1=list(snp=sumstats_locus_gene$variant_id,
       			beta=sumstats_locus_gene$beta,
                     varbeta=sumstats_locus_gene$standard_error^2,
                     type="cc"),
       	dataset2=list(snp=sumstats_locus_gene$variant_id,
       			beta=sumstats_locus_gene$slope,
                     varbeta=sumstats_locus_gene$slope_se^2,
                     sdY=1,
                     type="quant"))

	h4.pp <- as.data.frame(coloc_res_pval$summary, colnames = x)
	h4.pp <- as.data.frame(coloc_res_pval$summary)
	
	top_snp <- coloc_res_pval$results[which.max(coloc_res_pval$results$SNP.PP.H4),]$snp
	gene_info <- filter(sumstats_locus_gene, variant_id == top_snp) %>%
			select(variant_id, GWAS_beta = beta, eqtl_direction, phenotype_id)
	gene_info$h4.pp <- h4.pp[6,]
	gene_info$celltype <- celltype
	return(gene_info)
}

else{
      return (NULL)
    }

})

return(bind_rows(out))

}




coloc_results_all <- mclapply(1:nrow(loci), function(i){

  #Get coordinates from the GWAS locus
  chrom_locus <- loci$chrom[i]
  start <- loci$start[i] %>% as.numeric()
  end <- loci$end[i] %>% as.numeric()
  
  closest_gene_locus <- loci$MAPPED_GENE[i]
  GWAS_snp_name <- loci$GWAS_snp[i]
  GWAS_snp_pos_name <- loci$GWAS_snp_pos[i]	

  #Keep GWAS sumstats of SNPs in the locus
  sumstats_locus <- filter(sumstats,chromosome==gsub("chr", "", chrom_locus)) %>% 
    filter(base_pair_location>=start & base_pair_location<=end) %>% 
    dplyr::mutate(SNP_id_hg19=paste0("chr", chromosome,':',base_pair_location))

	vcf <- fread(paste0("~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/genotyping/HRC/", chrom_locus, "_plink_hg38.dbsnp_151.vcf.gz")) %>%
	select(., ID, REF, ALT)


granule_eqtl <- prepare_eqtl("Granule",chrom_locus,sumstats_locus, vcf)
microglia_eqtl <- prepare_eqtl("Microglia", chrom_locus,sumstats_locus, vcf)
purkinje_eqtl <- prepare_eqtl("Purkinje", chrom_locus,sumstats_locus, vcf)
oligo_eqtl <- prepare_eqtl("Oligodendrocytes", chrom_locus,sumstats_locus, vcf)
opc_eqtl <- prepare_eqtl("OPC", chrom_locus, sumstats_locus, vcf)
bergmann_eqtl <- prepare_eqtl("Bergmann", chrom_locus, sumstats_locus, vcf)
astrocytes_eqtl <- prepare_eqtl("Astrocytes", chrom_locus, sumstats_locus, vcf)
endocytes_eqtl <- prepare_eqtl("Endocytes", chrom_locus, sumstats_locus, vcf)
mli_1_eqtl <- prepare_eqtl("MLI_1", chrom_locus, sumstats_locus, vcf)
mli_2_eqtl <- prepare_eqtl("MLI_2", chrom_locus, sumstats_locus, vcf)
pericytes_eqtl <- prepare_eqtl("Pericytes", chrom_locus, sumstats_locus, vcf)
ubc_eqtl <- prepare_eqtl("UBC", chrom_locus, sumstats_locus, vcf)
golgi_eqtl <- prepare_eqtl("Golgi", chrom_locus, sumstats_locus, vcf)
pli_eqtl <- prepare_eqtl("PLI", chrom_locus, sumstats_locus, vcf)

granule_coloc <- run_coloc(granule_eqtl, "Granule", sumstats_locus)
microglia_coloc <- run_coloc(microglia_eqtl, "Microglia", sumstats_locus)
purkinje_coloc <- run_coloc(purkinje_eqtl, "Purkinje", sumstats_locus)
oligo_coloc <- run_coloc(oligo_eqtl, "Oligodendrocytes", sumstats_locus)
opc_coloc <- run_coloc(opc_eqtl, "OPC", sumstats_locus)
bergmann_coloc <- run_coloc(bergmann_eqtl, "Bergmann", sumstats_locus)
astrocytes_coloc <- run_coloc(astrocytes_eqtl, "Astrocytes", sumstats_locus)
endocytes_coloc <- run_coloc(endocytes_eqtl, "Endocytes", sumstats_locus)
mli_1_coloc <- run_coloc(mli_1_eqtl, "MLI_1", sumstats_locus)
mli_2_coloc <- run_coloc(mli_2_eqtl, "MLI_2", sumstats_locus)
pericytes_coloc <- run_coloc(pericytes_eqtl, "Pericytes", sumstats_locus)
ubc_coloc <- run_coloc(ubc_eqtl, "UBC", sumstats_locus)
golgi_coloc <- run_coloc(golgi_eqtl, "Golgi", sumstats_locus)
pli_coloc <- run_coloc(pli_eqtl, "PLI", sumstats_locus)


merged_eqtl <- rbind(granule_coloc, microglia_coloc, purkinje_coloc,
		oligo_coloc, opc_coloc, bergmann_coloc, astrocytes_coloc,
		endocytes_coloc, mli_1_coloc, mli_2_coloc, pericytes_coloc,
		ubc_coloc, golgi_coloc, pli_coloc) %>%
		mutate(closest_gene=closest_gene_locus,
		lead_GWAS_snp=GWAS_snp_name,
             GWAS_snp_pos=GWAS_snp_pos_name)

return(merged_eqtl)
}, mc.cores = Sys.getenv("SLURM_CPUS_PER_TASK"),mc.preschedule = FALSE) 

coloc_results_all <- bind_rows(coloc_results_all)
write.table(coloc_results_all, "~/projects/def-grouleau/cec/ET_cereb_splitseq/CEREB_eQTL/coloc/PRK_all_celltypes_coloc.txt")
