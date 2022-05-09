#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# setRepositories(ind=c(1:6))
# 
# BiocManager::install("biomaRt")
# 
# # install sleuth
# BiocManager::install("rhdf5")
# BiocManager::install("devtools")
# devtools::install_github("pachterlab/sleuth")

library(biomaRt)
library(sleuth)

#
# 0. user-defined variables
#
setwd("~/scratch/")
kallisto_dir = "/home/adrian/projects/reynisfjara/results/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/reynisfjara/metadata/reynisfjara_project_metadata\ -\ Sheet1.tsv"
results_dir = '/home/adrian/projects/reynisfjara/results/DEGs_sleuth'

kallisto_dir = "/Users/adrian/gd15/tmp/kallisto.100"
metadata_file = "/Users/adrian/gd15/hi/research/reynisfjara/metadata/reynisfjara_project_metadata\ -\ Sheet1.tsv"
results_dir = '/Users/adrian/gd15/hi/research/reynisfjara/results/sleuth'

#
# 1. generate gene to transcript mapping
#

# annotation defined from sleuth walkthrough
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'https://www.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
View(t2g)

#
# 2. read metadata
#
metadata = read.table(metadata_file, sep='\t', header=TRUE)
View(metadata)

#
# 3. work on different contrasts
#
comparisons = unique(metadata$mouse)

for (comparison in comparisons){
  print(comparison)

  # 3.1. work on WT contrast
  s2c = metadata[metadata$mouse == comparison, ]
  paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
  s2c$path = paths
  print(s2c)
  
  # prepare object for sleuth
  print('preparing sleuth object...')
  so = sleuth_prep(s2c, 
                   target_mapping=t2g, 
                   aggregation_column='ens_gene', 
                   read_bootstrap_tpm=TRUE, 
                   extra_bootstrap_summary=TRUE)
  
  # build full and partial models 
  print('building models...')
  so = sleuth_fit(so, ~time, 'full')
  so = sleuth_fit(so, ~1, 'reduced')
  
  # Wald test
  print('Wald testing...')
  wald = sleuth_wt(so, which_beta='time')
  wald_table = sleuth_results(wald, test='time', test_type='wt', show_all=FALSE, pval_aggregate=TRUE)
  wald_table = dplyr::filter(wald_table, qval <= 0.05)
  print(dim(wald_table))
  
  # LRT test
  print('LRT testing...')
  lrt = sleuth_lrt(so, 'reduced', 'full')
  lrt_table = sleuth_results(lrt, 'reduced:full', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
  lrt_table = dplyr::filter(lrt_table, qval <= 0.05)
  print(dim(lrt_table))
  
  # store into files
  print('storing...')
  write.csv(wald_table, file.path(results_dir, paste(comparison, 'wald', 'csv', sep='.')))
  write.csv(lrt_table, file.path(results_dir, paste(comparison, 'LRT', 'csv', sep='.')))

}



