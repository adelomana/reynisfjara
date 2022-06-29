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
library(crayon) # so the messages are blue

#
# 0. user-defined variables
#
setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/reynisfjara/results/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/reynisfjara/metadata/reynisfjara_project_metadata\ -\ Sheet1.tsv"
results_dir = '/home/adrian/projects/reynisfjara/results/DEGs_sleuth'

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
mice = unique(metadata$mouse)
View(metadata)

#
# 3. work on different contrasts
#

#
# 3.1. testing that WT t = 0 vs WT t = 72 are the same 
# 
cat(blue('comparing WT with time'), fill=TRUE)

# 3.1.1. prepare metadata including paths
rules = metadata$genotype == 'wt'
s2c = metadata[rules, ]
s2c$time_factor = paste(as.factor(s2c$time), 'h', sep='')

paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
s2c$path = paths
print(s2c)

# 3.1.2. prepare object for sleuth
cat(blue('preparing sleuth object...'), fill=TRUE)
so = sleuth_prep(s2c, 
                 target_mapping=t2g, 
                 aggregation_column='ens_gene', 
                 read_bootstrap_tpm=TRUE, 
                 extra_bootstrap_summary=TRUE)

# 3.1.3. build full and partial models 
cat(blue('building models...'), fill=TRUE)
so = sleuth_fit(so, ~time, 'full')
so = sleuth_fit(so, ~1, 'reduced')

# 3.1.4. Wald test: Harold Pimentel comments that this may give too many false positives: https://www.biostars.org/p/226022/
cat(blue('Wald testing...'), fill=TRUE)
wald = sleuth_wt(so, which_beta='time')
wald_table = sleuth_results(wald, test='time', test_type='wt', show_all=FALSE, pval_aggregate=TRUE)
wald_table = dplyr::filter(wald_table, qval <= 0.05)
print(dim(wald_table))

# 3.1.5. LRT test: should be the default approach from what sleuth developers think
cat(blue('LRT testing...', fill=TRUE))
lrt = sleuth_lrt(so, 'reduced', 'full')
lrt_table = sleuth_results(lrt, 'reduced:full', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
lrt_table = dplyr::filter(lrt_table, qval <= 0.05)
print(dim(lrt_table))

# 3.1.6. store into files
cat(blue('storing...'), fill=TRUE)
write.csv(wald_table, file.path(results_dir, paste('WT.t72overt0', 'wald', 'csv', sep='.')))
write.csv(lrt_table, file.path(results_dir, paste('WT.t72overt0', 'LRT', 'csv', sep='.')))

#
# 3.2. testing that MUT t = 72 vs MUT t = 0 are different
# 
for (mouse in mice[2:4]){ # 2:4 refers to the three mutants
  cat(blue(mouse), fill=TRUE)
  
  # 3.2.1. prepare metadata including paths
  rules = metadata$mouse == mouse
  s2c = metadata[rules, ]
  s2c$time_factor = paste(as.factor(s2c$time), 'h', sep='')
  
  paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
  s2c$path = paths
  print(s2c)
  
  # 3.2.2. prepare object for sleuth
  cat(blue('preparing sleuth object...'), fill=TRUE)
  so = sleuth_prep(s2c, 
                   target_mapping=t2g, 
                   aggregation_column='ens_gene', 
                   read_bootstrap_tpm=TRUE, 
                   extra_bootstrap_summary=TRUE)
  
  # 3.2.3. build full and partial models 
  cat(blue('building models...'), fill=TRUE)
  so = sleuth_fit(so, ~time, 'full')
  so = sleuth_fit(so, ~1, 'reduced')
  
  # 3.2.4. Wald test
  cat(blue('Wald testing...'), fill=TRUE)
  wald = sleuth_wt(so, which_beta='time')
  wald_table = sleuth_results(wald, test='time', test_type='wt', show_all=FALSE, pval_aggregate=TRUE)
  wald_table = dplyr::filter(wald_table, qval <= 0.05)
  print(dim(wald_table))
  
  # 3.2.5. LRT test
  cat(blue('LRT testing...'), fill=TRUE)
  lrt = sleuth_lrt(so, 'reduced', 'full')
  lrt_table = sleuth_results(lrt, 'reduced:full', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
  lrt_table = dplyr::filter(lrt_table, qval <= 0.05)
  print(dim(lrt_table))
  
  # 3.2.6. store into files
  cat(blue('storing...'), fill=TRUE)
  write.csv(wald_table, file.path(results_dir, paste(mouse, 't72overt0', 'wald', 'csv', sep='.')))
  write.csv(lrt_table, file.path(results_dir, paste(mouse, 't72overt0', 'LRT', 'csv', sep='.')))
}
