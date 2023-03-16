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
# 3.1. testing time + genotype + time:genotyppe as in 
# https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html
#
cat(blue('comparing WT with time'), fill=TRUE)

# 3.1.1. prepare metadata including paths
s2c = metadata
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
so = sleuth_fit(so, ~time_factor + genotype + time_factor:genotype, 'full')
so = sleuth_fit(so, ~time_factor + genotype, 'reduced')

# 3.1.5. LRT test: should be the default approach from what sleuth developers think
cat(blue('LRT testing...', fill=TRUE))
lrt = sleuth_lrt(so, 'reduced', 'full')
lrt_table = sleuth_results(lrt, 'reduced:full', 'lrt', show_all=FALSE, pval_aggregate=TRUE)
lrt_table = dplyr::filter(lrt_table, qval <= 0.05)
print(dim(lrt_table))

# 3.1.6. store into files
cat(blue('storing...'), fill=TRUE)
write.csv(lrt_table, file.path(results_dir, paste('interactions', 'LRT', 'csv', sep='.')))
