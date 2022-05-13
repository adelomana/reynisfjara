#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install("tictoc")

library(biomaRt)
library(sleuth)
library(crayon) # so the messages are blue
library(tictoc)
  
#
# 0. user-defined variables
#
setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/reynisfjara/results/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/reynisfjara/metadata/reynisfjara_project_metadata\ -\ Sheet1.tsv"
results_dir = '/home/adrian/projects/reynisfjara/results/tpm'

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
# 3. create a sleuth object
#
s2c = metadata
paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
s2c$path = paths
print(s2c)

# prepare object for sleuth
cat(blue('preparing sleuth object...'), fill=TRUE)
tic()
so = sleuth_prep(s2c, 
                 target_mapping=t2g, 
                 aggregation_column='ens_gene', 
                 read_bootstrap_tpm=TRUE, 
                 extra_bootstrap_summary=TRUE, 
                 gene_mode=TRUE)
toc()

#
# 4. store TPMs
#
cat(blue('storing'), fill=TRUE)
tpm_table = sleuth_to_matrix(so, 'obs_norm', 'tpm')
View(tpm_table)

write.csv(tpm_table, file.path(results_dir, 'sleuth_TPM_gene.csv'))