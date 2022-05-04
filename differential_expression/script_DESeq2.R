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

### annotation defined by ALO
# working_atributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name")
# ensembl96 = useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=96)
# t2g = getBM(attributes=working_atributes, mart=ensembl96)
# t2g = dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
# dim(t2g)
# View(t2g)

#
# 2. read metadata
#
metadata = read.table(metadata_file, sep='\t', header=TRUE)
View(metadata)

#
# 3. work on a full model
#

# define the metadata table
s2c = metadata
paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
s2c$path = paths
View(s2c)

# prepare object for sleuth
so = sleuth_prep(s2c, target_mapping=t2g, aggregation_column='ens_gene')

# build full and partial models
so <- sleuth_fit(so, ~genotype, 'reduced')
so <- sleuth_fit(so, ~genotype + time, 'full')
so <- sleuth_lrt(so, 'reduced', 'full')

# tables
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
print(dim(sleuth_table_gene))
head(sleuth_table_gene, 20)
write.csv(sleuth_table_gene, file.path(results_dir, paste('complex', '.csv', sep='')))

#
# 4. work with smaller contrasts
#
mice = unique(metadata$mouse)
induction_times = c(48, 72)

for (mice_index in 1:length(mice)) {
  for (time_index in 1:length(induction_times)) {
    
    working_mice = mice[mice_index]
    working_time = induction_times[time_index]

    s2c = metadata[which(metadata$mouse == working_mice & (metadata$time == 0 | metadata$time == working_time)), ]
    paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
    s2c$path = paths
    View(s2c)
    
    so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)
    so <- sleuth_fit(so, ~time, 'full')
    so <- sleuth_fit(so, ~1, 'reduced')
    so <- sleuth_lrt(so, 'reduced', 'full')
    
    sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
    sleuth_table <- dplyr::filter(sleuth_table, qval <= 0.05)
    print(dim(sleuth_table))
    head(sleuth_table, 20)
    write.csv(sleuth_table, file.path(results_dir, paste(paste(working_mice, working_time, sep='.'), '.csv', sep='')))
    
  }
}


















