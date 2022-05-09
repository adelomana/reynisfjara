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
# 3. work on different contrasts
#

# 3.1. work on WT contrast
s2c = metadata[metadata$genotype == 'wt', ]
paths = file.path(kallisto_dir, s2c$sample, "abundance.h5")
s2c$path = paths
View(s2c)

# prepare object for sleuth
design <- ~ time
so = sleuth_prep(s2c, 
                 full_model = design,
                 target_mapping=t2g, 
                 aggregation_column='ens_gene', 
                 read_bootstrap_tpm=TRUE, 
                 extra_bootstrap_summary=TRUE,
                 transformation_function = function(x) log2(x + 0.5))
so <- sleuth_fit(so)
models(so)

# build full and partial models ==> 6,008 DEGs
#so = sleuth_fit(so, ~time, 'full')
#so = sleuth_fit(so, ~1, 'reduced')
#so = sleuth_lrt(so, 'reduced', 'full')  ===> can i do the same with wt?

oe <- sleuth_wt(so, which_beta='time', which_model='full')


# tables
sleuth_table_gene <- sleuth_results(so, test='full', test_type = 'wt')

View(sleuth_table_gene)
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
print(dim(sleuth_table_gene))
head(sleuth_table_gene, 20)
write.csv(sleuth_table_gene, file.path(results_dir, paste('WT', '.csv', sep='')))








