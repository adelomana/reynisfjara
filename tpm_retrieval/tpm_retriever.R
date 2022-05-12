#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install("DESeq2")
# BiocManager::install("tximport")

#
# 0. load libraries
#
library(DESeq2)
library(tximport)
library(biomaRt)

#
# 1. user-defined variables
#
kallisto_dir = "/home/adrian/projects/reynisfjara/results/kallisto/kallisto.100"
results_dir = '/home/adrian/projects/reynisfjara/results/tpm/'

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
# 3. read files
#
dirnames = list.dirs(kallisto_dir, full.names=TRUE, recursive=FALSE)
dirnames

files = file.path(dirnames, 'abundance.h5')
files

labels = sapply(strsplit(files, split='/',fixed=TRUE), function(x) (x[9]))
labels

txi = tximport(files, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

#
# 4. find abundance
#
tpm = txi$abundance
colnames(tpm) = labels
dim(tpm)
View(tpm)

#
# 5. store
#
store = paste(results_dir, 'DESeq2_TPM_values.tsv', sep='')
write.table(tpm, file=store, quote=FALSE, sep='\t', col.names=NA)