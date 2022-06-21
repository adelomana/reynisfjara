#
# -1. packages installation
#
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install("DESeq2")

#
# 0. load libraries
#
library(biomaRt)
library(tximport)
library(DESeq2)
library(BiocParallel)
library(crayon) # so the messages are blue
 
#
# 1. user defined variables
#
register(MulticoreParam(20))

setwd("~/scratch/")

kallisto_dir = "/home/adrian/projects/reynisfjara/results/kallisto/kallisto.100"
metadata_file = "/home/adrian/projects/reynisfjara/metadata/reynisfjara_project_metadata - Sheet1.tsv"
results_dir = '/home/adrian/projects/reynisfjara/results/DEGs_DESeq2_reference_endpoint/'

#
# 2. generate gene to transcript mapping
#

# annotation defined from sleuth walkthrough
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'https://www.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                     "ensembl_gene_id",
                                     "gene_biotype", 
                                     "description", 
                                     "external_gene_name"), 
                      mart = mart)
t2g <- dplyr::rename(t2g, 
                     target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, 
                     ext_gene = external_gene_name)
View(t2g)

#
# 3. read metadata
#
metadata = read.table(metadata_file, sep='\t', header=TRUE)
View(metadata)

#
# 4. hypothesis testing
#
mice = unique(metadata$mouse)

# to be commented out
mouse = mice[1]

# f.1. slice metadata
working_metadata = metadata[metadata$mouse == mouse, ]
working_metadata$time_factor = paste(as.factor(working_metadata$time), 'h', sep='')
View(working_metadata)
  
 