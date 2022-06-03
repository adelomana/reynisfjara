library(biomaRt)

#
# 0. user-defined variables
#
setwd("~/scratch/")

annotation_file = "/home/adrian/projects/reynisfjara/results/annotation/annotation.csv"

#
# 1. generate gene to transcript mapping
#

# annotation defined from sleuth walkthrough
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'https://www.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", 'description', 'gene_biotype'), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
View(t2g)

write.csv(t2g, annotation_file)