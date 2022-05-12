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
results_dir = '/home/adrian/projects/reynisfjara/results/DEGs_DESeq2/'

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
mouse = mice[1]

for (mouse in mice){
  cat(blue(mouse))
  
  # f.1. slice metadata
  working_metadata = metadata[metadata$mouse == mouse, ]
  working_metadata$time_factor = paste(as.factor(working_metadata$time), 'h', sep='')
  View(working_metadata)
  
  # f.2. define files
  files = file.path(kallisto_dir, working_metadata$sample, "abundance.h5")
  cat(blue('files'), fill=TRUE)
  print(files)
  
  # f.3. read files
  txi = tximport(files, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)
  dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~time_factor) 
  cat(blue(paste('original dimensions', dim(dds)[1], dim(dds)[2])), fill=TRUE)
  
  # f.4. preliminary filter on poorly detected genes
  threshold = 10
  keep = rowMaxs(counts(dds)) >= threshold
  dds = dds[keep, ]
  cat(blue(paste('dimensions of filtered transcripts', dim(dds)[1], dim(dds)[2])), fill=TRUE)
  
  # f.5. run the model over variable time, actually "time_factor"
  dds = DESeq(dds, test="LRT", reduced=~1)
  
  # f.6. store original results
  res = results(dds, parallel=TRUE) 
  filt1 = res[which(res$pvalue < 0.05), ]
  filt2 = filt1[which(filt1$padj < 0.1), ]
  cat(blue(paste('DEGs found', dim(filt2)[1], sep=' ')), fill=TRUE)
  write.table(filt2, file=paste(results_dir, 'unformatted/unformatted_results_', mouse, '.tsv', sep=''), quote=FALSE, sep='\t')
  
  # store formatted results
  df = as.data.frame(filt2)
  df['common'] = rownames(df)
  
  annotation_table = t2g
  annotation_table['common'] = annotation_table$ens_gene
  annotation_table_unique = annotation_table[!duplicated(annotation_table$common), ]
  dn = merge(df, annotation_table_unique, by='common')
  # check for not missing DEGs because of annotation
  if (dim(df)[1] != dim(dn)[1]){
    print('ERROR: DEG missed on annotation step')
    stop()
  }
  
  up = dn[dn$log2FoldChange > 1, ]
  down = dn[dn$log2FoldChange < -1, ]
  
  cat(blue(paste('DEGs after log2FC > 1 are ', dim(up)[1], sep='')), fill=TRUE)
  cat(blue(paste('DEGs after log2FC < -1 are ', dim(down)[1], sep='')), fill=TRUE)
  
  sorted_up = up[order(up$log2FoldChange, decreasing=TRUE), ]
  sorted_down = down[order(down$log2FoldChange), ]
  
  store = paste(results_dir, mouse, '_up', '.tsv', sep='')
  write.table(sorted_up, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  store = paste(results_dir, mouse, '_down', '.tsv', sep='')
  write.table(sorted_down, file=store, quote=FALSE, sep='\t', row.names=FALSE)
  
  cat(blue('---'), fill=TRUE)
  
}


