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
results_dir = '/home/adrian/projects/reynisfjara/results/DEGs_DESeq2/start_end'

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
mice = unique(metadata$mouse)
View(metadata)

#
# 4. hypothesis testing
#

#
# 4.2. testing that WT t = 0 vs WT t = 72 are the same
#

# 4.2.2. slice metadata
rules = metadata$genotype == 'wt' & (metadata$time == 0 | metadata$time == 72)
working_metadata = metadata[rules, ]
working_metadata$time_factor = paste(as.factor(working_metadata$time), 'h', sep='')
print(working_metadata)

# 4.2.2. define files
files = file.path(kallisto_dir, working_metadata$sample, "abundance.h5")
cat(blue('files'), fill=TRUE)
print(files)

# 4.2.3. read files
txi = tximport(files, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)

# 4.2.4. define DESeq2 object
dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~time_factor)
dds$time_factor = relevel(dds$time_factor, ref="0h")
cat(blue(paste('original dimensions', dim(dds)[1], dim(dds)[2])), fill=TRUE)

# 4.2.5. preliminary filter on poorly detected genes
threshold = 10
keep = rowMaxs(counts(dds)) >= threshold
dds = dds[keep, ]
cat(blue(paste('dimensions of filtered transcripts', dim(dds)[1], dim(dds)[2])), fill=TRUE)

# 4.2.6. run the model
dds = DESeq(dds, test="LRT", reduced=~1)

# 4.2.7. retrieve results and filter
res = results(dds, parallel=TRUE) 
filt1 = res[which(res$pvalue < 0.05), ]
filt2 = filt1[which(filt1$padj < 0.1), ]
filt3 = filt2[which(abs(filt2$log2FoldChange) > 1), ]

# 4.2.8. print info and save table
print(dim(res))
print(dim(filt1))
print(dim(filt2))
print(dim(filt3))
cat(blue(paste('DEGs found', dim(filt3)[1], sep=' ')), fill=TRUE)
write.table(filt3, file=paste(results_dir, '/start_end_WTt72_over_WTt0', '.tsv', sep=''), quote=FALSE, sep='\t')

#
# 4.4. testing that MUT t = 72 vs MUT t = 0 are different
#
for (mouse in mice[2:4]){ # 2:4 refers to the three mutants
  
  # 4.4.1. slice metadata
  rules = metadata$mouse == mouse & (metadata$time == 0 | metadata$time == 72)
  working_metadata = metadata[rules, ]
  working_metadata$time_factor = paste(as.factor(working_metadata$time), 'h', sep='')
  print(working_metadata)
  
  # 4.4.2. define working files
  files = file.path(kallisto_dir, working_metadata$sample, "abundance.h5")
  cat(blue('files'), fill=TRUE)
  print(files)
  
  # 4.4.3. read files
  txi = tximport(files, type="kallisto", tx2gene=t2g, ignoreTxVersion=TRUE)
  
  # 4.4.4. define DESeq2 object
  dds = DESeqDataSetFromTximport(txi, colData=working_metadata, design=~time_factor) 
  dds$time_factor = relevel(dds$time_factor, ref="0h")
  cat(blue(paste('original dimensions', dim(dds)[1], dim(dds)[2])), fill=TRUE)
  
  # 4.4.5. preliminary filter on poorly detected genes
  threshold = 10
  keep = rowMaxs(counts(dds)) >= threshold
  dds = dds[keep, ]
  cat(blue(paste('dimensions of filtered transcripts', dim(dds)[1], dim(dds)[2])), fill=TRUE)
  
  # 4.4.6. run the model
  dds = DESeq(dds, test="LRT", reduced=~1)
  
  # 4.4.7. retrieve results and filter
  res = results(dds, parallel=TRUE) 
  filt1 = res[which(res$pvalue < 0.05), ]
  filt2 = filt1[which(filt1$padj < 0.1), ]
  filt3 = filt2[which(abs(filt2$log2FoldChange) > 1), ]
  
  # 4.4.8. print info and save table
  print(dim(res))
  print(dim(filt1))
  print(dim(filt2))
  print(dim(filt3))
  cat(blue(paste('DEGs found', dim(filt3)[1], sep=' ')), fill=TRUE)
  write.table(filt3, file=paste(results_dir, '/start_end_MUTt72_over_MUTt0_', mouse, '.tsv', sep=''), quote=FALSE, sep='\t')
}



