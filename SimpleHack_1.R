#installed GenomicDataCommons and DESeq2
#extract stable gene IDs from dataframe
library(GenomicDataCommons)
names <- DDR_epi_genes_names$Gene.stable.ID
#print(names)
raw_counts <- readHTSeqFile('c98b8334-836a-4c3c-b149-b08a7256aa68.htseq.counts.gz', sample = 'counts')

#extract only raw counts of genes in list
counts <- c()

for (i in 1:length(names)) {
  counts[i] <- raw_counts$counts[sapply(names[i], grepl, raw_counts$feature, ignore.case=TRUE)]
  counts_genes <- cbind(names, counts)
}

counts_genes <- as.data.frame(counts_genes)

#convert counts to log counts
counts_genes$log_counts <- log10(as.numeric(counts_genes$counts)+1)

#normalise log raw counts with DESeq2

