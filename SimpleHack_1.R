#installed GenomicDataCommons and DESeq2
#extract stable gene IDs from dataframe
library(GenomicDataCommons)
library(DESeq2)

#load curated gene list
DDR_epi_genes_names <- readRDS("DDR_epi_genes_names.rds")

#add Menin to gene list
Men1 <- data.frame("ENSG00000133895","MEN1","tumor suppressor associated with multiple endocrine neoplasia type 1", "GeneCards")
rownames(Men1) <- Men1[,1]
names(Men1) <- c("Gene.stable.ID","Gene.name","Gene.description","Sourcelist")
DDR_epi_genes_names <- rbind(DDR_epi_genes_names, Men1)

#load cancer datasets from TCGA website
#tcga biolinks package?
list <- c('c98b8334-836a-4c3c-b149-b08a7256aa68.htseq.counts.gz','21feac5d-641f-4802-a910-51b73431de3c.htseq.counts.gz')

#extract counts for genes of interest from datasets
counts_genes <- ensmbl <- DDR_epi_genes_names$Gene.stable.ID

for (i in 1:length(list)){
  raw_counts <- readHTSeqFile(list[i], sample = 'counts')
  
  counts <- c()
  
  for (j in 1:length(counts_genes)) {
    counts[j] <- raw_counts$counts[sapply(ensmbl[j], grepl, raw_counts$feature, ignore.case=TRUE)]
  }
  counts_genes <- cbind(counts_genes, as.numeric(counts))
}

#reformat counts matrix to make applicable to DESeq2
rownames(counts_genes) <- counts_genes[,1]
counts_genes <- counts_genes[,-1]
counts_genes <- as.matrix(counts_genes)
colnames(counts_genes) <- c('counts1','counts2')
counts_genes <- as.data.frame(counts_genes)
counts_genes_tmp <- apply(as.matrix.noquote(counts_genes),2,as.numeric)
rownames(counts_genes_tmp) <- rownames(counts_genes)
counts_genes <- counts_genes_tmp

#convert counts to log counts
#counts_genes$log_counts <- log10(as.numeric(counts_genes$counts)+1)

#normalise log raw counts with DESeq2:
#create DESeq2 object with dummy metaData object
metaData <- as.factor(c(rep("Ctr",1),rep("Irr",1)))
mycols <- data.frame(row.names = c('counts1','counts2'),metaData)

dds <- DESeqDataSetFromMatrix(counts_genes, DataFrame(metaData), ~ metaData)

#doesnt work:
#dds <- DESeq(dds)

#from paper/ github:
dds <- estimateSizeFactors(dds)
normalised_counts <- counts(dds, normalized = TRUE)
