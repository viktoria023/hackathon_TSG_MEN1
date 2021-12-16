library(GenomicDataCommons)
library(dplyr)
setwd("~/Desktop/Hackhathon")

#load curated gene list
DDR_epi_genes_names <- readRDS("DDR_epi_genes_names.rds")

#add Menin to gene list
Men1 <- data.frame("ENSG00000133895","MEN1","tumor suppressor associated with multiple endocrine neoplasia type 1", "GeneCards")
rownames(Men1) <- Men1[,1]
names(Men1) <- c("Gene.stable.ID","Gene.name","Gene.description","Sourcelist")
DDR_epi_genes_names <- rbind(DDR_epi_genes_names, Men1)

#read in filenames per cancer from list of lists and iterate through list
file_names <- readRDS("cancer_files_lists.rds")
cancers <- names(file_names)

setwd("~/Desktop/Hackhathon/cancer-data-all-files")

#k==chooses cancer type from list. 1 is already done, and I will continue with 2-4
k=1
list <- file_names[[k]]

#extract counts for genes of interest from datasets
counts_genes <- ensmbl <- DDR_epi_genes_names$Gene.stable.ID

#upload all files and make them as dataframe
for (i in 1:length(list)){
  #id <- paste0("raw_counts_",i)
  assign(paste("raw_counts",i, sep = "_"), as.data.frame(readHTSeqFile(list[i], sample = 'counts')))
  #names(get(paste("raw_counts",i, sep = "_"))) <- c("feature", id) # does not work
  if (i == 1){
    counts <- get(paste("raw_counts",i, sep = "_"))
    next
    }
  counts <- merge(counts, get(paste("raw_counts",i, sep = "_")), by = "feature", all = TRUE )
  }
colnames(counts) <- c("feature", 1:length(list))

counts$feature <- gsub("\\..*", "", counts$feature)

counts <- counts %>% 
  filter(feature %in% counts_genes)

# This is the dataframe above, it is called "counts"
# now we need to make it formatted for DESeq2

#reformat counts matrix to make applicable to DESeq2
rownames(counts_genes) <- counts_genes[,1]
counts_genes <- counts_genes[,-1]
counts_genes <- as.matrix(counts_genes)
colnames(counts_genes) <- NULL

counts_genes <- as.data.frame(counts_genes)
counts_genes_tmp <- apply(as.matrix.noquote(counts_genes),2,as.numeric)
rownames(counts_genes_tmp) <- rownames(counts_genes)
counts_genes <- counts_genes_tmp
counts_genes <- as.data.frame(counts_genes)

#save as a list of dataframes
cancer_data <- list(counts_genes)
savefile<-paste0(paste0(cancers[k]), "_raw_data.rds")
saveRDS(cancer_data, file = savefile)
