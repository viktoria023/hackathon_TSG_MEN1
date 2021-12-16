#Preliminaries------------
library(DESeq2)
library(dplyr)
#1. Load list of dataframes (each element is a tissue)
filename='tissues_datasets.rds'
df_list=readRDS(filename)

#Select a start column for the counts (the first columns may be gene names, IDs etc.)
sc=3

#MOR Normalisation with DSeq2--------------------
norm_df_list=list() #List of normalised counts dataframes (one dataframe per cancer type/tissue) 
for (k in c(1:length(df_list))){
  #First, create a metadata matrix for the current tissue
  cur_df=df_list[[k]]
  meta=data.frame(sample_id=colnames(cur_df[,sc:ncol(cur_df)]))
  #Create a DESeq2 object
  dds=DESeqDataSetFromMatrix(countData=cur_df[,sc:ncol(cur_df)],colData = meta,design=~sample_id)
  #EstimateSizeFactors
  dds <- estimateSizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  norm_df_list[[k]]=normalized_counts
}
saveRDS(norm_df_list,file='gtex_norm_df_list.rds')
#Correlation_analysis----------
#define the MEN1 ENSG ID for simple search through our dataframe
men1="ENSG00000133895.14"
#Define a list of vectors of correlations, once per datasets
corr_list=list()
#For each tissue/cancer_type, compute Spearman correlation between MEN1 exprresion and the remaining genes
#Then store result in a vector and append to the corr_list
for (k in c(1:length(norm_df_list))){
  cur_df=norm_df_list[[k]] #currant dataset (single tissue or cancer)
  cur_corr_vec=as.vector(matrix(0,nrow=nrow(cur_df))) #correlation vector for the current datset
  men1_data=as.vector(cur_df[men1,])
  #Now compute correlatio of men1 expression with each of the other query genes
  for (i in c(1:nrow(cur_df))){
    cur_gene_data=cur_df[i,]
    cur_corr=cor(men1_data,cur_gene_data,method='spearman')
    cur_corr_vec[i]=cur_corr
  }
  corr_list[[k]]=cur_corr_vec
}
#Name the correlation list with the names of the original list-of-dataframes
names(corr_list)=names(df_list)
#Save the correlation list for future analysis
saveRDS(corr_list,file='gtex_correlations.rds')


#Analysis of co-expression results----------------
#Generate correlation matrix (each row is a gene, each column a tissue)
M=data.frame(matrix(ncol=length(corr_list),nrow=nrow(df_list[[1]])))
colnames(M)=names(corr_list)
rownames(M)=rownames(df_list[[1]])
#Fill-in the correlation matrix
for (k in c(1:ncol(M))){
  cur_corr_vec=corr_list[[k]]
  M[,k]=cur_corr_vec
}
#Remove the row corresponding to Men1 correlation with itself
M=M[!row.names(M) %in% c(men1),]
saveRDS(M,file='correlation_matrix.rds')
#Select % of genes to extract for further analysis
p=0.1
ngenes=round(p*nrow(M))
#Initialise a new dataframe that will contain the p% highest correlated genes with MEN1 in each tissue
Mp=data.frame(matrix(nrow=nrow(M),ncol=ncol(M)))
#and a Map containing the gene ID for each element of Mp
Mp_id=data.frame(matrix(nrow=nrow(M),ncol=ncol(M)))

colnames(Mp)=colnames(M)
colnames(Mp_id)=colnames(M)
#Now select the ngenes most correlated genes in each tissue and the two new data frames
for (k in c(1:ncol(M))){
  sorted_ix=order(abs(M[,k]),decreasing=TRUE)
  
  Mp[,k]=M[sorted_ix,k]
  Mp_id[,k]=rownames(M)[sorted_ix]
}
#Now subset the first p% rows of the two dataframes
Mp=Mp[c(1:ngenes),]
Mp_id=Mp_id[c(1:ngenes),]

#Look for genes that are highly correlated across multiple tissues
frequency=as.data.frame(table(as.vector(as.matrix(Mp_id))))
frequency<-frequency[order(frequency[,2],decreasing=TRUE),]
saveRDS(frequency,file='candidates.rds')







  
