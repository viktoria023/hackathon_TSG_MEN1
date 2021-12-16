#Load_data
freq=readRDS('candidates.rds')
df_list=readRDS('gtex_norm_df_list.rds')
M=readRDS('correlation_matrix.rds') #correlation matrix

men1="ENSG00000133895.14"
#Plot data for the ith candidate gene (i=1 being the most widely correlated across tissues)
i=1
candidate=freq[i,1]

#check the 12 tissues where correlation of candidate with MEN1 is high
tissues_ix=order(abs(M[candidate,]),decreasing=TRUE)[1:12]
par(mfrow=c(3,4))
for (k in tissues_ix){
  cur_df=df_list[[k]]
  cur_men1_data=cur_df[men1,]
  cur_cand_data=cur_df[candidate,]
  plot(cur_men1_data,cur_cand_data,main=colnames(M)[k],xlab='MEN1',ylab='RAD9A',cex.main=1,col.main='red')
}
