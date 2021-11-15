###Hyper geometric tests for clusters
assigned_names = read.delim('/Users/vickyingham/Dropbox (LSTM)/Vicky Documents/Post Doc/Visitors/Henry Youd/AllGenes_with_added.txt',header=T)

###Read in your own data. The gene ID file should match the by.x name below
all_data = read.delim('', header=T)

##Change this to gambiae column
hyper_data = merge(all_data,assigned_names,by.x = 'AnGamTranscriptstableID',by.y = 'AnGamTranscriptstableID',all=T)
hyper_data = hyper_data[,-c((ncol(hyper_data)-4):(ncol(hyper_data)-1))]
hyper_data = as.matrix(hyper_data)

hyper_data[is.na(hyper_data[,ncol(hyper_data)])==T,ncol(hyper_data)] = 'Unknown'

assignments = unique(hyper_data[,ncol(hyper_data)])  
assignments = assignments[-which(assignments== 'Unknown')]

library(stringr)
hyper_tests = c()
total_genes = nrow(all_data)
for(j in 1:length(cl))
{
  p_val = c()
  cluster_members = as.data.frame(str_replace(as.character(cl[[j]]$NAME), '\\.', '-'))
  colnames(cluster_members) = 'gene'
  ###Change to gambiae column
  cluster_assignment = merge(cluster_members,hyper_data,by.x='gene',by.y='AnGamTranscriptstableID',all=F)
  for(i in 1:length(assignments))
  {
    total_on_array = length(which(hyper_data[,ncol(hyper_data)]==assignments[i]))
    total_in_cluster = length(which(cluster_assignment[,ncol(cluster_assignment)]==assignments[i]))
    p_val = c(p_val,p.adjust(1-phyper(total_in_cluster-1,total_on_array,(total_genes-total_on_array),nrow(cluster_members)),method='fdr'))
  }
  hyper_tests = rbind(hyper_tests,c(p_val,j))
}
colnames(hyper_tests) = c(assignments,'cluster')

write.table(hyper_tests,'/Users/vickyingham/Dropbox (LSTM)/Vicky Documents/Post Doc/Visitors/Henry Youd/Clusters/hyper_tests.txt',sep='\t',row.names=F)
