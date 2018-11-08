
exp_data=read.table('expression_matrix.csv',sep=',',row.names=1,header=F)
gene=read.table('rows_metadata.csv',sep=',',header=T)[,4]
rownames(exp_data)=gene


