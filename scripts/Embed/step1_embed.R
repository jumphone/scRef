
library(Seurat)
source('scRef.R')
load('CASE.RObj')

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

tag=cbind(names(pbmc@ident),pbmc@ident)
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)


exp_sc_mat=read.table('WT.txt',sep='\t',check.name=F,row.names=1,header=T)

out=.get_cor(exp_sc_mat, LocalRef, method='kendall',CPU=10, print_step=10)
tag=.get_tag_max(out)

W=(1-out)/2

getw=function(X){
  X=X/sum(X)
  return(X)
}

W=apply(W,2,getw)



