
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




