
library(Seurat)
source('scRef.R')
load('CASE.RObj')


TSNE_VEC=pbmc@dr$tsne@cell.embeddings
IDENT=as.character(pbmc@ident)



COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

tag=cbind(names(pbmc@ident),pbmc@ident)
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)


exp_sc_mat=read.table('WT.txt',sep='\t',check.name=F,row.names=1,header=T)

out=.get_cor(exp_sc_mat, LocalRef, method='kendall',CPU=10, print_step=10)
tag=.get_tag_max(out)




i=1
while(i<=length(tag[,1])){
    this_tag = as.character(tag[i,2])
    
    vec_index=which(IDENT==this_tag)
    
    
    exp_ref_mat[,vec_index]
    
    i=i+1}


library("RColorBrewer")
display.brewer.all()
brewer.pal(n=8,name='Set2')
#install.packages("wesanderson")
library(wesanderson)

