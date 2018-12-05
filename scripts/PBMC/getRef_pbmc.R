library(Seurat)
load('pbmc.RData')
source('scRef.R')
####################
  
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
#################
exp_ref_mat=read.table('PeripheralBlood_ref_human.txt',header=T,row.names=1,sep='\t',check.name=F) 

REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 

out=SCREF(exp_sc_mat, NewRef)
tag1=out$tag1
tag2=out$tag2
pbmc@meta.data$tag1=tag1[,2]
pbmc@meta.data$tag2=tag2[,2]
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='tag1')
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='tag2')

adjustedRandIndex(pbmc@ident,tag1[,2])
#0.6084547
adjustedRandIndex(pbmc@ident,tag2[,2])
#0.7172649
