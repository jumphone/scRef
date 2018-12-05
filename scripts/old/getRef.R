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

tag=SCREF(exp_sc_mat, NewRef)$tag2
pbmc@meta.data$scref=tag[,2]
TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='scref')

adjustedRandIndex(pbmc@ident,tag[,2])

