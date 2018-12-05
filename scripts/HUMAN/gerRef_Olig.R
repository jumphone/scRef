
library(Seurat)
load('pbmc.RData')
source('scRef.R')

########
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", 
    "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)
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
TSNEPlot(object = pbmc, do.label=T, group.by ='tag1', pt.size = 0.5)
TSNEPlot(object = pbmc, do.label=T, group.by ='tag2', pt.size = 0.5)


merged_label=as.character(pbmc@ident)
merged_label[which(merged_label=='CD4 T cells')]='T cell'
merged_label[which(merged_label=='CD8 T cells')]='T cell'
merged_label[which(merged_label=='CD14+ Monocyte')]='Monocyte'
merged_label[which(merged_label=='FCGR3A+ Monocytes')]='Monocyte'
#merged_label[which(merged_label=='Megakaryocytes')]='Monocyte'

library(mclust)
adjustedRandIndex(merged_label,tag1[,2])
#0.7612499

sout=.get_cor(exp_sc_mat, NewRef, method='spearman',CPU=4, print_step=10)
stag=.get_tag_max(sout)
adjustedRandIndex(merged_label,stag[,2])
#0.7854657

pout=.get_cor(exp_sc_mat, NewRef, method='pearson',CPU=4, print_step=10)
ptag=.get_tag_max(pout)
adjustedRandIndex(merged_label,ptag[,2])
#0.6203971

mout=.get_log_p_sc_given_ref(exp_sc_mat, NewRef, CPU=4, print_step=10)
mtag=.get_tag_max(mout)
adjustedRandIndex(merged_label,mtag[,2])
#0.7340889

adjustedRandIndex(merged_label,tag2[,2])
#0.8189932
