library(Seurat)
load('TSNE.RData')

TAG=read.table('Zeisel_semi.txt',header=T,sep='\t')
TAG[,2]=as.character(TAG[,2])
TAG[which(TAG[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocytes'
TAG[which(TAG[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocytes'
pbmc@meta.data$zhang2014 = TAG[,2]#as.character(pbmc@ident)

#Astrocytes
#Endothelial.Cells
#Microglia
#Neuron 
#Oligodendrocyte.Precursor.Cell
#Oligodendrocytes


dim(exp_sc_mat)
source('scRef.R')

exp_ref_mat=read.table('MCA_Brain_ref_mouse.txt',header=T,sep='\t',row.names=1)
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}

REF_TAG[which(REF_TAG=='Myelinating.oligodendrocyte.Brain.')]='Oligodendrocytes'
REF_TAG[which(REF_TAG=='Oligodendrocyte.precursor.cell.Brain.')]='OPC'

NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
tag=out$tag2

pbmc@meta.data$mca = as.character(tag[,2])

pbmc@meta.data$mca[which(pbmc@meta.data$mca=='Astrocyte')]='Astrocytes'
pbmc@meta.data$mca[which(pbmc@meta.data$mca=='Neuron.Brain.')]='Neuron'
pbmc@meta.data$mca[which(pbmc@meta.data$mca=='OPC')]='Oligodendrocyte.Precursor.Cell'
pbmc@meta.data$mca[which(pbmc@meta.data$mca=='Hypothalamic.ependymal.cell.Brain.')]='Hypothalamic.ependymal.Cell'
pbmc@meta.data$mca[which(pbmc@meta.data$mca=='Microglia.Brain.')]='Microglia'
pbmc@meta.data$mca[which(pbmc@meta.data$mca=='Oligodendrocytes')]='Oligodendrocytes'
pbmc@meta.data$mca[which(pbmc@meta.data$mca=='Schwann.cell.Brain.')]='Schwann.Cell'
#TSNEPlot(pbmc,group.by='mca')



exp_ref_mat=read.table('MBA_GRCm38.81.P60Cortex_noRep5_FRONTALonly_mouse.txt',header=T,sep='\t',row.names=1)
REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "\\.")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}

#REF_TAG[which(REF_TAG=='Myelinating.oligodendrocyte.Brain')]='Oligodendrocytes'
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 
out=SCREF(exp_sc_mat, NewRef)
tag=out$tag2
pbmc@meta.data$mba = as.character(tag[,2])

pbmc@meta.data$mba[which(pbmc@meta.data$mba=='Astrocyte')]='Astrocytes'
pbmc@meta.data$mba[which(pbmc@meta.data$mba=='Neuron.Brain.')]='Neuron'
pbmc@meta.data$mba[which(pbmc@meta.data$mba=='OPC')]='Oligodendrocyte.Precursor.Cell'
pbmc@meta.data$mba[which(pbmc@meta.data$mba=='Endothelial_Stalk')]='Endothelial.Cells'
pbmc@meta.data$mba[which(pbmc@meta.data$mba=='NeuronSlc17a7')]='Neuron'
pbmc@meta.data$mba[which(pbmc@meta.data$mba=='Oligodendrocyte')]='Oligodendrocytes'



ori=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori[,2]

pdf('MultiRef.pdf',width=7,height=7)
TSNEPlot(pbmc,group.by='ori')
TSNEPlot(pbmc,group.by='zhang2014')
TSNEPlot(pbmc,group.by='mca')
TSNEPlot(pbmc,group.by='mba')
dev.off()

save.image('MultiRef.RData')

NRM=which(!pbmc@meta.data$ori %in% c('pyramidal SS','pyramidal CA1','interneurons'))

adjustedRandIndex(pbmc@meta.data$ori[NRM],pbmc@meta.data$zhang2014[NRM])
#0.9133551
adjustedRandIndex(pbmc@meta.data$ori[NRM],pbmc@meta.data$mba[NRM])
#0.8134128
adjustedRandIndex(pbmc@meta.data$ori[NRM],pbmc@meta.data$mca[NRM])
#0.8738621


CON=which(pbmc@meta.data$zhang2014[NRM]==pbmc@meta.data$mba[NRM] & pbmc@meta.data$mba[NRM]==pbmc@meta.data$mca[NRM])

adjustedRandIndex(pbmc@meta.data$ori[NRM][CON],pbmc@meta.data$zhang2014[NRM][CON])
#0.9784506

