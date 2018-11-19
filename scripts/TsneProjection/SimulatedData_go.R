load('TSNE.RData')
library('Seurat')
source('scRef.R')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]

TSNEPlot(pbmc, group.by="ori")

table(pbmc@meta.data$ori)
ref_vec=pbmc@dr$tsne@cell.embeddings

#COL=rep('gray80',length(ref_vec[,1]))


OLIG=which(pbmc@meta.data$ori=='oligodendrocytes')

plot(ref_vec[,1],ref_vec[,2], pch=16, col='gray80', xlim=c(-30, 30), ylim=c(-40,35), xlab='tSNE_1',ylab='tSNE_2')
par(new=T)
plot(ref_vec[OLIG,1],ref_vec[OLIG,2], pch=16, col='red', xlim=c(-30, 30), ylim=c(-40,35), xlab='',ylab='')


