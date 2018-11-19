load('TSNE.RData')
library('Seurat')
source('scRef.R')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]

TSNEPlot(pbmc, group.by="ori")


