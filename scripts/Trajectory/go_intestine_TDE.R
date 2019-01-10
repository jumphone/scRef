library('Seurat')
exp_mat=read.table('GSE92332_atlas_UMIcounts_reheader.txt',sep='\t',header=T,row.names=1)
pbmc <- CreateSeuratObject(raw.data = exp_mat, min.cells = 0, min.genes = 0, project = "10X_PBMC")
tag=cbind(as.character(pbmc@meta.data$orig.ident),as.character(pbmc@meta.data$orig.ident))
source('scRef.R')
NewRef=.generate_ref(exp_mat,tag)
out=.get_log_p_sc_given_ref(pbmc@raw.data, NewRef)
result=.trajectory(out, plot_type='polygon', plot_size=1.7, label_dist=1.2, label_size=10, random_ratio=0.03)

png(filename = "TDE.png",width = 1024, height = 1024)
result$ggplot
dev.off()

library(igraph)

MST=.generate_mst(result$mat)
pdf('MST.pdf')
plot(MST)
dev.off()

save.image('TDE_intestine.RData')
