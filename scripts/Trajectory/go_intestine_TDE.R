library('Seurat')
exp_mat=read.table('GSE92332_atlas_UMIcounts_reheader.txt',sep='\t',header=T,row.names=1)
pbmc <- CreateSeuratObject(raw.data = exp_mat, min.cells = 0, min.genes = 0, project = "10X_PBMC")
tag=cbind(as.character(pbmc@meta.data$orig.ident),as.character(pbmc@meta.data$orig.ident))
source('scRef.R')
NewRef=.generate_ref(exp_mat,tag)



