library('Seurat')
library(dplyr)
load('this.RData')
ls()

exp_ref_mat=a

#########PCA########


pca=prcomp(exp_ref_mat,scale=T,center=T)
PC1=pca$rotation[,1]
PC2=pca$rotation[,2]
xmin=min(PC1)-0.1
xmax=max(PC1)+0.1
ymin=min(PC2)-0.1
ymax=max(PC2)+0.1
LABEL=colnames(exp_ref_mat)
plot(PC1,PC2,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16)
text(PC1, PC2, labels=LABEL, cex= 0.7, pos=c(1,2,3,4))


PC1=pca$rotation[,1]
PC3=pca$rotation[,3]
xmin=min(PC1)-0.1
xmax=max(PC1)+0.1
ymin=min(PC3)-0.1
ymax=max(PC3)+0.1
LABEL=colnames(exp_ref_mat)
plot(PC1,PC3,xlim=c(xmin,xmax),ylim=c(ymin,ymax),pch=16)
text(PC1, PC3, labels=LABEL, cex= 0.7, pos=c(1,2,3,4))


#########tSNE workflow########

pbmc.data <- exp_sc_mat
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 0, min.genes = 0)

allgene = rownames(exp_sc_mat)
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
pbmc <- RunPCA(object = pbmc, pcs.compute=150, pc.genes = allgene, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:150, do.fast = TRUE)

old_ident = pbmc@ident
pbmc@ident = as.factor(TAG[,2])
names(pbmc@ident)=names(old_ident)
pbmc_ori=pbmc
pbmc_ori@ident = as.factor(h[,5])
names(pbmc_ori@ident)=names(old_ident)


par(mfrow=c(1,2))
TSNEPlot(object = pbmc)
TSNEPlot(object = pbmc_ori)
