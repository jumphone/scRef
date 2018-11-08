
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


TAG=read.table('Zeisel_semi.txt',header=T,sep='\t')
TAG[,2]=as.character(TAG[,2])
TAG[which(TAG[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocytes'
TAG[which(TAG[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocytes'
TAG[which(TAG[,2]=='Oligodendrocyte.Precursor.Cell'),2]='OPC'

h=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')

old_ident = pbmc@ident
pbmc@ident = as.factor(TAG[,2])
names(pbmc@ident)=names(old_ident)
pbmc_ori=pbmc
pbmc_ori@ident = as.factor(h[,2])
names(pbmc_ori@ident)=names(old_ident)

TSNEPlot(object = pbmc)
TSNEPlot(object = pbmc_ori)

save.image(file='TSNE.RData')



#######Error analysis#########################

FLAG=read.table('ERR_FLAG.txt',sep='\t',header=F)
old_ident = pbmc@ident
pbmc@ident = as.factor(as.character(FLAG[,1]))
names(pbmc@ident)=names(old_ident)

TSNEPlot(object = pbmc)

##############OTHER 1######################

library(plotly)
pbmc_3 <- RunTSNE(object = pbmc, dims.use = 1:150, do.fast = TRUE, dim.embed = 3)
TSNE_VEC=pbmc_3@dr$tsne@cell.embeddings
p=plot_ly(x=TSNE_VEC[,1],y=TSNE_VEC[,2],z=TSNE_VEC[,3],color=as.factor(TAG[,2]))
htmlwidgets::saveWidget(p, "kendall.html")
p=plot_ly(x=TSNE_VEC[,1],y=TSNE_VEC[,2],z=TSNE_VEC[,3],color=as.factor(h[,2]))
htmlwidgets::saveWidget(p, "original.html")


##############Error analysis tSNE distance######################

TSNE_VEC=pbmc@dr$tsne@cell.embeddings
D=as.matrix(dist(TSNE_VEC))
OD=apply(D,2,order)

TOP=10
DIS=c()
i=1
while(i<=length(D[1,])){
    this_dis=mean(D[OD[,i][1:TOP],i])
    DIS=c(DIS,this_dis)
    i=i+1}

hist(DIS,breaks=100)
CUTOFF=quantile(DIS,0.95)

COL=rep(0,length(DIS))
COL[which(DIS>=CUTOFF)]=1
old_ident = pbmc@ident
pbmc@ident = as.factor(as.character(COL))
names(pbmc@ident)=names(old_ident)
TSNEPlot(object = pbmc)

HDIS=which(DIS >= CUTOFF)

FLAG=read.table('ERR_FLAG.txt',sep='\t',header=F)
ERR=which(FLAG[,1]==0)

length(which(HDIS %in% ERR))
length(ERR)
length(HDIS)


fisher.test(cbind(c(21,92),c(131,2761)))


######hierarchical cluster#########

library('Seurat')
load('TSNE.RData')

pbmc_3 <- RunTSNE(object = pbmc, dims.use = 1:150, do.fast = TRUE, dim.embed = 3)
TSNE_VEC=pbmc_3@dr$tsne@cell.embeddings

D=dist(TSNE_VEC)
H=hclust(D)
C=cutree(H,k=6)

old_ident = pbmc@ident
pbmc@ident = as.factor(C)
names(pbmc@ident)=names(old_ident)

TSNEPlot(object = pbmc)
















