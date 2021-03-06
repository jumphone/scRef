
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


######## Cluster & Unknown Cell types#########

library('Seurat')
load('TSNE.RData')

pbmc_3 <- RunTSNE(object = pbmc, dims.use = 1:150, do.fast = TRUE, dim.embed = 3)
TSNE_VEC=pbmc_3@dr$tsne@cell.embeddings


# Hierarchical
D=dist(TSNE_VEC)
H=hclust(D)
C=cutree(H,k=8)

# K-means
#C=kmeans(TSNE_VEC,centers=10)$cluster

old_ident = pbmc@ident
pbmc@ident = as.factor(C)
names(pbmc@ident)=names(old_ident)
TSNEPlot(object = pbmc, do.label=T)



TAG_cluster=cbind(rownames(TSNE_VEC),C)
TAG=read.table('Zeisel_semi.txt',header=T,sep='\t')

.compare_two_tag <- function(TAG1,TAG2){
    OUT=c()
    tag1_names=as.character(unique(TAG1[,2]))
    tag2_names=as.character(unique(TAG2[,2]))
    i=1
    while(i<=length(tag1_names)){
        tag1 = tag1_names[i]
        tag1_index = which(TAG1[,2]== tag1)
        j=1
        while(j<=length(tag2_names)){
            tag2 = tag2_names[j] 
            #print(tag2)
            tag2_index = which(TAG2[,2]== tag2)
            over = length(which(tag1_index %in% tag2_index))
            tag1_over = over/length(tag1_index)
            tag2_over = over/length(tag2_index)
            max_over = max(tag1_over, tag2_over)
            OUT=cbind(OUT, c(max_over, tag1, tag1_over ,tag2, tag2_over)) 
            j=j+1
            }
        i=i+1
        }
    OUT=t(OUT)
    #OUT=as.matrix(OUT)
    OUT[,1]=as.numeric(OUT[,1])
    OUT[,3]=as.numeric(OUT[,3])
    OUT[,5]=as.numeric(OUT[,5])
    colnames(OUT)=c('max_over','tag1','tag1_over','tag2','tag2_over')
    OUT=OUT[order(OUT[,1],decreasing=T),]
    return(OUT)
    }


OUT=.compare_two_tag(TAG_cluster,TAG)

write.table(OUT,'COM.txt',sep='\t',quote=F,col.names=T,row.names=F)


######## wokflow demo #########
load('TSNE.RData')
library('Seurat')
source('scRef.R')
ori=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
ken=read.table('Kendall.txt',header=T,sep='\t')
semi=read.table('Zeisel_semi.txt',header=T,sep='\t')

pbmc@meta.data$ori=ori[,2]
pbmc@meta.data$ken=ken[,2]
pbmc@meta.data$semi=semi[,2]

TSNEPlot(pbmc, group.by='ori')

MIC=which(ori[,2]=='microglia')
MICK=which(ori[,2]=='microglia' & ken[,2]=='Microglia')
MICS=which(ori[,2]=='microglia' & semi[,2]=='Microglia')

OP=which(ori[,2]=='oligodendrocytes')
OPK=which( (ken[,2]=='Newly.Formed.Oligodendrocyte'| ken[,2]=='Myelinating.Oligodendrocytes'))
OPS=which( (semi[,2]=='Newly.Formed.Oligodendrocyte'| semi[,2]=='Myelinating.Oligodendrocytes'))

AC=which(ori[,2]=='astrocytes_ependymal')
ACK=which( ken[,2]=='Astrocytes')
ACS=which( semi[,2]=='Astrocytes')

EN=which(ori[,2]=='endothelial-mural')
ENK=which(ken[,2]=='Endothelial.Cells')
ENS=which(semi[,2]=='Endothelial.Cells')

NU=which(ori[,2]=='interneurons' | ori[,2]=='pyramidal SS' | ori[,2]=='pyramidal CA1')
NUK=which(ken[,2]=='Neuron')
NUS=which(semi[,2]=='Neuron')



length(AC)
length(ACK)
length(which(AC %in% ACK))
length(ACS)
length(which(AC %in% ACS))

length(EN)
length(ENK)
length(which(EN %in% ENK))
length(ENS)
length(which(EN %in% ENS))

length(MIC)
length(MICK)
length(which(MIC %in% MICK))
length(MICS)
length(which(MIC %in% MICS))

length(OP)
length(OPK)
length(which(OP %in% OPK))
length(OPS)
length(which(OP %in% OPS))

length(NU)
length(NUK)
length(which(NU %in% NUK))
length(NUS)
length(which(NU %in% NUS))

TSNE_VEC=pbmc@dr$tsne@cell.embeddings

plot(TSNE_VEC[OP,1],TSNE_VEC[OP,2],pch=16)

