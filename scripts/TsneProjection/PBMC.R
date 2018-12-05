load('pbmc.RData')
source('scRef.R')

library(Seurat)

######human pbmc####
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

sc_tag=cbind(names(pbmc@ident),as.character(pbmc@ident))    
exp_raw_mat=as.matrix(pbmc@raw.data)[,COL]
exp_sc_mat=exp_raw_mat[,which(sc_tag[,2]=='5')] # NK cells




######mouse blood###
raw_mat=read.table('PeripheralBlood_com.txt.anno.txt.human',header=T,row.names=1,sep='\t')
pbmc <- CreateSeuratObject(raw.data = raw_mat, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(pbmc,do.label = TRUE, pt.size = 0.5)

#######tSNE plot alignment####
ref_vec=pbmc@dr$tsne@cell.embeddings
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

ref_tag=cbind(names(pbmc@ident),as.character(pbmc@ident))    
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]
LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )  

sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2
out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, min_cell=10, CPU=4, print_step=10)


XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))

plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70',cex=0.5)
par(new=T)
plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red',cex=0.5)



#######CCA####

# we follow the instruction in https://satijalab.org/seurat/immune_alignment.html
set.seed(123)
ctrl <- CreateSeuratObject(raw.data = exp_ref_mat, project = "SIM", min.cells = 5)
ctrl@meta.data$stim <- "CTRL"
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)

stim <- CreateSeuratObject(raw.data = exp_sc_mat, project = "SIM", min.cells = 5)
stim@meta.data$stim <- "STIM"
stim <- NormalizeData(stim)
stim <- ScaleData(stim, display.progress = F)

ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))

immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30, add.cell.id1='Allmouse', add.cell.id2='NKhuman')

immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", 
    dims.align = 1:10)


immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:10, 
    do.fast = T)

#pdf('simulationresult_CCA.pdf',width=4.5, height=5)
CEX=0.5
BLUE=rgb(0, 0, 255, 150, maxColorValue=255)
RED=rgb(255, 0, 0, 150, maxColorValue=255)
ALLVEC=immune.combined@dr$tsne@cell.embeddings
XLIM=c(min(ALLVEC[,1])-1,max(ALLVEC[,1])+1)
YLIM=c(min(ALLVEC[,2])-1,max(ALLVEC[,2])+1)
plot(ALLVEC, pch=16, col='grey70',xlim=XLIM,ylim=YLIM,cex=CEX)
par(new=T)
USE=which(immune.combined@meta.data$orig.ident=='NK.cell')
plot(ALLVEC[USE,], pch=16, col=BLUE,xlim=XLIM,ylim=YLIM,cex=CEX)
par(new=T)
plot(ALLVEC[which(immune.combined@ident=='NKhuman'),], pch=16, col=RED,xlim=XLIM,ylim=YLIM,cex=CEX)
#dev.off()



save.image('tsnealign_pbmc.RObj')


