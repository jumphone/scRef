load('TSNE.RData')
library('Seurat')
source('scRef.R')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]
#TSNEPlot(pbmc, group.by="ori")
table(pbmc@meta.data$ori)


####### Data preparation ##############
ref_vec=pbmc@dr$tsne@cell.embeddings

USE=which(pbmc@meta.data$ori=='astrocytes_ependymal')

plot(ref_vec[,1],ref_vec[,2], pch=16, col='gray80', xlim=c(-30, 30), ylim=c(-40,35), xlab='tSNE_1',ylab='tSNE_2')
par(new=T)
plot(ref_vec[USE,1],ref_vec[USE,2], pch=16, col='red', xlim=c(-30, 30), ylim=c(-40,35), xlab='',ylab='')


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

ref_tag=cbind(names(pbmc@ident), as.character(pbmc@meta.data$ori))    
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]


exp_sc_mat= exp_ref_mat[,USE]

getRanGene <- function(X){
    POS = which(X >0 )
    N=length(POS)/2
    KEEP = sample(x=POS, size=N )
    NEG = POS[which(!POS %in% KEEP)]
    X[NEG]=0
    return(X)
    }

set.seed(123)
sim_exp_sc_mat = apply(exp_sc_mat,2, getRanGene)

###### tSNE projection (default) ############

LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )
out=.get_cor(sim_exp_sc_mat, LocalRef, method='kendall',CPU=4, print_step=10)
tag=.get_tag_max(out)
LocalRef= .generate_ref(sim_exp_sc_mat, tag, min_cell = 10 )
out=.get_log_p_sc_given_ref(sim_exp_sc_mat, LocalRef, CPU=10, print_step=10)
tag=.get_tag_max(out)
sc_tag=tag

out =.vec_projection(exp_sc_mat=sim_exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, min_cell=10, CPU=4, print_step=10)

pdf('simulationresult_tSNEprojection.pdf',width=4.5, height=5)
CEX=0.7
BLUE=rgb(0, 0, 255, 150, maxColorValue=255)
RED=rgb(255, 0, 0, 150, maxColorValue=255)
plot(ref_vec,xlim=c(-30, 30), ylim=c(-40,35),pch=16,col='grey70', cex=CEX)
par(new=T)
plot(ref_vec[USE,1],ref_vec[USE,2],xlim=c(-30, 30), ylim=c(-40,35),pch=16,col=BLUE, cex=CEX, xlab='',ylab='')
par(new=T)
plot(out$vec,xlim=c(-30, 30), ylim=c(-40,35),pch=16,col=RED, cex=CEX)
dev.off()

###### tSNE projection (All) ############
out_all =.vec_projection(exp_sc_mat=sim_exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=1, alpha=0.5, random_size=3000, 
        random_seed=123, min_cell=10, CPU=4, print_step=10)

pdf('simulationresult_tSNEprojection_all.pdf',width=4.5, height=5)
CEX=0.7
BLUE=rgb(0, 0, 255, 150, maxColorValue=255)
RED=rgb(255, 0, 0, 150, maxColorValue=255)
plot(ref_vec,xlim=c(-30, 30), ylim=c(-40,35),pch=16,col='grey70', cex=CEX)
par(new=T)
plot(ref_vec[USE,1],ref_vec[USE,2],xlim=c(-30, 30), ylim=c(-40,35),pch=16,col=BLUE, cex=CEX, xlab='',ylab='')
par(new=T)
plot(out_all$vec,xlim=c(-30, 30), ylim=c(-40,35),pch=16,col=RED, cex=CEX)
dev.off()

###### CCA ############

# we follow the instruction in https://satijalab.org/seurat/immune_alignment.html
set.seed(123)
ctrl <- CreateSeuratObject(raw.data = exp_ref_mat, project = "SIM", min.cells = 5)
ctrl@meta.data$stim <- "CTRL"
ctrl <- NormalizeData(ctrl)
ctrl <- ScaleData(ctrl, display.progress = F)

stim <- CreateSeuratObject(raw.data = sim_exp_sc_mat, project = "SIM", min.cells = 5)
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

immune.combined <- RunCCA(ctrl, stim, genes.use = genes.use, num.cc = 30, add.cell.id1='All', add.cell.id2='Sim')

immune.combined <- AlignSubspace(immune.combined, reduction.type = "cca", grouping.var = "stim", 
    dims.align = 1:25)

# We have tried to use 10, 15, 20, and 25 aligned CCA subspaces to draw tSNE plot. 
# And we can see a clear "astrocytes_ependymal" cluster when using 10 aligned CCA subspaces.

immune.combined <- RunTSNE(immune.combined, reduction.use = "cca.aligned", dims.use = 1:10, 
    do.fast = T)

pdf('simulationresult_CCA.pdf',width=4.5, height=5)
CEX=0.7
BLUE=rgb(0, 0, 255, 150, maxColorValue=255)
RED=rgb(255, 0, 0, 150, maxColorValue=255)
ALLVEC=immune.combined@dr$tsne@cell.embeddings
XLIM=c(min(ALLVEC[,1])-1,max(ALLVEC[,1])+1)
YLIM=c(min(ALLVEC[,2])-1,max(ALLVEC[,2])+1)
plot(ALLVEC, pch=16, col='grey70',xlim=XLIM,ylim=YLIM,cex=CEX)
par(new=T)
plot(ALLVEC[USE,], pch=16, col=BLUE,xlim=XLIM,ylim=YLIM,cex=CEX)
par(new=T)
plot(ALLVEC[which(immune.combined@ident=='Sim'),], pch=16, col=RED,xlim=XLIM,ylim=YLIM,cex=CEX)
dev.off()


#save.image('tSNEprojection.RData')

