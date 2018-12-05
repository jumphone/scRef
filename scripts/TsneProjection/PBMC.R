load('pbmc.RData')
source('scRef.R')

library(Seurat)


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





raw_mat=read.table('PeripheralBlood_com.txt.anno.txt.human',header=T,row.names=1,sep='\t')
pbmc <- CreateSeuratObject(raw.data = raw_mat, min.cells = 3, min.genes = 200, 
    project = "10X_PBMC")

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))






