
library(Seurat)

exp_raw_data=read.table('MGH54_mat.txt',sep='\t',header=T,row.names=1)
pbmc <- CreateSeuratObject(raw.data = exp_raw_data, min.cells = 3, min.genes = 200,  project = "10X_PBMC")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)

save(pbmc,file='OligTumor.RData')


load('OligTumor.RData')
source('scRef.R')

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
#################
exp_ref_mat=read.table('Brain_ref_human.txt',header=T,row.names=1,sep='\t',check.name=F) 

REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 

out=SCREF(exp_sc_mat, NewRef)
tag1=out$tag1
tag2=out$tag2
pbmc@meta.data$tag1=tag1[,2]
pbmc@meta.data$tag2=tag2[,2]
#TSNEPlot(object = pbmc, do.label=T, group.by ='tag1', pt.size = 0.5)
TSNEPlot(object = pbmc, do.label=T, group.by ='tag2', pt.size = 0.5)



