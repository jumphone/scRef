
library(Seurat)
source('scRef.R')



exp_raw_data=read.table('GSE72056_melanoma_single_cell_revised_v2.txt',header=T,row.names=1)
exp_raw_data=exp_raw_data[3:length(exp_raw_data[,1]),]

pbmc <- CreateSeuratObject(raw.data = exp_raw_data, min.cells = 0, min.genes = 0, project = "10X_PBMC")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)    
