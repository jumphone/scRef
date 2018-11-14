
library(Seurat)
library(dplyr)
pbmc.data <- read.table('CASE.txt',sep='\t',check.name=F,row.names=1,header=T)
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "Natalie")
#save(pbmc,file='CASE_raw.RObj')


mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
pbmc <- AddMetaData(object = pbmc, metadata = pbmc@ident, col.name = "batch")
pbmc=FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(3500, 0.1))
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot = F, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff =0, y.cutoff = 0.5)
length(x=pbmc@var.genes) #3726
pbmc = ScaleData(object = pbmc,vars.to.regress = c("percent.mito", "nUMI", "batch"), genes.use=pbmc@var.genes)

PCNUM=40
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5,    genes.print = 5, pcs.compute=PCNUM, maxit = 500, weight.by.var = FALSE )

PCUSE=1:35
pbmc = RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)

TSNEPlot(pbmc,pt.size=0.5)
#save(pbmc,file='CASE_tsne.RObj')

TSNE_VEC=pbmc@dr$tsne@cell.embeddings
D=dist(TSNE_VEC)
H=hclust(D)
C=cutree(H,k=12)

old_ident = pbmc@ident
pbmc@ident = as.factor(C)
names(pbmc@ident)=names(old_ident)
TSNEPlot(object = pbmc, do.label=T)
