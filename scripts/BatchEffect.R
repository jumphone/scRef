

library('Seurat')

exp_sc_mat=read.table('./Benchmark/Zeisel_exp_sc_mat.txt',header=T,row.names=1,sep='\t')
exp_ref_mat=read.table('Reference_expression.txt',header=T,row.names=1,sep='\t')


sc=CreateSeuratObject(raw.data = exp_sc_mat, min.cells = 0, min.genes = 0)
ref=CreateSeuratObject(raw.data = exp_ref_mat, min.cells = 0, min.genes = 0)

combined <- MergeSeurat(object1 = sc, object2 = ref, add.cell.id1 = "sc", 
    add.cell.id2 = "ref")

pbmc=combined
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#pdf('batch_effect_ref.pdf',width=5,height=4.5)
PT=rep(0.5, length(pbmc@ident))
PT[which(pbmc@ident=='ref')]=3
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2, pt.size=PT)
#dev.off()

source('scRef.R')
out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
tag=.get_tag_max(out)
LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)

locref=CreateSeuratObject(raw.data = LocalRef, min.cells = 0, min.genes = 0)
combined_loc <- MergeSeurat(object1 = sc, object2 = locref, add.cell.id1 = "sc", 
    add.cell.id2 = "locref")
pbmc=combined_loc
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=F)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

#pdf('batch_effect_locref.pdf',width=5,height=4.5)
PT=rep(0.5, length(pbmc@ident))
PT[which(pbmc@ident=='locref')]=3
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2, pt.size=PT)
#dev.off()


