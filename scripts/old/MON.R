source("http://bioconductor.org/biocLite.R")
biocLite()
library(monocle)
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")
biocLite(c("DDRTree", "pheatmap"))

load('OligTumor.RData')

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]


expr_matrix <- exp_sc_mat

sample_name=names(pbmc@ident)
group=as.character(pbmc@ident)
sample_sheet=cbind(sample_name, group)
sample_sheet=as.data.frame(sample_sheet)
rownames(sample_sheet)=colnames(expr_matrix)


gene_short_name=rownames(expr_matrix)
gene_type=rownames(gene_name)
gene_annotation=cbind(gene_short_name,gene_type)
gene_annotation=as.data.frame(gene_annotation)
rownames(gene_annotation)=rownames(expr_matrix)


pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)

HSMM <- newCellDataSet(as.matrix(expr_matrix), phenoData = pd, featureData = fd,expressionFamily=negbinomial.size())

diff_test_res <- differentialGeneTest(cds,
    fullModelFormulaStr = "~Media")


HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 6, reduction_method = 'tSNE', verbose = T)
                
HSMM <- clusterCells(HSMM, num_clusters = 5)
plot_cell_clusters(HSMM, 1, 2)




diff_test_res <- differentialGeneTest(HSMM,  fullModelFormulaStr = "~Cluster")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)








