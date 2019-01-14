
load('TSNE.RData')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]
library(Seurat)
library(SingleCellExperiment)
library(scmap)
head(ann)


tmp=as.matrix(as.factor(pbmc@meta.data$ori))
#tmp[,1]=as.factor(tmp[,1])
rownames(tmp)=colnames(pbmc@data)
colnames(tmp)="cell_type1"
ann=tmp

set.seed(1)
yan=as.matrix(pbmc@data)
sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
rowData(sce)$feature_symbol <- rownames(sce)
logcounts(sce) <- normcounts(sce)
sce <- selectFeatures(sce, suppress_plot = FALSE)

set.seed(1)
sce <- indexCell(sce)
names(metadata(sce)$scmap_cell_index)
length(metadata(sce)$scmap_cell_index$subcentroids)
dim(metadata(sce)$scmap_cell_index$subcentroids[[1]])
metadata(sce)$scmap_cell_index$subcentroids[[1]][,1:5]


################

USE=which(pbmc@meta.data$ori=='astrocytes_ependymal')
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
sim_ann=as.matrix(ann[which(ann[,1]=='astrocytes_ependymal'),])
colnames(sim_ann)='cell_type1'

ssce=SingleCellExperiment(assays = list(normcounts = as.matrix(sim_exp_sc_mat)), colData = sim_ann)
rowData(ssce)$feature_symbol <- rownames(ssce)
logcounts(ssce) <- normcounts(ssce)
ssce <- selectFeatures(ssce, suppress_plot = FALSE)

##############

scmapCell_results <- scmapCell(
  ssce, 
  list(
    yan = metadata(sce)$scmap_cell_index
  )
)


scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, threshold = 0,w=1,
  list(
    as.character(colData(sce)$cell_type1)
  )
)

head(scmapCell_clusters$combined_labs)


plot(
  getSankey(
    colData(ssce)$cell_type1, 
    scmapCell_clusters$scmap_cluster_labs[,"yan"],
    plot_height = 400
  )
)

####################################






