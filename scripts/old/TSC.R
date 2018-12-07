
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TSCAN", version = "3.8")

library(TSCAN)
load('OligTumor.RData')

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

procdata <- preprocess(exp_sc_mat,clusternum=4)
lpsmclust <- exprmclust(procdata, clusternum =4)
plotmclust(lpsmclust,show_cell_names=F)

