library('Seurat')
exp_mat=read.table('GSE92332_atlas_UMIcounts_reheader.txt',sep='\t',header=T,row.names=1)
pbmc <- CreateSeuratObject(raw.data = exp_mat, min.cells = 0, min.genes = 0, project = "10X_PBMC")
tag=cbind(as.character(pbmc@meta.data$orig.ident),as.character(pbmc@meta.data$orig.ident))
source('scRef.R')
NewRef=.generate_ref(exp_mat,tag)
out=.get_log_p_sc_given_ref(pbmc@raw.data, NewRef)
result=.trajectory(out, plot_type='polygon', plot_size=1.7, label_dist=1.2, label_size=10, random_ratio=0.03)

png(filename = "TDE.png",width = 1024, height = 1024)
result$ggplot
dev.off()

library(igraph)

MST=.generate_mst(result$mat)
pdf('MST.pdf')
plot(MST)
dev.off()

save.image('TDE_intestine.RData')

###########TSCAN##########################

library(TSCAN)
COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
VAR_GENE=apply(exp_sc_mat,1,var)
USED_GENE=which(VAR_GENE>0)
used_data=exp_sc_mat[USED_GENE,]


procdata <- preprocess(used_data,clusternum=15, pseudocount = 1, minexpr_value = 1, minexpr_percent = 0.4, cvcutoff = 0.5)
lpsmclust <- exprmclust(procdata, clusternum=15 )
tmp=names(lpsmclust$clusterid)
lpsmclust$clusterid=pbmc@meta.data$orig.ident
names(lpsmclust$clusterid)=tmp

pdf('TSCAN_InstestineDEV.pdf',width=14,height=14)
plotmclust(lpsmclust,show_cell_names=F,cell_name_size =1)
dev.off()
