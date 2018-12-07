
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TSCAN", version = "3.8")

library(TSCAN)
load('OligTumor.RData')
library('Seurat')
source('scRef.R')




library(TSCAN)
library('Seurat')

#exp_raw_data=read.table('MGH36_mat.txt',header=T,row.names=1)
#exp_raw_data=read.table('MGH53_mat.txt',header=T,row.names=1)
#exp_raw_data=read.table('MGH54_mat.txt',header=T,row.names=1)
#exp_raw_data=read.table('MGH60_mat.txt',header=T,row.names=1)
#exp_raw_data=read.table('MGH93_mat.txt',header=T,row.names=1)
#exp_raw_data=read.table('MGH97_mat.txt',header=T,row.names=1)

pbmc <- CreateSeuratObject(raw.data =exp_raw_data, min.cells = 3, min.genes = 200, project = "10X_PBMC")

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]

exp_ref_mat=read.table('Reference_expression_human.txt',header=T,row.names=1,sep='\t',check.name=F)
a=SCREF(exp_sc_mat, exp_ref_mat)
out=a$out2
tag=a$tag2

rownames(out)[which(rownames(out)=='Newly Formed Oligodendrocyte')]='Oligodendrocytes'

tag[which(tag[,2]=="Newly Formed Oligodendrocyte"),2]='Oligodendrocytes'
tag[which(tag[,2]=="Oligodendrocyte Precursor Cell"),2]='OPC'

VAR_GENE=apply(exp_sc_mat,1,var)
USED_GENE=which(VAR_GENE>0)
used_data=exp_sc_mat[USED_GENE,]

#36
#procdata <- preprocess(used_data, clusternum=4,  pseudocount = 1, minexpr_value = 1, minexpr_percent = 0.4, cvcutoff = 0.5)
#53
#procdata <- preprocess(used_data, clusternum=4  )

#54
procdata <- preprocess(used_data,clusternum=4)
#60
#procdata <- preprocess(used_data, clusternum=4)
#93
#procdata <- preprocess(used_data, clusternum=4)
#97
#procdata <- preprocess(used_data, clusternum=4)

####
lpsmclust <- exprmclust(procdata, clusternum=4 )

tmp=names(lpsmclust$clusterid)
lpsmclust$clusterid=tag[,2]
names(lpsmclust$clusterid)=tmp

#pdf('TSCAN_MGH36.pdf',width=14,height=14)
#pdf('TSCAN_MGH53.pdf',width=14,height=14)
#pdf('TSCAN_MGH54.pdf',width=14,height=14)
#pdf('TSCAN_MGH60.pdf',width=14,height=14)
#pdf('TSCAN_MGH93.pdf',width=14,height=14)
#pdf('TSCAN_MGH97.pdf',width=14,height=14)
plotmclust(lpsmclust,show_cell_names=F,cell_name_size =1)
dev.off()

########################################################


########################################################

library(TSCAN)
pbmc=readRDS('GSE75330.RDS')
tmp=read.table('GSE75330_original_tag.txt',sep='\t')
tag=c()
for(one in names(pbmc@ident)){
tag=c(tag, as.character(tmp[which(as.character(tmp[,1])==one),2]))
}
#COP MFOL1 MFOL2  MOL1  MOL2  MOL3  MOL4  MOL5  MOL6 NFOL1 NFOL2   OPC   PPR
tag[which(tag %in% c('MFOL1','MFOL2'))]='MFOL'
tag[which(tag %in% c('MOL1','MOL2','MOL3','MOL4','MOL5','MOL6'))]='MOL'
tag[which(tag %in% c('NFOL1','NFOL2'))]='NFOL'


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

procdata <- preprocess(used_data, clusternum=6,  pseudocount = 1, minexpr_value = 1, minexpr_percent = 0.4, cvcutoff = 0.5)
lpsmclust <- exprmclust(procdata, clusternum =6)

tmp=names(lpsmclust$clusterid)
lpsmclust$clusterid=tag
names(lpsmclust$clusterid)=tmp
pdf('TSCAN_OligDEV.pdf',width=14,height=14)
plotmclust(lpsmclust,show_cell_names=F,cell_name_size =1)
dev.off()




