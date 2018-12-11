
library('Seurat')
source('scRef.R')

exp_raw_data= read.table('MGH54_mat.txt',sep='\t',header=T,row.names=1)
#exp_raw_data= read.table('GSE70630_OG_processed_data_v2_MGH53.txt',sep='\t',header=T,row.names=1)
#exp_raw_data= read.table('GSE70630_OG_processed_data_v2_MGH60.txt',sep='\t',header=T,row.names=1)
#exp_raw_data= read.table('GSE70630_OG_processed_data_v2_MGH93.txt',sep='\t',header=T,row.names=1)
#exp_raw_data= read.table('GSE70630_OG_processed_data_v2_MGH97.txt',sep='\t',header=T,row.names=1)
#exp_raw_data= read.table('GSE70630_OG_processed_data_v2_MGH36.txt',sep='\t',header=T,row.names=1)

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
tag[which(tag[,2]=='Newly Formed Oligodendrocyte'),2]='Oligodendrocytes'


######STEMNESS##
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",  scale.factor = 10000)

exp_mat=as.matrix(pbmc@data)
STEM=read.table('STEM.txt')
OC=read.table('OC.txt')
AC=read.table('AC.txt')

STEM_ROW=which(rownames(exp_mat) %in% STEM[,1])
OC_ROW=which(rownames(exp_mat) %in% OC[,1])
AC_ROW=which(rownames(exp_mat) %in% AC[,1])

STEM_SCORE=apply(exp_mat[STEM_ROW,], 2,mean)
OC_SCORE=apply(exp_mat[OC_ROW,], 2,mean)
AC_SCORE=apply(exp_mat[AC_ROW,], 2,mean)

AO_SCORE=cbind(OC_SCORE,AC_SCORE)
MAX_AO_SCORE=apply(AO_SCORE,1,max)

FINAL_STEM_SCORE=STEM_SCORE-MAX_AO_SCORE
BX_USE=which(tag[,2]!='Microglia')
boxplot(FINAL_STEM_SCORE[BX_USE]~tag[BX_USE,2],out=F)

AC_STEM=FINAL_STEM_SCORE[which(tag[,2]=='Astrocytes')]
OC_STEM=FINAL_STEM_SCORE[which(tag[,2]=='Oligodendrocytes')]
OPC_STEM=FINAL_STEM_SCORE[which(tag[,2]=='Oligodendrocyte Precursor Cell')]

ks.test(OPC_STEM,AC_STEM)
#0.00131
ks.test(OPC_STEM,OC_STEM)
#1.067e-13
######



rownames(out)[which(rownames(out)=='Newly Formed Oligodendrocyte')]='Oligodendrocytes'
result=.trajectory(out, plot_type='polygon', plot_size=1.7, cell_size=2,label_dist=1.2, label_size=10, random_ratio=0.03)

png(filename = "TraOligGliomaMGH54.png",width = 1024, height = 1024)
#png(filename = "TraOligGliomaMGH53.png",width = 1024, height = 1024)
#png(filename = "TraOligGliomaMGH60.png",width = 1024, height = 1024)
#png(filename = "TraOligGliomaMGH93.png",width = 1024, height = 1024)
#png(filename = "TraOligGliomaMGH97.png",width = 1024, height = 1024)
#png(filename = "TraOligGliomaMGH36.png",width = 1024, height = 1024)
result$ggplot
dev.off()

library(igraph)

MST=.generate_mst(result$mat)
pdf('MST.pdf')
plot(MST)
plot(MST,edge.width=2, vertex.size=1,vertex.label.color = "black")
dev.off()


pbmc=readRDS('GSE75330.RDS')
exp_ref_mat=read.table('GSE75330_mouse_combined_reference.txt',header=T,row.names=1,sep='\t')

# Reference-based annotation - multinomial
out=.get_log_p_sc_given_ref(pbmc@raw.data, exp_ref_mat)

# Construct trajectory
result=.trajectory(out, plot_type='polygon', plot_size=1.7, cell_size=2,label_dist=1.2, label_size=10, random_ratio=0.03)
png(filename = "TraOPC.png",width = 1024, height = 1024)
result$ggplot
dev.off()

MST=.generate_mst(result$mat)
plot(MST,edge.width=2, vertex.size=1,vertex.label.color = "black")

