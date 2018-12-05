
library('Seurat')
source('scRef.R')

exp_raw_data= read.table('GSE70630_OG_processed_data_v2_MGH54.txt',sep='\t',header=T,row.names=1)
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
plot(MST)

