
library(Seurat)
source('scRef.R')

exp_raw_data=read.table('GSE72056_melanoma_single_cell_revised_v2.txt.pure',header=T,row.names=1)

tag1=read.table('GSE72056_cnv.txt',header=T,row.names=1)
tag2=read.table('GSE72056_tag.txt',header=T,row.names=1)

tag1=t(tag1)
tag2=t(tag2)
tag1[which(tag1!='malignant')]=tag2[which(tag1!='malignant')]
tag1[which(tag1=='6')]='NK cell'
tag1[which(tag1=='0')]='unresolved'

pbmc <- CreateSeuratObject(raw.data = exp_raw_data, min.cells = 0, min.genes = 0, project = "10X_PBMC")


pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
    scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, 
    x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
    genes.print = 5)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:20, do.fast = TRUE)    
pbmc@meta.data$tag=tag1
saveRDS(pbmc, file = "GSE72056.RDS")

pdf('ORI.pdf',width=10,height=7)
TSNEPlot(object = pbmc, do.label=T,group.by='tag',cex=0.5)
dev.off()

######################

COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
exp_ref_mat=read.table('PeripheralBlood_ref_human.txt',header=T,row.names=1,sep='\t',check.name=F) 

REF_TAG=colnames(exp_ref_mat)
tmp=strsplit(REF_TAG, "_")
REF_TAG=c()
for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 

out=SCREF(exp_sc_mat, NewRef)

tag2=out$tag2
pbmc@meta.data$scref=tag2[,2]

pdf('SCREF.pdf',width=10,height=7)
TSNEPlot(object = pbmc, do.label=T, group.by ='scref', pt.size = 0.5)
dev.off()



merged_tag=pbmc@meta.data$tag
merged_tag[which(merged_tag=='NK cell')]='T & NK cell'
merged_tag[which(merged_tag=='T-cells')]='T & NK cell'





library(mclust)
#USED=which(pbmc@meta.data$tag!='malignant')

adjustedRandIndex(pbmc@meta.data$tag,tag2[,2])
#0.6123499
tmp=tag2
tmp[which(tag2[,2]=='NK cell'),2]='T & NK cell'
tmp[which(tag2[,2]=='T cell'),2]='T & NK cell'
adjustedRandIndex(merged_tag,tmp[,2])
#0.7508944


#Ken
adjustedRandIndex(pbmc@meta.data$tag, out$tag1[,2])
#0.5866281
tmp=out$tag1
tmp[which(tag2[,2]=='NK cell'),2]='T & NK cell'
tmp[which(tag2[,2]=='T cell'),2]='T & NK cell'
adjustedRandIndex(merged_tag,tmp[,2])
#0.7441506


sout=.get_cor(exp_sc_mat, NewRef, method='spearman',CPU=4, print_step=10)
stag=.get_tag_max(sout)
pbmc@meta.data$spea=stag[,2]
pdf('SPEA.pdf',width=10,height=7)
TSNEPlot(object = pbmc, do.label=T, group.by ='spea', pt.size = 0.5)
dev.off()
adjustedRandIndex(pbmc@meta.data$tag,stag[,2])
#0.6243437
tmp=stag
tmp[which(tag2[,2]=='NK cell'),2]='T & NK cell'
tmp[which(tag2[,2]=='T cell'),2]='T & NK cell'
adjustedRandIndex(merged_tag,tmp[,2])
#0.7297304


pout=.get_cor(exp_sc_mat, NewRef, method='pearson',CPU=4, print_step=10)
ptag=.get_tag_max(pout)
adjustedRandIndex(pbmc@meta.data$tag,ptag[,2])
#0.4173703



mout=.get_log_p_sc_given_ref(exp_sc_mat, NewRef, CPU=4, print_step=10)
mtag=.get_tag_max(mout)
adjustedRandIndex(pbmc@meta.data$tag,mtag[,2])
#0.1388988


