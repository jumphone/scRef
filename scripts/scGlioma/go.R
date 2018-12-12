
exp_data=read.table('GSE104276_all_pfc_2394_UMI_count_NOERCC.xls',header=T,row.names=1,sep='\t')

TYPE=read.table('NPCneuglia.txt',header=T,sep='\t')


exp_raw_data=c()

COLNAME=colnames(exp_data)
i=1
for(cell in TYPE[,1]){
exp_raw_data=cbind(exp_raw_data,exp_data[,which(COLNAME==as.character(cell))])
#print(cell)
print(i);i=i+1
}

TMP=cbind(as.character(TYPE[,2]),as.character(TYPE[,1]))
NEWH=paste0(TMP[,1],'_',TMP[,2])

rownames(exp_raw_data)=rownames(exp_data)
colnames(exp_raw_data)=NEWH

write.table(exp_raw_data,file='GSE104276_withLabel.txt',sep='\t',quote=F,row.names=T,col.names=T)


###########################


library('Seurat')
pbmc <- CreateSeuratObject(raw.data = exp_raw_data, min.cells = 3, min.genes = 200, project = "10X_PBMC")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, do.plot=F,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

PCNUM=100
pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM, pcs.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 1:5, 
    genes.print = 5)

PCUSE=1:50
pbmc <- RunTSNE(object = pbmc, dims.use = PCUSE, do.fast = TRUE)


TSNEPlot(object = pbmc)

saveRDS(pbmc, "GSE104276.RDS")
##################################


source('scRef.R')


pdf('GENE_CUTOFF.pdf')
SYMBOL='PTPRC' #Microglia
LOC=.get_gene_cutoff(SYMBOL, pbmc@data)
PTPRC=.use_gene_cutoff(SYMBOL, pbmc@data, LOC[1])

SYMBOL='PAX6' #NPCs
LOC=.get_gene_cutoff(SYMBOL, pbmc@data)
PAX6=.use_gene_cutoff(SYMBOL, pbmc@data, LOC[1])

SYMBOL='PDGFRA' #OPCs
LOC=.get_gene_cutoff(SYMBOL, pbmc@data, BW=0.1)
PDGFRA=.use_gene_cutoff(SYMBOL, pbmc@data, LOC[1])

SYMBOL='NEUROD2' #Excitatory Neurons
LOC=.get_gene_cutoff(SYMBOL, pbmc@data)
NEUROD2=.use_gene_cutoff(SYMBOL, pbmc@data, LOC[1])

SYMBOL='GAD1' #Interneurons
LOC=.get_gene_cutoff(SYMBOL, pbmc@data,BW=0.1)
GAD1=.use_gene_cutoff(SYMBOL, pbmc@data, LOC[1])

SYMBOL='AQP4' #Astricytes
LOC=.get_gene_cutoff(SYMBOL, pbmc@data)
AQP4=.use_gene_cutoff(SYMBOL, pbmc@data, LOC[1])
dev.off()


#####
TAG=rep('NA',length(AQP4))
TAG[which(PTPRC==1)]=paste0(TAG[which(PTPRC==1)],'_','PTPRC')
TAG[which(PAX6==1)]=paste0(TAG[which(PAX6==1)],'_','PAX6')
TAG[which(PDGFRA==1)]=paste0(TAG[which(PDGFRA==1)],'_','PDGFRA')
TAG[which(NEUROD2==1)]=paste0(TAG[which(NEUROD2==1)],'_','NEUROD2')
TAG[which(GAD1==1)]=paste0(TAG[which(GAD1==1)],'_','GAD1')
TAG[which(AQP4==1)]=paste0(TAG[which(AQP4==1)],'_','AQP4')




pbmc@meta.data$geneoriref=TAG
TSNEPlot(pbmc, group.by='geneoriref')
GeneRef=.generate_ref(exp_raw_data, cbind(TAG,TAG) , min_cell = 10 ) 

GeneRef=GeneRef[,which(colnames(GeneRef)!='NA')]

OLDNAME=colnames(GeneRef)
tmp=strsplit(OLDNAME, "_")
NEWNAME=c()
for(one in tmp[]){NEWNAME=c(NEWNAME, paste(one[2: length(one) ], collapse='_')  ) }
colnames(GeneRef)=NEWNAME


out=.get_log_p_sc_given_ref(exp_raw_data, GeneRef)
result=.trajectory(out, plot_type='polygon', plot_size=1.7, label_dist=1.2, label_size=10, random_ratio=0.03)


png(filename = "CortexDev.png",width = 1024, height = 1024)
result$ggplot
dev.off()


tag=.get_tag_max(out)

MST=.generate_mst(result$mat, min_cell=5)
plot(MST)

pdf('MST.pdf')
plot(MST)
dev.off()

pbmc@meta.data$generef=tag[,2]
pdf('TSNE.pdf',width=10,height=6)
TSNEPlot(pbmc, group.by='generef')
TSNEPlot(pbmc,group.by='orig.ident')
dev.off()

save.image('GSE104276_Done.RData')





