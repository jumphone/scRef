library('Seurat')
exp_mat=read.table('GSE92332_atlas_UMIcounts_reheader.txt',sep='\t',header=T,row.names=1)
pbmc <- CreateSeuratObject(raw.data = exp_mat, min.cells = 0, min.genes = 0, project = "10X_PBMC")
tag=cbind(as.character(pbmc@meta.data$orig.ident),as.character(pbmc@meta.data$orig.ident))
source('scRef.R')
NewRef=.generate_ref(exp_mat,tag)
write.table(NewRef, file='GSE92332_intestin_mouse_ref.txt',quote=F,row.names=T,col.names=T,sep='\t')


pbmc=readRDS('KO.RDS')

table(pbmc@meta.data$MCASmallIntest)

#Epithelial.cell
#Columnar.epithelium.Small.Intestine.
#Epithelium.of.small.intestinal.villi 


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
exp_sc_mat=exp_sc_mat[,which(pbmc@meta.data$MCASmallIntest %in% c('Epithelial.cell','Columnar.epithelium.Small.Intestine.','Epithelium.of.small.intestinal.villi'))]
out=SCREF(exp_sc_mat, NewRef, min_cell=1)
table(out$tag2[,2])
ko_out=out


pbmc=readRDS('WT.RDS')
table(pbmc@meta.data$MCASmallIntest)

#Epithelial.cell
#Columnar.epithelium.Small.Intestine.
#Epithelium.of.small.intestinal.villi 


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    }      
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
exp_sc_mat=exp_sc_mat[,which(pbmc@meta.data$MCASmallIntest %in% c('Epithelial.cell','Columnar.epithelium.Small.Intestine.','Epithelium.of.small.intestinal.villi'))]
out=SCREF(exp_sc_mat, NewRef, min_cell=1)
table(out$tag2[,2])
wt_out=out


table(wt_out$tag2[,2])


pie(table(wt_out$tag2[,2]))
