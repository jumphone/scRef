#http://dropviz.org/

anno=readRDS('annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS')

exp=readRDS('metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS')

COL=which(anno[,1]=='CB')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60Cerebellum_ALT_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='FC')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60Cortex_noRep5_FRONTALonly_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='PC')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60Cortex_noRep5_POSTERIORonly_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='ENT')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60EntoPeduncular_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='GP')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60GlobusPallidus_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='HC')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60Hippocampus_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='STR')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60Striatum_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='SN')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60SubstantiaNigra_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)

COL=which(anno[,1]=='TH')
out=exp[,COL]
colnames(out)= anno[COL,5]
write.table(out, file='GRCm38.81.P60Thalamus_mouse.txt', quote=F, sep='\t',row.names=T,col.names=T)









