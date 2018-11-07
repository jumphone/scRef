
a=read.table('rank.txt',sep='\t',header=T,row.names=1)
a=as.matrix(a)
library('gplots')

heatmap.2(a,Rowv=F,Colv=F, dendrogram='none',col=colorpanel(256, low="Blue", high="Red", mid="grey90"), scale="row", key=TRUE, keysize=1, symkey=FALSE, densadj = 0.1, density.info="none", trace="none",margin=c(20,20),cexRow=2,cexCol=2)
