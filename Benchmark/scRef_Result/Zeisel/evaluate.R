a=read.table('Spearman.txt',sep='\t',header=T)
b=read.table('../../Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
c=read.table('../../Zeisel_exp_sc_mat_cluster_merged.txt',sep='\t',header=T)

a[,2]=as.character(a[,2])
a[which(a[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocytes'
a[which(a[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocytes'

library('cidr')
library('NMI')


write.table(table(a[,2],b[,2]),file='TABLE.txt',sep='\t',quote=F,row.names=T,col.names=T)
#########Original###########
print('ori')
print('ARI')
adjustedRandIndex(a[,2],b[,2])
ids=1:length(a[,2])
pred=data.frame(ids=ids,label=a[,2])
real=data.frame(ids=ids,label=b[,2])
print('NMI')
NMI(real,pred)$value

#########Neuron Merged######
print('merge')
print('ARI')
adjustedRandIndex(a[,2],c[,2])
pred=data.frame(ids=ids,label=a[,2])
real=data.frame(ids=ids,label=c[,2])
print('NMI')
NMI(real,pred)$value


#########Neuron Removed#####
print('remove')
print('ARI')
r=which(!c[,2] %in% c('neurons'))
adjustedRandIndex(a[r,2],b[r,2])
ids=1:length(r)
pred=data.frame(ids=ids,label=a[r,2])
real=data.frame(ids=ids,label=b[r,2])
print('NMI')
NMI(real,pred)$value




P=a[,2]
R=c[,2]

getPR <-function(Rlabel,Plabel){
TP=length(which(R %in% Rlabel & P %in% Plabel))
PALL=length(which(P==Plabel))
RALL=length(which(R==Rlabel))
Pre=TP/PALL
Rec=TP/RALL
print('Pre')
print(Pre)
print('Rec')
print(Rec)}

Rlabel='astrocytes_ependymal'
Plabel='Astrocytes'

getPR('astrocytes_ependymal','Astrocytes')
getPR('endothelial-mural','Endothelial.Cells')
getPR('microglia','Microglia')
getPR('oligodendrocytes','Oligodendrocytes')
getPR('neurons','Neuron')









