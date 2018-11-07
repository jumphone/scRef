a=read.table('TAG.txt',sep='\t',header=T)
b=read.table('../../Hochgerner_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
c=read.table('../../Hochgerner_exp_sc_mat_cluster_merged.txt',sep='\t',header=T)

library('cidr')
library('NMI')



a[,2]=as.character(a[,2])
a[which(a[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocytes'
a[which(a[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocytes'
write.table(table(a[,2],b[,2]),file='TABLE.txt',sep='\t',quote=F,row.names=T,col.names=T)

#########Original###########
ids=1:length(a[,2])
pred=data.frame(ids=ids,label=a[,2])
real=data.frame(ids=ids,label=b[,2])
NMI(real,pred)$value
adjustedRandIndex(a[,2],b[,2])

#########Neuron Merged######
r=which(!c[,2] %in% c( 'OTHER'))
pred=data.frame(ids=ids[r],label=a[r,2])
real=data.frame(ids=ids[r],label=c[r,2])
NMI(real,pred)$value
adjustedRandIndex(a[r,2],c[r,2])



#table(a[,2],b[,2])
#write.table(table(a[,2],b[,2]),file='TABLE.txt',sep='\t',row.names=T,col.names=T,quote=F)



##########Neuron Removed#####
r=which(!c[,2] %in% c('Neuron', 'OTHER'))
#r=which(!c[,2] %in% c( 'OTHER'))
ids=1:length(r)
pred=data.frame(ids=ids,label=a[r,2])
real=data.frame(ids=ids,label=c[r,2])
NMI(real,pred)$value
adjustedRandIndex(a[r,2],c[r,2])

