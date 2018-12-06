library('mclust')
library('NMI')

a=read.table('TAG.txt',header=T,sep='\t')
b=read.table('../../Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
c=read.table('../../Zeisel_exp_sc_mat_cluster_merged.txt',sep='\t',header=T)



a[,2]=as.character(a[,2])
a[which(a[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocytes'
a[which(a[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocytes'


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
r=which(!c[,2] %in% c('OTHER'))
adjustedRandIndex(a[r,2],c[r,2])
pred=data.frame(ids=ids[r],label=a[r,2])
real=data.frame(ids=ids[r],label=c[r,2])
print('NMI')
NMI(real,pred)$value


#########Neuron Removed#####
print('remove')
print('ARI')
r=which(!c[,2] %in% c('neurons'))
adjustedRandIndex(a[r,2],b[r,2])
pred=data.frame(ids=ids[r],label=a[r,2])
real=data.frame(ids=ids[r],label=b[r,2])
print('NMI')
NMI(real,pred)$value
