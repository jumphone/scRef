library('mclust')
library('NMI')

TAG=read.table('TAG.txt',header=T,sep='\t')
TAG[,2]=as.character(TAG[,2])
TAG[which(TAG[,2]=='Myelinating Oligodendrocytes'),2]='Oligodendrocyte'
TAG[which(TAG[,2]=='Newly Formed Oligodendrocyte'),2]='Oligodendrocyte'

ORI=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
ORI[,2]=as.character(ORI[,2])
MER=read.table('Zeisel_exp_sc_mat_cluster_merged.txt',header=T,sep='\t')
MER[,2]=as.character(MER[,2])

r=which(MER[,2]!='neurons')

ids=c(1:length(ORI[,2]))
TAG_N=data.frame(cbind(ids,TAG[,2]))
ORI_N=data.frame(cbind(ids,ORI[,2]))
MER_N=data.frame(cbind(ids,MER[,2]))


O_ARI=adjustedRandIndex(TAG[,2],ORI[,2])
O_NMI=NMI(TAG_N,ORI_N)

M_ARI=adjustedRandIndex(TAG[,2],MER[,2])
M_NMI=NMI(TAG_N,MER_N)

print('O_ARI:')
print(O_ARI)
print('O_NMI:')
print(O_NMI)
print('M_ARI:')
print(M_ARI)
print('M_NMI:')
print(M_NMI)

TAG_N=data.frame(cbind(ids[r],TAG[r,2]))
ORI_N=data.frame(cbind(ids[r],ORI[r,2]))
MER_N=data.frame(cbind(ids[r],MER[r,2]))

R_ARI=adjustedRandIndex(TAG[r,2],MER[r,2])
R_NMI=NMI(TAG_N,MER_N)

print('R_ARI:')
print(R_ARI)
print('R_NMI:')
print(R_NMI)
