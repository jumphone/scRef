
a=read.table('Kendall.txt',header=T,sep='\t')
a[,2]=as.character(a[,2])
a[which(a[,2]=='Oligodendrocyte.Precursor.Cell'),2]='Oligodendrocyte'
a[which(a[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocyte'
a[which(a[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocyte'
K=a


a=read.table('Multinomial.txt',header=T,sep='\t')
a[,2]=as.character(a[,2])
a[which(a[,2]=='Oligodendrocyte.Precursor.Cell'),2]='Oligodendrocyte'
a[which(a[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocyte'
a[which(a[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocyte'
M=a

a=read.table('Pearson.txt',header=T,sep='\t')
a[,2]=as.character(a[,2])
a[which(a[,2]=='Oligodendrocyte.Precursor.Cell'),2]='Oligodendrocyte'
a[which(a[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocyte'
a[which(a[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocyte'
P=a

a=read.table('Spearman.txt',header=T,sep='\t')
a[,2]=as.character(a[,2])
a[which(a[,2]=='Oligodendrocyte.Precursor.Cell'),2]='Oligodendrocyte'
a[which(a[,2]=='Myelinating.Oligodendrocytes'),2]='Oligodendrocyte'
a[which(a[,2]=='Newly.Formed.Oligodendrocyte'),2]='Oligodendrocyte'
S=a


R=read.table('../../Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
R[,2]=as.character(R[,2])

R[which(R[,2]=='astrocytes_ependymal'),2]='Astrocytes'
R[which(R[,2]=='endothelial-mural'),2]='Endothelial.Cells'
R[which(R[,2]=='interneurons'),2]='Neuron'
R[which(R[,2]=='microglia'),2]='Microglia'
R[which(R[,2]=='oligodendrocytes'),2]='Oligodendrocyte'
R[which(R[,2]=='pyramidal CA1'),2]='Neuron'
R[which(R[,2]=='pyramidal SS'),2]='Neuron'



a=rep(0,length(R[,2]))
a[which(K[,2]==R[,2])]=1
KP=a

a=rep(0,length(R[,2]))
a[which(M[,2]==R[,2])]=1
MP=a

a=rep(0,length(R[,2]))
a[which(S[,2]==R[,2])]=1
SP=a

a=rep(0,length(R[,2]))
a[which(P[,2]==R[,2])]=1
PP=a

ERR=which(PP==0 & MP==0 &SP==0 & KP==0)
ONE=which(PP==1 | MP==1 | SP==1 | KP==1)

FLAG=rep(0,length(R[,2]))
FLAG[ONE]=1
write.table(FLAG, 'ERR_FLAG.txt', sep='\t',quote=F,row.names=F,col.names=F)



load('this.RData')


getGN=function(X){
    
	X[which(X>0)]=1
	return(sum(X))
}

GN=apply(exp_sc_mat,2,getGN)
#MEAN=apply(exp_sc_mat, 2, mean)


#boxplot(MEAN~FLAG,xlab='0: Wrong in all four methods. 1: Correct in at least one method', ylab='MEAN of UMI for one cell')
boxplot(GN~FLAG,xlab='0: Wrong in all four methods. 1: Correct in at least one method', ylab='# Sequenced gene for one cell')










