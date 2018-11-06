
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

Kendall_correct=which(KP==1)
Multinomial_correct=which(MP==1)
Spearman_correct=which(SP==1)
Pearson_correct=which(PP==1)

library(VennDiagram)

venn.diagram(x=list(Kendall=Kendall_correct, Multinomial=Multinomial_correct,
Spearman=Spearman_correct, Pearson=Pearson_correct), 
"VENN.png", height = 750, width = 750, 
resolution =300, imagetype="png", 
alpha=c(0.6, 0.6, 0.6, 0.6),lwd=0.5, cex=0.5,cat.cex=0.5)



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



###########Other CON########################

CON=which(K[,2]==P[,2] & P[,2]==S[,2] & S[,2] == M[,2])

#CON precision
2521/length(CON)
#0.9714836

#KEN precision
length(which(K[,2]==R[,2]))/3005
#0.9341098


a=which(K[,2]=='Astrocytes')
Astrocytes=a[which(a %in% CON)]

a=which(K[,2]=='Endothelial.Cells')
Endothelial.Cells=a[which(a %in% CON)]

a=which(K[,2]=='Microglia')
Microglia=a[which(a %in% CON)]

a=which(K[,2]=='Neuron')
Neuron=a[which(a %in% CON)]

a=which(K[,2]=='Oligodendrocyte')
Oligodendrocyte=a[which(a %in% CON)]


Astrocytes = apply(exp_sc_mat[,Astrocytes],1,sum)
Endothelial.Cells = apply(exp_sc_mat[,Endothelial.Cells],1,sum)
Microglia=apply(exp_sc_mat[,Microglia],1,sum)
Neuron=apply(exp_sc_mat[,Neuron],1,sum)
Oligodendrocyte=apply(exp_sc_mat[,Oligodendrocyte],1,sum)
NewRef=cbind(Astrocytes,Endothelial.Cells,Microglia,Neuron,Oligodendrocyte)

source('scRef.R')
out=.get_cor(exp_sc_mat, NewRef, method='kendall',CPU=4, print_step=10 )
tag=.get_tag_max(out)
length(which(tag[,2]==R[,2]))/3005
#0.9580699

source('scRef.R')
out=.get_cor(exp_sc_mat, NewRef, method='pearson',CPU=4, print_step=10 )
tag=.get_tag_max(out)
tag[CON,2]=K[CON,2]
length(which(tag[,2]==R[,2]))/3005
#0.909817

source('scRef.R')
out=.get_cor(exp_sc_mat, NewRef, method='spearman',CPU=4, print_step=10 )
tag=.get_tag_max(out)
length(which(tag[,2]==R[,2]))/3005
#0.9554077

source('scRef.R')
out=.get_log_p_sc_given_ref(exp_sc_mat, NewRef,CPU=4, print_step=10 )
tag=.get_tag_max(out)
length(which(tag[,2]==R[,2]))/3005
#0.9600666


