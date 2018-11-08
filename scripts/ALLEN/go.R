
exp_data=read.table('exp_mat.txt',sep='\t',row.names=1,header=F)
meta_data=read.table('meta.txt',sep='\t',header=T)

library(pcaPP)

################################################################
Stage=meta_data[,9]
OUT=c()
i=1
while(i<=length(exp_data[,1])){
this_gene=rownames(exp_data)[i]
this_exp=t(exp_data[i,])
this_cor=cor.fk(this_exp, Stage)
OUT=cbind(OUT,c(this_gene,this_cor))
i=i+1
print(i)
}
OUT=t(OUT)
OUT[which(OUT[,2]=='NaN'),2]='0'
hist(as.numeric(OUT[,2]),breaks=100)
save(OUT,file='COR.RData')
################################################################
COR=as.numeric(OUT[,2])
STEM = which(COR < -0.5)
STEM_GENES=OUT[STEM,1]
save(STEM_GENES,file='STEM.RData')
################################################################





