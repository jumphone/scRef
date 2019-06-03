
##################################
RESULT=c()
setwd('F:/SCREF_COM')
source('evaluate.R')

result <- evaluate('./OUT/scRef_True_Labels.csv', './OUT/scRef_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT/scPred_True_Labels.csv', './OUT/scPred_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT/CHETAH_True_Labels.csv', './OUT/CHETAH_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT/scID_True_Labels.csv', './OUT/scID_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT/scmapcell_True_Labels.csv', './OUT/scmapcell_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT/scmapcluster_True_Labels.csv', './OUT/scmapcluster_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT/singleCellNet_True_Labels.csv', './OUT/singleCellNet_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT/True_Labels_CaSTLe.csv', './OUT/Pred_Labels_CaSTLe.csv')
RESULT=cbind(RESULT,result$F1)


result <- evaluate('./OUT/SingleR_True_Labels.csv', './OUT/SingleR_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)


colnames(RESULT)=c('scRef','scPred','CHETAH','scID', 'scmapcell', 'scmapcluster', 'singleCellNet','CaSTLe','SingleR')

library(gplots)
rRESULT=t(apply(RESULT,1,rank))
rownames(rRESULT)=rownames(RESULT)
colnames(rRESULT)=colnames(RESULT)

SSS=apply(RESULT,2,sum)
library(gplots)
heatmap.2(rRESULT[,order(SSS)],scale=c("none"),key=F,cellnote=round(RESULT[,order(SSS)],2),
notecol='black',dendrogram='row',Colv=F,Rowv=T, trace='none',col=colorRampPalette(c('grey90','pink','red3')),margins=c(5,5))
######################################################################################################



##################################
RESULT=c()
setwd('F:/SCREF_COM')
source('evaluate.R')


result <- evaluate('./OUT_REV/scRef_True_Labels.csv', './OUT_REV/scRef_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV/scPred_True_Labels.csv', './OUT_REV/scPred_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV/CHETAH_True_Labels.csv', './OUT_REV/CHETAH_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV/scID_True_Labels.csv', './OUT_REV/scID_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV/scmapcell_True_Labels.csv', './OUT_REV/scmapcell_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV/scmapcluster_True_Labels.csv', './OUT_REV/scmapcluster_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV/singleCellNet_True_Labels.csv', './OUT_REV/singleCellNet_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV/True_Labels_CaSTLe.csv', './OUT_REV/Pred_Labels_CaSTLe.csv')
RESULT=cbind(RESULT,result$F1)


result <- evaluate('./OUT_REV/SingleR_True_Labels.csv', './OUT_REV/SingleR_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

colnames(RESULT)=c('scRef','scPred','CHETAH','scID', 'scmapcell', 'scmapcluster', 'singleCellNet','CaSTLe','SingleR')

library(gplots)
rRESULT=t(apply(RESULT,1,rank))
rownames(rRESULT)=rownames(RESULT)
colnames(rRESULT)=colnames(RESULT)

SSS=apply(RESULT,2,sum)
library(gplots)
heatmap.2(rRESULT[,order(SSS)],scale=c("none"),key=F,cellnote=round(RESULT[,order(SSS)],2),
notecol='black',dendrogram='row',Colv=F,Rowv=T, trace='none',col=colorRampPalette(c('grey90','pink','red3')),margins=c(5,5))


##############################



##################################
RESULT=c()
setwd('F:/SCREF_COM')
source('evaluate.R')


result <- evaluate('./OUT_REV_NOI/scRef_True_Labels.csv', './OUT_REV_NOI/scRef_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV_NOI/scPred_True_Labels.csv', './OUT_REV_NOI/scPred_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV_NOI/CHETAH_True_Labels.csv', './OUT_REV_NOI/CHETAH_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV_NOI/scID_True_Labels.csv', './OUT_REV_NOI/scID_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV_NOI/scmapcell_True_Labels.csv', './OUT_REV_NOI/scmapcell_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV_NOI/scmapcluster_True_Labels.csv', './OUT_REV_NOI/scmapcluster_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV_NOI/singleCellNet_True_Labels.csv', './OUT_REV_NOI/singleCellNet_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

result <- evaluate('./OUT_REV_NOI/True_Labels_CaSTLe.csv', './OUT_REV_NOI/Pred_Labels_CaSTLe.csv')
RESULT=cbind(RESULT,result$F1)


result <- evaluate('./OUT_REV_NOI/SingleR_True_Labels.csv', './OUT_REV_NOI/SingleR_Pred_Labels.csv')
RESULT=cbind(RESULT,result$F1)

colnames(RESULT)=c('scRef','scPred','CHETAH','scID', 'scmapcell', 'scmapcluster', 'singleCellNet','CaSTLe','SingleR')

library(gplots)
rRESULT=t(apply(RESULT,1,rank))
rownames(rRESULT)=rownames(RESULT)
colnames(rRESULT)=colnames(RESULT)

SSS=apply(RESULT,2,sum)
library(gplots)
heatmap.2(rRESULT[,order(SSS)],scale=c("none"),key=F,cellnote=round(RESULT[,order(SSS)],2),
notecol='black',dendrogram='row',Colv=F,Rowv=T, trace='none',col=colorRampPalette(c('grey90','pink','red3')),margins=c(5,5))



