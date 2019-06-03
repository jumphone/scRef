
source('evaluate.R')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV_NOI/'
GeneOrderPath = NULL
NumGenes = NULL


#scPred Seurat 2
###################################
source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_scPred.R')
run_scPred(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/scPred_True_Labels.csv', './OUT_REV_NOI/scPred_Pred_Labels.csv')
result$F1
result$MedF1

#scRef
###################################
source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_scRef.R')
run_scRef(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/scRef_True_Labels.csv', './OUT_REV_NOI/scRef_Pred_Labels.csv')
result$F1
result$MedF1

#CHETAH R=3.6 Seurat3
###################################
setwd('F:/SCREF_COM')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV_NOI/'
GeneOrderPath = NULL
NumGenes = NULL
#

source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_CHETAH.R')
run_CHETAH(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/CHETAH_True_Labels.csv', './OUT_REV_NOI/CHETAH_Pred_Labels.csv')
result$F1
result$MedF1


#scID  Seurat 3
###################################

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV_NOI/'
GeneOrderPath = NULL
NumGenes = NULL

source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_scID.R')
run_scID(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/scID_True_Labels.csv', './OUT_REV_NOI/scID_Pred_Labels.csv')
result$F1
result$MedF1


#scmap  R 3.6
###################################

source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_scmap.R')
run_scmap(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/scmapcell_True_Labels.csv', './OUT_REV_NOI/scmapcell_Pred_Labels.csv')
result$F1
result$MedF1

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/scmapcluster_True_Labels.csv', './OUT_REV_NOI/scmapcluster_Pred_Labels.csv')
result$F1
result$MedF1





#singleCellNet  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV_NOI/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_singleCellNet.R')
run_singleCellNet(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/singleCellNet_True_Labels.csv', './OUT_REV_NOI/singleCellNet_Pred_Labels.csv')
result$F1
result$MedF1






#run_CaSTLe  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='./Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV_NOI/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_CaSTLe.R')
run_CaSTLe(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/True_Labels_CaSTLe.csv', './OUT_REV_NOI/Pred_Labels_CaSTLe.csv')
result$F1
result$MedF1


#run_SingleR  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='./Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV_NOI/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/revnoi/run_SingleR.R')
run_SingleR(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV_NOI/SingleR_True_Labels.csv', './OUT_REV_NOI/SingleR_Pred_Labels.csv')
result$F1
result$MedF1









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
heatmap.2(RESULT,scale=c("none"),dendrogram='row',Colv=T,trace='none',col=colorRampPalette(c('blue','yellow','yellow3','red3')),margins=c(5,5))

rRESULT=t(apply(RESULT,1,rank))
rownames(rRESULT)=rownames(RESULT)
colnames(rRESULT)=colnames(RESULT)

SSS=apply(RESULT,2,sum)
library(gplots)
heatmap.2(rRESULT[,order(SSS)],scale=c("none"),key=F,cellnote=round(RESULT[,order(SSS)],2),
notecol='black',dendrogram='row',Colv=F,Rowv=T, trace='none',col=colorRampPalette(c('grey90','pink','red3')),margins=c(5,5))

