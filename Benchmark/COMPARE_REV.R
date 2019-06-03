#a=read.table('Zeisel_exp_sc_mat.txt',header=T,row.names=1)
#write.table(t(a),'Zeisel_exp_sc_mat.txt.csv',quote=F,sep=',',row.names=T,col.names=T)
#a=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,row.names=1,sep='\t')
#write.table(a,'Zeisel_exp_sc_mat_cluster_original.txt.csv',quote=F,sep=',',row.names=T,col.names=T)
#install.packages('devtools')
#devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))
#install.packages("./", repos=NULL, type="source",dep=T)

#source('Cross_Validation.R')
#Cross_Validation('./Zeisel_exp_sc_mat_cluster_original.txt.csv',1,'./OUT/')
#setwd('../')

source('evaluate.R')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV/'
GeneOrderPath = NULL
NumGenes = NULL


#scPred Seurat 2
###################################
source('./scRNAseq_Benchmark-master/Scripts/rev/run_scPred.R')
run_scPred(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/scPred_True_Labels.csv', './OUT_REV/scPred_Pred_Labels.csv')
result$F1
result$MedF1

#scRef
###################################
source('./scRNAseq_Benchmark-master/Scripts/rev/run_scRef.R')
run_scRef(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/scRef_True_Labels.csv', './OUT_REV/scRef_Pred_Labels.csv')
result$F1
result$MedF1

#CHETAH R=3.6 Seurat3
###################################
setwd('F:/SCREF_COM')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV/'
GeneOrderPath = NULL
NumGenes = NULL
#

source('./scRNAseq_Benchmark-master/Scripts/rev/run_CHETAH.R')
run_CHETAH(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/CHETAH_True_Labels.csv', './OUT_REV/CHETAH_Pred_Labels.csv')
result$F1
result$MedF1


#scID  Seurat 3
###################################

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV/'
GeneOrderPath = NULL
NumGenes = NULL

source('./scRNAseq_Benchmark-master/Scripts/rev/run_scID.R')
run_scID(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/scID_True_Labels.csv', './OUT_REV/scID_Pred_Labels.csv')
result$F1
result$MedF1


#scmap  R 3.6
###################################

source('./scRNAseq_Benchmark-master/Scripts/rev/run_scmap.R')
run_scmap(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/scmapcell_True_Labels.csv', './OUT_REV/scmapcell_Pred_Labels.csv')
result$F1
result$MedF1

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/scmapcluster_True_Labels.csv', './OUT_REV/scmapcluster_Pred_Labels.csv')
result$F1
result$MedF1





#singleCellNet  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/rev/run_singleCellNet.R')
run_singleCellNet(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/singleCellNet_True_Labels.csv', './OUT_REV/singleCellNet_Pred_Labels.csv')
result$F1
result$MedF1






#run_CaSTLe  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='./Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/rev/run_CaSTLe.R')
run_CaSTLe(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/True_Labels_CaSTLe.csv', './OUT_REV/Pred_Labels_CaSTLe.csv')
result$F1
result$MedF1


#run_SingleR  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='./Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT_REV/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/rev/run_SingleR.R')
run_SingleR(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT_REV/SingleR_True_Labels.csv', './OUT_REV/SingleR_Pred_Labels.csv')
result$F1
result$MedF1









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
heatmap.2(RESULT,scale=c("none"),dendrogram='row',Colv=T,trace='none',col=colorRampPalette(c('blue','yellow','yellow3','red3')),margins=c(5,5))

rRESULT=t(apply(RESULT,1,rank))
rownames(rRESULT)=rownames(RESULT)
colnames(rRESULT)=colnames(RESULT)


heatmap.2(rRESULT,scale=c("row"),dendrogram='row',Colv=T,trace='none',col=colorRampPalette(c('blue','red')),margins=c(5,5))


