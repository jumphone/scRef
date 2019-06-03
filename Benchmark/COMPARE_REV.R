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
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL


#scPred Seurat 2
###################################
source('./scRNAseq_Benchmark-master/Scripts/run_scPred.R')
run_scPred(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/scPred_True_Labels.csv', './OUT/scPred_Pred_Labels.csv')
result$F1
result$MedF1

#scRef
###################################
source('./scRNAseq_Benchmark-master/Scripts/run_scRef.R')
run_scRef(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/scRef_True_Labels.csv', './OUT/scRef_Pred_Labels.csv')
result$F1
result$MedF1

#CHETAH R=3.6 Seurat3
###################################
setwd('F:/SCREF_COM')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL
#

source('./scRNAseq_Benchmark-master/Scripts/run_CHETAH.R')
run_CHETAH(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/CHETAH_True_Labels.csv', './OUT/CHETAH_Pred_Labels.csv')
result$F1
result$MedF1

#scID  Seurat 3
###################################

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL

source('./scRNAseq_Benchmark-master/Scripts/run_scID.R')
run_scID(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/scID_True_Labels.csv', './OUT/scID_Pred_Labels.csv')
result$F1
result$MedF1


#scmap  R 3.6
###################################

source('./scRNAseq_Benchmark-master/Scripts/run_scmap.R')
run_scmap(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/scmapcell_True_Labels.csv', './OUT/scmapcell_Pred_Labels.csv')
result$F1
result$MedF1

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/scmapcluster_True_Labels.csv', './OUT/scmapcluster_Pred_Labels.csv')
result$F1
result$MedF1


#SCINA  R 3.6 $CRASHED
###################################
setwd('F:/SCREF_COM')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/run_SCINA.R')
run_SCINA(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/SCINA_True_Labels.csv', './OUT/SCINA_Pred_Labels.csv')
result$F1
result$MedF1



#singleCellNet  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/run_singleCellNet.R')
run_singleCellNet(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/singleCellNet_True_Labels.csv', './OUT/singleCellNet_Pred_Labels.csv')
result$F1
result$MedF1



#garnett_CV  R 3.6 
#:set fileformat=unix 
#:w
###################################
setwd('F:/SCREF_COM')


DataPath='./Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/run_Garnett_CV.R')
run_Garnett_CV(DataPath, LabelsPath, CV_RDataPath, OutputDir, MarkerPath='MK.txt',Human=F)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/Garnett_CV_True_Labels.csv', './OUT/Garnett_CV_Pred_Labels.csv')
result$F1
result$MedF1



#run_CaSTLe  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='./Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/run_CaSTLe.R')
run_CaSTLe(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/True_Labels_CaSTLe.csv', './OUT/Pred_Labels_CaSTLe.csv')
result$F1
result$MedF1


#run_SingleR  R 3.6 
###################################
setwd('F:/SCREF_COM')

DataPath='./Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
OutputDir='./OUT/'
GeneOrderPath = NULL
NumGenes = NULL


source('./scRNAseq_Benchmark-master/Scripts/run_SingleR.R')
run_SingleR(DataPath, LabelsPath, CV_RDataPath, OutputDir)

setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/SingleR_True_Labels.csv', './OUT/SingleR_Pred_Labels.csv')
result$F1
result$MedF1









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
heatmap.2(RESULT,scale=c("none"),dendrogram='row',Colv=T,trace='none',col=colorRampPalette(c('blue','green','yellow','yellow3','red3')),margins=c(5,5))
heatmap.2(RESULT,scale=c("row"),dendrogram='row',Colv=T,trace='none',col=colorRampPalette(c('blue','red')),margins=c(5,5))


