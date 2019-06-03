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


#scPred
###################################
source('./scRNAseq_Benchmark-master/Scripts/run_scPred.R')
run_scPred(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/scPred_True_Labels.csv', './OUT/scPred_Pred_Labels.csv')
result$F1

#scRef
###################################
source('./scRNAseq_Benchmark-master/Scripts/run_scRef.R')
run_scRef(DataPath, LabelsPath, CV_RDataPath, OutputDir)
###################################
setwd('F:/SCREF_COM')
source('evaluate.R')
result <- evaluate('./OUT/scRef_True_Labels.csv', './OUT/scRef_Pred_Labels.csv')
result$F1

#CHETAH R=3.6
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


#scID V1.1 R3.6 
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
result$F1










