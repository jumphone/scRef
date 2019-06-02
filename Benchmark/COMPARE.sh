a=read.table('Zeisel_exp_sc_mat.txt',header=T,row.names=1)
write.table(t(a),'Zeisel_exp_sc_mat.txt.csv',quote=F,sep=',',row.names=T,col.names=T)
a=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,row.names=1,sep='\t')
write.table(a,'Zeisel_exp_sc_mat_cluster_original.txt.csv',quote=F,sep=',',row.names=T,col.names=T)

devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))

install.packages("./", repos=NULL, type="source",dep=T)

source('Cross_Validation.R')
Cross_Validation('./Zeisel_exp_sc_mat_cluster_original.txt.csv',1,'./OUT/')
setwd('../')

source('./scRNAseq_Benchmark-master/Scripts/run_scPred.R')
run_scPred('Zeisel_exp_sc_mat.txt.csv','./Zeisel_exp_sc_mat_cluster_original.txt.csv','./OUT/CV_folds.RData','~/OUT/')


DataPath='Zeisel_exp_sc_mat.txt.csv'
LabelsPath='./Zeisel_exp_sc_mat_cluster_original.txt.csv'
CV_RDataPath='./OUT/CV_folds.RData'
GeneOrderPath = NULL
NumGenes = NULL


  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  print(dim(Data))
  print(dim(Labels))

  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }

