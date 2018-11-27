#!/home/zhangfeng/bin/bin/Rscript
#######################################################################################
# Description: Run SAFE and all the methods compared in SAFE for single cell data
# Date: 2018-10-22 14:40
# douyaguang
# package required: 
#    SAFEclustering:  for SAFEclustering algorithm
#    NMI:             calculate NMI
#    cidr:            calculate ARI
#######################################################################################

usage <- "Usage: 
    program  <input_data_path>  <output_data_path>  <core_number>
        input_data_path:     the input single cell data path

        output_data_path:    the output cluster result file
"
args <- commandArgs(T)
if (TRUE %in% (c("-h", "-H", "--help") %in% args))  {
	message(usage)
	quit(save = "default", status = 0, runLast = TRUE)
}
input_path <- args[1]
output_path <- args[2]
cores_number <- 10
if (length(args) > 2) cores_number <- as.numeric(args[3])
time_info_file <- paste(output_path, "time", sep = ".")
cat("", file = time_info_file)

# SAFEclustering can not be installed directly by github, many dependent R packages
# must be installed via BioConductor or github manually, after all the dependent R
# packages have been installed, the the SAFEclustering can be install by instruction 
# from github.
library("SAFEclustering")

home <- "/home/disk/deepexp/"
# default seed in individual_clustering toturial
seed <- 123 
methods <- c("SC3", "CIDR", "Seurat", "tSNE", "SAFE")

## loading dataset
single_cell_data <- read.table(input_path, row.names = 1, header = T, sep = "\t", check.names = F)

## run SC3 CIDR t-SNE+k-means methods compared by SAFE, The methods exec order is
## SC3 -> CIDR -> Seurat -> tSNE for all the methods assigned TRUE.

time_start <- as.numeric(Sys.time())

single_cell_cluster_results <- individual_clustering(inputTags = single_cell_data,
	datatype = "count", mt_filter = FALSE, nGene_filter = FALSE, 
	SC3 = TRUE, gene_filter = FALSE, 
	CIDR = TRUE, nPC.cidr = NULL, 
	Seurat = TRUE, nPC.seurat = NULL, resolution = 0.9, 
	tSNE = TRUE, dimensions = 3, perplexity = 30, 
	SEED = seed,
	time_info_file = time_info_file)

time_end <- as.numeric(Sys.time())
time_used <- time_end - time_start
cat("other_methods:\t", time_used, "\n", sep = "", file = time_info_file, append = T)

# SAFE program part. Ensemble learning based on the clustering result of the above methods.
# Graph-partition algorithm is applied while perform ensemble learning. The Graph-partition
# exec has already downloaded and compiled in the path:
# "../5CompareMethods/SAFEclustering/single_cell_clustering"
program_dir <- paste(home, "5CompareMethods/SAFEclustering/single_cell_clustering", sep = "/")

time_start <- as.numeric(Sys.time())

single_cell_cluster_ensemble <- SAFE(cluster_results = single_cell_cluster_results, 
	program.dir = program_dir, MCLA = TRUE, CSPA = TRUE, HGPA = TRUE, SEED = seed)

time_end <- as.numeric(Sys.time())
time_used <- time_end - time_start
cat("SAFE:\t", time_used, "\n", sep = "", file = time_info_file, append = T)

## combine all the results together
comb_results <- rbind(single_cell_cluster_results, single_cell_cluster_ensemble$optimal_clustering)
rownames(comb_results) <- methods
colnames(comb_results) <- colnames(single_cell_data)
write.table(file = output_path, t(comb_results), quote = F, sep = "\t")
