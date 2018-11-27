#!/home/zhangfeng/bin/bin/Rscript
#######################################################################################
# Description: Run pagoda2 for single cell expression data
# Date: 2018-10-24 22:20
# douyaguang
# package required: 
#    pagoda2:         basic pagoda2 algorithm package
#    NMI:             calculate NMI
#    cidr:            calculate ARI
#    igraph:          for cluster methods used by pagoda2
#######################################################################################

usage <- "Usage: 
    program  <input_data_path>  <output_data_path>  <core_number>
        input_data_path:     the input single cell data path

        output_data_path:    the output cluster result file prefix, 
                             the output will be two files, 
                             one is cluster result by infomap method that ends with \"_infomap\"
                             another one is cluster result by multilevel method that ends with
                             \"_multilevel\"

        core_number:         running core numbers (default: 10)
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
cat("core_number:\t", cores_number, "\n", sep = "", file = time_info_file, append = F)

seed <- 123
set.seed(seed)

library(pagoda2)
library(Matrix)

single_cell_data <- read.table(input_path, row.names = 1, header = T, sep = "\t", check.names = F)

time_start <- as.numeric(Sys.time())

counts <- gene.vs.molecule.cell.filter(single_cell_data, min.cell.size = 500)
r <- Pagoda2$new(as.matrix(counts), log.scale = TRUE, n.cores = cores_number)
r$adjustVariance(plot = F, gam.k = 10)
r$calculatePcaReduction(nPcs = 50, n.odgenes = 3e3)
r$makeKnnGraph(k = 40, type = 'PCA', center = T, distance = "cosine")

time_end <- as.numeric(Sys.time())
time_used <- time_end - time_start
cat(time_used, sep = "\n", file = time_info_file, append = T)

library(igraph)
# infomap cluster method

time_start <- as.numeric(Sys.time())

r$getKnnClusters(method = infomap.community, type = "PCA")
infomap_result <- r$clusters$PCA$community

time_end <- as.numeric(Sys.time())
time_used <- time_end - time_start
cat(time_used, sep = "\n", file = time_info_file, append = T)


# multilevel cluster method
time_start <- as.numeric(Sys.time())

r$getKnnClusters(method = multilevel.community, type = "PCA")
multilevel_result <- r$clusters$PCA$community

time_end <- as.numeric(Sys.time())
time_used <- time_end - time_start
cat(time_used, sep = "\n", file = time_info_file, append = T)

comb_results <- rbind(infomap_result, multilevel_result)
rownames(comb_results) <- c("pogoda2_infomap", "pogoda2_multilevel")
write.table(file = output_path, t(comb_results), quote = F, sep = "\t")
