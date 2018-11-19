load('TSNE.RData')
library('Seurat')
source('scRef.R')
ori_label=read.table('Zeisel_exp_sc_mat_cluster_original.txt',header=T,sep='\t')
pbmc@meta.data$ori=ori_label[,2]
#TSNEPlot(pbmc, group.by="ori")
table(pbmc@meta.data$ori)
ref_vec=pbmc@dr$tsne@cell.embeddings

OLIG=which(pbmc@meta.data$ori=='oligodendrocytes')

plot(ref_vec[,1],ref_vec[,2], pch=16, col='gray80', xlim=c(-30, 30), ylim=c(-40,35), xlab='tSNE_1',ylab='tSNE_2')
par(new=T)
plot(ref_vec[OLIG,1],ref_vec[OLIG,2], pch=16, col='red', xlim=c(-30, 30), ylim=c(-40,35), xlab='',ylab='')


COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

ref_tag=cbind(names(pbmc@ident), as.character(pbmc@meta.data$ori))    
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]


exp_sc_mat= exp_ref_mat[,OLIG]

getRanGene <- function(X){
    POS = which(X >0 )
    N=length(POS)/2
    KEEP = sample(x=POS, size=N )
    NEG = POS[which(!POS %in% KEEP)]
    X[NEG]=0
    return(X)
    }

set.seed(123)
sim_exp_sc_mat = apply(exp_sc_mat,2, getRanGene)


LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )
out=.get_cor(sim_exp_sc_mat, LocalRef, method='kendall',CPU=4, print_step=10)
tag=.get_tag_max(out)
LocalRef= .generate_ref(sim_exp_sc_mat, tag, min_cell = 10 )
out=.get_log_p_sc_given_ref(sim_exp_sc_mat, LocalRef, CPU=10, print_step=10)
tag=.get_tag_max(out)
sc_tag=tag

out =.vec_projection(exp_sc_mat=sim_exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
        method='kendall', nearest_cell=3, alpha=0.5, random_size=30, 
        random_seed=123, min_cell=10, CPU=4, print_step=10)

pdf('simulationresult_tSNEprojection.pdf',width=7, height=7)
CEX=0.5
plot(ref_vec,xlim=c(-30, 30), ylim=c(-40,35),pch=16,col='grey70', cex=CEX)
par(new=T)
plot(ref_vec[OLIG,1],ref_vec[OLIG,2],xlim=c(-30, 30), ylim=c(-40,35),pch=16,col='blue', cex=CEX, xlab='',ylab='')
par(new=T)
plot(out$vec,xlim=c(-30, 30), ylim=c(-40,35),pch=16,col='red', cex=CEX)
dev.off()







