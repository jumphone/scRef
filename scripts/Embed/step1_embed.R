
library(Seurat)
source('scRef.R')
load('CASE.RObj')


TSNE_VEC=pbmc@dr$tsne@cell.embeddings
IDENT=as.character(pbmc@ident)



COL=c()
i=1
while(i <=length(pbmc@ident)){
    this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
    COL=c(COL,this_col)
    i=i+1
    } 

tag=cbind(names(pbmc@ident),pbmc@ident)
exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]
exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)


exp_sc_mat=read.table('WT.txt',sep='\t',check.name=F,row.names=1,header=T)

out=.get_cor(exp_sc_mat, LocalRef, method='kendall',CPU=10, print_step=10)
tag=.get_tag_max(out)

library(parallel)

set.seed(123456)

print_step=10
N=5
RANDOM=50

SINGLE = function(i){   
    library('pcaPP')
    .get_dis= function(this_sc, this_ref){
        exp_sc_mat=this_sc
        exp_ref_mat=this_ref
        exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
        exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
        gene_sc=rownames(exp_sc_mat)
        gene_ref=rownames(exp_ref_mat)
        gene_over= gene_sc[which(gene_sc %in% gene_ref)]
        exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
        exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
        colname_sc=colnames(exp_sc_mat)
        colname_ref=colnames(exp_ref_mat)
        log_p_sc_given_ref_list=c()
        exp_sc = as.array(exp_sc_mat[,1])
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            log_p_sc_given_ref=cor.fk(exp_sc,exp_ref)
            log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref) 
            j=j+1
            }
        return(log_p_sc_given_ref_list)
        }
      
    this_tag = as.character(tag[i,2])
    
    vec_index=which(IDENT==this_tag)
    vec_index=sample(vec_index, size=RANDOM, replace = FALSE) 
    
    this_vec = TSNE_VEC[vec_index,]
    this_ref= exp_ref_mat[,vec_index]
    this_sc = cbind(exp_sc_mat[,i],exp_sc_mat[,i])
    rownames(this_sc) = rownames(exp_sc_mat)
    colnames(this_sc)= c('rep1','rep2')
    this_out = .get_dis(this_sc, this_ref);
    this_out_rank=rank(-this_out)
    
    this_weight = rep(0,length(this_out))
    this_weight[which(this_out_rank <=N)]=1
    this_weight=this_weight/sum(this_weight)
    
    v1=sum(this_weight * this_vec[,1])
    v2=sum(this_weight * this_vec[,2])
    this_out_vec=c(v1,v2)
    if(i%%print_step==1){print(i)}
    return(this_out_vec)
    }

CPU=10
RUN = mclapply(1:length(tag[,1]), SINGLE, mc.cores=CPU)

VEC_OUT = c()
for(this_out_vec in RUN){
    VEC_OUT=cbind(VEC_OUT, this_out_vec)
    }
VEC_OUT=t(VEC_OUT)
rownames(VEC_OUT)=colnames(exp_sc_mat)


                                 
#library("RColorBrewer")
#display.brewer.all()
#brewer.pal(n=8,name='Set2')
#install.packages("wesanderson")

library(wesanderson)
colpal=wes_palette(name="Moonrise3", n=24, type = c("continuous"))
COLOR=colpal[as.factor(IDENT)]

XLIM=c(min(TSNE_VEC[,1]),max(TSNE_VEC[,1]))
YLIM=c(min(TSNE_VEC[,2]),max(TSNE_VEC[,2]))

pdf('tSNEinjection.pdf',width=15,height=10)
plot(TSNE_VEC, col=COLOR, pch=16,cex=0.5,xlim=XLIM, ylim=YLIM)
par(new=TRUE)
plot(VEC_OUT, col='red', pch=16,cex=0.5, xlim=XLIM, ylim=YLIM)
dev.off()



