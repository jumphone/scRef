
load('TSNE.RData')
library(Seurat)


TSNE_VEC=pbmc@dr$tsne@cell.embeddings
TAG=tag


library(ClusterR)

MAXC=5
OUT_TAG=TAG
tag_list=TAG[,2]
uniq_tag_list=unique(tag_list)

i=1
while(i<=length(uniq_tag_list))
    this_uniq_tag=uniq_tag_list[i]
    used_cell=which(tag_list==this_uniq_tag)
    this_tsne_vec=TSNE_VEC[used_cell,]
    this_tag=tag_list[used_cell]
    
    this_opt=Optimal_Clusters_KMeans(this_tsne_vec,max_clusters =MAXC,plot_clusters =FALSE)


    i=i+1
    }



