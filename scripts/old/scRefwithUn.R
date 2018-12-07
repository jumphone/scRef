
load('TSNE.RData')
library(Seurat)

TAG=tag

    TSNE_VEC=pbmc@dr$tsne@cell.embeddings


.add_clustering <-function(TAG, TSNE_VEC, MAXC=10, random_seed=123){

    random_seed=random_seed
    TSNE_VEC=TSNE_VEC
    TAG=as.matrix(TAG)
    TAG[,1]=as.character(TAG[,1])
    TAG[,2]=as.character(TAG[,2])
    library(factoextra)
    set.seed(random_seed)
    MAXC=MAXC
    OUT_TAG=TAG
    tag_list=TAG[,2]
    uniq_tag_list=unique(tag_list)

    i=1
    while(i<=length(uniq_tag_list)){
        this_uniq_tag=uniq_tag_list[i]
        print(this_uniq_tag)
        used_cell=which(tag_list==this_uniq_tag)
        this_tsne_vec=TSNE_VEC[used_cell,]
        this_tag=tag_list[used_cell]
    
        this_opt=fviz_nbclust(this_tsne_vec, kmeans, k.max = MAXC,method = "silhouette")
        this_k = as.numeric(this_opt$data[which(this_opt$data[,2] == max(this_opt$data[,2])),1])
    
        if(this_k==2){this_k=1}
        this_out=kmeans(this_tsne_vec,centers=this_k)
    
        this_new_tag=as.character(this_out$cluster)
    
        j=1
        while(j<=length(used_cell)){
            OUT_TAG[used_cell,2][j]=paste0(OUT_TAG[used_cell,2][j],'_',this_new_tag[j])
            j=j+1}
    
        i=i+1
        }
    return(OUT_TAG)
    }



pbmc@meta.data$new=OUT_TAG[,2]
TSNEPlot(pbmc,group.by='new')


