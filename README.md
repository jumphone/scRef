# scRef
A reference-based toolkit for single-cell analysis

# dbRef
A web resource of references. <a href='/Reference/'>Visit dbRef</a>

# Table of content

* [Citation](#Citation)
* [Download](#Download)
* [Requirement](#Requirement)
* [Input format](#Input-format)
* [Usage](#Usage)
    * [1. Single-round classification](#1-single-round-classification)
    * [2. scRef-two-round classification (RBC)](#2-SCREF-two-round-classification-RBC)
    * [3. Compare the results of clustering method and scRef](#3-Compare-the-results-of-clustering-method-and-scRef)
    * [4. scRef & Seurat](#4-scRef--Seurat)
    * [5. scRef-tSNE plot alignment (TPA)](#5-scRef-tSNE-plot-alignment-TPA)
    * [6. scRef-trajectory detection (TDE)](#6-scRef-trajectory-detection-TDE)
* [License](#License)


# Citation

Feng Zhang, Yaguang Dou, Rohit R Rao, Q. Richard Lu, Weidong Tian; SCREF: a reference-based solution for single-cell analysis, Coming Soon

# Download

### Whole repository:

Download: https://github.com/jumphone/scRef/archive/master.zip

### Only "scRef.R":

For Linux and Mac: 

curl https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R > scRef.R

For Windows (in cmd):

curl https://raw.githubusercontent.com/jumphone/scRef/master/windows_version/scRef.R > scRef.R

# Requirement

R: 3.5.0

Denpendent R packages: pcaPP, MASS, ggplot2, igraph, pastecs

# Input format

Single-cell and reference expression matrix must be tab-delimited with header and row names.

For reference expression matrix, we recommend using TPM, FPKM, RPKM, RPK, or UMI matrix. 

For single-cell expression matrix, we recommend using UMI matrix.

# Usage

# 1. Single-round classification

### Code:

    #install.packages('pcaPP')    
    source('scRef.R')
    
    # Load data
    exp_sc_mat=read.table('Zeisel_exp_sc_mat.txt',header=T,row.names=1,sep='\t')
    exp_ref_mat=read.table('Reference_expression.txt',header=T,row.names=1,sep='\t')
    
    # Kendall
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Kendall.txt',quote=F,row.names=F,col.names=T,sep='\t')
    
    # Spearman
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='spearman',CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Spearman.txt',quote=F,row.names=F,col.names=T,sep='\t')

    # Multinomial
    out=.get_log_p_sc_given_ref(exp_sc_mat, exp_ref_mat, CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Multinomial.txt',quote=F,row.names=F,col.names=T,sep='\t')

    # Pearson
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='pearson',CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Pearson.txt',quote=F,row.names=F,col.names=T,sep='\t')

# 2. scRef two-round classification (RBC)

   Users can directly use "SCREF" funtion to do semi-supervised (two-round) annotation.
   
### Code:

    #install.packages('pcaPP') 
    source('scRef.R')
    
    # Load data
    exp_sc_mat=read.table('Zeisel_exp_sc_mat.txt',header=T,row.names=1,sep='\t')
    exp_ref_mat=read.table('Reference_expression.txt',header=T,row.names=1,sep='\t')
   
    # SCREF annotation
    tag=SCREF(exp_sc_mat, exp_ref_mat, CPU=4, print_step=10)$tag2
    
    # Output results
    write.table(tag,file='Semi.txt',quote=F,row.names=F,col.names=T,sep='\t')
    
### OR:
    
    source('scRef.R')
    
    # Load data
    exp_sc_mat=read.table('Zeisel_exp_sc_mat.txt',header=T,row.names=1,sep='\t')
    exp_ref_mat=read.table('Reference_expression.txt',header=T,row.names=1,sep='\t')
   
    # First-round annotation - Kendall
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
    tag=.get_tag_max(out)  
    
    # Build local reference
    LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)  
    
    # Second-round annotation - Multinomial
    out=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=4, print_step=10)
    tag=.get_tag_max(out)
    
    # Output results
    write.table(tag,file='Semi.txt',quote=F,row.names=F,col.names=T,sep='\t')

# 3. Compare the results of clustering method and scRef 

### Step 1.

Use “Seurat” package in R to do the tSNE:

    pbmc = RunTSNE(object = pbmc, dims.use = 1:150, do.fast = TRUE, dim.embed = 2)
    pbmc_3 = RunTSNE(object = pbmc, dims.use = 1:150, do.fast = TRUE, dim.embed = 3)

### Step 2.

Use the tSNE vectors to do hierarchical (or kmeans clustering by using “kmeans” in R) clustering. Users can choose a proper “k” according to the tSNE plot: 

    TSNE_VEC=pbmc_3@dr$tsne@cell.embeddings
    D=dist(TSNE_VEC)
    H=hclust(D)
    C=cutree(H, k=8) 
    pbmc@meta.data$C=C
    TSNEPlot(object = pbmc, do.label=T, group.by ='C')

### Step 3.

Compare the clusters given by hierarchical clustering and the annotation given by scRef. Users can choose a threshold for the “max_over” to get final result (we suggest using a value around 0.8). Users should be careful about the annotation result when there are multiple hits, because the multiple hits might be caused by unknown cell types. 

    source(“scRef.R”)
    TAG1=cbind(rownames(TSNE_VEC),C)
    TAG2=read.table('Zeisel_semi.txt', header=T, sep='\t')
    OUT= .compare_two_tag(TAG1, TAG2)
    write.table(OUT, 'COMPARE.txt', sep='\t', quote=F, col.names=T, row.names=F)
 
### Format of “OUT” (output of “.compare_two_tag”):

    Column 1: Max(Column3, Column5)
    Column 2: a cluster label of TAG1
    Column 3: number of overlapped cells divided by number of cells with TAG1’s label 
    Column 4: a cluster label of TAG2
    Column 5: number of overlapped cells divided by number of cells with TAG2’s label 
    

# 4. scRef & Seurat

### Original labels for human Peripheral Blood Mononuclear Cells (PBMC) 

Source: https://satijalab.org/seurat/pbmc3k_tutorial.html

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage4_Original.png" width="520">
</a>

### scRef labels (Reference: HumanTissue_Bulk_GtexV7, [Download](/Reference/HumanTissue_Bulk_GtexV7/)) 

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage4_scRef.png" width="570">
</a>

### scRef labels (Reference: HumanBrainDev_Bulk_AllenBrain, [Download](/Reference/HumanBrainDev_Bulk_AllenBrain/)) 

Please note that this brain reference may not be suitable for PBMC. 

Smaller number indicates earlier development stage.

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage4_scRef_dev.png" width="450">
</a>

### Code:

    library(Seurat)
    load('pbmc.RData')
    source('scRef.R')
    
    # To get Seurat object "pbmc", please follow the instruction of Seurat:
    # https://satijalab.org/seurat/pbmc3k_tutorial.html 
        
    COL=c()
    i=1
    while(i <=length(pbmc@ident)){
        this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
        COL=c(COL,this_col)
	    i=i+1
        }      
    exp_sc_mat=as.matrix(pbmc@raw.data)[,COL]
    
    #############GTEX############
    exp_ref_mat=read.table('GTEx_v7_median_tpm_human.txt',header=T,row.names=1,sep='\t',check.name=F)    
    tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
    pbmc@meta.data$scref=tag[,2]
    TSNEPlot(object = pbmc, do.label=T, label.size=2.2, group.by ='scref')
    
    ##############DEV#################
    exp_ref_mat=read.table('exp_ref_mat_human_brain_dev',header=T,row.names=1,sep='\t',check.name=F)
    tag=SCREF(exp_sc_mat, exp_ref_mat)$tag2
    pbmc@meta.data$scref=tag[,2]
    COLOR=heat.colors(n=length(table(pbmc@meta.data$scref))+2)
    TSNEPlot(object = pbmc, colors.use=COLOR, group.by ='scref')
    
    ##################################      

# 5. scRef tSNE plot alignment (TPA)

Visualize scRNA-seq data of different batches in a single dimension reduction (e.g. tSNE) plot .

Users can use this function to align cells to any given dimension reduction plot (e.g. tSNE plot, PCA plot, and UMAP plot, etc.)

Cell types of the low-quality data should be covered by the high-quality data.

### Workflow: 

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage5.png" width="420">
</a>

### Result:

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage5_TsneProjection.png" width="420">
</a>

### Code:

    library(Seurat)
    load('pbmc.RData')
    source('scRef.R')
     
    ref_vec=pbmc@dr$tsne@cell.embeddings
    COL=c()
    i=1
    while(i <=length(pbmc@ident)){
        this_col=which(colnames(pbmc@raw.data)==names(pbmc@ident)[i])
        COL=c(COL,this_col)
        i=i+1
        } 
	
    ref_tag=cbind(names(pbmc@ident),as.character(pbmc@ident))    
    exp_ref_mat=as.matrix(pbmc@raw.data)[,COL]
    exp_sc_mat=exp_ref_mat[,which(ref_tag[,2]=='5')]        
    
    # Generate local reference of high-quality data
    LocalRef= .generate_ref(exp_ref_mat, ref_tag, min_cell = 10 )    
    
    # Semi-supervised annotation
    sc_tag=SCREF(exp_sc_mat, LocalRef)$tag2
    
    # tSNE projection
    out =.vec_projection(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec, 
            method='kendall', nearest_cell=3, alpha=0.5, random_size=30, 
            random_seed=123, min_cell=10, CPU=4, print_step=10)
    
    XLIM=c(min(ref_vec[,1]-1),max(ref_vec[,1]+1))
    YLIM=c(min(ref_vec[,2]-1),max(ref_vec[,2]+1))

    plot(ref_vec,xlim=XLIM,ylim=YLIM,pch=16,col='grey70')
    par(new=T)
    plot(out$vec,xlim=XLIM,ylim=YLIM,pch=16,col='red')


# 6. scRef trajectory detection (TDE)

Drawing trajectory based on the results of scRef

Note: the input should have at least three cell types

### Demo data (GSE75330_Oligodendrocyte): 

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage6_original.png" width="420">
</a>

### Our result (GSE75330_Oligodendrocyte): 

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage6.png" width="420">
</a>

### Minimum Spanning Tree (MST) (GSE75330_Oligodendrocyte): 

Distance Calculation:

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage6_distance.png" width="240">
</a>

Result:

<a href='https://github.com/jumphone/scRef'>
<img src="/source/Usage6_MST.png" width="420">
</a>

### Code:
     
    #install.packages('MASS')
    #install.packages('ggplot2')
    
    library(Seurat)
    source('scRef.R')
    
    pbmc=readRDS('GSE75330.RDS')
    exp_ref_mat=read.table('GSE75330_mouse_combined_reference.txt',header=T,row.names=1,sep='\t')
    
    # Reference-based annotation - multinomial
    out=.get_log_p_sc_given_ref(pbmc@raw.data, exp_ref_mat)
    
    # Construct trajectory
    result=.trajectory(out, plot_type='polygon', plot_size=1.7, label_dist=1.2, label_size=10, random_ratio=0.03)

    png(filename = "TraOPC.png",width = 1024, height = 1024)
    result$ggplot
    dev.off()
    
    # MST:
    #install.packages('igraph')
    library(igraph)
    
    MST=.generate_mst(result$mat)
    pdf('MST.pdf')
    plot(MST)
    dev.off()
    
# License

    Copyright (c) 2019 Zhang, Feng

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.





