<img src="/Logo.png" width="80%">

# scREF

Reference-based single-cell annotation for single-cell RNA-seq data

# Citation:

Feng Zhang, Yaguang Dou, Weidong Tian; Reference-based single-cell annotation for single-cell RNA-seq data, Coming Soon

# Requirement

R: 3.5.0

Denpendent R packages: pcaPP

# Input format

Single-cell and reference expression matrix must be tab-delimited with header and row names.

For reference expression matrix, we recommend using TPM, FPKM, RPKM, or UMI matrix. 

For single-cell expression matrix, we recommend using UMI matrix.

# Usage

# 1. Reference-based annotation
    
    source('scRef.R')
    
    # Load data
    exp_sc_mat=read.table('Zeisel_exp_sc_mat.txt',header=T,row.names=1,'\t')
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
    out=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Multinomial.txt',quote=F,row.names=F,col.names=T,sep='\t')

    # Pearson
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='pearson',CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Pearson.txt',quote=F,row.names=F,col.names=T,sep='\t')

# 2. Semi-supervised annotation

    source('scRef.R')
    
    # First-round annotation - Kendall
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
    tag=.get_tag_max(out)
    
    # Build local reference
    LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)
    
    # Second-round annotation - Multinomial
    out=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Semi.txt',quote=F,row.names=F,col.names=T,sep='\t')

# 3. Compare tags or Combine the results of clustering method and scRef 
   
    source('scRef.R')
    
    # Compare tags (column 1: cell_id; column 2: labels)
    OUT=.compare_two_tag(TAG1, TAG2)
    write.table(OUT, 'COMPARE.txt', sep='\t', quote=F, col.names=T, row.names=F)

# 4. ScRef & Seurat

Original labels for Peripheral Blood Mononuclear Cells (PBMC):

https://satijalab.org/seurat/pbmc3k_tutorial.html

<img src="/source/Usage4_Original.png" width="500">

scRef labels (Reference: Gtex_v7_Tissue, [Download](/Reference/Gtex_v7_Tissue/)):

<img src="/source/Usage4_scRef.png" width="570">

scRef labels (Reference: BrainDev_AllenBrain, [Download](/Reference/BrainDev_AllenBrain/)):

Please note that this reference may not be suitable for PBMC. 

Smaller number indicates earlier development stage.

<img src="/source/Usage4_scRef_dev.png" width="450">

    library(Seurat)
    source('scRef.R')
    
    load('pbmc.RData')
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
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=10, print_step=10)
    tag=.get_tag_max(out)
    LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)
    out=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=10, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='TAG_GTEX.txt',quote=F,row.names=F,col.names=T,sep='\t')  
    TAG=tag
    old_ident = pbmc@ident
    pbmc@ident = as.factor(TAG[,2])
    names(pbmc@ident)=names(old_ident)
    pbmc_ori=pbmc
    TSNEPlot(object = pbmc, do.label=T, label.size=2.2)
    
    ##############DEV#################
    exp_ref_mat=read.table('GTEx_v7_median_tpm_human.txt',header=T,row.names=1,sep='\t',check.name=F)
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=10, print_step=10)
    tag=.get_tag_max(out)
    LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=10)
    out=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=10, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='TAG_DEV.txt',quote=F,row.names=F,col.names=T,sep='\t')
    TAG=tag
    old_ident = pbmc@ident
    pbmc@ident = as.factor(TAG[,2])
    names(pbmc@ident)=names(old_ident)
    pbmc_ori=pbmc
    COLOR=heat.colors(n=length(table(pbmc@ident))+2)
    TSNEPlot(object = pbmc, colors.use=COLOR)
    
    
    
     
# MIT License

    Copyright (c) 2018 Zhang, Feng

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





