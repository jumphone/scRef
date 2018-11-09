# scREF: reference-based single-cell annotation for single-cell RNA-seq data

![image](/Logo.png "image")

# Citation:

Feng Zhang, Yaguang Dou, Weidong Tian; Reference-based single-cell annotation for single-cell RNA-seq data, Coming Soon

# Requirement

R: 3.5.0

Denpendent R packages: pcaPP

# Usage

# 1. Plain reference-based annotation
    
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

# 2. Semi-supervised reference-based annotation

    source('scRef.R')
    
    # First-round annotation - Kendall
    out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
    tag=.get_tag_max(out)
    
    # Build local reference
    LocalRef=.generate_ref(exp_sc_mat, tag, min_cell=1)
    
    # Second-round annotation - Multinomial
    out=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=4, print_step=10)
    tag=.get_tag_max(out)
    write.table(tag,file='Semi.txt',quote=F,row.names=F,col.names=T,sep='\t')

# 3. Compare tags or Combine the results of tSNE and scRef 
   
    source('scRef.R')
    
    #Compare tags (column 1: cell_id; column 2: labels)
    OUT=.compare_two_tag(TAG1, TAG2)
    write.table(OUT, 'COMPARE.txt', sep='\t', quote=F, col.names=T, row.names=F)



