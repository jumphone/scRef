# Notice:

When you use MCA as a reference, we recommend combining all single-gene-expression-based subtypes.

For example, merging "Astrocyte_Pla2g7 high(Brain)", "Astrocyte_Atp1b2 high(Brain)", and "Astrocyte_Mfe8 high(Brain)" into "Astrocyte".

You can use the following scripts to generate a new reference by combining those subtypes:

    source('scRef.R')
    exp_ref_mat=read.table('MCA_Brain_ref_mouse.txt',header=T,sep='\t',row.names=1)
    REF_TAG=colnames(exp_ref_mat)
    tmp=strsplit(REF_TAG, "_")
    REF_TAG=c()
    for(one in tmp){REF_TAG=c(REF_TAG, one[1])}
    NewRef=.generate_ref(exp_ref_mat, cbind(REF_TAG,REF_TAG), min_cell=1) 

We provide the combined file which includes all main cell types in MCA:

MCA_combined_mouse.txt.zip and MCA_combined_human.txt.zip



# Data Type:

Mouse, UMI

# Data Source:

This reference is built from Mouse Cell Atlas (MCA) single-cell RNA-seq Data

http://bis.zju.edu.cn/MCA/

https://figshare.com/articles/MCA_DGE_Data/5435866

# Workflow:

* For each gene in each cell type, we add all UMI up.

* Detailed scripts are in: https://github.com/jumphone/scRef/tree/master/scripts/MCA

* We have generated a human reference by mapping mouse gene to human gene.



