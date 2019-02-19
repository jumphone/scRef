   
# 5. Using aggregated single-cell data as reference

Download single-cell references: [link](https://github.com/jumphone/scRef/tree/master/Reference)

If you are using aggregated single-cell data as reference, please check the number of expressed genes of each column (reference), and you can remove some columns with less expressed genes than a given cutoff:
    
    source('scRef.R')
    NUMBER=.check_pos(exp_ref_mat)
    CUTOFF=15000
    exp_ref_mat=exp_ref_mat[,which(NUMBER < CUTOFF)] 

or you can trim the number of expressed genes (keep top genes) of all columns to a specific number:

    source('scRef.R')
    exp_ref_mat=.trim_pos(exp_ref_mat, 12000) 
