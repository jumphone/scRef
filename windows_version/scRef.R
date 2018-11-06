#library('pcaPP')
########################################
#Author: Feng Zhang
#Email: 15110700005@fudan.edu.cn
#######################################
.get_log_p_sc_given_ref <- function(exp_sc_mat, exp_ref_mat, CPU=4, print_step=10){
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name   
    ##################
    library(parallel)
    ##################
    Refprob <- function(exp_sc, exp_ref){
    log_p_sc_given_ref = dmultinom(x=exp_sc,log=T,prob=exp_ref)
    return(log_p_sc_given_ref)}
    #################
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    #Step 2. calculate prob
    SINGLE <- function(i){
        exp_sc = as.array(exp_sc_mat[,i])
        log_p_sc_given_ref_list=c()
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            #####
            exp_ref[which(exp_ref==0)]=min(exp_ref[which(exp_ref>0)])
            #####
            log_p_sc_given_ref=Refprob(exp_sc,exp_ref)
            log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
            j=j+1}
        ################################
        #if(i%%print_step==1){cat(i);cat('\n')}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    cl= makeCluster(CPU)
    RUN = parLapply(cl=cl,1:length(colname_sc), SINGLE)
    LOG_P_SC_GIVEN_REF = c()
    for(log_p_sc_given_ref_list in RUN){
        LOG_P_SC_GIVEN_REF=cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)}
    #######################################    
    rownames(LOG_P_SC_GIVEN_REF)=colname_ref
    colnames(LOG_P_SC_GIVEN_REF)=colname_sc
    return(LOG_P_SC_GIVEN_REF)
    }

.get_cor  <- function(exp_sc_mat, exp_ref_mat, method='pearson',CPU=4, print_step=10){
    #method = "pearson", "kendall", "spearman"
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
    ##################
    library(parallel)
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    #Step 2. calculate prob
    SINGLE <- function(i){
        library('pcaPP')
        exp_sc = as.array(exp_sc_mat[,i])
        log_p_sc_given_ref_list=c()
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            #####
            exp_ref[which(exp_ref==0)]=min(exp_ref[which(exp_ref>0)])
            #####
            #if(method=='rococo'){log_p_sc_given_ref=rococo(exp_sc,exp_ref)} else
            if(method=='kendall'){log_p_sc_given_ref=cor.fk(exp_sc,exp_ref)}
            else{
            log_p_sc_given_ref=cor(exp_sc,exp_ref, method=method)}
            log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
            j=j+1}
        ################################
        #if(i%%print_step==1){cat(i);cat('\n')}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    cl= makeCluster(CPU)
    RUN = parLapply(cl=cl,1:length(colname_sc), SINGLE)
    LOG_P_SC_GIVEN_REF = c()
    for(log_p_sc_given_ref_list in RUN){
        LOG_P_SC_GIVEN_REF=cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)}
    #######################################
    rownames(LOG_P_SC_GIVEN_REF)=colname_ref
    colnames(LOG_P_SC_GIVEN_REF)=colname_sc
    return(LOG_P_SC_GIVEN_REF)
    }


.get_p_ref_given_sc <- function(LOG_P_SC_GIVEN_REF){
    #######
    ref2sc <- function(log_p_sc_given_ref_list){
        weight_list = exp(log_p_sc_given_ref_list - max(log_p_sc_given_ref_list))
        p_ref_given_sc = weight_list / sum(weight_list)
        return(p_ref_given_sc)
        }
    #######
    P_REF_GIVEN_SC = apply(LOG_P_SC_GIVEN_REF, 2, ref2sc)
    colnames(P_REF_GIVEN_SC)=colnames(LOG_P_SC_GIVEN_REF)
    rownames(P_REF_GIVEN_SC)=rownames(LOG_P_SC_GIVEN_REF)
    return(P_REF_GIVEN_SC)
    }


.get_tag_max <- function(P_REF_GIVEN_SC){
    RN=rownames(P_REF_GIVEN_SC)
    CN=colnames(P_REF_GIVEN_SC)
    TAG=cbind(CN,rep('NA',length(CN)))
    i=1
    while(i<=length(CN)){
        this_rn_index=which(P_REF_GIVEN_SC[,i] == max(P_REF_GIVEN_SC[,i]))[1]
        TAG[i,2]=RN[this_rn_index]
        i=i+1
        }
    colnames(TAG)=c('cell_id','tag')
    return(TAG)
    }


.get_tag_min <- function(P_REF_GIVEN_SC){
    RN=rownames(P_REF_GIVEN_SC)
    CN=colnames(P_REF_GIVEN_SC)
    TAG=cbind(CN,rep('NA',length(CN)))
    i=1
    while(i<=length(CN)){
        this_rn_index=which(P_REF_GIVEN_SC[,i] == min(P_REF_GIVEN_SC[,i]))[1]
        TAG[i,2]=RN[this_rn_index]
        i=i+1
        }
    colnames(TAG)=c('cell_id','tag')
    return(TAG)
    }


.tag_iteration <- function(exp_sc_mat, TAG, method='multinomial', CPU=4, print_step=10){
	NewRef=c()
	TAG[,2]=as.character(TAG[,2])
    refnames=names(table(TAG[,2]))
    for(one in refnames){
    	this_col=which(TAG[,2]==one)
    	if(length(this_col)>=1){
        this_new_ref=apply(exp_sc_mat[,this_col],1,sum)        
        NewRef=cbind(NewRef,this_new_ref)}
    }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=refnames
    if(method=='multinomial'){
    OUT=.get_log_p_sc_given_ref(exp_sc_mat, NewRef, CPU=CPU, print_step=print_step)}
    else {OUT=.get_cor(exp_sc_mat, NewRef, method=method, CPU=CPU, print_step=print_step)}}
    TAG=.get_tag_max(OUT)
    return(TAG)
    }


