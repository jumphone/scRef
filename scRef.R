#library('pcaPP')
########################################
#Author: Feng Zhang
#Email: 15110700005@fudan.edu.cn
#######################################

.get_log_p_sc_given_ref <- function(exp_sc_mat, exp_ref_mat, CPU=4, print_step=10){
    delta = 0.01;
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
            exp_ref[which(exp_ref==0)]= delta * min(exp_ref[which(exp_ref>0)])
            #####
            log_p_sc_given_ref=Refprob(exp_sc,exp_ref)
            log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
            j=j+1}
        ################################
        if(i%%print_step==1){print(i)}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF = c()
    for(log_p_sc_given_ref_list in RUN){
        LOG_P_SC_GIVEN_REF=cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)}
    #######################################    
    rownames(LOG_P_SC_GIVEN_REF)=colname_ref
    colnames(LOG_P_SC_GIVEN_REF)=colname_sc
    return(LOG_P_SC_GIVEN_REF)
    }

.get_cor  <- function(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10){
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
            #exp_ref[which(exp_ref==0)]=min(exp_ref[which(exp_ref>0)])
            #####
            #if(method=='rococo'){log_p_sc_given_ref=rococo(exp_sc,exp_ref)} else
            if(method=='kendall'){log_p_sc_given_ref=cor.fk(exp_sc,exp_ref)}
            else{
            log_p_sc_given_ref=cor(exp_sc,exp_ref, method=method)}
            log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
            j=j+1}
        ################################
        if(i%%print_step==1){print(i)}
        return(log_p_sc_given_ref_list)
        }
    #######################################
    RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
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


.generate_ref <- function(exp_sc_mat, TAG, min_cell=1, refnames=FALSE){
    NewRef=c()
    TAG[,2]=as.character(TAG[,2])
    if(refnames==FALSE){
        refnames=names(table(TAG[,2]))}
        else{refnames=refnames}
    outnames=c()
    for(one in refnames){
        this_col=which(TAG[,2]==one)
        if(length(this_col)>= min_cell){
            outnames=c(outnames,one)
            if(length(this_col) >1){
                this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
                }
                else{this_new_ref = exp_sc_mat[,this_col]}
            NewRef=cbind(NewRef,this_new_ref)
            }
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }



.compare_two_tag <- function(TAG1, TAG2){
    OUT=c()
    tag1_names=as.character(unique(TAG1[,2]))
    tag2_names=as.character(unique(TAG2[,2]))
    i=1
    while(i<=length(tag1_names)){
        tag1 = tag1_names[i]
        tag1_index = which(TAG1[,2]== tag1)
        j=1
        while(j<=length(tag2_names)){
            tag2 = tag2_names[j] 
            #print(tag2)
            tag2_index = which(TAG2[,2]== tag2)
            over = length(which(tag1_index %in% tag2_index))
            tag1_over = over/length(tag1_index)
            tag2_over = over/length(tag2_index)
            max_over = max(tag1_over, tag2_over)
            OUT=cbind(OUT, c(max_over, tag1, tag1_over ,tag2, tag2_over)) 
            j=j+1
            }
        i=i+1
        }
    OUT=t(OUT)
    #OUT=as.matrix(OUT)
    OUT[,1]=as.numeric(OUT[,1])
    OUT[,3]=as.numeric(OUT[,3])
    OUT[,5]=as.numeric(OUT[,5])
    colnames(OUT)=c('max_over','tag1','tag1_over','tag2','tag2_over')
    OUT=OUT[order(OUT[,1],decreasing=T),]
    return(OUT)
    }




.vec_projection <- function(exp_sc_mat, sc_tag, exp_ref_mat, ref_tag, ref_vec,  method='kendall', nearest_cell=3, random_size=30, random_seed=123, alpha=0.5, min_cell=10, CPU=4, print_step=10){

    delta = 0.01;
    alpha = alpha;
    library(parallel);
    random_seed=random_seed;
    set.seed(random_seed);
    sc_cell_name=colnames(exp_sc_mat);
    tag=sc_tag;
    ref_tag=ref_tag;
    ref_vec=ref_vec;
    exp_ref_mat=exp_ref_mat;
    exp_sc_mat=exp_sc_mat;
    method=method;
    nearest_cell=nearest_cell;
    random_size=random_size;
    min_cell=min_cell;
    CPU=CPU;
    print_step=print_step;
    
    SINGLE = function(i){   
        library('pcaPP')
        Refprob <- function(exp_sc, exp_ref){
        	exp_ref[which(exp_ref==0)] = delta * min(exp_ref[which(exp_ref > 0)]) 
            log_p_sc_given_ref = dmultinom(x=exp_sc,log=T,prob=exp_ref)
            return(log_p_sc_given_ref)}
        .get_dis= function(this_sc, this_ref, method2=method2){
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
                if(method=='multinomial'){
                    log_p_sc_given_ref=Refprob(exp_sc,exp_ref)
                    } else if(method=='kendall'){
                	log_p_sc_given_ref=cor.fk(exp_sc,exp_ref)
                    } else {
                    log_p_sc_given_ref=cor(exp_sc,exp_ref, method=method)
                    }
                
                log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref) 
                j=j+1
                }
            return(log_p_sc_given_ref_list)
            }
        
        this_tag=as.character(tag[i,2])
        #print(this_tag)
        vec_index=which(as.character(ref_tag[,2])==this_tag)
        this_R=min(random_size, length(vec_index))
        vec_index=sample(vec_index, size= this_R, replace = FALSE)         
        this_vec = ref_vec[vec_index,]
        this_ref= exp_ref_mat[,vec_index]
        this_sc = cbind(exp_sc_mat[,i],exp_sc_mat[,i])
        rownames(this_sc) = rownames(exp_sc_mat)
        colnames(this_sc)= c('rep1','rep2')
        this_out = .get_dis(this_sc, this_ref, method2=method2)
        this_out_rank=rank(-this_out)
        used_index=which(this_out_rank <= nearest_cell)
        
        this_weight = rep(0,length(this_out))
        this_weight[used_index] = 1 #(1-this_out[used_index])/2
        this_weight = this_out_rank * this_weight
        this_weight[used_index] = alpha**this_weight[used_index] 
        this_weight=this_weight/sum(this_weight)

        this_out_vec = t(as.matrix(this_vec)) %*% as.matrix(this_weight)
        this_out_exp='none'
        #this_out_exp = as.matrix(this_ref) %*% as.matrix(this_weight)
        #names(this_out_exp) = rownames(this_ref)
        this_out=list(out_vec=this_out_vec, out_exp=this_out_exp)
        
        if(i%%print_step==1){print(i)}
        return(this_out)
        }
    #windows
    #cl= makeCluster(CPU)
    #RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    #unix
    RUN = mclapply(1:length(exp_sc_mat[1,]), SINGLE, mc.cores=CPU)

    OUT_VEC = c()
    #OUT_EXP = c()
    for(this_out in RUN){
        OUT_VEC = cbind(OUT_VEC, this_out$out_vec)
        #OUT_EXP = cbind(OUT_EXP, this_out$out_exp)
        }
    OUT_VEC = t(OUT_VEC)
    rownames(OUT_VEC) = sc_cell_name
    colnames(OUT_VEC) = colnames(ref_vec) 
    #rownames(OUT_EXP) = names(this_out$out_exp)    
    #colnames(OUT_EXP) = sc_cell_name
    #OUT=list(vec=OUT_VEC, exp=OUT_EXP, tag=tag)
    OUT=list(vec=OUT_VEC, tag=tag)
    return(OUT)
    }


SCREF <- function(exp_sc_mat, exp_ref_mat, method1='kendall', method2='multinomial', min_cell=10, CPU=4, print_step=10){
    print('First-round annotation:')
    print(method1)
    if(method1!='multinomail'){
        out1=.get_cor(exp_sc_mat, exp_ref_mat, method=method1,CPU=CPU, print_step=print_step)
        } else {
        out1=.get_log_p_sc_given_ref(exp_sc_mat, exp_ref_mat, CPU=CPU, print_step=print_step)
        }
    tag1=.get_tag_max(out1)

    print('Build local reference')
    LocalRef=.generate_ref(exp_sc_mat, tag1, min_cell=min_cell)

    print('Second-round annotation:')
    print(method2)
    if(method2!='multinomial'){
        out2=.get_cor(exp_sc_mat, LocalRef, method=method2,CPU=CPU, print_step=print_step)
        } else {
        out2=.get_log_p_sc_given_ref(exp_sc_mat, LocalRef, CPU=CPU, print_step=print_step)
        }
    tag2=.get_tag_max(out2)
    print('Finish!')
    return(tag2)
    }



.trajectory = function(sim_mat){

    library(MASS)
    library(ggplot2)

    input_value=sim_mat
    target_size=3
    label_size=3
    random_seed=123
    random_ratio=0.05
    plot_size=1.5

    set.seed(random_seed)

    CN=length(colnames(input_value))
    N=length(rownames(input_value))

    A=pi/2
    S=2*pi/N
    VEC=c()
    i=1
    while(i<=N){
        A=A+S
        x=cos(A)
        y=sin(A)
        VEC=cbind(VEC,c(x,y))
        i=i+1
    }
    VEC=t(VEC)
    rownames(VEC)=rownames(input_value)
    colnames(VEC)=c('X','Y')

    tmp = apply(input_value,2,scale)
    tmp[which(tmp<0)]=0
    tmp = apply(tmp,2,pnorm)
    tmp[which(tmp==0.5)]=0
    colnames(tmp)=colnames(input_value)
    rownames(tmp)=c(rownames(input_value))
    tmp=t(tmp)
    this_vec=tmp %*% VEC 

    r_this_vec=this_vec

    r_this_vec[,1]=this_vec[,1]+rnorm(CN)*random_ratio
    r_this_vec[,2]=this_vec[,2]+rnorm(CN)*random_ratio

    df=data.frame(r_this_vec); colnames(df) = c("x","y")
    target_df=data.frame(VEC); colnames(target_df) = c("x","y")

    output=ggplot(data=df,aes(x,y)) +  
      geom_point(colour='grey50') +
      stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black') + 
      scale_fill_continuous(low="green",high="red")  + 
      guides(alpha="none") +
      xlim(-plot_size, plot_size) +
      ylim(-plot_size, plot_size) +
      geom_point(data=target_df, aes(x, y), colour="royalblue",size=target_size) +
      geom_text(data=target_df,aes(label=rownames(target_df)), size=label_size)

    return(output)
    }
