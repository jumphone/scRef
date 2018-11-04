source('scRef.R')
exp_sc_mat=read.table('Zeisel_exp_sc_mat.txt',header=T,row.names=1,'\t',check.name=F)
exp_ref_mat=read.table('Reference_expression.txt',header=T,row.names=1,sep='\t',check.name=F)

out=.get_cor(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10)
tag=.get_tag_max(out)
cat('finished!!!')
write.table(tag,file='TAG.txt',quote=F,row.names=F,col.names=T,sep='\t')


