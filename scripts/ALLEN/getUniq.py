import sys

fa=open('rows_metadata.csv')
fa.readline()
fi=open('expression_matrix.csv')

fo=open('exp_mat.txt','w')
line1=fa.readline().replace('"','')
line2=fi.readline().replace('"','')
old=set()
while line1 !='':
    seq1=line1.rstrip().split(',')
    seq2=line2.rstrip().split(',')
    if seq1[3] not in old:
        old.add(seq1[3])
        fo.write(seq1[3]+'\t'+'\t'.join(seq2[1:])+'\n')
    line1=fa.readline()
    line2=fi.readline() 




