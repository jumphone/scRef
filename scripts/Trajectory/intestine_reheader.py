fi=open('GSE92332_atlas_UMIcounts.txt')
fo=open('GSE92332_atlas_UMIcounts_reheader.txt','w')
header=fi.readline().rstrip().split('\t')
nheader=[]
for one in header:
    seq=one.split('_')
    newone=seq[-1]+'_'+seq[-2]+'_'+seq[-3]
    nheader.append(newone)

fo.write('\t'.join(nheader)+'\n')
for line in fi:
    fo.write(line)





