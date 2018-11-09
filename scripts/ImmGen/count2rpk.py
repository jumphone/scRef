fa=open('MGI_HGNC_homologene.rpt')
L={}
fa.readline()
for line in fa:
    seq=line.rstrip().split('\t')
    try:
        l= abs(float(seq[11])-float(seq[12]))
        L[seq[1]]=l
    except Exception as e:
        pass


fi=open('GSE109125_Gene_count_table.csv.uniq')
fo=open('GSE109125_Gene_count_table.csv.uniq.rpk','w')


fo.write(fi.readline())
for line in fi:
    seq=line.rstrip().split('\t')
    if seq[0] in L:
        l=L[seq[0]]
        fo.write(seq[0])
        for one in seq[1:]:
            rpk=float(one)/l * 1000
            fo.write('\t'+str(rpk))
        fo.write('\n')



