import sys
fi=open(sys.argv[1])

#fi=open('GSE89567_IDH_A_processed_data.txt')

header=fi.readline().upper().replace('-','_').rstrip().split('\t')
P=set()
tmp={}
i=1
while i<len(header):

    p=header[i].split('_')[0]
    if p[0] in ['0','1','2','3','4','5','6','7','8','9']:
        p='MGH'+p
        header[i]='MGH'+header[i]

        
    P.add(p)
    if p in tmp:
        tmp[p].append(i)
    else:
        tmp[p]=[i]
    i=i+1

FO={}
for p in P:
    FO[p]='./Patients/'+p+'_mat.txt'
    new_header=['gene']
    for one in tmp[p]:
        new_header.append(header[one])
    open(FO[p],'w').write('\t'.join(new_header)+'\n')


for line in fi:
    seq=line.rstrip().split('\t')

    for p in P:
        fo=open(FO[p],'a')
        fo.write(seq[0].rstrip().lstrip())
        for one in tmp[p]:
            try:
                this_exp=str(float(seq[one]))
            except Exception as e:
                this_exp='0.0'
            if this_exp=='nan':
                this_exp='0.0'
            fo.write('\t'+this_exp)
        fo.write('\n')




