import sys
fi=open(sys.argv[1])
fo=open(sys.argv[1]+'.degree','w')
fo.write('GENE\tDEGREE\n')

DE={}

for line in fi:
    seq=line.rstrip().split('\t') 
    if seq[0] in DE:
        DE[seq[0]]+=1
    else:
        DE[seq[0]]=1
    if seq[1] in DE:
        DE[seq[1]]+=1
    else:
        DE[seq[1]]=1

for one in DE:
    fo.write(one+'\t'+str(DE[one])+'\n')
    


