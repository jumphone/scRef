import sys
fi=open(sys.argv[1])
fa=open(sys.argv[2])
fo=open(sys.argv[1]+'.net','w')
#fa.readline()
DE={}
MAX=1
for line in fa:
    seq=line.rstrip().split('\t')
    seq[0]=seq[0].upper()
    try:
        DE[seq[0]]=float(seq[1])
        MAX=max(MAX,float(seq[1]))
    except Exception as e:
        pass

fo.write(fi.readline())
for line in fi:
    seq=line.rstrip().split('\t')
    seq[0]=seq[0].upper()
    print(seq[0])
    if seq[0] in DE:
         this_de= DE[seq[0]] 
    else:
         this_de = 1.0
    j=0
    while j< this_de:
         fo.write(seq[0]+'_'+str(j)+'\t')
         fo.write('\t'.join(seq[1:])+'\n')
         #i=1
         #while i<len(seq):
         #    fo.write('\t'+str(this_de * float(seq[i])))
         #    i=i+1
         #fo.write('\n')
         j=j+1
    #else:
    #     fo.write(line)

