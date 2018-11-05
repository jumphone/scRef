import sys
import numpy as np
fi=open(sys.argv[1])
fa=open(sys.argv[2])
fo=open(sys.argv[1]+'.net','w')
#fa.readline()
G={}


fo.write(fi.readline())
for line in fi:
    seq=line.rstrip().split('\t')
    seq[0]=seq[0].upper()
    G[seq[0]]=seq[1:]
    tmp_len=len(seq[1:])



for line in fa:
    seq=line.rstrip().split('\t')
    tag=seq[0]
    tmp_exp=[0]*tmp_len
    i=0
    for g in seq[1:]:    
        g=g.upper()	
        if g in G:
            tmp=[]
            for exp in G[g]:
            	tmp.append(float(exp))
            tmp_exp=np.array(tmp_exp)+np.array(tmp)
            i=i+1
    if i!=0:
    	tmp_exp=tmp_exp/float(i)
    fo.write(tag)
    for one in tmp_exp:
    	fo.write('\t'+str(one))
    fo.write('\n')



    


            
    


