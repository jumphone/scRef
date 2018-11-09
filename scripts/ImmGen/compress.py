fi=open('GSE109125_Gene_count_table.csv.uniq.rpk')
fo=open('GSE109125_Gene_count_table.csv.uniq.rpk.pure','w')
#fi.readline()
header=fi.readline().rstrip().split('\t')
h=[]
for one in header:
    h.append(one.split('#')[0])
header=h

G={}
i=1
while i<len(header):
    if header[i] in G:
        G[header[i]].append(i)
    else:
        G[header[i]]=[i]
    i=i+1


new_header=[]
for one in G:
    new_header.append(one)
fo.write('\t'.join(new_header)+'\n')
for line in fi:
    seq=line.rstrip().split('\t') 
    fo.write(seq[0])
    for one in new_header:
        tmp=[]
        for exp in G[one]:
            tmp.append(float(seq[exp]))
        tmp=sum(tmp)/float(len(tmp))
        fo.write('\t'+str(tmp))
    fo.write('\n')
