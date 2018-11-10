import sys

fa=open('MCA_CellAssignments.csv')
name2anno={}
fa.readline()
for line in fa:
    seq=line.rstrip().split(',')

    name=seq[1]
    #print(name)
    anno=seq[6]
    name2anno[name]=anno



fi=open(sys.argv[1])
header=fi.readline().rstrip().replace('"','').split()
new_header=[]
used=[]
i=0
while i< len(header):
    if header[i] in name2anno:
        new_header.append(name2anno[header[i]])
        used.append(i)
    i=i+1

fo=open(sys.argv[1]+'.anno.txt','w')
fo.write('gene\t'+'\t'.join(new_header)+'\n')
for line in fi:
    seq=line.replace('"','').rstrip().split()
    fo.write(seq[0])
    for one in used:
            fo.write('\t'+seq[one+1])
    fo.write('\n')






