import sys

fa=open('annotations_facs.csv')
name2anno={}
fa.readline()
for line in fa:
    seq=line.rstrip().split(',')

    name=seq[2]
    #print(name)
    anno=seq[-3]+'_'+seq[3]+'_'+seq[-4]+'_'+seq[6]
    name2anno[name]=anno


fi=open(sys.argv[1])
header=fi.readline().rstrip().replace('"','').split(',')[1:]
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
    seq=line.replace('"','').rstrip().split(',')
    fo.write(seq[0])
    for one in used:
            fo.write('\t'+seq[one+1])
    fo.write('\n')






