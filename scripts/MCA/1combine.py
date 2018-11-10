import sys
output=sys.argv[1]
input_list=(sys.argv[2:])


EXP={}
header=[]
for input_file in input_list:
    fi=open(input_file)
    header=header+fi.readline().replace('"','').rstrip().split()
    for line in fi:
        seq=line.replace('"','').rstrip().split()
        if seq[0] in EXP:
            EXP[seq[0]]=EXP[seq[0]]+seq[1:]
        else:
        	EXP[seq[0]]=seq[1:]
    fi.close()

fo=open(output,'w')
fo.write('\t'.join(header)+'\n')
for gene in EXP:
    if len(EXP[gene])==len(header):
        fo.write(gene+'\t'+'\t'.join(EXP[gene])+'\n')

fo.close()








    



