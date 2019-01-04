import sys
output=sys.argv[1]
input_list=(sys.argv[2:])


EXP={}
header=[]
for input_file in input_list:
    fi=open(input_file)
    header=header+fi.readline().replace('"','').rstrip().split('\t')
    for line in fi:
        seq=line.replace('"','').rstrip().split('\t')
        if seq[0] in EXP:
            EXP[seq[0]]=EXP[seq[0]]+seq[1:]
        else:
        	EXP[seq[0]]=seq[1:]
    fi.close()
print(header)
fo=open(output,'w')
fo.write('\t'.join(header)+'\n')
for gene in EXP:
    #print(len(EXP[gene]))
    #print(len(header))
    if len(EXP[gene])==len(header):
	#print('ok')
        fo.write(gene+'\t'+'\t'.join(EXP[gene])+'\n')

fo.close()








    



