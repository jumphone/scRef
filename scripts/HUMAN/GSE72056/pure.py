fi=open('GSE72056_melanoma_single_cell_revised_v2.txt')
header=fi.readline()
fi.readline()

fo=open('GSE72056_cnv.txt','w')
fo.write(header)
fo.write(fi.readline().replace('2','malignant').replace('1','non-malignant').replace('0','unresolved'))

fo=open('GSE72056_tag.txt','w')
fo.write(header)
fo.write(fi.readline().replace('1','T-cells').replace('2','B-cells').replace('3','Macrophages').replace('4','Endothelial').replace('5','CAFs'))

fo=open('GSE72056_melanoma_single_cell_revised_v2.txt.pure','w')
old=set()
for line in fi:
    seq=line.rstrip('\t')
    if seq[0] not in old:
    	fo.write(line)
    	old.add(seq[0])
