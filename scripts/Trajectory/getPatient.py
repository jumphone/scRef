import sys
fi=open('GSE70630_OG_processed_data_v2.txt')
fo=open('GSE70630_OG_processed_data_v2_MGH54.txt','w')
#fo=open('GSE70630_OG_processed_data_v2_MGH36.txt','w')
#fo=open('GSE70630_OG_processed_data_v2_MGH53.txt','w')
#fo=open('GSE70630_OG_processed_data_v2_MGH60.txt','w')
#fo=open('GSE70630_OG_processed_data_v2_MGH93.txt','w')
#fo=open('GSE70630_OG_processed_data_v2_MGH97.txt','w')
header=fi.readline().rstrip().split('\t')

tmp=[]
i=0
while i<len(header):
    #print(header[i])
    if 'MGH54' in header[i] or header[i][:2]=='54':
    #if 'MGH36' in header[i] or header[i][:2]=='36':
    #if 'MGH53' in header[i] or header[i][:2]=='53':
    #if 'MGH60' in header[i] or header[i][:2]=='60':
    #if 'MGH93' in header[i] or header[i][:2]=='93':
    #if 'MGH97' in header[i] or header[i][:2]=='97':
        tmp.append(i)
    i=i+1


new_header=['gene']
for one in tmp:
    #print header[one]
    if header[one][0]!='M':
	    header[one]='MGH'+header[one] 
    new_header.append(header[one])
fo.write('\t'.join(new_header)+'\n')

for line in fi:
    seq=line.rstrip().split('\t')

    fo.write(seq[0].replace("'",''))
    for one in tmp:
        try:
            this_exp=str(float(seq[one]))
        except Exception as e:
            this_exp='0.0'
        if this_exp=='nan':
            this_exp='0.0'
        fo.write('\t'+this_exp)
    fo.write('\n')





