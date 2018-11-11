workspace='./FACS/FACS/'

TAG='Aorta'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Bladder'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Brain_Myeloid'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Brain_Non-Myeloid'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Diaphragm'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Fat'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Heart'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Kidney'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Large_Intestine'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Limb_Muscle'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Liver'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Lung'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Mammary_Gland'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Marrow'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Pancreas'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Skin'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Spleen'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Thymus'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Tongue'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Trachea'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1getAnno.py $workspace$TAG\-counts.csv
python 2compress.py $workspace$TAG\-counts.csv.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt


