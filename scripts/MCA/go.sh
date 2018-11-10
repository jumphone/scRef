workspace='./MCA_BatchRemove_dge/rmbatch_dge/'

TAG='Bladder'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='BoneMarrow'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
f3=$workspace$TAG\3_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2 $f3
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='BoneMarrowcKit'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
f3=$workspace$TAG\3_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2 $f3
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Brain'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='CJ7.EB.WT'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='CJ7.EB14.Ezh2.1'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='EB.Ezh2'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='EB.WT'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='EmbryonicMesenchymeE14.5'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='EmbryonicStemCells'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalBrain'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalFemaleGonad'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalIntestine'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalKidney'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalLiverE14.1'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalLung'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalMaleGonad'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalPancreas'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='FetalStomach'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Kidney'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Liver'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Lung'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
f3=$workspace$TAG\3_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2 $f3
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MammaryGland.Involution.CD45.'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MammaryGland.Involution'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MammaryGland.Lactation'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MammaryGland.Pregnancy'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MammaryGland.Virgin.CD45.'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MammaryGland.Virgin'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
f3=$workspace$TAG\3_rm.batch_dge.txt
f4=$workspace$TAG\4_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2 $f3 $f4
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MesenchymalStemCells'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='MesenchymalStemCellsPrimary'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Muscle'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='NeonatalCalvaria'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='NeonatalHeart'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='NeonatalMuscle'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='NeonatalPancreas'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='NeonatalRib'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
f3=$workspace$TAG\3_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2 $f3
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='NeonatalSkin'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='NeontalBrain'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Ovary'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Pancreas'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='PeripheralBlood'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
f3=$workspace$TAG\3_rm.batch_dge.txt
f4=$workspace$TAG\4_rm.batch_dge.txt
f5=$workspace$TAG\5_rm.batch_dge.txt
f6=$workspace$TAG\6_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2 $f3 $f4 $f5 $f6
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Prostate'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='SmallIntestine.CD45'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='SmallIntestine'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
f3=$workspace$TAG\3_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2 $f3
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Spleen'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Stomach'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Testis'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Thymus'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='TrophoblastStemCells'
f1=$workspace$TAG\_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt

TAG='Uterus'
f1=$workspace$TAG\1_rm.batch_dge.txt
f2=$workspace$TAG\2_rm.batch_dge.txt
python 1combine.py $workspace$TAG\_com.txt $f1 $f2
python 2getAnno.py $workspace$TAG\_com.txt
python 3compress.py $workspace$TAG\_com.txt.anno.txt ./OUT/$TAG\_ref_mouse.txt
python mouse2human.py ./OUT/$TAG\_ref_mouse.txt ./OUT/$TAG\_ref_human.txt








