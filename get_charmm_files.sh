
for x in `seq 0 2`
do
for dir in `ls -d *rmsd*/`; do tmp=${dir##*./}; base=${tmp%%_rmsd*};  cd $dir; for file in `ls -v model.00${x}.*`; do name=${file%%-rmsd*}; num=`proplister -p r_piper_charmm_energy $file |  grep -v piper | grep -v "\-\-"` ; echo  "$name $num" >> ../${base}_model${x}_charmm.txt ; done; cd ../; done
done

