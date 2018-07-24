#!/bin/bash

mkdir results

for dir in `ls -d */ | grep -v results`
do 
cd $dir/
echo ${dir%%/}_work
for subdir in `ls -d */ | grep -v piper`
do
if [ "$(ls -A $subdir)" ]
then
cd $subdir/
#mv maefiles/ ${dir%%/}_maefiles
#mv rmsd-maefiles/ ${dir%%/}_rmsd-maefiles
#mv ${dir%%/}_maefiles ../../results/
#mv ${dir%%/}_rmsd-maefiles ../../results
cp rmsd.txt ../../results/${dir%%/}_rmsd-maefiles

fi
done
cd ../../
done


