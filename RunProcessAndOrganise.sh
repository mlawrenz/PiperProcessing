#!/bin/bash

mkdir results

for dir in `ls -d */ | grep -v results | grep -v test`
do 
cd $dir/
echo ${dir%%/}_work
for subdir in `ls -d */ | grep -v piper`
do
if [ "$(ls -A $subdir)" ]
then
cd $subdir/
name=${dir%%__*}
ref=`ls ../../../Brd4_aligned_proteins/clean*${name}*pdb`
echo $ref
$SCHRODINGER/run $PIPER/processPiperNew.py --rmsd-ref $ref --chain A
cp -r  rmsd-maefiles/ ../../results/${dir%%/}_rmsd-maefiles
cp report.txt ../../results/${dir%%/}_rmsd-maefiles
cp rmsd.txt ../../results/${dir%%/}_rmsd-maefiles
fi
done
cd ../../
done

#cd results/
#for file in `find . -name "*rmsd*"txt`; do cat $file | sort -nk2 | head; echo "---------------------"; echo $file; echo "------------------------"; done
