#!/bin/bash

rec_atom="C1"
lig_residue='JQ1'
rec_residue='V32'
chain='A'

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
if [[ $ref == *3MX4* ]];
then
echo "on 3MX4 atom N2"
lig_atom="N2"
fi
if [[ $ref == *5T35* ]];
then
echo "on 5T35 atom NBL"
lig_atom="NBL"
fi
if [[ $ref == *5U2C* ]];
then
echo "on 5U2C atom N1"
lig_atom="N1"
fi
$SCHRODINGER/run $PIPER/processPiper.py --rec_residue $rec_residue --lig_residue $lig_residue --lig_atom $lig_atom --rec_atom $rec_atom  --rmsd-ref $ref --chain ${chain}
cp -r  rmsd-maefiles/ ../../results/${dir%%/}_rmsd-maefiles
cp report.txt ../../results/${dir%%/}_rmsd-maefiles
cp rmsd.txt ../../results/${dir%%/}_rmsd-maefiles
fi
done
cd ../../
done

