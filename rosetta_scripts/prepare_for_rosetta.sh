file=$1
lig1=$2
lig2=$3
cat $file | grep -v ANISOU | grep -v REMARK | grep -v TITLE | grep -v HEADER | grep -v CRYS | grep -v EXPD | grep -v CONECT > test.pdb
mv test.pdb clean-$file

grep "$lig1" clean-$file  > ${lig1}.pdb
babel  ${lig1}.pdb  ${lig1}.mol2
if [ ! -z $lig2 ]
then
grep "$lig2" clean-$file  > ${lig2}.pdb
babel  ${lig2}.pdb  ${lig2}.mol2
fi

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/apps/Rosetta/main/source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/

/opt/apps/Rosetta/main//source/scripts/python/public/molfile_to_params.py -n ${lig1}  -p ${lig1} --conformers-in-one-file ${lig1}.mol2
if [ ! -z $lig2 ]
then
/opt/apps/Rosetta/main//source/scripts/python/public/molfile_to_params.py -n ${lig2}  -p ${lig2} --conformers-in-one-file ${lig2}.mol2
fi


if [ ! -z $lig2 ]
then
grep -v ${lig1} clean-$file | grep -v TER | grep -v END | grep -v ${lig2} > rosetta-clean-$file
cat ${lig1}.pdb >>  rosetta-clean-$file
cat ${lig2}.pdb >>  rosetta-clean-$file
else
grep -v ${lig1} clean-$file | grep -v TER | grep -v END > rosetta-clean-$file
cat ${lig1}.pdb >>  rosetta-clean-$file
fi






