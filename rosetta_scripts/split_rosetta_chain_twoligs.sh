file=$1
chain1=$2
lig1=$3
chain2=$4
lig2=$5

$SCHRODINGER/run split_structure.py -k  -many_files $file ${file%%.pdb*}.split.mae
for maefile in `ls *${file%%.pdb*}*.mae`; do structconvert -imae $maefile -opdb ${maefile%%.mae*}.pdb; done
grep -v CONECT ${file%%.pdb*}.split_chain${chain1}.pdb | grep -v END > tmp
grep $lig1 *X*pdb >> tmp
mv tmp  ${file%%.pdb*}.split_chain${chain1}.pdb

grep -v CONECT ${file%%.pdb*}.split_chain${chain2}.pdb | grep -v END > tmp
grep $lig2 *X*pdb >> tmp
mv tmp  ${file%%.pdb*}.split_chain${chain2}.pdb

sed '1d' < $file > ${file%%.pdb*}_split.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}_split.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chain${chain1}.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chain${chain2}.pdb

mv score_only.sc ${file%%.pdb*}_score_only.sc
