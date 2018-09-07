file=$1

$SCHRODINGER/run split_structure.py -k  -many_files $file ${file%%.pdb*}.split.mae
for maefile in `ls *${file%%.pdb*}*.mae`; do structconvert -imae $maefile -opdb ${maefile%%.mae*}.pdb; done
grep -v CONECT ${file%%.pdb*}.split_chainD.pdb | grep -v END > tmp
grep V32 *X*pdb >> tmp
mv tmp  ${file%%.pdb*}.split_chainD.pdb

grep -v CONECT ${file%%.pdb*}.split_chainA.pdb | grep -v END > tmp
grep JQ1 *X*pdb >> tmp
mv tmp  ${file%%.pdb*}.split_chainA.pdb

sed '1d' < $file > ${file%%.pdb*}_split.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}_split.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chainD.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chainA.pdb

mv score_only.sc ${file%%.pdb*}_score_only.sc
