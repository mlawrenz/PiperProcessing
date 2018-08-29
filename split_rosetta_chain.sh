file=$1

$SCHRODINGER/run split_structure.py -k  -many_files $file ${file%%.pdb*}.split.mae
for maefile in `ls *${file%%.pdb*}*.mae`; do structconvert -imae $maefile -opdb ${maefile%%.mae*}.pdb; done
grep -v CONECT ${file%%.pdb*}.split_chainD.pdb | grep -v END > tmp
grep v32 *X*pdb >> tmp
mv tmp  ${file%%.pdb*}.split_chainD.pdb

grep -v CONECT ${file%%.pdb*}.split_chainA.pdb | grep -v END > tmp
grep JQ1 *X*pdb >> tmp
mv tmp  ${file%%.pdb*}.split_chainA.pdb

sed '1d' < $file > ${file%%.pdb*}_split.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}_split.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chainD.pdb
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chainA.pdb

#grep SCORE score_only.sc  | awk '{print $2, $4, $5, $10, $11, $12, $13, $14, $26}' >> rosetta_score_model.000.summary.txt


