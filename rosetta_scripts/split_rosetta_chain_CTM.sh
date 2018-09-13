file=$1
chain1=$2
chain2=$3
ctm=$4

ctm_chain=`grep $ctm $file | awk '{print $5}' | head -10 | tail -n 1`

$SCHRODINGER/run split_structure.py -k  -many_files $file ${file%%.pdb*}.split.mae
for maefile in `ls *${file%%.pdb*}*.mae`; do structconvert -imae $maefile -opdb ${maefile%%.mae*}.pdb; done
./run_rosetta_score_only.sh ${file%%.pdb*}_split.pdb ${ctm}.params
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chain${chain1}.pdb ${ctm}.params
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chain${chain2}.pdb ${ctm}.params
./run_rosetta_score_only.sh ${file%%.pdb*}.split_chain${ctm_chain}.pdb ${ctm}.params

mv score_only.sc ${file%%.pdb*}_score_only.sc
