
for file in `ls ../complex-model.000.*mae`; do base=`basename $file`; mkdir ${base%%.mae*}/; structconvert -imae $file -opdb ${base%%.mae*}/${base%%.mae*}.pdb; done
for file in `ls *.pdb`; do mkdir ${file%%.pdb*}/; mv $file ${file%%.pdb*}/; cd ${file%%.pdb*}/; done

for dir in `ls -d complex*/`
do
cd $dir
file=`ls *pdb`
cp ../prepare_for_rosetta.sh .
./prepare_for_rosetta.sh $file V32 JQ1
cd ../
done
