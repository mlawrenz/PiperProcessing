ref=$1

for file in `ls prime*out*mae*`; do structconvert -imae $file -opdb ${file%%.mae}.pdb; done
mkdir rmsd-pdbfiles
mv prime*pdb rmsd-pdbfiles

cd rmsd-pdbfiles
ls -v *pdb >> list
$SCHRODINGER/run $PIPER/calculate_piper_rmsd.py -r ../$ref -l list --postprocess
cd ../
$SCHRODINGER/run $PIPER/add_rmsd_property.py --rmsd rmsd-pdbfiles/rmsd.txt

