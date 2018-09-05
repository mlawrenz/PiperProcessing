for dir in `ls -d complex*min`; do cd $dir; name=rosetta-${dir%%/*}; cp ../run_rosetta.sh .;  qsub -N ${name}  -cwd -q cpu -pe mpi 10 ./run_rosetta.sh rosetta-clean-${dir%%/*}.pdb; cd ../; done


