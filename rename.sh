for dir in `cat dirs_no_prime`; 
# works for all dirs
do   tmp=${dir%%3_*}; tmp2=${dir##*__};
tmp3=${tmp2%%_rmsd-*}; cd $dir/; for file in `ls *mae*`; do
file_base=${file%%-out*}; mv $file
${tmp}_${tmp3}_calc_prime_mini_${file_base##*complex-}.mae; done; cd ../; done^C

