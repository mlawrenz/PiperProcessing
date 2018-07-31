#!/bin/bash


#clean
for file in `ls *mae.mae`; do if [ -e ${file%%.mae}.mae ]; then rm $file; fi;
done


for file in `ls prep*model.000.*mae`
do 
minfile=prime_mini_${file%%.mae*}-out.maegz
mcfile=prime_hybridmc_${file%%.mae*}-out.maegz
#proplister -p r_psp_Prime_Energy $minfile -o tmp
#val=`tail -n 1 tmp`
#echo $file $val >> prime_mini_model.000.summary.txt 
#proplister -p r_psp_Prime_Energy $mcfile -o tmp
#val=`tail -n 1 tmp`
#echo $file $val >> prime_hybridmc_model.000.summary.txt 
#proplister -p r_user_CA_RMSD $mcfile -o tmp
#val=`tail -n 1 tmp`
#echo $file $val >> rmsd_model.000.summary.txt
proplister -p s_piper_cluster_size  $mcfile -o tmp
val=`tail -n 1 tmp`
echo $file $val >> cluster_model.000.summary.txt
proplister -p r_piper_charmm_energy $mcfile -o tmp
val=`tail -n 1 tmp`
echo $file $val >> charmm_model.000.summary.txt
done
