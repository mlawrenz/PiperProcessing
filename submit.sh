system=$1


for file in `ls prep*mae`
do 
sed "s/XXX/${file%%.mae*}/g" < prime_mini_${system}.inp > prime_mini_${file%%.mae*}.inp
sed "s/XXX/${file%%.mae*}/g" < prime_hybridmc_${system}.inp > prime_hybridmc_${file%%.mae*}.inp
done

total=`ls prep*mae | wc -l`
for file in `ls prime_mini_prep-complex-model.000.*`
do
prime -HOST cpu_only:1 $file
done

for file in `ls prime_hybridmc_prep-complex-model.000.*`
do
prime -HOST cpu_only:1 $file
done
