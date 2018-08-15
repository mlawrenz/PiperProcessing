type=$1

#for file in `ls *mini*model.${type}*mae* | grep -v calc | grep -v chain`
#do
#if [ ! -e ${file%%.mae*}.split_chainC.mae ]
#then
#echo "need split $dir $file" 
#$SCHRODINGER/run split_structure.py -k  -many_files $file ${file%%.mae*}.split.mae
#fi
#done


for file in `ls *mini*model.${type}*mae* | grep -v calc`
do 
str='maegz'
if [[ $file =~ .maegz ]]
then
if [ ! -e  ${file%%.maegz}.mae ]
then
echo "convert $file"
#structconvert -imae $file -omae ${file%%.maegz}.mae
fi
rm $file
fi
file=${file%%.maegz}.mae
sed "s/XXX/${file%%.mae*}/g" < ../prime_calcenergy.inp > prime_calcenergy_${file%%.mae*}.inp
done




rm prime*calc*prime*calc*

for input in `ls prime*calc*model.${type}*inp`
do 
logfile=`ls ${input%%.inp*}*.log`
if [ -z $logfile ]
then
echo "run PRIME for $dir $input"
prime -HOST cpu_only:1 $input
else
echo "have $logfile"
fi
done
