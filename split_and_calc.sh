type=000

#split prime minimized 

for file in `ls *mini*model.${type}*mae* | grep -v calc | grep -v chain`
do
if [ ! -e ${file%%.mae*}.split_chainC.mae ]
then
echo "need split $dir $file" 
$SCHRODINGER/run split_structure.py -k  -many_files $file ${file%%.mae*}.split.mae
fi
done


for file in `ls *mini*model.${type}*mae* | grep -v calc | grep chain`
do 
file=${file%%.mae*}.mae
sed "s/XXX/${file%%.mae*}/g" < $PIPER/prime_calcenergy.inp > prime_calcenergy_${file%%.mae*}.inp
done




for input in `ls prime*calc*model.${type}*inp`
do 
prime -HOST cpu_only:1 $input
done
