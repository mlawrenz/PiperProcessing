ref=-87.420

for x in `seq 1 150`
do
dir=x${x}/
cd $dir
high=`awk '{print $2, $23}' score.sc  | sort -nk1 | grep complex  | awk '{print $1}' | tail -n 1`
low=`awk '{print $2, $23}' score.sc  | sort -nk1 | grep complex  | awk '{print $1}' | head -1`
low_file=`awk '{print $2, $23}' score.sc  | sort -nk1 | grep complex |  awk '{print $2}'  | head -1`

avg=`awk '{print $2, $23}' score.sc  | sort -nk1 | grep complex  | awk '{sum+=$1} END {print sum/NR}'`
echo $dir $low $low_file $high $avg
mkdir low_score
cp ${low_file}.pdb low_score
cd ../
done

for file in `ls -v clean-complex-model.000.*pdb`; do grep ${file%%.pdb*}
score.sc  | sort -nk 2 | awk '{print $23, $2}' | head -1 >>
rosetta_model.000.summary.txt; done

