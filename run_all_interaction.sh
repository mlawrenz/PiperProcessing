
for input in `ls prep*complex*mae`
do
if [ ! -e ${input%%.mae*}.interaction.csv ]
then
$SCHRODINGER/run protein_interaction_analysis.py -structure $input -group1 C -group2 D -outfile ${input%%.mae*}.interaction.csv -hbond_max_dist 4.0 -allowable_vdw_overlap 3.0 -salt_bridge_max_dist 5.0 -max_stack_dist 5.0
echo "done $input"
fi
done

