input=$1
addside=$2


program=/opt/apps/Schrodinger/default/utilities/prepwizard
if [ -z "$side" ]
then
output=prep-${input%%.mae}.mae
$program -f 3 -fix -noepik -r 0.3 -HOST cpu_only:1 -NOLOCAL $input $output 
else
echo "RUNNING WITH ADD SIDECHAINS"
output=prep-sidech-${input%%.mae}.mae
$program -f 3 -epik_pHt 0.0 -r 0.3 -fillsidechains -HOST cpu_only:1 -NOLOCAL $input $output

fi


