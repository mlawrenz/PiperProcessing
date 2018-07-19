
for file in `ls *pdb`
do 
cat $file | grep -v ANISOU | grep -v REMARK | grep -v TITLE | grep -v HEADER | grep -v CRYS | grep -v EXPD | grep -v CONECT > test.pdb
mv test.pdb clean-$file
done
