for file in `ls model*.000.*min*pdb`
do 
name=${file%%.min*}; 
tar=`ls ../*tar*`
mkdir cluster-members-${name}/;
piper_cluster_members  --output-dir  cluster-members-${name}/ $tar $file; 
done
