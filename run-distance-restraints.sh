receptor=$1
ligand=$2
restraints=distance_restraints.json
 
run_piper --rec $receptor --lig $ligand --restraint-set $restraints -np 12 --mpi --piper-license ~/piper-1.1.0_b1/piper_package/piper_license  --piper-base ~/piper-1.1.0_b1/piper_package

