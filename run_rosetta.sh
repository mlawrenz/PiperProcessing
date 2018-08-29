structure=$1
lig1=JQ1.params
lig2=v32.params

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/apps/Rosetta/main/source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/

time `/opt/apps/Rosetta/main/source/bin/relax.default.linuxgccrelease -database /opt/apps/Rosetta/main/database/  -relax:constrain_relax_to_start_coords -relax:coord_constrain_sidechains -relax:ramp_constraints false -ex1 -ex2 -use_input_sc -flip_HNQ -no_optH false -s $structure -extra_res_fa $lig1 $lig2 -j 10 -nstruct 10 >& rosetta_${structure%%.pdb*}.log` 


























