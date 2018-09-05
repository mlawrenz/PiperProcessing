structure=$1
lig1=JQ1.params
lig2=V32.params

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/opt/apps/Rosetta/main/source/build/external/release/linux/3.10/64/x86/gcc/4.8/default/

/opt/apps/Rosetta/main/source/bin/score_jd2.linuxgccrelease  -s $structure -out:file:scorefile score_only.sc -ex1 -ex2 -use_input_sc -extra_res_fa $lig1 $lig2  >& rosetta_dock_score_only_${structure%%.pdb*}.log



























