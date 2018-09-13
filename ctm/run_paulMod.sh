whA_prot_file=whA_protOnly_Brd4_5T35.pdb
whA_bindingConf_file=prep-5T35-fixbond_AT1_CTM_JQ1bindconf.sdf
whB_prot_file=whB_protOnly_VHL_5T35.pdb
whB_bindingConf_file=prep-5T35-fixbond_AT1_CTM_VHLbindconf.sdf
whA_linkers_file=autoSetupOut_whA_linkers.sdf
whB_linkerAtoms_file=autoSetupOut_whB_attach.sdf
smirksFile=autoSetupOut_smirks.smk



program='/mnt/public/pnovick_public/ctm/alignComps_rotTrans.py'
#program='/home/mlawrenz/paulsComplexGeneratingScript/alignComps_rotTrans/alignComps_rotTrans.py'

python $program --whA_prot_file $whA_prot_file \
--whA_bindingConf_file $whA_bindingConf_file --whB_prot_file $whB_prot_file --whB_bindingConf_file $whB_bindingConf_file \
--whA_linkers_file $whA_linkers_file --whB_linkerAtoms_file $whB_linkerAtoms_file --smirksFile $smirksFile --outMolsFile \
outMolsFile.oeb --outComplexFile outComplexFile.oeb --outProteinFile outProteinFile.oeb    --outPDBdir outPDBDir/ \
--ligsOnly --nproc 12 --maxProtLigClashes 30 --maxProtProtClashes 30 --whAconfs_maxConfs 3000 --whBconfs_rms .01 \


