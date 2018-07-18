export SCHRODINGER=/opt/apps/Schrodinger/suite2017-4
export name=MD_model03_2ns

cd ~/$name

$SCHRODINGER/run ~mlawrenz/AnalyzeMD/ExtractMD.py -c ${name}-out.cms -t ${name}_trj/ -a "protein or ligand"

cd analysis

export reference=solute-${name}_reference.pdb
export dcd=aln-solute-${name}.dcd

$SCHRODINGER/run ~mlawrenz/AnalyzeMD/AnalyzeMD.py -r $reference -t $dcd rmsd_calc rmsd-all

$SCHRODINGER/run ~mlawrenz/AnalyzeMD/AnalyzeMD.py -r solute-MD_model03_2ns_reference.pdb -t aln-solute-MD_model03_2ns.dcd hbonds --selection '*.B' '*.A'  -o interchain
